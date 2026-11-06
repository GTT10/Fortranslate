program test_regular_held_godunov_3d
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_mass_properties
  use gas_transport_mod, only: gas_transport_species
  use transport_database_mod, only: load_h2o2_full_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use dual_energy_mod, only: dual_energy_policy
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_sdc_3d,advance_amr_patch_tree_sdc_subcycle_3d
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_sdc_3d,advance_mpi_amr_tree_sdc_subcycle_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d
  use amr_patch_tree_3d_mod, only: initialize_amr_patch_tree_3d,composite_integrals_amr_patch_tree_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_sdc_amr_patch_tree_hydro_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: compute_sdc_mpi_amr_tree_hydro_3d
  implicit none
  type :: result
    type(amr_tree_state_3d) :: gas,held,work
    type(amr_tree_flux_3d), allocatable :: flux(:),aux(:)
  end type
  type(result) :: serial,distributed(2),saved
  type :: native_state
    type(amr_tree_state_3d) :: gas,temp,held
    type(amr_tree_flux_3d), allocatable :: flux(:),aux(:),physical(:),physical_aux(:)
    real(dp), allocatable :: boundary(:),source(:)
    real(dp) :: prefix,correction,theta,work_change
  end type
  type(native_state) :: native_serial,native_mpi,native_saved
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction) :: reactions(0)
  type(dual_energy_policy) :: policy
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: boxes(1),empty(0)
  type(amr_tree_state_3d) :: gas(2),temp,held,source,local_gas(2),local_temp,local_held,local_source
  type(amr_tree_state_3d) :: packet,future,local_packet,local_future,first_rate
  type(nasa7_species), allocatable :: species(:)
  integer, allocatable :: owners(:)
  real(dp), allocatable :: y(:)
  real(dp) :: spacing(3),dt,viscosity,hh(3),lhs,scale,error,nonzero,mach(6),tin(6),pressure(6)
  integer :: rank,nrank,status,nv,kind,model,a,node,n(3),i,j,k,c,pass
  character(len=16) :: faces(6)
  logical :: good
  call MPI_Init(status);call MPI_Comm_rank(MPI_COMM_WORLD,rank,status);call MPI_Comm_size(MPI_COMM_WORLD,nrank,status)
  call require(nrank==3,'two owners and idle')
  call load_h2o2_full_thermo(species,good);call require(good,'thermo')
  call load_h2o2_full_transport(transport,good);call require(good,'transport')
  policy%eta2=2;policy%allow_negative_energy=.true.;policy%update_total_from_internal=.false.
  ! Select the pinned low-energy-floor branch that preserves admissible held
  ! H when eta2 refuses the total-energy projection, including formation H<0.
  policy%allow_small_energy=.false.
  nv=5+size(species);allocate(y(size(species)));y=0
  do c=1,size(species)
    if(species(c)%name=='N2') y(c)=1
  enddo
  spacing=.001_dp;dt=2e-8_dp
  do kind=1,3
    if(kind==1.or.kind==3) then
      call initialize_amr_patch_tree_3d([4,4,4],[integer ::],empty,tree,good);owners=[0]
    else
      call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),good)
      call require(good,'fine patch')
      call initialize_amr_patch_tree_3d([4,4,4],[1],boxes,tree,good);owners=[0,1]
    endif
    call require(good,'tree')
    if(kind==3) then
      y=0
      do c=1,size(species)
        if(species(c)%name=='H2O') y(c)=1
      enddo
    endif
    call initialize()
    if(kind==3) call require(all(held%fields(1)%values<0),'ordinary signed formation energy exercised')
    do a=1,2
      call localize(gas(a),local_gas(a))
    enddo
    call localize(temp,local_temp);call localize(held,local_held);call localize(source,local_source)
    do model=1,2
      viscosity=merge(0._dp,.1_dp,model==1)
      faces='periodic'
      do pass=1,3
        faces='periodic';mach=-1;tin=600;pressure=0
        if(pass==2) faces(1:2)='slip_wall'
        if(pass==3) then
          faces(1)='inflow';faces(2)='outflow';mach(1)=.2_dp;pressure(2)=80000
          ! Keep the prescribed pressure close to this mixture's interior
          ! pressure: all four characteristic halo layers must remain valid.
          if(kind==3) pressure(2)=125000
        endif
        call compute_sdc_amr_patch_tree_hydro_3d(tree,species,gas(1),temp,source,spacing,dt,faces, &
          serial%gas,serial%flux,good,difmag=viscosity,inflow_mach=mach,inflow_temperature=tin, &
          outlet_pressure=pressure,internal_energy=held, &
          internal_energy_rate=serial%held,internal_energy_integrated_flux=serial%aux, &
          internal_energy_pressure_work=serial%work)
        if(.not.good) print *, 'Serial spatial failure kind/model/pass:',kind,model,pass
        call require(good,'serial Godunov actualH spatial')
        do a=1,2
          call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local_gas(a),local_temp,local_source,spacing,dt,faces, &
            distributed(a)%gas,distributed(a)%flux,MPI_COMM_WORLD,good,difmag=viscosity, &
            inflow_mach=mach,inflow_temperature=tin,outlet_pressure=pressure,internal_energy=local_held, &
            internal_energy_rate=distributed(a)%held,internal_energy_integrated_flux=distributed(a)%aux, &
            internal_energy_pressure_work=distributed(a)%work)
          call require(good,'MPI Godunov actualH spatial')
        enddo
        nonzero=0
        do node=1,size(owners)
          if(owners(node)/=rank) then
            call require(.not.allocated(distributed(1)%held%fields(node)%values),'owner-local held output')
            cycle
          endif
          call compare_node(serial,distributed(1),node,'serial MPI exact paired family')
          if(any(distributed(1)%held%fields(node)%values/=distributed(2)%held%fields(node)%values)) &
            print *, 'Held mismatch kind/model/pass/node/max:',kind,model,pass,node, &
              maxval(abs(distributed(1)%held%fields(node)%values-distributed(2)%held%fields(node)%values))
          call require(all(distributed(1)%held%fields(node)%values==distributed(2)%held%fields(node)%values), &
            'actualH residual independent of bulk energy')
          call require(all(distributed(1)%work%fields(node)%values==distributed(2)%work%fields(node)%values), &
            'actualH pressure work independent of bulk energy')
          call compare_faces(distributed(1)%aux(node),distributed(2)%aux(node),'actualH face independence')
          n=tree%nodes(node)%shape;hh=spacing/real(2**tree%nodes(node)%level,dp)
          error=0;scale=1
          do k=1,n(3)
          do j=1,n(2)
          do i=1,n(1)
            lhs=distributed(1)%held%fields(node)%values(1,i,j,k)+distributed(1)%work%fields(node)%values(1,i,j,k)
            lhs=lhs+(distributed(1)%aux(node)%x(1,i+1,j,k)-distributed(1)%aux(node)%x(1,i,j,k))/(dt*hh(1))+ &
              (distributed(1)%aux(node)%y(1,i,j+1,k)-distributed(1)%aux(node)%y(1,i,j,k))/(dt*hh(2))+ &
              (distributed(1)%aux(node)%z(1,i,j,k+1)-distributed(1)%aux(node)%z(1,i,j,k))/(dt*hh(3))
            error=max(error,abs(lhs));scale=max(scale,abs(distributed(1)%held%fields(node)%values(1,i,j,k)))
          enddo
          enddo
          enddo
          call require(error<=128*epsilon(1._dp)*scale,'UEINT residual uses paired face divergence minus pdivu')
          nonzero=max(nonzero,maxval(abs(distributed(1)%work%fields(node)%values)))
        enddo
        call MPI_Allreduce(MPI_IN_PLACE,nonzero,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
        call require(nonzero>1,'nonzero pressure work exercised')
      enddo
    enddo
    saved=distributed(1)
    if(rank==owners(size(owners))) local_held%fields(size(owners))%values(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
    call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local_gas(1),local_temp,local_source,spacing,dt,faces, &
      distributed(1)%gas,distributed(1)%flux,MPI_COMM_WORLD,good,internal_energy=local_held, &
      internal_energy_rate=distributed(1)%held,internal_energy_integrated_flux=distributed(1)%aux, &
      internal_energy_pressure_work=distributed(1)%work)
    call require(.not.good,'bad actualH refuses collectively')
    do node=1,size(owners)
      if(owners(node)==rank) call compare_node(saved,distributed(1),node,'all spatial outputs rollback')
    enddo
    call localize(held,local_held)
    if(rank==2) then
      call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local_gas(1),local_temp,local_source,spacing,dt,faces, &
        distributed(1)%gas,distributed(1)%flux,MPI_COMM_WORLD,good,internal_energy=local_held, &
        internal_energy_rate=distributed(1)%held,internal_energy_integrated_flux=distributed(1)%aux)
    else
      call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local_gas(1),local_temp,local_source,spacing,dt,faces, &
        distributed(1)%gas,distributed(1)%flux,MPI_COMM_WORLD,good,internal_energy=local_held, &
        internal_energy_rate=distributed(1)%held,internal_energy_integrated_flux=distributed(1)%aux, &
        internal_energy_pressure_work=distributed(1)%work)
    endif
    call require(.not.good,'pressure-work presence mismatch refuses collectively')
    do node=1,size(owners)
      if(owners(node)==rank) call compare_node(saved,distributed(1),node,'presence rollback')
    enddo
    if(kind==2) then
      ! A native fine interval consumes full gas/T/H parent endpoints.
      ! Alter only future H to expose accidental projection from gas E.
      call pack(packet);future=packet;future%fields(1)%values(nv+2,:,:,:)=1.001_dp*held%fields(1)%values(1,:,:,:)
      call localize(packet,local_packet);call localize(future,local_future)
      faces='periodic'
      do pass=1,2
        call compute_sdc_amr_patch_tree_hydro_3d(tree,species,gas(1),temp,source,spacing,dt,faces, &
          serial%gas,serial%flux,good,target_depth=1,coarse_start=packet,coarse_end=future, &
          coarse_fractions=[.5_dp*real(pass-1,dp)],source_start=source,source_end=source, &
          source_fractions=[.5_dp*real(pass-1,dp)],internal_energy=held, &
          internal_energy_rate=serial%held,internal_energy_integrated_flux=serial%aux, &
          internal_energy_pressure_work=serial%work)
        call require(good,'serial fine actualH time interpolation')
        call compute_sdc_mpi_amr_tree_hydro_3d(tree,owners,species,local_gas(1),local_temp,local_source, &
          spacing,dt,faces,distributed(1)%gas,distributed(1)%flux,MPI_COMM_WORLD,good,target_depth=1, &
          coarse_start=local_packet,coarse_end=local_future,coarse_fractions=[.5_dp*real(pass-1,dp)], &
          source_start=local_source,source_end=local_source,source_fractions=[.5_dp*real(pass-1,dp)], &
          internal_energy=local_held,internal_energy_rate=distributed(1)%held, &
          internal_energy_integrated_flux=distributed(1)%aux,internal_energy_pressure_work=distributed(1)%work)
        call require(good,'MPI fine actualH time interpolation')
        if(rank==owners(2)) call compare_node(serial,distributed(1),2,'serial MPI parent H chronology')
        if(rank==owners(1)) call require(all(distributed(1)%held%fields(1)%values==0).and. &
          all(distributed(1)%work%fields(1)%values==0),'inactive parent has zero residual and pressure work')
        if(pass==1) first_rate=serial%held
      enddo
      call require(any(serial%held%fields(2)%values/=first_rate%fields(2)%values),'future parent H affects fine residual')
    endif
    call check_native()
  enddo
  if(rank==0) print *, 'PASS: regular held Godunov serial/MPI spatial families and rollback'
  call MPI_Finalize(status)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if(condition) return
    print *, 'FAIL rank',rank,': ',message
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
  subroutine initialize()
    real(dp) :: x(3),thermal,v(3),mw,rg,cp,cv,gamma,enthalpy,e,entropy,bulk_e
    do a=1,2
      if(allocated(gas(a)%fields)) deallocate(gas(a)%fields)
      allocate(gas(a)%fields(size(tree%nodes)))
    enddo
    if(allocated(temp%fields)) deallocate(temp%fields,held%fields,source%fields)
    allocate(temp%fields(size(tree%nodes)),held%fields(size(tree%nodes)),source%fields(size(tree%nodes)))
    call mixture_mass_properties(species,y,1000._dp,mw,rg,cp,cv,gamma,enthalpy,bulk_e,entropy,good)
    call require(good,'bulk thermo')
    do node=1,size(tree%nodes)
      n=tree%nodes(node)%shape;hh=spacing/real(2**tree%nodes(node)%level,dp)
      allocate(gas(1)%fields(node)%values(nv,n(1),n(2),n(3)),gas(2)%fields(node)%values(nv,n(1),n(2),n(3)), &
        temp%fields(node)%values(1,n(1),n(2),n(3)),held%fields(node)%values(1,n(1),n(2),n(3)), &
        source%fields(node)%values(nv,n(1),n(2),n(3)))
      source%fields(node)%values=0
      do k=1,n(3)
      do j=1,n(2)
      do i=1,n(1)
        x=(real([i,j,k],dp)-.5_dp)*hh
        thermal=700+25*sin(2*acos(-1._dp)*x(1)/.004_dp)
        v=[8*sin(2*acos(-1._dp)*x(1)/.004_dp),2*sin(2*acos(-1._dp)*x(2)/.004_dp),0._dp]
        call mixture_mass_properties(species,y,thermal,mw,rg,cp,cv,gamma,enthalpy,e,entropy,good)
        call require(good,'actual thermo')
        do a=1,2
          gas(a)%fields(node)%values(:,i,j,k)=0
          gas(a)%fields(node)%values(1,i,j,k)=.4_dp;gas(a)%fields(node)%values(2:4,i,j,k)=.4_dp*v
          gas(a)%fields(node)%values(5,i,j,k)=.4_dp*merge(bulk_e,e,a==1)+.2_dp*sum(v*v)
          gas(a)%fields(node)%values(6:,i,j,k)=.4_dp*y
        enddo
        temp%fields(node)%values(1,i,j,k)=thermal;held%fields(node)%values(1,i,j,k)=.4_dp*e
      enddo
      enddo
      enddo
    enddo
  end subroutine
  subroutine localize(full,local)
    type(amr_tree_state_3d), intent(in) :: full
    type(amr_tree_state_3d), intent(inout) :: local
    integer :: peer
    if(allocated(local%fields)) deallocate(local%fields)
    allocate(local%fields(size(owners)))
    do peer=1,size(owners)
      if(owners(peer)==rank) local%fields(peer)%values=full%fields(peer)%values
    enddo
  end subroutine
  subroutine pack(full)
    type(amr_tree_state_3d), intent(inout) :: full
    integer :: peer,dims(3)
    if(allocated(full%fields)) deallocate(full%fields)
    allocate(full%fields(size(owners)))
    do peer=1,size(owners)
      dims=tree%nodes(peer)%shape
      allocate(full%fields(peer)%values(nv+2,dims(1),dims(2),dims(3)))
      full%fields(peer)%values(:nv,:,:,:)=gas(1)%fields(peer)%values
      full%fields(peer)%values(nv+1,:,:,:)=temp%fields(peer)%values(1,:,:,:)
      full%fields(peer)%values(nv+2,:,:,:)=held%fields(peer)%values(1,:,:,:)
    enddo
  end subroutine
  subroutine check_native()
    integer :: cycle,peer
    real(dp) :: tolerance,inventory_change,reset_residual,work_oracle,work_scale,old_inventory(1),new_inventory(1),hboundary
    type(amr_tree_state_3d) :: full_held
    faces='periodic'
    ! The one-root interval has an independent spatial oracle. Native
    ! Godunov must retain the full H family and subtract pdivu exactly once.
    if(kind/=2) then
      call compute_sdc_amr_patch_tree_hydro_3d(tree,species,gas(1),temp,source,spacing,dt,faces, &
        serial%gas,serial%flux,good,difmag=0._dp,internal_energy=held, &
        internal_energy_rate=serial%held,internal_energy_integrated_flux=serial%aux, &
        internal_energy_pressure_work=serial%work)
      call require(good,'native spatial oracle')
    endif
    do cycle=1,3
      if(cycle==3.and.kind/=1) cycle
      policy%allow_small_energy=.false.;policy%allow_negative_energy=.true.
      policy%update_total_from_internal=.false.
      if(cycle==3) then
        ! This pinned branch repairs bulk energy from held H. Its nonzero
        ! initial/final corrections must enter the accepted ledger once.
        policy%allow_small_energy=.true.;policy%allow_negative_energy=.false.
        policy%update_total_from_internal=.true.
      endif
      native_serial%gas=gas(1);native_serial%temp=temp;native_serial%held=held
      native_mpi%gas=local_gas(1);native_mpi%temp=local_temp;native_mpi%held=local_held
      if(allocated(native_serial%boundary)) deallocate(native_serial%boundary,native_serial%source, &
        native_mpi%boundary,native_mpi%source)
      allocate(native_serial%boundary(nv),native_serial%source(nv),native_mpi%boundary(nv),native_mpi%source(nv))
      native_serial%prefix=2;native_mpi%prefix=2
      native_serial%work_change=19;native_mpi%work_change=19
      native_serial%boundary=3;native_serial%source=4;native_serial%correction=5;native_serial%theta=6
      native_mpi%boundary=3;native_mpi%source=4;native_mpi%correction=5;native_mpi%theta=6
      call native_step(native_serial,.true.,cycle,good);call require(good,'serial native held Godunov accepted')
      call native_step(native_mpi,.false.,cycle,good);call require(good,'MPI native held Godunov accepted')
      do peer=1,size(owners)
        if(owners(peer)/=rank) cycle
        tolerance=128*epsilon(1._dp)*max(1._dp,maxval(abs(native_serial%held%fields(peer)%values)))
        call require(maxval(abs(native_serial%held%fields(peer)%values-native_mpi%held%fields(peer)%values))<=tolerance, &
          'serial MPI native held Godunov endpoint')
        call require(maxval(abs(native_serial%gas%fields(peer)%values-native_mpi%gas%fields(peer)%values))<= &
          128*epsilon(1._dp)*max(1._dp,maxval(abs(native_serial%gas%fields(peer)%values))), &
          'serial MPI native Godunov gas')
        call require(maxval(abs(native_serial%temp%fields(peer)%values-native_mpi%temp%fields(peer)%values))<= &
          128*epsilon(1._dp)*maxval(abs(native_serial%temp%fields(peer)%values)),'serial MPI native Godunov EOS')
        call compare_faces(native_mpi%aux(peer),native_mpi%physical_aux(peer),'Godunov H register retains full family')
        call compare_faces(native_mpi%flux(peer),native_mpi%physical(peer),'Godunov gas register retains full family')
        if(kind/=2) then
          if(maxval(abs(native_mpi%held%fields(peer)%values-(held%fields(peer)%values+ &
            dt*serial%held%fields(peer)%values)))>tolerance) print *, 'Native H oracle kind/cycle/error/tol/prefix:', &
              kind,cycle,maxval(abs(native_mpi%held%fields(peer)%values-(held%fields(peer)%values+ &
              dt*serial%held%fields(peer)%values))),tolerance,native_mpi%prefix
          call require(maxval(abs(native_mpi%held%fields(peer)%values-(held%fields(peer)%values+ &
            dt*serial%held%fields(peer)%values)))<=tolerance,'Godunov UEINT uses full spatial residual')
          call compare_faces(serial%aux(peer),native_mpi%aux(peer),'full Godunov H spatial face oracle')
          if(cycle==3) then
            inventory_change=product(spacing)*(sum(native_mpi%gas%fields(peer)%values(5,:,:,:))- &
              sum(gas(1)%fields(peer)%values(5,:,:,:)))
            reset_residual=(native_mpi%prefix-2)-(inventory_change-native_mpi%boundary(5)-native_mpi%source(5))
            call require(abs(native_mpi%prefix-2)>1e-4_dp,'nonzero Godunov reset ledger exercised')
            call require(abs(reset_residual)<=256*epsilon(1._dp)*max(1._dp,abs(inventory_change)), &
              'accepted Godunov reset ledger excludes discarded predictor resets')
          endif
        endif
      enddo
      call require(abs(native_mpi%prefix-native_serial%prefix)<=128*epsilon(1._dp)* &
        max(1._dp,abs(native_serial%prefix)),'serial MPI Godunov reset history')
      call require(native_mpi%work_change==native_serial%work_change,'serial MPI accepted pressure-work history')
      if(kind/=2) then
        work_oracle=-dt*product(spacing)*sum(serial%work%fields(1)%values)
        work_scale=dt*product(spacing)*sum(abs(serial%work%fields(1)%values))
        call require(abs(native_mpi%work_change-work_oracle)<=256*epsilon(1._dp)*max(1e-30_dp,work_scale), &
          'pressure-work output uses final Godunov volume source once')
      endif
      ! Independent composite H inventory uses the accepted scalar boundary
      ! family, including fine boundary replacements after synchronization.
      full_held=native_serial%held
      call composite_integrals_amr_patch_tree_3d(tree,held,spacing,old_inventory,good)
      call require(good,'old composite H inventory')
      call composite_integrals_amr_patch_tree_3d(tree,full_held,spacing,new_inventory,good)
      call require(good,'new composite H inventory')
      n=tree%nodes(1)%shape
      hboundary=(sum(native_serial%physical_aux(1)%x(1,1,:,:))- &
        sum(native_serial%physical_aux(1)%x(1,n(1)+1,:,:)))*spacing(2)*spacing(3)+ &
        (sum(native_serial%physical_aux(1)%y(1,:,1,:))-sum(native_serial%physical_aux(1)%y(1,:,n(2)+1,:)))* &
        spacing(1)*spacing(3)+(sum(native_serial%physical_aux(1)%z(1,:,:,1))- &
        sum(native_serial%physical_aux(1)%z(1,:,:,n(3)+1)))*spacing(1)*spacing(2)
      call require(abs(new_inventory(1)-old_inventory(1)-hboundary-native_serial%work_change)<= &
        512*epsilon(1._dp)*max(1e-30_dp,abs(old_inventory(1)),abs(new_inventory(1))), &
        'composite H inventory closes with boundary and accepted pressure work')
    enddo
    if(kind==1) then
      native_saved=native_mpi
      if(rank==2) then
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,native_mpi%gas,native_mpi%temp, &
          spacing,dt,faces,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          native_mpi%boundary,native_mpi%source,native_mpi%correction,native_mpi%theta,MPI_COMM_WORLD,good, &
          internal_energy=native_mpi%held,energy_correction=native_mpi%prefix,energy_policy=policy)
      else
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,native_mpi%gas,native_mpi%temp, &
          spacing,dt,faces,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp, &
          native_mpi%boundary,native_mpi%source,native_mpi%correction,native_mpi%theta,MPI_COMM_WORLD,good, &
          internal_energy=native_mpi%held,energy_correction=native_mpi%prefix,energy_policy=policy, &
          internal_energy_pressure_work_change=native_mpi%work_change)
      endif
      call require(.not.good,'native pressure-work presence mismatch refuses collectively')
      call require(native_mpi%prefix==native_saved%prefix.and.native_mpi%work_change==native_saved%work_change.and. &
        all(native_mpi%boundary==native_saved%boundary).and.all(native_mpi%source==native_saved%source), &
        'native pressure-work disagreement preserves all ledgers')
      if(rank==owners(1)) call require(all(native_mpi%gas%fields(1)%values==native_saved%gas%fields(1)%values).and. &
        all(native_mpi%held%fields(1)%values==native_saved%held%fields(1)%values), 'native pressure-work gas/H rollback')
    endif
  end subroutine
  subroutine native_step(t,is_serial,cycle,accepted)
    type(native_state), intent(inout) :: t
    logical, intent(in) :: is_serial
    integer, intent(in) :: cycle
    logical, intent(out) :: accepted
    if(is_serial) then
      if(cycle==1) then
        call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,t%gas,t%temp,spacing,dt,faces,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,t%boundary,t%source, &
          t%correction,t%theta,accepted,difmag=0._dp,mol_mode=.false.,internal_energy=t%held,energy_correction=t%prefix, &
          energy_policy=policy,integrated_flux=t%flux,internal_energy_integrated_flux=t%aux, &
          physical_integrated_flux=t%physical,physical_internal_energy_integrated_flux=t%physical_aux, &
          internal_energy_pressure_work_change=t%work_change)
      else
        call advance_amr_patch_tree_sdc_subcycle_3d(tree,species,reactions,transport,t%gas,t%temp,spacing,dt,faces,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,t%boundary,t%source, &
          t%correction,t%theta,accepted,difmag=0._dp,internal_energy=t%held,energy_correction=t%prefix, &
          energy_policy=policy,integrated_flux=t%flux,internal_energy_integrated_flux=t%aux, &
          physical_integrated_flux=t%physical,physical_internal_energy_integrated_flux=t%physical_aux, &
          internal_energy_pressure_work_change=t%work_change)
      endif
    else
      if(cycle==1) then
        call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,t%gas,t%temp,spacing,dt,faces,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,t%boundary,t%source, &
          t%correction,t%theta,MPI_COMM_WORLD,accepted,difmag=0._dp,mol_mode=.false., &
          internal_energy=t%held,energy_correction=t%prefix,energy_policy=policy, &
          integrated_flux=t%flux,internal_energy_integrated_flux=t%aux, &
          physical_integrated_flux=t%physical,physical_internal_energy_integrated_flux=t%physical_aux, &
          internal_energy_pressure_work_change=t%work_change)
      else
        call advance_mpi_amr_tree_sdc_subcycle_3d(tree,owners,species,reactions,transport,t%gas,t%temp,spacing,dt,faces,2,.false., &
          [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-6_dp,1e-14_dp,t%boundary,t%source, &
          t%correction,t%theta,MPI_COMM_WORLD,accepted,difmag=0._dp, &
          internal_energy=t%held,energy_correction=t%prefix,energy_policy=policy, &
          integrated_flux=t%flux,internal_energy_integrated_flux=t%aux, &
          physical_integrated_flux=t%physical,physical_internal_energy_integrated_flux=t%physical_aux, &
          internal_energy_pressure_work_change=t%work_change)
      endif
    endif
  end subroutine
  subroutine compare_faces(one,two,message)
    type(amr_tree_flux_3d), intent(in) :: one,two
    character(len=*), intent(in) :: message
    call require(all(one%x==two%x).and.all(one%y==two%y).and.all(one%z==two%z),message)
  end subroutine
  subroutine compare_node(one,two,peer,message)
    type(result), intent(in) :: one,two
    integer, intent(in) :: peer
    character(len=*), intent(in) :: message
    call require(all(one%gas%fields(peer)%values==two%gas%fields(peer)%values),message)
    call require(all(one%held%fields(peer)%values==two%held%fields(peer)%values),message)
    call require(all(one%work%fields(peer)%values==two%work%fields(peer)%values),message)
    call compare_faces(one%flux(peer),two%flux(peer),message)
    call compare_faces(one%aux(peer),two%aux(peer),message)
  end subroutine
end program
