program test_mpi_mol_reaction_3d
  use iso_fortran_env, only: int64
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use mixture_thermo_mod, only: mixture_mass_properties
  use constant_volume_reactor_mod, only: advance_forced_reactor_interval
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d,build_eb_tree_extended_geometry_3d, &
    average_down_reactive_eb_tree_3d,synchronize_reactive_eb_tree_flux_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(0)
  type(amr_patch_tree_3d) :: tree
  type(eb_geometry_3d) :: geometry(1),extended(1)
  type(amr_tree_state_3d) :: initial,initial_t,state,temp,saved,saved_t
  type(eb_reactive_tree_3d) :: serial,serial_initial
  real(dp), allocatable :: y(:),q(:),u(:),reference(:),b(:),w(:),ch(:),source(:),aux(:)
  real(dp) :: mw,r,cp,cv,gamma,h,e,entropy,t,c,reference_t,theta,error
  integer :: rank,status,nv,s,iteration,steps,i,j,k,calls,global_calls,bad_species,aux_calls,serial_aux_calls
  integer :: pulse_species
  real(dp) :: pulse_energy,pulse_peak,global_peak
  real(dp) :: aux_force(3)=[2e8_dp,-3e8_dp,4e8_dp]
  real(dp) :: force(3)=[1e9_dp,-2e9_dp,3e9_dp],weight,volume
  real(dp), allocatable :: serial_b(:),serial_w(:),serial_ch(:),serial_source(:)
  logical :: ok,fail_late
  character(len=16) :: faces(6)
  call MPI_Init(status);call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  call load_h2o2_elementary_mechanism(reactions,ok);call require(ok,'reactions')
  nv=reactive_nvar(size(species))
  allocate(y(size(species)),q(reactive_nprim(size(species))),u(nv),reference(size(species)), &
    b(nv),w(nv),ch(nv),source(nv),aux(nv))
  y=0
  do s=1,size(species)
    select case(species(s)%name)
    case('H2');y(s)=.1_dp
    case('O2');y(s)=.2_dp
    case('N2');y(s)=.699_dp
    case('H');y(s)=.001_dp
    end select
  end do
  call mixture_mass_properties(species,y,1200._dp,mw,r,cp,cv,gamma,h,e,entropy,ok)
  call require(ok,'properties')
  q=0;q(1)=1;q(5)=r*1200;q(6:)=y
  call reactive_primitive_to_conserved(species,q,u,t,c,ok);call require(ok,'state')
  reference=u(6:);reference_t=t
  call advance_forced_reactor_interval(species,reactions,u(5),0*y,0._dp,1e-10_dp,1e-8_dp,1e-12_dp, &
    reference,reference_t,steps,ok)
  call require(ok.and.maxval(abs(reference-u(6:)))>1e-10_dp,'nontrivial homogeneous reactor')
  call initialize_amr_patch_tree_3d([2,2,2],[integer::],boxes,tree,ok);call require(ok,'tree')
  call build_axis_plane_eb_geometry_3d(14,14,14,-3._dp,4._dp,-3._dp,4._dp,-3._dp,4._dp, &
    'x',-10._dp,extended(1),ok);call require(ok,'extended regular geometry')
  call extract_eb_geometry_3d(extended(1),[7,7,7],[8,8,8],geometry(1),ok);call require(ok,'geometry')
  allocate(initial%fields(1),initial_t%fields(1))
  if(rank==0) then
    allocate(initial%fields(1)%values(nv,2,2,2),initial_t%fields(1)%values(1,2,2,2))
    do k=1,2
    do j=1,2
    do i=1,2
      initial%fields(1)%values(:,i,j,k)=u
    end do
    end do
    end do
    initial_t%fields(1)%values=t
  end if
  faces='periodic';fail_late=.false.
  do iteration=1,2
    state=initial;temp=initial_t;calls=0
    call step(iteration,ok);call require(ok,'reacting MOL step')
    call MPI_Allreduce(calls,global_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(global_calls==8*iteration,'MOL reaction iteration count')
    error=0
    if(rank==0) then
      do k=1,2
      do j=1,2
      do i=1,2
        error=max(error,maxval(abs(state%fields(1)%values(6:,i,j,k)-reference)))
        error=max(error,abs(temp%fields(1)%values(1,i,j,k)-reference_t)/reference_t)
      end do
      end do
      end do
      call require(error<1e-11_dp,'homogeneous reactor agreement')
      call require(maxval(abs(state%fields(1)%values(:,1,1,1)-u-b-w-ch-source))<1e-9_dp,'chemical ledger')
      print *, 'MOL iterations, reactor calls, maximum error:',iteration,global_calls,error
    end if
  end do
  ! Analytic homogeneous forcing distinguishes the first old/new pass from
  ! later reacting MOL passes, which retain old external forcing only in the
  ! pinned source. Compare the independently owned serial and MPI APIs too.
  allocate(serial_b(nv),serial_w(nv),serial_ch(nv),serial_source(nv))
  volume=8*geometry(1)%dx*geometry(1)%dy*geometry(1)%dz
  if(rank==0) then
    call initialize_reactive_eb_tree_3d(species,tree,geometry(1),initial%fields(1)%values, &
      initial_t%fields(1)%values(1,:,:,:),serial_initial,ok)
    call require(ok,'serial homogeneous forced tree')
  endif
  do iteration=1,3
    state=initial;temp=initial_t;calls=0;aux_calls=0
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
      1e-10_dp,iteration,.true.,[.false.,.false.,.false.,.false.],force, &
      1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok,evaluate_reaction=reactor,mol_mode=.true., &
      auxiliary_source=positive_source,auxiliary_change=aux)
    call require(ok.and.aux_calls==2,'homogeneous forced reacting MOL')
    call MPI_Allreduce(calls,global_calls,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,status)
    call require(global_calls==8*iteration,'forced MOL reaction iteration count')
    weight=1
    if(iteration>1) weight=.5_dp
    call require(maxval(abs(source(2:4)-volume*weight*1e-10_dp*force))<1e-11_dp, &
      'MOL external momentum ledger follows pinned old-only later passes')
    call require(maxval(abs(aux(2:4)-volume*weight*1e-10_dp*aux_force))<1e-11_dp, &
      'MOL auxiliary momentum ledger follows old-only later passes')
    if(rank==0) then
      do k=1,2
      do j=1,2
      do i=1,2
        call require(maxval(abs(state%fields(1)%values(2:4,i,j,k)-u(2:4)- &
          weight*1e-10_dp*(force+aux_force)))<1e-11_dp, &
          'MOL homogeneous momentum analytic update')
      enddo
      enddo
      enddo
      serial=serial_initial;calls=0;serial_aux_calls=0
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-10_dp,iteration,.true., &
        [.false.,.false.,.false.,.false.],force,1e-8_dp,1e-12_dp,serial,extended,faces, &
        serial_b,serial_w,serial_ch,serial_source,theta,ok,evaluate_reaction=reactor,mol_mode=.true., &
        auxiliary_source=serial_positive_source,auxiliary_change=aux)
      call require(ok.and.calls==8*iteration.and.serial_aux_calls==2,'serial forced reacting MOL')
      call require(maxval(abs(aux(2:4)-volume*weight*1e-10_dp*aux_force))<1e-11_dp, &
        'serial MOL old-only auxiliary momentum ledger')
      call require(maxval(abs(serial_source(2:4)-volume*weight*1e-10_dp*force))<1e-11_dp, &
        'serial MOL analytic momentum ledger')
      error=maxval(abs(state%fields(1)%values-serial%nodes(1)%state)/max(1._dp,abs(serial%nodes(1)%state)))
      error=max(error,maxval(abs(temp%fields(1)%values(1,:,:,:)-serial%nodes(1)%temperature)/ &
        serial%nodes(1)%temperature))
      call require(error<1e-11_dp,'serial and MPI forced MOL state and temperature agreement')
      print *, 'Forced MOL iterations and serial/MPI maximum scaled error:',iteration,error
    endif
  enddo
  ! A new-time-only pulse belongs to the first reacting predictor, then
  ! disappears from the pinned later MOL RHS. The corrected spatial flux
  ! must not export that unfunded transient species. Exercise serial/MPI
  ! limiting, rather than just the homogeneous safe-pass path above.
  pulse_species=0
  do s=1,size(species)
    if(species(s)%name=='H2O') pulse_species=s
  enddo
  call require(pulse_species>0,'pulse species')
  call require(u(5+pulse_species)==0,'pulse starts from absent species')
  block
    real(dp) :: pulse_y(size(species))
    pulse_y=0;pulse_y(pulse_species)=1
    call mixture_mass_properties(species,pulse_y,1200._dp,mw,r,cp,cv,gamma,h,pulse_energy,entropy,ok)
    call require(ok,'pulse thermodynamics')
  end block
  state=initial;temp=initial_t;aux_calls=0;pulse_peak=0;calls=0
  call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
    1e-8_dp,2,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
    1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok,evaluate_reaction=pulse_reactor,mol_mode=.true., &
    auxiliary_source=pulse_source,auxiliary_change=aux)
  call require(ok.and.aux_calls==2,'MOL transient pulse correction')
  call MPI_Allreduce(pulse_peak,global_peak,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
  call require(global_peak>0,'first MOL pulse really entered the reacted predictor')
  call require(all(aux==0),'later MOL pulse has no funded auxiliary inventory')
  if(rank==0) then
    call require(all(state%fields(1)%values(5+pulse_species,:,:,:)==0),'unfunded pulse retained zero species inventory')
    serial=serial_initial;serial_aux_calls=0;pulse_peak=0
    call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-8_dp,2,.true., &
      [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial,extended,faces, &
      serial_b,serial_w,serial_ch,serial_source,theta,ok,evaluate_reaction=pulse_reactor,mol_mode=.true., &
      auxiliary_source=serial_pulse_source,auxiliary_change=aux)
    call require(ok.and.serial_aux_calls==2.and.pulse_peak>0,'serial MOL transient pulse correction')
    call require(all(serial%nodes(1)%state(5+pulse_species,:,:,:)==0),'serial pulse keeps zero species inventory')
    error=maxval(abs(state%fields(1)%values-serial%nodes(1)%state)/max(1._dp,abs(serial%nodes(1)%state)))
    call require(error<1e-11_dp,'active MOL flux correction serial/MPI state agreement')
    print *, 'MOL transient pulse peak, final tracer zero, serial/MPI error:',global_peak,error
  endif
  call check_disabled_chemistry_rates()
  call check_nested_flux_registration()
  call check_partial_registered_reflux()
  ! Even a sub-tolerance negative species must fail before publication: the
  ! restart format requires genuinely nonnegative partial densities.
  state=initial;temp=initial_t;saved=state;saved_t=temp
  bad_species=5+findloc(y,0._dp,dim=1);aux_calls=0;aux=31
  b=17;w=19;ch=23;source=29
  call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
    1e-10_dp,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
    1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok, &
    auxiliary_source=negative_source,auxiliary_change=aux,mol_mode=.true.)
  call require(aux_calls==2,'negative source reached both stages')
  call require(.not.ok,'sub-tolerance negative species rejection')
  call require(all(aux==31),'negative-species auxiliary rollback')
  call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(source==29),'negative-species ledger rollback')
  if(rank==0) then
    call require(all(state%fields(1)%values==saved%fields(1)%values),'negative-species state rollback')
    call require(all(temp%fields(1)%values==saved_t%fields(1)%values),'negative-species temperature rollback')
  end if
  saved=state;saved_t=temp;calls=0;fail_late=.true.;b=17;w=19;ch=23;source=29
  call step(2,ok);call require(.not.ok,'late reactor failure')
  call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(source==29),'ledger rollback')
  if(rank==0) then
    call require(all(state%fields(1)%values==saved%fields(1)%values),'state rollback')
    call require(all(temp%fields(1)%values==saved_t%fields(1)%values),'temperature rollback')
    print *, 'PASS: reacting MOL and late-stage transactional rollback'
  end if
  call MPI_Finalize(status)
contains
  subroutine check_nested_flux_registration()
    type(amr_patch_3d) :: patch(1)
    type(amr_patch_tree_3d) :: nested
    type(eb_geometry_3d) :: finest,padded,physical_geometry(2)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(eb_reactive_tree_3d) :: initialized,candidate
    type(amr_tree_state_3d) :: local,local_t
    type(eb_integrated_flux_3d), allocatable :: registered(:),physical(:),serial_registered(:),serial_physical(:)
    real(dp) :: root_u(nv,2,2,2),root_t(2,2,2),register_weight,local_error,global_error,denom
    integer(int64) :: counts(2),serial_counts(2)
    integer :: owners(2),node,i,j,k,pass
    owners=[0,1]
    call initialize_amr_patch_3d(2,2,2,1,2,1,2,1,2,2,patch(1),ok)
    call require(ok,'registration full refinement patch')
    call initialize_amr_patch_tree_3d([2,2,2],[1],patch,nested,ok)
    call require(ok,'registration two-level topology')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-10._dp,finest,ok)
    call require(ok,'registration nested regular geometry')
    do k=1,2
    do j=1,2
    do i=1,2
      root_u(:,i,j,k)=u
    enddo
    enddo
    enddo
    root_t=t
    call initialize_reactive_eb_tree_3d(species,nested,finest,root_u,root_t,initialized,ok)
    call require(ok,'registration nested initial state')
    do node=1,2
      physical_geometry(node)=initialized%nodes(node)%geometry
    enddo
    call build_axis_plane_eb_geometry_3d(28,28,28,-3._dp,4._dp,-3._dp,4._dp,-3._dp,4._dp,'x',-10._dp,padded,ok)
    call require(ok,'registration padded finest geometry')
    call build_eb_tree_extended_geometry_3d(nested,physical_geometry,padded,wide,ok,halo_width=6)
    call require(ok,'registration nested halo geometry')
    allocate(local%fields(2),local_t%fields(2))
    denom=1e-10_dp*q(5)
    do pass=1,3,2
      candidate=initialized
      do node=1,2
        if(owners(node)/=rank) cycle
        local%fields(node)%values=initialized%nodes(node)%state
        allocate(local_t%fields(node)%values(1,nested%nodes(node)%shape(1), &
          nested%nodes(node)%shape(2),nested%nodes(node)%shape(3)))
        local_t%fields(node)%values(1,:,:,:)=initialized%nodes(node)%temperature
      enddo
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical_geometry, &
        local,local_t,1e-10_dp,pass,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-8_dp,1e-12_dp,wide,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
        registered_flux=registered,integrated_flux=physical,level_steps=counts)
      call require(ok.and.all(counts==[1_int64,2_int64]),'nested registration closes both fine intervals')
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-10_dp,pass,.false., &
        [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,candidate,wide,faces, &
        serial_b,serial_w,serial_ch,serial_source,theta,ok,mol_mode=.true., &
        registered_flux=serial_registered,integrated_flux=serial_physical,level_steps=serial_counts)
      call require(ok.and.all(serial_counts==counts),'serial nested registration intervals')
      register_weight=merge(1._dp,.5_dp,pass==1);local_error=0
      do node=1,2
        if(owners(node)/=rank) cycle
        local_error=max(local_error,maxval(abs(registered(node)%x(2,:,:,:)-register_weight*denom))/denom, &
          maxval(abs(registered(node)%y(3,:,:,:)-register_weight*denom))/denom, &
          maxval(abs(registered(node)%z(4,:,:,:)-register_weight*denom))/denom, &
          maxval(abs(physical(node)%x(2,:,:,:)-denom))/denom, &
          maxval(abs(physical(node)%y(3,:,:,:)-denom))/denom,maxval(abs(physical(node)%z(4,:,:,:)-denom))/denom)
        call require(all(registered(node)%x==serial_registered(node)%x).and. &
          all(registered(node)%y==serial_registered(node)%y).and.all(registered(node)%z==serial_registered(node)%z), &
          'nested serial/MPI registered flux exact agreement')
        call require(all(physical(node)%x==serial_physical(node)%x).and. &
          all(physical(node)%y==serial_physical(node)%y).and.all(physical(node)%z==serial_physical(node)%z), &
          'nested serial/MPI physical flux exact agreement')
        deallocate(local_t%fields(node)%values)
      enddo
      call MPI_Allreduce(local_error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
      call require(global_error<1e-12_dp,'nested registered and physical analytic pressure flux')
      if(rank==0) print *, 'Nested MOL registration iterations, pressure flux error:',pass,global_error
    enddo
    do node=1,2
      if(owners(node)/=rank) cycle
      registered(node)%x=13;registered(node)%y=13;registered(node)%z=13
      allocate(local_t%fields(node)%values(1,nested%nodes(node)%shape(1), &
        nested%nodes(node)%shape(2),nested%nodes(node)%shape(3)))
      local_t%fields(node)%values(1,:,:,:)=initialized%nodes(node)%temperature
    enddo
    if(rank==0) then
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical_geometry, &
        local,local_t,1e-10_dp,3,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-8_dp,1e-12_dp,wide,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok,mol_mode=.true., &
        registered_flux=registered)
    else
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,owners,physical_geometry, &
        local,local_t,1e-10_dp,3,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-8_dp,1e-12_dp,wide,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok,mol_mode=.true.)
    endif
    call require(.not.ok,'rank-disagreeing registration presence rejects collectively')
    do node=1,2
      if(owners(node)/=rank) cycle
      call require(all(registered(node)%x==13).and.all(registered(node)%y==13).and.all(registered(node)%z==13), &
        'registration presence rejection retains output')
    enddo
  end subroutine check_nested_flux_registration
  ! Recover the unsynchronized FV endpoint from the measured physical face
  ! flux, then apply the registered flux as an independent reflux oracle.
  ! A partially refined pressure gradient makes the two flux choices observable.
  subroutine check_partial_registered_reflux()
    type(amr_patch_3d) :: patch(1)
    type(amr_patch_tree_3d) :: nested
    type(eb_geometry_3d) :: finest,padded,g(2)
    type(eb_geometry_3d), allocatable :: wide(:)
    type(eb_reactive_tree_3d) :: initialized,actual,expected,physical_control
    type(amr_tree_state_3d) :: local,local_t
    type(eb_integrated_flux_3d), allocatable :: registered(:),physical(:),oracle_flux(:)
    real(dp) :: root_u(nv,2,2,2),root_t(2,2,2),cell(nv),pt,cs,dt,error,control_gap,mpi_error,global_error
    real(dp) :: oracle_q(size(q))
    integer :: node,x,y,z,d(3),mode,pass
    call initialize_amr_patch_3d(2,2,2,1,1,1,2,1,2,2,patch(1),ok)
    call require(ok,'partial reflux patch')
    call initialize_amr_patch_tree_3d([2,2,2],[1],patch,nested,ok)
    call require(ok,'partial reflux topology')
    call build_axis_plane_eb_geometry_3d(4,4,4,0._dp,1._dp,0._dp,1._dp,0._dp,1._dp,'x',-10._dp,finest,ok)
    call require(ok,'partial reflux geometry')
    do z=1,2
    do y=1,2
    do x=1,2
      root_u(:,x,y,z)=u
    enddo
    enddo
    enddo
    root_t=t
    call initialize_reactive_eb_tree_3d(species,nested,finest,root_u,root_t,initialized,ok)
    call require(ok,'partial reflux initialization')
    d=nested%nodes(2)%shape
    do z=1,d(3)
    do y=1,d(2)
    do x=1,d(1)
      oracle_q=q;oracle_q(5)=q(5)*(1+.05_dp*real(x,dp))
      call reactive_primitive_to_conserved(species,oracle_q,cell,pt,cs,ok)
      call require(ok,'partial pressure perturbation')
      initialized%nodes(2)%state(:,x,y,z)=cell;initialized%nodes(2)%temperature(x,y,z)=pt
    enddo
    enddo
    enddo
    call average_down_reactive_eb_tree_3d(species,initialized,ok)
    call require(ok,'partial restricted initial state')
    do node=1,2
      g(node)=initialized%nodes(node)%geometry
    enddo
    call build_axis_plane_eb_geometry_3d(28,28,28,-3._dp,4._dp,-3._dp,4._dp,-3._dp,4._dp,'x',-10._dp,padded,ok)
    call require(ok,'partial padded geometry')
    call build_eb_tree_extended_geometry_3d(nested,g,padded,wide,ok,halo_width=6)
    call require(ok,'partial prepared halo geometry')
    allocate(local%fields(2),local_t%fields(2));dt=1e-8_dp
    do mode=1,2
    do pass=1,3,2
      actual=initialized;expected=initialized
      if(mode==1) then
        call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,dt,pass,.false., &
          [.false.,.false.,.false.,.false.],0*force,1e-8_dp,1e-12_dp,actual,wide,faces, &
          serial_b,serial_w,serial_ch,serial_source,theta,ok,mol_mode=.true., &
          registered_flux=registered,integrated_flux=physical)
      else
        call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,dt,pass,.false., &
          [.false.,.false.,.false.,.false.],0*force,1e-8_dp,1e-12_dp,actual,wide,faces, &
          serial_b,serial_w,serial_ch,serial_source,theta,ok,mol_mode=.true., &
          registered_flux=registered,integrated_flux=physical)
      endif
      call require(ok,'partial serial native registered reflux')
      do node=1,2
        d=nested%nodes(node)%shape
        do z=1,d(3)
        do y=1,d(2)
        do x=1,d(1)
          expected%nodes(node)%state(:,x,y,z)=initialized%nodes(node)%state(:,x,y,z)- &
            (physical(node)%x(:,x,y,z)-physical(node)%x(:,x-1,y,z))/g(node)%dx- &
            (physical(node)%y(:,x,y,z)-physical(node)%y(:,x,y-1,z))/g(node)%dy- &
            (physical(node)%z(:,x,y,z)-physical(node)%z(:,x,y,z-1))/g(node)%dz
        enddo
        enddo
        enddo
        call recover_reactive_temperatures_3d(species,expected%nodes(node)%state,initialized%nodes(node)%temperature, &
          d(1),d(2),d(3),expected%nodes(node)%temperature,ok)
        call require(ok,'partial FV endpoint EOS')
      enddo
      physical_control=expected
      if(pass>1) then
        oracle_flux=physical
        call synchronize_reactive_eb_tree_flux_3d(species,physical_control,oracle_flux,ok,boundary_faces=faces)
        call require(ok,'partial physical flux control')
        oracle_flux=registered
        call synchronize_reactive_eb_tree_flux_3d(species,expected,oracle_flux,ok,boundary_faces=faces)
        call require(ok,'partial native register oracle')
      else
        ! Full-weight physical output already contains canonical reflux faces.
        call average_down_reactive_eb_tree_3d(species,expected,ok)
        call require(ok,'partial physical endpoint restriction')
      endif
      error=0;control_gap=0
      do node=1,2
        error=max(error,maxval(abs(actual%nodes(node)%state-expected%nodes(node)%state)/ &
          max(1._dp,abs(expected%nodes(node)%state))))
        control_gap=max(control_gap,maxval(abs(actual%nodes(node)%state-physical_control%nodes(node)%state)/ &
          max(1._dp,abs(actual%nodes(node)%state))))
        if(rank/=node-1) cycle
        local%fields(node)%values=initialized%nodes(node)%state
        d=nested%nodes(node)%shape
        if(allocated(local_t%fields(node)%values)) deallocate(local_t%fields(node)%values)
        allocate(local_t%fields(node)%values(1,d(1),d(2),d(3)))
        local_t%fields(node)%values(1,:,:,:)=initialized%nodes(node)%temperature
      enddo
      call require(error<2e-12_dp,'partial FV/native-register endpoint oracle')
      if(pass>1) call require(control_gap>1e-8_dp,'partial physical/native flux difference is observed')
      if(mode==1) then
        call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,nested,[0,1],g,local,local_t,wide, &
          dt,pass,.false.,[.false.,.false.,.false.,.false.],0*force,1e-8_dp,1e-12_dp,faces, &
          b,w,ch,source,theta,MPI_COMM_WORLD,ok,mol_mode=.true.)
      else
        call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,nested,[0,1],g,local,local_t, &
          dt,pass,.false.,[.false.,.false.,.false.,.false.],0*force,1e-8_dp,1e-12_dp,wide,faces, &
          b,w,ch,source,theta,MPI_COMM_WORLD,ok,mol_mode=.true.)
      endif
      call require(ok,'partial MPI native registered reflux without diagnostic outputs')
      mpi_error=0
      do node=1,2
        if(rank/=node-1) cycle
        mpi_error=max(mpi_error,maxval(abs(local%fields(node)%values-expected%nodes(node)%state)/ &
          max(1._dp,abs(expected%nodes(node)%state))))
      enddo
      call MPI_Allreduce(mpi_error,global_error,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,status)
      call require(global_error<2e-12_dp,'partial MPI FV/native-register endpoint oracle')
      if(rank==0) print *, 'Partial MOL mode/iterations, native oracle, physical-control gap:',mode,pass,error,control_gap
    enddo
    enddo
  end subroutine check_partial_registered_reflux
  subroutine check_disabled_chemistry_rates()
    type(amr_tree_state_3d) :: ir,heat,serial_ir,serial_heat
    type(eb_integrated_flux_3d), allocatable :: registered(:),physical(:),serial_registered(:),serial_physical(:)
    real(dp) :: register_weight,registered_error,physical_error
    integer :: mode,pass,budget_iters
    allocate(ir%fields(1),heat%fields(1))
    if(rank==0) then
      allocate(ir%fields(1)%values(nv-4,2,2,2),heat%fields(1)%values(1,2,2,2))
      serial_ir=ir;serial_heat=heat
    endif
    allocate(registered(1),physical(1))
    if(rank==0) then
      allocate(registered(1)%x(nv,0:2,2,2),registered(1)%y(nv,2,0:2,2),registered(1)%z(nv,2,2,0:2))
      registered(1)%x=13;registered(1)%y=13;registered(1)%z=13
      physical=registered;serial_registered=registered;serial_physical=physical
    endif
    do budget_iters=1,3,2
    register_weight=merge(1._dp,.5_dp,budget_iters==1)
    do mode=1,2
    do pass=1,2
      state=initial;temp=initial_t;aux_calls=0;aux=31
      b=17;w=19;ch=23;source=29
      if(rank==0) then
        ir%fields(1)%values=1e20_dp;heat%fields(1)%values=37
        registered(1)%x=13;registered(1)%y=13;registered(1)%z=13
        physical(1)%x=17;physical(1)%y=17;physical(1)%z=17
      endif
      ! A very large stored rate would poison the predictor if reused.
      ! Pass 2 rejects only after both source stages, testing private reset.
      bad_species=5+findloc(y,0._dp,dim=1)
      if(mode==1) then
        if(pass==1) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
            1e-10_dp,budget_iters,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok, &
            reaction_source=ir,heat_release=heat,mol_mode=.true., &
            registered_flux=registered,integrated_flux=physical)
        else
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
            1e-10_dp,budget_iters,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok, &
            reaction_source=ir,heat_release=heat,mol_mode=.true., &
            registered_flux=registered,integrated_flux=physical, &
            auxiliary_source=negative_source,auxiliary_change=aux)
        endif
      else
        if(pass==1) then
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,[0],geometry,state,temp, &
            1e-10_dp,budget_iters,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,extended,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok, &
            reaction_source=ir,heat_release=heat,mol_mode=.true., &
            registered_flux=registered,integrated_flux=physical)
        else
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,tree,[0],geometry,state,temp, &
            1e-10_dp,budget_iters,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,extended,faces,b,w,ch,source,theta,MPI_COMM_WORLD,ok, &
            reaction_source=ir,heat_release=heat,mol_mode=.true., &
            registered_flux=registered,integrated_flux=physical, &
            auxiliary_source=negative_interval_source,auxiliary_change=aux)
        endif
      endif
      if(pass==1) then
        call require(ok,'disabled chemistry clears stored MPI rates')
        call require(all(ch==0),'disabled chemistry has no chemical inventory')
      else
        call require(.not.ok.and.aux_calls==2,'disabled chemistry late source rejection')
        call require(all(b==17).and.all(w==19).and.all(ch==23).and.all(source==29).and.all(aux==31), &
          'disabled chemistry ledger rollback')
      endif
      if(rank==0) then
        call require(all(state%fields(1)%values==initial%fields(1)%values), &
          'disabled chemistry does not use stored rates')
        call require(all(temp%fields(1)%values==initial_t%fields(1)%values),'disabled chemistry temperature')
        if(pass==1) then
          call require(all(ir%fields(1)%values==0).and.all(heat%fields(1)%values==0),'disabled MPI IR and heat reset')
          registered_error=max(maxval(abs(registered(1)%x(2,:,:,:)-register_weight*1e-10_dp*q(5))), &
            maxval(abs(registered(1)%y(3,:,:,:)-register_weight*1e-10_dp*q(5))), &
            maxval(abs(registered(1)%z(4,:,:,:)-register_weight*1e-10_dp*q(5))))/(1e-10_dp*q(5))
          physical_error=max(maxval(abs(physical(1)%x(2,:,:,:)-1e-10_dp*q(5))), &
            maxval(abs(physical(1)%y(3,:,:,:)-1e-10_dp*q(5))), &
            maxval(abs(physical(1)%z(4,:,:,:)-1e-10_dp*q(5))))/(1e-10_dp*q(5))
          call require(registered_error<1e-12_dp.and.physical_error<1e-12_dp, &
            'pinned MOL registration separated from physical flux')
        else
          call require(all(ir%fields(1)%values==1e20_dp).and.all(heat%fields(1)%values==37), &
            'disabled MPI IR and heat rollback')
          call require(all(registered(1)%x==13).and.all(registered(1)%y==13).and.all(registered(1)%z==13), &
            'registered flux rollback')
          call require(all(physical(1)%x==17).and.all(physical(1)%y==17).and.all(physical(1)%z==17), &
            'physical flux rollback')
        endif
        serial=serial_initial;serial_ir%fields(1)%values=1e20_dp;serial_heat%fields(1)%values=37
        if(mode==1) then
          call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-10_dp,budget_iters,.false., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial,extended,faces, &
            serial_b,serial_w,serial_ch,serial_source,theta,ok,reaction_source=serial_ir,heat_release=serial_heat, &
            mol_mode=.true.,registered_flux=serial_registered,integrated_flux=serial_physical)
        else
          call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-10_dp,budget_iters,.false., &
            [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial,extended,faces, &
            serial_b,serial_w,serial_ch,serial_source,theta,ok,reaction_source=serial_ir,heat_release=serial_heat, &
            mol_mode=.true.,registered_flux=serial_registered,integrated_flux=serial_physical)
        endif
        call require(ok,'disabled chemistry serial common/recursive')
        call require(all(serial_ir%fields(1)%values==0).and.all(serial_heat%fields(1)%values==0), &
          'disabled serial IR and heat reset')
        call require(all(serial%nodes(1)%state==initial%fields(1)%values),'disabled chemistry serial state')
        if(pass==1) then
          call require(all(serial_registered(1)%x==registered(1)%x).and. &
            all(serial_registered(1)%y==registered(1)%y).and.all(serial_registered(1)%z==registered(1)%z), &
            'serial/MPI registered flux agreement')
          call require(all(serial_physical(1)%x==physical(1)%x).and. &
            all(serial_physical(1)%y==physical(1)%y).and.all(serial_physical(1)%z==physical(1)%z), &
            'serial/MPI physical flux agreement')
        endif
        if(pass==2) then
          serial=serial_initial;serial_ir%fields(1)%values=1e20_dp;serial_heat%fields(1)%values=37
          serial_b=17;serial_w=19;serial_ch=23;serial_source=29;aux=31;serial_aux_calls=0
          serial_registered(1)%x=13;serial_registered(1)%y=13;serial_registered(1)%z=13
          serial_physical(1)%x=17;serial_physical(1)%y=17;serial_physical(1)%z=17
          if(mode==1) then
            call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-10_dp,budget_iters,.false., &
              [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial,extended,faces, &
              serial_b,serial_w,serial_ch,serial_source,theta,ok,reaction_source=serial_ir,heat_release=serial_heat, &
              mol_mode=.true.,registered_flux=serial_registered,integrated_flux=serial_physical, &
              auxiliary_source=serial_negative_source,auxiliary_change=aux)
          else
            call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-10_dp,budget_iters,.false., &
              [.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial,extended,faces, &
              serial_b,serial_w,serial_ch,serial_source,theta,ok,reaction_source=serial_ir,heat_release=serial_heat, &
              mol_mode=.true.,registered_flux=serial_registered,integrated_flux=serial_physical, &
              auxiliary_source=serial_negative_interval_source,auxiliary_change=aux)
          endif
          call require(.not.ok.and.serial_aux_calls==2,'disabled serial late rejection')
          call require(all(serial_ir%fields(1)%values==1e20_dp).and.all(serial_heat%fields(1)%values==37), &
            'disabled serial IR and heat rollback')
          call require(all(serial%nodes(1)%state==initial%fields(1)%values),'disabled serial state rollback')
          call require(all(serial_b==17).and.all(serial_w==19).and.all(serial_ch==23).and. &
            all(serial_source==29).and.all(aux==31),'disabled serial ledger rollback')
          call require(all(serial_registered(1)%x==13).and.all(serial_registered(1)%y==13).and. &
            all(serial_registered(1)%z==13),'serial registered flux rollback')
          call require(all(serial_physical(1)%x==17).and.all(serial_physical(1)%y==17).and. &
            all(serial_physical(1)%z==17),'serial physical flux rollback')
        endif
      endif
    enddo
    enddo
    enddo
    if(rank==0) print *, 'Disabled chemistry: common/recursive IR reset and MPI late rejection rollback'
    if(rank==0) print *, 'Pinned MOL registration weights 1 and 1/2; physical full flux; serial/MPI exact agreement'
  end subroutine check_disabled_chemistry_rates
  subroutine serial_negative_source(gas,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: gas
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    serial_aux_calls=serial_aux_calls+1
    rate%fields(1)%values=0
    if(stage==1) rate%fields(1)%values(bad_species,:,:,:)=-1e-40_dp
    good=.true.
  end subroutine serial_negative_source
  subroutine serial_negative_interval_source(gas,interval,stage,depth,time,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: gas
    real(dp), intent(in) :: interval,time
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    call serial_negative_source(gas,interval,stage,depth,rate,good)
  end subroutine serial_negative_interval_source
  subroutine negative_interval_source(gas,temperature,interval,stage,depth,time,rate,comm,good)
    real(dp), intent(in) :: time,interval
    type(amr_tree_state_3d), intent(in) :: gas,temperature
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    call negative_source(gas,temperature,interval,stage,depth,rate,comm,good)
  end subroutine negative_interval_source
  subroutine step(iterations,good)
    integer, intent(in) :: iterations
    logical, intent(out) :: good
    call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,tree,[0],geometry,state,temp,extended, &
      1e-10_dp,iterations,.true.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
      1e-8_dp,1e-12_dp,faces,b,w,ch,source,theta,MPI_COMM_WORLD,good,evaluate_reaction=reactor,mol_mode=.true.)
  end subroutine
  subroutine positive_source(gas,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: gas,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    integer :: c
    aux_calls=aux_calls+1
    rate=gas
    if(rank==0) then
      rate%fields(1)%values=0
      do c=1,3
        rate%fields(1)%values(c+1,:,:,:)=aux_force(c)
      enddo
    endif
    good=stage==0.or.stage==1
  end subroutine
  subroutine serial_positive_source(gas,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: gas
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    integer :: c
    serial_aux_calls=serial_aux_calls+1
    rate%fields(1)%values=0
    do c=1,3
      rate%fields(1)%values(c+1,:,:,:)=aux_force(c)
    enddo
    good=stage==0.or.stage==1
  end subroutine
  subroutine pulse_source(gas,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: gas,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    aux_calls=aux_calls+1;rate=gas
    if(rank==0) then
      rate%fields(1)%values=0
      if(stage==1) then
        rate%fields(1)%values(1,1,:,:)=1e5_dp
        rate%fields(1)%values(5+pulse_species,1,:,:)=1e5_dp
        rate%fields(1)%values(5,1,:,:)=1e5_dp*pulse_energy
      endif
    endif
    good=stage==0.or.stage==1
  end subroutine
  subroutine serial_pulse_source(gas,interval,stage,depth,rate,good)
    type(eb_reactive_tree_3d), intent(in) :: gas
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    logical, intent(out) :: good
    serial_aux_calls=serial_aux_calls+1;rate%fields(1)%values=0
    if(stage==1) then
      rate%fields(1)%values(1,1,:,:)=1e5_dp
      rate%fields(1)%values(5+pulse_species,1,:,:)=1e5_dp
      rate%fields(1)%values(5,1,:,:)=1e5_dp*pulse_energy
    endif
    good=stage==0.or.stage==1
  end subroutine
  subroutine pulse_reactor(sp,rx,e0,external_y,external_e,duration,rtol,atol,densities,temp,count,good,maximum_steps)
    type(nasa7_species), intent(in) :: sp(:)
    type(elementary_reaction), intent(in) :: rx(:)
    real(dp), intent(in) :: e0,external_y(:),external_e,duration,rtol,atol
    real(dp), intent(inout) :: densities(:),temp
    integer, intent(out) :: count
    logical, intent(out) :: good
    integer, optional, intent(in) :: maximum_steps
    type(elementary_reaction) :: no_reactions(0)
    call advance_forced_reactor_interval(sp,no_reactions,e0,external_y,external_e,duration,rtol,atol, &
      densities,temp,count,good,maximum_steps=maximum_steps)
    if(good) pulse_peak=max(pulse_peak,densities(pulse_species))
  end subroutine
  subroutine negative_source(gas,temperature,interval,stage,depth,rate,comm,good)
    type(amr_tree_state_3d), intent(in) :: gas,temperature
    real(dp), intent(in) :: interval
    integer, intent(in) :: stage,depth
    type(amr_tree_state_3d), intent(inout) :: rate
    type(MPI_Comm), intent(in) :: comm
    logical, intent(out) :: good
    aux_calls=aux_calls+1
    rate=gas
    if(rank==0) then
      rate%fields(1)%values=0
      if(stage>0) then
        rate%fields(1)%values(1,:,:,:)=-1e-8_dp
        rate%fields(1)%values(bad_species,:,:,:)=-1e-8_dp
      end if
    end if
    good=.true.
  end subroutine negative_source
  subroutine reactor(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps)
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    real(dp), intent(in) :: initial_energy_density,external_species(:),external_energy,duration,rtol,atol
    real(dp), intent(inout) :: species_density(:),temperature
    integer, intent(out) :: completed_steps
    logical, intent(out) :: ok
    integer, optional, intent(in) :: maximum_steps
    calls=calls+1
    if(fail_late.and.calls==9) then
      species_density=-1;temperature=-1;completed_steps=0;ok=.false.;return
    end if
    call advance_forced_reactor_interval(species,reactions,initial_energy_density,external_species,external_energy, &
      duration,rtol,atol,species_density,temperature,completed_steps,ok,maximum_steps=maximum_steps)
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: code
    if(condition) return
    print *, 'FAIL ',label
    call MPI_Abort(MPI_COMM_WORLD,1,code)
  end subroutine
end program
