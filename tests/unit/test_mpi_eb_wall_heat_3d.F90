program test_mpi_eb_wall_heat_3d
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use mixture_thermo_mod, only: mixture_mass_properties
  use reactive_1d_mod, only: reactive_nvar,reactive_nprim,reactive_mass_fraction_component,reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d,extract_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    build_eb_tree_extended_geometry_3d,compute_eb_tree_wall_transport_3d, &
    advance_reactive_eb_tree_sdc_3d,advance_reactive_eb_tree_sdc_subcycle_3d
  use eb_reactive_reflux_3d_mod, only: eb_integrated_flux_3d
  use eb_wall_heat_source_3d_mod, only: compute_eb_wall_heat_source_3d,compute_eb_wall_transport_source_3d
  use mpi_eb_wall_heat_3d_mod, only: compute_mpi_eb_wall_heat_3d,compute_mpi_eb_wall_transport_3d
  use mpi_eb_tree_transport_3d_mod, only: compute_sdc_mpi_eb_tree_transport_3d
  use mpi_eb_tree_full_3d_mod, only: advance_mpi_eb_tree_sdc_3d,advance_mpi_eb_tree_sdc_subcycle_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  use elementary_kinetics_mod, only: elementary_reaction
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:)
  type(amr_patch_3d) :: boxes(2)
  type(amr_patch_tree_3d) :: topology
  type(eb_reactive_tree_3d) :: serial,visc_serial,serial_result
  type(eb_geometry_3d) :: padded,finest,geometry(3),padded_six
  type(eb_geometry_3d), allocatable :: extended(:),extended_six(:)
  type(amr_tree_state_3d) :: local,lt,rate,stage,baseline
  type(amr_tree_state_3d) :: evolved,evolved_t
  type(amr_tree_state_3d) :: visc_state,visc_rate,serial_rate,saved_rate
  type(eb_integrated_flux_3d), allocatable :: flux(:)
  real(dp), allocatable :: q(:),y(:),u(:),root(:,:,:,:),rt(:,:,:),wide(:,:,:,:),wt(:,:,:),heat(:,:,:)
  real(dp) :: mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,t,sound,expected,ref_frequency,power,frequency
  real(dp) :: theta,wall_energy,stage_frequency,control
  real(dp), allocatable :: before(:),after(:),boundary(:),wall(:),chemical(:),external(:)
  real(dp), allocatable :: visc_reference(:),visc_exchange(:),visc_source(:,:,:,:)
  real(dp) :: residual,serial_frequency
  real(dp), allocatable :: serial_exchange(:)
  character(len=16) :: faces(6)
  integer :: status,rank,ranks,nv,s,node,owners(3),n(3),mode,scheme,j
  logical :: ok
  call MPI_Init(status)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,status)
  call MPI_Comm_size(MPI_COMM_WORLD,ranks,status)
  call load_h2o2_elementary_thermo(species,ok);call require(ok,'thermo')
  call load_h2o2_elementary_transport(transport,ok);call require(ok,'transport')
  nv=reactive_nvar(size(species))
  allocate(reactions(0),before(nv),after(nv),boundary(nv),wall(nv),chemical(nv),external(nv))
  allocate(visc_reference(nv),visc_exchange(nv),visc_source(nv,8,8,8),serial_exchange(nv))
  faces='outflow'
  allocate(q(reactive_nprim(size(species))),y(size(species)),u(nv),root(nv,4,4,4),rt(4,4,4))
  y=0
  do s=1,size(species)
    if(species(s)%name=='N2') y(s)=1
  end do
  call mixture_mass_properties(species,y,500._dp,mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,ok)
  call require(ok,'properties')
  q=0;q(1)=1;q(5)=rgas*500
  do s=1,size(species)
    q(reactive_mass_fraction_component(s))=y(s)
  end do
  call reactive_primitive_to_conserved(species,q,u,t,sound,ok);call require(ok,'state')
  do s=1,nv
    root(s,:,:,:)=u(s)
  end do
  rt=500
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,boxes(1),ok);call require(ok,'left patch')
  call initialize_amr_patch_3d(4,4,4,3,4,1,4,1,4,2,boxes(2),ok);call require(ok,'right patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1],boxes,topology,ok);call require(ok,'tree')
  call build_axis_plane_eb_geometry_3d(16,16,16,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp, &
    'x',.43_dp,padded,ok);call require(ok,'padded geometry')
  call extract_eb_geometry_3d(padded,[5,5,5],[12,12,12],finest,ok);call require(ok,'finest')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,rt,serial,ok);call require(ok,'serial tree')
  do node=1,3
    geometry(node)=serial%nodes(node)%geometry
  end do
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded,extended,ok,halo_width=2)
  call require(ok,'extended tree geometry')
  call build_axis_plane_eb_geometry_3d(32,32,32,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
    'x',.43_dp,padded_six,ok);call require(ok,'hydro padded geometry')
  call build_eb_tree_extended_geometry_3d(topology,geometry,padded_six,extended_six,ok)
  call require(ok,'hydro extended tree geometry')
  allocate(wide(nv,16,16,16),wt(16,16,16),heat(8,8,8))
  do s=1,nv
    wide(s,:,:,:)=u(s)
  end do
  wt=500
  call compute_eb_wall_heat_source_3d(species,transport,wide,wt,padded,[5,5,5],[12,12,12], &
    400._dp,heat,expected,ref_frequency,ok);call require(ok,'whole fine reference')
  wide(3,:,:,:)=wide(3,:,:,:)+2*wide(1,:,:,:)
  wide(5,:,:,:)=wide(5,:,:,:)+2*wide(1,:,:,:)
  call compute_eb_wall_transport_source_3d(species,transport,wide,wt,padded,[5,5,5],[12,12,12], &
    .false.,.true.,0._dp,visc_source,visc_reference,ref_frequency,ok)
  call require(ok.and.visc_reference(3)<0,'whole fine stationary-wall drag')
  call require(visc_reference(1)==0.and.all(visc_reference(5:)==0),'stationary wall zero mass and energy')
  visc_serial=serial
  do node=1,3
    visc_serial%nodes(node)%state(3,:,:,:)=visc_serial%nodes(node)%state(3,:,:,:)+2*visc_serial%nodes(node)%state(1,:,:,:)
    visc_serial%nodes(node)%state(5,:,:,:)=visc_serial%nodes(node)%state(5,:,:,:)+2*visc_serial%nodes(node)%state(1,:,:,:)
  end do
  call compute_eb_tree_wall_transport_3d(species,transport,serial,extended,.true.,.false.,400._dp, &
    serial_rate,serial_exchange,serial_frequency,ok)
  call require(ok,'serial composite wall source')
  saved_rate=serial_rate;serial_exchange=123;serial_frequency=456
  call compute_eb_tree_wall_transport_3d(species,transport,serial,extended,.true.,.false.,-1._dp, &
    serial_rate,serial_exchange,serial_frequency,ok)
  call require(.not.ok.and.all(serial_exchange==123).and.serial_frequency==456,'serial source ledger rollback')
  do node=1,3
    call require(all(serial_rate%fields(node)%values==saved_rate%fields(node)%values),'serial source rate rollback')
  end do
  call compute_eb_tree_wall_transport_3d(species,transport,serial,extended,.true.,.false.,400._dp, &
    serial_rate,serial_exchange,serial_frequency,ok)
  call require(ok,'serial source after rejection')
  call check_held_wall()
  do mode=1,2
    owners=[0,0,mod(1,ranks)]
    if(mode==2) owners=0
    allocate(local%fields(3),lt%fields(3))
    do node=1,3
      if(owners(node)/=rank) cycle
      n=topology%nodes(node)%shape
      local%fields(node)%values=serial%nodes(node)%state
      allocate(lt%fields(node)%values(1,n(1),n(2),n(3)))
      lt%fields(node)%values(1,:,:,:)=serial%nodes(node)%temperature
    end do
    call compute_mpi_eb_wall_heat_3d(species,transport,topology,owners,geometry,extended,local,lt,400._dp, &
      rate,power,frequency,MPI_COMM_WORLD,ok)
    call require(ok,'distributed heat')
    call require(abs(power-expected)<1e-12_dp*abs(expected),'composite heat excludes covered coarse wall')
    call require(frequency>0.and.frequency==serial_frequency,'exact serial MPI wall frequency')
    call require(abs(power-serial_exchange(5))<1e-12_dp*abs(power),'serial MPI composite heat')
    do node=1,3
      if(owners(node)/=rank) cycle
      call require(all(rate%fields(node)%values==serial_rate%fields(node)%values),'exact serial MPI wall rates')
    end do
    visc_state=local
    do node=1,3
      if(owners(node)/=rank) cycle
      visc_state%fields(node)%values(3,:,:,:)=visc_state%fields(node)%values(3,:,:,:)+ &
        2*visc_state%fields(node)%values(1,:,:,:)
      visc_state%fields(node)%values(5,:,:,:)=visc_state%fields(node)%values(5,:,:,:)+ &
        2*visc_state%fields(node)%values(1,:,:,:)
    end do
    call compute_mpi_eb_wall_transport_3d(species,transport,topology,owners,geometry,extended,visc_state,lt, &
      .false.,.true.,0._dp,visc_rate,visc_exchange,stage_frequency,MPI_COMM_WORLD,ok)
    call require(ok,'distributed stationary-wall drag')
    call require(maxval(abs(visc_exchange-visc_reference))<1e-12_dp*maxval(abs(visc_reference)), &
      'composite traction excludes covered coarse wall')
    call require(visc_exchange(1)==0.and.all(visc_exchange(5:)==0),'MPI zero mass and energy')
    if(rank==0) print *, 'Stationary wall force:',visc_exchange(2:4)
    call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1e-7_dp, &
      .false.,.true.,.false.,.false.,local,lt,baseline,flux,theta,MPI_COMM_WORLD,ok)
    call require(ok,'baseline MPI stage')
    call compute_sdc_mpi_eb_tree_transport_3d(species,transport,topology,owners,geometry,1e-7_dp, &
      .false.,.true.,.false.,.false.,local,lt,stage,flux,theta,MPI_COMM_WORLD,ok, &
      embedded_wall_temperature=400._dp,extended_geometry=extended,embedded_wall_energy=wall_energy, &
      embedded_frequency=stage_frequency)
    call require(ok,'heated MPI stage')
    call require(abs(wall_energy-1e-7_dp*power)<1e-20_dp,'integrated stage ledger')
    do node=1,3
      if(owners(node)/=rank) cycle
      call require(maxval(abs(stage%fields(node)%values-baseline%fields(node)%values-rate%fields(node)%values)) &
        <1e-10_dp,'source routed once')
    end do
    control=400
    if(rank==ranks-1) control=-1
    power=123;frequency=456
    call compute_mpi_eb_wall_heat_3d(species,transport,topology,owners,geometry,extended,local,lt,control, &
      rate,power,frequency,MPI_COMM_WORLD,ok)
    call require(.not.ok.and.power==123.and.frequency==456,'collective failure rollback')
    if(mode==1) then
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,before,MPI_COMM_WORLD,ok)
      call require(ok,'initial inventory')
      do scheme=1,4
        evolved=local;evolved_t=lt
        boundary=0;wall=0;chemical=0;external=0;theta=1
        if(mod(scheme,2)==1) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
            extended_six,1e-7_dp,2,.false.,[.false.,.true.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
            embedded_wall_temperature=400._dp,mol_mode=scheme==3)
        else
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
            1e-7_dp,2,.false.,[.false.,.true.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,extended_six,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
            embedded_wall_temperature=400._dp,mol_mode=scheme==4)
        end if
        call require(ok,'evolved isothermal wall')
        call compare_serial_wall(.true.,scheme)
        call integrate_mpi_eb_tree_3d(topology,owners,geometry,evolved,after,MPI_COMM_WORLD,ok)
        call require(ok,'final inventory')
        residual=maxval(abs(after-before-boundary-wall-chemical-external)/max(1._dp,abs(before)))
        call require(residual<1e-12_dp,'time-integrated conservation')
        call require(wall(5)<0.and.abs(wall(5)-1e-7_dp*expected)<1e-5_dp*abs(1e-7_dp*expected), &
          'time-integrated cooling ledger')
        if(rank==0) print *, 'Wall scheme, heat, budget residual:',scheme,wall(5),residual
      end do
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,visc_state,before,MPI_COMM_WORLD,ok)
      call require(ok,'initial viscous inventory')
      do scheme=1,4
        evolved=visc_state;evolved_t=lt
        boundary=0;wall=0;chemical=0;external=0;theta=1
        if(mod(scheme,2)==1) then
          call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
            extended_six,1e-7_dp,2,.false.,[.true.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
            embedded_noslip=.true.,mol_mode=scheme==3)
        else
          call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
            1e-7_dp,2,.false.,[.true.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
            1e-8_dp,1e-12_dp,extended_six,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
            embedded_noslip=.true.,mol_mode=scheme==4)
        end if
        call require(ok,'evolved adiabatic stationary wall')
        call compare_serial_wall(.false.,scheme)
        call integrate_mpi_eb_tree_3d(topology,owners,geometry,evolved,after,MPI_COMM_WORLD,ok)
        call require(ok,'final viscous inventory')
        residual=maxval(abs(after-before-boundary-wall-chemical-external)/max(1._dp,abs(before)))
        call require(residual<1e-12_dp,'viscous time-integrated conservation')
        call require(wall(5)==0,'stationary adiabatic wall energy')
        call require(wall(3)<0.and.abs(wall(3)-1e-7_dp*visc_reference(3)) &
          <1e-5_dp*abs(1e-7_dp*visc_reference(3)),'time-integrated tangential drag ledger')
        if(rank==0) print *, 'Wall scheme, tangential impulse, budget residual:',scheme,wall(3),residual
      end do
    end if
    if(mode==1) then
      evolved=visc_state;evolved_t=lt
      do node=1,3
        if(owners(node)/=rank) cycle
        do j=1,topology%nodes(node)%shape(2)
          evolved%fields(node)%values(:,:,j,:)=evolved%fields(node)%values(:,:,j,:)* &
            (1+.01_dp*(geometry(node)%y_lower+(j-.5_dp)*geometry(node)%dy))
        end do
      end do
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,evolved,before,MPI_COMM_WORLD,ok)
      call require(ok,'MOL inviscid initial inventory')
      call advance_mpi_eb_tree_sdc_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
        extended_six,1e-7_dp,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-8_dp,1e-12_dp,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok,mol_mode=.true.)
      call require(ok,'MOL inviscid nonuniform step')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,evolved,after,MPI_COMM_WORLD,ok)
      call require(ok,'MOL inviscid final inventory')
      residual=maxval(abs(after-before-boundary-wall-chemical-external)/max(1._dp,abs(before)))
      call require(residual<1e-12_dp,'MOL inviscid conservation')
      call require(any(after/=before),'MOL nontrivial update')
      if(rank==0) print *, 'MOL inviscid nonuniform budget residual:',residual
      before=after
      call advance_mpi_eb_tree_sdc_subcycle_3d(species,reactions,transport,topology,owners,geometry,evolved,evolved_t, &
        1e-7_dp,2,.false.,[.false.,.false.,.false.,.false.],[0._dp,0._dp,0._dp], &
        1e-8_dp,1e-12_dp,extended_six,faces,boundary,wall,chemical,external,theta,MPI_COMM_WORLD,ok, &
        redistribution='state_redist',mol_mode=.true.)
      call require(ok,'recursive redistributed MOL step')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,evolved,after,MPI_COMM_WORLD,ok)
      call require(ok,'recursive MOL final inventory')
      residual=maxval(abs(after-before-boundary-wall-chemical-external)/max(1._dp,abs(before)))
      call require(residual<1e-12_dp,'recursive redistributed MOL conservation')
      if(rank==0) print *, 'Recursive redistributed MOL budget residual:',residual
    end if
    deallocate(local%fields,lt%fields)
  end do
  if(rank==0) print *, 'PASS: cross-owner donors, composite wall heat, SDC source and collective rejection',expected
  call MPI_Finalize(status)
contains
  subroutine check_held_wall()
    type(eb_reactive_tree_3d) :: trial
    type(amr_tree_state_3d) :: stored,held,gas,temperature,reference,pending,previous
    real(dp) :: held_exchange(nv),reference_exchange(nv),previous_exchange(nv),held_frequency,reference_frequency
    real(dp) :: previous_frequency,held_power
    integer :: layout,p,dims(3),partition(3)
    call mixture_mass_properties(species,y,500._dp,mw,rgas,cp,cv,gamma,enthalpy,energy,entropy,ok)
    call require(ok,'held wall independent NASA')
    trial=visc_serial;allocate(stored%fields(3))
    do p=1,3
      dims=topology%nodes(p)%shape
      allocate(stored%fields(p)%values(1,dims(1),dims(2),dims(3)),source=energy)
      trial%nodes(p)%state(5,:,:,:)=-1.e20_dp
    enddo
    call compute_eb_tree_wall_transport_3d(species,transport,trial,extended,.true.,.true.,400._dp, &
      reference,reference_exchange,reference_frequency,ok,internal_energy=stored)
    call require(ok,'serial actual held wall family')
    call require(abs(reference_exchange(5)-expected)<1.e-12_dp*abs(expected), &
      'independent composite held heat ledger')
    call require(maxval(abs(reference_exchange(2:4)-visc_reference(2:4)))< &
      1.e-12_dp*maxval(abs(visc_reference)),'independent composite held traction')
    do layout=1,2
      partition=[0,0,mod(1,ranks)]
      if(layout==2) partition=0
      allocate(held%fields(3),gas%fields(3),temperature%fields(3))
      do p=1,3
        if(partition(p)/=rank) cycle
        dims=topology%nodes(p)%shape
        held%fields(p)=stored%fields(p);gas%fields(p)%values=trial%nodes(p)%state
        allocate(temperature%fields(p)%values(1,dims(1),dims(2),dims(3)),source=500._dp)
      enddo
      call compute_mpi_eb_wall_transport_3d(species,transport,topology,partition,geometry,extended,gas,temperature, &
        .true.,.true.,400._dp,pending,held_exchange,held_frequency,MPI_COMM_WORLD,ok,internal_energy=held)
      call require(ok,'distributed actual held thermal/viscous wall')
      call require(all(held_exchange==reference_exchange).and.held_frequency==reference_frequency, &
        'exact composite held wall ledgers/frequency')
      do p=1,3
        if(partition(p)/=rank) cycle
        call require(all(pending%fields(p)%values==reference%fields(p)%values),'exact owner held wall rates')
      enddo
      previous=pending;previous_exchange=held_exchange;previous_frequency=held_frequency
      if(partition(3)==rank) held%fields(3)%values(1,1,1,1)=ieee_value(0._dp,ieee_quiet_nan)
      call compute_mpi_eb_wall_transport_3d(species,transport,topology,partition,geometry,extended,gas,temperature, &
        .true.,.true.,400._dp,pending,held_exchange,held_frequency,MPI_COMM_WORLD,ok,internal_energy=held)
      call require(.not.ok.and.all(held_exchange==previous_exchange).and.held_frequency==previous_frequency, &
        'remote held wall donor refusal ledgers')
      do p=1,3
        if(partition(p)/=rank) cycle
        call require(all(pending%fields(p)%values==previous%fields(p)%values),'remote held wall rate rollback')
      enddo
      if(partition(3)==rank) held%fields(3)=stored%fields(3)
      call compute_mpi_eb_wall_heat_3d(species,transport,topology,partition,geometry,extended,gas,temperature,400._dp, &
        pending,held_power,held_frequency,MPI_COMM_WORLD,ok,internal_energy=held)
      call require(ok.and.held_power==reference_exchange(5),'held heat-only MPI API')
      deallocate(held%fields,gas%fields,temperature%fields)
    enddo
  end subroutine check_held_wall
  subroutine compare_serial_wall(heating,scheme)
    logical, intent(in) :: heating
    integer, intent(in) :: scheme
    real(dp), allocatable :: target_temperature
    real(dp) :: b(nv),w(nv),ch(nv),ext(nv),limiter,error
    logical :: good,flags(4)
    integer :: p
    flags=[.not.heating,heating,.false.,.false.]
    if(heating) then
      serial_result=serial
      allocate(target_temperature,source=400._dp)
    else
      serial_result=visc_serial
    end if
    if(mod(scheme,2)==1) then
      call advance_reactive_eb_tree_sdc_3d(species,reactions,transport,1e-7_dp,2,.false.,flags, &
        [0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial_result,extended_six,faces,b,w,ch,ext,limiter,good, &
        embedded_wall_temperature=target_temperature,embedded_noslip=.not.heating,mol_mode=scheme==3)
    else
      call advance_reactive_eb_tree_sdc_subcycle_3d(species,reactions,transport,1e-7_dp,2,.false.,flags, &
        [0._dp,0._dp,0._dp],1e-8_dp,1e-12_dp,serial_result,extended_six,faces,b,w,ch,ext,limiter,good, &
        embedded_wall_temperature=target_temperature,embedded_noslip=.not.heating,mol_mode=scheme==4)
    end if
    call require(good,'serial embedded wall evolution')
    do p=1,3
      if(owners(p)/=rank) cycle
      call require(all(serial_result%nodes(p)%state==evolved%fields(p)%values),'exact serial MPI heated/viscous state')
      call require(all(serial_result%nodes(p)%temperature==evolved_t%fields(p)%values(1,:,:,:)), &
        'exact serial MPI heated/viscous temperature')
    end do
    error=max(maxval(abs(b-boundary)/max(1._dp,abs(boundary))),maxval(abs(w-wall)/max(1._dp,abs(wall))), &
      maxval(abs(ch-chemical)/max(1._dp,abs(chemical))),maxval(abs(ext-external)/max(1._dp,abs(external))))
    call require(error<1e-12_dp.and.limiter==theta,'serial MPI wall ledgers and limiter')
  end subroutine compare_serial_wall
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    integer :: code
    if(condition) return
    print *, 'FAIL: ',message
    call MPI_Abort(MPI_COMM_WORLD,1,code)
  end subroutine
end program
