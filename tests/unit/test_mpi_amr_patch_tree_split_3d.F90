program test_mpi_amr_patch_tree_split_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d
  use transport_database_mod, only: gas_transport_species,load_h2o2_elementary_transport
  use amr_patch_tree_split_3d_mod, only: advance_amr_patch_tree_split_3d,advance_amr_patch_tree_sdc_3d
  use mpi_amr_patch_tree_split_3d_mod, only: advance_mpi_amr_tree_split_3d,advance_mpi_amr_tree_sdc_3d
  use elementary_kinetics_mod, only: elementary_reaction
  use h2o2_elementary_mechanism_mod, only: load_h2o2_elementary_mechanism
  use amr_patch_tree_chemistry_3d_mod, only: advance_amr_patch_tree_chemistry_3d,advance_amr_patch_tree_strang_3d
  use mpi_amr_patch_tree_chemistry_3d_mod, only: advance_mpi_amr_tree_chemistry_3d,advance_mpi_amr_tree_strang_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(gas_transport_species), allocatable :: transport(:)
  type(elementary_reaction), allocatable :: reactions(:),bad(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  real(dp) :: pressure(6),mach(6),pin(6)
  character(len=32) :: faces(6),reconstruction,integrator,arg,solver
  real(dp) :: coefficient
  logical :: ppm
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth
  logical :: ok,enabled,bad_enabled
  real(dp) :: theta,serial_theta,wall_t(6),wall_v(3,6)
  call get_command_argument(1,arg)
  ppm=trim(arg)=='ppm'
  coefficient=0;solver='hllc'
  if(ppm) then
    coefficient=.037_dp;solver='pelec'
  end if
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call load_h2o2_elementary_mechanism(reactions,ok)
  call require(ok,'mechanism')
  call load_h2o2_elementary_transport(transport,ok)
  call require(ok,'transport')
  call initialize_amr_patch_3d(2,2,2,1,1,1,1,1,1,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(2,2,2,2,2,2,2,2,2,3,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([2,2,2],[1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(global%fields(3),reference%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(global%fields(node)%values(12,n(1),n(2),n(3)),reference%fields(node)%values(1,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=202650*(1+.02_dp*node)
          q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
          q(6)=.025_dp+.003_dp*i;q(7)=1e-5_dp;q(9)=.22_dp;q(12)=1-q(6)-q(7)-q(9)
          call reactive_primitive_to_conserved(species,q,global%fields(node)%values(:,i,j,k),exact_t,sound,ok)
          call require(ok,'initial conserved state')
        end do
      end do
    end do
    call recover_reactive_temperatures_3d(species,global%fields(node)%values, &
      spread(spread(spread(800.0_dp,1,n(1)),2,n(2)),3,n(3)),n(1),n(2),n(3),reference%fields(node)%values(1,:,:,:),ok)
    call require(ok,'serial recovery oracle')
  end do
  spacing=[.001_dp,.0015_dp,.002_dp]
  call compute_amr_patch_tree_cfl_3d(tree,species,global,reference,spacing,.15_dp,expected_dt,ok)
  call require(ok,'serial CFL')
  inlet=spread(global%fields(1)%values(:,1,1,1),2,6)
  inlet_t=reference%fields(1)%values(1,1,1,1)
  if(trim(arg)=='sdc') then
    call check_sdc()
    if(rank==0) print *, 'PASS: coupled MPI SDC parity and rollback'
    call MPI_Finalize(ierr)
    stop
  end if
  do mode=1,6
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if(ppm) reconstruction='characteristic_ppm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='slip_wall'
    end if
    if (mode==3) faces(1)='inflow'
    integrator='explicit'
    if (mod(mode,2)==0) integrator='implicit'
    wall_t=1100;wall_v=0;enabled=mode/=3
    if (mod(mode,2)==0) then
      faces='outflow';faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
    end if
    pressure=0;mach=-1;pin=0
    if (mode==4) pressure([1,6])=203650
    if(mode>=5) then
      faces(1)='inflow';mach(1)=.2_dp
      if(mode==6) pin(1)=210000
    end if
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    oracle=global;oracle_t=reference;dt=1e-10_dp
    call advance_amr_patch_tree_split_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces, &
      reconstruction,'mc',solver,enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      serial_budget,serial_largest,serial_theta,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6,difmag=coefficient)
    call require(ok,'serial full operator split')
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc',solver,enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6,difmag=coefficient)
    call require(ok,'distributed full operator split')
    call compare_state()
    call require(all(budget==serial_budget).and.largest==serial_largest.and.theta==serial_theta,'exact coupled diagnostics')
    saved=local;saved_t=temperature
    if(ppm) then
      if(rank==nranks-1) coefficient=.05_dp
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc',solver,enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6,difmag=coefficient)
      call require(.not.ok,'coefficient mismatch rolls back chemistry and transport')
      call unchanged()
      coefficient=.037_dp
    end if
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc','bad',enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6)
    call require(.not.ok,'hydro failure rolls back preceding chemistry and transport')
    call unchanged()
    bad_enabled=enabled
    if (rank==nranks-1) bad_enabled=.not.enabled
    call advance_mpi_amr_tree_split_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces, &
      reconstruction,'mc',solver,bad_enabled,enabled,enabled,enabled,enabled,1e-6_dp,1e-14_dp, &
      budget,largest,theta,MPI_COMM_WORLD,ok,integrator,inlet,inlet_t,wall_t,wall_v, &
      outlet_pressure=pressure,inflow_mach=mach,inflow_pressure=pin,subcycling=mode==6)
    call require(.not.ok,'inconsistent enabled operators rejected')
    call unchanged()
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: distributed reaction transport hydro split, disabled operators, serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine check_sdc()
    type(amr_tree_state_3d) :: ir,serial_ir,retained_ir
    real(dp) :: external(12),serial_external(12),force(3),visc
    integer :: peer,dims(3)
    owners=mod([0,1,0],nranks);dt=1e-10_dp;faces='periodic';force=[20000._dp,-30000._dp,40000._dp]
    allocate(local%fields(3),temperature%fields(3),ir%fields(3),serial_ir%fields(3))
    do peer=1,3
      dims=tree%nodes(peer)%shape
      allocate(serial_ir%fields(peer)%values(8,dims(1),dims(2),dims(3)))
      serial_ir%fields(peer)%values=0
      if(owners(peer)/=rank) cycle
      local%fields(peer)=global%fields(peer);temperature%fields(peer)=reference%fields(peer)
      ir%fields(peer)=serial_ir%fields(peer)
    end do
    oracle=global;oracle_t=reference
    call advance_amr_patch_tree_sdc_3d(tree,species,reactions,transport,oracle,oracle_t,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,serial_budget,serial_external,serial_largest,serial_theta, &
      ok,reaction_source=serial_ir)
    call require(ok,'serial coupled SDC reference')
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok,reaction_source=ir)
    call require(ok,'MPI coupled SDC iteration')
    call compare_state()
    call require(all(budget==serial_budget).and.all(external==serial_external),'exact MPI SDC boundary/source inventory')
    call require(largest==serial_largest.and.theta==serial_theta,'exact MPI SDC diagnostics')
    do peer=1,3
      call require(allocated(ir%fields(peer)%values).eqv.(owners(peer)==rank),'owner-local SDC estimate')
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==serial_ir%fields(peer)%values),'exact MPI SDC reaction estimates')
    end do
    saved=local;saved_t=temperature;retained_ir=ir;visc=.1_dp
    if(rank==nranks-1) visc=-1
    call advance_mpi_amr_tree_sdc_3d(tree,owners,species,reactions,transport,local,temperature,spacing,dt,faces,2,.true., &
      [.true.,.true.,.true.,.true.],force,1e-6_dp,1e-14_dp,budget,external,largest,theta,MPI_COMM_WORLD,ok, &
      reaction_source=ir,difmag=visc)
    call require(.not.ok,'idle-rank invalid hydro after transport rejects whole SDC step')
    call unchanged()
    call require(all(external==serial_external),'source inventory rollback')
    do peer=1,3
      if(owners(peer)/=rank) cycle
      call require(all(ir%fields(peer)%values==retained_ir%fields(peer)%values),'reaction estimate rollback')
    end do
  end subroutine check_sdc
  subroutine compare_state()
    integer :: peer
    do peer=1,3
      call require(allocated(local%fields(peer)%values).eqv.(owners(peer)==rank),'owned state only')
      call require(allocated(temperature%fields(peer)%values).eqv.(owners(peer)==rank),'owned temperature only')
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==oracle%fields(peer)%values),'exact serial conserved state')
      call require(all(temperature%fields(peer)%values==oracle_t%fields(peer)%values),'exact serial temperature')
    end do
  end subroutine
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest.and.theta==serial_theta,'diagnostic rollback')
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if (condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program
