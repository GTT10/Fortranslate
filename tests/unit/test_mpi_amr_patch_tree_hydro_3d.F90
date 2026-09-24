program test_mpi_amr_patch_tree_hydro_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d,advance_amr_patch_tree_hydro_3d
  use mpi_amr_patch_tree_hydro_3d_mod, only: advance_mpi_amr_tree_hydro_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  real(dp) :: pressure(6),mach(6),pin(6)
  character(len=32) :: faces(6),reconstruction
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(4,4,4,2,2,2,2,2,2,3,patches(2),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(global%fields(3),reference%fields(3))
  do node=1,3
    n=tree%nodes(node)%shape
    allocate(global%fields(node)%values(12,n(1),n(2),n(3)),reference%fields(node)%values(1,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          q=0;q(1)=.3_dp*(1+.02_dp*(i+j+k));q(5)=101325*(1+.02_dp*node)
          q(2:4)=[10.0_dp+i,20.0_dp+j,-4.0_dp+.5_dp*k]
          q(6)=.025_dp+.003_dp*i;q(9)=.22_dp;q(12)=1-q(6)-q(9)
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
  do mode=1,8
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='slip_wall'
    end if
    if (mode==4) faces(1)='inflow'
    pressure=0;mach=-1;pin=0
    if (mode>=5) pressure([1,6])=102325
    if(mode>=7) then
      pressure=0;faces(1)='inflow';mach(1)=.2_dp
      if(mode==8) pin(1)=120000
    end if
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    oracle=global;oracle_t=reference;dt=.1_dp*expected_dt
    call advance_amr_patch_tree_hydro_3d(tree,species,oracle,oracle_t,spacing,dt,faces,reconstruction,'mc','hllc', &
      serial_budget,serial_largest,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'serial RK2 oracle')
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc','hllc', &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'distributed RK2')
    do node=1,3
      call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owned state only')
      call require(allocated(temperature%fields(node)%values).eqv.(owners(node)==rank),'owned temperature only')
      if (owners(node)/=rank) cycle
      call require(all(local%fields(node)%values==oracle%fields(node)%values),'exact serial conserved state')
      call require(all(temperature%fields(node)%values==oracle_t%fields(node)%values),'exact serial temperature')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest,'exact boundary and reflux diagnostics')
    saved=local;saved_t=temperature
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc','bad', &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'flux solver failure rejects collectively')
    call unchanged()
    if (rank==nranks-1) dt=2*dt
    call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc','hllc', &
      budget,largest,MPI_COMM_WORLD,ok,inflow_state=inlet,inflow_temperature=inlet_t,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'mismatched timestep rejected')
    call unchanged()
    if (mode>=5.and.mode<=6) then
      dt=.1_dp*expected_dt
      if (rank==nranks-1) pressure(1)=103325
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc','hllc', &
        budget,largest,MPI_COMM_WORLD,ok,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'mismatched pressure rejected')
      call unchanged()
      pressure(1)=102325
      if (rank==nranks-1) pressure(3)=102325
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces,reconstruction,'mc','hllc', &
        budget,largest,MPI_COMM_WORLD,ok,outlet_pressure=pressure, &
      inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'rank-local wall pressure rejected')
      call unchanged()
    end if
    if(mode>=7) then
      dt=.1_dp*expected_dt
      if(rank==nranks-1) mach(1)=.3_dp
      call advance_mpi_amr_tree_hydro_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
        reconstruction,'mc','hllc',budget,largest,MPI_COMM_WORLD,ok,inflow_temperature=inlet_t, &
        inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'rank-disagreeing Mach rejected')
      call unchanged()
    end if
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: distributed tree RK2, PCM/PLM, physical boundaries, exact serial parity and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest,'diagnostic rollback')
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
