program test_mpi_amr_patch_tree_thermo_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: recover_mpi_amr_tree_temperatures_3d,compute_mpi_amr_tree_cfl_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
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
  do mode=1,3
    owners=[0,1,0]
    if (mode==2) owners=0
    if (mode==3) owners=mod([0,1,2],nranks)
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
      temperature%fields(node)%values=800
    end do
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,local,temperature,MPI_COMM_WORLD,ok,depth=1)
    call require(ok,'single depth recovery')
    do node=1,3
      if (owners(node)/=rank) cycle
      if (node==2) then
        call require(all(temperature%fields(node)%values==reference%fields(node)%values),'selected depth exact')
      else
        call require(all(temperature%fields(node)%values==800),'unselected depths unchanged')
      end if
      temperature%fields(node)%values=800
    end do
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,local,temperature,MPI_COMM_WORLD,ok)
    call require(ok,'all-depth recovery')
    do node=1,3
      call require(allocated(temperature%fields(node)%values).eqv.(owners(node)==rank),'owner-only temperatures')
      if (owners(node)/=rank) cycle
      call require(all(temperature%fields(node)%values==reference%fields(node)%values),'exact serial temperature')
      call require(all(local%fields(node)%values==global%fields(node)%values),'conserved source unchanged')
    end do
    call compute_mpi_amr_tree_cfl_3d(tree,owners,species,local,temperature,spacing,.15_dp,dt,MPI_COMM_WORLD,ok)
    call require(ok.and.dt==expected_dt,'exact global CFL')
    saved=temperature
    if (rank==owners(3)) local%fields(3)%values(1,1,1,1)=-1
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,local,temperature,MPI_COMM_WORLD,ok)
    call require(.not.ok,'failed EOS on one owner rejects every rank')
    call unchanged()
    call compute_mpi_amr_tree_cfl_3d(tree,owners,species,local,temperature,spacing,.15_dp,dt,MPI_COMM_WORLD,ok)
    call require(.not.ok.and.dt==expected_dt,'invalid EOS preserves CFL output')
    if (rank==owners(3)) local%fields(3)=global%fields(3)
    depth=0
    if (rank==nranks-1) depth=1
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,local,temperature,MPI_COMM_WORLD,ok,depth=depth)
    call require(.not.ok,'inconsistent selected depth')
    call unchanged()
    cfl=.15_dp
    if (rank==nranks-1) cfl=.2_dp
    call compute_mpi_amr_tree_cfl_3d(tree,owners,species,local,temperature,spacing,cfl,dt,MPI_COMM_WORLD,ok)
    call require(.not.ok.and.dt==expected_dt,'inconsistent CFL controls')
    deallocate(local%fields,temperature%fields)
  end do
  if (rank==0) print *, 'PASS: owner EOS recovery, selected depth, global CFL and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(temperature%fields(peer)%values==saved%fields(peer)%values),'temperature rollback')
    end do
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
