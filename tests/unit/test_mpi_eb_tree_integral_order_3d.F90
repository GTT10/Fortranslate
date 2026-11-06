program test_mpi_eb_tree_integral_order_3d
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use h2o2_full_thermo_mod, only: load_h2o2_full_thermo
  use mixture_thermo_mod, only: mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_level_set_eb_geometry_3d
  use eb_reactive_patch_tree_3d_mod, only: eb_reactive_tree_3d,initialize_reactive_eb_tree_3d, &
    integrate_reactive_eb_tree_3d
  use mpi_eb_tree_transfer_3d_mod, only: integrate_mpi_eb_tree_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  type(amr_patch_3d) :: patches(1)
  type(amr_patch_tree_3d) :: topology
  type(eb_geometry_3d) :: finest,geometry(2)
  type(eb_reactive_tree_3d) :: tree
  type(amr_tree_state_3d) :: local
  real(dp) :: q(15),cell(15),root(15,4,4,4),root_t(4,4,4),phi(0:8,0:8,0:8),t,sound,r
  real(dp) :: serial(15),parallel(15),oracle(15),reference(15),large,tiny,volume
  real(dp), volatile :: coarse_total,fine_total,flat_total
  integer :: rank,nproc,ierr,p,layout,mode,owners(2)
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  call require(nproc==3,'three ranks include an idle rank')
  call load_h2o2_full_thermo(species,ok);call require(ok,'thermo')
  q=0;q(1)=1;q(15)=1
  r=mixture_specific_gas_constant(species,q(6:),ok);call require(ok,'gas constant')
  q(5)=r*1000
  call reactive_primitive_to_conserved(species,q,cell,t,sound,ok);call require(ok,'initial gas')
  do p=1,15
    root(p,:,:,:)=cell(p)
  enddo
  root_t=t;phi=1
  call build_level_set_eb_geometry_3d(phi,[0._dp,0._dp,0._dp],[4._dp,4._dp,4._dp],finest,ok)
  call require(ok,'geometry')
  call initialize_amr_patch_3d(4,4,4,1,2,1,4,1,4,2,patches(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patches,topology,ok)
  call require(ok,'tree')
  call initialize_reactive_eb_tree_3d(species,topology,finest,root,root_t,tree,ok)
  call require(ok,'initial fields')
  do p=1,2
    geometry(p)=tree%nodes(p)%geometry
    tree%nodes(p)%state=0
  enddo
  volume=geometry(1)%dx*geometry(1)%dy*geometry(1)%dz
  call require(volume==1._dp,'exact root cell volume')
  volume=geometry(2)%dx*geometry(2)%dy*geometry(2)%dz
  call require(volume==.125_dp,'exact fine cell volume')
  ! These signed fields represent source/chemical ledgers, not physical states.
  ! A covered coarse sentinel must not contribute to the composite integral.
  tree%nodes(1)%state(:,1:2,:,:)=scale(1._dp,60)
  large=scale(1._dp,54);tiny=scale(1._dp,-80)
  do mode=1,2
    if(mode==2) tiny=-tiny
    tree%nodes(1)%state(13,3,1,1)=large*tiny
    tree%nodes(2)%state(13,1,1,1)=-8*large*tiny
    tree%nodes(2)%state(13,2,1,1)=8*tiny
    tree%nodes(1)%state(5,3,1,1)=2
    tree%nodes(2)%state(5,1,1,1)=32
    tree%nodes(2)%state(5,2,1,1)=64
    ! Independent three-term floating-point oracle: fixed per-node grouping.
    ! VOLATILE prevents constant folding/reassociation of the intended order.
    coarse_total=large*tiny
    fine_total=-large*tiny
    fine_total=fine_total+tiny
    oracle=0;oracle(13)=coarse_total+fine_total;oracle(5)=14
    flat_total=coarse_total
    flat_total=flat_total-large*tiny
    flat_total=flat_total+tiny
    call require(flat_total/=oracle(13),'cancellation control distinguishes old global-cell grouping')
    do layout=1,2
      owners=[0,2]
      if(layout==2) owners=[2,0]
      if(allocated(local%fields)) deallocate(local%fields)
      allocate(local%fields(2))
      do p=1,2
        if(owners(p)==rank) local%fields(p)%values=tree%nodes(p)%state
      enddo
      call require(rank/=1.or.(.not.allocated(local%fields(1)%values).and. &
        .not.allocated(local%fields(2)%values)),'idle rank owns no field arrays')
      serial=-17;parallel=-19
      call integrate_reactive_eb_tree_3d(tree,serial,ok)
      call require(ok,'serial source integral')
      call integrate_mpi_eb_tree_3d(topology,owners,geometry,local,parallel,MPI_COMM_WORLD,ok)
      call require(ok,'MPI source integral')
      if(rank==0) print *, 'ORDER',mode,layout,serial(13),parallel(13),oracle(13),flat_total
      call require(all(parallel==oracle),'MPI agrees with independent canonical grouping and masking')
      call require(all(serial==oracle),'serial agrees with independent canonical grouping and masking')
      call require(all(serial==parallel),'bit-exact serial/MPI integral, no tolerance relaxation')
      if(layout==1) reference=parallel
      call require(all(parallel==reference),'swapped owner mapping does not change integral')
    enddo
  enddo
  if(rank==0) print *, 'PASS: canonical node-ordered cancellation integral with swapped owners and idle rank'
  call MPI_Finalize(ierr)
contains
  subroutine require(condition,message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    logical :: agreed
    integer :: status
    call MPI_Allreduce(condition,agreed,1,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,status)
    if(status/=MPI_SUCCESS.or..not.agreed) then
      print *, 'FAIL:',rank,trim(message)
      call MPI_Abort(MPI_COMM_WORLD,1,status)
    endif
  end subroutine require
end program test_mpi_eb_tree_integral_order_3d
