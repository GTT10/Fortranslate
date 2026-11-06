program test_mpi_eb_interface_source_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,build_axis_plane_eb_geometry_3d
  use eb_reactive_redistribution_3d_mod, only: state_redist_register_3d, &
    initialize_state_redist_tree_register_3d,assemble_state_redist_interface_source_3d, &
    assemble_state_redist_interface_rereflux_3d,initialize_flux_rereflux_register_3d
  use mpi_eb_tree_redistribution_3d_mod, only: assemble_state_redist_mpi_interface_source_3d, &
    assemble_state_redist_mpi_interface_rereflux_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(1)
  type(eb_geometry_3d) :: geometry(2)
  type(state_redist_register_3d) :: reference_registers(2),local_registers(2)
  type(amr_tree_state_3d) :: source,reference,kept,increment,serial_increment,kept_increment
  character(len=16) :: faces(6)
  integer :: rank,nranks,ierr,mode,layout,p,owners(2),depth
  real(dp) :: expected,denominator
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(4,4,4,2,3,2,3,2,3,2,patches(1),ok)
  call require(ok,'patch')
  call initialize_amr_patch_tree_3d([4,4,4],[1],patches,tree,ok)
  call require(ok,'tree')
  call build_axis_plane_eb_geometry_3d(16,16,16,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp,-1.5_dp,2.5_dp, &
    'y',.1_dp,geometry(1),ok)
  call require(ok,'coarse geometry')
  call build_axis_plane_eb_geometry_3d(16,16,16,-.5_dp,1.5_dp,-.5_dp,1.5_dp,-.5_dp,1.5_dp, &
    'y',.1_dp,geometry(2),ok)
  call require(ok,'fine geometry')
  block
    real(dp) :: mismatch(1,4,4,4)
    type(state_redist_register_3d) :: ledger,saved
    faces='outflow';mismatch=0
    mismatch(1,2,1,2)=.25_dp ! Coarse/fine boundary cut cell.
    mismatch(1,2,2,2)=100._dp ! Covered coarse state is not interface inventory.
    mismatch(1,1,4,4)=200._dp ! A diagonal neighbor is not a face boundary.
    call initialize_flux_rereflux_register_3d(tree,1,geometry(1),mismatch,ledger,ok,faces)
    call require(ok,'flux mismatch register')
    call require(ledger%coarse_delta(1,8,7,8)==.25_dp.and.sum(ledger%coarse_delta)==.25_dp, &
      'only coarse face-boundary mismatch enters re-reflux')
    saved=ledger
    mismatch(1,2,1,2)=ieee_value(0._dp,ieee_quiet_nan)
    call initialize_flux_rereflux_register_3d(tree,1,geometry(1),mismatch,ledger,ok,faces)
    call require(.not.ok,'invalid flux mismatch rejected')
    call require(all(ledger%coarse_delta==saved%coarse_delta),'flux mismatch register rollback')
  end block
  do mode=0,1
    faces='outflow'
    if(mode==1) faces(1:2)='periodic'
    do p=1,2
      call initialize_state_redist_tree_register_3d(tree,p,1,6,reference_registers(p),ok,faces,.5_dp)
      call require(ok,'reference registers')
    enddo
    reference_registers(1)%coarse_delta(1,8,7,8)=.25_dp
    reference_registers(1)%coarse_delta(1,7,8,8)=.5_dp
    reference_registers(2)%fine_delta(1,7:8,6,7:8)=2
    reference_registers(2)%fine_delta(1,7:8,5,7:8)=100
    denominator=sum(geometry(2)%volume_fraction(7:8,5:6,7:8))
    expected=.25_dp+8._dp/denominator
    call assemble_state_redist_interface_source_3d(tree,geometry,reference_registers,1,0,reference,ok,faces)
    call require(ok,'serial reference')
    call assemble_state_redist_interface_rereflux_3d(tree,geometry,reference_registers,1,0,serial_increment,ok,faces)
    call require(ok,'serial rereflux reference')
    do layout=1,3
      owners=[0,1]
      if(layout==2) owners=[1,0]
      if(layout==3) owners=0
      block
        type(state_redist_register_3d) :: empty
        do p=1,2
          local_registers(p)=empty
          if(owners(p)==rank) local_registers(p)=reference_registers(p)
        enddo
      end block
      call assemble_state_redist_mpi_interface_source_3d(tree,owners,geometry,local_registers,1,0, &
        source,MPI_COMM_WORLD,ok,faces)
      call require(ok,'owner-local interface assembly')
      call require(allocated(source%fields(1)%values).eqv.(owners(1)==rank),'coarse owner allocation')
      call require(.not.allocated(source%fields(2)%values),'fine source absent')
      if(owners(1)==rank) then
        call require(all(source%fields(1)%values==reference%fields(1)%values),'exact full halo serial agreement')
        call require(abs(source%fields(1)%values(1,8,7,8)-expected)<1.e-14_dp,'independent kappa denominator')
      endif
      call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,geometry,local_registers,1,0, &
        increment,MPI_COMM_WORLD,ok,faces)
      call require(ok,'MPI owner-local re-reflux')
      do p=1,2
        call require(allocated(increment%fields(p)%values).eqv.(owners(p)==rank),'re-reflux owner allocation')
        if(owners(p)/=rank) cycle
        call require(all(increment%fields(p)%values==serial_increment%fields(p)%values),'exact serial/MPI re-reflux')
      enddo
      kept_increment=increment
      kept=source
      if(owners(2)==rank) local_registers(2)%fine_delta(1,7,6,7)=ieee_value(0._dp,ieee_quiet_nan)
      call assemble_state_redist_mpi_interface_source_3d(tree,owners,geometry,local_registers,1,0, &
        source,MPI_COMM_WORLD,ok,faces)
      call rejected('late fine NaN collective rollback')
      call assemble_state_redist_mpi_interface_rereflux_3d(tree,owners,geometry,local_registers,1,0, &
        increment,MPI_COMM_WORLD,ok,faces)
      call require(.not.ok,'invalid re-reflux rejected')
      do p=1,2
        if(owners(p)/=rank) cycle
        call require(all(increment%fields(p)%values==kept_increment%fields(p)%values),'re-reflux collective rollback')
      enddo
      if(owners(2)==rank) then
        local_registers(2)=reference_registers(2)
        local_registers(2)%fine_delta=0
        local_registers(2)%fine_delta(1,7,6,7)=huge(1._dp)*.75_dp
      endif
      if(owners(1)==rank) local_registers(1)%coarse_delta(1,8,7,8)=huge(1._dp)*.95_dp
      call assemble_state_redist_mpi_interface_source_3d(tree,owners,geometry,local_registers,1,0, &
        source,MPI_COMM_WORLD,ok,faces)
      call rejected('cross-owner coarse plus fine overflow rollback')
      do p=1,2
        if(owners(p)==rank) local_registers(p)=reference_registers(p)
      enddo
      depth=0
      if(rank==0) depth=1
      call assemble_state_redist_mpi_interface_source_3d(tree,owners,geometry,local_registers,1,depth, &
        source,MPI_COMM_WORLD,ok,faces)
      call rejected('inconsistent depth rollback')
    enddo
  enddo
  if(rank==0) print *, 'PASS: MPI sparse interface source, full-block kappa, periodic halos and collective rollback'
  call MPI_Finalize(ierr)
contains
  subroutine rejected(label)
    character(len=*), intent(in) :: label
    integer :: p
    call require(.not.ok,label)
    do p=1,2
      call require(allocated(source%fields(p)%values).eqv.allocated(kept%fields(p)%values),label)
      if(.not.allocated(source%fields(p)%values)) cycle
      call require(all(source%fields(p)%values==kept%fields(p)%values),label)
    enddo
  end subroutine
  subroutine require(condition,label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    integer :: status
    if(condition) return
    print *, 'FAIL rank ',rank,': '//label
    call MPI_Abort(MPI_COMM_WORLD,1,status)
  end subroutine
end program test_mpi_eb_interface_source_3d
