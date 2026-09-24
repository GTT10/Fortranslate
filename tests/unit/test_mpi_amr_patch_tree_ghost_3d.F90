program test_mpi_amr_patch_tree_ghost_3d
  use mpi_f08
  use, intrinsic :: ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_ghost_3d_mod, only: fill_amr_patch_tree_ghosts_3d
  use mpi_amr_patch_tree_ghost_3d_mod, only: fill_mpi_amr_tree_ghosts_3d
  implicit none
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(3)
  type(amr_tree_state_3d) :: global,local,ghosts,saved,future,local_future
  real(dp), allocatable :: oracle(:,:,:,:)
  real(dp) :: inlet(7,6),velocity(3,6),rho,fractions(2),bad_fractions(2)
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,scheme,owners(4),active,blend,levels
  character(len=32) :: faces(6)
  character(len=:), allocatable :: limiter
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call require(nranks>=2,'rank count')
  call initialize_amr_patch_3d(4,4,4,1,2,1,2,1,2,2,patches(1),ok)
  call require(ok,'child')
  call initialize_amr_patch_3d(4,4,4,3,4,1,2,1,2,2,patches(2),ok)
  call require(ok,'peer')
  call initialize_amr_patch_3d(4,4,4,1,1,1,1,1,1,3,patches(3),ok)
  call require(ok,'grandchild')
  call initialize_amr_patch_tree_3d([4,4,4],[1,1,2],patches,tree,ok)
  call require(ok,'tree')
  allocate(global%fields(4))
  do node=1,4
    n=tree%nodes(node)%shape
    allocate(global%fields(node)%values(7,n(1),n(2),n(3)))
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          rho=2+.01_dp*(i+2*j+3*k+node)
          global%fields(node)%values(:,i,j,k)=[rho,.2_dp*i,.3_dp*j,.4_dp*k,1000.0_dp+node,.3_dp*rho,.7_dp*rho]
        end do
      end do
    end do
  end do
  future=global
  do node=1,4
    if (tree%nodes(node)%level==2) then
      deallocate(future%fields(node)%values)
    else
      future%fields(node)%values=1.25_dp*future%fields(node)%values
    end if
  end do
  inlet=spread([3.0_dp,4.0_dp,5.0_dp,6.0_dp,2000.0_dp,1.0_dp,2.0_dp],2,6)
  do mode=1,6
    owners=[0,1,0,1]
    if (mode==4.or.mode==6) owners=0
    allocate(local%fields(4),local_future%fields(4))
    do node=1,4
      if (owners(node)==rank) then
        local%fields(node)=global%fields(node)
        local_future%fields(node)=future%fields(node)
      end if
    end do
    faces='periodic';velocity=0
    if (mode>=2) faces=[character(len=32)::'slip_wall','outflow','inflow','no_slip_wall','isothermal_wall','inflow']
    if (mode>=3) then
      velocity(1,4)=.2_dp;velocity(2,5)=.3_dp
    end if
    if(mode>=5) then
      faces='periodic';faces(1:2)='no_slip_wall';velocity=0
      velocity(2,1)=1;velocity(2,2)=3
    end if
    do scheme=1,3
      if (allocated(limiter)) deallocate(limiter)
      if (scheme==2) limiter='mc'
      if (scheme==3) limiter='minmod'
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
        viscous_wall=mode>=3,inflow_state=inlet,wall_velocity=velocity,coarse_limiter=limiter,species_count=2)
      call require(ok,'collective halo fill')
      do node=1,4
        call require(allocated(ghosts%fields(node)%values).eqv.(owners(node)==rank),'only owned halos')
        if (owners(node)/=rank) cycle
        call fill_amr_patch_tree_ghosts_3d(tree,global,node,faces,oracle,ok,viscous_wall=mode>=3, &
          inflow_state=inlet,wall_velocity=velocity,coarse_limiter=limiter,species_count=2)
        call require(ok,'serial oracle')
        call require(all(shape(oracle)==shape(ghosts%fields(node)%values)),'halo shape')
        call require(all(oracle==ghosts%fields(node)%values),'exact serial halo parity')
      end do
      block
        call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
          viscous_wall=mode>=3,inflow_state=inlet,wall_velocity=velocity, &
          coarse_limiter=limiter,species_count=2,ghost_width=7)
        call require(ok,'wide collective halo fill')
        do node=1,4
          call require(allocated(ghosts%fields(node)%values).eqv.(owners(node)==rank),'only owned wide halos')
          if(owners(node)/=rank) cycle
          call fill_amr_patch_tree_ghosts_3d(tree,global,node,faces,oracle,ok, &
            viscous_wall=mode>=3,inflow_state=inlet,wall_velocity=velocity, &
            coarse_limiter=limiter,species_count=2,ghost_width=7)
          call require(ok.and.all(oracle==ghosts%fields(node)%values),'wide exact serial halo parity')
        end do
        saved=ghosts
        call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
          viscous_wall=mode>=3,inflow_state=inlet,wall_velocity=velocity, &
          coarse_limiter=limiter,species_count=2,ghost_width=7+rank)
        call require(.not.ok,'rank-disagreeing halo width rejected')
        do node=1,4
          if(owners(node)/=rank) cycle
          call require(all(ghosts%fields(node)%values==saved%fields(node)%values),'wide halo mismatch rollback')
        end do
      end block
    end do
    do scheme=1,3
      if (allocated(limiter)) deallocate(limiter)
      if (scheme==2) limiter='mc'
      if (scheme==3) limiter='minmod'
    do active=-1,2
      levels=2
      if (active>=0) levels=active
      do blend=0,2
        fractions=real(blend,dp)*.5_dp
        if (blend==1) fractions=[.25_dp,.75_dp]
        call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
          viscous_wall=mode>=3,inflow_state=inlet,wall_velocity=velocity,coarse_limiter=limiter,species_count=2, &
          coarse_future=local_future,coarse_fractions=fractions(:levels),depth=active)
        call require(ok,'collective temporal halo fill')
        do node=1,4
          call require(allocated(ghosts%fields(node)%values).eqv. &
            (owners(node)==rank.and.(active<0.or.tree%nodes(node)%level==active)),'only active owned temporal halos')
          if (.not.allocated(ghosts%fields(node)%values)) cycle
          call fill_amr_patch_tree_ghosts_3d(tree,global,node,faces,oracle,ok,viscous_wall=mode>=3, &
            inflow_state=inlet,wall_velocity=velocity,coarse_limiter=limiter,species_count=2, &
            coarse_future=future,coarse_fractions=fractions(:tree%nodes(node)%level))
          call require(ok,'serial temporal oracle')
          call require(all(oracle==ghosts%fields(node)%values),'exact temporal serial halo parity')
        end do
      end do
    end do
    end do
    saved=ghosts
    bad_fractions=fractions
    if (rank==nranks-1) bad_fractions(1)=.5_dp
    call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
      inflow_state=inlet,coarse_future=local_future,coarse_fractions=bad_fractions,depth=2)
    call require(.not.ok,'different valid fractions rejected')
    call unchanged()
    if (rank==0) then
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
        inflow_state=inlet,coarse_future=local_future,coarse_fractions=fractions,depth=2)
    else
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok,inflow_state=inlet,depth=2)
    end if
    call require(.not.ok,'different endpoint presence rejected')
    call unchanged()
    if (rank==0) deallocate(local_future%fields(1)%values)
    call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
      inflow_state=inlet,coarse_future=local_future,coarse_fractions=fractions,depth=2)
    call require(.not.ok,'missing coarse endpoint rejected')
    call unchanged()
    if (rank==0) local_future%fields(1)=future%fields(1)
    saved=ghosts
    if (rank==0) local%fields(1)%values(1,1,1,1)=ieee_value(0.0_dp,ieee_quiet_nan)
    call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
      inflow_state=inlet,coarse_limiter=limiter)
    call require(.not.ok,'nonfinite donor rejected collectively')
    call unchanged()
    if (rank==0) local%fields(1)%values=global%fields(1)%values
    if (mode>=3) then
      if (rank==0) local%fields(1)%values(1,1,1,1)=-1.0_dp
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok, &
        viscous_wall=.true.,inflow_state=inlet,wall_velocity=velocity)
      call require(.not.ok,'wall transform failure after exchange rejected collectively')
      call unchanged()
      if (rank==0) local%fields(1)%values=global%fields(1)%values
    end if
    if (rank==nranks-1) faces(1)='invalid' 
    call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,7,faces,ghosts,MPI_COMM_WORLD,ok,inflow_state=inlet)
    call require(.not.ok,'invalid boundary on empty rank rejected')
    call unchanged()
    deallocate(local%fields,local_future%fields)
  end do
  call mach_checks()
  if (rank==0) print *, 'PASS: MPI PCM/PLM halos match serial, physical walls, empty owners and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine mach_checks()
    type(nasa7_species), allocatable :: thermo(:)
    real(dp) :: q(12),sound,ma(6),tin(6),pin(6),bad_ma(6),bad_t(6),bad_p(6)
    integer :: face,check
    call load_h2o2_elementary_thermo(thermo,ok)
    call require(ok,'Mach MPI thermo')
    do node=1,4
      deallocate(global%fields(node)%values)
      n=tree%nodes(node)%shape
      allocate(global%fields(node)%values(13,n(1),n(2),n(3)))
      do k=1,n(3)
        do j=1,n(2)
          do i=1,n(1)
            q=0;q(1)=.3_dp+.01_dp*i;q(5)=101325+100*j
            q(6)=.03_dp;q(9)=.23_dp;q(12)=.74_dp;q(2:4)=[real(i,dp),real(j,dp),real(k,dp)]
            call reactive_primitive_to_conserved(thermo,q,global%fields(node)%values(:12,i,j,k), &
              global%fields(node)%values(13,i,j,k),sound,ok)
            call require(ok,'Mach MPI initial gas')
          end do
        end do
      end do
    end do
    future=global
    do node=1,4
      if(tree%nodes(node)%level==2) then
        deallocate(future%fields(node)%values)
      else
        future%fields(node)%values(:12,:,:,:)=1.2_dp*future%fields(node)%values(:12,:,:,:)
      end if
    end do
    allocate(local%fields(4),local_future%fields(4))
    owners=[0,1,0,1];fractions=[.25_dp,.75_dp]
    do node=1,4
      if(owners(node)/=rank) cycle
      local%fields(node)=global%fields(node);local_future%fields(node)=future%fields(node)
    end do
    do face=1,6
      faces='outflow';faces(face)='inflow';ma=-1;ma(face)=.2_dp;tin=700;pin=0
      if(modulo(face,2)==0) pin(face)=120000
      call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,13,faces,ghosts,MPI_COMM_WORLD,ok, &
        species=thermo,inflow_mach=ma,inflow_temperature=tin,inflow_pressure=pin,ghost_width=7, &
        coarse_limiter='mc',species_count=7,coarse_future=local_future,coarse_fractions=fractions)
      call require(ok,'Mach MPI interpolated wide halo')
      do node=1,4
        call require(allocated(ghosts%fields(node)%values).eqv.(owners(node)==rank),'Mach owned halos')
        if(owners(node)/=rank) cycle
        call fill_amr_patch_tree_ghosts_3d(tree,global,node,faces,oracle,ok, &
          species=thermo,inflow_mach=ma,inflow_temperature=tin,inflow_pressure=pin,ghost_width=7, &
          coarse_limiter='mc',species_count=7,coarse_future=future, &
          coarse_fractions=fractions(:tree%nodes(node)%level))
        call require(ok,'Mach serial interpolated halo')
        call require(all(oracle==ghosts%fields(node)%values),'Mach exact distributed halo parity')
      end do
      saved=ghosts
      do check=1,4
        bad_ma=ma;bad_t=tin;bad_p=pin
        if(rank==nranks-1) then
          if(check==1) bad_ma(face)=.3_dp
          if(check==2) bad_t(face)=701
          if(check==3) bad_p(face)=130000
          if(check==4) bad_t(face)=-1
        end if
        call fill_mpi_amr_tree_ghosts_3d(tree,owners,local,13,faces,ghosts,MPI_COMM_WORLD,ok, &
          species=thermo,inflow_mach=bad_ma,inflow_temperature=bad_t,inflow_pressure=bad_p,ghost_width=7)
        call require(.not.ok,'Mach rank mismatch or invalid input rejected')
        call unchanged()
      end do
    end do
  end subroutine mach_checks

  subroutine unchanged()
    integer :: peer
    do peer=1,4
      call require(allocated(ghosts%fields(peer)%values).eqv.allocated(saved%fields(peer)%values),'rollback allocation')
      if (.not.allocated(saved%fields(peer)%values)) cycle
      call require(all(ghosts%fields(peer)%values==saved%fields(peer)%values),'rollback values')
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
