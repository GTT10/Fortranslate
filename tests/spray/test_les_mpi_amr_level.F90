program test_les_mpi_amr_level
  use mpi_f08
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use amr_patch_tree_hydro_3d_mod, only: compute_amr_patch_tree_cfl_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_level_3d
  use les_mpi_amr_level_mod, only: advance_mpi_amr_tree_les_level_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: theta,serial_theta,saved_theta,wall_t(6),wall_v(3,6)
  character(len=:), allocatable :: coarse_limiter
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  type(amr_tree_state_3d) :: coarse_start,coarse_end,start_local,end_local
  type(amr_tree_flux_3d), allocatable :: flux(:),oracle_flux(:),saved_flux(:)
  real(dp) :: start_fraction(2),end_fraction(2),wrong_fraction(2)
  real(dp) :: q(12),sound,exact_t,spacing(3),dt,expected_dt,cfl
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  character(len=32) :: faces(6),reconstruction,model
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3),depth,m
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
  do m=1,4
    select case(m)
    case(1)
      model="smagorinsky"
    case(2)
      model="wale"
    case(3)
      model="vreman"
    case(4)
      model="dynamic"
    end select
  do mode=1,4
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';reconstruction='pcm'
    wall_t=1100;wall_v=0
    if (allocated(coarse_limiter)) deallocate(coarse_limiter)
    if (mod(mode,2)==0) coarse_limiter='mc'
    if (mod(mode,2)==0) reconstruction='characteristic_plm'
    if (mode>=3) then
      faces='outflow';faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
    end if
    if (mode==4) faces(1)='inflow'
    do depth=0,2
      allocate(local%fields(3),temperature%fields(3),coarse_start%fields(3),coarse_end%fields(3), &
        start_local%fields(3),end_local%fields(3))
      do node=1,3
        n=tree%nodes(node)%shape
        if (tree%nodes(node)%level<depth) then
          allocate(coarse_start%fields(node)%values(13,n(1),n(2),n(3)))
          coarse_start%fields(node)%values(:12,:,:,:)=global%fields(node)%values
          coarse_start%fields(node)%values(13,:,:,:)=reference%fields(node)%values(1,:,:,:)
          coarse_end%fields(node)=coarse_start%fields(node)
          coarse_end%fields(node)%values(:12,:,:,:)=1.02_dp*coarse_end%fields(node)%values(:12,:,:,:)
          if (owners(node)==rank) then
            start_local%fields(node)=coarse_start%fields(node)
            end_local%fields(node)=coarse_end%fields(node)
          end if
        end if
        if (owners(node)/=rank) cycle
        local%fields(node)=global%fields(node)
        temperature%fields(node)=reference%fields(node)
      end do
      start_fraction=[.1_dp,.25_dp];end_fraction=[.6_dp,.75_dp]
      oracle=global;oracle_t=reference;dt=1e-9_dp
      call advance_amr_patch_tree_les_level_3d(tree,species,oracle,oracle_t,depth,spacing,dt,faces, &
        .17_dp,.7_dp,.7_dp,oracle_flux,serial_theta,ok,coarse_start,coarse_end, &
        start_fraction(:depth),end_fraction(:depth), &
        inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=wall_v,coarse_limiter=coarse_limiter,model=model)
      call require(ok,'serial level RK2')
      call advance_mpi_amr_tree_les_level_3d(tree,owners,species,local,temperature,depth,spacing,dt,faces, &
        .17_dp,.7_dp,.7_dp,flux,theta,MPI_COMM_WORLD,ok,start_local,end_local,start_fraction(:depth), &
        end_fraction(:depth), &
        inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=wall_v,coarse_limiter=coarse_limiter,model=model)
      call require(ok,'MPI level RK2')
      call require(theta==serial_theta,'exact diffusivity diagnostic')
      do node=1,3
        call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owned state only')
        call require(allocated(flux(node)%x).eqv.(owners(node)==rank.and.tree%nodes(node)%level==depth),'active flux only')
        if (owners(node)/=rank) cycle
        call require(all(local%fields(node)%values==oracle%fields(node)%values),'exact serial state')
        call require(all(temperature%fields(node)%values==oracle_t%fields(node)%values),'exact serial temperature')
        if (tree%nodes(node)%level/=depth) then
          call require(all(local%fields(node)%values==global%fields(node)%values),'inactive state unchanged')
          call require(all(temperature%fields(node)%values==reference%fields(node)%values),'inactive temperature unchanged')
          cycle
        end if
        call require(all(flux(node)%x==oracle_flux(node)%x),'exact x integrated flux')
        call require(all(flux(node)%y==oracle_flux(node)%y),'exact y integrated flux')
        call require(all(flux(node)%z==oracle_flux(node)%z),'exact z integrated flux')
      end do
      saved=local;saved_t=temperature;saved_flux=flux;saved_theta=theta
      call advance_mpi_amr_tree_les_level_3d(tree,owners,species,local,temperature,depth,spacing,1._dp,faces, &
        .17_dp,.7_dp,.7_dp,flux,theta,MPI_COMM_WORLD,ok,start_local,end_local,start_fraction(:depth), &
        end_fraction(:depth), &
        inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=wall_v,coarse_limiter=coarse_limiter,model=model)
      call require(.not.ok,'excessive interval rejected')
      call unchanged()
      if (depth>0) then
        wrong_fraction=end_fraction
        if (rank==nranks-1) wrong_fraction(1)=.7_dp
        call advance_mpi_amr_tree_les_level_3d(tree,owners,species,local,temperature,depth,spacing,dt,faces, &
          .17_dp,.7_dp,.7_dp,flux,theta,MPI_COMM_WORLD,ok,start_local,end_local, &
          start_fraction(:depth),wrong_fraction(:depth),inflow_state=inlet,inflow_temperature=inlet_t, &
        wall_temperature=wall_t,wall_velocity=wall_v,coarse_limiter=coarse_limiter,model=model)
        call require(.not.ok,'second-stage endpoint disagreement rejected')
        call unchanged()
      end if
      deallocate(local%fields,temperature%fields,coarse_start%fields,coarse_end%fields,start_local%fields,end_local%fields)
    end do
  end do
  end do
  if (rank==0) print *, 'PASS: MPI level LES, temporal donors, exact state/T/flux, inactive levels and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged()
    integer :: peer
    call require(theta==saved_theta,"diffusivity rollback")
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
      call require(allocated(flux(peer)%x).eqv.allocated(saved_flux(peer)%x),'flux allocation rollback')
      if (.not.allocated(saved_flux(peer)%x)) cycle
      call require(all(flux(peer)%x==saved_flux(peer)%x),'x flux rollback')
      call require(all(flux(peer)%y==saved_flux(peer)%y),'y flux rollback')
      call require(all(flux(peer)%z==saved_flux(peer)%z),'z flux rollback')
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
