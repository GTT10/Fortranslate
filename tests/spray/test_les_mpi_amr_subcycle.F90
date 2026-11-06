program test_les_mpi_amr_subcycle
  use mpi_f08
  use, intrinsic :: iso_fortran_env, only: int64
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use thermo_database_mod, only: load_h2o2_elementary_thermo
  use reactive_1d_mod, only: reactive_primitive_to_conserved
  use reactive_3d_mod, only: recover_reactive_temperatures_3d
  use rotating_frame_3d_mod, only: rotating_frame_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d
  use les_amr_subcycle_mod, only: advance_amr_patch_tree_subcycle_les_3d
  use les_amr_patch_tree_mod, only: advance_amr_patch_tree_les_level_3d,convert_amr_les_tree_energy_3d
  use les_mpi_amr_level_mod, only: advance_mpi_amr_tree_les_level_3d
  use amr_patch_tree_reflux_3d_mod, only: amr_tree_flux_3d
  use les_mpi_amr_subcycle_mod, only: advance_mpi_amr_tree_subcycle_les_3d
  implicit none
  type(nasa7_species), allocatable :: species(:)
  real(dp) :: diffusivity,serial_diffusivity,wall_t(6),wall_v(3,6)
  character(len=:), allocatable :: coarse_limiter
  type(amr_patch_tree_3d) :: tree
  type(amr_patch_3d) :: patches(2)
  type(amr_tree_state_3d) :: global,reference,local,temperature,saved,oracle,oracle_t,saved_t
  type(amr_tree_state_3d) :: coarse_plain,coarse_rotated
  type(rotating_frame_3d) :: frame,trial_frame
  type(amr_tree_flux_3d), allocatable :: plain_flux(:),rotated_flux(:)
  real(dp) :: q(12),sound,exact_t,spacing(3),dt
  real(dp) :: budget(12),serial_budget(12),largest,serial_largest,inlet(12,6),inlet_t(6)
  real(dp) :: pressure(6),mach(6),pin(6)
  character(len=32) :: faces(6)
  character(len=16), parameter :: models(4)=[character(len=16) :: 'smagorinsky','wale','vreman','dynamic']
  integer :: rank,nranks,ierr,node,n(3),i,j,k,mode,owners(3)
  integer(int64) :: steps(3),serial_steps(3)
  logical :: ok
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nranks,ierr)
  call load_h2o2_elementary_thermo(species,ok)
  call require(ok,'thermo')
  wall_t=900;wall_v=0;wall_v(1,3:4)=2
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
  inlet=spread(global%fields(1)%values(:,1,1,1),2,6)
  inlet_t=reference%fields(1)%values(1,1,1,1)
  do mode=1,5
    owners=[0,1,0]
    if (mode==3) owners=0
    if (mode==4) owners=mod([0,1,2],nranks)
    faces='periodic';wall_v=0
    if (allocated(coarse_limiter)) deallocate(coarse_limiter)
    if (mod(mode,2)==0) coarse_limiter='mc'
    if (mode>=3) then
      faces='outflow';faces(3:4)='isothermal_wall';wall_v(1,3:4)=2
    end if
    if (mode==4) faces(1)='inflow'
    pressure=0;mach=-1;pin=0
    if (mode==4) then
      faces(1)='outflow';pressure(1)=102325
    end if
    if(mode==5) then
      faces(1)='inflow';mach(1)=.2_dp;coarse_limiter='mc'
    end if
    allocate(local%fields(3),temperature%fields(3))
    do node=1,3
      if (owners(node)/=rank) cycle
      local%fields(node)=global%fields(node)
      temperature%fields(node)=reference%fields(node)
    end do
    oracle=global;oracle_t=reference;dt=1.e-9_dp
    call advance_amr_patch_tree_subcycle_les_3d(tree,species,oracle,oracle_t,spacing,dt,faces, &
      .17_dp,.7_dp,.7_dp,serial_budget,serial_largest,serial_diffusivity,ok, &
      inlet,inlet_t,wall_t,wall_v,coarse_limiter,serial_steps, &
      outlet_pressure=pressure,model=trim(models(min(mode,4))),inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'serial RK2 oracle')
    call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
      .17_dp,.7_dp,.7_dp,budget,largest,diffusivity,MPI_COMM_WORLD,ok, &
      inlet,inlet_t,wall_t,wall_v,coarse_limiter,steps, &
      outlet_pressure=pressure,model=trim(models(min(mode,4))),inflow_mach=mach,inflow_pressure=pin)
    call require(ok,'distributed recursive RK2')
    call require(all(steps==[1_int64,4_int64,36_int64]).and.all(steps==serial_steps),'exact recursive step counts')
    do node=1,3
      call require(allocated(local%fields(node)%values).eqv.(owners(node)==rank),'owned state only')
      call require(allocated(temperature%fields(node)%values).eqv.(owners(node)==rank),'owned temperature only')
      if (owners(node)/=rank) cycle
      call require(all(local%fields(node)%values==oracle%fields(node)%values),'exact serial conserved state')
      call require(all(temperature%fields(node)%values==oracle_t%fields(node)%values),'exact serial temperature')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest.and.diffusivity==serial_diffusivity, &
      'exact boundary and reflux diagnostics')
    saved=local;saved_t=temperature
    dt=1._dp
    call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
      .17_dp,.7_dp,.7_dp,budget,largest,diffusivity,MPI_COMM_WORLD,ok, &
      inlet,inlet_t,wall_t,wall_v,coarse_limiter,steps, &
      outlet_pressure=pressure,model=trim(models(min(mode,4))),inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'excessive timestep rejects collectively')
    call unchanged()
    dt=1.e-9_dp
    if (rank==nranks-1) dt=2*dt
    call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
      .17_dp,.7_dp,.7_dp,budget,largest,diffusivity,MPI_COMM_WORLD,ok, &
      inlet,inlet_t,wall_t,wall_v,coarse_limiter,steps, &
      outlet_pressure=pressure,model=trim(models(min(mode,4))),inflow_mach=mach,inflow_pressure=pin)
    call require(.not.ok,'mismatched timestep rejected')
    call unchanged()
    if (mode==4) then
      dt=1e-9_dp
      if (rank==nranks-1) pressure(1)=103325
      call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
        .17_dp,.7_dp,.7_dp,budget,largest,diffusivity,MPI_COMM_WORLD,ok, &
        wall_temperature=wall_t,wall_velocity=wall_v,outlet_pressure=pressure,level_steps=steps, &
      model=trim(models(min(mode,4))),inflow_mach=mach,inflow_pressure=pin)
      call require(.not.ok,'mismatched pressure rejected')
      call unchanged()
    end if
    deallocate(local%fields,temperature%fields)
  end do
  frame%axis=3;frame%angular_speed=400._dp;frame%axis_point=[.001_dp,-.002_dp,0._dp]
  owners=[0,1,0];faces='outflow';dt=1e-9_dp
  oracle=global;oracle_t=reference
  allocate(local%fields(3),temperature%fields(3))
  do node=1,3
    if(owners(node)/=rank) cycle
    local%fields(node)=global%fields(node)
    temperature%fields(node)=reference%fields(node)
  enddo
  call convert_amr_les_tree_energy_3d(tree,oracle,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok)
  call require(ok,'serial rotating recursive LES storage')
  call convert_amr_les_tree_energy_3d(tree,local,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok,owners,rank)
  call require(ok,'MPI rotating recursive LES storage')
  call advance_amr_patch_tree_subcycle_les_3d(tree,species,oracle,oracle_t,spacing,dt,faces, &
    .17_dp,.7_dp,.7_dp,serial_budget,serial_largest,serial_diffusivity,ok, &
    frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'serial rotating recursive LES oracle')
  call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
    .17_dp,.7_dp,.7_dp,budget,largest,diffusivity,MPI_COMM_WORLD,ok, &
    frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'MPI rotating recursive LES step')
  do node=1,3
    if(owners(node)/=rank) cycle
    call require(all(local%fields(node)%values==oracle%fields(node)%values), &
      'MPI rotating recursive LES state parity')
    call require(all(temperature%fields(node)%values==oracle_t%fields(node)%values), &
      'MPI rotating recursive LES temperature parity')
  enddo
  call require(all(budget==serial_budget),'MPI rotating recursive LES budget parity')
  oracle=global;oracle_t=reference
  deallocate(local%fields,temperature%fields)
  allocate(local%fields(3),temperature%fields(3))
  do node=1,3
    if(owners(node)/=rank) cycle
    local%fields(node)=global%fields(node)
    temperature%fields(node)=reference%fields(node)
  enddo
  call convert_amr_les_tree_energy_3d(tree,local,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok,owners,rank)
  call require(ok,'MPI rotating LES level storage')
  call advance_amr_patch_tree_les_level_3d(tree,species,oracle,oracle_t,0,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,plain_flux,serial_diffusivity,ok)
  call require(ok,'MPI rotating LES level oracle')
  call advance_mpi_amr_tree_les_level_3d(tree,owners,species,local,temperature,0,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,rotated_flux,diffusivity,MPI_COMM_WORLD,ok,frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'MPI rotating LES level step')
  call convert_amr_les_tree_energy_3d(tree,local,spacing,frame,[0._dp,0._dp,0._dp],.false.,ok,owners,rank)
  call require(ok,'MPI rotating LES level recovery')
  if(rank==0) then
    call require(maxval(abs(local%fields(1)%values-oracle%fields(1)%values)/ &
      max(1._dp,abs(oracle%fields(1)%values)))<2e-12_dp,'MPI rotating LES level state parity')
    call require(maxval(abs(rotated_flux(1)%x-plain_flux(1)%x))<1e-12_dp,'MPI rotating LES level flux parity')
  endif
  oracle=global;oracle_t=reference
  do node=1,3
    if(owners(node)/=rank) cycle
    local%fields(node)=global%fields(node)
    temperature%fields(node)=reference%fields(node)
  enddo
  call convert_amr_les_tree_energy_3d(tree,local,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok,owners,rank)
  call require(ok,'MPI rotating LES fine level storage')
  allocate(coarse_plain%fields(3),coarse_rotated%fields(3))
  allocate(coarse_plain%fields(1)%values(13,4,4,4))
  coarse_plain%fields(1)%values(:12,:,:,:)=global%fields(1)%values
  coarse_plain%fields(1)%values(13,:,:,:)=reference%fields(1)%values(1,:,:,:)
  if(rank==0) coarse_rotated%fields(1)=coarse_plain%fields(1)
  call convert_amr_les_tree_energy_3d(tree,coarse_rotated,spacing,frame,[0._dp,0._dp,0._dp],.true.,ok, &
    owners,rank,maximum_depth=0)
  call require(ok,'MPI rotating LES coarse endpoint storage')
  call advance_amr_patch_tree_les_level_3d(tree,species,oracle,oracle_t,1,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,plain_flux,serial_diffusivity,ok,coarse_start=coarse_plain,coarse_end=coarse_plain, &
    fractions_start=[0._dp],fractions_end=[1._dp])
  call require(ok,'MPI rotating LES fine level oracle')
  call advance_mpi_amr_tree_les_level_3d(tree,owners,species,local,temperature,1,spacing,1e-9_dp,faces, &
    .17_dp,.7_dp,.7_dp,rotated_flux,diffusivity,MPI_COMM_WORLD,ok, &
    coarse_start=coarse_rotated,coarse_end=coarse_rotated, &
    fractions_start=[0._dp],fractions_end=[1._dp],frame=frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(ok,'MPI rotating LES fine level step')
  call convert_amr_les_tree_energy_3d(tree,local,spacing,frame,[0._dp,0._dp,0._dp],.false.,ok,owners,rank)
  call require(ok,'MPI rotating LES fine level recovery')
  do node=1,3
    if(tree%nodes(node)%level/=1.or.owners(node)/=rank) cycle
    call require(maxval(abs(local%fields(node)%values-oracle%fields(node)%values)/ &
      max(1._dp,abs(oracle%fields(node)%values)))<2e-12_dp,'MPI rotating LES fine level state parity')
    call require(maxval(abs(rotated_flux(node)%x-plain_flux(node)%x))<1e-12_dp, &
      'MPI rotating LES fine level flux parity')
  enddo
  saved=local;saved_t=temperature
  budget=-17;largest=-19;diffusivity=-29;steps=-23
  trial_frame=frame
  if(rank==nranks-1) trial_frame%angular_speed=trial_frame%angular_speed+1._dp
  call advance_mpi_amr_tree_subcycle_les_3d(tree,owners,species,local,temperature,spacing,dt,faces, &
    .17_dp,.7_dp,.7_dp,budget,largest,diffusivity,MPI_COMM_WORLD,ok,level_steps=steps, &
    frame=trial_frame,domain_lower=[0._dp,0._dp,0._dp])
  call require(.not.ok.and.all(budget==-17).and.largest==-19.and.diffusivity==-29.and.all(steps==-23), &
    'MPI rotating recursive LES frame disagreement rollback')
  do node=1,3
    if(owners(node)/=rank) cycle
    call require(all(local%fields(node)%values==saved%fields(node)%values), &
      'MPI rotating recursive LES disagreement state rollback')
    call require(all(temperature%fields(node)%values==saved_t%fields(node)%values), &
      'MPI rotating recursive LES disagreement temperature rollback')
  enddo
  if (rank==0) print *, 'PASS: MPI recursive LES: four models, 1/4/36 steps, exact state/budget, empty ranks and rollback'
  call MPI_Finalize(ierr)
contains
  subroutine unchanged()
    integer :: peer
    do peer=1,3
      if (owners(peer)/=rank) cycle
      call require(all(local%fields(peer)%values==saved%fields(peer)%values),'state rollback')
      call require(all(temperature%fields(peer)%values==saved_t%fields(peer)%values),'temperature rollback')
    end do
    call require(all(budget==serial_budget).and.largest==serial_largest.and.diffusivity==serial_diffusivity,'diagnostic rollback')
    call require(all(steps==serial_steps),'counter rollback')
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
