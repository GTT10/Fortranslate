module spray_mpi_amr_application_mod
  use mpi_f08
  use precision_mod, only: dp
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use simulation_config_reactive_3d_mod, only: reactive_3d_config,resolved_reactive_boundaries_3d
  use reactive_entropy_wave_3d_problem_mod, only: initialize_reactive_problem_3d
  use amr_hierarchy_3d_mod, only: amr_patch_3d,initialize_amr_patch_3d
  use amr_patch_tree_3d_mod, only: amr_patch_tree_3d,amr_tree_state_3d,initialize_amr_patch_tree_3d, &
    amr_tree_integer_layout_3d
  use mpi_amr_patch_tree_thermo_3d_mod, only: compute_mpi_amr_tree_cfl_3d,recover_mpi_amr_tree_temperatures_3d
  use mpi_amr_patch_tree_transport_3d_mod, only: compute_mpi_amr_tree_transport_timestep_3d
  use mpi_amr_patch_tree_restriction_3d_mod, only: average_down_mpi_amr_tree_3d
  use mpi_amr_patch_tree_csv_3d_mod, only: write_mpi_amr_tree_csv_3d
  use spray_parcel_mod
  use spray_mpi_mod, only: inject_mpi_cone_parcels
  use spray_mpi_amr_tree_mod, only: migrate_amr_tree_spray_parcels
  use spray_mpi_amr_reactive_mod, only: advance_mpi_amr_spray_reactive
  use spray_mpi_amr_remap_mod, only: remap_mpi_amr_spray
  implicit none
  private
  public :: run_mpi_amr_spray_application
contains
  subroutine run_mpi_amr_spray_application(input_path,topology_path,spray_path,config,species,reactions, &
      transport,mole_fractions,comm,chemistry_integrator)
    character(len=*), intent(in) :: input_path,topology_path,spray_path
    character(len=*), optional, intent(in) :: chemistry_integrator
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), intent(in) :: mole_fractions(:)
    type(MPI_Comm), intent(in) :: comm
    type(reactive_3d_config) :: local
    type(amr_patch_tree_3d) :: tree,target
    type(amr_tree_state_3d) :: state,temperature
    type(spray_parcel), allocatable :: parcels(:)
    type(spray_liquid), allocatable :: liquids(:)
    type(spray_options) :: options
    real(dp), allocatable :: masses(:,:),x(:),y(:),z(:),mf(:),ledger(:),global_ledger(:)
    integer, allocatable :: owners(:),assigned(:),ratios(:)
    integer(int64), allocatable :: origin(:,:),domain(:,:)
    integer(int64) :: last_id
    real(dp) :: lower(3),spacing(3),h(3),density,time,dt,transport_dt,diffusivity,correction,theta,stop_time
    integer :: rank,nranks,ierr,unit,ios,node,n(3),i,p,nv,step,attempt,particle_steps,number
    integer :: component_count,parcel_count,maximum_parcels
    character(len=24) :: vapor_names(16)
    real(dp) :: fractions(16),liquid_density(16),liquid_cp(16),reference_temperature(16), &
      reference_pressure(16),latent_heat(16)
    real(dp) :: total_mass,center(3),axis(3),speed,half_angle,diameter,drop_temperature,regrid_time,local_mass
    logical :: ok,subcycling,regridded
    character(len=32) :: faces(6)
    character(len=1024) :: regrid_topology,particle_path
    character(len=24) :: suffix
    character(len=:), allocatable :: header
    namelist /amr_spray/ component_count,vapor_names,fractions,liquid_density,liquid_cp,reference_temperature, &
      reference_pressure,latent_heat,parcel_count,total_mass,center,axis,speed,half_angle,diameter,drop_temperature, &
      maximum_parcels,options,subcycling,regrid_time,regrid_topology
    call MPI_Comm_rank(comm,rank,ierr)
    call MPI_Comm_size(comm,nranks,ierr)
    faces=resolved_reactive_boundaries_3d(config)
    if (any(faces/='periodic')) call fatal('AMR spray requires periodic gas boundaries')
    ! Inputs are shared files, as for the existing MPI tree applications.
    component_count=1;vapor_names='';vapor_names(1)='H2O';fractions=0;fractions(1)=1
    liquid_density=1000;liquid_cp=4184;reference_temperature=373.15_dp
    reference_pressure=101325;latent_heat=2.257e6_dp
    parcel_count=1;total_mass=1e-10_dp;center=0;axis=[1._dp,0._dp,0._dp]
    speed=10;half_angle=0;diameter=1e-5_dp;drop_temperature=350;maximum_parcels=100000
    subcycling=.false.;regrid_time=-1;regrid_topology='';regridded=.false.
    open(newunit=unit,file=spray_path,status='old',action='read',iostat=ios)
    if (ios/=0) call fatal('Cannot open AMR spray input')
    read(unit,nml=amr_spray,iostat=ios)
    close(unit)
    if (ios/=0) call fatal('Invalid amr_spray namelist')
    if (component_count<1.or.component_count>16) call fatal('Invalid liquid component count')
    if (.not.ieee_is_finite(regrid_time)) call fatal('Invalid regrid time')
    if ((regrid_time>=0).neqv.(len_trim(regrid_topology)>0)) call fatal('Regrid needs time and topology')
    call check_output(trim(config%output_file))
    write(suffix,'(i0)') rank
    call check_output(trim(config%output_file)//'.particles.'//trim(suffix)//'.csv')
    options%maximum_parcels=maximum_parcels
    if (.not.valid_spray_options(options)) call fatal('Invalid spray options')
    allocate(liquids(component_count),parcels(0),masses(component_count,0))
    do i=1,component_count
      call initialize_spray_liquid(species,trim(vapor_names(i)),liquid_density(i),liquid_cp(i), &
        reference_temperature(i),reference_pressure(i),latent_heat(i),liquids(i),ok)
      if (.not.ok) call fatal('Invalid liquid properties/species')
    end do
    call read_topology(topology_path,tree)
    if (regrid_time>=0) call read_topology(trim(regrid_topology),target)
    owners=[(modulo(node-1,nranks),node=1,size(tree%nodes))]
    call amr_tree_integer_layout_3d(tree,origin,domain,ratios,ok)
    if (.not.ok) call fatal('Invalid tree coordinate layout')
    lower=[config%x_lower,config%y_lower,config%z_lower]
    spacing=([config%x_upper,config%y_upper,config%z_upper]-lower)/real(tree%nodes(1)%shape,dp)
    nv=reactive_nvar(size(species))
    allocate(state%fields(size(tree%nodes)),temperature%fields(size(tree%nodes)),mf(size(species)))
    allocate(ledger(nv),global_ledger(nv))
    do node=1,size(tree%nodes)
      if (owners(node)/=rank) cycle
      n=tree%nodes(node)%shape
      h=spacing*real(domain(:,1),dp)/real(domain(:,tree%nodes(node)%level+1),dp)
      x=[(lower(1)+(real(origin(1,node),dp)+i-.5_dp)*h(1),i=1,n(1))]
      y=[(lower(2)+(real(origin(2,node),dp)+i-.5_dp)*h(2),i=1,n(2))]
      z=[(lower(3)+(real(origin(3,node),dp)+i-.5_dp)*h(3),i=1,n(3))]
      allocate(state%fields(node)%values(nv,n(1),n(2),n(3)),temperature%fields(node)%values(1,n(1),n(2),n(3)))
      local=config;local%nx=n(1);local%ny=n(2);local%nz=n(3)
      call initialize_reactive_problem_3d(species,local,x,y,z,state%fields(node)%values, &
        temperature%fields(node)%values(1,:,:,:),density,mf,ok,mole_fractions)
      if (.not.ok) call fatal('Gas initialization failed')
    end do
    call average_down_mpi_amr_tree_3d(tree,owners,state,nv,comm,ok)
    if (.not.ok) call fatal('Initial parent synchronization failed')
    call recover_mpi_amr_tree_temperatures_3d(tree,owners,species,state,temperature,comm,ok)
    if (.not.ok) call fatal('Initial gas recovery failed')
    number=0;local_mass=0;last_id=0
    if (rank==0) then
      number=parcel_count;local_mass=total_mass
    end if
    call inject_mpi_cone_parcels(comm,liquids(1),local_mass,number,center,axis,speed,half_angle,diameter,drop_temperature, &
      maximum_parcels,last_id,parcels,ledger,ok,liquids,fractions(:component_count),masses)
    if (.not.ok) call fatal('Initial injection failed')
    call MPI_Allreduce(ledger,global_ledger,nv,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    if (ierr/=MPI_SUCCESS) call fatal('Injection ledger reduction failed')
    call migrate_amr_tree_spray_parcels(comm,tree,owners,parcels,lower,spacing,[.true.,.true.,.true.], &
      maximum_parcels,ok,masses)
    if (.not.ok) call fatal('Initial particle migration failed')
    time=0;step=0
    do
      if (.not.regridded.and.regrid_time>=0) then
        if (time==regrid_time) then
          assigned=[(modulo(node-1,nranks),node=1,size(target%nodes))]
          call remap_mpi_amr_spray(comm,species,tree,owners,state,temperature,parcels,masses,target,assigned, &
            lower,spacing,[.true.,.true.,.true.],maximum_parcels,ok)
          if (.not.ok) call fatal('Coupled regrid failed')
          regridded=.true.
        end if
      end if
      if (time>=config%final_time) exit
      if (step>=config%maximum_steps) call fatal('Maximum steps reached')
      stop_time=config%final_time
      if (.not.regridded.and.regrid_time>time) stop_time=min(stop_time,regrid_time)
      ! Finest-level bounds also safely bound recursive gas stepping.
      call compute_mpi_amr_tree_cfl_3d(tree,owners,species,state,temperature,spacing,config%cfl,dt,comm,ok)
      if (.not.ok) call fatal('Gas CFL calculation failed')
      dt=min(dt,stop_time-time)
      if (config%transport_enabled) then
        call compute_mpi_amr_tree_transport_timestep_3d(tree,owners,species,transport,state,temperature,spacing, &
          config%transport_cfl,config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          transport_dt,diffusivity,comm,ok)
        if (.not.ok) call fatal('Transport timestep failed')
        dt=min(dt,transport_dt)
      end if
      do attempt=1,20
        if (time+dt<=time) call fatal('Timestep cannot advance time')
        correction=0;theta=1
        call advance_mpi_amr_spray_reactive(comm,tree,owners,species,reactions,transport,liquids,options, &
          state,temperature,parcels,masses,lower,spacing,dt,config%reconstruction,config%limiter,config%riemann_solver, &
          config%chemistry_enabled,config%transport_enabled.and.config%viscosity_enabled, &
          config%transport_enabled.and.config%thermal_conduction_enabled, &
          config%transport_enabled.and.config%species_diffusion_enabled, &
          config%transport_enabled.and.config%barodiffusion_enabled, &
          config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance,particle_steps,correction,theta, &
          ok,last_id,chemistry_integrator,subcycling)
        if (ok) exit
        dt=.5_dp*dt
      end do
      if (.not.ok) call fatal('Coupled step failed after retries')
      if (dt==stop_time-time) then
        time=stop_time
      else
        time=time+dt
      end if
      step=step+1
    end do
    call write_mpi_amr_tree_csv_3d(trim(config%output_file),tree,owners,species,state,temperature,lower,spacing,time,comm,ok)
    if (.not.ok) call fatal('Gas CSV output failed')
    write(suffix,'(i0)') rank
    particle_path=trim(config%output_file)//'.particles.'//trim(suffix)//'.csv'
    open(newunit=unit,file=trim(particle_path),status='replace',action='write',iostat=ios)
    if (ios/=0) call fatal('Cannot open particle CSV')
    header='id,parent_id,material_id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do i=1,component_count
      write(suffix,'(i0)') i
      header=header//',component_mass_'//trim(suffix)
    end do
    write(unit,'(a)',iostat=ios) header
    if (ios/=0) call fatal('Cannot write particle CSV header')
    do p=1,size(parcels)
      write(unit,'(i0,2(",",i0),",",*(es25.16e3,:,","))',iostat=ios) parcels(p)%id,parcels(p)%parent_id, &
        parcels(p)%material_id,parcels(p)%position,parcels(p)%velocity,parcels(p)%mass,parcels(p)%temperature, &
        parcels(p)%multiplicity,masses(:,p)
      if (ios/=0) call fatal('Cannot write particle CSV')
    end do
    close(unit,iostat=ios)
    if (ios/=0) call fatal('Cannot close particle CSV')
    if (rank==0) then
      write(*,'(a,i0,a,es24.16,a,i0)') 'Steps: ',step,' Final time: ',time,' Last particle ID: ',last_id
      write(*,'(a,l1)') 'Regridded: ',regridded
      write(*,'(a,*(es24.16,1x))') 'Injected inventory: ',global_ledger
    end if
  contains
    subroutine check_output(path)
      character(len=*), intent(in) :: path
      if (selected_reactor_paths_alias_lexically(path,input_path).or. &
          selected_reactor_paths_alias_lexically(path,topology_path).or. &
          selected_reactor_paths_alias_lexically(path,spray_path)) call fatal('Output aliases input')
      if (len_trim(regrid_topology)>0) then
        if (selected_reactor_paths_alias_lexically(path,trim(regrid_topology))) call fatal('Output aliases regrid input')
      end if
    end subroutine
    subroutine fatal(message)
      character(len=*), intent(in) :: message
      integer :: status
      write(*,'(a,i0,2a)') 'AMR spray rank ',rank,': ',message
      call MPI_Abort(comm,2,status)
      error stop 2
    end subroutine
    subroutine read_topology(path,result)
      character(len=*), intent(in) :: path
      type(amr_patch_tree_3d), intent(out) :: result
      type(amr_patch_3d), allocatable :: patches(:)
      integer, allocatable :: parents(:),shapes(:,:)
      integer :: input,ios,nchild,peer,parent,dims(3),lower(3),upper(3),ratio
      logical :: valid
      open(newunit=input,file=path,status='old',action='read',iostat=ios)
      if (ios/=0) call fatal('Cannot open tree topology')
      read(input,*,iostat=ios) nchild
      if (ios/=0) call fatal('Invalid tree child count')
      if (nchild<0.or.nchild==huge(nchild)) call fatal('Invalid tree child count')
      allocate(patches(nchild),parents(nchild),shapes(3,nchild+1),stat=ios)
      if (ios/=0) call fatal('Cannot allocate tree topology')
      shapes(:,1)=[config%nx,config%ny,config%nz]
      do peer=1,nchild
        read(input,*,iostat=ios) parent,lower(1),upper(1),lower(2),upper(2),lower(3),upper(3),ratio
        if (ios/=0) call fatal('Invalid tree patch record')
        if (parent<1.or.parent>peer) call fatal('Tree parents must precede children')
        parents(peer)=parent;dims=shapes(:,parent)
        call initialize_amr_patch_3d(dims(1),dims(2),dims(3),lower(1),upper(1),lower(2),upper(2), &
          lower(3),upper(3),ratio,patches(peer),valid)
        if (.not.valid) call fatal('Invalid tree patch geometry')
        if (any(int(upper-lower+1,int64)*int(ratio,int64)>int(huge(1),int64))) call fatal('Tree patch extent overflow')
        shapes(:,peer+1)=[patches(peer)%fine_nx(),patches(peer)%fine_ny(),patches(peer)%fine_nz()]
      end do
      close(input)
      call initialize_amr_patch_tree_3d(shapes(:,1),parents,patches,result,valid)
      if (.not.valid) call fatal('Invalid tree topology')
    end subroutine
  end subroutine
end module
