module spray_eb_application_mod
  use precision_mod, only: dp
  use, intrinsic :: iso_c_binding, only: c_char,c_size_t,c_ptr,c_associated,c_null_char
  use simulation_config_selected_reactor_mod, only: selected_reactor_paths_alias_lexically
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species
  use reactive_1d_mod, only: reactive_nvar
  use simulation_config_reactive_eb_3d_mod, only: reactive_eb_3d_config,build_configured_eb_geometry_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d
  use reactive_eb_3d_driver_mod, only: initialize_reactive_eb_density_sheet_3d,write_reactive_eb_3d_csv
  use reactive_eb_cfl_3d_mod, only: compute_reactive_eb_cfl_timestep_3d
  use eb_reactive_transport_3d_mod, only: reactive_eb_transport_timestep_3d
  use spray_parcel_mod
  use spray_mixture_mod, only: liquid_mixture_properties
  use spray_eb_reactive_mod, only: advance_eb_spray_reactive
  implicit none
  private
  public :: run_eb_spray_application
  interface
    function system_getcwd(buffer,length) bind(C,name='getcwd') result(pointer)
      import c_char,c_size_t,c_ptr
      character(kind=c_char), intent(out) :: buffer(*)
      integer(c_size_t), value :: length
      type(c_ptr) :: pointer
    end function
  end interface
contains
  subroutine run_eb_spray_application(input_path,spray_path,config,species,reactions,transport,mole_fractions,integrator)
    character(len=*), intent(in) :: input_path,spray_path
    type(reactive_eb_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    real(dp), optional, intent(in) :: mole_fractions(:)
    character(len=*), optional, intent(in) :: integrator
    integer :: nc=1,number=1
    character(len=24) :: vapor_names(16)='H2O'
    real(dp) :: densities(16)=1000,heat_capacities(16)=4184,reference_temperatures(16)=373.15_dp
    real(dp) :: reference_pressures(16)=101325,latent_heats(16)=2.257e6_dp,fractions(16)=0
    real(dp) :: total_mass=1e-10_dp,center(3)=0,axis(3)=[1._dp,0._dp,0._dp],speed=1,half_angle=0
    real(dp) :: diameter=1e-5_dp,liquid_temperature=350,maximum_dt=1e-8_dp
    type(spray_options) :: options
    namelist /eb_spray/ nc,number,vapor_names,densities,heat_capacities,reference_temperatures, &
      reference_pressures,latent_heats,fractions,total_mass,center,axis,speed,half_angle,diameter, &
      liquid_temperature,maximum_dt,options
    type(eb_geometry_3d) :: g
    type(spray_liquid), allocatable :: liquids(:)
    type(spray_liquid) :: effective
    type(spray_parcel), allocatable :: drops(:)
    real(dp), allocatable :: phi(:,:,:),state(:,:,:,:),temperature(:,:,:),masses(:,:)
    real(dp), allocatable :: boundary(:),wall(:),chemical(:),sb(:),sw(:),sc(:)
    real(dp) :: time,dt,hydro_dt,transport_dt,diffusivity,theta,rho,cp,e
    integer :: unit,status,close_status,nv,n,p,step,attempt,particle_steps
    logical :: ok
    character(len=1024) :: message
    character(len=24) :: component_label
    character(len=:), allocatable :: header
    ! Configuration and namelist objects are local to one invocation.
    nc=1;number=1;fractions=0;fractions(1)=1
    vapor_names='H2O';densities=1000;heat_capacities=4184;reference_temperatures=373.15_dp
    reference_pressures=101325;latent_heats=2.257e6_dp;total_mass=1e-10_dp;center=0
    axis=[1._dp,0._dp,0._dp];speed=1;half_angle=0;diameter=1e-5_dp;liquid_temperature=350
    maximum_dt=1e-8_dp;options=spray_options()
    call check_output(config%output_file)
    call check_output(trim(config%output_file)//'.particles.csv')
    if (config%geometry_type/='level_set') error stop 'EB spray requires nodal level-set geometry'
    if (config%checkpoint_interval_steps/=0.or.len_trim(config%restart_file)>0.or.config%stop_after_checkpoint) &
      error stop 'EB spray checkpoint is not yet supported'
    open(newunit=unit,file=spray_path,status='old',action='read',iostat=status)
    if(status/=0) error stop 'Cannot open EB spray input'
    read(unit,nml=eb_spray,iostat=status)
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'Invalid EB spray namelist'
    if(nc<1.or.nc>16.or.number<1.or.number>125000) error stop 'Invalid EB spray population or components'
    if(.not.ieee_is_finite(maximum_dt)) error stop 'Invalid EB spray maximum_dt'
    if(maximum_dt<=0.or.maximum_dt>1e6_dp) error stop 'Invalid EB spray maximum_dt'
    if(.not.all(ieee_is_finite(fractions)).or.any(fractions<0)) error stop 'Invalid liquid fractions'
    if(abs(sum(fractions(:nc))-1)>1e-12_dp.or.any(fractions(nc+1:)/=0)) error stop 'Liquid fractions must sum to one'
    if(.not.valid_spray_options(options).or.options%khrt_children) error stop 'Unsupported EB spray options'
    call build_configured_eb_geometry_3d(config,g,ok,phi)
    if(.not.ok) error stop 'Cannot build EB spray geometry'
    allocate(liquids(nc))
    do n=1,nc
      call initialize_spray_liquid(species,vapor_names(n),densities(n),heat_capacities(n),reference_temperatures(n), &
        reference_pressures(n),latent_heats(n),liquids(n),ok)
      if(.not.ok) error stop 'Invalid EB spray liquid table'
    end do
    call liquid_mixture_properties(liquids,fractions(:nc),liquid_temperature,rho,cp,e,ok)
    if(.not.ok) error stop 'Invalid injected liquid mixture'
    effective=liquids(1);effective%density=rho;effective%heat_capacity=cp
    effective%reference_temperature=liquid_temperature;effective%reference_energy=e
    call inject_cone_parcels(effective,total_mass,number,1_int64,center,axis,speed,half_angle,diameter, &
      liquid_temperature,drops,ok)
    if(.not.ok) error stop 'Invalid EB injection'
    allocate(masses(nc,number))
    do p=1,number
      masses(:,p)=drops(p)%mass*fractions(:nc)
    end do
    nv=reactive_nvar(size(species))
    allocate(state(nv,g%nx,g%ny,g%nz),temperature(g%nx,g%ny,g%nz))
    allocate(boundary(nv),wall(nv),chemical(nv),sb(nv),sw(nv),sc(nv))
    call initialize_reactive_eb_density_sheet_3d(species,config,g,state,temperature,ok,mole_fractions)
    if(.not.ok) error stop 'Invalid EB initial gas'
    time=0;step=0;boundary=0;wall=0;chemical=0
    do while(time<config%final_time)
      if(step>=config%maximum_steps) error stop 'EB spray maximum_steps exhausted'
      call compute_reactive_eb_cfl_timestep_3d(species,state,temperature,g,config%cfl,hydro_dt,ok)
      if(.not.ok) error stop 'Invalid EB gas CFL bound'
      dt=min(hydro_dt,maximum_dt,config%final_time-time)
      if(config%transport_enabled) then
        call reactive_eb_transport_timestep_3d(species,transport,state,temperature,g,config%transport_cfl, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          transport_dt,diffusivity,ok)
        if(.not.ok) error stop 'Invalid EB transport bound'
        dt=min(dt,transport_dt)
      end if
      do attempt=1,30
        ok=.false.
        if(dt<=0.or.time+dt<=time) exit
        theta=1;sb=0;sw=0;sc=0
        call advance_eb_spray_reactive(species,reactions,transport,liquids,options,g,phi,state,temperature, &
          drops,masses,dt,config%riemann_solver,config%redistribution,config%chemistry_enabled,config%transport_enabled, &
          config%viscosity_enabled,config%thermal_conduction_enabled,config%species_diffusion_enabled, &
          config%barodiffusion_enabled,config%chemistry_relative_tolerance,config%chemistry_absolute_tolerance, &
          particle_steps,theta,sb,sw,sc,ok,config%state_redist_target_volume_fraction,integrator)
        if(ok) exit
        dt=.5_dp*dt
      end do
      if(.not.ok) error stop 'EB spray coupled step rejected'
      time=time+dt;step=step+1;boundary=boundary+sb;wall=wall+sw;chemical=chemical+sc
    end do
    call write_reactive_eb_3d_csv(config%output_file,species,g,state,temperature,time,ok,message)
    if(.not.ok) error stop 'Cannot write EB spray gas CSV'
    open(newunit=unit,file=trim(config%output_file)//'.particles.csv',status='replace',action='write',iostat=status)
    if(status/=0) error stop 'Cannot write EB spray particle CSV'
    header='id,x,y,z,u,v,w,mass,temperature,multiplicity'
    do n=1,nc
      write(component_label,'(i0)') n
      header=header//',component_mass_'//trim(component_label)
    end do
    write(unit,'(a)',iostat=status) header
    do p=1,size(drops)
      if(status/=0) exit
      write(unit,'(i0,*(",",es25.16,:))',iostat=status) drops(p)%id,drops(p)%position,drops(p)%velocity, &
        drops(p)%mass,drops(p)%temperature,drops(p)%multiplicity,masses(:,p)
    end do
    close(unit,iostat=close_status)
    if(status/=0.or.close_status/=0) error stop 'EB particle CSV write failed'
    write(*,'(a,i0,a,es25.16)') 'EB spray steps: ',step,' time: ',time
    write(*,'(a,*(es25.16,1x))') 'Boundary inventory: ',boundary
    write(*,'(a,*(es25.16,1x))') 'Wall inventory: ',wall
    write(*,'(a,*(es25.16,1x))') 'Chemistry inventory: ',chemical
  contains
    function absolute_path(path) result(value)
      character(len=*), intent(in) :: path
      character(len=:), allocatable :: value
      character(kind=c_char) :: buffer(4096)
      type(c_ptr) :: pointer
      integer :: index
      if (len_trim(path)==0) error stop 'Empty EB path'
      if (path(1:1)=='/') then
        value=trim(path)
      else
        pointer=system_getcwd(buffer,int(size(buffer),c_size_t))
        if (.not.c_associated(pointer)) error stop 'Cannot resolve EB working directory'
        value=''
        do index=1,size(buffer)
          if (buffer(index)==c_null_char) exit
          value=value//buffer(index)
        end do
        value=value//'/'//trim(path)
      end if
      if (len(value)>1024) error stop 'EB path too long'
    end function
    subroutine check_output(path)
      character(len=*), intent(in) :: path
      if(selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(input_path)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(spray_path)).or. &
         selected_reactor_paths_alias_lexically(absolute_path(path),absolute_path(config%level_set_file))) &
        error stop 'EB output aliases input'
    end subroutine
  end subroutine
end module
