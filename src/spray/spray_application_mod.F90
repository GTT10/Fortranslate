module spray_application_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use constants_mod, only: pelef_version
  use nasa7_thermo_mod, only: nasa7_species
  use elementary_kinetics_mod, only: elementary_reaction
  use gas_transport_mod, only: gas_transport_species, compatible_transport_database
  use mixture_thermo_mod, only: mass_fractions_from_mole_fractions, mixture_specific_gas_constant
  use reactive_1d_mod, only: reactive_nvar, reactive_nprim, reactive_primitive_to_conserved
  use reactive_3d_mod, only: compute_reactive_cfl_timestep_3d
  use spray_parcel_mod
  use spray_coupling_3d_mod, only: spray_system_totals
  use spray_reactive_3d_mod, only: advance_spray_reactive_3d, spray_cvode_available
  use spray_checkpoint_mod, only: write_spray_checkpoint, read_spray_checkpoint
  implicit none
  private
  public :: run_spray_application
contains
  subroutine run_spray_application(input_path,species,reactions,transport,mechanism_id)
    character(len=*), intent(in) :: input_path,mechanism_id
    type(nasa7_species), intent(in) :: species(:)
    type(elementary_reaction), intent(in) :: reactions(:)
    type(gas_transport_species), intent(in) :: transport(:)
    integer :: cells(3),parcel_count,maximum_steps,stop_after_steps,maximum_substeps
    integer :: unit,ios,i,j,k,s,steps,event,retry,history_unit,gas_unit,parcel_unit,nvar
    real(dp) :: domain_length(3),gas_temperature,gas_pressure,gas_velocity(3),shear_velocity
    real(dp) :: mole_amounts(32),liquid_density,liquid_cp,reference_temperature,reference_pressure,latent_heat
    real(dp) :: diameter_shape,diameter_min,diameter_max
    real(dp) :: drop_temperature,drop_diameter,initial_liquid_mass,injection_center(3),injection_axis(3)
    real(dp) :: injection_speed,cone_angle,mass_flow_rate,injection_start,injection_end,injection_interval
    real(dp) :: final_time,maximum_dt,cfl,rtol,atol,smagorinsky_constant,turbulent_prandtl,turbulent_schmidt
    real(dp) :: max_fractional_change,particle_cfl,minimum_diameter,time,dt,cfl_dt,next_event,tolerance,amount
    real(dp) :: spacing(3),sound,local_temperature,gas_constant,totals5(5),volume
    real(dp), allocatable :: state(:,:,:,:),temperature(:,:,:),x(:),y(:),primitive(:)
    real(dp), allocatable :: reference_totals(:),injected_totals(:),totals(:),expected_initial(:)
    type(spray_parcel), allocatable :: parcels(:),new_parcels(:)
    type(spray_liquid) :: liquid
    type(spray_options) :: options
    logical :: chemistry_enabled,transport_enabled,drag_enabled,heat_enabled,evaporation_enabled,abramzon_sirignano
    logical :: ok,exists
    character(len=24) :: mole_names(32),vapor_name,chemistry_method,hydro_reconstruction,riemann_solver
    character(len=1024) :: output_prefix,checkpoint_file,restart_file,message
    character(len=16384) :: context
    character(len=1050) :: paths(3)
    integer(int64) :: first_id
    namelist /spray/ cells,domain_length,hydro_reconstruction,riemann_solver,gas_temperature,gas_pressure, &
      gas_velocity,shear_velocity,mole_names,mole_amounts, &
      vapor_name,liquid_density,liquid_cp,reference_temperature,reference_pressure,latent_heat,drop_temperature, &
      drop_diameter,diameter_shape,diameter_min,diameter_max,initial_liquid_mass,parcel_count, &
      injection_center,injection_axis,injection_speed,cone_angle, &
      mass_flow_rate,injection_start,injection_end,injection_interval,final_time,maximum_dt,cfl,maximum_steps, &
      chemistry_enabled,transport_enabled,chemistry_method,rtol,atol,drag_enabled,heat_enabled,evaporation_enabled, &
      abramzon_sirignano,max_fractional_change,particle_cfl,minimum_diameter,maximum_substeps,smagorinsky_constant, &
      turbulent_prandtl,turbulent_schmidt,output_prefix,checkpoint_file,restart_file,stop_after_steps

    cells=[4,4,4]; domain_length=0.004_dp
    hydro_reconstruction='characteristic_plm'; riemann_solver='hllc'
    gas_temperature=1200.0_dp; gas_pressure=101325.0_dp; gas_velocity=0.0_dp; shear_velocity=0.0_dp
    mole_names=''; mole_names(1:3)=[character(len=24)::'H2','O2','N2']
    mole_amounts=0.0_dp; mole_amounts(1:3)=[2.0_dp,1.0_dp,3.76_dp]
    vapor_name='H2O'; liquid_density=1000.0_dp; liquid_cp=4184.0_dp
    reference_temperature=373.15_dp; reference_pressure=101325.0_dp; latent_heat=2.257e6_dp
    drop_temperature=330.0_dp; drop_diameter=20.0e-6_dp; initial_liquid_mass=1.0e-10_dp; parcel_count=4
    diameter_shape=0.0_dp; diameter_min=1.0e-12_dp; diameter_max=1.0_dp
    injection_center=0.0015_dp; injection_axis=[0.0_dp,0.0_dp,1.0_dp]; injection_speed=10.0_dp; cone_angle=0.1_dp
    mass_flow_rate=0.0_dp; injection_start=0.0_dp; injection_end=1.0e-5_dp; injection_interval=1.0e-6_dp
    final_time=1.0e-5_dp; maximum_dt=1.0e-7_dp; cfl=0.3_dp; maximum_steps=100000; stop_after_steps=0
    chemistry_enabled=.true.; transport_enabled=.true.; chemistry_method='implicit'; rtol=1.0e-7_dp; atol=1.0e-12_dp
    drag_enabled=.true.; heat_enabled=.true.; evaporation_enabled=.true.; abramzon_sirignano=.true.
    max_fractional_change=0.05_dp; particle_cfl=0.25_dp; minimum_diameter=1.0e-9_dp; maximum_substeps=100000
    smagorinsky_constant=0.0_dp; turbulent_prandtl=0.7_dp; turbulent_schmidt=0.7_dp
    output_prefix='spray'; checkpoint_file=''; restart_file=''
    open(newunit=unit,file=trim(input_path),status='old',action='read',iostat=ios)
    if (ios/=0) call fail('Cannot open spray input')
    read(unit,nml=spray,iostat=ios,iomsg=message)
    close(unit)
    if (ios/=0) call fail('Invalid spray namelist: '//trim(message))
    nvar=reactive_nvar(size(species))
    if (nvar<7) call fail('Spray requires 2--32 gas species')
    if (.not. compatible_transport_database(species,transport)) call fail('Incompatible gas transport database')
    if (trim(hydro_reconstruction)/='pcm' .and. trim(hydro_reconstruction)/='characteristic_plm') &
      call fail('Invalid hydro reconstruction')
    if (trim(hydro_reconstruction)=='characteristic_plm' .and. any(cells<3)) &
      call fail('Characteristic PLM requires at least three cells per direction')
    if (trim(riemann_solver)/='hllc' .and. trim(riemann_solver)/='rusanov' .and. trim(riemann_solver)/='pelec') &
      call fail('Invalid Riemann solver')
    if (any(cells<2) .or. any(cells>64)) call fail('Spray cell counts must be in [2,64]')
    if (.not. all(ieee_is_finite([domain_length,gas_temperature,gas_pressure,gas_velocity,shear_velocity, &
        mole_amounts,liquid_density,liquid_cp,reference_temperature,reference_pressure,latent_heat,drop_temperature, &
        drop_diameter,diameter_shape,diameter_min,diameter_max,initial_liquid_mass, &
        injection_center,injection_axis,injection_speed,cone_angle,mass_flow_rate, &
        injection_start,injection_end,injection_interval,final_time,maximum_dt,cfl,rtol,atol,smagorinsky_constant, &
        turbulent_prandtl,turbulent_schmidt,max_fractional_change,particle_cfl,minimum_diameter]))) &
      call fail('Nonfinite spray input')
    if (any(domain_length<1.0e-6_dp) .or. any(domain_length>1.0e3_dp)) call fail('Invalid domain length')
    if (gas_temperature<1.0_dp .or. gas_temperature>1.0e5_dp .or. gas_pressure<1.0_dp .or. &
        gas_pressure>1.0e10_dp .or. maxval(abs(gas_velocity))>1.0e5_dp .or. &
        abs(shear_velocity)>1.0e5_dp) call fail('Invalid gas state')
    if (maximum_dt<1.0e-15_dp .or. maximum_dt>1.0_dp .or. final_time<=0.0_dp .or. &
        final_time>1000.0_dp .or. cfl<=0.0_dp .or. cfl>0.5_dp .or. maximum_steps<1 .or. &
        maximum_steps>10000000 .or. stop_after_steps<0) call fail('Invalid time controls')
    if (rtol<=0.0_dp .or. rtol>0.01_dp .or. atol<=0.0_dp .or. atol>0.001_dp) call fail('Invalid chemistry tolerances')
    if (smagorinsky_constant<0.0_dp .or. smagorinsky_constant>0.5_dp .or. &
        min(turbulent_prandtl,turbulent_schmidt)<0.1_dp .or. &
        max(turbulent_prandtl,turbulent_schmidt)>10.0_dp) call fail('Invalid LES controls')
    if (trim(chemistry_method)/='implicit' .and. trim(chemistry_method)/='explicit' .and. &
        trim(chemistry_method)/='cvode') call fail('Unknown chemistry backend')
    if (trim(chemistry_method)=='cvode' .and. .not. spray_cvode_available()) &
      call fail('CVODE spray backend not built; enable PELEF_ENABLE_SUNDIALS')
    if (initial_liquid_mass<0.0_dp .or. initial_liquid_mass>1.0_dp .or. &
        parcel_count<1 .or. parcel_count>100000) call fail('Invalid initial injection')
    if (mass_flow_rate<0.0_dp .or. mass_flow_rate>1.0_dp .or. injection_start<0.0_dp .or. &
        injection_end<injection_start .or. injection_end>1000.0_dp .or. &
        injection_interval<1.0e-15_dp .or. injection_interval>1000.0_dp) call fail('Invalid injection schedule')
    if (mass_flow_rate>0.0_dp) then
      if (injection_end<=injection_start) call fail('Nonzero injection requires a nonempty time window')
      if ((injection_end-injection_start)/injection_interval>1000000.0_dp) call fail('Too many injection events')
    end if
    if (any(injection_center<0.0_dp) .or. any(injection_center>=domain_length)) &
      call fail('Injection center must lie in the periodic domain')
    call initialize_spray_liquid(species,trim(vapor_name),liquid_density,liquid_cp,reference_temperature, &
      reference_pressure,latent_heat,liquid,ok)
    if (.not. ok) call fail('Invalid liquid/vapor thermodynamics')
    options%drag=drag_enabled; options%heat=heat_enabled; options%evaporation=evaporation_enabled
    options%abramzon_sirignano=abramzon_sirignano; options%max_fractional_change=max_fractional_change
    options%particle_cfl=particle_cfl; options%minimum_diameter=minimum_diameter; options%maximum_substeps=maximum_substeps
    if (.not. valid_spray_options(options)) call fail('Invalid parcel integration controls')
    ! Validate injection geometry even for a zero-mass/gas-only run.
    call inject_cone_parcels(liquid,1.0e-12_dp,1,1_int64,injection_center,injection_axis,injection_speed, &
      cone_angle,drop_diameter,drop_temperature,new_parcels,ok,diameter_shape,diameter_min,diameter_max)
    if (.not. ok) call fail('Invalid parcel injection geometry or material state')
    deallocate(new_parcels)
    allocate(x(size(species)),y(size(species)),primitive(reactive_nprim(size(species))))
    x=0.0_dp
    if (any(mole_amounts<0.0_dp) .or. maxval(mole_amounts)>1.0e10_dp) call fail('Invalid gas mole amounts')
    do i=1,size(mole_names)
      if (len_trim(mole_names(i))==0) then
        if (mole_amounts(i)/=0.0_dp) call fail('Mole amount lacks a species name')
        cycle
      end if
      s=0
      do j=1,size(species)
        if (trim(species(j)%name)==trim(mole_names(i))) s=j
      end do
      if (s==0) call fail('Unknown gas species: '//trim(mole_names(i)))
      if (count(mole_names==mole_names(i))/=1) call fail('Duplicate input species name')
      x(s)=mole_amounts(i)
    end do
    if (sum(x)<=0.0_dp) call fail('Gas mole amounts sum to zero')
    x=x/sum(x)
    call mass_fractions_from_mole_fractions(species,x,y,ok)
    if (.not. ok) call fail('Invalid gas composition')
    gas_constant=mixture_specific_gas_constant(species,y,ok)
    if (.not. ok) call fail('Invalid gas molecular weight')
    spacing=domain_length/real(cells,dp); volume=product(spacing)
    allocate(state(nvar,cells(1),cells(2),cells(3)),temperature(cells(1),cells(2),cells(3)))
    allocate(reference_totals(nvar),injected_totals(nvar),totals(nvar))
    primitive(1)=gas_pressure/(gas_constant*gas_temperature); primitive(2:4)=gas_velocity
    primitive(5)=gas_pressure; primitive(6:)=y
    do k=1,cells(3)
      do j=1,cells(2)
        primitive(2)=gas_velocity(1)+shear_velocity*sin(2.0_dp*acos(-1.0_dp)*(real(j,dp)-0.5_dp)/real(cells(2),dp))
        do i=1,cells(1)
          call reactive_primitive_to_conserved(species,primitive,state(:,i,j,k),local_temperature,sound,ok)
          if (.not. ok) call fail('Initial gas state is outside the mechanism thermodynamic range')
          temperature(i,j,k)=local_temperature
        end do
      end do
    end do
    allocate(parcels(0)); time=0.0_dp; steps=0; event=0; injected_totals=0.0_dp
    if (initial_liquid_mass>0.0_dp) then
      call inject_cone_parcels(liquid,initial_liquid_mass,parcel_count,1_int64,injection_center,injection_axis, &
        injection_speed,cone_angle,drop_diameter,drop_temperature,parcels,ok,diameter_shape,diameter_min,diameter_max)
      if (.not. ok) call fail('Initial parcel injection failed')
    end if
    call spray_system_totals(state,spacing,parcels,liquid,reference_totals)
    expected_initial=reference_totals
    ! The checkpoint contract includes every physical, injection and step-control
    ! input. Paths, requested stop time/count and output scheduling are not physics.
    write(context,'(a,1x,a,1x,*(i0,1x))',iostat=ios) &
      'spray-v1/'//pelef_version,trim(mechanism_id),cells,parcel_count,maximum_substeps
    if (ios/=0) call fail('Cannot encode checkpoint context')
    do i=1,size(species)
      context=trim(context)//' '//trim(species(i)%name)
    end do
    call append_reals([domain_length,gas_temperature,gas_pressure,gas_velocity,shear_velocity,x,liquid_density,liquid_cp, &
      reference_temperature,reference_pressure,latent_heat,drop_temperature,drop_diameter, &
      diameter_shape,diameter_min,diameter_max,initial_liquid_mass, &
      injection_center,injection_axis,injection_speed,cone_angle,mass_flow_rate,injection_start,injection_end, &
      injection_interval,maximum_dt,cfl,rtol,atol,smagorinsky_constant,turbulent_prandtl,turbulent_schmidt, &
      max_fractional_change,particle_cfl,minimum_diameter])
    context=trim(context)//' '//trim(vapor_name)//' '//trim(chemistry_method)// &
      ' '//trim(hydro_reconstruction)//' '//trim(riemann_solver)
    block
      character(len=32) :: flags
      write(flags,'(7l2)') chemistry_enabled,transport_enabled,drag_enabled,heat_enabled,evaporation_enabled, &
        abramzon_sirignano,spray_cvode_available()
      context=trim(context)//' '//trim(flags)
    end block
    if (len_trim(restart_file)>0) then
      call read_spray_checkpoint(trim(restart_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message)
      if (.not. ok) call fail(trim(message))
      if (time>final_time .or. steps>maximum_steps .or. event>1000000) call fail('Restart exceeds requested limits')
      if (any(abs(reference_totals-expected_initial)>1.0e-11_dp* &
          max(abs(expected_initial),1.0e-30_dp))) call fail('Restart initial inventory mismatch')
      if (mass_flow_rate==0.0_dp .and. event/=0) call fail('Restart has an impossible injection event')
      amount=mass_flow_rate*min(real(event,dp)*injection_interval,injection_end-injection_start)
      if (abs(injected_totals(1)-amount)>1.0e-11_dp*max(amount,1.0e-30_dp)) &
        call fail('Restart injection mass ledger mismatch')
      if (int(size(parcels),int64)/=int(event,int64)*int(parcel_count,int64)+ &
          int(merge(parcel_count,0,initial_liquid_mass>0.0_dp),int64)) call fail('Restart parcel/event count mismatch')
      tolerance=32.0_dp*epsilon(1.0_dp)*max(time,maximum_dt,injection_interval)
      if (mass_flow_rate>0.0_dp) then
        next_event=injection_start+real(event,dp)*injection_interval
        if (event>0) then
          if (next_event-injection_interval>=injection_end .or. next_event-injection_interval>time+tolerance) &
            call fail('Restart injection event is in the future')
        end if
        if (next_event<injection_end .and. next_event<time-tolerance) call fail('Restart missed an injection event')
      end if
      call spray_system_totals(state,spacing,parcels,liquid,totals)
      if (abs(totals(1)-reference_totals(1)-injected_totals(1))> &
          1.0e-9_dp*max(reference_totals(1)+injected_totals(1),1.0e-30_dp)) call fail('Restart mass ledger mismatch')
      if (abs(totals(5)-reference_totals(5)-injected_totals(5))> &
          1.0e-9_dp*max(abs(reference_totals(5))+abs(injected_totals(5)),1.0e-20_dp)) &
        call fail('Restart energy ledger mismatch')
      if (any(parcels%id<1_int64)) call fail('Invalid restart parcel IDs')
      do i=1,size(parcels)
        if (parcels(i)%id/=int(i,int64)) call fail('Restart parcel IDs are not canonical')
        if (any(parcels(i)%position<0.0_dp) .or. any(parcels(i)%position>=domain_length)) &
          call fail('Restart parcel outside periodic domain')
      end do
    end if
    if (len_trim(output_prefix)==0) call fail('Empty output prefix')
    paths=[character(len=1050)::trim(output_prefix)//'.history.csv', &
      trim(output_prefix)//'.gas.csv',trim(output_prefix)//'.parcels.csv']
    do i=1,3
      inquire(file=trim(paths(i)),exist=exists)
      if (exists) call fail('Output exists; refusing to overwrite: '//trim(paths(i)))
    end do
    if (len_trim(checkpoint_file)>0) then
      inquire(file=trim(checkpoint_file),exist=exists)
      if (exists) call fail('Checkpoint destination already exists')
      if (any(paths==checkpoint_file)) call fail('Checkpoint path aliases an output')
    end if
    open(newunit=history_unit,file=trim(paths(1)),status='new',action='write',iostat=ios)
    if (ios/=0) call fail('Cannot create history output')
    write(history_unit,'(a)',iostat=ios) 'time,steps,gas_mass,liquid_mass,system_mass,px,py,pz,total_energy,'// &
      'injected_mass,mass_error,px_error,py_error,pz_error,energy_error,tmin,tmax,live_parcels,smd'
    if (ios/=0) call fail('Cannot write history header')
    call write_history()
    do while(time<final_time)
      if (steps>=maximum_steps) call fail('Maximum step count reached before final time')
      if (stop_after_steps>0) then
        if (steps>=stop_after_steps) exit
      end if
      tolerance=16.0_dp*epsilon(1.0_dp)*max(time,maximum_dt,injection_interval)
      next_event=huge(1.0_dp)
      if (mass_flow_rate>0.0_dp) then
        next_event=injection_start+real(event,dp)*injection_interval
        if (next_event>=injection_end) next_event=huge(1.0_dp)
        if (next_event<=time+tolerance) then
          amount=mass_flow_rate*min(injection_interval,injection_end-next_event)
          if (size(parcels)+parcel_count>1000000) call fail('Parcel storage limit reached')
          first_id=1_int64
          if (size(parcels)>0) first_id=parcels(size(parcels))%id+1_int64
          call inject_cone_parcels(liquid,amount,parcel_count,first_id,injection_center,injection_axis,injection_speed, &
            cone_angle,drop_diameter,drop_temperature,new_parcels,ok,diameter_shape,diameter_min,diameter_max)
          if (.not. ok) call fail('Scheduled parcel injection failed')
          do i=1,size(new_parcels)
            totals5=parcel_totals(new_parcels(i),liquid)
            injected_totals(1:5)=injected_totals(1:5)+totals5
            injected_totals(5+liquid%vapor_index)=injected_totals(5+liquid%vapor_index)+totals5(1)
          end do
          parcels=[parcels,new_parcels]; deallocate(new_parcels); event=event+1
          next_event=injection_start+real(event,dp)*injection_interval
          if (next_event>=injection_end) next_event=huge(1.0_dp)
        end if
      end if
      call compute_reactive_cfl_timestep_3d(species,state,temperature,cells(1),cells(2),cells(3), &
        spacing(1),spacing(2),spacing(3),cfl,cfl_dt,ok)
      if (.not. ok) call fail('Cannot compute gas CFL step')
      dt=min(maximum_dt,cfl_dt,final_time-time)
      if (next_event>time) dt=min(dt,next_event-time)
      ! Avoid a spurious vanishing final step from accumulated roundoff.
      if (abs(final_time-(time+dt)) <= 8.0_dp*epsilon(1.0_dp)*max(abs(final_time),abs(time))) dt=final_time-time
      do retry=1,30
        call advance_spray_reactive_3d(species,reactions,transport,liquid,options,state,temperature,parcels, &
          [0.0_dp,0.0_dp,0.0_dp],spacing,dt,chemistry_enabled,transport_enabled,trim(chemistry_method),rtol,atol, &
          smagorinsky_constant,turbulent_prandtl,turbulent_schmidt,ok,message,hydro_reconstruction,riemann_solver)
        if (ok) exit
        dt=0.5_dp*dt
        if (time+dt<=time .or. dt<1.0e-18_dp) exit
      end do
      if (.not. ok) call fail('Coupled step failed with rollback: '//trim(message))
      time=time+dt; steps=steps+1
      call write_history()
    end do
    close(history_unit,iostat=ios)
    if (ios/=0) call fail('History close failed')
    if (len_trim(checkpoint_file)>0) then
      call write_spray_checkpoint(trim(checkpoint_file),trim(context),species,state,temperature,parcels,time,steps,event, &
        reference_totals,injected_totals,ok,message)
      if (.not. ok) call fail(trim(message))
    end if
    call write_snapshots()
    write(*,'(a,es14.6,a,i0,a,i0)') 'Spray run completed: t=',time,', steps=',steps,', parcels=',size(parcels)
  contains
    subroutine fail(text)
      character(len=*), intent(in) :: text
      write(*,'(a)') trim(text)
      error stop 2
    end subroutine
    subroutine append_reals(values)
      real(dp), intent(in) :: values(:)
      character(len=8192) :: line
      write(line,'(*(es26.17e3,1x))',iostat=ios) values
      if (ios/=0 .or. len_trim(context)+len_trim(line)+1>len(context)) call fail('Checkpoint context too long')
      context=trim(context)//' '//trim(line)
    end subroutine
    subroutine write_history()
      real(dp) :: gas_mass,liquid_mass,error5(5),d2,d3,d,multiplicity,smd
      integer :: p,live
      call spray_system_totals(state,spacing,parcels,liquid,totals)
      gas_mass=sum(state(1,:,:,:))*volume; liquid_mass=0.0_dp; d2=0.0_dp; d3=0.0_dp; live=0
      do p=1,size(parcels)
        liquid_mass=liquid_mass+parcels(p)%mass*parcels(p)%multiplicity
        if (parcels(p)%mass<=0.0_dp) cycle
        live=live+1; d=parcel_diameter(parcels(p),liquid); multiplicity=parcels(p)%multiplicity
        d2=d2+multiplicity*d*d; d3=d3+multiplicity*d*d*d
      end do
      smd=0.0_dp
      if (d2>0.0_dp) smd=d3/d2
      error5=totals(1:5)-reference_totals(1:5)-injected_totals(1:5)
      write(history_unit,'(es26.17e3,a,i0,15(a,es26.17e3),a,i0,a,es26.17e3)',iostat=ios) time,',',steps, &
        ',',gas_mass,',',liquid_mass,',',totals(1),',',totals(2),',',totals(3),',',totals(4),',',totals(5), &
        ',',injected_totals(1),',',error5(1),',',error5(2),',',error5(3),',',error5(4),',',error5(5), &
        ',',minval(temperature),',',maxval(temperature),',',live,',',smd
      if (ios/=0) call fail('History write failed')
    end subroutine
    subroutine write_snapshots()
      integer :: p
      character(len=2048) :: header
      open(newunit=gas_unit,file=trim(paths(2)),status='new',action='write',iostat=ios)
      if (ios/=0) call fail('Cannot create gas snapshot')
      header='x,y,z,temperature,rho,rhou,rhov,rhow,rhoE'
      do p=1,size(species)
        header=trim(header)//',rhoY_'//trim(species(p)%name)
      end do
      write(gas_unit,'(a)',iostat=ios) trim(header)
      if (ios/=0) call fail('Gas header write failed')
      do k=1,cells(3)
        do j=1,cells(2)
          do i=1,cells(1)
            write(gas_unit,'(*(es26.17e3,:,a))',iostat=ios) &
              (real(i,dp)-0.5_dp)*spacing(1),',',(real(j,dp)-0.5_dp)*spacing(2),',', &
              (real(k,dp)-0.5_dp)*spacing(3),',',temperature(i,j,k),',',(state(p,i,j,k),',',p=1,nvar-1),state(nvar,i,j,k)
            if (ios/=0) call fail('Gas snapshot write failed')
          end do
        end do
      end do
      close(gas_unit,iostat=ios)
      if (ios/=0) call fail('Gas snapshot close failed')
      open(newunit=parcel_unit,file=trim(paths(3)),status='new',action='write',iostat=ios)
      if (ios/=0) call fail('Cannot create parcel snapshot')
      write(parcel_unit,'(a)',iostat=ios) 'id,x,y,z,u,v,w,mass,temperature,multiplicity,diameter'
      if (ios/=0) call fail('Parcel header write failed')
      do p=1,size(parcels)
        write(parcel_unit,'(i0,10(a,es26.17e3))',iostat=ios) parcels(p)%id,',',parcels(p)%position(1), &
          ',',parcels(p)%position(2),',',parcels(p)%position(3),',',parcels(p)%velocity(1), &
          ',',parcels(p)%velocity(2),',',parcels(p)%velocity(3),',',parcels(p)%mass,',',parcels(p)%temperature, &
          ',',parcels(p)%multiplicity,',',parcel_diameter(parcels(p),liquid)
        if (ios/=0) call fail('Parcel snapshot write failed')
      end do
      close(parcel_unit,iostat=ios)
      if (ios/=0) call fail('Parcel snapshot close failed')
    end subroutine
  end subroutine run_spray_application
end module spray_application_mod
