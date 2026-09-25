module simulation_config_reactive_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use selected_composition_mod, only: &
    selected_composition_max_species, selected_composition_name_length, &
    validate_selected_composition_fields, resolve_selected_composition
  implicit none
  private

  integer, parameter, public :: reactive_3d_max_species = &
    selected_composition_max_species

  type, public :: reactive_3d_config
    integer :: nx = 12
    integer :: ny = 12
    integer :: nz = 12
    integer :: maximum_steps = 100000
    real(dp) :: x_lower = 0.0_dp
    real(dp) :: x_upper = 0.01_dp
    real(dp) :: y_lower = 0.0_dp
    real(dp) :: y_upper = 0.01_dp
    real(dp) :: z_lower = 0.0_dp
    real(dp) :: z_upper = 0.01_dp
    real(dp) :: final_time = 1.0e-6_dp
    real(dp) :: cfl = 0.30_dp
    character(len=32) :: problem = "entropy_wave"
    character(len=32) :: reconstruction = "pcm"
    character(len=32) :: limiter = "mc"
    character(len=32) :: thermo_model = "full_h2o2"
    character(len=32) :: riemann_solver = "rusanov"
    character(len=32) :: boundary_condition = "periodic"
    character(len=32) :: boundary_faces(6) = ""
    real(dp) :: wall_temperature(6)=300.0_dp
    real(dp) :: inflow_mach(6)=-1.0_dp
    real(dp) :: difmag=-1.0_dp ! -1 selects the solver policy; nonnegative overrides PPM.
    real(dp) :: outlet_pressure(6)=0.0_dp
    real(dp) :: wall_velocity(3,6)=0.0_dp
    real(dp) :: inflow_temperature(6)=1000.0_dp, inflow_pressure(6)=101325.0_dp
    real(dp) :: inflow_velocity(3,6)=0.0_dp
    integer :: inflow_composition_count=0
    character(len=selected_composition_name_length) :: inflow_composition_species(reactive_3d_max_species)=''
    real(dp) :: inflow_composition_mole_fractions(reactive_3d_max_species)=0.0_dp
    logical :: chemistry_enabled = .false.
    real(dp) :: chemistry_relative_tolerance = 2.0e-7_dp
    real(dp) :: chemistry_absolute_tolerance = 1.0e-12_dp
    logical :: transport_enabled = .false.
    logical :: viscosity_enabled = .true.
    logical :: thermal_conduction_enabled = .true.
    logical :: species_diffusion_enabled = .true.
    logical :: barodiffusion_enabled = .true.
    real(dp) :: transport_cfl = 0.35_dp
    real(dp) :: initial_temperature = 1000.0_dp
    real(dp) :: initial_pressure = 101325.0_dp
    real(dp) :: initial_velocity_x = 300.0_dp
    real(dp) :: initial_velocity_y = 200.0_dp
    real(dp) :: initial_velocity_z = -100.0_dp
    real(dp) :: density_wave_amplitude = 0.08_dp
    integer :: wave_number_x = 1
    integer :: wave_number_y = 1
    integer :: wave_number_z = 1
    real(dp) :: hotspot_temperature_rise = 250.0_dp
    real(dp) :: hotspot_center_x = 0.005_dp
    real(dp) :: hotspot_center_y = 0.005_dp
    real(dp) :: hotspot_center_z = 0.005_dp
    real(dp) :: hotspot_width = 0.0012_dp
    real(dp) :: x_h2 = 0.29570_dp
    real(dp) :: x_h = 1.0e-5_dp
    real(dp) :: x_o = 1.0e-5_dp
    real(dp) :: x_o2 = 0.14784_dp
    real(dp) :: x_oh = 1.0e-5_dp
    real(dp) :: x_h2o = 0.0_dp
    real(dp) :: x_ho2 = 0.0_dp
    real(dp) :: x_h2o2 = 0.0_dp
    real(dp) :: x_ar = 0.0_dp
    real(dp) :: x_n2 = 0.55643_dp
    integer :: composition_count = 0
    character(len=selected_composition_name_length) :: &
      composition_species(reactive_3d_max_species) = ""
    real(dp) :: &
      composition_mole_fractions(reactive_3d_max_species) = 0.0_dp
    character(len=512) :: output_file = "reactive_entropy_wave_3d.csv"
  end type reactive_3d_config

  public :: fill_pressure_outlet_ghosts_3d
  public :: characteristic_pressure_ghost_3d
  public :: valid_wall_velocity_3d
  public :: fill_reactive_boundary_ghosts_3d
  public :: resolved_reactive_boundaries_3d, valid_reactive_boundaries_3d
  public :: read_reactive_3d_configuration
  public :: validate_reactive_3d_configuration
  public :: reactive_3d_mole_fractions
  public :: resolve_reactive_3d_selected_composition

  public :: fill_mach_inlet_boundary_ghosts_3d
  public :: mach_inlet_ghost_3d
  public :: map_reactive_boundary_donor_3d
contains

  ! Frozen-composition acoustic pressure extrapolation, as in pinned
  ! PeleC EB-ConvergingNozzle/prob.H. distance_ratio=2 is the first ghost center.
  ! Pinned EB-ConvergingNozzle inlet: interior composition, prescribed T/Mach,
  ! extrapolated pressure by default, or an explicit pressure target.
  subroutine mach_inlet_ghost_3d(species,interior,temperature,face,inlet_temperature,mach,ghost,ghost_temperature,ok, &
      inlet_pressure,velocity_fluctuation,total_velocity)
    use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive,reactive_primitive_to_conserved
    use mixture_thermo_mod, only: mixture_mass_properties
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: interior(:),temperature,inlet_temperature,mach
    integer, intent(in) :: face
    real(dp), intent(inout) :: ghost(:),ghost_temperature
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inlet_pressure,velocity_fluctuation(3),total_velocity(3)
    real(dp) :: q(reactive_nvar(size(species))),candidate(reactive_nvar(size(species)))
    real(dp) :: recovered,sound,pressure,mw,r,cp,cv,gamma,h,e,entropy,new_t,new_sound
    integer :: axis
    logical :: valid
    ok=.false.
    if(face<1.or.face>6.or.size(interior)/=size(q).or.size(ghost)/=size(q)) return
    if(.not.all(ieee_is_finite([inlet_temperature,mach]))) return
    if(inlet_temperature<=0.or.mach<0.or.mach>100) return
    if(present(velocity_fluctuation).and.present(total_velocity)) return
    if(present(velocity_fluctuation)) then
      if(.not.all(ieee_is_finite(velocity_fluctuation))) return
    end if
    if(present(total_velocity)) then
      if(.not.all(ieee_is_finite(total_velocity))) return
    end if
    call reactive_conserved_to_primitive(species,interior,temperature,q,recovered,sound,valid)
    if(.not.valid) return
    pressure=q(5)
    if(present(inlet_pressure)) pressure=inlet_pressure
    if(.not.ieee_is_finite(pressure)) return
    if(pressure<=0.or.pressure>1e15_dp) return
    call mixture_mass_properties(species,q(6:),inlet_temperature,mw,r,cp,cv,gamma,h,e,entropy,valid)
    if(.not.valid) return
    q(1)=pressure/(r*inlet_temperature);q(5)=pressure;q(2:4)=0
    axis=(face+1)/2
    q(axis+1)=merge(1._dp,-1._dp,modulo(face,2)==1)*mach*sqrt(gamma*r*inlet_temperature)
    if(present(velocity_fluctuation)) q(2:4)=q(2:4)+velocity_fluctuation
    if(present(total_velocity)) q(2:4)=total_velocity
    call reactive_primitive_to_conserved(species,q,candidate,new_t,new_sound,valid)
    if(.not.valid) return
    ghost=candidate;ghost_temperature=new_t;ok=.true.
  end subroutine mach_inlet_ghost_3d

  ! Mach >= 0 enables this inlet; negative entries retain ordinary face behavior.
  subroutine fill_mach_inlet_boundary_ghosts_3d(species,u,t,faces,mach,ghost,gt,ok, &
      inflow_state,inflow_temperature,inflow_pressure,ghost_width,viscous_wall,wall_velocity)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: u(:,:,:,:),t(:,:,:),mach(6)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(inout) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: ok
    real(dp), optional, intent(in) :: inflow_state(:,:),inflow_temperature(:),inflow_pressure(6),wall_velocity(:,:)
    integer, optional, intent(in) :: ghost_width
    logical, optional, intent(in) :: viscous_wall
    real(dp), allocatable :: candidate(:,:,:,:),temps(:,:,:)
    real(dp) :: pressure(6),shift(3),momentum(3),rho
    character(len=32) :: ordinary(6)
    integer :: n(3),dims(3),width,i,j,k,q(3),inlet,axis
    logical :: valid,reflect(3),viscous
    ok=.false.;n=shape(t);width=2;pressure=0;viscous=.false.
    if(present(ghost_width)) width=ghost_width
    if(present(inflow_pressure)) pressure=inflow_pressure
    if(present(viscous_wall)) viscous=viscous_wall
    if(width<1.or.width>128.or.any(n<2).or.any(n>huge(1)-256)) return
    dims=n+2*width
    if(any(shape(u)/=[5+size(species),n]).or.any(shape(ghost)/=[5+size(species),dims])) return
    if(any(shape(gt)/=dims)) return
    if(.not.all(ieee_is_finite(mach)).or..not.all(ieee_is_finite(pressure))) return
    if(any(mach>100).or.any(mach>=0.and.faces/='inflow')) return
    if(any(pressure<0).or.any(pressure>1e15_dp).or.any(pressure>0.and.mach<0)) return
    if(any(mach>=0)) then
      if(.not.present(inflow_temperature)) return
      if(size(inflow_temperature)/=6) return
      if(.not.all(ieee_is_finite(inflow_temperature))) return
      if(any(inflow_temperature<=0.and.mach>=0)) return
    end if
    ordinary=faces
    where(mach>=0) ordinary='outflow'
    allocate(candidate,mold=ghost);allocate(temps,mold=gt)
    call fill_reactive_boundary_ghosts_3d(u,t,ordinary,candidate,temps,valid,viscous_wall=viscous, &
      inflow_state=inflow_state,inflow_temperature=inflow_temperature,wall_velocity=wall_velocity,ghost_width=width)
    if(.not.valid) return
    do k=1,dims(3)
      do j=1,dims(2)
        do i=1,dims(1)
          call map_reactive_boundary_donor_3d([i,j,k]-width,n,faces,viscous,q,reflect,shift,inlet,valid,wall_velocity)
          if(.not.valid) return
          if(inlet==0) cycle
          if(mach(inlet)<0) cycle
          if(pressure(inlet)>0) then
            call mach_inlet_ghost_3d(species,u(:,q(1),q(2),q(3)),t(q(1),q(2),q(3)),inlet, &
              inflow_temperature(inlet),mach(inlet),candidate(:,i,j,k),temps(i,j,k),valid,inlet_pressure=pressure(inlet))
          else
            call mach_inlet_ghost_3d(species,u(:,q(1),q(2),q(3)),t(q(1),q(2),q(3)),inlet, &
              inflow_temperature(inlet),mach(inlet),candidate(:,i,j,k),temps(i,j,k),valid)
          end if
          if(.not.valid) return
          momentum=candidate(2:4,i,j,k)
          do axis=1,3
            if(reflect(axis)) candidate(axis+1,i,j,k)=-candidate(axis+1,i,j,k)
          end do
          if(any(abs(shift)>0)) then
            rho=candidate(1,i,j,k)
            candidate(2:4,i,j,k)=candidate(2:4,i,j,k)+rho*shift
            candidate(5,i,j,k)=candidate(5,i,j,k)+(sum(candidate(2:4,i,j,k)**2)-sum(momentum**2))/(2*rho)
          end if
        end do
      end do
    end do
    if(.not.all(ieee_is_finite(candidate)).or..not.all(ieee_is_finite(temps))) return
    ghost=candidate;gt=temps;ok=.true.
  end subroutine fill_mach_inlet_boundary_ghosts_3d

  subroutine characteristic_pressure_ghost_3d(species,interior,temperature,face,pressure,distance_ratio, &
      ghost,ghost_temperature,ok)
    use reactive_1d_mod, only: reactive_nvar,reactive_conserved_to_primitive,reactive_primitive_to_conserved
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: interior(:),temperature,pressure,distance_ratio
    integer, intent(in) :: face
    real(dp), intent(inout) :: ghost(:),ghost_temperature
    logical, intent(out) :: ok
    real(dp) :: q(reactive_nvar(size(species))),candidate(reactive_nvar(size(species)))
    real(dp) :: recovered,sound,new_t,new_sound,dpres,rho,sign_normal
    integer :: axis
    ok=.false.
    if(face<1.or.face>6) return
    if(size(interior)/=size(q).or.size(ghost)/=size(q)) return
    if(.not.all(ieee_is_finite([pressure,distance_ratio]))) return
    if(pressure<=0.or.pressure>1e15_dp.or.distance_ratio<0.or.distance_ratio>1e6_dp) return
    call reactive_conserved_to_primitive(species,interior,temperature,q,recovered,sound,ok)
    if(.not.ok) return
    ok=.false.;rho=q(1);dpres=distance_ratio*(pressure-q(5));axis=(face+1)/2
    sign_normal=merge(-1._dp,1._dp,modulo(face,2)==1)
    q(1)=rho+dpres/(sound*sound)
    q(axis+1)=q(axis+1)-sign_normal*dpres/(rho*sound)
    q(5)=q(5)+dpres
    call reactive_primitive_to_conserved(species,q,candidate,new_t,new_sound,ok)
    if(.not.ok) return
    ghost=candidate;ghost_temperature=new_t
  end subroutine

  ! Apply after ordinary ghost filling; face order defines corner composition.
  ! boundary_cells gives the first/last physical-domain cell indices in the
  ! padded buffer. Unspecified bounds use the patch interior and width=2.
  subroutine fill_pressure_outlet_ghosts_3d(species,n,faces,pressure,ghost,gt,ok,ghost_width,boundary_cells)
    type(nasa7_species), intent(in) :: species(:)
    integer, intent(in) :: n(3)
    character(len=*), intent(in) :: faces(6)
    real(dp), optional, intent(in) :: pressure(6)
    real(dp), intent(inout) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: ok
    integer, optional, intent(in) :: ghost_width,boundary_cells(2,3)
    integer :: face,axis,layer,i,j,k,idx(3),donor(3),dims(3),width,bounds(2,3),layers
    real(dp) :: state(size(ghost,1)),temperature
    ok=.false.;width=2
    if(present(ghost_width)) width=ghost_width
    if(width<1.or.width>128) return
    if(any(n>huge(1)-2*width)) return
    dims=n+2*width
    bounds(1,:)=width+1;bounds(2,:)=width+n
    if(present(boundary_cells)) bounds=boundary_cells
    if(any(bounds<1).or.any(bounds(1,:)>dims).or.any(bounds(2,:)>dims)) return
    if(any(bounds(2,:)<bounds(1,:))) return
    if(any(n<1).or.any(shape(gt)/=dims)) return
    if(any(shape(ghost)/=[5+size(species),dims])) return
    if(.not.present(pressure)) then
      ok=.true.;return
    end if
    if(.not.all(ieee_is_finite(pressure))) return
    if(any(pressure<0).or.any(pressure>1e15_dp)) return
    if(any(pressure>0.and.faces/='outflow')) return
    do face=1,6
      if(pressure(face)<=0) cycle
      axis=(face+1)/2
      do k=1,dims(3)
      do j=1,dims(2)
      do i=1,dims(1)
        donor=[i,j,k]
        if(donor(axis)/=1) cycle
        donor(axis)=bounds(1,axis);layers=donor(axis)-1
        if(modulo(face,2)==0) then
          donor(axis)=bounds(2,axis);layers=dims(axis)-donor(axis)
        end if
        state=ghost(:,donor(1),donor(2),donor(3));temperature=gt(donor(1),donor(2),donor(3))
        do layer=1,layers
          idx=donor;idx(axis)=donor(axis)-layer
          if(modulo(face,2)==0) idx(axis)=donor(axis)+layer
          call characteristic_pressure_ghost_3d(species,state,temperature,face,pressure(face),real(2*layer,dp), &
            ghost(:,idx(1),idx(2),idx(3)),gt(idx(1),idx(2),idx(3)),ok)
          if(.not.ok) return
        end do
      end do
      end do
      end do
    end do
    ok=.true.
  end subroutine

  pure subroutine fill_reactive_boundary_ghosts_3d(u,t,faces,ghost,gt,ok,viscous_wall, &
      inflow_state,inflow_temperature,wall_velocity,ghost_width)
    real(dp), intent(in) :: u(:,:,:,:),t(:,:,:)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(out) :: ghost(:,:,:,:),gt(:,:,:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: inflow_state(:,:),inflow_temperature(:)
    logical, intent(in), optional :: viscous_wall
    real(dp), intent(in), optional :: wall_velocity(:,:)
    integer, optional, intent(in) :: ghost_width
    integer :: width
    real(dp) :: shift(3), old_momentum(3),rho
    logical :: use_viscous_wall
    integer :: n(3),dims(3),i,j,k,a,idx(3),q(3),inlet
    logical :: reflect(3)
    ok=.false.; ghost=0.0_dp; gt=0.0_dp
    use_viscous_wall=.false.
    if (present(viscous_wall)) use_viscous_wall=viscous_wall
    if (present(wall_velocity)) then
      if (.not.valid_wall_velocity_3d(faces,wall_velocity)) return
    end if
    width=2
    if (present(ghost_width)) width=ghost_width
    if (width<1.or.width>128) return
    n=shape(t)
    if (any(n>huge(1)-2*width)) return
    dims=n+2*width
    if (any(faces=='inflow')) then
      if (.not.present(inflow_state).or..not.present(inflow_temperature)) return
      if (any(shape(inflow_state)/=[size(u,1),6]).or.size(inflow_temperature)/=6) return
      if (.not.all(ieee_is_finite(inflow_state)).or..not.all(ieee_is_finite(inflow_temperature))) return
    end if
    if (any(n<2) .or. size(u,1)<5) return
    if (any(shape(u)/=[size(u,1),n]) .or. any(shape(ghost)/=[size(u,1),dims])) return
    if (any(shape(gt)/=dims) .or. .not.valid_reactive_boundaries_3d(faces)) return
      do k=1,dims(3)
        do j=1,dims(2)
          do i=1,dims(1)
            idx=[i,j,k]-width
            call map_reactive_boundary_donor_3d(idx,n,faces,use_viscous_wall,q,reflect,shift,inlet,ok,wall_velocity)
            if (.not.ok) return
            ok=.false.
            ghost(:,i,j,k)=u(:,q(1),q(2),q(3)); gt(i,j,k)=t(q(1),q(2),q(3))
            if (inlet>0) then
              ghost(:,i,j,k)=inflow_state(:,inlet); gt(i,j,k)=inflow_temperature(inlet)
            end if
            old_momentum=ghost(2:4,i,j,k)
            do a=1,3
              if (reflect(a)) ghost(a+1,i,j,k)=-ghost(a+1,i,j,k)
            end do
            if (any(abs(shift)>0.0_dp)) then
              rho=ghost(1,i,j,k)
              if (rho<=0.0_dp) return
              ghost(2:4,i,j,k)=ghost(2:4,i,j,k)+rho*shift
              ghost(5,i,j,k)=ghost(5,i,j,k)+ &
                (sum(ghost(2:4,i,j,k)**2)-sum(old_momentum**2))/(2.0_dp*rho)
              if (.not.all(ieee_is_finite(ghost(:,i,j,k)))) return
            end if
          end do
        end do
      end do
    ok=.true.
  end subroutine fill_reactive_boundary_ghosts_3d

  pure subroutine map_reactive_boundary_donor_3d(idx,n,faces,viscous,q,reflect,shift,inlet,ok,wall_velocity)
    integer, intent(in) :: idx(3),n(3)
    character(len=*), intent(in) :: faces(6)
    logical, intent(in) :: viscous
    integer, intent(out) :: q(3),inlet
    logical, intent(out) :: reflect(3),ok
    real(dp), intent(out) :: shift(3)
    real(dp), optional, intent(in) :: wall_velocity(:,:)
    integer :: a,face,bounces,bounce,wall_faces(128)
    ok=.false.;q=idx;reflect=.false.;shift=0;inlet=0
    if(any(n<1).or.any(n>huge(1)-256)) return
    if(any(idx<1-128).or.any(idx>n+128)) return
    if(.not.valid_reactive_boundaries_3d(faces)) return
    if(present(wall_velocity)) then
      if(.not.valid_wall_velocity_3d(faces,wall_velocity)) return
    end if
    do a=1,3
      bounces=0
      do while (q(a)<1.or.q(a)>n(a))
        face=2*a
        if (q(a)<1) face=face-1
        select case(trim(faces(face)))
        case('periodic')
          q(a)=modulo(q(a)-1,n(a))+1
        case('inflow')
          inlet=face
          q(a)=max(1,min(n(a),q(a)))
        case('outflow')
          q(a)=max(1,min(n(a),q(a)))
        case('slip_wall','no_slip_wall','isothermal_wall')
          bounces=bounces+1
          if (bounces>size(wall_faces)) return
          wall_faces(bounces)=face
          if (q(a)<1) then
            q(a)=1-q(a)
          else
            q(a)=n(a)+1-(q(a)-n(a))
          end if
        end select
      end do
      ! Apply the encountered walls from the interior donor outward.
      ! Reversing this order matters for distinct moving wall speeds.
      do bounce=bounces,1,-1
        face=wall_faces(bounce)
        ! Euler reflects normal velocity only; diffusion imposes no-slip.
        if (faces(face)/='slip_wall'.and.viscous) then
          reflect=.not.reflect
          shift=-shift
          if (present(wall_velocity)) shift=shift+2.0_dp*wall_velocity(:,face)
        else
          reflect(a)=.not.reflect(a)
          shift(a)=-shift(a)
        end if
      end do
    end do
    ok=.true.
  end subroutine map_reactive_boundary_donor_3d

  pure logical function valid_wall_velocity_3d(faces,velocity) result(ok)
    character(len=*), intent(in) :: faces(6)
    real(dp), intent(in) :: velocity(:,:)
    integer, parameter :: axis(6)=[1,1,2,2,3,3]
    integer :: face
    ok=.false.
    if (any(shape(velocity)/=[3,6])) return
    if (.not.all(ieee_is_finite(velocity))) return
    do face=1,6
      if (abs(velocity(axis(face),face))>0.0_dp) return
      if (faces(face)/='no_slip_wall'.and.faces(face)/='isothermal_wall') then
        if (any(abs(velocity(:,face))>0.0_dp)) return
      end if
    end do
    ok=.true.
  end function valid_wall_velocity_3d

  pure function resolved_reactive_boundaries_3d(config) result(faces)
    type(reactive_3d_config), intent(in) :: config
    character(len=32) :: faces(6)
    integer :: f
    faces = config%boundary_faces
    do f = 1, 6
      if (len_trim(faces(f)) == 0) faces(f) = config%boundary_condition
    end do
  end function

  pure logical function valid_reactive_boundaries_3d(faces) result(ok)
    character(len=*), intent(in) :: faces(6)
    integer :: f
    ok = .false.
    do f = 1, 6
      select case(trim(faces(f)))
      case("periodic", "outflow", "slip_wall", "no_slip_wall", "inflow", "isothermal_wall")
      case default
        return
      end select
    end do
    do f = 1, 6, 2
      if ((faces(f) == "periodic") .neqv. (faces(f+1) == "periodic")) return
    end do
    ok = .true.
  end function

  subroutine read_reactive_3d_configuration( &
      path, config, ok, message, allow_selected, allow_pressure_outlet, allow_mach_inlet, allow_characteristic_ppm, allow_difmag)
    character(len=*), intent(in) :: path
    type(reactive_3d_config), intent(out) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected, allow_pressure_outlet, allow_mach_inlet
    logical, intent(in), optional :: allow_characteristic_ppm,allow_difmag

    integer :: nx, ny, nz, maximum_steps, unit, io_status
    integer :: wave_number_x, wave_number_y, wave_number_z
    real(dp) :: x_lower, x_upper, y_lower, y_upper, z_lower, z_upper
    real(dp) :: final_time, cfl, initial_temperature, initial_pressure,difmag
    real(dp) :: chemistry_relative_tolerance, chemistry_absolute_tolerance
    real(dp) :: transport_cfl
    real(dp) :: initial_velocity_x, initial_velocity_y, initial_velocity_z
    real(dp) :: density_wave_amplitude
    real(dp) :: hotspot_temperature_rise, hotspot_center_x
    real(dp) :: hotspot_center_y, hotspot_center_z, hotspot_width
    real(dp) :: x_h2, x_h, x_o, x_o2, x_oh, x_h2o
    real(dp) :: x_ho2, x_h2o2, x_ar, x_n2
    real(dp) :: composition_mole_fractions(reactive_3d_max_species)
    integer :: composition_count,inflow_composition_count
    real(dp) :: wall_temperature(6),wall_velocity(3,6),outlet_pressure(6),inflow_mach(6)
    real(dp) :: inflow_temperature(6),inflow_pressure(6),inflow_velocity(3,6)
    real(dp) :: inflow_composition_mole_fractions(reactive_3d_max_species)
    character(len=selected_composition_name_length) :: inflow_composition_species(reactive_3d_max_species)
    character(len=selected_composition_name_length) :: &
      composition_species(reactive_3d_max_species)
    character(len=32) :: problem, reconstruction, limiter
    character(len=32) :: thermo_model, riemann_solver
    character(len=32) :: boundary_condition, boundary_faces(6)
    character(len=512) :: output_file
    logical :: chemistry_enabled, transport_enabled, viscosity_enabled
    logical :: thermal_conduction_enabled, species_diffusion_enabled
    logical :: barodiffusion_enabled
    logical :: selected_allowed

    namelist /reactive_3d/ &
      nx, ny, nz, maximum_steps, &
      x_lower, x_upper, y_lower, y_upper, z_lower, z_upper, &
      final_time, cfl, problem, reconstruction, limiter, difmag, &
      thermo_model, riemann_solver, &
      boundary_condition, boundary_faces, wall_temperature, wall_velocity, inflow_temperature, inflow_pressure, inflow_velocity, &
      outlet_pressure, inflow_mach, inflow_composition_count,inflow_composition_species,inflow_composition_mole_fractions, &
      chemistry_enabled, chemistry_relative_tolerance, &
      chemistry_absolute_tolerance, transport_enabled, viscosity_enabled, &
      thermal_conduction_enabled, species_diffusion_enabled, &
      barodiffusion_enabled, transport_cfl, &
      initial_temperature, initial_pressure, &
      initial_velocity_x, initial_velocity_y, initial_velocity_z, &
      density_wave_amplitude, wave_number_x, wave_number_y, wave_number_z, &
      hotspot_temperature_rise, hotspot_center_x, hotspot_center_y, &
      hotspot_center_z, hotspot_width, &
      x_h2, x_h, x_o, x_o2, x_oh, x_h2o, x_ho2, x_h2o2, x_ar, x_n2, &
      composition_count, composition_species, composition_mole_fractions, &
      output_file

    config = reactive_3d_config()
    nx = config%nx
    ny = config%ny
    nz = config%nz
    maximum_steps = config%maximum_steps
    x_lower = config%x_lower
    x_upper = config%x_upper
    y_lower = config%y_lower
    y_upper = config%y_upper
    z_lower = config%z_lower
    z_upper = config%z_upper
    final_time = config%final_time
    cfl = config%cfl
    problem = config%problem
    difmag=config%difmag
    reconstruction = config%reconstruction
    limiter = config%limiter
    thermo_model = config%thermo_model
    riemann_solver = config%riemann_solver
    boundary_condition = config%boundary_condition
    inflow_temperature = config%inflow_temperature
    inflow_pressure = config%inflow_pressure
    inflow_velocity = config%inflow_velocity
    inflow_composition_count = config%inflow_composition_count
    inflow_composition_species = config%inflow_composition_species
    inflow_composition_mole_fractions = config%inflow_composition_mole_fractions
    wall_temperature=config%wall_temperature
    inflow_mach=config%inflow_mach
    outlet_pressure=config%outlet_pressure
    wall_velocity=config%wall_velocity
    boundary_faces = config%boundary_faces
    chemistry_enabled = config%chemistry_enabled
    chemistry_relative_tolerance = config%chemistry_relative_tolerance
    chemistry_absolute_tolerance = config%chemistry_absolute_tolerance
    transport_enabled = config%transport_enabled
    viscosity_enabled = config%viscosity_enabled
    thermal_conduction_enabled = config%thermal_conduction_enabled
    species_diffusion_enabled = config%species_diffusion_enabled
    barodiffusion_enabled = config%barodiffusion_enabled
    transport_cfl = config%transport_cfl
    initial_temperature = config%initial_temperature
    initial_pressure = config%initial_pressure
    initial_velocity_x = config%initial_velocity_x
    initial_velocity_y = config%initial_velocity_y
    initial_velocity_z = config%initial_velocity_z
    density_wave_amplitude = config%density_wave_amplitude
    wave_number_x = config%wave_number_x
    wave_number_y = config%wave_number_y
    wave_number_z = config%wave_number_z
    hotspot_temperature_rise = config%hotspot_temperature_rise
    hotspot_center_x = config%hotspot_center_x
    hotspot_center_y = config%hotspot_center_y
    hotspot_center_z = config%hotspot_center_z
    hotspot_width = config%hotspot_width
    x_h2 = config%x_h2
    x_h = config%x_h
    x_o = config%x_o
    x_o2 = config%x_o2
    x_oh = config%x_oh
    x_h2o = config%x_h2o
    x_ho2 = config%x_ho2
    x_h2o2 = config%x_h2o2
    x_ar = config%x_ar
    x_n2 = config%x_n2
    composition_count = config%composition_count
    composition_species = config%composition_species
    composition_mole_fractions = config%composition_mole_fractions
    output_file = config%output_file
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected

    open(newunit=unit, file=trim(path), status="old", action="read", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not open input file:", trim(path)
      return
    end if
    read(unit, nml=reactive_3d, iostat=io_status)
    close(unit)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') &
        "Could not read &reactive_3d from:", trim(path)
      return
    end if

    config%nx = nx
    config%ny = ny
    config%nz = nz
    config%maximum_steps = maximum_steps
    config%x_lower = x_lower
    config%x_upper = x_upper
    config%y_lower = y_lower
    config%y_upper = y_upper
    config%z_lower = z_lower
    config%z_upper = z_upper
    config%final_time = final_time
    config%cfl = cfl
    config%problem = trim(problem)
    config%difmag=difmag
    config%reconstruction = trim(reconstruction)
    config%limiter = trim(limiter)
    config%thermo_model = trim(thermo_model)
    config%riemann_solver = trim(riemann_solver)
    config%boundary_condition = trim(boundary_condition)
    config%inflow_temperature = inflow_temperature
    config%inflow_pressure = inflow_pressure
    config%inflow_velocity = inflow_velocity
    config%inflow_composition_count = inflow_composition_count
    config%inflow_composition_species = inflow_composition_species
    config%inflow_composition_mole_fractions = inflow_composition_mole_fractions
    config%wall_temperature=wall_temperature
    config%inflow_mach=inflow_mach
    config%outlet_pressure=outlet_pressure
    config%wall_velocity=wall_velocity
    config%boundary_faces = boundary_faces
    config%chemistry_enabled = chemistry_enabled
    config%chemistry_relative_tolerance = chemistry_relative_tolerance
    config%chemistry_absolute_tolerance = chemistry_absolute_tolerance
    config%transport_enabled = transport_enabled
    config%viscosity_enabled = viscosity_enabled
    config%thermal_conduction_enabled = thermal_conduction_enabled
    config%species_diffusion_enabled = species_diffusion_enabled
    config%barodiffusion_enabled = barodiffusion_enabled
    config%transport_cfl = transport_cfl
    config%initial_temperature = initial_temperature
    config%initial_pressure = initial_pressure
    config%initial_velocity_x = initial_velocity_x
    config%initial_velocity_y = initial_velocity_y
    config%initial_velocity_z = initial_velocity_z
    config%density_wave_amplitude = density_wave_amplitude
    config%wave_number_x = wave_number_x
    config%wave_number_y = wave_number_y
    config%wave_number_z = wave_number_z
    config%hotspot_temperature_rise = hotspot_temperature_rise
    config%hotspot_center_x = hotspot_center_x
    config%hotspot_center_y = hotspot_center_y
    config%hotspot_center_z = hotspot_center_z
    config%hotspot_width = hotspot_width
    config%x_h2 = x_h2
    config%x_h = x_h
    config%x_o = x_o
    config%x_o2 = x_o2
    config%x_oh = x_oh
    config%x_h2o = x_h2o
    config%x_ho2 = x_ho2
    config%x_h2o2 = x_h2o2
    config%x_ar = x_ar
    config%x_n2 = x_n2
    config%composition_count = composition_count
    config%composition_species = composition_species
    config%composition_mole_fractions = composition_mole_fractions
    config%output_file = trim(output_file)
    call validate_reactive_3d_configuration( &
      config, ok, message, allow_selected=selected_allowed)
    if (.not.ok) return
    if (config%difmag>=0) then
      ok=.false.
      if (present(allow_difmag)) ok=allow_difmag
      if (.not.ok) then
        message='Artificial-viscosity coefficient input is not connected in this application';return
      end if
    end if
    if (trim(config%reconstruction) == "characteristic_ppm") then
      ok = .false.
      if (present(allow_characteristic_ppm)) ok = allow_characteristic_ppm
      if (.not. ok) then
        message = 'Characteristic PPM input is not connected in this application'
        return
      end if
    end if
    if (any(config%inflow_mach>=0.0_dp)) then
      ok=.false.
      if (present(allow_mach_inlet)) ok=allow_mach_inlet
      if (.not.ok) then
        message='Mach inlet input is not connected in this application'; return
      end if
    end if
    if (any(config%outlet_pressure>0.0_dp)) then
      if (present(allow_pressure_outlet)) then
        if (allow_pressure_outlet) return
      end if
      ok=.false.
      message='Pressure outlet input is not connected in this application'
    end if
  end subroutine read_reactive_3d_configuration

  pure subroutine validate_reactive_3d_configuration( &
      config, ok, message, allow_selected)
    type(reactive_3d_config), intent(in) :: config
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    logical, intent(in), optional :: allow_selected

    real(dp) :: mole_sum
    logical :: selected_allowed, selected_model
    integer :: face
    character(len=32) :: faces(6)

    ok = .false.
    message = ""
    selected_allowed = .false.
    if (present(allow_selected)) selected_allowed = allow_selected
    selected_model = trim(config%thermo_model) == "selected"
    faces=resolved_reactive_boundaries_3d(config)
    if (.not.ieee_is_finite(config%difmag)) then
      message='Artificial-viscosity coefficient must be finite';return
    end if
    if (config%difmag<0.and.config%difmag/=-1) then
      message='Artificial-viscosity coefficient must be nonnegative or -1 for automatic';return
    end if
    if (config%difmag>=0.and.config%reconstruction/='characteristic_ppm') then
      message='Artificial-viscosity coefficient input requires characteristic_ppm';return
    end if
    if (.not.all(ieee_is_finite(config%inflow_mach))) then
      message='Inflow Mach numbers must be finite'; return
    end if
    if (any(config%inflow_mach>100).or.any(config%inflow_mach>=0.and.faces/='inflow')) then
      message='Active Mach numbers must be at most 100 and assigned to inflow faces'; return
    end if
    if (.not.all(ieee_is_finite(config%outlet_pressure))) then
      message='Outlet pressure must be finite'; return
    end if
    if (any(config%outlet_pressure<0.0_dp).or.any(config%outlet_pressure>1e15_dp)) then
      message='Outlet pressure must be between zero and 1e15 Pa'; return
    end if
    if (any(config%outlet_pressure>0.0_dp.and.faces/='outflow')) then
      message='Positive outlet pressure requires an outflow face'; return
    end if
    if (.not.valid_wall_velocity_3d(faces,config%wall_velocity)) then
      message='Wall velocity must be finite, tangential, and assigned only to no-slip walls'; return
    end if
    if (any(abs(config%wall_velocity)>0.0_dp)) then
      if (.not.config%transport_enabled.or..not.config%viscosity_enabled) then
        message='Moving walls require molecular viscosity'; return
      end if
    end if
    do face=1,6
      if (faces(face)=='isothermal_wall') then
        if (.not.ieee_is_finite(config%wall_temperature(face))) then
          message='Wall temperature must be finite'; return
        end if
        if (config%wall_temperature(face)<=0.0_dp) then
          message='Wall temperature must be positive'; return
        end if
        if (.not.config%transport_enabled.or..not.config%thermal_conduction_enabled) then
          message='Isothermal walls require molecular heat conduction'; return
        end if
      end if
      if (faces(face)/='inflow') cycle
      if (.not.all(ieee_is_finite([config%inflow_temperature(face),config%inflow_pressure(face), &
          config%inflow_velocity(:,face)]))) then
        message='Inflow thermodynamic state and velocity must be finite'; return
      end if
      if (config%inflow_mach(face)>=0) then
        if (config%inflow_temperature(face)<=0.or.config%inflow_pressure(face)<0.or. &
            config%inflow_pressure(face)>1e15_dp) then
          message='Mach inlet requires positive temperature and pressure between zero and 1e15 Pa'; return
        end if
        if (any(config%inflow_velocity(:,face)/=0)) then
          message='Mach inlet sets inward normal velocity; inflow_velocity must be zero'; return
        end if
      else if (min(config%inflow_temperature(face),config%inflow_pressure(face))<=0.0_dp) then
        message='Inflow temperature and pressure must be positive'; return
      end if
    end do
    if (any(faces=='inflow').and.config%inflow_composition_count/=0) then
      call validate_selected_composition_fields('Inflow',config%inflow_composition_count, &
        config%inflow_composition_species,config%inflow_composition_mole_fractions,ok,message)
      if (.not.ok) return
      ok=.false.
    end if
    if (config%nx < 4 .or. config%ny < 4 .or. config%nz < 4) then
      message = "Reactive 3D grid extents must all be at least 4"
    else if (config%maximum_steps <= 0) then
      message = "Reactive 3D maximum_steps must be positive"
    else if (config%x_upper <= config%x_lower .or. &
             config%y_upper <= config%y_lower .or. &
             config%z_upper <= config%z_lower) then
      message = "Reactive 3D domain maxima must exceed minima"
    else if (.not. all(ieee_is_finite([ &
        config%x_lower, config%x_upper, config%y_lower, config%y_upper, &
        config%z_lower, config%z_upper, config%final_time, config%cfl, &
        config%transport_cfl, &
        config%chemistry_relative_tolerance, &
        config%chemistry_absolute_tolerance, &
        config%initial_temperature, config%initial_pressure, &
        config%initial_velocity_x, config%initial_velocity_y, &
        config%initial_velocity_z, config%density_wave_amplitude, &
        config%hotspot_temperature_rise, config%hotspot_center_x, &
        config%hotspot_center_y, config%hotspot_center_z, &
        config%hotspot_width]))) then
      message = "Reactive 3D scalar inputs must be finite"
    else if (config%final_time <= 0.0_dp) then
      message = "Reactive 3D final_time must be positive"
    else if (config%cfl <= 0.0_dp .or. config%cfl > 1.0_dp) then
      message = "Reactive 3D cfl must be in (0, 1]"
    else if (config%transport_cfl <= 0.0_dp .or. &
             config%transport_cfl > 0.5_dp) then
      message = "Reactive 3D transport_cfl must be in (0, 0.5]"
    else if (config%transport_enabled .and. &
             .not. (config%viscosity_enabled .or. &
               config%thermal_conduction_enabled .or. &
               config%species_diffusion_enabled)) then
      message = "Reactive 3D transport requires an enabled process"
    else if (config%transport_enabled .and. &
             config%barodiffusion_enabled .and. &
             .not. config%species_diffusion_enabled) then
      message = "Reactive 3D barodiffusion requires species diffusion"
    else if (config%chemistry_relative_tolerance <= 0.0_dp .or. &
             config%chemistry_absolute_tolerance <= 0.0_dp) then
      message = "Reactive 3D chemistry tolerances must be positive"
    else if (config%initial_temperature <= 0.0_dp .or. &
             config%initial_pressure <= 0.0_dp) then
      message = "Reactive 3D initial temperature and pressure must be positive"
    else if (trim(config%problem) == "entropy_wave" .and. &
             abs(config%density_wave_amplitude) >= 1.0_dp) then
      message = "Reactive 3D relative density amplitude must be below one"
    else if (trim(config%problem) == "entropy_wave" .and. &
             config%wave_number_x == 0 .and. &
             config%wave_number_y == 0 .and. &
             config%wave_number_z == 0) then
      message = "Reactive 3D requires at least one nonzero wave number"
    else if (trim(config%problem) /= "entropy_wave" .and. &
             trim(config%problem) /= "uniform_reactor" .and. &
             trim(config%problem) /= "reactive_hotspot") then
      message = "Unknown reactive 3D problem"
    else if (trim(config%problem) == "reactive_hotspot" .and. &
             (config%hotspot_temperature_rise < 0.0_dp .or. &
              config%hotspot_width <= 0.0_dp .or. &
              config%hotspot_center_x < config%x_lower .or. &
              config%hotspot_center_x > config%x_upper .or. &
              config%hotspot_center_y < config%y_lower .or. &
              config%hotspot_center_y > config%y_upper .or. &
              config%hotspot_center_z < config%z_lower .or. &
              config%hotspot_center_z > config%z_upper)) then
      message = "Invalid reactive 3D hotspot geometry"
    else if (trim(config%thermo_model) /= "elementary" .and. &
             trim(config%thermo_model) /= "full_h2o2" .and. &
             .not. (selected_model .and. selected_allowed)) then
      message = &
        "Reactive 3D thermo_model must be elementary, full_h2o2, or selected"
    else if (trim(config%reconstruction) /= "pcm" .and. &
             trim(config%reconstruction) /= "characteristic_plm" .and. &
             trim(config%reconstruction) /= "characteristic_ppm") then
      message = "Reactive 3D reconstruction must be pcm, characteristic_plm, or characteristic_ppm"
    else if (trim(config%reconstruction) == "characteristic_plm" .and. &
             min(config%nx, config%ny, config%nz) < 3) then
      message = "Reactive 3D characteristic_plm requires nx, ny, nz >= 3"
    else if (trim(config%limiter) /= "minmod" .and. &
             trim(config%limiter) /= "mc") then
      message = "Reactive 3D limiter must be minmod or mc"
    else if (trim(config%riemann_solver) /= "rusanov" .and. &
             trim(config%riemann_solver) /= "hllc" .and. &
             trim(config%riemann_solver) /= "pelec") then
      message = "Reactive 3D Riemann solver must be rusanov, hllc, or pelec"
    else if (.not. valid_reactive_boundaries_3d(resolved_reactive_boundaries_3d(config))) then
      message = "Unknown reactive 3D face or unpaired periodic boundary"
    else if (len_trim(config%output_file) == 0) then
      message = "Reactive 3D output_file must not be empty"
    else
      if (selected_model) then
        call validate_selected_composition_fields( &
          "Reactive 3D", config%composition_count, &
          config%composition_species, config%composition_mole_fractions, &
          ok, message)
        return
      end if
      if (.not. all(ieee_is_finite([ &
          config%x_h2, config%x_h, config%x_o, config%x_o2, config%x_oh, &
          config%x_h2o, config%x_ho2, config%x_h2o2, config%x_ar, &
          config%x_n2]))) then
        message = "Reactive 3D mole fractions must be finite"
        return
      end if
      if (min(config%x_h2, config%x_h, config%x_o, config%x_o2, &
              config%x_oh, config%x_h2o, config%x_ho2, &
              config%x_h2o2, config%x_ar, config%x_n2) < 0.0_dp) then
        message = "Reactive 3D mole fractions must be nonnegative"
        return
      end if
      mole_sum = config%x_h2 + config%x_h + config%x_o + config%x_o2 + &
        config%x_oh + config%x_h2o + config%x_ho2 + config%x_h2o2 + &
        config%x_ar + config%x_n2
      if (abs(mole_sum - 1.0_dp) > 5.0e-10_dp) then
        message = "Reactive 3D mole fractions must sum to one"
      else if (trim(config%thermo_model) == "elementary" .and. &
               max(config%x_ho2, config%x_h2o2, config%x_ar) > &
                 5.0e-14_dp) then
        message = "Elementary thermo requires zero HO2/H2O2/AR fractions"
      else
        ok = .true.
      end if
    end if
  end subroutine validate_reactive_3d_configuration

  subroutine reactive_3d_mole_fractions(config, nspecies, mole_fractions, ok)
    type(reactive_3d_config), intent(in) :: config
    integer, intent(in) :: nspecies
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok

    mole_fractions = 0.0_dp
    ok = size(mole_fractions) == nspecies
    if (.not. ok) return
    select case (trim(config%thermo_model))
    case ("elementary")
      if (nspecies /= 7) then
        ok = .false.
        return
      end if
      mole_fractions = [ &
        config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_n2]
    case ("full_h2o2")
      if (nspecies /= 10) then
        ok = .false.
        return
      end if
      mole_fractions = [ &
        config%x_h2, config%x_h, config%x_o, config%x_o2, &
        config%x_oh, config%x_h2o, config%x_ho2, config%x_h2o2, &
        config%x_ar, config%x_n2]
    case default
      ok = .false.
      return
    end select
    ok = minval(mole_fractions) >= 0.0_dp .and. &
      abs(sum(mole_fractions) - 1.0_dp) <= 5.0e-10_dp
  end subroutine reactive_3d_mole_fractions


  subroutine resolve_reactive_3d_selected_composition( &
      config, species, mole_fractions, ok, message)
    type(reactive_3d_config), intent(in) :: config
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(out) :: mole_fractions(:)
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    mole_fractions = 0.0_dp
    if (trim(config%thermo_model) /= "selected") then
      ok = .false.
      message = &
        "Reactive 3D selected composition requires thermo_model='selected'"
      return
    end if
    call resolve_selected_composition( &
      "Reactive 3D", config%composition_count, config%composition_species, &
      config%composition_mole_fractions, species, mole_fractions, ok, message)
  end subroutine resolve_reactive_3d_selected_composition

end module simulation_config_reactive_3d_mod
