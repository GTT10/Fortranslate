program test_reactive_eb_3d_config
  use, intrinsic :: ieee_arithmetic, only: &
    ieee_value, ieee_quiet_nan
  use precision_mod, only: dp
  use simulation_config_reactive_eb_3d_mod, only: &
    reactive_eb_3d_config, read_reactive_eb_3d_configuration, &
    validate_reactive_eb_3d_configuration,build_configured_eb_geometry_3d
  use simulation_config_reactive_3d_mod, only: map_reactive_boundary_donor_3d
  use eb_geometry_3d_mod, only: eb_geometry_3d,extract_eb_geometry_3d
  implicit none

  character(len=*), parameter :: path = "reactive_eb_3d_config_test.nml"
  type(reactive_eb_3d_config) :: config
  character(len=256) :: message
  logical :: ok
  integer :: unit

  call verify_extended_geometry()
  config=reactive_eb_3d_config()
  config%boundary_faces='outflow';config%boundary_faces(3)='inflow'
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(.not.ok,'physical boundary requires CTU')
  config%hydro_method='plm_ctu';config%riemann_solver='pelec'
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(ok,'physical inflow and outflow config')
  config%boundary_faces(1)='legacy'
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(.not.ok,'mixed legacy and physical boundaries reject')
  config%boundary_faces(1)='outflow';config%inflow_density(3)=-1
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(.not.ok,'negative inlet density rejects')
  config%inflow_density=.31_dp;config%boundary_faces='isothermal_wall'
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(.not.ok,'isothermal boundary requires transport')
  config%transport_enabled=.true.
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(ok,'isothermal boundary input')
  config%wall_temperature(1)=-1
  call validate_reactive_eb_3d_configuration(config,ok,message)
  call require(.not.ok,'invalid wall temperature rejects')


  config = reactive_eb_3d_config()
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "default public 3D EB configuration")
  config%chemistry_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "reacting public 3D EB configuration")

  config%plane_axis = "y"
  config%plane_position = 0.49375_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "y-normal public 3D EB configuration")
  config%plane_axis = "z"
  config%plane_position = 0.4916666666666667_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "z-normal public 3D EB configuration")

  config = reactive_eb_3d_config()
  config%state_redist_target_volume_fraction = 0.0_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "invalid public StateRedist target")
  config = reactive_eb_3d_config()
  config%state_redist_target_volume_fraction = 0.05_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "inactive public StateRedist target")
  config%redistribution = "flux_redist"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "FluxRedist ignores StateRedist target activation")
  config = reactive_eb_3d_config()
  config%redistribution = "raw"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "unstabilized public route rejection")
  config = reactive_eb_3d_config()
  config%plane_position = 0.4_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "face-aligned public plane rejection")
  config = reactive_eb_3d_config()
  config%plane_position = 0.95_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "receiver-free public plane rejection")
  config = reactive_eb_3d_config()
  config%initial_pressure = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonfinite public configuration rejection")
  config = reactive_eb_3d_config()
  config%chemistry_relative_tolerance = 0.0_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonpositive chemistry tolerance rejection")
  config = reactive_eb_3d_config()
  config%chemistry_absolute_tolerance = &
    ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonfinite chemistry tolerance rejection")
  config = reactive_eb_3d_config()
  config%transport_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "transport public 3D EB configuration")
  config%redistribution = "flux_redist"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "transport requires StateRedist")
  config = reactive_eb_3d_config()
  config%transport_enabled = .true.
  config%viscosity_enabled = .false.
  config%thermal_conduction_enabled = .false.
  config%species_diffusion_enabled = .false.
  config%barodiffusion_enabled = .false.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "transport process rejection")
  config%thermal_conduction_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "thermal-only transport configuration")
  config%barodiffusion_enabled = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "barodiffusion requires species diffusion")
  config = reactive_eb_3d_config()
  config%transport_cfl = 0.0_dp
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonpositive transport CFL rejection")
  config%transport_cfl = ieee_value(0.0_dp, ieee_quiet_nan)
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "nonfinite transport CFL rejection")
  config = reactive_eb_3d_config()
  config%checkpoint_interval_steps = -1
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "negative checkpoint interval rejection")
  config = reactive_eb_3d_config()
  config%stop_after_checkpoint = .true.
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "checkpoint stop without interval rejection")
  config%checkpoint_interval_steps = 1
  config%checkpoint_file = ""
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "empty checkpoint path rejection")
  config%checkpoint_file = config%output_file
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "checkpoint output collision rejection")
  config%checkpoint_file = "./reactive_eb_3d.csv"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "lexical checkpoint output alias rejection")
  config%checkpoint_file = "scratch/../reactive_eb_3d.csv"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "parent checkpoint output alias rejection")
  config = reactive_eb_3d_config()
  config%restart_file = config%output_file
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "restart output collision rejection")
  config%restart_file = "./reactive_eb_3d.csv"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(.not. ok, "lexical restart output alias rejection")
  config = reactive_eb_3d_config()
  config%checkpoint_interval_steps = 2
  config%checkpoint_file = "resume.chk"
  config%restart_file = "resume.chk"
  call validate_reactive_eb_3d_configuration(config, ok, message)
  call require(ok, "same checkpoint and restart path")

  open(newunit=unit, file=path, status="replace", action="write")
  write(unit, '(a)') "&reactive_eb_3d"
  write(unit, '(a)') " nx = 12"
  write(unit, '(a)') " plane_position = 0.3291666666666667"
  write(unit, '(a)') " redistribution = 'state_redist'"
  write(unit, '(a)') " chemistry_enabled = .true."
  write(unit, '(a)') " chemistry_relative_tolerance = 1.0e-8"
  write(unit, '(a)') " chemistry_absolute_tolerance = 1.0e-13"
  write(unit, '(a)') " transport_enabled = .true."
  write(unit, '(a)') " viscosity_enabled = .false."
  write(unit, '(a)') " thermal_conduction_enabled = .true."
  write(unit, '(a)') " species_diffusion_enabled = .false."
  write(unit, '(a)') " barodiffusion_enabled = .false."
  write(unit, '(a)') " transport_cfl = 0.2"
  write(unit, '(a)') " checkpoint_interval_steps = 3"
  write(unit, '(a)') " stop_after_checkpoint = .true."
  write(unit, '(a)') " checkpoint_file = 'configured_eb_3d.chk'"
  write(unit, '(a)') " restart_file = 'restart_eb_3d.chk'"
  write(unit, '(a)') " output_file = 'configured_eb_3d.csv'"
  write(unit, '(a)') "/"
  close(unit)
  call read_reactive_eb_3d_configuration(path, config, ok, message)
  call require(ok .and. config%difmag==.1_dp .and. config%nx == 12 .and. &
    trim(config%redistribution) == "state_redist" .and. &
    config%chemistry_enabled .and. &
    config%chemistry_relative_tolerance == 1.0e-8_dp .and. &
    config%chemistry_absolute_tolerance == 1.0e-13_dp .and. &
    config%transport_enabled .and. .not. config%viscosity_enabled .and. &
    config%thermal_conduction_enabled .and. &
    .not. config%species_diffusion_enabled .and. &
    .not. config%barodiffusion_enabled .and. &
    config%transport_cfl == 0.2_dp .and. &
    config%checkpoint_interval_steps == 3 .and. &
    config%stop_after_checkpoint .and. &
    trim(config%checkpoint_file) == "configured_eb_3d.chk" .and. &
    trim(config%restart_file) == "restart_eb_3d.chk" .and. &
    trim(config%output_file) == "configured_eb_3d.csv", &
    "public 3D EB namelist read")
  open(newunit=unit, file=path, status="old")
  close(unit, status="delete")

  write(*, '(a)') "test_reactive_eb_3d_config: PASS"

contains
  subroutine verify_extended_geometry()
    type(reactive_eb_3d_config) :: cfg
    type(eb_geometry_3d) :: base,wide,cropped,saved
    real(dp), allocatable :: nodes(:,:,:),expanded(:,:,:)
    real(dp) :: section(0:4,0:4),field(0:4,0:4,0:4),error,shift(3)
    logical :: reflect(3)
    integer :: boundary_mode,donor(3),inlet
    integer :: i,j,k,axis,kind,handle,point(3)
    character(len=1), parameter :: axes(3)=['x','y','z']
    logical :: valid
    cfg=reactive_eb_3d_config();cfg%nx=4;cfg%ny=4;cfg%nz=4
    cfg%plane_position=.395_dp
    cfg%section_file='extended_eb_section.dat';cfg%level_set_file='extended_eb_volume.dat'
    do k=0,4
    do j=0,4
    do i=0,4
      field(i,j,k)=real(i,dp)+.5_dp*j+.25_dp*k-2.3_dp
      section(i,j)=real(i,dp)+.5_dp*j-2.3_dp
    end do
    end do
    end do
    open(newunit=handle,file=trim(cfg%section_file),status='replace')
    write(handle,*) 4,4
    write(handle,*) section
    close(handle)
    open(newunit=handle,file=trim(cfg%level_set_file),status='replace')
    write(handle,*) 4,4,4
    write(handle,*) field
    close(handle)
    do boundary_mode=0,2
    cfg%boundary_faces='legacy'
    if(boundary_mode==1) cfg%boundary_faces='slip_wall'
    if(boundary_mode==2) then
      cfg%boundary_faces='outflow';cfg%boundary_faces(1:6:2)='slip_wall'
    end if
    do kind=1,3
      select case(kind)
      case(1);cfg%geometry_type='plane'
      case(2);cfg%geometry_type='extruded'
      case(3);cfg%geometry_type='level_set'
      end select
      do axis=1,3
        if(kind==3.and.axis>1) cycle
        cfg%plane_axis=axes(axis);cfg%extrusion_axis=axes(axis)
        if(kind==3) then
          call build_configured_eb_geometry_3d(cfg,base,valid,nodal_field=nodes)
          call require(valid,'base nodal geometry')
          call build_configured_eb_geometry_3d(cfg,wide,valid,nodal_field=expanded,halo=6)
        else
          call build_configured_eb_geometry_3d(cfg,base,valid)
          call require(valid,'base plane/extrusion geometry')
          call build_configured_eb_geometry_3d(cfg,wide,valid,halo=6)
        end if
        call require(valid,'six-cell physical geometry extension')
        call require(all([wide%nx,wide%ny,wide%nz]==16),'extended dimensions')
        call extract_eb_geometry_3d(wide,[7,7,7],[10,10,10],cropped,valid)
        call require(valid,'physical geometry recovery')
        error=max(maxval(abs(cropped%volume_fraction-base%volume_fraction)), &
          maxval(abs(cropped%x_face_fraction-base%x_face_fraction)), &
          maxval(abs(cropped%y_face_fraction-base%y_face_fraction)), &
          maxval(abs(cropped%z_face_fraction-base%z_face_fraction)), &
          maxval(abs(cropped%cell_centroid_x-base%cell_centroid_x)), &
          maxval(abs(cropped%cell_centroid_y-base%cell_centroid_y)), &
          maxval(abs(cropped%cell_centroid_z-base%cell_centroid_z)), &
          maxval(abs(cropped%boundary_area-base%boundary_area)))
        call require(error<2e-12_dp.and.all(cropped%cell_type==base%cell_type), &
          'extension preserves physical cut geometry')
        if(kind==3) then
          do k=0,16
          do j=0,16
          do i=0,16
            select case(boundary_mode)
            case(0);point=max(0,min(4,[i,j,k]-6))
            case(1)
              point=modulo([i,j,k]-6,8)
              point=min(point,8-point)
            case(2);point=min(4,abs([i,j,k]-6))
            end select
            call require(expanded(i,j,k)==nodes(point(1),point(2),point(3)), &
              'domain-face nodal extension includes edges and corners')
          end do
          end do
          end do
        end if
        if(boundary_mode==1) then
          do k=1,16
          do j=1,16
          do i=1,16
            call map_reactive_boundary_donor_3d([i,j,k]-6,[4,4,4],cfg%boundary_faces,.false., &
              donor,reflect,shift,inlet,valid)
            call require(valid,'wall state donor mapping')
            call require(abs(wide%volume_fraction(i,j,k)-base%volume_fraction(donor(1),donor(2),donor(3)))<2e-12_dp, &
              'reflected geometry matches repeatedly reflected state donors')
          end do
          end do
          end do
        end if
      end do
    end do
    end do
    saved=wide
    call build_configured_eb_geometry_3d(cfg,wide,valid,halo=-1)
    call require(.not.valid.and.all(wide%volume_fraction==saved%volume_fraction),'negative halo rollback')
    open(newunit=handle,file=trim(cfg%section_file),status='old');close(handle,status='delete')
    open(newunit=handle,file=trim(cfg%level_set_file),status='old');close(handle,status='delete')
  end subroutine verify_extended_geometry


  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label

    if (.not. condition) then
      write(*, '(a,1x,a)') trim(label), trim(message)
      error stop label
    end if
  end subroutine require

end program test_reactive_eb_3d_config
