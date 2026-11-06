module reactive_eb_cfl_3d_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use precision_mod, only: dp
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_conserved_to_primitive
  use eb_geometry_3d_mod, only: eb_geometry_3d, eb_covered_cell_3d
  implicit none
  private

  public :: compute_reactive_eb_cfl_timestep_3d
  public :: compute_reactive_eb_euler_cfl_timestep_3d,reactive_eb_signal_speeds_3d

contains

  subroutine compute_reactive_eb_cfl_timestep_3d( &
      species, state, temperature, geometry, cfl, dt, ok,internal_energy_density)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(in) :: cfl
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok

    real(dp), allocatable :: primitive(:),held
    real(dp) :: local_temperature, sound_speed, rate, maximum_rate
    logical :: local_ok
    integer :: i, j, k, active_cells, nvar

    dt = 0.0_dp
    ok = .false.
    nvar = reactive_nvar(size(species))
    if (nvar <= 0 .or. .not. geometry%is_valid() .or. &
        .not. ieee_is_finite(cfl) .or. cfl <= 0.0_dp .or. &
        cfl > 1.0_dp) return
    if (size(state, 1) /= nvar .or. &
        size(state, 2) /= geometry%nx .or. &
        size(state, 3) /= geometry%ny .or. &
        size(state, 4) /= geometry%nz .or. &
        any(shape(temperature) /= &
          [geometry%nx, geometry%ny, geometry%nz])) return

    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=shape(temperature))) return
      allocate(held)
    endif
    allocate(primitive(reactive_nprim(size(species))))
    maximum_rate = 0.0_dp
    active_cells = 0
    do k = 1, geometry%nz
      do j = 1, geometry%ny
        do i = 1, geometry%nx
          if (geometry%cell_type(i, j, k) == eb_covered_cell_3d) cycle
          active_cells = active_cells + 1
          if(present(internal_energy_density)) held=internal_energy_density(i,j,k)
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), &
            primitive, local_temperature, sound_speed, local_ok,internal_energy_density=held)
          if (.not. local_ok) return
          rate = (abs(primitive(2)) + sound_speed) / geometry%dx + &
            (abs(primitive(3)) + sound_speed) / geometry%dy + &
            (abs(primitive(4)) + sound_speed) / geometry%dz
          maximum_rate = max(maximum_rate, rate)
        end do
      end do
    end do
    if (active_cells == 0 .or. .not. ieee_is_finite(maximum_rate) .or. &
        maximum_rate <= 0.0_dp) return

    dt = cfl / maximum_rate
    ok = ieee_is_finite(dt) .and. dt > 0.0_dp
    if (.not. ok) dt = 0.0_dp
  end subroutine compute_reactive_eb_cfl_timestep_3d

  ! Maximum Cartesian signal speeds plus an orientation-independent wall bound.
  ! A hierarchy can combine these four values across patches/ranks and pass the
  ! result as donor_speeds to include remote face donors in the timestep rate.
  subroutine reactive_eb_signal_speeds_3d(species,state,temperature,geometry,speeds,ok,internal_energy_density)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:)
    real(dp), optional, intent(in) :: internal_energy_density(:,:,:)
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: speeds(4)
    logical, intent(out) :: ok
    real(dp) :: primitive(reactive_nprim(size(species))),local_temperature,sound_speed,candidate(4)
    real(dp), allocatable :: held
    integer :: i,j,k,active
    logical :: valid
    speeds=0._dp;ok=.false.;candidate=0._dp;active=0
    if (.not.geometry%is_valid()) return
    if (size(state,1)/=reactive_nvar(size(species))) return
    if (any(shape(state(1,:,:,:))/=[geometry%nx,geometry%ny,geometry%nz])) return
    if (any(shape(temperature)/=[geometry%nx,geometry%ny,geometry%nz])) return
    if(present(internal_energy_density)) then
      if(any(shape(internal_energy_density)/=shape(temperature))) return
      allocate(held)
    endif
    do k=1,geometry%nz
      do j=1,geometry%ny
        do i=1,geometry%nx
          if (geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
          if(present(internal_energy_density)) held=internal_energy_density(i,j,k)
          call reactive_conserved_to_primitive(species,state(:,i,j,k),temperature(i,j,k), &
            primitive,local_temperature,sound_speed,valid,internal_energy_density=held)
          if (.not.valid) return
          candidate(1:3)=max(candidate(1:3),abs(primitive(2:4))+sound_speed)
          candidate(4)=max(candidate(4),sum(abs(primitive(2:4)))+sound_speed)
          active=active+1
        end do
      end do
    end do
    if (active==0.or..not.all(ieee_is_finite(candidate))) return
    if (any(candidate<=0._dp)) return
    speeds=candidate;ok=.true.
  end subroutine reactive_eb_signal_speeds_3d

  ! Plain Euler small-cell estimate: dt = CFL / max(sum(A * signal)/(kappa*V)).
  ! Both Cartesian faces and the embedded wall enter the rate. Using full
  ! face sums is deliberately twice the usual Cartesian acoustic rate for a
  ! uniform regular grid. This is a hydro estimate, not a positivity guarantee
  ! for nonlinear Riemann solvers, reflux, chemistry or molecular transport.
  subroutine compute_reactive_eb_euler_cfl_timestep_3d( &
      species,state,temperature,geometry,cfl,dt,ok,donor_speeds,internal_energy_density)
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: state(:,:,:,:),temperature(:,:,:),cfl
    type(eb_geometry_3d), intent(in) :: geometry
    real(dp), intent(out) :: dt
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: donor_speeds(4)
    real(dp), intent(in), optional :: internal_energy_density(:,:,:)
    real(dp) :: speeds(4),rate,maximum_rate,volume
    integer :: i,j,k
    logical :: valid
    dt=0._dp;ok=.false.
    if (.not.ieee_is_finite(cfl)) return
    if (cfl<=0._dp.or.cfl>1._dp) return
    call reactive_eb_signal_speeds_3d(species,state,temperature,geometry,speeds,valid,internal_energy_density)
    if (.not.valid) return
    if (present(donor_speeds)) then
      if (.not.all(ieee_is_finite(donor_speeds))) return
      if (any(donor_speeds<0._dp)) return
      speeds=max(speeds,donor_speeds)
    end if
    volume=geometry%dx*geometry%dy*geometry%dz
    maximum_rate=0._dp
    do k=1,geometry%nz
      do j=1,geometry%ny
        do i=1,geometry%nx
          if (geometry%cell_type(i,j,k)==eb_covered_cell_3d) cycle
          rate=(speeds(1)*(geometry%x_face_fraction(i-1,j,k)+geometry%x_face_fraction(i,j,k))/geometry%dx &
            +speeds(2)*(geometry%y_face_fraction(i,j-1,k)+geometry%y_face_fraction(i,j,k))/geometry%dy &
            +speeds(3)*(geometry%z_face_fraction(i,j,k-1)+geometry%z_face_fraction(i,j,k))/geometry%dz &
            +speeds(4)*geometry%boundary_area(i,j,k)/volume)/geometry%volume_fraction(i,j,k)
          if (.not.ieee_is_finite(rate).or.rate<=0._dp) return
          maximum_rate=max(maximum_rate,rate)
        end do
      end do
    end do
    if (maximum_rate<=0._dp) return
    dt=cfl/maximum_rate
    ok=ieee_is_finite(dt).and.dt>0._dp
    if (.not.ok) dt=0._dp
  end subroutine compute_reactive_eb_euler_cfl_timestep_3d

end module reactive_eb_cfl_3d_mod
