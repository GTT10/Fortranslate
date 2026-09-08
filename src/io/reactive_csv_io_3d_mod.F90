module reactive_csv_io_3d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: imx, imy, imz, iet
  use nasa7_thermo_mod, only: nasa7_species
  use reactive_1d_mod, only: &
    reactive_nvar, reactive_nprim, reactive_species_component, &
    reactive_mass_fraction_component, reactive_conserved_to_primitive
  implicit none
  private

  public :: write_reactive_3d_csv

contains

  subroutine write_reactive_3d_csv( &
      path, species, x, y, z, state, temperature, nx, ny, nz, time, &
      ok, message)
    character(len=*), intent(in) :: path
    type(nasa7_species), intent(in) :: species(:)
    real(dp), intent(in) :: x(:), y(:), z(:)
    real(dp), intent(in) :: state(:, :, :, :), temperature(:, :, :)
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: time
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    real(dp), allocatable :: primitive(:)
    real(dp) :: local_temperature, sound_speed
    logical :: cell_ok
    integer :: unit, io_status, i, j, k, species_index, nspecies

    ok = .false.
    message = ""
    nspecies = size(species)
    if (size(x) /= nx .or. size(y) /= ny .or. size(z) /= nz .or. &
        size(state, 1) /= reactive_nvar(nspecies) .or. &
        size(state, 2) /= nx .or. size(state, 3) /= ny .or. &
        size(state, 4) /= nz .or. size(temperature, 1) /= nx .or. &
        size(temperature, 2) /= ny .or. size(temperature, 3) /= nz) then
      message = "Reactive 3D CSV array shapes are inconsistent"
      return
    end if

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=io_status)
    if (io_status /= 0) then
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if

    write(unit, '(a)', advance='no') &
      "time,x,y,z,rho,u,v,w,pressure,temperature,rhoE,rhou,rhov,rhow"
    do species_index = 1, nspecies
      write(unit, '(a)', advance='no') &
        ",Y_" // trim(species(species_index)%name)
    end do
    do species_index = 1, nspecies
      write(unit, '(a)', advance='no') &
        ",rhoY_" // trim(species(species_index)%name)
    end do
    write(unit, '(a)') ""

    allocate(primitive(reactive_nprim(nspecies)))
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call reactive_conserved_to_primitive( &
            species, state(:, i, j, k), temperature(i, j, k), primitive, &
            local_temperature, sound_speed, cell_ok)
          if (.not. cell_ok) then
            close(unit)
            write(message, '(a,i0,a,i0,a,i0,a)') &
              "Non-physical reactive state while writing cell (", i, ",", &
              j, ",", k, ")"
            return
          end if
          write(unit, '(*(es25.16e3,:,","))') &
            time, x(i), y(j), z(k), primitive(1), primitive(2), &
            primitive(3), primitive(4), primitive(5), local_temperature, &
            state(iet, i, j, k), state(imx, i, j, k), &
            state(imy, i, j, k), state(imz, i, j, k), &
            (primitive(reactive_mass_fraction_component(species_index)), &
              species_index = 1, nspecies), &
            (state(reactive_species_component(species_index), i, j, k), &
              species_index = 1, nspecies)
        end do
      end do
    end do
    close(unit)
    ok = .true.
  end subroutine write_reactive_3d_csv

end module reactive_csv_io_3d_mod
