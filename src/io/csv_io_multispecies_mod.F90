module csv_io_multispecies_mod
  use precision_mod, only: dp
  use state_indices_mod, only: irho, imx, iet, iei, item, ncons, nprim, qrho, qu, qp
  use state_conversion_mod, only: conserved_to_primitive
  use multispecies_state_mod, only: &
    max_supported_species, multispecies_nvar, species_component, &
    mass_fractions_from_state
  use eos_ideal_mod, only: ideal_gas_specific_internal_energy, &
    ideal_gas_sound_speed
  implicit none
  private

  public :: write_multispecies_solution_csv

contains

  subroutine write_multispecies_solution_csv( &
      path, x, conserved, nx, nspecies, gamma, ok, message)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nx, nspecies
    real(dp), intent(in) :: x(nx), conserved(:, 0:), gamma
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: unit, io_status, i, species, nvar
    real(dp) :: primitive(nprim), internal_energy, sound_speed
    real(dp) :: mass_fractions(max_supported_species)
    logical :: cell_ok

    ok = .false.
    message = ""
    nvar = multispecies_nvar(nspecies)
    if (nvar == 0 .or. size(conserved, 1) /= nvar) then
      message = "Invalid multispecies state layout"
      return
    end if

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=io_status)
    if (io_status /= 0) then
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if

    write(unit, '(a)', advance="no") &
      "x,rho,u,p,specific_internal_energy,internal_energy_density," // &
      "temperature,total_energy_density,momentum_density,sound_speed"
    do species = 1, nspecies
      write(unit, '(a,i0,a,i0)', advance="no") &
        ",rhoY", species, ",Y", species
    end do
    write(unit, '(a)') ""

    do i = 1, nx
      call conserved_to_primitive( &
        conserved(1:ncons, i), gamma, primitive, cell_ok)
      if (.not. cell_ok) then
        close(unit)
        write(message, '(a,i0)') "Non-physical state while writing cell ", i
        return
      end if
      call mass_fractions_from_state( &
        conserved(:, i), nspecies, mass_fractions, cell_ok)
      if (.not. cell_ok) then
        close(unit)
        write(message, '(a,i0)') &
          "Invalid species state while writing cell ", i
        return
      end if

      internal_energy = ideal_gas_specific_internal_energy( &
        primitive(qrho), primitive(qp), gamma)
      sound_speed = ideal_gas_sound_speed( &
        primitive(qrho), primitive(qp), gamma)

      write(unit, '(*(es25.16e3,:,","))') &
        x(i), primitive(qrho), primitive(qu), primitive(qp), &
        internal_energy, conserved(iei, i), conserved(item, i), &
        conserved(iet, i), conserved(imx, i), sound_speed, &
        (conserved(species_component(species), i), &
         mass_fractions(species), species = 1, nspecies)
    end do

    close(unit)
    ok = .true.
  end subroutine write_multispecies_solution_csv

end module csv_io_multispecies_mod
