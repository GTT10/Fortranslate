module csv_io_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    irho, imx, iet, ncons, nprim, qrho, qu, qp
  use state_conversion_mod, only: conserved_to_primitive
  use eos_ideal_mod, only: ideal_gas_specific_internal_energy, ideal_gas_sound_speed
  implicit none
  private

  public :: write_solution_csv

contains

  subroutine write_solution_csv(path, x, conserved, nx, gamma, ok, message)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nx
    real(dp), intent(in) :: x(nx)
    real(dp), intent(in) :: conserved(ncons, 0:nx + 1)
    real(dp), intent(in) :: gamma
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: unit, io_status, i
    real(dp) :: primitive(nprim), internal_energy, sound_speed
    logical :: cell_ok

    open(newunit=unit, file=trim(path), status="replace", action="write", iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if

    write(unit, '(a)') &
      "x,rho,u,p,specific_internal_energy,total_energy_density,momentum_density,sound_speed"

    do i = 1, nx
      call conserved_to_primitive(conserved(:, i), gamma, primitive, cell_ok)
      if (.not. cell_ok) then
        close(unit)
        ok = .false.
        write(message, '(a,i0)') "Non-physical state while writing cell ", i
        return
      end if

      internal_energy = ideal_gas_specific_internal_energy( &
        primitive(qrho), primitive(qp), gamma)
      sound_speed = ideal_gas_sound_speed(primitive(qrho), primitive(qp), gamma)

      write(unit, '(es24.16,7(",",es24.16))') &
        x(i), primitive(qrho), primitive(qu), primitive(qp), internal_energy, &
        conserved(iet, i), conserved(imx, i), sound_speed
    end do

    close(unit)
    ok = .true.
    message = ""
  end subroutine write_solution_csv

end module csv_io_mod
