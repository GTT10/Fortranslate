module csv_io_2d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    imx, imy, iet, ncons, nprim, qrho, qu, qv, qp
  use state_conversion_mod, only: conserved_to_primitive
  implicit none
  private

  public :: write_solution_csv_2d

contains

  subroutine write_solution_csv_2d( &
      path, x, y, conserved, nx, ny, gamma, ok, message)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nx, ny
    real(dp), intent(in) :: x(nx), y(ny)
    real(dp), intent(in) :: conserved(ncons, nx, ny)
    real(dp), intent(in) :: gamma
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: unit, io_status, i, j
    real(dp) :: primitive(nprim)
    logical :: cell_ok

    open(newunit=unit, file=trim(path), status="replace", action="write", &
      iostat=io_status)
    if (io_status /= 0) then
      ok = .false.
      write(message, '(a,1x,a)') "Could not create output file:", trim(path)
      return
    end if

    write(unit, '(a)') &
      "x,y,rho,u,v,p,total_energy_density,momentum_x_density,momentum_y_density"

    do j = 1, ny
      do i = 1, nx
        call conserved_to_primitive( &
          conserved(:, i, j), gamma, primitive, cell_ok)
        if (.not. cell_ok) then
          close(unit)
          ok = .false.
          write(message, '(a,i0,a,i0)') &
            "Non-physical state while writing cell (", i, ",", j
          return
        end if

        write(unit, '(es24.16,8(",",es24.16))') &
          x(i), y(j), primitive(qrho), primitive(qu), primitive(qv), &
          primitive(qp), conserved(iet, i, j), conserved(imx, i, j), &
          conserved(imy, i, j)
      end do
    end do

    close(unit)
    ok = .true.
    message = ""
  end subroutine write_solution_csv_2d

end module csv_io_2d_mod
