module csv_io_3d_mod
  use precision_mod, only: dp
  use state_indices_mod, only: &
    imx, imy, imz, iet, ncons, nprim, qrho, qu, qv, qw, qp
  use state_conversion_mod, only: conserved_to_primitive
  implicit none
  private

  public :: write_solution_csv_3d

contains

  subroutine write_solution_csv_3d( &
      path, x, y, z, conserved, nx, ny, nz, gamma, ok, message)
    character(len=*), intent(in) :: path
    integer, intent(in) :: nx, ny, nz
    real(dp), intent(in) :: x(nx), y(ny), z(nz)
    real(dp), intent(in) :: conserved(ncons, nx, ny, nz)
    real(dp), intent(in) :: gamma
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message

    integer :: unit, io_status, i, j, k
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
      "x,y,z,rho,u,v,w,p,total_energy_density," // &
      "momentum_x_density,momentum_y_density,momentum_z_density"
    do k = 1, nz
      do j = 1, ny
        do i = 1, nx
          call conserved_to_primitive( &
            conserved(:, i, j, k), gamma, primitive, cell_ok)
          if (.not. cell_ok) then
            close(unit)
            ok = .false.
            write(message, '(a,i0,a,i0,a,i0,a)') &
              "Non-physical state while writing cell (", i, ",", j, &
              ",", k, ")"
            return
          end if
          write(unit, '(es24.16,11(",",es24.16))') &
            x(i), y(j), z(k), primitive(qrho), primitive(qu), &
            primitive(qv), primitive(qw), primitive(qp), &
            conserved(iet, i, j, k), conserved(imx, i, j, k), &
            conserved(imy, i, j, k), conserved(imz, i, j, k)
        end do
      end do
    end do

    close(unit)
    ok = .true.
    message = ""
  end subroutine write_solution_csv_3d

end module csv_io_3d_mod
