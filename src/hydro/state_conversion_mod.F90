module state_conversion_mod
  use precision_mod, only: dp
  use constants_mod, only: density_floor, pressure_floor
  use state_indices_mod, only: &
    irho, imx, imy, imz, iet, ncons, qrho, qu, qv, qw, qp, nprim
  use eos_ideal_mod, only: ideal_gas_pressure, ideal_gas_total_energy_density
  implicit none
  private

  public :: primitive_to_conserved
  public :: conserved_to_primitive
  public :: state_is_physical

contains

  pure subroutine primitive_to_conserved(primitive, gamma, conserved, ok)
    real(dp), intent(in) :: primitive(nprim)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: conserved(ncons)
    logical, intent(out) :: ok
    real(dp) :: rho, total_energy_density

    conserved = 0.0_dp
    rho = primitive(qrho)

    if (rho <= density_floor .or. primitive(qp) <= pressure_floor) then
      ok = .false.
      return
    end if

    total_energy_density = ideal_gas_total_energy_density( &
      rho, primitive(qu), primitive(qv), primitive(qw), primitive(qp), gamma)
    if (total_energy_density <= 0.0_dp) then
      ok = .false.
      return
    end if

    conserved(irho) = rho
    conserved(imx) = rho * primitive(qu)
    conserved(imy) = rho * primitive(qv)
    conserved(imz) = rho * primitive(qw)
    conserved(iet) = total_energy_density
    ok = .true.
  end subroutine primitive_to_conserved

  pure subroutine conserved_to_primitive(conserved, gamma, primitive, ok)
    real(dp), intent(in) :: conserved(ncons)
    real(dp), intent(in) :: gamma
    real(dp), intent(out) :: primitive(nprim)
    logical, intent(out) :: ok
    real(dp) :: rho, pressure

    primitive = 0.0_dp
    rho = conserved(irho)
    if (rho <= density_floor) then
      ok = .false.
      return
    end if

    pressure = ideal_gas_pressure(conserved, gamma)
    if (pressure <= pressure_floor) then
      ok = .false.
      return
    end if

    primitive(qrho) = rho
    primitive(qu) = conserved(imx) / rho
    primitive(qv) = conserved(imy) / rho
    primitive(qw) = conserved(imz) / rho
    primitive(qp) = pressure
    ok = .true.
  end subroutine conserved_to_primitive

  pure logical function state_is_physical(conserved, gamma) result(is_physical)
    real(dp), intent(in) :: conserved(ncons)
    real(dp), intent(in) :: gamma
    real(dp) :: primitive(nprim)

    call conserved_to_primitive(conserved, gamma, primitive, is_physical)
  end function state_is_physical

end module state_conversion_mod
