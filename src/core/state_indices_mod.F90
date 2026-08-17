module state_indices_mod
  implicit none
  private

  ! Conserved variables used by the Phase-1 Euler solver.
  integer, parameter, public :: irho = 1
  integer, parameter, public :: imx  = 2
  integer, parameter, public :: imy  = 3
  integer, parameter, public :: imz  = 4
  integer, parameter, public :: iet  = 5
  integer, parameter, public :: ncons = 5

  ! Reserved PeleF base-state positions. Internal energy density and
  ! temperature are derived quantities in Phase 1 and are not advanced.
  integer, parameter, public :: iei  = 6
  integer, parameter, public :: item = 7
  integer, parameter, public :: nbase = 7

  ! Primitive variables.
  integer, parameter, public :: qrho = 1
  integer, parameter, public :: qu   = 2
  integer, parameter, public :: qv   = 3
  integer, parameter, public :: qw   = 4
  integer, parameter, public :: qp   = 5
  integer, parameter, public :: nprim = 5
end module state_indices_mod
