module constants_mod
  use precision_mod, only: dp
  implicit none
  private
  real(dp), parameter, public :: default_gamma = 1.4_dp
  real(dp), parameter, public :: density_floor = 1.0e-12_dp
  real(dp), parameter, public :: pressure_floor = 1.0e-12_dp
  real(dp), parameter, public :: tiny_speed = 1.0e-14_dp
  character(len=*), parameter, public :: pelef_version = "0.187.0"
end module constants_mod
