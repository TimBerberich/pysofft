!> ------
!! @brief Defines pi
module math_constants
  use precision
  implicit none
  real(kind=dp), parameter :: pi = 4._dp*atan2(1._dp,1._dp)
  integer(kind=idp), parameter :: fubar = 123
end module math_constants
