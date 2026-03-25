!> ------
!! @brief Defines kind aliases for both fortran and c
module precision
  use, intrinsic :: iso_c_binding, only : c_int, c_int64_t, c_double, c_double_complex, c_bool,c_int32_t
  use, intrinsic :: iso_fortran_env, only : real64, int32, int64
  implicit none
  private

  ! Public kind parameters
  public :: dp,qp,isp,idp,c_dp,c_isp,c_idp, c_cdp,c_b

  integer, parameter :: dp = selected_real_kind(15,307)
  integer, parameter :: qp = selected_real_kind(33, 4931)
  integer, parameter :: isp = selected_int_kind(9)
  integer, parameter :: idp = selected_int_kind(18)
  
  ! Real kind for floating-point computations
  integer, parameter :: c_dp = c_double
  ! Integer kind for indices / counters
  integer, parameter :: c_idp = c_int64_t
  integer, parameter :: c_isp = c_int32_t
  ! Logical kind interoperable with C _Bool
  integer, parameter :: c_b = c_bool
  ! Complex kind interoperable with C double complex
  integer, parameter :: c_cdp = c_double_complex

end module precision
