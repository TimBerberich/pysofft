program test
  use precision
  use soft
  integer(kind = dp),parameter :: bandw = 3
  real :: start, end
  complex(kind=dp) :: so3func(2*bandw,2*bandw,2*bandw)
  complex(kind = dp), allocatable :: coeff_out(:) 

  call init(bandw,.True.,.False.)
  empty_coeff = 1
  allocate(coeff_out(size(empty_coeff,1)))
  
  call cpu_time(start)
  do i=1,1000
     call inverse_wigner_trf_cmplx(empty_coeff,so3func)
     call forward_wigner_trf_cmplx(so3func,coeff_out)
     !call inverse_soft_precomputed_cmplx(empty_coeff,so3func)
  end do
  print *, sum(coeff_out)-size(coeff_out,1)
  call cpu_time(end)
  print * , end-start
end program test
