program test
  use precision
  use soft
  implicit none
  integer(kind = dp),parameter :: bandw = 32
  integer(kind=dp) :: i,j
  real :: start, end
  complex(kind=dp) :: so3func(2*bandw,2*bandw,2*bandw),so3func2(2*bandw,2*bandw,2*bandw),so3func3(2*bandw,2*bandw,2*bandw)
  complex(kind = dp), allocatable :: coeff_out(:) 

  call init(bandw,.True.,.False.)
  empty_coeff = 0
  call random_seed()
  call random_number(empty_coeff%re)
  call random_number(empty_coeff%im)
  
  allocate(coeff_out(size(empty_coeff,1)))
  coeff_out=0
  
  call cpu_time(start)
  do i=1,1
     !call inverse_wigner_trf_cmplx(empty_coeff,so3func)
     !call forward_wigner_trf_cmplx(so3func,coeff_out)
     call inverse_soft_cmplx(empty_coeff,so3func)
     call forward_soft_cmplx(so3func,coeff_out)

     !so3func=0
     !call random_number(so3func%re)
     !call random_number(so3func%im)
     !so3func3=0
     
     !call dfftw_execute_dft(plan_c2c_forward,so3func,so3func2)
     !call dfftw_execute_dft(plan_c2c_backward,so3func2,so3func3)
     !so3func3 = so3func3/((2*bw)**2)

  end do
  print *, sum(abs(empty_coeff-coeff_out))
  !write(*,'(F16.6)') so3func
  !print *
  !write(*,'(F16.6)') coeff_out

  !write(*,'(F16.6)') so3func3
  !print *, sum(coeff_out)-size(coeff_out,1)
  call cpu_time(end)
  print * , end-start
end program test
