program test
  use precision
  use soft
  implicit none
  integer(kind = dp),parameter :: bandw = 64
  integer(kind=dp) :: i,j
  real :: start, end
  complex(kind=dp) :: so3func(2*bandw,2*bandw,2*bandw),so3func2(2*bandw,2*bandw,2*bandw),so3func3(2*bandw,2*bandw,2*bandw)
  complex(kind = dp), allocatable :: coeff_out(:),coeff_in(:)

  !call init(bandw,.True.,.False.)
  call init(bandw,.False.,.True.,64)

  allocate(coeff_in(total_num_coeffs(bw)))
  coeff_in = 0
  call random_seed()
  call random_number(coeff_in%re)
  call random_number(coeff_in%im)

  call random_number(fft_r2c_out)
  !call random_number(so3func)
  
  allocate(coeff_out(total_num_coeffs(bw)))
  coeff_out=0

  fft_c2c_out = 0
  
  call cpu_time(start)
  do i=1,7*20
     call inverse_wigner_trf_cmplx(coeff_in,fft_c2c_out)
     call forward_wigner_trf_cmplx(fft_c2c_out,coeff_out)
     !call inverse_soft_cmplx(empty_coeff,so3func)
     !call forward_soft_cmplx(so3func,coeff_out)

     !so3func=0
     !call random_number(so3func%re)
     !call random_number(so3func%im)
     !so3func3=0
     
     !call dfftw_execute_dft(plan_c2c_forward,so3func,so3func2)
     !call dfftw_execute_dft(plan_c2c_backward,so3func2,so3func3)
     !so3func3 = so3func3/((2*bw)**2)
     !print*,shape(so3func),shape(fft_c2r_in)
     !call fft(so3func,so3func2)
     !print *, fft_c2r_in
     !call irfft(fft_c2r_in,so3func2)
     !print * ,'next'
     !print *, fft_c2r_in
  end do
  !print *, sum(abs(so3func-so3func2))
  print *, sum(abs(coeff_in-coeff_out))
  !write(*,'(F16.6)') so3func
  !print *
  !write(*,'(F16.6)') coeff_out

  !write(*,'(F16.6)') so3func3
  !print *, sum(coeff_out)-size(coeff_out,1)
  call cpu_time(end)
  print * , (end-start)/(7.0_dp*20.0_dp)
end program test
