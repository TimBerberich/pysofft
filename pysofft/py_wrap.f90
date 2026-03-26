!> --------
!! @brief Hack: Object oriented Fortran into f2py  
!!
!! f2py compatible wrappers for class `so3ft` in `softclass` module.
!! 
!! The problem is that custom types as well as newer object oriented features of Fortran are not supported by f2py.
!! We thus have to exclued the relevate code parts from beeing wraped by f2py using the `skip::` commandline option.
!!
!! Instead this module provides wrappers that expose a pointer to `so3ft` instances as integer to python.
!! Any call to a method first converts the provided integer into a `so3ft` pointer which allows access to the class methods from python since the wrapper function are f2py compatible.
module py
  !! Contains versions of the type bound procedures of so3ft that can be wrapped with f2py.
  use precision
  use softclass, only: so3ft,so3ft_ptr,kostelec_recurrence,risbo_recurrence
  use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  
contains
  function load_fftw_wisdom(file_path) result(error_code)
    ! error_code == 0 means something went wrong in fftw
    character(len=*),intent(in) :: file_path
    integer(kind = idp) :: i
    integer(C_INT) :: error_code
    character(kind=c_char), dimension(LEN(file_path)+1) :: c_path
    ! convert fortran to C string.
    do i = 1, LEN(file_path)
       c_path(i) = file_path(i:i)
    end do
    c_path(LEN(file_path)+1) = C_NULL_CHAR
    error_code = fftw_import_wisdom_from_filename(c_path)
  end function load_fftw_wisdom
  function save_fftw_wisdom(file_path) result(error_code)
    ! error_code == 0 means something went wrong in fftw
    character(len=*),intent(in) :: file_path
    integer(kind = idp) :: i
    integer(C_INT) :: error_code
    character(kind=c_char), dimension(LEN(file_path)+1) :: c_path
    ! convert fortran to C string.
    do i = 1, LEN(file_path)
       c_path(i) = file_path(i:i)
    end do
    c_path(LEN(file_path)+1) = C_NULL_CHAR
    error_code = fftw_export_wisdom_to_filename(c_path)
  end function save_fftw_wisdom

  function get_kostelec_recurrence_val(f2py_bug) result(id)
    logical, intent(in) :: f2py_bug
    logical :: f2py_bug_dummy
    !f2py logical :: f2py_bug = 0
    integer :: id
    f2py_bug_dummy = f2py_bug ! only there to suppress warning during f2py compilation.
    id = kostelec_recurrence
  end function get_kostelec_recurrence_val
  function get_risbo_recurrence_val(f2py_bug) result(id)
    logical, intent(in) :: f2py_bug
    logical :: f2py_bug_dummy
    !f2py logical :: f2py_bug = 0
    integer :: id
    f2py_bug_dummy = f2py_bug ! only there to suppress warning during f2py compilation.
    id = risbo_recurrence
  end function get_risbo_recurrence_val
  subroutine OMP_set_num_threads_(nthreads)
    integer(kind = isp), intent(in) :: nthreads
    call OMP_set_num_threads(nthreads)
  end subroutine OMP_set_num_threads_
  function OMP_get_max_threads_(f2py_bug) result(nthreads)
    !! This function does not need an input ...
    !! but f2py generates an error if it does not have one ...
    logical, intent(in) :: f2py_bug
    logical :: f2py_bug_dummy
    !f2py logical :: f2py_bug = 0
    integer(kind = isp) :: nthreads
    f2py_bug_dummy = f2py_bug ! only there to suppress warning during f2py compilation.
    nthreads = OMP_get_max_threads()
  end function OMP_get_max_threads_
  
  function py_init_soft(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags) result(self)
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp), intent(in) :: lmax
    integer(kind = idp), intent(in) :: recurrence_type
    logical, intent(in) :: init_ffts
    logical, intent(in) :: precompute_wigners
    integer(kind = isp), intent(in) :: fftw_flags
    type(so3ft_ptr) :: self_ptr
    integer(kind = idp) :: self
    ALLOCATE(self_ptr%p)
    self_ptr%p = so3ft(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
    self = TRANSFER(self_ptr,self)
  end function py_init_soft

  subroutine int_to_soft_pointer(self_int,self_ptr,self)
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr),intent(inout) :: self_ptr
    type(so3ft),intent(inout),pointer :: self
    self_ptr = TRANSFER(self_int,self_ptr)
    self => self_ptr%p
  end subroutine int_to_soft_pointer
  function py_get_bw(self_int) result(bw)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = idp) :: bw
    call int_to_soft_pointer(self_int,self_ptr,self)
    bw = self%bw
  end function py_get_bw
  function py_get_lmax(self_int) result(lmax)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = idp) :: lmax
    call int_to_soft_pointer(self_int,self_ptr,self)
    lmax = self%lmax
  end function py_get_lmax
  function py_wigners_are_precomputed(self_int) result(precomputed)
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    logical :: precomputed
    call int_to_soft_pointer(self_int,self_ptr,self)
    precomputed = (allocated(self%wigner_d) .AND. allocated(self%wigner_d_trsp))
  end function py_wigners_are_precomputed
  function py_get_recurrence_type(self_int) result(recurrence_type)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = idp) :: recurrence_type
    call int_to_soft_pointer(self_int,self_ptr,self)
    recurrence_type = self%recurrence_type
  end function py_get_recurrence_type
  subroutine py_set_lmax(self_int,lmax)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    integer(kind = idp),intent(in) :: lmax
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    self%lmax = lmax
  end subroutine py_set_lmax
  function py_get_fftw_flags(self_int) result(fftw_flags)
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = idp) :: fftw_flags
    call int_to_soft_pointer(self_int,self_ptr,self)
    fftw_flags = self%fftw_flags
  end function py_get_fftw_flags
  subroutine py_reset(self_int,bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp), intent(in) :: lmax
    !f2py integer :: lmax = bw - 1
    logical, intent(in) :: init_ffts
    !f2py logical :: init_ffts = 0
    logical, intent(in) :: precompute_wigners
    !f2py logical :: precompute_wigners = 0
    integer(kind = idp), intent(in) :: recurrence_type
    !f2py integer :: recurrence_type = 0
    ! kistelec_recurrence=0, risbo_recurrence=1
    integer(kind = isp), intent(in) :: fftw_flags
    !f2py integer :: fftw_flags = 0
    ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%init(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
  end subroutine py_reset
  subroutine py_destroy(self_int)
    integer(kind = idp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%destroy()

    DEALLOCATE(self)
  end subroutine py_destroy
  subroutine py_init_fft(self_int,use_real_fft)
    integer(kind = idp),intent(in) :: self_int
    logical, intent(in) :: use_real_fft
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%init_fft(use_real_fft)
  end subroutine py_init_fft

  subroutine py_inverse_wigner_trf_cmplx(self_int,coeff,so3func,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%inverse_wigner_trf_cmplx(coeff, so3func, use_mp)
  end subroutine py_inverse_wigner_trf_cmplx
  subroutine py_forward_wigner_trf_cmplx(self_int,so3func,coeff,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%forward_wigner_trf_cmplx(so3func,coeff, use_mp)
  end subroutine py_forward_wigner_trf_cmplx
  subroutine py_inverse_wigner_trf_real(self_int,coeff,so3func,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%inverse_wigner_trf_real(coeff, so3func, use_mp)
  end subroutine py_inverse_wigner_trf_real
  subroutine py_forward_wigner_trf_real(self_int,so3func,coeff,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%forward_wigner_trf_real(so3func,coeff, use_mp)
  end subroutine py_forward_wigner_trf_real


  subroutine py_isoft(self_int,coeff,so3func,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%isoft(coeff,so3func,use_mp)
  end subroutine py_isoft
  subroutine py_soft(self_int,so3func,coeff,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%soft(so3func,coeff,use_mp)
  end subroutine py_soft
  subroutine py_irsoft(self_int,coeff,so3func,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(in) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%irsoft(coeff,so3func,use_mp)
  end subroutine py_irsoft
  subroutine py_rsoft(self_int,so3func,coeff,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%rsoft(so3func,coeff,use_mp)
  end subroutine py_rsoft
  subroutine py_isoft_many(self_int,coeffs,so3funcs,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(in) :: coeffs(:,:)
    complex(kind = dp), intent(inout) :: so3funcs(:,:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%isoft_many(coeffs,so3funcs,use_mp)
  end subroutine py_isoft_many
  subroutine py_soft_many(self_int,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeffs(:,:)
    complex(kind = dp), intent(inout) :: so3funcs(:,:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%soft_many(so3funcs,coeffs,use_mp)
  end subroutine py_soft_many
  subroutine py_irsoft_many(self_int,coeffs,so3funcs,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(in) :: coeffs(:,:)
    real(kind = dp), intent(inout) :: so3funcs(:,:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%irsoft_many(coeffs,so3funcs,use_mp)
  end subroutine py_irsoft_many
  subroutine py_rsoft_many(self_int,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeffs(:,:)
    real(kind = dp), intent(inout) :: so3funcs(:,:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%rsoft_many(so3funcs,coeffs,use_mp)
  end subroutine py_rsoft_many

  function py_integrate_over_so3_cmplx(self_int,f) result(integral)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp),intent(in) :: f(:,:,:)
    complex(kind = dp) :: integral
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    integral =  self%integrate_over_so3_cmplx(f)
  end function py_integrate_over_so3_cmplx
  function py_integrate_over_so3_real(self_int,f) result(integral)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    real(kind = dp),intent(in) :: f(:,:,:)
    real(kind = dp) :: integral
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    integral =  self%integrate_over_so3_real(f)
  end function py_integrate_over_so3_real

  subroutine py_inverse_wigner_trf_corr_cmplx(self_int,f_lm,g_lm,fft_array,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%inverse_wigner_trf_corr_cmplx(f_lm,g_lm,fft_array,use_mp)
  end subroutine py_inverse_wigner_trf_corr_cmplx
  subroutine py_inverse_wigner_trf_corr_real(self_int,f_ml,g_ml,fft_array,use_mp)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%inverse_wigner_trf_corr_real(f_ml,g_ml,fft_array,use_mp)
  end subroutine py_inverse_wigner_trf_corr_real
  subroutine py_cross_correlation_ylm_cmplx(self_int,f_lm,g_lm,cc,use_mp)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: cc(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%cross_correlation_ylm_cmplx(f_lm,g_lm,cc,use_mp)
  end subroutine py_cross_correlation_ylm_cmplx
  subroutine py_cross_correlation_ylm_cmplx_3d(self_int,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp),target,intent(in) :: f_lms(:,:),g_lms(:,:)
    complex(kind = dp), intent(inout) :: cc(:,:,:)
    real(kind = dp), intent(in) :: radial_sampling_points(:)
    integer(kind = idp), intent(in) :: radial_limits(:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%cross_correlation_ylm_cmplx_3d(f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
  end subroutine py_cross_correlation_ylm_cmplx_3d
  subroutine py_cross_correlation_ylm_real(self_int,f_lm,g_lm,cc,use_mp)
    !f2py threadsafe
    integer(kind = idp),intent(in) :: self_int
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    real(kind = dp), intent(inout) :: cc(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%cross_correlation_ylm_real(f_lm,g_lm,cc,use_mp)
  end subroutine py_cross_correlation_ylm_real
  subroutine py_cross_correlation_ylm_real_3d(self_int,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
    !f2py threadsafe
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp),target,intent(in) :: f_lms(:,:),g_lms(:,:)
    real(kind = dp), intent(inout) :: cc(:,:,:)
    real(kind = dp), intent(in) :: radial_sampling_points(:)
    integer(kind = idp), intent(in) :: radial_limits(:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%cross_correlation_ylm_real_3d(f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
  end subroutine py_cross_correlation_ylm_real_3d

  subroutine py_fft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = idp), intent(in) :: self_int
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%fft(f1,f2)
  end subroutine py_fft
  subroutine py_ifft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = idp), intent(in) :: self_int
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%ifft(f1,f2)
  end subroutine py_ifft
  subroutine py_rfft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = idp), intent(in) :: self_int
    real(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%rfft(f1,f2)
  end subroutine py_rfft
  subroutine py_irfft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = idp), intent(in) :: self_int
    complex(kind = dp), intent(in) :: f1(:,:,:)
    real(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%irfft(f1,f2)
  end subroutine py_irfft

end module py
