!> --------
!! @brief Main Module containing the so3ft class which implements all transforms.
module softclass
  use precision
  use utils
  use make_wigner
  use omp_lib
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  
  ! f2py has issues if there is no variable declaration before the type declaration
  ! for unknown reasons so make sure there is at least one variable beeing declared
  integer :: risbo_recurrence = 1
  integer :: kostelec_recurrence = 0
  type :: so3ft
     integer(kind = idp) :: wigner_size
     integer(kind = idp) :: bw = 0
     integer(kind = idp) :: Lmax = 0
     integer(kind = idp) :: recurrence_type = 0
     real(kind = dp), allocatable :: wigner_d(:,:)
     real(kind = dp), allocatable :: wigner_d_trsp(:)
     real(kind = dp), allocatable :: wigner_dlml(:,:)
     real(kind = dp), allocatable :: trig_samples(:,:)
     real(kind = dp), allocatable :: trig_samples_risbo(:,:)
     real(kind = dp), allocatable :: sqrts_risbo(:)
     real(kind = dp), allocatable :: legendre_weights(:)
     real(kind = c_dp), pointer :: fft_r2c_out(:,:,:) => null()
     type(c_ptr) :: fft_r2c_out_p = c_null_ptr
     complex(kind = c_cdp), pointer :: fft_c2c_in(:,:,:) => null()
     type(c_ptr) :: fft_c2c_in_p = c_null_ptr
     complex(kind = c_cdp), pointer :: fft_c2r_in(:,:,:) => null()
     type(c_ptr) :: fft_c2r_in_p = c_null_ptr
     complex(kind = c_cdp), pointer :: fft_c2c_out(:,:,:) => null()
     type(c_ptr) :: fft_c2c_out_p = c_null_ptr
     integer(kind=dp) :: plan_c2c_forward,plan_c2c_backward,plan_r2c_forward,plan_c2r_backward
     type(c_ptr) :: plan_r2c_forward_p = c_null_ptr
     type(c_ptr) :: plan_c2r_backward_p = c_null_ptr
     type(c_ptr) :: plan_c2c_forward_p = c_null_ptr
     type(c_ptr) :: plan_c2c_backward_p = c_null_ptr
     integer(kind = isp) :: fftw_flags   ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
     logical :: plans_allocated_c  = .FALSE.
     logical :: plans_allocated_r = .FALSE.
   contains
     procedure :: destroy
     procedure :: destroy_fft
     procedure :: alloc_fft_arrays
     procedure :: dealloc_fft_arrays
     procedure :: init
     procedure :: init_kostelec_recurrence
     procedure :: init_risbo_recurrence
     procedure :: init_fft
     procedure :: init_wigners
     procedure :: check_recurrence_type

     ! transforms

     !! inverse wigner complex
     procedure :: inverse_wigner_loop_body_cmplx_kostelec_alloc
     procedure :: inverse_wigner_loop_body_cmplx_kostelec
     procedure :: inverse_wigner_trf_cmplx_kostelec
     procedure :: inverse_wigner_loop_body_cmplx_risbo
     procedure :: inverse_wigner_trf_cmplx_risbo
     procedure :: inverse_wigner_trf_cmplx
     !! forward wigner complex
     procedure :: forward_wigner_loop_body_cmplx_kostelec_alloc
     procedure :: forward_wigner_loop_body_cmplx_kostelec
     procedure :: forward_wigner_trf_cmplx_kostelec
     procedure :: forward_wigner_loop_body_cmplx_risbo
     procedure :: forward_wigner_trf_cmplx_risbo
     procedure :: forward_wigner_trf_cmplx
     !! inverse wigner real
     procedure :: inverse_wigner_loop_body_real_kostelec_alloc
     procedure :: inverse_wigner_loop_body_real_kostelec
     procedure :: inverse_wigner_trf_real_kostelec
     procedure :: inverse_wigner_loop_body_real_risbo
     procedure :: inverse_wigner_trf_real_risbo
     procedure :: inverse_wigner_trf_real
     !! forward wigner real
     procedure :: forward_wigner_loop_body_real_kostelec_alloc
     procedure :: forward_wigner_loop_body_real_kostelec
     procedure :: forward_wigner_trf_real_kostelec
     procedure :: forward_wigner_loop_body_real_risbo
     procedure :: forward_wigner_trf_real_risbo
     procedure :: forward_wigner_trf_real
     !! Full SO(3) transforms
     procedure :: isoft
     procedure :: isoft_
     procedure :: soft
     procedure :: soft_
     procedure :: rsoft
     procedure :: rsoft_
     procedure :: irsoft
     procedure :: irsoft_
     procedure :: isoft_many
     procedure :: soft_many
     procedure :: rsoft_many
     procedure :: irsoft_many

     ! derived functions
     procedure :: integrate_over_so3_cmplx
     procedure :: integrate_over_so3_real
     procedure :: inverse_wigner_loop_body_corr_cmplx_kostelec_alloc
     procedure :: inverse_wigner_loop_body_corr_cmplx_kostelec
     procedure :: inverse_wigner_trf_corr_cmplx_kostelec
     procedure :: inverse_wigner_loop_body_corr_cmplx_risbo
     procedure :: inverse_wigner_trf_corr_cmplx_risbo
     procedure :: inverse_wigner_trf_corr_cmplx
     procedure :: cross_correlation_ylm_cmplx
     procedure :: cross_correlation_ylm_cmplx_
     procedure :: cross_correlation_ylm_cmplx_3d
     procedure :: inverse_wigner_loop_body_corr_real_kostelec_alloc
     procedure :: inverse_wigner_loop_body_corr_real_kostelec
     procedure :: inverse_wigner_trf_corr_real_kostelec
     procedure :: inverse_wigner_loop_body_corr_real_risbo
     procedure :: inverse_wigner_trf_corr_real_risbo
     procedure :: inverse_wigner_trf_corr_real
     procedure :: cross_correlation_ylm_real_
     procedure :: cross_correlation_ylm_real
     procedure :: cross_correlation_ylm_real_3d

     ! handles for used fftw routines
     procedure :: fft
     procedure :: ifft
     procedure :: rfft
     procedure :: irfft
  end type so3ft
  interface so3ft
     module procedure :: init_soft
  end interface so3ft
  type :: so3ft_ptr
     type(so3ft),pointer :: p=>Null()
  end type so3ft_ptr
  abstract interface
     subroutine inverse_wigner_kostelec_interface(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
       import :: dp  ! Otherwise dp is undefined in _abstract interface blocks
       import :: so3ft
       class(so3ft),intent(in),target :: self
       complex(kind = dp), intent(in) :: coeff(:)
       complex(kind = dp), intent(inout) :: so3func(:,:,:)
       integer(kind=dp), intent(in) :: m1,m2
       real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
     end subroutine inverse_wigner_kostelec_interface
     subroutine forward_wigner_kostelec_interface(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
       import :: dp  ! Otherwise dp is undefined in _abstract interface blocks
       import :: so3ft
       class(so3ft),intent(in),target :: self
       complex(kind = dp), intent(inout) :: coeff(:)
       complex(kind = dp), intent(in) :: so3func(:,:,:)
       integer(kind=dp), intent(in) :: m1,m2
       real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
     end subroutine forward_wigner_kostelec_interface
     subroutine wigner_corr_interface(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
       import :: dp  ! Otherwise dp is undefined in _abstract interface blocks
       import :: so3ft
       class(so3ft),intent(in),target :: self
       complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
       complex(kind = dp), intent(inout) :: so3func(:,:,:)
       integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:),nm1_slice(:)
       real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
     end subroutine wigner_corr_interface
     subroutine wigner_corr_real_interface(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
       import :: so3ft
       import :: dp  ! Otherwise dp is undefined in _abstract interface blocks
       class(so3ft),intent(in),target :: self
       complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
       complex(kind = dp), intent(inout) :: so3func(:,:,:)
       integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:)
       real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
     end subroutine wigner_corr_real_interface
  end interface
  
contains
  ! Instance functions
  subroutine destroy_fft(self,apply_to_real_fft)
    class(so3ft),intent(inout) :: self
    logical, intent(in) :: apply_to_real_fft
    
    
    if (self%plans_allocated_c .AND. (.NOT. apply_to_real_fft)) then
       call self%dealloc_fft_arrays(.False.)
       call fftw_destroy_plan(self%plan_c2c_forward_p)
       call fftw_destroy_plan(self%plan_c2c_backward_p)
       self%plans_allocated_c = .FALSE.
    end if

    if (self%plans_allocated_r .AND. apply_to_real_fft) then
       call self%dealloc_fft_arrays(.True.)
       call fftw_destroy_plan(self%plan_r2c_forward_p)
       call fftw_destroy_plan(self%plan_c2r_backward_p)
       self%plans_allocated_r = .FALSE.
    end if
    
  end subroutine destroy_fft
  subroutine init_wigners(self)
    class(so3ft),intent(inout) :: self
    integer(kind = idp) :: m1, m2, d_shape(2),l_slice(2),trsp_slice(2),wig_len,lid,bw
    bw = self%bw
    d_shape = wigner_d_shape(bw)
    allocate(self%wigner_d(d_shape(1),d_shape(2)))
    if (self%recurrence_type == kostelec_recurrence) then
       call genwig_all_kostelec_preallocated(bw,self%wigner_d,.True.)
    else
       call genwig_all_risbo_preallocated(bw,self%wigner_d,.True.)
    end if
    allocate(self%wigner_d_trsp(d_shape(1)*d_shape(2)))
    
    do m1=0,bw-1
       do m2=m1, bw-1
          l_slice = mnl_to_flat_l_slice(m1,m2,bw)
          trsp_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
          wig_len = 2*bw*(bw-max(abs(m1),abs(m2)))          
          self%wigner_d_trsp(trsp_slice(1):trsp_slice(2)) = reshape(transpose(self%wigner_d(:,l_slice(1):l_slice(2))),[wig_len])
          ! multiply untransposed wigners with legendre weights
          do lid=l_slice(1), l_slice(2)
             self%wigner_d(:,lid) = self%wigner_d(:,lid)*self%legendre_weights
          end do
       end do
    end do
  end subroutine init_wigners
  subroutine init_fft(self,use_real_fft)
    logical, intent(in) :: use_real_fft
    class(so3ft),intent(inout) :: self
    integer(kind = c_isp) :: fft_rank,fft_n(2),fft_howmany,fft_inembed(2),fft_istride,fft_idist,fft_onembed(2),fft_ostride,fft_odist,bw2

    ! Destroy FFTs if they already exist.
    call self%destroy_fft(use_real_fft)
    
    bw2 = int(2_idp*self%bw,kind = isp)
    fft_rank = 2

    !wrong result but fast <=> np.fft2(,axes = (1,0))
    !fft_n = [bw2,bw2]
    !fft_howmany = bw2
    !fft_istride = 1
    !fft_ostride = 1
    !fft_inembed = fft_n
    !fft_onembed = fft_n
    !fft_idist = product(fft_n)
    !fft_odist = product(fft_n)
    
    !slower correct result !wrong result but fast <=> np.fft2(,axes = (2,1))
    fft_n = [bw2,bw2]
    fft_howmany = bw2
    fft_istride = bw2
    fft_ostride = bw2
    fft_inembed = fft_n
    fft_onembed = fft_n
    fft_idist = 1
    fft_odist = 1 ! product(fft_n)
    
    if (use_real_fft) then
       
       fft_onembed = [bw2/2+1,bw2]
       call self%alloc_fft_arrays(.TRUE.)
       self%plan_r2c_forward_p = fftw_plan_many_dft_r2c(&
            & fft_rank,fft_n,fft_howmany,&
            & self%fft_r2c_out,fft_inembed,fft_istride,fft_idist,&
            & self%fft_c2r_in,fft_onembed,fft_ostride,fft_odist,&
            & self%fftw_flags)       
       self%plan_c2r_backward_p = fftw_plan_many_dft_c2r(&
            & fft_rank,fft_n,fft_howmany,&
            & self%fft_c2r_in,fft_onembed,fft_ostride,fft_odist,&
            & self%fft_r2c_out,fft_inembed,fft_istride,fft_idist,&
            & self%fftw_flags)
       self%plans_allocated_r = .TRUE.
       ! fftw_plan_many_dft(rank,n,howmany,in,inembed,istride,idist
       ! rank (int) = dimension of fft
       ! n (array 1d) = shape of a single fft, array has to be of length rank
       ! howmany (int) = number of individual ffts to be computed
       ! in (array dim = rank+1) = input array
       ! inembed (array dim = rank) = n  (specifyies the subarray size of in to be used for each fft)
       ! istride (int) = distance between sucessive lements of the input data (here 1)
       ! idist (int) = distance between howmany different ffts (here prod(n))
       ! out (array dim = rank+1) = output array
       ! onembed (array dim = rank) = n (specifies the subarray size of the out to be used for each fft)
       ! ostride (int) = distance between sucessive elements of the output data set (here 1)
       ! odist (int) = distance betwene start of each of the howmany outputs (here prod(n))
       
    else
       fft_onembed = fft_n
       call self%alloc_fft_arrays(.FALSE.)
       self%plan_c2c_forward_p = fftw_plan_many_dft(&
            & fft_rank,fft_n,fft_howmany,&
            & self%fft_c2c_in,fft_inembed,fft_istride,fft_idist,&
            & self%fft_c2c_out,fft_onembed,fft_ostride,fft_odist,&
            & FFTW_FORWARD,self%fftw_flags)
       self%plan_c2c_backward_p = fftw_plan_many_dft(&
            & fft_rank,fft_n,fft_howmany,&
            & self%fft_c2c_out,fft_onembed,fft_ostride,fft_odist,&
            & self%fft_c2c_in,fft_inembed,fft_istride,fft_idist,&
            & FFTW_BACKWARD,self%fftw_flags)
       self%plans_allocated_c = .TRUE.
       ! fftw_plan_many_dft(rank,n,howmany,in,inembed,istride,idist
       ! rank (int) = dimension of fft
       ! n (array 1d) = shape of a single fft, array has to be of length rank
       ! howmany (int) = number of individual ffts to be computed
       ! in (array dim = rank+1) = input array
       ! inembed (array dim = rank) = n  (specifyies the subarray size of in to be used for each fft)
       ! istride (int) = distance between sucessive lements of the input data (here 1)
       ! idist (int) = distance between howmany different ffts (here prod(n))
       ! out (array dim = rank+1) = output array
       ! onembed (array dim = rank) = n (specifies the subarray size of the out to be used for each fft)
       ! ostride (int) = distance between sucessive elements of the output data set (here 1)
       ! odist (int) = distance betwene start of each of the howmany outputs (here prod(n))
    end if
  end subroutine init_fft
  subroutine alloc_fft_arrays(self,use_real_fft)
    logical, intent(in) :: use_real_fft
    class(so3ft), intent(inout) :: self
    integer(kind = idp) :: bw2
    bw2=2_idp*self%bw

    if (use_real_fft) then
       self%fft_r2c_out_p = fftw_alloc_real(int(bw2 * bw2 * bw2, C_SIZE_T))
       call c_f_pointer(self%fft_r2c_out_p, self%fft_r2c_out, [bw2,bw2,bw2])
       self%fft_c2r_in_p = fftw_alloc_complex(int(bw2 * (bw2/2+1) * bw2, C_SIZE_T))
       call c_f_pointer(self%fft_c2r_in_p, self%fft_c2r_in, [bw2,bw2/2+1,bw2])
    else
       self%fft_c2c_in_p = fftw_alloc_complex(int(bw2 * bw2 * bw2, C_SIZE_T))
       call c_f_pointer(self%fft_c2c_in_p, self%fft_c2c_in, [bw2,bw2,bw2])
       self%fft_c2c_out_p = fftw_alloc_complex(int(bw2 * bw2 * bw2, C_SIZE_T))
       call c_f_pointer(self%fft_c2c_out_p, self%fft_c2c_out, [bw2,bw2,bw2])
    end if
  end subroutine alloc_fft_arrays
  subroutine dealloc_fft_arrays(self,use_real_fft)
    logical, intent(in) :: use_real_fft
    class(so3ft), intent(inout) :: self
    
    if (use_real_fft) then
       if (c_associated(self%fft_c2r_in_p))then
          call fftw_free(self%fft_c2r_in_p)
          self%fft_c2r_in_p = c_null_ptr
       end if
       if (c_associated(self%fft_r2c_out_p))then
          call fftw_free(self%fft_r2c_out_p)
          self%fft_r2c_out_p = c_null_ptr
       end if
    else
       if (c_associated(self%fft_c2c_in_p))then
          call fftw_free(self%fft_c2c_in_p)
          self%fft_c2c_in_p = c_null_ptr
       end if
       if (c_associated(self%fft_c2c_out_p))then
          call fftw_free(self%fft_c2c_out_p)
          self%fft_c2c_out_p = c_null_ptr
       end if
    end if
  end subroutine dealloc_fft_arrays
  subroutine check_recurrence_type(self)
    class(so3ft),intent(inout) :: self
    if ((self%recurrence_type /= kostelec_recurrence) .AND. (self%recurrence_type /= risbo_recurrence)) then
        print *, "Error: unknown recurrence type", self%recurrence_type
        stop
     end if
  end subroutine check_recurrence_type
  subroutine init(self,bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp), intent(in) :: lmax
    integer(kind = idp), intent(in) :: recurrence_type
    logical, intent(in) :: init_ffts,precompute_wigners
    integer(kind = isp), intent(in) :: fftw_flags
    class(so3ft),intent(inout) :: self

    call self%destroy()

    self%bw = bw
    self%lmax = lmax
    self%fftw_flags = fftw_flags
    self%recurrence_type = recurrence_type
    call self%check_recurrence_type()
    allocate(self%legendre_weights(2*bw))
    self%legendre_weights = legendre_quadrature_weights(bw)
    self%wigner_size = size_wigner_d(bw)
    if (precompute_wigners) then
       call self%init_wigners()
    else
       if (self%recurrence_type == kostelec_recurrence) then
          call self%init_kostelec_recurrence()
       else
          call self%init_risbo_recurrence()
       end if
    end if
    if (init_ffts) then
       call self%init_fft(.FALSE.)
       call self%init_fft(.TRUE.)
    end if
  end subroutine init
  subroutine init_kostelec_recurrence(self)
    class(so3ft),intent(inout) :: self
    integer(kind=dp) :: bw
    bw = self%bw
    allocate(self%trig_samples(2*bw,3))
    self%trig_samples = create_trig_samples(bw)
    allocate(self%wigner_dlml(2_idp*bw,triangular_size(bw)))
    self%wigner_dlml = compute_all_dlml_l_contiguous(bw, self%trig_samples(:,2), self%trig_samples(:,3),.TRUE.)
  end subroutine init_kostelec_recurrence
  subroutine init_risbo_recurrence(self)
    class(so3ft),intent(inout) :: self
    integer(kind=dp) :: bw
    bw = self%bw
    allocate(self%trig_samples_risbo(2*bw,2))
    self%trig_samples_risbo = create_trig_samples_risbo(bw)
    allocate(self%sqrts_risbo(2*bw-1))
    self%sqrts_risbo = create_sqrts_risbo(bw)
  end subroutine init_risbo_recurrence
  
  function init_soft(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags) result(self)
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp), intent(in) :: lmax
    integer(kind = idp), intent(in) :: recurrence_type
    logical, intent(in) :: init_ffts
    logical, intent(in) :: precompute_wigners
    integer(kind = isp), intent(in) :: fftw_flags
    type(so3ft) :: self
    call self%init(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
  end function init_soft
  subroutine destroy(self)
    class(so3ft),intent(inout) :: self
    ! Reset soft module variables
    if (allocated(self%wigner_d)) then
       deallocate(self%wigner_d)
    end if
    if (allocated(self%wigner_d_trsp)) then
       deallocate(self%wigner_d_trsp)
    end if
    if (allocated(self%wigner_dlml)) then
       deallocate(self%wigner_dlml)
    end if
    if (allocated(self%legendre_weights)) then
       deallocate(self%legendre_weights)
    end if
    if (allocated(self%trig_samples)) then
       deallocate(self%trig_samples)
    end if
    if (allocated(self%sqrts_risbo)) then
       deallocate(self%sqrts_risbo)
    end if
    if (allocated(self%trig_samples_risbo)) then
       deallocate(self%trig_samples_risbo)
    end if
    call self%destroy_fft(.False.)
    call self%destroy_fft(.True.)
    self%bw=0
    self%lmax=0
  end subroutine destroy
  
  ! inverse  wigner complex
  subroutine inverse_wigner_loop_body_cmplx_kostelec_alloc(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    complex(kind = dp) :: coeff_part(self%bw-m2)
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_slice(2),bw,bw2,t_slice(2)
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    
    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    
    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))

    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    c_slice = coeff_slice_mnl(m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))
    so3func(:,m1+1,m2+1) = matmul(coeff_part,wig_mat)
    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    c_slice = coeff_slice_mnl(-m2,-m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))
    so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
    
    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       c_slice = coeff_slice_mnl(m2,m1,bw)
       coeff_part = coeff(c_slice(1):c_slice(2))*(sym_const_m1*sym_const_m2)
       so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice_mnl(-m1,-m2,bw)
       coeff_part = coeff(c_slice(1):c_slice(2))*(sym_const_m1*sym_const_m2)
       so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    c_slice = coeff_slice_mnl(m1,-m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice_mnl(-m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice_mnl(m2,-m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)


    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice_mnl(-m2,m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
  end subroutine inverse_wigner_loop_body_cmplx_kostelec_alloc
  subroutine inverse_wigner_loop_body_cmplx_kostelec(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp) :: wig_mat(self%bw-m2,2*self%bw)
    complex(kind = dp) :: coeff_part(self%bw-m2)
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_slice(2),bw,bw2,dlml_id

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genwig_l2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id),.True.)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    c_slice = coeff_slice_mnl(m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))
    so3func(:,m1+1,m2+1) = matmul(coeff_part,wig_mat)
    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    c_slice = coeff_slice_mnl(-m2,-m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))
    so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
    
    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       c_slice = coeff_slice_mnl(m2,m1,bw)
       coeff_part = coeff(c_slice(1):c_slice(2))*(sym_const_m1*sym_const_m2)
       so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice_mnl(-m1,-m2,bw)
       coeff_part = coeff(c_slice(1):c_slice(2))*(sym_const_m1*sym_const_m2)
       so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    c_slice = coeff_slice_mnl(m1,-m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice_mnl(-m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice_mnl(m2,-m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)


    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice_mnl(-m2,m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
  end subroutine inverse_wigner_loop_body_cmplx_kostelec
  subroutine inverse_wigner_trf_cmplx_kostelec(self,coeff,so3func,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    integer(kind = idp) :: i,m1,m2,L
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(inverse_wigner_kostelec_interface),pointer :: loop_body => Null()

    ! Select loop_body
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_cmplx_kostelec_alloc
    else
       loop_body => inverse_wigner_loop_body_cmplx_kostelec
    end if
       
    ! initiallizing some constants
    L = self%lmax
    do i=0,L
       sym_array(i+1) = (-1.0)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(i,m1,m2,sym_const_m1) SHARED(so3func,coeff,sym_array)
       !$OMP DO
       do i=1,((L+1)*(L+2))/2
          call flat_to_triangular_index(m1,m2,i,L)
          !print *,i,m1,m2
          sym_const_m1 = (-1.0)**m1
          call loop_body(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
         end do
       !$OMP END DO
       !$OMP END PARALLEL 
    else
       do m1=0, L
          sym_const_m1 = (-1.0)**m1
          do m2=m1, L
             call loop_body(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
          end do
        end do
    end if
  end subroutine inverse_wigner_trf_cmplx_kostelec
  subroutine inverse_wigner_loop_body_cmplx_risbo(self,coeff,so3func,dlmn,l,m1,m2,sym_const_l,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: l,m1,m2
    real(kind=dp),intent(in) :: sym_const_l,sym_const_m1
    real(kind=dp),intent(in) :: dlmn(:)
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_loc,c_slice(2),bw,bw2,o,i

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    c_slice = coeff_slice_lmn(l,m1,m2)
    c_loc = c_slice(1)
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !$OMP SIMD 
    do i=1,bw2
       so3func(i,m1+1,m2+1) = so3func(i,m1+1,m2+1) + coeff(c_loc)*dlmn(i)
    end do
    !$OMP END SIMD

    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    s_ids = order_to_ids(-m2,-m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + coeff(c_loc+1)*dlmn(i)
    end do
    !$OMP END SIMD
    
    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + (sym_const_m1*sym_const_m2*coeff(c_loc+2_idp))*dlmn(i)
       end do
       !$OMP END SIMD

       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + (sym_const_m1*sym_const_m2*coeff(c_loc+3_idp))*dlmn(i)
       end do
       !$OMP END SIMD
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    o = MERGE(2_idp,0_idp,(m1/=m2))
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+o+2_idp))*dlmn(i)
    end do
    !$OMP END SIMD

    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m2*coeff(c_loc+o+3_idp))*dlmn(i)
    end do
    !$OMP END SIMD
        
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i,s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i,s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+6_idp))*dlmn(i)
    end do
    !$OMP END SIMD

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i,s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i,s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m2*coeff(c_loc + 7_idp))*dlmn(i)
    end do
    !$OMP END SIMD
  end subroutine inverse_wigner_loop_body_cmplx_risbo
  subroutine inverse_wigner_trf_cmplx_risbo(self,coeff,so3func,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    integer(kind = idp) :: m1,m2,l,mnid
    real(kind=dp) :: sym_const_m1,sym_const_l
    real(kind=dp), allocatable :: dl(:,:)
    ! zeroing so3func needed since it will be populated by +=
    so3func = 0
    
    allocate(dl(2*self%bw,(self%bw*(self%bw+1))/2))
    if (use_mp) then       
       ! non-fft part of the SO(3) fourier transform
       !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,dl,sym_const_l,l) SHARED(so3func,coeff)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_idp,((l+1_idp)*(l+2_idp))/2_idp
             call flat_to_triangular_index(m1,m2,mnid,l)
             sym_const_m1 = (-1.0)**m1 
             call inverse_wigner_loop_body_cmplx_risbo(self,coeff,so3func,dl(:,mnid),l,m1,m2,sym_const_l,sym_const_m1)
          end do
          !$OMP END DO
       end do
       !$OMP END PARALLEL       
    else
       ! non-fft part of the SO(3) fourier transform
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_idp)
                call inverse_wigner_loop_body_cmplx_risbo(self,coeff,so3func,dl(:,mnid),l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
    deallocate(dl)
  end subroutine inverse_wigner_trf_cmplx_risbo
  subroutine inverse_wigner_trf_cmplx(self,coeff,so3func,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (self%recurrence_type == kostelec_recurrence) then
       call self%inverse_wigner_trf_cmplx_kostelec(coeff, so3func, use_mp)
    else
       if (allocated(self%wigner_d_trsp)) then
          ! genwig_all_risbo_preallocated has been used compute wigners and store them in mnl order
          ! compatible with faster kostelec transfom
          call self%inverse_wigner_trf_cmplx_kostelec(coeff, so3func, use_mp)
       else
          call self%inverse_wigner_trf_cmplx_risbo(coeff,so3func,use_mp)
       end if
    end if
  end subroutine inverse_wigner_trf_cmplx
  ! forward  wigner complex
  subroutine forward_wigner_loop_body_cmplx_kostelec_alloc(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2     
    integer(kind=dp) :: s_ids(2),c_slice(2),bw,bw2,l_slice(2)
    complex(kind=dp) :: so3func_part(2*self%bw)

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    l_slice = mnl_to_flat_l_slice(m1, m2, bw)
    wig_mat(1:bw2,1:bw-m2) => self%wigner_d(:,l_slice(1):l_slice(2))
    !call get_wigner_matrix(m1,m2,wig_mat,wig_mat_arr)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    c_slice = coeff_slice_mnl(m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(so3func(:,m1+1,m2+1),wig_mat)
    !return
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!

    s_ids = order_to_ids(-m2,-m1,bw)
    !print * , -m2,-m1
    !write(*,'(F16.5, F16.5)', advance='yes') so3func(:,s_ids(1),s_ids(2))
    c_slice = coeff_slice_mnl(-m2,-m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(so3func(:,s_ids(1),s_ids(2)),wig_mat)
    !return    
    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!       
       s_ids = order_to_ids(m2,m1,bw)
       c_slice = coeff_slice_mnl(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(so3func(:,s_ids(1),s_ids(2)),wig_mat)

       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice_mnl(-m1,-m2,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(so3func(:,s_ids(1),s_ids(2)),wig_mat)
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    c_slice = coeff_slice_mnl(m1,-m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice_mnl(-m1,m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m2
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice_mnl(m2,-m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice_mnl(-m2,m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m2
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
  end subroutine forward_wigner_loop_body_cmplx_kostelec_alloc
  subroutine forward_wigner_loop_body_cmplx_kostelec(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp) :: sym_const_m2     
    real(kind = dp) :: wig_mat(2*self%bw,self%bw-m2) ! possible source of optimization can make it allocatable
    integer(kind=dp) :: s_ids(2),c_slice(2),bw,bw2,dlml_id

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genwig_l2(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id),.True.)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    c_slice = coeff_slice_mnl(m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(:,m1+1,m2+1),wig_mat)
    !return
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!

    s_ids = order_to_ids(-m2,-m1,bw)
    !print * , -m2,-m1
    !write(*,'(F16.5, F16.5)', advance='yes') so3func(:,s_ids(1),s_ids(2))
    c_slice = coeff_slice_mnl(-m2,-m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)
    !return    
    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!       
       s_ids = order_to_ids(m2,m1,bw)
       c_slice = coeff_slice_mnl(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(self%legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)

       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice_mnl(-m1,-m2,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(self%legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    c_slice = coeff_slice_mnl(m1,-m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice_mnl(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m2*sym_array(m2+1:))
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice_mnl(m2,-m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice_mnl(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m2*sym_array(m2+1:))
  end subroutine forward_wigner_loop_body_cmplx_kostelec
  subroutine forward_wigner_trf_cmplx_kostelec(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: i,m1,m2,L
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(forward_wigner_kostelec_interface),pointer :: loop_body => Null()

    ! Select loop_body
    if (allocated(self%wigner_d)) then
       loop_body => forward_wigner_loop_body_cmplx_kostelec_alloc
    else
       loop_body => forward_wigner_loop_body_cmplx_kostelec
    end if
       
    ! initiallizing some constants
    L = self%lmax
    do i=0,L
       sym_array(i+1) = (-1)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(i,m1,m2,sym_const_m1) SHARED(so3func,coeff,sym_array)
       !$OMP DO
       do i=1,(L+1)*(L+2)/2
          call flat_to_triangular_index(m1,m2,i,L)
          sym_const_m1 = (-1.0)**m1
          call loop_body(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
       end do
       !$OMP END DO
       !$OMP END PARALLEL 
    else
       do m1=0, L
          sym_const_m1 = (-1.0)**m1
          do m2=m1, L
             call loop_body(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
          end do
       end do
    end if
  end subroutine forward_wigner_trf_cmplx_kostelec
  subroutine forward_wigner_loop_body_cmplx_risbo(self,so3func,coeff,dlmn,l,m1,m2,sym_const_l,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    real(kind = dp), intent(in) :: dlmn(:)
    integer(kind=dp), intent(in) ::l,m1,m2
    real(kind=dp),intent(in) :: sym_const_l,sym_const_m1
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_slice(2),c_loc,bw,bw2,o

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    c_slice = coeff_slice_lmn(l,m1,m2)
    c_loc = c_slice(1)

    ! This method assiumes 0<=m1<=m2<=l
    ! which also means m = max(abs(m1),abs(m2)) = m2

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!        
    coeff(c_loc) = DOT_PRODUCT(dlmn,so3func(:,m1+1,m2+1))
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    s_ids = order_to_ids(-m2,-m1,bw)
    coeff(c_loc+1_idp) = DOT_PRODUCT(dlmn,so3func(:,s_ids(1),s_ids(2)))
   

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids = order_to_ids(m2,m1,bw)
       coeff(c_loc + 2_idp) = DOT_PRODUCT(dlmn,so3func(:,s_ids(1),s_ids(2)))*sym_const_m1*sym_const_m2
       
       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)       
       coeff(c_loc + 3_idp) = DOT_PRODUCT(dlmn,so3func(:,s_ids(1),s_ids(2)))*sym_const_m1*sym_const_m2
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    o=MERGE(2_idp,0_idp,(m1/=m2))
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    coeff(c_loc + o + 2_idp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_m1*sym_const_l
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    coeff(c_loc + o + 3_idp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_m2*sym_const_l

    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    coeff(c_loc + 6_idp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2))) *sym_const_m1*sym_const_l
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    coeff(c_loc + 7_idp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_m2*sym_const_l

  end subroutine forward_wigner_loop_body_cmplx_risbo
  subroutine forward_wigner_trf_cmplx_risbo(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: m1,m2,l,mnid
    real(kind=dp) :: sym_const_m1,sym_const_l,dl_tmp(2_idp*self%bw)
    real(kind=dp), allocatable :: dl(:,:)

    allocate(dl(2*self%bw,(self%bw*(self%bw+1))/2))
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,dl_tmp,dl,sym_const_l,l) SHARED(so3func,coeff)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_idp,((l+1_idp)*(l+2_idp))/2_idp
             call flat_to_triangular_index(m1,m2,mnid,l)
             sym_const_m1 = (-1.0)**m1
             dl_tmp = dl(:,mnid)*self%legendre_weights
             call forward_wigner_loop_body_cmplx_risbo(self,so3func,coeff,dl_tmp,l,m1,m2,sym_const_l,sym_const_m1)
          end do
          !$OMP END DO
       end do
       !$OMP END PARALLEL
    else
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_idp)
                dl_tmp = dl(:,mnid)*self%legendre_weights
                call forward_wigner_loop_body_cmplx_risbo(self,so3func,coeff,dl_tmp,l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
    deallocate(dl)
  end subroutine forward_wigner_trf_cmplx_risbo
  subroutine forward_wigner_trf_cmplx(self,so3func,coeff,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (self%recurrence_type == kostelec_recurrence) then
       call self%forward_wigner_trf_cmplx_kostelec(so3func,coeff, use_mp)
    else
       if (allocated(self%wigner_d)) then
          ! genwig_all_risbo_preallocated has been used compute wigners and store them in mnl order
          ! compatible with faster kostelec transfom
          call self%forward_wigner_trf_cmplx_kostelec(so3func,coeff,use_mp)
       else
          call self%forward_wigner_trf_cmplx_risbo(so3func,coeff,use_mp)
       end if
    end if
  end subroutine forward_wigner_trf_cmplx
  ! inverse  wigner real
  subroutine inverse_wigner_loop_body_real_kostelec_alloc(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw,bw2,t_slice(2)

    sym_const_m2 = (-1._dp)**m2
    bw = self%bw
    bw2 = 2_idp*bw

    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice_mnl(m1,m2,bw)
    so3func(:,m1+1,m2+1) = matmul(coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)),wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    c_nm2nm1_slice = coeff_slice_mnl(-m2,-m1,bw)
    so3func(:,s_ids(1),s_ids(2))=CONJG(matmul(coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)),wig_mat))

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    c_pm1nm2_slice = coeff_slice_mnl(m1,-m2,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)),wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice_mnl(m2,-m1,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)),wig_mat)    
  end subroutine inverse_wigner_loop_body_real_kostelec_alloc
  subroutine inverse_wigner_loop_body_real_kostelec(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    !real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2
    real(kind = dp) :: wig_mat(self%bw-m2,2*self%bw)
    complex(kind = dp):: coeff_part(self%bw-m2)
    integer(kind=dp) :: s_ids(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw,bw2,dlml_id

    sym_const_m2 = (-1._dp)**m2
    bw = self%bw
    bw2 = 2_idp*bw

    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genwig_l2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id),.True.)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice_mnl(m1,m2,bw)
    so3func(:,m1+1,m2+1) = matmul(coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)),wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    if (m1/=m2) then
       !! m2,m1 !!
       s_ids = order_to_ids(m2,m1,bw)
       c_nm2nm1_slice = coeff_slice_mnl(-m2,-m1,bw)
       so3func(:,s_ids(1),s_ids(2))=CONJG(matmul(coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)),wig_mat))
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! pre@vents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    c_pm1nm2_slice = coeff_slice_mnl(m1,-m2,bw)
    coeff_part = coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice_mnl(m2,-m1,bw)
    coeff_part = coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)    
  end subroutine inverse_wigner_loop_body_real_kostelec
  subroutine inverse_wigner_trf_real_kostelec(self,coeff,so3func,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(in) :: coeff(:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: i,m1,m2,L,s_ids(2),s_ids_sym(2),bw
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(inverse_wigner_kostelec_interface),pointer :: loop_body

    bw  = self%bw    
    ! Select loop_body
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_real_kostelec_alloc
    else
       loop_body => inverse_wigner_loop_body_real_kostelec
    end if
    
    ! initiallizing some constants
    L = self%lmax
    do i=0,L
       sym_array(i+1) = (-1)**i 
    end do

    if (use_mp) then
       ! non-fft part of the SO(3) fourier transform

       !$OMP PARALLEL PRIVATE(i,m1,m2,sym_const_m1,s_ids,s_ids_sym) SHARED(so3func,coeff,sym_array)
       !$OMP DO
       do i=1,(L+1)*(L+2)/2
          call flat_to_triangular_index(m1,m2,i,L)
          sym_const_m1 = (-1.0)**m1
          call loop_body(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
       end do
       !$OMP END DO

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*       
       !$OMP DO
       do m2=1,L
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          so3func(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(:,s_ids(1),s_ids(2)))
       end do
       !$OMP END DO
       !$OMP END PARALLEL 
     else
        ! non-fft part of the SO(3) fourier transform        
       do m1=0, L
          sym_const_m1 = (-1.0)**m1
          do m2=m1, L
             call loop_body(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,L
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          so3func(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(:,s_ids(1),s_ids(2)))
       end do
    end if
  end subroutine inverse_wigner_trf_real_kostelec
  subroutine inverse_wigner_loop_body_real_risbo(self,coeff,so3func,dlmn,l,m1,m2,sym_const_l,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: l,m1,m2
    real(kind=dp),intent(in) :: sym_const_l,sym_const_m1
    real(kind=dp),intent(in) :: dlmn(:)
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_loc,c_slice(2),bw,bw2,o,i

    sym_const_m2 = (-1._dp)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    c_slice = coeff_slice_lmn(l,m1,m2)
    c_loc = c_slice(1)    
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    !$OMP SIMD 
    do i=1,bw2
       so3func(i,m1+1,m2+1) = so3func(i,m1+1,m2+1) + coeff(c_loc)*dlmn(i)
    end do
    !$OMP END SIMD
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    if (m1/=m2) then
       !! m2,m1 !!
       s_ids = order_to_ids(m2,m1,bw)
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + CONJG(coeff(c_loc+1)*dlmn(i))
       end do
       !$OMP END SIMD
    end if


    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! pre@vents sign swaps on 0 ids which are already covered
    o = MERGE(2_idp,0_idp,(m1/=m2))
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+o+2_idp))*dlmn(i)
    end do
    !$OMP END SIMD
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i,s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i,s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+6_idp))*dlmn(i)
    end do
    !$OMP END SIMD       
  end subroutine inverse_wigner_loop_body_real_risbo
  subroutine inverse_wigner_trf_real_risbo(self,coeff,so3func,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(in) :: coeff(:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: j,i,m1,m2,l,mnid,s_ids(2),s_ids_sym(2),bw,bw2
    real(kind=dp) :: sym_const_m1,sym_const_l
    real(kind = dp), allocatable :: dl(:,:) 

    bw  = self%bw
    bw2  = 2_idp*bw
    ! I had issues with segfaults due to stack overflow in the mp routines
    ! so lets make dl allocatable
    allocate(dl(2*self%bw,(self%bw*(self%bw+1))/2))
    
    if (use_mp) then
       ! non-fft part of the SO(3) fourier transform
       !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,s_ids,s_ids_sym,i,j,dl,sym_const_l,l) SHARED(so3func,coeff,bw,bw2)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_idp,((l+1_idp)*(l+2_idp))/2_idp
             call flat_to_triangular_index(m1,m2,mnid,l)
             sym_const_m1 = (-1.0)**m1 
             call inverse_wigner_loop_body_real_risbo(self,coeff,so3func,dl(:,mnid),l,m1,m2,sym_const_l,sym_const_m1)
          end do
          !$OMP END DO
       end do
       !$OMP END PARALLEL
       
       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       !$OMP PARALLEL PRIVATE(m2,s_ids,s_ids_sym,i,j) SHARED(so3func,bw,bw2)
       !$OMP DO
       do j=1,(self%lmax+1_idp)*bw2
          m2 = ((j-1_idp)/bw2)
          i = MOD(j,bw2)+1_idp
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          so3func(i,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(i,s_ids(1),s_ids(2)))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       ! non-fft part of the SO(3) fourier transform
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_idp)
                call inverse_wigner_loop_body_real_risbo(self,coeff,so3func,dl(:,mnid),l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,self%lmax
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          do i=1,bw2
             so3func(i,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(i,s_ids(1),s_ids(2)))
          end do
       end do
    end if
    deallocate(dl)
  end subroutine inverse_wigner_trf_real_risbo
  subroutine inverse_wigner_trf_real(self,coeff,so3func,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (self%recurrence_type == kostelec_recurrence) then
       call self%inverse_wigner_trf_real_kostelec(coeff, so3func, use_mp)
    else
       if (allocated(self%wigner_d_trsp)) then
          ! genwig_all_risbo_preallocated has been used compute wigners and store them in mnl order
          ! compatible with faster kostelec transfom
          call self%inverse_wigner_trf_real_kostelec(coeff, so3func, use_mp)
       else
          call self%inverse_wigner_trf_real_risbo(coeff,so3func,use_mp)
       end if
    end if
  end subroutine inverse_wigner_trf_real
  ! forward  wigner real
  subroutine forward_wigner_loop_body_real_kostelec_alloc(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    complex(kind = dp), intent(inout) :: coeff(:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2
    complex(kind = dp) :: so3func_part(2*self%bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw,bw2,l_slice(2)

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw

    ! load wigner d matrices
    l_slice = mnl_to_flat_l_slice(m1, m2, bw)
    wig_mat(1:bw2,1:bw-m2) => self%wigner_d(:,l_slice(1):l_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice_mnl(m1,m2,bw)
    coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)) = matmul(so3func(:,m1+1,m2+1),wig_mat)

    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    so3func_part = get_so3func_part_halfcomplex(bw,-m2,-m1,so3func)
    c_nm2nm1_slice = coeff_slice_mnl(-m2,-m1,bw)
    coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)) = matmul(so3func_part,wig_mat)

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    !! Due to real input both transforms have already been computed  !!
    !! Use symmetry $$D^l_{m1,m2}(\alpha,\beta,\gamma) =(-1)^{m2-m1} (D^l_{-m1,-m2}(\alpha,\beta,\gamma))^*$$
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       c_slice = coeff_slice_mnl(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)))
       
       !! -m1,-m2 !!
       c_slice = coeff_slice_mnl(-m1,-m2,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)))
    end if
    

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    c_pm1nm2_slice = coeff_slice_mnl(m1,-m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m1,m2 !!
    !! Use real symmetry
    c_slice = coeff_slice_mnl(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice_mnl(m2,-m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m2,m1 !!
    c_slice = coeff_slice_mnl(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)))
  end subroutine forward_wigner_loop_body_real_kostelec_alloc
  subroutine forward_wigner_loop_body_real_kostelec(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    complex(kind = dp), intent(inout) :: coeff(:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp) :: sym_const_m2
    real(kind = dp) :: wig_mat(2*self%bw,self%bw-m2)
    complex(kind = dp) :: so3func_part(2*self%bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw,bw2,dlml_id

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    
    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genwig_l2(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id),.True.)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice_mnl(m1,m2,bw)
    coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)) = matmul(self%legendre_weights*so3func(:,m1+1,m2+1),wig_mat)
    
    
    if (m2==0 .AND. m2==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    so3func_part = get_so3func_part_halfcomplex(bw,-m2,-m1,so3func)
    c_nm2nm1_slice = coeff_slice_mnl(-m2,-m1,bw)
    coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)) = matmul(self%legendre_weights*so3func_part,wig_mat)


    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    !! Due to real input both transforms have already been computed  !!
    !! Use symmetry $$D^l_{m1,m2}(\alpha,\beta,\gamma) =(-1)^{m2-m1} (D^l_{-m1,-m2}(\alpha,\beta,\gamma))^*$$
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       c_slice = coeff_slice_mnl(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)))
       
       !! -m1,-m2 !!
       c_slice = coeff_slice_mnl(-m1,-m2,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)))
    end if
    

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    c_pm1nm2_slice = coeff_slice_mnl(m1,-m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)) = matmul(self%legendre_weights*so3func_part,wig_mat)*sym_array(m2+1:)
   
    !! -m1,m2 !!
    !! Use real symmetry
    c_slice = coeff_slice_mnl(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice_mnl(m2,-m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)) = matmul(self%legendre_weights*so3func_part,wig_mat)*sym_array(m2+1:)
    
    
    !! -m2,m1 !!
    c_slice = coeff_slice_mnl(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)))
  end subroutine forward_wigner_loop_body_real_kostelec
  subroutine forward_wigner_trf_real_kostelec(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: i,m1,m2,L
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(forward_wigner_kostelec_interface),pointer :: loop_body

    ! Select loop_body
    if (allocated(self%wigner_d)) then
       loop_body => forward_wigner_loop_body_real_kostelec_alloc
    else
       loop_body => forward_wigner_loop_body_real_kostelec
    end if
    
    ! initiallizing some constants
    L = self%lmax
    do i=0,L
       sym_array(i+1) = (-1)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(i,m1,m2,sym_const_m1) SHARED(so3func,coeff,sym_array)
       !$OMP DO
       do i=1,(L+1)*(L+2)/2
          call flat_to_triangular_index(m1,m2,i,L)
          sym_const_m1 = (-1.0)**m1             
          call loop_body(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
       end do
       !$OMP END DO
       !$OMP END PARALLEL 
    else
       do m1=0, L
          sym_const_m1 = (-1.0)**m1
          do m2=m1, L
             call loop_body(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
          end do
       end do
    end if
  end subroutine forward_wigner_trf_real_kostelec
  subroutine forward_wigner_loop_body_real_risbo(self,so3func,coeff,dlmn,l,m1,m2,sym_const_l,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    real(kind = dp), intent(in) :: dlmn(:)
    integer(kind=dp), intent(in) ::l,m1,m2
    real(kind=dp),intent(in) :: sym_const_l,sym_const_m1
    real(kind = dp) :: sym_const_m2
    integer(kind=dp) :: s_ids(2),c_slice(2),c_loc,bw,bw2,o

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_idp*bw
    c_slice = coeff_slice_lmn(l,m1,m2)
    c_loc = c_slice(1)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    coeff(c_loc) = DOT_PRODUCT(dlmn,so3func(:,m1+1,m2+1))    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    coeff(c_loc+1_idp) = DOT_PRODUCT(dlmn,CONJG(so3func(:,s_ids(1),s_ids(2))))



    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    !! Due to real input both transforms have already been computed  !!
    !! Use symmetry $$D^l_{m1,m2}(\alpha,\beta,\gamma) =(-1)^{m2-m1} (D^l_{-m1,-m2}(\alpha,\beta,\gamma))^*$$
    if (m1/=m2) then
       !!  m2,m1  !!
       coeff(c_loc + 2_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 1_idp))
       
       !! -m1,-m2 !!
       coeff(c_loc + 3_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc))
    end if
    
    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    o=MERGE(2_idp,0_idp,(m1/=m2))
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    coeff(c_loc + o + 2_idp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_l*sym_const_m1
   
    !! -m1,m2 !!
    !! Use real symmetry
    coeff(c_loc + o + 3_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + o + 2_idp)) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    coeff(c_loc + 6_idp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_l*sym_const_m1    
    
    !! -m2,m1 !!
    coeff(c_loc + 7_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 6_idp))
  end subroutine forward_wigner_loop_body_real_risbo
  subroutine forward_wigner_trf_real_risbo(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: m1,m2,l,mnid
    real(kind=dp) :: sym_const_m1,sym_const_l,dl_tmp(2_idp*self%bw)
    real(kind=dp), allocatable :: dl(:,:)

    ! I had issues with segfaults due to stack overflow in the mp routines
    ! so lets make dl allocatable
    allocate(dl(2*self%bw,(self%bw*(self%bw+1))/2))
    
    if (use_mp) then
      !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,dl_tmp,dl,sym_const_l,l) SHARED(so3func,coeff)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_idp,((l+1_idp)*(l+2_idp))/2_idp
             call flat_to_triangular_index(m1,m2,mnid,l)
             sym_const_m1 = (-1.0)**m1
             dl_tmp = dl(:,mnid)*self%legendre_weights
             call forward_wigner_loop_body_real_risbo(self,so3func,coeff,dl_tmp,l,m1,m2,sym_const_l,sym_const_m1)
          end do
          !$OMP END DO
       end do
       !$OMP END PARALLEL
    else       
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_idp)
                dl_tmp = dl(:,mnid)*self%legendre_weights
                call forward_wigner_loop_body_real_risbo(self,so3func,coeff,dl_tmp,l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
    deallocate(dl)
  end subroutine forward_wigner_trf_real_risbo
  subroutine forward_wigner_trf_real(self,so3func,coeff,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (self%recurrence_type == kostelec_recurrence) then
       call self%forward_wigner_trf_real_kostelec(so3func,coeff,use_mp)
    else
       if (allocated(self%wigner_d)) then
          ! genwig_all_risbo_preallocated has been used compute wigners and store them in mnl order
          ! compatible with faster kostelec transfom
          call self%forward_wigner_trf_real_kostelec(so3func,coeff,use_mp)
       else
          call self%forward_wigner_trf_real_risbo(so3func,coeff,use_mp)
       end if
    end if
  end subroutine forward_wigner_trf_real

  ! SO(3) FFTs
  subroutine isoft_(self,coeff,so3func,fft_array,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical,intent(in) :: use_mp

    fft_array=0._dp
    call self%inverse_wigner_trf_cmplx(coeff,fft_array,use_mp)
    call fftw_execute_dft(self%plan_c2c_forward_p,fft_array,so3func)
    so3func = so3func * (1._dp/(2._dp*pi))
  end subroutine isoft_
  subroutine isoft(self,coeff,so3func,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    call self%isoft_(coeff,so3func,self%fft_c2c_in,use_mp)
  end subroutine isoft
  subroutine soft_(self,so3func,coeff,fft_array,use_mp)    
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    fft_array = 0._dp
    call fftw_execute_dft(self%plan_c2c_backward_p,so3func,fft_array)
    fft_array = fft_array * (2._dp*pi/real(2_idp*self%bw,kind=dp)**2) ! * [(2*pi/(2*bw) * 1/sqrt(2*pi))] * [2*pi/(2*bw) * 1/sqrt(2*pi)]
    ! 2*pi/(2*bw) corrects the integration range of a single fourier transform while 1/(sqrt(2*pi)) is due to the used normalization
    ! of the Wigner d matrices, that is  (D^l_m,n)^{ortho} = \sqrt{(2l+1)/8\pi^2}, where \sqrt{(2l+1)/2} is already contained in
    ! the small wigner matrices d^l_mn
    call self%forward_wigner_trf_cmplx(fft_array,coeff,use_mp)
  end subroutine soft_
  subroutine soft(self,so3func,coeff,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    
    call self%soft_(so3func,coeff,self%fft_c2c_out,use_mp)
  end subroutine soft
  subroutine irsoft_(self,coeff,so3func,fft_array,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(in) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical,intent(in) :: use_mp

    fft_array=0._dp
    call self%inverse_wigner_trf_real(coeff,fft_array,use_mp)
    fft_array = CONJG(fft_array) ! to correct for the fact that we have to compute the forward not the backward fft.
    call fftw_execute_dft_c2r(self%plan_c2r_backward_p,fft_array,so3func)
    so3func = so3func * (1._dp/(2._dp*pi))
  end subroutine irsoft_
  subroutine irsoft(self,coeff,so3func,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(in) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    call self%irsoft_(coeff,so3func,self%fft_c2r_in,use_mp)
  end subroutine irsoft
  subroutine rsoft_(self,so3func,coeff,fft_array,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    fft_array=0._dp
    call fftw_execute_dft_r2c(self%plan_r2c_forward_p,so3func,fft_array)
    fft_array = fft_array * (2._dp*pi/real(2_idp*self%bw,kind=dp)**2) ! * 1/(2*bw) * 1/(2*bw)
    fft_array = CONJG(fft_array) ! to correct for the fact that we have to compute the backward not the forward fft.
    call self%forward_wigner_trf_real(fft_array,coeff,use_mp)    
  end subroutine rsoft_
  subroutine rsoft(self,so3func,coeff,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    call self%rsoft_(so3func,coeff,self%fft_c2r_in,use_mp)
  end subroutine rsoft
  subroutine soft_many(self,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    complex(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(inout) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    complex(kind=dp) :: fft_c2c_in(self%bw*2,self%bw*2,self%bw*2)
    complex(kind=dp), allocatable :: fft_c2c_in_2(:,:,:)
    integer(kind=dp) :: n,i
    n = size(so3funcs,4)

    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    
    if (use_mp) then
       ! We can not use self%fft_c2c_in in openMP since it is not thread private and openMP does not support 
       ! declaration of derived type arguments as private, so we create a local copy, fft_c2c_in.
       
       !$OMP PARALLEL PRIVATE(i,fft_c2c_in_2) SHARED(so3funcs,coeffs,n)
       !! allocatable array is precaution to not cause stack overflows
       !! for large self%bw
       allocate(fft_c2c_in_2(self%bw*2,self%bw*2,self%bw*2))
       !$OMP DO
       do i=1,n
          call self%soft_(so3funcs(:,:,:,i),coeffs(:,i),fft_c2c_in_2,.False.)
       end do
       !$OMP END DO
       deallocate(fft_c2c_in_2)
       !$OMP END PARALLEL
    else
       do i=1,n
          call self%soft_(so3funcs(:,:,:,i),coeffs(:,i),fft_c2c_in,.FALSE.)
       end do
    end if
       
  end subroutine soft_many
  subroutine isoft_many(self,coeffs,so3funcs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    complex(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(in) :: coeffs(:,:)
    complex(kind=dp) :: fft_c2c_in(self%bw*2,self%bw*2,self%bw*2)
    complex(kind=dp), allocatable :: fft_c2c_in_2(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind=dp) :: n,i
    n = size(coeffs,2)

    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    
    if (use_mp) then
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,fft_c2c_in_2)
       !! allocatable array is precaution to not cause stack overflows
       !! for large self%bw
       allocate(fft_c2c_in_2(self%bw*2,self%bw*2,self%bw*2))
       !$OMP DO
       do i=1,n
          call self%isoft_(coeffs(:,i),so3funcs(:,:,:,i),fft_c2c_in_2,.FALSE.)       
       end do
       !$OMP END DO
       deallocate(fft_c2c_in_2)
       !$OMP END PARALLEL
    else
       do i=1,n
          call self%isoft_(coeffs(:,i),so3funcs(:,:,:,i),fft_c2c_in,.FALSE.)       
       end do
    end if
  end subroutine isoft_many
  subroutine rsoft_many(self,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    real(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(inout) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    complex(kind=dp) :: fft_c2r_in(self%bw*2,self%bw+1,self%bw*2)
    complex(kind=dp), allocatable :: fft_c2r_in_2(:,:,:)
    integer(kind=dp) :: n,i
    n = size(so3funcs,4)

    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(i,fft_c2r_in_2) SHARED(so3funcs,coeffs,n)
       !! allocatable array is precaution to not cause stack overflows
       !! for large self%bw
       allocate(fft_c2r_in_2(self%bw*2,self%bw+1,self%bw*2))
       !$OMP DO
       do i=1,n
          call self%rsoft_(so3funcs(:,:,:,i),coeffs(:,i),fft_c2r_in_2,.FALSE.)
       end do
       !$OMP END DO
       deallocate(fft_c2r_in_2)
       !$OMP END PARALLEL
    else
       do i=1,n
          call self%rsoft_(so3funcs(:,:,:,i),coeffs(:,i),fft_c2r_in,.FALSE.)
          !call self%rsoft(so3funcs(:,:,:,i),coeffs(:,i),.FALSE.)       
       end do
    end if
  end subroutine rsoft_many
  subroutine irsoft_many(self,coeffs,so3funcs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    real(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(in) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    complex(kind=dp) :: fft_c2r_in(self%bw*2,self%bw+1,self%bw*2)
    complex(kind=dp),allocatable :: fft_c2r_in_2(:,:,:)
    integer(kind=dp) :: n,i
    n = size(coeffs,2)

    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(i,fft_c2r_in_2) SHARED(so3funcs,coeffs,n)
       allocate(fft_c2r_in_2(self%bw*2,self%bw+1,self%bw*2))
       !$OMP DO
       do i=1,n
          call self%irsoft_(coeffs(:,i),so3funcs(:,:,:,i),fft_c2r_in_2,.FALSE.)       
       end do
       !$OMP END DO
       deallocate(fft_c2r_in_2)
       !$OMP END PARALLEL
    else
       do i=1,n
          call self%irsoft_(coeffs(:,i),so3funcs(:,:,:,i),fft_c2r_in,.FALSE.)       
       end do
    end if
  end subroutine irsoft_many

  ! Numerical integration over SO(3)
  function integrate_over_so3_cmplx(self,f) result(integral)
    ! Integrates f(\alpha,\gamma,\beta) over SO(3)
    ! If f = 1 then the result sould be 8*pi**2
    ! In this case SUM(SUM(f,1),1) = (2*bw)**2 and SUM(legendre_weights)=2
    ! => The normalization constant has to be (pi/bw)**2
    ! Maybe the legendre weights should be divided by 2, it is strange for their sum to be 2
    class(so3ft),intent(in) :: self
    complex(kind = dp),intent(in) :: f(:,:,:)
    complex(kind = dp) :: integral
    integral = SUM(SUM(SUM(f,1),1)*self%legendre_weights)
    integral = integral*(pi/real(self%bw,kind=dp))**2
  end function integrate_over_so3_cmplx
  function integrate_over_so3_real(self,f) result(integral)
    ! Integrates f(\alpha,\gamma,\beta) over SO(3)
    ! If f = 1 then the result sould be 8*pi**2
    ! In this case SUM(SUM(f,1),1) = (2*bw)**2 and SUM(legendre_weights)=2
    ! => The normalization constant has to be (pi/bw)**2
    ! Maybe the legendre weights should be divided by 2, it is strange for their sum to be 2
    class(so3ft),intent(in) :: self
    real(kind = dp),intent(in) :: f(:,:,:)
    real(kind = dp) :: integral
    integral = SUM(SUM(SUM(f,1),1)*self%legendre_weights)
    integral = integral*(pi/real(self%bw,kind=dp))**2
  end function integrate_over_so3_real

  ! correlation wigner inverse complex
  subroutine inverse_wigner_loop_body_corr_cmplx_kostelec_alloc(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:),nm1_slice(:)
    real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn(self%bw-m2)
    integer(kind=dp) :: s_ids(2),bw,bw2,nm2_slice(2),pm2_slice(2),l_start,t_slice(2)

    sym_const_m2 = (-1._dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_idp
    bw = self%bw
    bw2 = 2_idp*bw
    pm2_slice = MLc_slice(m2,bw)
    nm2_slice = MLc_slice(-m2,bw) 
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1m2
    so3func(:,m1+1,m2+1) = matmul(cc_lmn,wig_mat)
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(nm2_slice(1):nm2_slice(2)) * g_ml(nm1_slice(1):nm1_slice(2)) * sym_const_m1m2
    so3func(:,s_ids(1),s_ids(2)) = matmul(cc_lmn,wig_mat)

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)) * sym_const_m1m2
       so3func(:,s_ids(1),s_ids(2)) = matmul(cc_lmn,wig_mat)*sym_const_m1m2       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       cc_lmn = wig_norm(l_start:) * f_ml(nm1_slice(1):nm1_slice(2)) * g_ml(nm2_slice(1):nm2_slice(2)) * sym_const_m1m2
       so3func(:,s_ids(1),s_ids(2)) = matmul(cc_lmn,wig_mat) * sym_const_m1m2
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(nm2_slice(1):nm2_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)
    
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(nm1_slice(1):nm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(cc_lmn,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(nm1_slice(1):nm1_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(nm2_slice(1):nm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(cc_lmn,wig_mat)
  end subroutine inverse_wigner_loop_body_corr_cmplx_kostelec_alloc
  subroutine inverse_wigner_loop_body_corr_cmplx_kostelec(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:),nm1_slice(:)
    real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
    real(kind = dp) :: wig_mat(self%bw-m2,2*self%bw)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn(self%bw-m2)
    integer(kind=dp) :: s_ids(2),bw,bw2,nm2_slice(2),pm2_slice(2),l_start,dlml_id

    sym_const_m2 = (-1._dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_idp
    bw = self%bw
    bw2 = 2_idp*bw
    pm2_slice = MLc_slice(m2,bw)
    nm2_slice = MLc_slice(-m2,bw) 
    
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genwig_l2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id),.True.)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    !print *, 'wig_norm',wig_norm(l_start:)
    !print *, 'sym_const',sym_const_m1m2
    !print *, 'flm',f_ml(pm1_slice(1):pm1_slice(2))
    !print *, 'glm', g_ml(pm2_slice(1):pm2_slice(2))
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1m2
    !print *, 'cc_lmn',cc_lmn
    so3func(:,m1+1,m2+1) = matmul(cc_lmn,wig_mat)
    
    if (m1==0 .AND. m2==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(nm2_slice(1):nm2_slice(2)) * g_ml(nm1_slice(1):nm1_slice(2)) * sym_const_m1m2
    so3func(:,s_ids(1),s_ids(2)) = matmul(cc_lmn,wig_mat)

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)) * sym_const_m1m2
       so3func(:,s_ids(1),s_ids(2)) = matmul(cc_lmn,wig_mat)*sym_const_m1m2       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       cc_lmn = wig_norm(l_start:) * f_ml(nm1_slice(1):nm1_slice(2)) * g_ml(nm2_slice(1):nm2_slice(2)) * sym_const_m1m2
       so3func(:,s_ids(1),s_ids(2)) = matmul(cc_lmn,wig_mat) * sym_const_m1m2
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(nm2_slice(1):nm2_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)
    
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(nm1_slice(1):nm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(cc_lmn,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(nm1_slice(1):nm1_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(nm2_slice(1):nm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)) * sym_const_m1m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(cc_lmn,wig_mat)
  end subroutine inverse_wigner_loop_body_corr_cmplx_kostelec
  subroutine inverse_wigner_trf_corr_cmplx_kostelec(self,f_lm,g_lm,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    complex(kind = dp) ::  f_ml(size(f_lm)),g_ml(size(f_lm))
    integer(kind = idp) :: i,n,m1,m2,m,l,lmax,bw,pm1_slice(2),nm1_slice(2),pm1_tmp(2),nm1_tmp(2)   
    real(kind=dp) :: sym_const_m1,sym_array(self%bw),wig_norm(self%bw)
    procedure(wigner_corr_interface),pointer :: loop_body
    
    bw = self%bw
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_corr_cmplx_kostelec_alloc
    else
       loop_body => inverse_wigner_loop_body_corr_cmplx_kostelec
    end if
    

    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    fft_array = 0._dp
    ! initiallizing some constants
    lmax = self%lmax
    do i=0,lmax
       sym_array(i+1) = (-1.0)**i
       wig_norm(i+1_idp) = 2._dp*pi*SQRT(2._dp/real(2_idp*i+1_idp,kind = dp))
    end do
    if (use_mp) then
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,l,m,m1,m2,sym_const_m1,pm1_slice,nm1_slice)

       ! "Transpose" input coefficient layout to be adapted to the wigner memory layout (l contiguous)   
       !$OMP DO 
       do i=1,bw**2
          call flat_to_pyramid_index(l,m,i)
          f_ml(MLc(m,l,bw)) = f_lm(LMc(l,m))
          g_ml(MLc(m,l,bw)) = CONJG(g_lm(LMc(l,m)))
       end do
       !$OMP END DO
       
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2
       !$OMP DO
       do i=1,(lmax+1)*(lmax+2)/2
          call flat_to_triangular_index(m1,m2,i,lmax)
          sym_const_m1 = (-1.0)**m1
          pm1_slice = MLc_slice(m1,bw)
          nm1_slice = MLc_slice(-m1,bw)
          pm1_slice(1) = pm1_slice(1) + m2-m1
          nm1_slice(1) = nm1_slice(1) + m2-m1
          
          call loop_body(self,f_ml,g_ml,fft_array,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       ! "Transpose" input coefficient layout to be adapted to the wigner memory layout (l contiguous)
       do l=0,bw-1
          do m=-l,l
             f_ml(MLc(m,l,bw)) = f_lm(LMc(l,m))
             g_ml(MLc(m,l,bw)) = CONJG(g_lm(LMc(l,m)))
          end do
       end do
       
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2       
       do m1=0, lmax
          sym_const_m1 = (-1.0)**m1
          pm1_slice = MLc_slice(m1,bw)
          pm1_tmp = pm1_slice
          nm1_slice = MLc_slice(-m1,bw)
          nm1_tmp = nm1_slice
          do m2=m1, lmax
             n = m2 - m1
             pm1_tmp(1) = pm1_slice(1) + n
             nm1_tmp(1) = nm1_slice(1) + n 
             call loop_body(self,f_ml,g_ml,fft_array,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_tmp,nm1_tmp)
          end do
       end do
    end if
  end subroutine inverse_wigner_trf_corr_cmplx_kostelec
  subroutine inverse_wigner_loop_body_corr_cmplx_risbo(self,f_lm,g_lm,so3func,dlmn,l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: l,m1,m2
    real(kind=dp),intent(in) :: wig_norm,sym_const_l,sym_const_m1
    real(kind=dp),intent(in) :: dlmn(:)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn
    integer(kind=dp) :: s_ids(2),bw,bw2,l_start,i,nm1,nm2,pm1,pm2

    sym_const_m2 = (-1._dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_idp
    bw = self%bw
    bw2 = 2_idp*bw
    pm1 = LMc(l,m1)
    nm1 = LMc(l,-m1)
    pm2 = LMc(l,m2)
    nm2 = LMc(l,-m2)
    
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !print *, 'wig_norm',wig_norm
    !print *, 'sym_const',sym_const_m1m2
    !print *, 'flm',f_lm(pm1)
    !print *, 'glm', CONJG(g_lm(pm2))
    cc_lmn = wig_norm * f_lm(pm1) * CONJG(g_lm(pm2)) * sym_const_m1m2
    !print *, 'cc_lmn',cc_lmn
    !$OMP SIMD 
    do i=1,bw2
       so3func(i,m1+1,m2+1) = so3func(i,m1+1,m2+1) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
    
    if (m1==0 .AND. m2==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    s_ids=order_to_ids(-m2,-m1,bw)
    cc_lmn = wig_norm * f_lm(nm2) * CONJG(g_lm(nm1)) * sym_const_m1m2
    !$OMP SIMD
    do i=1,bw2
       so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (m1/=m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       cc_lmn = wig_norm * f_lm(pm2) * CONJG(g_lm(pm1))
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
       end do
       !$OMP END SIMD

       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       cc_lmn = wig_norm * f_lm(nm1) * CONJG(g_lm(nm2))
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
       end do
       !$OMP END SIMD
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    cc_lmn = wig_norm * f_lm(pm1) * CONJG(g_lm(nm2)) * sym_const_m2 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
    
    
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    cc_lmn = wig_norm * f_lm(nm1) * CONJG(g_lm(pm2)) * sym_const_m1 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD

    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm * f_lm(pm2) * CONJG(g_lm(nm1)) * sym_const_m2 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    cc_lmn = wig_norm * f_lm(nm2) * CONJG(g_lm(pm1)) * sym_const_m1 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
  end subroutine inverse_wigner_loop_body_corr_cmplx_risbo
  subroutine inverse_wigner_trf_corr_cmplx_risbo(self,f_lm,g_lm,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: mnid,m1,m2,l,bw
    real(kind=dp) :: sym_const_m1,sym_const_l,wig_norm,dl(2*self%bw,(self%bw*(self%bw+1))/2)
        
    bw = self%bw
 
    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    fft_array = 0._dp

    if (use_mp) then
       ! non-fft part of the SO(3) fourier transform
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_idp*l+1_idp,kind = dp))
          !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1) SHARED(fft_array,f_lm,g_lm,dl,sym_const_l,l)
          !$OMP DO
          do mnid=1_idp,((l+1_idp)*(l+2_idp))/2_idp
             call flat_to_triangular_index(m1,m2,mnid,l)
             sym_const_m1 = (-1.0)**m1 
             call inverse_wigner_loop_body_corr_cmplx_risbo(self,f_lm,g_lm,fft_array,dl(:,mnid),l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
          end do
          !$OMP END DO
          !$OMP END PARALLEL
       end do
    else
       ! non-fft part of the SO(3) fourier transform
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_idp*l+1_idp,kind = dp))
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_idp)
                call inverse_wigner_loop_body_corr_cmplx_risbo(self,f_lm,g_lm,fft_array,dl(:,mnid),l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
  end subroutine inverse_wigner_trf_corr_cmplx_risbo
  subroutine inverse_wigner_trf_corr_cmplx(self,f_lm,g_lm,fft_array,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    if (self%recurrence_type == kostelec_recurrence) then
       call self%inverse_wigner_trf_corr_cmplx_kostelec(f_lm,g_lm, fft_array, use_mp)
    else
       if (allocated(self%wigner_d_trsp)) then
          ! genwig_all_risbo_preallocated has been used compute wigners and store them in mnl order
          ! compatible with faster kostelec transfom
          call self%inverse_wigner_trf_corr_cmplx_kostelec(f_lm,g_lm, fft_array, use_mp)
       else
          call self%inverse_wigner_trf_corr_cmplx_risbo(f_lm,g_lm, fft_array, use_mp)
       end if
    end if
  end subroutine inverse_wigner_trf_corr_cmplx
  ! cross correlations ylm complex
  subroutine cross_correlation_ylm_cmplx_(self,f_lm,g_lm,cc,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: cc(:,:,:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp

    ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2
    call self%inverse_wigner_trf_corr_cmplx(f_lm, g_lm,fft_array,use_mp)
    ! Compute inverse fft
    call fftw_execute_dft(self%plan_c2c_forward_p,fft_array,cc)
    cc = cc * (1/(2._dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)   
  end subroutine cross_correlation_ylm_cmplx_
  subroutine cross_correlation_ylm_cmplx(self,f_lm,g_lm,cc,use_mp)
    ! let f,g be two square integrable functions on the 2 sphere 
    ! Define CC: SO(3) ---> \mathbb{R},   R |---> <f,g \circ R> = \int_{S^2} dx f(x)*\overline{g(Rx)}
    ! This function calculates CC(R) for all R defined in make_SO3_grid
    ! Spherical harmonic coefficient arrays , g_lm,f_lm, are assumed to follow the indexing conventions of the SHTNS package,
    ! i.e. those given by the pure functions LMc LMr from above.
    !
    ! arguments :
    !   f_lm: f_{l,m} spherical harmonic coefficients of f 
    !   g_lm: g_{l,m} spherical harmonic coefficients of g 
    !            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
    !
    !
    ! output :
    !   C_values: array of shape (2*bw,2*bw,2*bw)
    !   The maximum in C corresponds to the Rotation that best maps f to g
    !
    !   The Idis of the maximum in C i_max,j_max,k_max correspond to the euler anglees beta,alpha,gamma in this order.
    !   Somehow there is still a bug. The resulting euler angles need to be modified as follows to yield correct results
    !   alpha,beta,gamma -----> 2*pi - alpha, beta , 2*pi-gamma
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: cc(:,:,:)
    logical, intent(in) :: use_mp

    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    call self%cross_correlation_ylm_cmplx_(f_lm,g_lm,cc,self%fft_c2c_in,use_mp)
  end subroutine cross_correlation_ylm_cmplx
  subroutine cross_correlation_ylm_cmplx_3d(self,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
    ! let f,g be two square integrable functions on the $\mathbb{R}^3$ in spherical coordinates 
    ! Define CC: SO(3) ---> \mathbb{R},   R |---> <f,g \circ R> = \int_R dr r^2 \int_{S^2} dw f(r,w)*\overline{g(r,Rw)}
    ! This function calculates CC(R) for all R defined in make_SO3_grid
    ! Spherical harmonic coefficient arrays , g_lm,f_lm, are assumed to follow the indexing conventions of the SHTNS package,
    ! i.e. those given by the pure functions LMc LMr from above.
    !
    ! arguments :
    !   f_lm: f_{l,m} spherical harmonic coefficients of f 
    !   g_lm: g_{l,m} spherical harmonic coefficients of g 
    !            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
    !
    !
    ! output :
    !   C_values: array of shape (2*bw,2*bw,2*bw)
    !   The maximum in C corresponds to the Rotation that best maps f to g
    !
    !   The Idis of the maximum in C i_max,j_max,k_max correspond to the euler anglees beta,alpha,gamma in this order.
    !   Somehow there is still a bug. The resulting euler angles need to be modified as follows to yield correct results
    !   alpha,beta,gamma -----> 2*pi - alpha, beta , 2*pi-gamma

    !f2py threadsafe
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lms(:,:),g_lms(:,:)
    complex(kind = dp), intent(inout) :: cc(:,:,:)
    real(kind = dp), intent(in) :: radial_sampling_points(:)
    integer(kind = idp), intent(in) :: radial_limits(:)
    logical, intent(in) :: use_mp
    complex(kind = dp), allocatable :: tmp_corr(:,:,:),fft_c2c_in(:,:,:)
    logical :: dimensions_wrong
    real(kind = dp) :: inv_radial_range,radial_step
    integer(kind=dp) :: rid
    
    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    
    cc = 0._dp
    radial_step = radial_sampling_points(2)-radial_sampling_points(1)
    inv_radial_range = 1._dp/(radial_sampling_points(radial_limits(2))-radial_sampling_points(radial_limits(1)) + radial_step)

    !make sure array dimensions are correct
    dimensions_wrong = (size(radial_sampling_points,1)/=size(f_lms,2)) .OR. (size(radial_sampling_points,1)/=size(g_lms,2))
    if (dimensions_wrong) then
       print *, "Size of last dimension of f_lms and g_lms have to be equal to length of radial_sampling points, but:"
       print *, "shape(f_lm) = ", shape(f_lms)
       print *, "shape(g_lm) = ", shape(g_lms)
       print *, "shape(radial_sampling_points) = ", shape(radial_sampling_points)
       ERROR STOP "Wrong input array dimensions!"
    end if
    ! make sure radial_limits stay within bounds
    if ((radial_limits(1)<1) .or. (radial_limits(2)>size(radial_sampling_points,1))) then
       print *,"radial_limits out of bounds for radial_sampling points of length",size(radial_sampling_points,1)
       ERROR STOP "radial_limits contains out of bounds indices"
    end if

    if (use_mp) then
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rid,tmp_corr,fft_c2c_in)
       !! Using automatic arrays for fft_c2c_in and tmp_corr caused
       !! segfault due to stack overflow for cc of shape (64,64,64)
       !! and 8 parallel processes. 
       allocate(tmp_corr(self%bw*2,self%bw*2,self%bw*2))
       allocate(fft_c2c_in(self%bw*2,self%bw*2,self%bw*2))
       !$OMP DO REDUCTION(+:cc) SCHEDULE(STATIC)
       do rid=radial_limits(1),radial_limits(2)
          call self%cross_correlation_ylm_cmplx_(f_lms(:,rid),g_lms(:,rid),tmp_corr,fft_c2c_in,.FALSE.)          
          cc = cc + tmp_corr*radial_sampling_points(rid)**2
       end do
       !$OMP END DO
       deallocate(tmp_corr)
       deallocate(fft_c2c_in)
       !$OMP END PARALLEL
    else
       do rid=radial_limits(1),radial_limits(2)
          !print *, rid,radial_limits(2)
          call self%cross_correlation_ylm_cmplx(f_lms(:,rid),g_lms(:,rid),self%fft_c2c_out,.FALSE.)
          cc = cc+self%fft_c2c_out*radial_sampling_points(rid)**2
       end do
    end if
    cc = cc*inv_radial_range
  end subroutine cross_correlation_ylm_cmplx_3d

  ! correlation wigner inverse real
  subroutine inverse_wigner_loop_body_corr_real_kostelec_alloc(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:)
    real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn(self%bw-m2)
    integer(kind=dp) :: s_ids(2),bw,bw2,pm2_slice(2),l_start,t_slice(2)

    sym_const_m2 = (-1._dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_idp
    bw = self%bw
    bw2 = 2_idp*bw
    pm2_slice = MLr_slice(m2,bw)

    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * CONJG(g_ml(pm2_slice(1):pm2_slice(2))) * sym_const_m1m2
    so3func(:,m1+1,m2+1) = matmul(cc_lmn,wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    if (m1/=m2) then
       !! m2,m1 !!
       s_ids = order_to_ids(m2,m1,bw)
       cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * CONJG(g_ml(pm1_slice(1):pm1_slice(2)))
       
       so3func(:,s_ids(1),s_ids(2))=matmul(cc_lmn,wig_mat)
    end if

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)) * sym_const_m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)    
  end subroutine inverse_wigner_loop_body_corr_real_kostelec_alloc
  subroutine inverse_wigner_loop_body_corr_real_kostelec(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:)
    real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
    real(kind = dp) :: wig_mat(self%bw-m2,2*self%bw)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn(self%bw-m2)
    integer(kind=dp) :: s_ids(2),bw,bw2,pm2_slice(2),l_start,dlml_id

    sym_const_m2 = (-1._dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_idp
    bw = self%bw
    bw2 = 2_idp*bw
    pm2_slice = MLr_slice(m2,bw)

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genwig_l2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id),.True.)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * CONJG(g_ml(pm2_slice(1):pm2_slice(2))) * sym_const_m1m2
    so3func(:,m1+1,m2+1) = matmul(cc_lmn,wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    if (m1/=m2) then
       !! m2,m1 !!
       s_ids = order_to_ids(m2,m1,bw)
       cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * CONJG(g_ml(pm1_slice(1):pm1_slice(2)))
       
       so3func(:,s_ids(1),s_ids(2))=matmul(cc_lmn,wig_mat)
    end if
    
    !! -m2,-m1 !!
    !s_ids = order_to_ids(m2,m1,bw)
    !cc_lmn = wig_norm(l_start:) * CONJG(f_ml(pm2_slice(1):pm2_slice(2))) * g_ml(pm1_slice(1):pm1_slice(2))
    !so3func(:,s_ids(1),s_ids(2))=CONJG(matmul(cc_lmn,wig_mat))

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0 .or. m2==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)) * sym_const_m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)    
  end subroutine inverse_wigner_loop_body_corr_real_kostelec
  subroutine inverse_wigner_trf_corr_real_kostelec(self,f_ml,g_ml,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: i,n,m,m1,m2,lmax,bw,pm1_slice(2),pm1_tmp(2),s_ids(2),s_ids_sym(2)
    real(kind=dp) :: sym_const_m1,sym_array(self%bw),wig_norm(self%bw)
    procedure(wigner_corr_real_interface),pointer :: loop_body

    bw = self%bw
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_corr_real_kostelec_alloc
    else
       loop_body => inverse_wigner_loop_body_corr_real_kostelec
    end if

    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    fft_array = 0._dp
    
    ! initiallizing some constants
    lmax = self%lmax
    do i=0,lmax
       sym_array(i+1) = (-1.0)**i
       wig_norm(i+1) = 2._dp*pi*SQRT(2._dp/real(2_idp*i+1_idp,kind = dp))
    end do

    if (use_mp) then
       
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(i,m,m1,m2,sym_const_m1,pm1_slice,pm1_tmp,n)
       
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2
       !$OMP DO
       do i=1,(lmax+1)*(lmax+2)/2
          call flat_to_triangular_index(m1,m2,i,lmax)
          sym_const_m1 = (-1.0)**m1
          pm1_slice = MLr_slice(m1, bw)
          pm1_slice(1) = pm1_slice(1) + m2-m1
          call loop_body(self,f_ml,g_ml,fft_array,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
       end do
       !$OMP END DO

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       !$OMP DO
       do m2=1,lmax
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          fft_array(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(fft_array(:,s_ids(1),s_ids(2)))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2       
       do m1=0, lmax
          sym_const_m1 = (-1.0)**m1
          pm1_slice= MLr_slice(m1,bw)
          pm1_tmp = pm1_slice
          !print * , pm1_tmp
          do m2=m1, lmax
             n = m2-m1
             pm1_tmp(1) = pm1_slice(1)+n
             call loop_body(self,f_ml,g_ml,fft_array,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_tmp)
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,lmax
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          fft_array(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(fft_array(:,s_ids(1),s_ids(2)))
       end do
    end if
  end subroutine inverse_wigner_trf_corr_real_kostelec
  subroutine inverse_wigner_loop_body_corr_real_risbo(self,f_ml,g_ml,so3func,dlmn,l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: l,m1,m2
    real(kind=dp),intent(in) :: wig_norm,sym_const_l,sym_const_m1
    real(kind=dp),intent(in) :: dlmn(:)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn
    integer(kind=dp) :: s_ids(2),bw,bw2,l_start,i,pm1,pm2


    sym_const_m2 = (-1._dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_idp
    bw = self%bw
    bw2 = 2_idp*bw
    pm1 = MLr(m1,l,bw)
    pm2 = MLr(m2,l,bw)
        
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    cc_lmn = wig_norm * f_ml(pm1) * CONJG(g_ml(pm2)) * sym_const_m1m2
    !print *, "wig_norm",wig_norm
    !print *, "sym_const_m1m2",sym_const_m1m2
    !print *, "f_ml",f_ml(pm1),pm1
    !print *, "g_ml",CONJG(g_ml(pm2)),pm2
    !print *, "cc_lmn",cc_lmn
    !$OMP SIMD 
    do i=1,bw2
       so3func(i,m1+1,m2+1) = so3func(i,m1+1,m2+1) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD

    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    if (m1/=m2) then
       !! m2,m1 !!
       s_ids = order_to_ids(m2,m1,bw)
       cc_lmn = wig_norm * f_ml(pm2) * CONJG(g_ml(pm1))
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) +cc_lmn*dlmn(i)
       end do
       !$OMP END SIMD
    end if
    

    !! branch for sgn(m1)!=sgn(m2)                                   !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(l+m1)*d_{m1,-m2}^l(\pi-\beta)      !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause:                                                   !!
    !!    a constant sign swap by (-1)^M , M=m1 or m2                !!
    !!    a l dependent sign swap by (-1)^l                          !!
    !!    an inversion of the \beta coordinate axis                  !!
    if (m1==0) return       ! prevents sign swaps on 0 ids which are already covered
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    cc_lmn = wig_norm * f_ml(pm1) * g_ml(pm2) * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm * f_ml(pm2) * g_ml(pm1) * sym_const_m1m2 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_idp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
  end subroutine inverse_wigner_loop_body_corr_real_risbo
  subroutine inverse_wigner_trf_corr_real_risbo(self,f_ml,g_ml,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = idp) :: mnid,m1,m2,l,bw,bw2,i,j,s_ids(2),s_ids_sym(2)
    real(kind=dp) :: sym_const_m1,sym_const_l,wig_norm,dl(2*self%bw,(self%bw*(self%bw+1))/2)

    bw = self%bw
    bw2 = 2_idp*bw
    
    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    fft_array = 0._dp

    if (use_mp) then       
       ! non-fft part of the SO(3) fourier transform
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_idp*l+1_idp,kind = dp))
          !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1) SHARED(fft_array,f_ml,g_ml,dl,sym_const_l,l,wig_norm)
          !$OMP DO
          do mnid=1_idp,((l+1_idp)*(l+2_idp))/2_idp
             call flat_to_triangular_index(m1,m2,mnid,l)
             sym_const_m1 = (-1.0)**m1 
             call inverse_wigner_loop_body_corr_real_risbo(self,f_ml,g_ml,fft_array,dl(:,mnid),l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
          end do
          !$OMP END DO
          !$OMP END PARALLEL
       end do
       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       !$OMP PARALLEL PRIVATE(m2,s_ids,s_ids_sym,i,j) SHARED(fft_array,bw)
       !$OMP DO
       do j=1,(self%lmax+1_idp)*bw2
          m2 = ((j-1_idp)/bw2)
          i = MOD(j,bw2)+1_idp
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          fft_array(i,s_ids_sym(1),s_ids_sym(2)) = CONJG(fft_array(i,s_ids(1),s_ids(2)))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2       
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_idp*l+1_idp,kind = dp))
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_idp)
                call inverse_wigner_loop_body_corr_real_risbo(self,f_ml,g_ml,fft_array,dl(:,mnid),l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
             end do
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,self%lmax
          s_ids = order_to_ids(0_idp,m2,bw)
          s_ids_sym = order_to_ids(0_idp,-m2,bw)
          !$OMP SIMD
          do i=1,bw2
             fft_array(i,s_ids_sym(1),s_ids_sym(2)) = CONJG(fft_array(i,s_ids(1),s_ids(2)))
          end do
          !$OMP END SIMD
       end do

    end if
  end subroutine inverse_wigner_trf_corr_real_risbo
  subroutine inverse_wigner_trf_corr_real(self,f_ml,g_ml,fft_array,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    if (self%recurrence_type == kostelec_recurrence) then
       call self%inverse_wigner_trf_corr_real_kostelec(f_ml,g_ml, fft_array, use_mp)
    else
       if (allocated(self%wigner_d_trsp)) then
          ! genwig_all_risbo_preallocated has been used compute wigners and store them in mnl order
          ! compatible with faster kostelec transfom
          call self%inverse_wigner_trf_corr_real_kostelec(f_ml,g_ml, fft_array, use_mp)
       else
          call self%inverse_wigner_trf_corr_real_risbo(f_ml,g_ml, fft_array, use_mp)
       end if
    end if
  end subroutine inverse_wigner_trf_corr_real
  ! cross correlations ylm complex
  subroutine cross_correlation_ylm_real_(self,f_ml,g_ml,cc,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    real(kind = dp), intent(inout) :: cc(:,:,:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp

    ! non FFT wigner part of the correlation computation.
    call self%inverse_wigner_trf_corr_real(f_ml,g_ml,fft_array,use_mp)
    ! Compute fft.
    fft_array = CONJG(fft_array) ! to correct for the fact that we have to compute the forward not the backward fft.
    call fftw_execute_dft_c2r(self%plan_c2r_backward_p,fft_array,cc)
    cc = cc * (1/(2._dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)   
  end subroutine cross_correlation_ylm_real_
  subroutine cross_correlation_ylm_real(self,f_ml,g_ml,cc,use_mp)
    ! let f,g be two square integrable functions on the 2 sphere 
    ! Define CC: SO(3) ---> \mathbb{R},   R |---> <f,g \circ R> = \int_{S^2} dx f(x)*\overline{g(Rx)}
    ! This function calculates CC(R) for all R defined in make_SO3_grid
    ! Spherical harmonic coefficient arrays , g_lm,f_lm, are assumed to follow the indexing conventions of the SHTNS package,
    ! i.e. those given by the pure functions LMc LMr from above.
    !
    ! arguments :
    !   f_lm: f_{l,m} spherical harmonic coefficients of f 
    !   g_lm: g_{l,m} spherical harmonic coefficients of g 
    !            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
    !
    !
    ! output :
    !   C_values: array of shape (2*bw,2*bw,2*bw)
    !   The maximum in C corresponds to the Rotation that best maps f to g
    !
    !   The Idis of the maximum in C i_max,j_max,k_max correspond to the euler anglees beta,alpha,gamma in this order.
    !   Somehow there is still a bug. The resulting euler angles need to be modified as follows to yield correct results
    !   alpha,beta,gamma -----> 2*pi - alpha, beta , 2*pi-gamma
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    real(kind = dp), intent(inout) :: cc(:,:,:)
    logical, intent(in) :: use_mp

    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    call self%cross_correlation_ylm_real_(f_ml,g_ml,cc,self%fft_c2r_in,use_mp)
  end subroutine cross_correlation_ylm_real
  subroutine cross_correlation_ylm_real_3d(self,f_mls,g_mls,cc,radial_sampling_points,radial_limits,use_mp)
    ! let f,g be two square integrable functions on the $\mathbb{R}^3$ in spherical coordinates 
    ! Define CC: SO(3) ---> \mathbb{R},   R |---> <f,g \circ R> = \int_R dr r^2 \int_{S^2} dw f(r,w)*\overline{g(r,Rw)}
    ! This function calculates CC(R) for all R defined in make_SO3_grid
    ! Spherical harmonic coefficient arrays , g_lm,f_lm, are assumed to follow the indexing conventions of the SHTNS package,
    ! i.e. those given by the pure functions LMc LMr from above.
    !
    ! arguments :
    !   f_lm: f_{l,m} spherical harmonic coefficients of f 
    !   g_lm: g_{l,m} spherical harmonic coefficients of g 
    !            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
    !
    !
    ! output :
    !   C_values: array of shape (2*bw,2*bw,2*bw)
    !   The maximum in C corresponds to the Rotation that best maps f to g
    !
    !   The Idis of the maximum in C i_max,j_max,k_max correspond to the euler anglees beta,alpha,gamma in this order.
    !   Somehow there is still a bug. The resulting euler angles need to be modified as follows to yield correct results
    !   alpha,beta,gamma -----> 2*pi - alpha, beta , 2*pi-gamma

    !f2py threadsafe
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_mls(:,:),g_mls(:,:)
    real(kind = dp), intent(inout) :: cc(:,:,:)
    real(kind = dp), intent(in) :: radial_sampling_points(:)
    integer(kind = idp), intent(in) :: radial_limits(:)
    logical, intent(in) :: use_mp
    real(kind = dp), allocatable :: tmp_corr(:,:,:)
    complex(kind = dp), allocatable :: fft_c2r_in(:,:,:)
    logical :: dimensions_wrong
    real(kind = dp) :: inv_radial_range,radial_step
    integer(kind=dp) :: rid


    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    
    cc = 0._dp
    radial_step = radial_sampling_points(2)-radial_sampling_points(1)
    inv_radial_range = 1._dp/(radial_sampling_points(radial_limits(2))-radial_sampling_points(radial_limits(1)) + radial_step)
    !print *, "step",radial_step
    !print *, "inv_range",inv_radial_range

    !make sure array dimensions are correct
    dimensions_wrong = (size(radial_sampling_points,1)/=size(f_mls,2)) .OR. (size(radial_sampling_points,1)/=size(g_mls,2))
    if (dimensions_wrong) then
       print *, "Size of last dimension of f_lms and g_lms have to be equal to length of radial_sampling points, but:"
       print *, "shape(f_lm) = ", shape(f_mls)
       print *, "shape(g_lm) = ", shape(g_mls)
       print *, "shape(radial_sampling_points) = ", shape(radial_sampling_points)
       ERROR STOP "Wrong input array dimensions!"
    end if
    ! make sure radial_limits stay within bounds
    if ((radial_limits(1)<1) .or. (radial_limits(2)>size(radial_sampling_points,1))) then
       print *,"radial_limits out of bounds for radial_sampling points of length",size(radial_sampling_points,1)
       ERROR STOP "radial_limits contains out of bounds indices"
    end if
    
    if (use_mp) then
       !$OMP PARALLEL DEFAULT(SHARED) PRIVATE(rid,tmp_corr,fft_c2r_in)
       !! Using automatic arrays for fft_c2r_out and tmp_corr caused
       !! segfault due to stack overflow for cc of shape (64,64,64)
       !! and 8 parallel processes. 
       allocate(tmp_corr(self%bw*2,self%bw*2,self%bw*2))
       allocate(fft_c2r_in(self%bw*2,self%bw+1,self%bw*2))
       !$OMP DO REDUCTION(+:cc) SCHEDULE(STATIC)
       do rid=radial_limits(1),radial_limits(2)
          call self%cross_correlation_ylm_real_(f_mls(:,rid),g_mls(:,rid),tmp_corr,fft_c2r_in,.FALSE.)
          cc = cc + tmp_corr*radial_sampling_points(rid)**2
       end do
       !$OMP END DO
       deallocate(tmp_corr)
       deallocate(fft_c2r_in)
       !$OMP END PARALLEL
    else
       do rid=radial_limits(1),radial_limits(2)
          !print *, 'rid',rid
          !print *, 'f_mls',f_mls(:,rid)
          !print *, 'g_mls',g_mls(:,rid)
          !print *, 'fft_array', SHAPE(self%fft_r2c_out)
          !print *, 'fft_array', radial_sampling_points(rid)**2
          call self%cross_correlation_ylm_real(f_mls(:,rid),g_mls(:,rid),self%fft_r2c_out,.FALSE.)
          cc = cc+self%fft_r2c_out*radial_sampling_points(rid)**2
       end do
    end if
    cc = cc*inv_radial_range
  end subroutine cross_correlation_ylm_real_3d

  ! Standalone 2D ffts used in SO(3)FFT
  subroutine fft(self,f1,f2)
    class(so3ft),intent(in) :: self
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft(self%plan_c2c_forward,f1,f2)
    f2 = f2 * (1._dp/real(Size(f1,1),kind=dp)) 
  end subroutine fft
  subroutine ifft(self,f1,f2)
    class(so3ft),intent(in) :: self
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft(self%plan_c2c_backward,f1,f2)
    f2 = f2 * (1._dp/real(Size(f1,1),kind=dp)) 
  end subroutine ifft
  subroutine rfft(self,f1,f2)
    class(so3ft),intent(in) :: self
    real(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft_r2c(self%plan_r2c_forward,f1,f2)
    f2 = f2 * (1._dp/real(Size(f1,1),kind=dp)) 
  end subroutine rfft
  subroutine irfft(self,f1,f2)
    class(so3ft),intent(in) :: self
    complex(kind = dp), intent(in) :: f1(:,:,:)
    real(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft_c2r(self%plan_c2r_backward,f1,f2)
    f2 = f2 * (1._dp/(real(Size(f1,1),kind=dp)))
  end subroutine irfft
  
  ! Class functions
  function get_so3func_part_halfcomplex(bw,m1,m2,so3func) result(so3func_part)
    complex(kind = dp),intent(in) :: so3func(:,:,:)
    integer(kind = idp),intent(in) :: m1,m2,bw
    complex(kind = dp) :: so3func_part(size(so3func,1))
    integer(kind = idp) :: s_ids(2)

    if (m1<0) then
       s_ids = order_to_ids(-m1,-m2,bw)
       so3func_part = CONJG(so3func(:,s_ids(1),s_ids(2)))
    else
       s_ids = order_to_ids(m1,m2,bw)
       so3func_part = so3func(:,s_ids(1),s_ids(2))
    end if
       
  end function get_so3func_part_halfcomplex
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
  
end module softclass
