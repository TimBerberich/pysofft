!> --------
!! @brief Module containing code to compute Wigner-D matrices 
module make_wigner
  use precision
  use math_constants
  use utils
  implicit none
contains
  !> ------
  !! @brief Number of all non-symmetric Wigner-d matrix elements, $l<\\mathrm{bw}$.
  !!
  !! $\\frac{\\mathrm{bw}^2(2+3\\mathrm{bw}+\\mathrm{bw}^2)}{3}$
  function size_wigner_d(bw) result(wsize)
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp) :: wsize
    wsize = bw*bw * int(2_idp + 3_idp*bw + bw*bw,idp)/3_idp
  end function size_wigner_d
  
  !> ------
  !! @brief Shape of all non-symmetric wigner-d matrix elements, $l<\\mathrm{bw}$.
  !!
  !! $\\left[2\\,\\mathrm{bw},\\frac{\\mathrm{bw}(2+3\\mathrm{bw}+\\mathrm{bw}^2)}{6}\\right]$
  function wigner_d_shape(bw) result(d_shape)
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp) :: d_shape(2)
    d_shape = [2_idp*bw,(bw*(2_idp + 3_idp*bw + bw*bw))/6_idp]
  end function wigner_d_shape

  !> ------
  !! @brief Returns n uniformly sampled angles in $(0,\\pi)$
  !!
  !! $\frac{(2 i + 1) \pi}{2 n}$ for $i=0,\\ldots,n-1$
  !! Note: These are angles used in Chebyshev nodes of the first kind.
  function create_beta_samples(n) result(betas)
    integer(kind = idp) :: n,i
    real(kind = dp) :: betas(n),factor
    factor = pi / real(2_idp * n, dp)
    betas = [(real(2_idp * i - 1_idp, dp) * factor, i = 1, n)]
  end function create_beta_samples

  !> -------
  !! @brief Returns trig samples for Wigner-d computations via kostelec.
  !! 
  !! Does the following:
  !! ```fortran
  !!  betas = create_beta_samples(2*bw)
  !!  trig_samples(:,1) = cos(betas)                      ! <= Chebyshev nodes
  !!  trig_samples(:,2) = cos(betas/2_idp)*sin(betas/2_idp) ! <= needed for d^l_ml computation
  !!  trig_samples(:,3) = cos(betas/2_idp)**2              ! <= needed for d^l_ml computation
  !! ```
  function create_trig_samples(bw) result(trig_samples)
    integer(kind=dp) :: bw
    real(kind = dp) :: trig_samples(2*bw,3),betas(2*bw)
    betas = create_beta_samples(2*bw)
    trig_samples(:,1) = cos(betas)                      ! <= Chebyshev nodes
    trig_samples(:,2) = cos(betas/2_idp)*sin(betas/2_idp) ! <= needed for d^l_ml computation
    trig_samples(:,3) = cos(betas/2_idp)**2              ! <= needed for d^l_ml computation
  end function create_trig_samples

  !> -------
  !! @brief Returns trig samples for Wigner-d computations via risbo.
  !! 
  !! Does the following:
  !! ```fortran
  !!  betas = create_beta_samples(2*bw)
  !!  trig_samples(:,1) = cos(betas/2)
  !!  trig_samples(:,2) = sin(betas/2)
  !! ```
  function create_trig_samples_risbo(bw) result(trig_samples)
    integer(kind=dp) :: bw
    real(kind = dp) :: trig_samples(2*bw,2),betas(2*bw)
    betas = create_beta_samples(2*bw)
    trig_samples(:,1) = cos(betas/2)
    trig_samples(:,2) = sin(betas/2)
  end function create_trig_samples_risbo

  !> -------
  !! @brief Returns square root samples for Wigner-d computations via risbo.
  !!
  !! Returns $\sqrt{l}$ for $0 \leq l \leq (2\mathrm{bw}-2)$.
  function create_sqrts_risbo(bw) result(sqrts)
    integer(kind=dp), intent(in) :: bw
    real(kind = dp) :: sqrts(2*bw-1)
    integer(kind=dp) :: i
    do i=0,2*bw-2
       sqrts(i+1) = SQRT(real(i,dp))
    end do
  end function create_sqrts_risbo
  
  !> ------
  !! @brief Returns n uniformly sampled angles in $[0,2 \\pi)$
  !!
  !! $\frac{i 2 \\pi}{n}$ for $i=0,\ldots,n-1$
  function create_alpha_gamma_samples(n) result(alpha)
    ! Note: These are angles nodes of an FFT
    integer(kind = idp) :: n,i
    real(kind = dp) :: alpha(n),factor
    factor = 2._dp*pi / real(n, dp)
    alpha = [(real(i-1_idp, dp) * factor, i = 1, n)]
  end function create_alpha_gamma_samples

  !> -----
  !! @brief Computes $d^l\\_{m\\_1,l}(\beta)$ via recursion over $m\\_1$
  !!
  !! The analytic formula for the normalized d^l\\_{m\\_1,l}(\beta) is
  !!
  !! $$
  !!  {d^l\\_{m,l}}(\\beta) = \\sqrt{\\frac{2l+1}{2}} \\sqrt{\\frac{2l!}{(l+m)!(l-m)!}} \\cos(\\frac{\\beta}{2})^{l+m} \\sin(\\frac{\\beta}{2})^{l-m}
  !! $$
  !!
  !! This function computes $d^l\\_{m\\_1,l}(\beta)$ at fixed $m\\_1$ for all orders $l$ based on $d^l\\_{m\\_1-1,l}(\beta)$, using the following recursion relations
  !! 
  !! $$
  !!    \begin{aligned}
  !!      {d^{l+1}\\_{m,l+1}}(\\beta) &= {d^l\\_{m-1,l}}\\sqrt{\\frac{(2l+2)(2l+3)}{(l+m)(l+m+1)}} \\cos(\\frac{\\beta}{2})^2 \\\\
  !!      {d^{l+1}\\_{m,l+1}}(\\beta) &= {d^l\\_{m,l}}\\sqrt{\\frac{(2l+2)(2l+3)}{(l+m+1)(l-m+2)}} \\cos(\\frac{\\beta}{2}) \\sin(\\frac{\\beta}{2})
  !!    \end{aligned}
  !! $$
  !! 
  !! The recursion can be started at $m\\_1=0$ by providing teh input array `dlml` with `dlml(1) = SQRT(0.5_idp)` (i.e.: $d^0\\_{0,0}(\beta) = \\sqrt{\frac{1}{2}}$ ). 
  !!
  !! 
  !! /// info|Recursion scheme\n
  !! ![dlml\\_l]{../images/dlml.svg#only-light}(width=400 align=left)\n
  !! ![dlml\\_l]{../images/dlml_white.svg#only-dark}(width=400 align=left)\n
  !! The diagram shows the recursion scheme for `bw=6`. Diagonal arrows correspond to usage of the first recursion relation, i.e. $d^l\\_{m_1-1,l}\rightarrow d^{l+1}\\_{m_1,l}$ , and horizontal steps imply usage of the second recursion relation, i.e. $d^l\\_{m_1,l}\rightarrow d^{l+1}\\_{m_1,l}$.\n
  !! ///
  !!
  !! For a detailed usage example see the source of [compute_all_dlml_l_contiguous](namespacemake__wigner.md#function-compute_all_dlml_l_contiguous).
  subroutine dlml_recursion_l_contiguous(dlml,m,bw,sincos,cos2,normalized)
    !! dlml for all l at fixed m. Has to be of size [Size(sincos,1),bw]
    real(kind = dp),intent(inout) :: dlml(:,:)
    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: m,bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    integer(kind = idp) :: i,l,swap_id,swap_id2,nc1,nc2

    !Normalization defining constants
    nc1 = MERGE(2_idp,1_idp,normalized)
    nc2 = nc1+1_idp
    
    swap_id = m+1_idp
    swap_id2 = Merge(m,swap_id,(m==0))
    do i=1,min(bw-m,swap_id2)
       !(l,m-1 -> l+1,m)
       l=m+(i-1_idp)-1_idp
       dlml(:,i) = dlml(:,i) * SQRT(real((2_idp*l+nc1)*(2_idp*l+nc2),kind=dp)/real((l+m)*(l+m+1_idp),kind=dp))*cos2
    end do
    do i=swap_id,bw-m-1_idp
       !(l -> l+1)
       l=m+(i-1_idp)
       dlml(:,i+1) = dlml(:,i) * SQRT(real((2_idp*l+nc1)*(2_idp*l+nc2),kind=dp)/real((l+m+1_idp)*(l-m+1_idp),kind=dp))*sincos
    end do
  end subroutine dlml_recursion_l_contiguous
  
  !> ----------
  !! @brief Computes all $d^l\\_{m\\_1,m\\_2}(\beta)$ for $m\\_2=l<\\mathrm{bw}$.
  !!
  !! Computes all Wigner little d coefficients of the form d_lml(beta) with l<bw and stores the values in an l contiguous way.
  !! That is $d^l\\_{m,l}$ is stored at index $=$ triangular_to_flat_index(m,l,bw)
  function compute_all_dlml_l_contiguous(bw,sincos,cos2,normalized) result(dlml)
    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(SIZE(cos2,1),(bw*(bw+1))/2_idp),dlml_tmp(Size(cos2,1),bw)
    integer(kind=dp) :: m,lm_id,next_lm_id
    
    ! Compute dlml values for m=0
    dlml(:,1) = MERGE(SQRT(0.5_dp),1._dp,normalized)     ! In agreement with our normalization $\sqrt{\frac{2l+1}{2}}$ at l=0
    ! Setup workspace
    dlml_tmp(:,1) = dlml(:,1)
    
    ! Do stair recursion alternating between (l -> l+1) and (l,m -> l+1,m+1 ) steps.
    do m=0,bw-1
       call dlml_recursion_l_contiguous(dlml_tmp, m, bw, sincos, cos2,normalized)
       lm_id = triangular_to_flat_index(m,m,bw)
       next_lm_id = triangular_to_flat_index(m,bw-1,bw)
       dlml(:,lm_id:next_lm_id) = dlml_tmp(:,:bw-m)
    end do
  end function compute_all_dlml_l_contiguous

  !> -----
  !! @brief Compute $d^l\\_{m,l}$ for specific $l,m$.
  !!
  !! This is a convenience function that computes the small wigner-d matrix elements $d^l_{m,l}$ for specific $l,m$.
  !! It does so using the l-contiguous recurrence relation `dlml_recursion_l_contuguous`.
  !! DO NOT use for perfomance critical computations!  
  function compute_dlml(l,m,sincos,cos2,normalized) result(dlml)
    !! Indices $l,m$ at which to compute d^l_{m,l}
    integer(kind = idp),intent(in) :: l,m
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(size(cos2,1)), dlml_tmp(size(cos2,1),l+1)
    integer(kind = idp) :: mm

    dlml_tmp = MERGE(SQRT(0.5_dp),1._dp,normalized) ! In agreement with our normalization $\sqrt{\frac{2l+1}{2}}$ at l=0
    
    ! Do stair recursion alternating between (l -> l+1) and (l,m -> l+1,m+1 ) steps.
    do mm=0,m
       call dlml_recursion_l_contiguous(dlml_tmp, mm, l+1, sincos, cos2,normalized)
    end do
    dlml = dlml_tmp(:,l+1-m)
  end function compute_dlml
  
  !> --------
  !! @brief three term recurrence relation
  !!
  !! This function implements the three term recurrence relation used by Kostelec
  !! For a usage example usage see the source code of [genwig_l2](namespacemake__wigner.md#function-genwig_l2).
  !!
  !! @param[in,out]   workspace   (SIZE(beta,1),3), storing the elements for the three term recurrence
  !! @param[in]   cos   (SIZE(beta,1),), containing $\cos(\beta)$ values
  !! @param[in]   l  degree to compute.
  !! @param[in]   m1  first order
  !! @param[in]   m2  second order
  !! @param[in]   normalized  Whether or not to compute the L2 normalized wigners.
  subroutine wig_l_recurrence_kostelec(workspace,cos,l,m1,m2,normalized)
    !
    integer(kind = idp), intent(in) :: l,m1,m2
    real(kind = dp),intent(in) :: cos(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp),intent(inout) :: workspace(:,:)
    integer(kind = idp) :: dl_id,dl_1id,o,i
    real(kind = dp) :: c_common,c_l,c_l_1
    
    i = l-m2+1_idp
    dl_1id = MODULO(i-1_idp,3)+1_idp
    dl_id = MODULO(i,3)+1_idp
    o = MODULO(i+1_idp,3)+1_idp
    
       
    c_common = real(l+1_idp,kind = dp)*SQRT(1._dp/real( ((l+1_idp)**2-m1**2)*((l+1_idp)**2-m2**2),kind = dp))
    c_common = MERGE(c_common*SQRT( real(2_idp*l+3_idp,kind=dp)),c_common,normalized)
    
    c_l = real(2_idp*l+1_idp,kind = dp)
    c_l = MERGE(SQRT(c_l),c_l,normalized)

    if (l==0) then
       workspace(:,o) = c_common*c_l*cos*workspace(:,dl_id)
    else if (l==m2) then
       workspace(:,o) = c_common*c_l*( cos-(real(m1*m2,kind=dp)/real(l*(l+1_idp),kind = dp)) )*workspace(:,dl_id)
    else
       workspace(:,o) = c_common*c_l*( cos-(real(m1*m2,kind=dp)/real(l*(l+1_idp),kind = dp)) )*workspace(:,dl_id)          
       c_l_1 = SQRT(real((l**2-m1**2)*(l**2-m2**2),kind = dp))/real(l,kind = dp)
       c_l_1 = MERGE(c_l_1/SQRT(real(2_idp*l-1_idp,kind = dp)),c_l_1,normalized)
       workspace(:,o) = workspace(:,o) - c_common*c_l_1*workspace(:,dl_1id)
    end if
    
  end subroutine wig_l_recurrence_kostelec

  !> --------------
  !!  @brief Computes $d^l\\_{m\\_1,m\\_2}(\beta)$ for $m\\_2 \\leq l \\leq \\mathrm{bw}$ using kostelec recursion.
  !!
  !!  Given orders 0<=m1<=m2<=bw, and a bandwidth bw, this function will
  !!  generate all the Wigner little d functions whose orders
  !!  are (m1, m2) and degrees are j = max(|m1|, |m2|) = m2 through j = bw - 1
  !!  using the 3-term recurrence. 
  !!
  !!  When normalization = True is used all of these Wigners
  !!  will have L2 norm = 1
  !!  
  !!  arguments: m1, m2 = orders of the functions
  !!             bw = bandwidth
  !!             trig_samples = array containing cos(beta),cos(beta/2) and sine(beta) values
  !!               
  !!             result = array to store result, length (bw-m)*n   (as (n,bw-m) matrix)
  !!                      where m = max( |m1|, |m2| ); 
  !!  
  !!  The routine won't be efficient, but it'll work.
  !!  NOTE: compared to the python/c version the arrays are transposed:
  !!        As mentioned before, here we compute the small wigner d matrices as
  !!        (n,bw-m) where as in the C/python version they have the shape (bw-m,n)
  function genwig_l2(m1,m2,bw,cos,dlml,normalized) result(wigners_m1m2)
    integer(kind = idp),intent(in) :: m1,m2,bw
    real(kind = dp),intent(in) :: cos(:),dlml(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp) :: wigners_m1m2(Size(cos,1),(bw-max(abs(m1),abs(m2))))
    
    real(kind = dp) :: workspace(Size(cos,1),3)
    integer(kind = idp) :: i,l_start,l,n_samples,o

    if (.NOT.(0<=m1 .AND. m1<=m2 .AND. m2<=bw )) then
       print *, 'Invalid arguments:  0<=m1<=m2<=bw is not sattisfied.'
    end if
    
    l_start = m2 !max(abs(m1),abs(m2))
    n_samples = size(cos,1)
    workspace(:,1) = 0._dp
    workspace(:,2) = dlml
    wigners_m1m2(:,1) = workspace(:,2)
    
    !print * , 'lstart', l_start,'m1m2',m1,m2
    do i=1, bw-l_start-1
       l=l_start+(i-1_idp)
       call wig_l_recurrence_kostelec(workspace,cos,l,m1,m2,normalized)
       o = MODULO(i+1_idp,3)+1_idp
       wigners_m1m2(:,i+1) = workspace(:,o)
    end do
  end function genwig_l2

  !> ---------
  !! @brief Same as genwig_l2 but transposed 
  !!
  !! This routine is faster than calling genwig_l2 and then transposing.
  function genwig_l2_trsp(m1,m2,bw,cos,dlml,normalized) result(wigners_m1m2)
    integer(kind = idp), intent(in) :: m1,m2,bw
    real(kind = dp), intent(in) :: cos(:),dlml(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp) :: wigners_m1m2((bw-max(abs(m1),abs(m2))),Size(cos,1))
    
    real(kind = dp) :: workspace(Size(cos,1),3)
    integer(kind = idp) :: i,l_start,l,n_samples,o

    
    if (.NOT.(0<=m1 .AND. m1<=m2 .AND. m2<=bw )) then
       print *, 'Invalid arguments:  0<=m1<=m2<=bw is not sattisfied.'
    end if
    
    l_start = m2 !max(abs(m1),abs(m2))
    n_samples = size(cos,1)
    workspace(:,1) = 0._dp
    workspace(:,2) = dlml
    wigners_m1m2(1,:) = workspace(:,2)
    
    !print * , 'lstart', l_start,'m1m2',m1,m2
    do i=1, bw-l_start-1
       l=l_start+(i-1_idp)
       call wig_l_recurrence_kostelec(workspace,cos,l,m1,m2,normalized)
       o = MODULO(i+1_idp,3)+1_idp
       wigners_m1m2(i+1,:) = workspace(:,o)
    end do
  end function genwig_l2_trsp

  !> ---------
  !! @brief Computes the small Wigner d matrix $d^l\\_{m,n}(\beta)$ for fixed $l$.
  !!
  !! DO NOT use for performace criticlal tasks!
  !!
  !! Convenience function to compute the small Wigner d matrix $d^l\\_{m,n}(\beta)$ for fixed $l$.
  !! It uses the three recurrence realtion from `wig_l_recurrence`.
  !! This implies that it has to compute ALL small wigner $d^k\\_{m,n}$ matrices with $k<l$ as well.
  !!
  !! IMPORTANT: Pi-betas must be equal to reversing the betas array.
  !!            Otherwise the routine will not give correct results
  function wigner_dl_kostelec(l,betas,normalized) result(dl)
    !! compute the small wigner-d matrix $d^l_{m_1,m_2}(\beta)$ for fixed l
    integer(kind = idp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp) :: dl(SIZE(betas,1),2*l+1,2*l+1)
    real(kind = dp) :: trig(SIZE(betas,1),3)
    real(kind = dp) :: workspace(SIZE(betas,1),3)
    real(kind=dp) :: dlml_workspace(SIZE(betas,1),l+1),sym_const_m,sym_const_n,sym_const_l
    integer(kind = idp) :: m,n,i,mid,nid,neg_mid,neg_nid,o,bw,nbeta

    dl = 0
    sym_const_l = (-1._dp)**l
    trig(:,1) = cos(betas)                      
    trig(:,2) = cos(betas/2_idp)*sin(betas/2_idp) 
    trig(:,3) = cos(betas/2_idp)**2
    nbeta = SIZE(betas,1)
    bw = l+1_idp

    dlml_workspace(:,1) = MERGE(SQRT(0.5_dp),1._dp,normalized)
    do m=0,l
       mid = l+1_idp+m
       neg_mid = l+1_idp-m
       sym_const_m = (-1._dp)**m
       call dlml_recursion_l_contiguous(dlml_workspace, m, bw, trig(:,2), trig(:,3),normalized)
       do n=m,l
          nid = l+1_idp+n
          neg_nid = l+1_idp-n
          sym_const_n = (-1._dp)**n
          !norm = SQRT(2._dp/real(2_idp*l+1_idp,kind=dp))
          
          workspace(:,1) = 0
          workspace(:,2) = dlml_workspace(:,1+n-m)
          workspace(:,3) = workspace(:,2)
          
          do i=n,l-1
             call wig_l_recurrence_kostelec(workspace,trig(:,1),i,m,n,normalized)
          end do
          
          o = MODULO(l-n+1_idp,3)+1_idp

          ! populate the small wigner d-matrix using all of the available symmetries
          dl(:,nid,mid)=workspace(:,o)          
          if ((m==0 .AND. n==0)) cycle
          dl(:,neg_mid,neg_nid) = workspace(:,o)
          if (.NOT. m==n) then
             dl(:,mid,nid) = (sym_const_m*sym_const_n)*workspace(:,o)
             dl(:,neg_nid,neg_mid) = (sym_const_m*sym_const_n)*workspace(:,o)
          end if
          if (m==0 .or. n==0) cycle
          dl(:,neg_nid,mid)=workspace(nbeta:1:-1,o)*(sym_const_l*sym_const_m)
          dl(:,nid,neg_mid)=workspace(nbeta:1:-1,o)*(sym_const_l*sym_const_n)
          if (m==n) cycle
          dl(:,neg_mid,nid)=workspace(nbeta:1:-1,o)*(sym_const_l*sym_const_m)
          dl(:,mid,neg_nid)=workspace(nbeta:1:-1,o)*(sym_const_l*sym_const_n)
          
       end do
    end do
  end function wigner_dl_kostelec
  
  !> -----
  !! @brief Same as genwig_all but writes into a given array
  !!
  !! Computes all non-symmetry equivalent Wigner-d matrix values and stores them in the input array `wigners`.
  subroutine genwig_all_kostelec_preallocated(bw,wigners,normalized)
    ! Same as genwig_all but writes the wigner d matricies into the provided wigners array
    ! instead of creating a new array
    !
    ! Computes all independent small wigner d-matrix elements
    ! upto a n order of l=bw-1.
    !
    ! Important: Different from The original C implementation and
    ! the previouse numba implementation within pysofft, this routine
    ! does compute d_(0,m) values instead of d(m,0) in the old versions.
    ! This allows one to loop over all values of m1,m2 in a simple double
    ! loop ,i.e.
    ! do m1=0, bw-1
    !   do m2=m1, bw-1
    ! instead of looping over all cases where m2=0 or m1=m2 separately.
    
    integer(kind = idp), intent(in) :: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp), intent(inout) :: wigners(:,:)
    real(kind = dp) :: trig_samples(2*bw,3)
    real(kind=dp) :: dlml_workspace(2_idp*bw,bw)
    ! wigner should have size total_n_wigners(bw) but seems like f2py does not support array size asignments trhough pure functions ...
    integer(kind = idp) :: m1,m2,slize(2)
    
    if (.NOT. size(wigners)==size_wigner_d(bw)) then
       print '(A, I3)','Fortran Error: Provided wigners array has wrong length, it shoud be of shape:', size_wigner_d(bw) 
       stop
    end if
    
    wigners = 0._dp
    dlml_workspace(:,1) = MERGE(SQRT(0.5_dp),1._dp,normalized)
    trig_samples = create_trig_samples(bw)    
    do m1=0,bw-1
       call dlml_recursion_l_contiguous(dlml_workspace, m1, bw, trig_samples(:,2), trig_samples(:,3),normalized)
       do m2=m1, bw-1
          slize = mnl_to_flat_l_slice(m1,m2,bw)
          !print*, m1,m2
          !print*, genwig_l2(m1,m2,bw,trig_samples)
          wigners(:,slize(1):slize(2)) = genwig_l2(m1,m2,bw,trig_samples(:,1),dlml_workspace(:,m2-m1+1),normalized)
       end do
    end do
  end subroutine genwig_all_kostelec_preallocated

  !> -------
  !! @brief Computes all $d^l\\_{m,n}$ for $0 \\leq m \\leq n \leq l < bw$.
  !! 
  !! Computes all independent small Wigner d-matrix elements d^l\\_{m,n}(\beta), i.e. those  
  !! for $0 \\leq m \\leq n \\leq l <bw$.
  !! This routine uses the three term recurrence `wig_l_recurrence_kostelec`.
  function genwig_all_kostelec(bw,normalized) result(wigners)
    ! Computes all independent small wigner d-matrix elements
    ! upto a n order of l=bw-1.
    !
    ! Important: Different from The original C implementation and
    ! the previouse numba implementation within pysofft, this routine
    ! does compute d_(0,m) values instead of d(m,0) in the old versions.
    ! This allows one to loop over all values of m1,m2 in a simple double
    ! loop ,i.e.
    ! do m1=0, bw-1
    !   do m2=m1, bw-1
    ! instead of looping over all cases where m2=0 or m1=m2 separately.
    
    integer(kind = idp),intent(in):: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp) :: wigners(bw*2_idp,(bw * (2_idp + 3_idp*bw + bw*bw))/6_idp)
    call genwig_all_kostelec_preallocated(bw,wigners,normalized)
  end function genwig_all_kostelec

  !> ----------
  !!
  !! @brief $d^{l-1}\\_{m,n}(\beta) \rightarrow d^l\\_{m,n}(\beta)$
  !!
  !! Given a wigner little-d matrix of degree L-1 evaluated
  !! at some angles beta, this function constructs the wigner little-d
  !! matrix of degree L, evaluated at the same betas.
  !! This algorithm is an implementation of the one given by T. Risbo in
  !!
  !! >    "Fourier transform summation of Legendre series and D-Functions"
  !!      T. Risbo
  !!      Journal of Geodesy
  !!      1996
  !!      volume 70: p. 383 - 396
  !!
  !! Example: To compute the wigner d-matrix of degree l you can use this reccurrence
  !!     as follows
  !! ```fortran
  !! real(kind=8) :: betas(1)
  !! real(kind=8) :: dl(Size(betas,1),1)
  !! real(kind=8) :: ls_sqrt(l+1)
  !! ! Setup
  !! betas(1) = 1.23
  !! dl = 0._dp
  !! cos_betas=COS(betas/2)
  !! sin_betas=SIN(betas/2)
  !! do l=0,2*l+1
  !!   ls_sqrt(l+1) = SQRT(l)
  !! end do
  !!
  !! ! Use reccurrence
  !! do l=1,L
  !!   dl = wigner_d_recurrence(dl,l,cos_beta,sin_beta,ls_sqrt)
  !! end do
  !! ```
  !!
  function wigner_mn_recurrence_risbo(d_l_1,l,cos_beta,sin_beta,sqrts) result(d_l)
    real(kind = dp), intent(in) :: cos_beta(:),sin_beta(:),sqrts(:),d_l_1(:,:)
    integer(kind = idp), intent(in)  :: l
    real(kind = dp) :: d_l(Size(d_l_1,1),(2*l+1)*(2*l+1)),temp(Size(d_l_1,1),(2*l)*(2*l)),inv_deg
    integer(kind = idp) :: deg,i,j,ij,ippj,ijpp,ippjpp,pdeg
    
    if (l==0) then
       d_l = 1._dp
    else if (l==1) then
       d_l(:,1) = cos_beta**2
       d_l(:,2) = sqrts(3) * cos_beta * sin_beta
       d_l(:,3) = sin_beta**2

       d_l(:,4) = -d_l(:,2)
       d_l(:,5) =  d_l(:,1)-d_l(:,3)
       d_l(:,6) =  d_l(:,2)

       d_l(:,7) =  d_l(:,3)
       d_l(:,8) = -d_l(:,2)
       d_l(:,9) =  d_l(:,1)

    else
       temp = 0._dp
       deg = 2*l-1
       pdeg = deg+1
       inv_deg = (1._dp/real(deg,kind=dp))
       do i=0,deg-1
          do j=0,deg-1
             ij     = i*pdeg + j+1
             ippj   = (i+1)*pdeg + j+1
             ijpp   = i*pdeg + (j+1)+1
             ippjpp = (i+1)*pdeg + (j+1)+1
             
             temp(:,ij)     = temp(:,ij)     + inv_deg*sqrts(deg-i+1)*sqrts(deg-j+1) *d_l_1(:,i*deg+j+1)*cos_beta
             temp(:,ippj)   = temp(:,ippj)   - inv_deg*sqrts(i+2)    *sqrts(deg-j+1) *d_l_1(:,i*deg+j+1)*sin_beta
             temp(:,ijpp)   = temp(:,ijpp)   + inv_deg*sqrts(deg-i+1)*sqrts(j+2)     *d_l_1(:,i*deg+j+1)*sin_beta
             temp(:,ippjpp) = temp(:,ippjpp) + inv_deg*sqrts(i+2)    *sqrts(j+2)     *d_l_1(:,i*deg+j+1)*cos_beta
          end do
       end do
       d_l  = 0._dp
       deg = 2*l
       pdeg = deg+1
       inv_deg = (1._dp/real(deg,kind=dp))
       do i=0,deg-1
          do j=0,deg-1
             ij     = i*pdeg + j+1
             ippj   = (i+1)*pdeg + j+1
             ijpp   = i*pdeg + (j+1)+1
             ippjpp = (i+1)*pdeg + (j+1)+1
             d_l(:,ij)     = d_l(:,ij)     + inv_deg*sqrts(deg-i+1)*sqrts(deg-j+1) *temp(:,i*deg+j+1)*cos_beta
             d_l(:,ippj)   = d_l(:,ippj)   - inv_deg*sqrts(i+2)    *sqrts(deg-j+1) *temp(:,i*deg+j+1)*sin_beta
             d_l(:,ijpp)   = d_l(:,ijpp)   + inv_deg*sqrts(deg-i+1)*sqrts(j+2)     *temp(:,i*deg+j+1)*sin_beta
             d_l(:,ippjpp) = d_l(:,ippjpp) + inv_deg*sqrts(i+2)    *sqrts(j+2)     *temp(:,i*deg+j+1)*cos_beta
          end do
       end do
    end if
  end function wigner_mn_recurrence_risbo

  !> ---------
  !! @brief Computes the small Wigner d matrix $d^l\\_{m,n}(\beta)$ for fixed $l$.
  !!
  !! This function uses the recurrence by Risbo et. al. from function `wigner_mn_recurrence_risbo`.
  function wigner_dl_risbo(l,betas,normalized) result(dl)
    integer(kind = idp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical, intent(in) :: normalized
    real(kind = dp) :: dl_tmp(Size(betas,1),(2*l+1)*(2*l+1)), dl(SIZE(betas,1),2*l+1,2*l+1)
    real(kind=dp) :: cos_b(Size(betas,1)),sin_b(Size(betas,1)),ls_sqrt(2*l+1)
    integer(kind = idp) :: i,j
    ! Setup
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    do i=0_idp,2_idp*l
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl_tmp = 0._dp
    do i=0_idp,l
       dl_tmp(:,1:(2*i+1)*(2*i+1)) = wigner_mn_recurrence_risbo(dl_tmp(:,1:(2*i-1)*(2*i-1)), i, cos_b, sin_b, ls_sqrt)
    end do

    ! Save transposed array to output
    do i=1,2*l+1
       do j=1,2*l+1
          dl(:,j,i) = dl_tmp(:,(i-1_idp)*(2*l+1_idp)+j)
       end do
    end do

    if (normalized) then
       dl = dl * SQRT(real(2*l+1_idp,kind=dp)/2._dp)
    end if
  end function wigner_dl_risbo

  !>--------
  !! @brief Symmetry reduced risbo recurrence for $d^l\\_{m,n}$
  !!
  !! Risbo recurrence for symmetry reduce small wigner-d matrices.
  function wigner_recurrence_risbo_reduced(dl_1,l,cos_beta,sin_beta,sqrts,normalized) result(dl)
    real(kind = dp), intent(in) :: cos_beta(:),sin_beta(:),sqrts(:),dl_1(:,:)
    integer(kind = idp), intent(in)  :: l
    logical, intent(in) :: normalized
    real(kind = dp) :: dl(Size(cos_beta,1),((l+1)*(l+2))/2),dl_halve(Size(cos_beta,1),(l*(l+1))/2),inv_deg
    real(kind = dp) :: l_sym,lpp_sym,n_sym,norm_const,norm_const2
    integer(kind = idp) :: ll,m,n,mn,nbeta
    procedure(triangular_to_flat_index), pointer :: tri

    tri => triangular_to_flat_index

    if (l==0) then
       dl=MERGE(SQRT(0.5_dp),1._dp,normalized)
    else if (l==1) then
       norm_const=SQRT(1.5_dp)
       dl(:,1)=MERGE((cos_beta**2-sin_beta**2)*norm_const    ,cos_beta**2-sin_beta**2,normalized)
       dl(:,2)=MERGE(sqrts(3)*cos_beta*sin_beta*norm_const ,sqrts(3)*cos_beta*sin_beta,normalized)
       dl(:,3)=MERGE(cos_beta**2*norm_const                ,cos_beta**2,normalized)
    else
       !The following is going to become a bit tediouse ...

       l_sym   = real((-1_idp)**l,dp)
       lpp_sym = real((-1_idp)**(l+1_idp),dp)
       nbeta = SIZE(cos_beta,1)
       !!!!!!!!!!!!!!!!!!
       ! Step l -> l+1/2!
       !!!!!!!!!!!!!!!!!!
       inv_deg = 1._dp/real(2_idp*l-1_idp,dp)
       ll=l-1_idp
       norm_const = MERGE(inv_deg*SQRT(real(2_idp*l,dp)/real(2_idp*l-1,dp)),inv_deg,normalized)
       !! 0 <= m <= n <= l-2
       do m=0,l-2_idp
          !!! m=n
          mn = tri(m,m,l)
          !print *,m,m,mn
          dl_halve(:,mn)= ( dl_1(:,mn)                  *real(l+m,dp)                *cos_beta &
                        &  -dl_1(:,tri(m,m+1_idp,l))     *SQRT(real((l+m)*(ll-m),dp)) *sin_beta &
                        &  -dl_1(:,tri(m,m+1_idp,l))     *SQRT(real((ll-m)*(l+m),dp)) *sin_beta &
                        &  +dl_1(:,tri(m+1_idp,m+1_idp,l))*real(ll-m,dp)               *cos_beta &
                        & )*norm_const
          !print *, dl_1(1,:)
          !print *, tri(m,m+1_idp,l),tri(m+1_idp,m+1_idp,l),dl_halve(1,mn)
          
          do n=m+1_idp,l-2_idp
             mn = tri(m,n,l)
             !print *,m,n,mn
             dl_halve(:,mn)= ( dl_1(:,mn)                  *SQRT(real((l+m)*(l+n),dp))  *cos_beta &
                           &  -dl_1(:,tri(m,n+1_idp,l))     *SQRT(real((l+m)*(ll-n),dp)) *sin_beta &
                           &  +dl_1(:,tri(m+1_idp,n,l))     *SQRT(real((ll-m)*(l+n),dp)) *sin_beta &
                           &  +dl_1(:,tri(m+1_idp,n+1_idp,l))*SQRT(real((ll-m)*(ll-n),dp))*cos_beta &
                           & )*norm_const
          end do
          !!! n=l-1
          mn = tri(m,l-1_idp,l)
          !print *,m,l-1_idp,mn
          dl_halve(:,mn)=( dl_1(:,mn)                   *SQRT(real((l+m)*(l+l-1_idp),dp))  *cos_beta &
                        &  +dl_1(:,tri(m+1_idp,l-1_idp,l))*SQRT(real((ll-m)*(l+l-1_idp),dp)) *sin_beta &
                        & )*norm_const

       end do
       !! m=n=l-1       
       mn = tri(l-1_idp,l-1_idp,l)
       !print *,l-1_idp,l-1_idp,mn
       dl_halve(:,mn)= dl_1(:,mn) * real((l+l-1_idp),dp)*cos_beta * norm_const


       
       !!!!!!!!!!!!!!!!!!!!
       ! Step l+1/2 -> l+1!
       !!!!!!!!!!!!!!!!!!!!
       inv_deg = 1._dp/real((2_idp*l),dp)
       norm_const = MERGE(inv_deg*SQRT(real(2_idp*l+1,dp)/real(2_idp*l,dp)),inv_deg,normalized)
       norm_const2 = MERGE(0.5_dp*SQRT(real(2_idp*l+1,dp)/real(2_idp*l,dp)),0.5_dp,normalized)
       !! m=n=0
       mn = tri(0_idp,0_idp,l+1)
       dl(:,mn)  = (  2._dp  *dl_halve(:,tri(0_idp,0_idp,l))         *cos_beta &
                 &  - lpp_sym*dl_halve(nbeta:1:-1,tri(0_idp,0_idp,l))*sin_beta &
                 &  + l_sym  *dl_halve(nbeta:1:-1,tri(0_idp,0_idp,l))*sin_beta &
                 & )*norm_const2
       !! m=0 1<=n<=l-1
       do n=1,l-1
          mn = tri(0_idp,n,l+1_idp)
          n_sym = real((-1_idp)**(-n),dp)

          dl(:,mn)  = ( l_sym*n_sym  *dl_halve(nbeta:1:-1,tri(0_idp,n-1_idp,l))*SQRT(real(l*(l+n),dp))*cos_beta &
                    &  -lpp_sym*n_sym*dl_halve(nbeta:1:-1,tri(0_idp,n,l))     *SQRT(real(l*(l-n),dp))*sin_beta &
                    &  +              dl_halve(:,tri(0_idp,n-1_idp,l))         *SQRT(real(l*(l+n),dp))*sin_beta &
                    &  +              dl_halve(:,tri(0_idp,n,l))              *SQRT(real(l*(l-n),dp))*cos_beta &
                    & )*norm_const
       end do
       !! m=0 n=l
       mn = tri(0_idp,l,l+1_idp)
       dl(:,mn)  = ( dl_halve(nbeta:1:-1,tri(0_idp,l-1_idp,l))*cos_beta &
                 &  +dl_halve(:,tri(0_idp,l-1_idp,l))         *sin_beta &
                 & )*norm_const2*SQRT(2._dp)
       !! 1<=m<=n<=l-1
       do m=1,l-1
          !!! m=n
          mn = tri(m,m,l+1_idp)
          dl(:,mn)= ( dl_halve(:,tri(m-1_idp,m-1_idp,l))*real(l+m,dp)               *cos_beta &
                  &  -dl_halve(:,tri(m-1_idp,m,l))     *SQRT(real((l+m)*(l-m),dp)) *sin_beta &
                  &  -dl_halve(:,tri(m-1_idp,m,l))     *SQRT(real((l-m)*(l+m),dp)) *sin_beta &
                  &  +dl_halve(:,tri(m,m,l))          *real(l-m,dp)               *cos_beta &
                  & )*norm_const
          do n=m+1,l-1
             mn = tri(m,n,l+1_idp)
             dl(:,mn)= ( dl_halve(:,tri(m-1_idp,n-1_idp,l))*SQRT(real((l+m)*(l+n),dp))*cos_beta &
                     &  -dl_halve(:,tri(m-1_idp,n,l))     *SQRT(real((l+m)*(l-n),dp))*sin_beta &
                     &  +dl_halve(:,tri(m,n-1_idp,l))     *SQRT(real((l-m)*(l+n),dp))*sin_beta &
                     &  +dl_halve(:,tri(m,n,l))          *SQRT(real((l-m)*(l-n),dp))*cos_beta &
                     & )*norm_const
          end do
          !!! n=l
          mn = tri(m,l,l+1_idp)
          dl(:,mn)= ( dl_halve(:,tri(m-1_idp,l-1_idp,l))*SQRT(real((l+m)*(l+l),dp))*cos_beta &
                  &  +dl_halve(:,tri(m,l-1_idp,l))     *SQRT(real((l-m)*(l+l),dp))*sin_beta &
                  & )*norm_const
       end do
       !! m=n=l
       mn = tri(l,l,l+1_idp)
       dl(:,mn)= dl_halve(:,tri(l-1_idp,l-1_idp,l)) * real((l+l),dp)*cos_beta * norm_const
    end if
  end function wigner_recurrence_risbo_reduced

  !> ---------
  !! @brief Computes symmetry reduced part of the small Wigner d matrix $d^l\\_{m,n}(\beta)$ for fixed $l$.
  !!
  !! This function uses the recurrence by Risbo et. al. from function `wigner_mn_recurrence_risbo_reduced`.
  !! It only computes the non symmetry equivalent part that is its values for  0<=m<=n<=l.
  function wigner_dl_risbo_reduced(l,betas,normalized) result(dl)
    integer(kind = idp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical, intent(in) :: normalized
    real(kind = dp) :: dl(Size(betas,1),((l+1)*(l+2))/2)
    real(kind=dp) :: cos_b(Size(betas,1)),sin_b(Size(betas,1)),ls_sqrt(2*l+1)
    integer(kind = idp) :: i,ndl,ndlmm
    ! Setup
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    do i=0_idp,2_idp*l
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl=0
    do i=0_idp,l
       ndl   = ((i+1)*(i+2))/2
       ndlmm = (i*(i+1))/2
       dl(:,1:ndl) = wigner_recurrence_risbo_reduced(dl(:,1:ndlmm), i, cos_b, sin_b, ls_sqrt,normalized)
    end do
  end function wigner_dl_risbo_reduced

  !> --------
  !! @brief Computes all $d^l\\_{m,n}$ for $0 \\leq m \\leq n \leq l < bw$.
  !! 
  !! Computes all independent small Wigner d-matrix elements d^l\\_{m,n}(\beta), i.e. those  
  !! for $0 \\leq m \\leq n \\leq l <bw$.
  !! This routine uses risbo recurrence `wigner_recurrence_risbo_reduced`.
  !! It however uses the mnl ordering same as `genwig_all` to store the computed values 
  subroutine genwig_all_risbo_preallocated(bw,wigners,normalized)
    integer(kind = idp),intent(in) :: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp), intent(inout) :: wigners(:,:)
    real(kind = dp) :: dl(2_idp*bw,(bw*(bw+1_idp))/2)
    real(kind = dp) :: betas(2_idp*bw),cos_b(2_idp*bw),sin_b(2_idp*bw),ls_sqrt(2_idp*bw)
    integer(kind=dp) :: ndl,ndlmm,i,l,mn,mnl,m,n

    betas = create_beta_samples(2_idp*bw)
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    do i=0_idp,2_idp*bw
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl=0
    do l=0_idp,bw-1_idp
       ndl   = ((l+1)*(l+2))/2
       ndlmm = (l*(l+1))/2
       dl(:,1:ndl) = wigner_recurrence_risbo_reduced(dl(:,1:ndlmm), l, cos_b, sin_b, ls_sqrt,normalized)
       do m=0_idp,l
          do n=m,l
             mn = triangular_to_flat_index(m,n,l+1_idp)
             mnl = mnl_to_flat_index(m,n,l,bw)
             !$OMP SIMD
             do i=1,2_idp*bw
                wigners(i,mnl) = dl(i,mn)
             end do
             !$OMP END SIMD
          end do
       end do
    end do
  end subroutine genwig_all_risbo_preallocated

  !> --------
  !! @brief Computes all $d^l\\_{m,n}$ for $0 \\leq m \\leq n \leq l < bw$.
  !! 
  !! Computes all independent small Wigner d-matrix elements d^l\\_{m,n}(\beta), i.e. those  
  !! for $0 \\leq m \\leq n \\leq l <bw$.
  !! This routine uses risbo recurrence `wigner_recurrence_risbo_reduced`.
  !! It however uses the mnl ordering same as `genwig_all` to store the computed values
  function genwig_all_risbo(bw,normalized) result(wigners)
    
    integer(kind = idp),intent(in):: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp) :: wigners(bw*2_idp,(bw * (2_idp + 3_idp*bw + bw*bw))/6_idp)

    call genwig_all_risbo_preallocated(bw,wigners,normalized)
  end function genwig_all_risbo

  
  !> --------
  !! @brief Computes all $D^l\\_{m,n}(\alpha,\beta,\gamma))$ for $0 \\leq m \\leq n \leq l < bw$.
  !! 
  !! Computes all independent Wigner D-matrix elements D^l\\_{m,n}(\alpha,\beta,\gamma), i.e. those  
  !! for $0 \\leq m \\leq n \\leq l <bw$.
  !! This routine uses risbo recurrence `wigner_recurrence_risbo_reduced`.
  !! It however uses the mnl ordering same as `genwig_all` to store the computed values
  function full_big_d_risbo(bw,alpha,beta,gamma,normalized) result(wigners)
    integer(kind = idp),intent(in):: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp), intent(in) :: alpha,beta,gamma
    complex(kind=dp) :: wigners((4_idp*(bw*bw*bw)-bw)/3_idp)
    complex(kind= dp) :: aarg,garg
    real(kind = dp) :: dl(2,(bw*(bw+1_idp))/2)
    real(kind = dp) :: betas(2),cos_b(2),sin_b(2),ls_sqrt(2_idp*bw),two_pi_inv
    integer(kind=dp) :: ndl,ndlmm,i,l,mn,mnl,m,n,loc,sym_const_n,sym_const_m,sym_const_l
    
    betas(1)=beta
    betas(2)=pi-beta
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    two_pi_inv = 1._dp/(2._dp*pi)
    do i=0_idp,2_idp*bw
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl=0
    do l=0_idp,bw-1_idp
       ndl   = ((l+1)*(l+2))/2
       ndlmm = (l*(l+1))/2
       sym_const_l = (-1._dp)**l
       dl(:,1:ndl) = wigner_recurrence_risbo_reduced(dl(:,1:ndlmm), l, cos_b, sin_b, ls_sqrt,normalized)
       do m=0_idp,l
          sym_const_m = (-1._dp)**m
          do n=m,l
             sym_const_n = (-1._dp)**n
             mn = triangular_to_flat_index(m,n,l+1_idp)

             !populate coefficients
             loc = coeff_location_mnl(m, n, l, bw)
             aarg = EXP(cmplx(0._dp,real(-m,kind = dp)*alpha,kind = dp))
             garg = EXP(cmplx(0._dp,real(-n,kind = dp)*gamma,kind = dp))
             wigners(loc) = aarg * garg * dl(1,mn)
             
             if ((m==0 .AND. n==0)) cycle
             loc = coeff_location_mnl(-n, -m, l, bw)
             aarg = EXP(cmplx(0._dp,real(n,kind = dp)*alpha,kind = dp))
             garg = EXP(cmplx(0._dp,real(m,kind = dp)*gamma,kind = dp))
             wigners(loc) = aarg * garg * dl(1,mn)
             
             if (.NOT. m==n) then
                loc = coeff_location_mnl(n, m, l, bw)
                aarg = EXP(cmplx(0._dp,real(-n,kind = dp)*alpha,kind = dp))
                garg = EXP(cmplx(0._dp,real(-m,kind = dp)*gamma,kind = dp))
                wigners(loc) = (sym_const_m*sym_const_n)*aarg * garg * dl(1,mn)
                
                loc = coeff_location_mnl(-m, -n, l, bw)
                aarg = EXP(cmplx(0._dp,real(m,kind = dp)*alpha,kind = dp))
                garg = EXP(cmplx(0._dp,real(n,kind = dp)*gamma,kind = dp))
                wigners(loc) = (sym_const_m*sym_const_n)*aarg * garg * dl(1,mn)
             end if
             
             if (m==0 .or. n==0) cycle
             loc = coeff_location_mnl(-m, n, l, bw)
             aarg = EXP(cmplx(0._dp,real(m,kind = dp)*alpha,kind = dp))
             garg = EXP(cmplx(0._dp,real(-n,kind = dp)*gamma,kind = dp))
             wigners(loc) = (sym_const_m*sym_const_l)*aarg * garg * dl(2,mn)
             
             loc = coeff_location_mnl(m, -n, l, bw)
             aarg = EXP(cmplx(0._dp,real(-m,kind = dp)*alpha,kind = dp))
             garg = EXP(cmplx(0._dp,real(n,kind = dp)*gamma,kind = dp))
             wigners(loc) = (sym_const_n*sym_const_l) * aarg * garg * dl(2,mn)
             
             if (m==n) cycle
             loc = coeff_location_mnl(-n, m, l, bw)
             aarg = EXP(cmplx(0._dp,real(n,kind = dp)*alpha,kind = dp))
             garg = EXP(cmplx(0._dp,real(-m,kind = dp)*gamma,kind = dp))
             wigners(loc) = (sym_const_m*sym_const_l) * aarg * garg * dl(2,mn)
             
             loc = coeff_location_mnl(n, -m, l, bw)
             aarg = EXP(cmplx(0._dp,real(-n,kind = dp)*alpha,kind = dp))
             garg = EXP(cmplx(0._dp,real(m,kind = dp)*gamma,kind = dp))
             wigners(loc) = (sym_const_n*sym_const_l) * aarg * garg * dl(2,mn)
             
          end do
       end do
    end do
    
    if (normalized) then
       wigners = wigners * two_pi_inv
    end if
  end function full_big_d_risbo
  
  !> ------
  !! @brief Convert symmetry reduced Wigner-d to full Wigner-d
  !!
  !! Converts symmetry reduced $d^l\\_{m,n}$, i.e. with $0\\leq m \\leq n \leq l$, to
  !! the full matrix with $-l \leq m,n \leq l$.
  function sym_reduced_to_full_wigner(sym_dl,l) result(dl)
        real(kind = dp), intent(in) :: sym_dl(:,:)
        integer(kind = idp), intent(in) :: l
        real(kind = dp) :: dl(SIZE(sym_dl,1),2*l+1,2*l+1)
        !real(kind = dp), intent(out) :: dl(SIZE(sym_dl,1),int(SQRT(9_dp/4_idp-2_idp*real(SIZE(sym_dl,2),dp))-3_idp/2_idp) ,int(SQRT(9_dp/4_idp - 2_idp*real(SIZE(sym_dl,2),dp))-3_idp/2_idp)) ! would get rid of the degree l as input but f2py cant handle this ...
        integer(kind = idp) :: m,n,mn,nbeta,mid,nid,neg_mid,neg_nid
        real(kind=dp) :: sym_const_m,sym_const_n,sym_const_l
        
        sym_const_l=(-1._dp)**l
        nbeta=SIZE(sym_dl,1)
        dl=0
        do m=0,l
           mid = l+1_idp+m
           neg_mid = l+1_idp-m
           sym_const_m = (-1._dp)**m

           do n=m,l
              nid = l+1_idp+n
              neg_nid = l+1_idp-n
              sym_const_n = (-1._dp)**n

              mn = triangular_to_flat_index(m,n,l+1_idp)
              
              ! populate the small wigner d-matrix using all of the available symmetries
              dl(:,nid,mid)=sym_dl(:,mn)          
              if ((m==0 .AND. n==0)) cycle
              dl(:,neg_mid,neg_nid) = sym_dl(:,mn)
              if (.NOT. m==n) then
                 dl(:,mid,nid) = (sym_const_m*sym_const_n)*sym_dl(:,mn)
                 dl(:,neg_nid,neg_mid) = (sym_const_m*sym_const_n)*sym_dl(:,mn)
              end if
              if (m==0 .or. n==0) cycle
              dl(:,neg_nid,mid)=sym_dl(nbeta:1:-1,mn)*(sym_const_l*sym_const_m)
              dl(:,nid,neg_mid)=sym_dl(nbeta:1:-1,mn)*(sym_const_l*sym_const_n)
              if (m==n) cycle
              dl(:,neg_mid,nid)=sym_dl(nbeta:1:-1,mn)*(sym_const_l*sym_const_m)
              dl(:,mid,neg_nid)=sym_dl(nbeta:1:-1,mn)*(sym_const_l*sym_const_n)
           end do
        end do
      end function sym_reduced_to_full_wigner
end module make_wigner
