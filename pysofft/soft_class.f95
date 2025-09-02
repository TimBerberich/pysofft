! Define constants
module precision
implicit none
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)
end module precision
module math_constants
  use precision
  implicit none
  real(kind=dp), parameter :: pi = 4.0_dp*atan2(1.0_dp,1.0_dp)
end module math_constants
module utils
  use precision
  use math_constants
  implicit none
contains 
  function num_coeffs(m1,m2,bw) result(num_coeff)
    !  For orders m1, m2, and bandwidth bw, returns how many coefficients
    !  will be computed.
    !  
    !  let m = Max( |m1|, |m2| )
    !  
    !  The number of coefficients is = bw - m
    integer(kind=dp) :: m1,m2,bw,num_coeff
    if (abs(m1)>=bw .OR. abs(m2)>=bw) then
       print *, "Out of bounds error (n_coeffs): |m1|>bw or |m2|>bw"
       stop
    end if
    num_coeff = bw - max(abs(m1),abs(m2))
  end function num_coeffs
  function total_num_coeffs(bw) result(num_coeff)
    !  For bandwidth bw, returns the TOTAL NUMBER of coefficients
    !  that will be computed when doing a full SO3 forward transform.
    !
    !  total number = 1/3 * bw * (4 * bw^2 - 1 ) 
    !
    !  note that for integer bw, 4 * bw^3 - bw ,is always a multiple of 3,
    !  so integer division won't be messy.   
    integer(kind=8) bw,num_coeff
    num_coeff = (4_dp*(bw*bw*bw)-bw)/3_dp
  end function total_num_coeffs

  function order_to_ids(m1,m2,bw) result(ids)
    ! For order m1, m2 and bandwidth bw, returns the slice in the 3-D,
    ! fft'd array of the data necessary for the Wigner-(m1,m2)
    ! transform. I.e. The slice of the fft'd data
    ! needed to compute f_{m1,m2}^l for all legal l. 
    !
    ! Note that the array is of size 2*bw x 2*bw x 2*bw, so I
    ! need to multiply by that 2*bw (I'm sampling at twice
    ! the bandwidth)
    integer(kind=dp) :: m1,m2,bw,bw2
    integer(kind=dp) :: ids(2)
    if (abs(m1)>=bw .OR. abs(m2)>=bw) then
       print *, "Out of bounds error (sample_slice): |m1|>bw or |m2|>bw"
       stop
    end if
    bw2=2*bw
    ids(1) = (1_dp-sign(1_dp,m1))/2_dp*bw2 + m1+1_dp
    ids(2) = (1_dp-sign(1_dp,m2))/2_dp*bw2 + m2+1_dp
  end function order_to_ids
  function coeff_slice_legacy(m1,m2,bw) result(slice)
    ! So sample_slice tells me the location of the samples 
    ! needed to do the order (m1,m2)-Wigner-transform on the
    ! just fft'd data.
    !
    ! This function, coeff_slice, tells me where in the coefficient
    ! array I should place the evaluated inner-products
    ! for bandwidth bw of the transform.
    integer(kind = dp) :: m1,m2,bw,sqrbw,m,size,slice(2)
    integer(kind = dp) :: start,a1,a2,lsum
    a1 = abs(m1)
    a2 = abs(m2)
    m = max(a1,a2)
    if (m>=bw) then
       print *, "Out of bounds error (coeff_slice): |m1|>bw or |m2|>bw"
       stop
    end if
    sqrbw = bw*bw
    size = bw-m       
    
    ! lsum is the number of coefficients with m1,n2,bw 
    ! where n2 ranges from (0 to m2-1) or (m2 to 1) depending on the sign of m2
    ! NOTE the below code coputes this without branching using
    ! (1_dp-sign(1_dp,a2-a1))/2_dp  <=>  if a1>=a2 then 1_dp else 0_dp
    ! (1_dp-sign(1_dp,a1-a2-1))/2_dp  <=>  if a2>a1 then 1_dp else 0_dp
    lsum = (1_dp-sign(1_dp,a2-a1))/2_dp * a2*(bw-a1) &     
         &+ (1_dp-sign(1_dp,a1-a2-1_dp))/2_dp &
         &* (a2*(2_dp*bw-a2+sign(1_dp,m2))-a1*(a1+sign(1_dp,m2)))/2_dp
    lsum = lsum*sign(1_dp,m2)
    
    
    ! no branching version
    ! (1_dp-sign(1_dp,m1))/2_dp  <=>  if m1>=0 then 1_dp else 0_dp
    start = total_num_coeffs(bw)*(1_dp-sign(1_dp,m1))/2_dp &
         &+ sqrbw*(m1+(1_dp-sign(1_dp,m2))/2_dp) &
         &- (m1-sign(1_dp,m2))*m1*(2_dp*m1 - sign(1_dp,m2))/6_dp
    
    start = start + lsum
    start = start + 1_dp ! +1 due to fortran indexing starting at 1   
    slice(1) = start
    slice(2) = start + (size-1)
  end function coeff_slice_legacy
  function coeff_location(m1,m2,l,bw) result(id)
    ! This function returns the index of the coefficient array corresponding to f_{m1,m2}^l.
    ! Note that this index has to lie within the slice given by coefLoc_so3.
    integer(kind = dp) :: m1,m2,bw,l,m,slice(2),id
    slice = coeff_slice(m1,m2,bw)
    m =  max(abs(m1),abs(m2))
    if (l<m .OR. l>bw) then
       print *, "Out of bounds error (coeff_location): l does not satisfy max(|m1|,|m2|) <= l < bw "
       stop
    end if
    id = slice(1)+l-m
  end function coeff_location
  function coeff_slice(n1,n2,bw) result(slice)
    ! Returns the slice of coefficients f_(n1,n2)^l for all possible l valaues.
    ! It is chosen such that in loops of the follwing type memory access
    ! for the coefficients is continuouse without jumps.
    !
    ! do m1=0,bw-1
    !    do m2=m1,bw-1
    !      f_{m1,m2}
    !      if m1==0 and m2==0, cycle  
    !      f_{-m2,-m1}
    !      if m1/=m2
    !         f_{m2,m1}
    !         f_{-m1,-m2}
    !      if m1==0 or m2==0
    !      f_{m1,-m2}
    !      f_{-m1,m2}
    !      if m1==m2, cycle
    !      f_{m2,-m1}
    !      f_{-m2,m1}
    !
    integer(kind=dp), intent(in) :: n1,n2,bw
    integer(kind=dp) :: m(2),slice(2),an1,an2,temp1,temp2,temp3
    logical :: swapped,n1_neg,n2_neg
    an1 = ABS(n1)
    an2 = ABS(n2)
    m = MERGE([an1,an2],[an2,an1],an1<an2)

    temp1 = 4*(bw*m(1)*(bw-m(1)+1))-2*(bw**2+m(1)**2)-bw + m(1)*(4*m(1)**2+2)/3
    slice(1) = MERGE(temp1,0_dp,m(1)>0)
    
    temp1 = bw+2*(2*bw-m(2))*(m(2)-1)
    temp2 = 4*((bw-m(1))+(m(2)-1 - m(1))*(2*bw - m(1) - m(2)))
    slice(1) = slice(1)+MERGE(MERGE(temp1,temp2,m(1)==0),0_dp,m(1)/=m(2))

    swapped = (m(1)/=ABS(n1))
    n1_neg = (n1<0)
    n2_neg = (n2<0)
  
    temp1 = MERGE(1_dp,MERGE(2_dp,3_dp,(n2_neg .AND.  n1/=0_dp ) .OR. n2==0_dp )&
         &, n1_neg .AND. (n2_neg .OR. n2==0_dp  ) )
    
    temp2 = MERGE(MERGE(MERGE(1_dp,3_dp,swapped),MERGE(7_dp,5_dp,swapped),n2_neg) &
         &,MERGE(MERGE(6_dp,4_dp,swapped),2_dp,n2_neg),n1_neg)
  
    temp3 = MERGE(temp1,temp2, (m(1)==0) .OR. (m(2)==0) .OR. (m(1)==m(2)) )
    slice(1) = slice(1) + MERGE(temp3*(bw-m(2)),0_dp,(m(1)/=n1) .OR. (m(2)/=n2))
    slice(2) = slice(1) + bw-m(2)
    slice(1) = slice(1) + 1_dp ! 1 indexing
  end function coeff_slice
  
  function wigLen_so3(m1,m2,bw) result(wigLen)
    ! Returns the number of small wigner d values d_{m1,m2}^l(\beta)
    ! for a given choice of m1,m2 and bw.
    ! Since beta is sampled at twice the bandwith the number of values is
    !
    ! number of values = (number of possible l) * 2*bw
    integer(kind = dp) :: m1,m2,m,bw,wigLen
    m = max(abs(m1),abs(m2))
    if (m>=bw) then
       print *, "Out of bounds error (wigLen_so3): |m1|>bw or |m2|>bw"
       stop
    end if
    wigLen = (bw-m)*2*bw 
  end function wigLen_so3
  function euler_shape(bw) result(eshape)
    integer(kind=dp), intent(in) :: bw
    integer(kind=dp) :: eshape(3)
    eshape = 2*bw
  end function euler_shape
  
  function legendre_quadrature_weights(bw) result(weights)
    ! This function computes the even legendre weights
    ! used in the quadrature of the SOFT algorithm
    ! For a fixed bandwith (bw) the weights are the
    ! unique solutions to the problem:
    ! $$\sum_{k=0}^{2\text{bw}-1} w_{\text{bw}}(k) P_m(\cos(\beta_k)) = \delta_{0,m} \quad \forall 0\leq m \leq \text{bw}$$
    ! where the sampling angles are given by $$\beta_k = \frac{\pi(2k+1)}{4\text{bw}} $$
    ! Their closed form expression is given by:
    ! $$ w_{\text{bw}}(k) = \frac{2}{\text{bw}}\sin(\frac{\pi(2k+1)}{4\text{bw}}) \sum_{j=0}^{\text{bw}-1}\frac{1}{2j+1}\sin((2k+1)(2j+1)\frac{\pi}{4\text{bw}})$$
    ! both of these Formulas can be found in equations (2.13) and (2.14) of
    ! :cite: P.J. Kostelec and D.N. Rockmore, J Fourier Anal Appl (2008) 14: 145–179
    ! The proof of the closed form is contained in>
    ! :cite: Driscoll, J.R., Healy, D.:  Proc. 34th IEEE FOCS (1989), pp. 344–349. Adv. in Appl. Math., vol. 15, pp. 202–250 (1994)
    integer(kind = dp) :: bw,j,k
    real(kind = dp) :: weights(2*bw),tempsum,k_odd,j_odd,xi
    xi = pi/real(4_dp*bw,kind=dp)
    
    do k=0,2*bw-1
       k_odd = real(2_dp*k+1_dp,kind=dp)
       tempsum = 0
       do j=0,bw-1
          j_odd = real(2_dp*j+1_dp,kind=dp)
          tempsum = tempsum + (1.0_dp/j_odd)*SIN(k_odd*j_odd*xi)
       end do
       tempsum = tempsum*( (2.0_dp/real(bw,kind=dp)) * SIN(k_odd*xi) )
       weights(k+1_dp)=tempsum
    end do
  end function legendre_quadrature_weights

  function get_empty_coeff(bw) result(coeff)
    integer(kind = dp) :: bw
    complex(kind = dp) :: coeff((4_dp*(bw*bw*bw)-bw)/3_dp)
    coeff = 0.0
  end function get_empty_coeff
  function get_empty_so3func_cmplx(bw) result(so3func)
    integer(kind = dp) :: bw
    complex(kind = dp) :: so3func(2*bw,2*bw,2*bw)
    so3func = 0.0
  end function get_empty_so3func_cmplx
  function get_empty_so3func_halfcmplx(bw) result(so3func)
    integer(kind = dp) :: bw
    complex(kind = dp) :: so3func(2*bw,bw+1,2*bw)
    so3func = 0.0
  end function get_empty_so3func_halfcmplx
  function get_empty_so3func_real(bw) result(so3func)
    integer(kind = dp) ::  bw
    real(kind = dp) :: so3func(2*bw,2*bw,2*bw)
    so3func = 0.0
  end function get_empty_so3func_real
  subroutine enforce_real_sym(coeff,bw)
    complex(kind=dp) ,intent(inout) :: coeff(:)
    integer(kind=dp) ,intent(in) :: bw
    real(kind=dp) :: sym_const_m1,sym_const_m2
    integer(kind=dp) :: m1,m2,c_slice(2),c_slice_sym(2)

    c_slice = coeff_slice(0_dp,0_dp,bw)
    coeff(c_slice(1):c_slice(2)) = coeff(c_slice(1):c_slice(2))%re
    do m1=-bw+1,bw-1
       sym_const_m1 = (-1.0)**m1
       do m2=-bw+1,bw-1
          sym_const_m2 = (-1.0)**m2
          c_slice = coeff_slice(m1,m2,bw)
          c_slice_sym = coeff_slice(-m1,-m2,bw)
          coeff(c_slice(1):c_slice(2)) = sym_const_m1*sym_const_m2*CONJG(coeff(c_slice_sym(1):c_slice_sym(2)))       
       end do
    end do
  end subroutine enforce_real_sym

  pure function triangular_size(bw) result(tri_size)
    !! Consider the triangular index 0<=i<=j<bw
    !! This function returns the total number of possible pairs i,j
    integer(kind = dp),intent(in) :: bw
    integer(kind = dp) :: tri_size
    tri_size = (bw*(bw+1))/2_dp
  end function triangular_size
  
  pure function pyramid_size(bw) result(pyr_size)
    !! consider the pyramid index 0<=i<bw and -i<=j<=i
    !!This function returns the total number of possible pairs i,j
    integer(kind = dp),intent(in) ::bw
    integer(kind = dp) :: pyr_size
    pyr_size = bw**2
  end function pyramid_size
  
  subroutine flat_to_pyramid_index(i,j,k)
    ! Converts a running index k=0 to bw*(bw+1)/2-1 into
    ! a triangular double index i,j with  0<=i<kw and i<=j<=bw
    ! This allows to reformulate Triangular loops as simple loops via
    ! do k=1, (N+1)**2
    !    call flat_to_pyramid_index(i,j,k)
    ! is the same as
    ! do i=0,N
    !   do j=-i,i
    integer(kind = dp), intent(inout) :: i,j
    integer(kind = dp), intent(in) :: k
    i = int(SQRT(real(k-1,dp)),kind = dp)
    j = -i+(k-1)-i**2 
  end subroutine flat_to_pyramid_index
  subroutine flat_to_triangular_index(i,j,k,N)
    ! Converts a running index k=0 to bw*(bw+1)/2-1 into
    ! a triangular double index i,j with  0<=i<kw and i<=j<=bw
    ! This allows to reformulate Triangular loops as simple loops via
    ! do k=1, (N+1)*(N+2)/2
    !    call flat_to_triangular_index(i,j,k,bw)
    ! is the same as
    ! do i=0,N
    !   do j=i,N
    integer(kind = dp), intent(inout) :: i,j
    integer(kind = dp), intent(in) :: k,N
    i = int(-SQRT(real((2_dp*N+3_dp)**2-8_dp*(k-1_dp),kind = dp))+real(2_dp*N+3_dp,kind=dp),kind=dp)/2_dp
    j =  (k-1_dp)-i*(2_dp*N+1_dp-i)/2_dp
    !i = (-int(SQRT(real((2_dp*N+3_dp)**2-8_dp*(k-1_dp),kind = dp)),kind = dp)+2_dp*N+3_dp)/2_dp
    !j =  (k-1_dp)-i*(2_dp*N+1_dp-i)/2_dp
  end subroutine flat_to_triangular_index

  function triangular_to_flat_index(i,j,bw) result(id)
    !! Consider the triangular index 0<=i<=j<bw 
    !! This function returns the slice of of a flattened array that corresponds
    !! to all valid j for a fixed i.
    !! This function allows to store a triangular array in a j contiguous way
    !! Designed for the following loops
    !! do i=0,bw-1
    !!    do j=i,bw-1
    integer(kind = dp), intent(in) :: i,j,bw
    integer(kind = dp) :: id
    id = (i*(2_dp*bw-i+1_dp))/2_dp + j-i + 1_dp
  end function triangular_to_flat_index
  function triangular_to_flat_index_reversed(j,i) result(id)
    !! Consider the triangular index 0<=i<=j<bw 
    !! This function returns the slice of of a flattened array that corresponds
    !! to all valid i for a fixed j.
    !! This function allows to store a triangular array in a i contiguous way
    integer(kind = dp), intent(in) :: j,i
    integer(kind = dp) :: id
    id = (j*(j+1_dp))/2_dp + i +1_dp
  end function triangular_to_flat_index_reversed

  function lmn_to_flat_index(l,m,n) result(id)
    !! Considers all possible indexes 0<=m<=n<=l
    !! such that the following loop is cintiguous in memory
    !! do i=0,l
    !!    do j=0,m
    !!       do k=j,n
    integer(kind = dp),intent(in) :: l,m,n
    integer(kind = dp) :: id
    if (.NOT. (0<=m .AND. m<=n .AND. n<=l)) then
       print *, "Invalid arguments: l,m,n have to satisfy 0<=m<=n<=l"
    end if
    id = (l*(l+1)*(l+2))/6 + (m*(2*l+3-m))/2 + n-m
  end function lmn_to_flat_index
  function mnl_to_flat_index(m,n,l,bw) result(id)
    !! Considers all possible indexes 0<=m<=n<=l<=bw
    !! such that the following loop is cintiguous in memory
    !! do j=0,m
    !!    do k=j,n
    !!       do i=k,l
    integer(kind = dp),intent(in) :: l,m,n,bw
    integer(kind = dp) :: id
    if (.NOT. (0<=m .AND. m<=n .AND. n<=l .AND. l<bw )) then
       print *, "Invalid arguments: l,m,n have to satisfy 0<=m<=n<=l<bw"
    end if
    id = (m*(3*bw**2+m**2+2-3*bw*(m-2)-3*m))/6+((n-m)*(2*bw+1-m-n))/2 + l-n+1
  end function mnl_to_flat_index
  function mnl_to_flat_l_slice(m,n,bw) result(slice)
    integer(kind = dp),intent(in) :: m,n,bw
    integer(kind = dp) :: slice(2)
    slice(1) = mnl_to_flat_index(m,n,n,bw)
    slice(2) = slice(1) + bw - 1_dp - n 
  end function mnl_to_flat_l_slice
  function mnl_to_flat_l_slice_padded(m,n,bw,pad_size) result(slice)
    integer(kind = dp),intent(in) :: m,n,bw,pad_size
    integer(kind = dp) :: slice(2)
    slice(1) = (mnl_to_flat_index(m,n,n,bw)-1_dp)*pad_size+1_dp
    slice(2) = slice(1) + pad_size*(bw - n) - 1_dp 
  end function mnl_to_flat_l_slice_padded
  
  pure function LMc(l,m) result(index)
    ! Index for the spherical harmonic coefficient Y_lm for complex data.
    ! This is the same convention used by the software SHTNS
    ! https://nschaeff.bitbucket.io/shtns/index.html
    ! Complex data => Coefficients are ordered by l and m contiguous
    ! Loops of the following kind are contiguous in memory with this indexing
    ! do l in [0,1,..., bw-1]
    !   do m in [-(bw-1),...,bw-1]
    integer(kind = dp),intent(in) :: l,m
    integer(kind = dp) :: index
    index = l*(l+1_dp)+m+1_dp
  end function LMc

  pure function MLr(m,l,bw) result(index)
    ! Index for the spherical harmonic coefficient Y_lm for real data.
    ! This is the same convention used by the software SHTNS
    ! https://nschaeff.bitbucket.io/shtns/index.html
    ! Real data => Coefficients are ordered by m and l contiguous
    ! Loops of the following kind are contiguous in memory with this indexing
    ! do m in [0,1,..., bw-1]
    !   do l in [m,...,bw-1]
    integer(kind = dp),intent(in) :: l,m,bw
    integer(kind = dp) :: index
    index = (m*(2_dp*bw-1_dp-m))/2_dp + l + 1_dp
  end function MLr
  pure function MLr_slice(m,bw) result(slice)
    integer(kind = dp), intent(in) :: m,bw
    integer(kind = dp) :: slice(2)

    slice(1) = MLr(m,0_dp,bw)
    slice(2) = slice(1) + bw-m-1_dp
  end function MLr_slice

  pure function MLc(m,l,bw) result(index)
    ! Index for the spherical harmonic coefficient Y_lm for complex data.
    ! Custom vesion that is ordered by m and l contiguous,
    ! while positive and negative m are interleaved.
    ! Loops of the following kind are contiguous in memory with this indexing
    ! do m in [0,1,-1,2,-2,..., bw-1,-(bw-1)]
    !   do l in [|m|,...,bw-1]
    integer(kind = dp),intent(in) :: l,m,bw
    integer(kind = dp) :: index
    index = MAX(0, ABS(m)*(2_dp*bw - ABS(m)) - bw) + MERGE(0_dp,bw-abs(m),m>=0) + l + 1_dp
  end function MLc
  pure function MLc_slice(m,bw) result(slice)
    ! Index for the spherical harmonic coefficient Y_lm for complex data.
    ! Custom vesion that is ordered by m and l contiguous,
    ! while positive and negative m are interleaved.
    ! Loops of the following kind are contiguous in memory with this indexing
    ! do m in [0,1,-1,2,-2,..., bw-1,-(bw-1)]
    !   do l in [|m|,...,bw-1]
    integer(kind = dp),intent(in) :: m,bw
    integer(kind = dp) :: slice(2)
    
    slice(1) = MLc(m,0_dp,bw)+abs(m)
    slice(2) = slice(1) + bw-abs(m)-1_dp
  end function MLc_slice

end module utils

! Start of routine code
module make_wigner
  use precision
  use math_constants
  use utils
  implicit none
contains
  function size_wigner_d(bw) result(wsize)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: wsize
    wsize = bw*bw * int(2_dp + 3_dp*bw + bw*bw,dp)/3_dp
  end function size_wigner_d
  function wigner_d_shape(bw) result(d_shape)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: d_shape(2)
    d_shape = [2_dp*bw,(bw*(2_dp + 3_dp*bw + bw*bw))/6_dp]
  end function wigner_d_shape
  
  function create_beta_samples(n) result(betas)
    ! returns n uniformly sampled angles in (0,pi)
    ! Note: These are angles used in Chebyshev nodes of the first kind.
    integer(kind = dp) :: n,i
    real(kind = dp) :: betas(n),factor
    factor = pi / real(2_dp * n, dp)
    betas = [(real(2_dp * i - 1_dp, dp) * factor, i = 1, n)]
  end function create_beta_samples
  function create_trig_samples(bw) result(trig_samples)
    integer(kind=dp) :: bw
    real(kind = dp) :: trig_samples(2*bw,3),betas(2*bw)
    betas = create_beta_samples(2*bw)
    trig_samples(:,1) = cos(betas)                      ! <= Chebyshev nodes
    trig_samples(:,2) = cos(betas/2_dp)*sin(betas/2_dp) ! <= needed for d^l_ml computation
    trig_samples(:,3) = cos(betas/2_dp)**2              ! <= needed for d^l_ml computation
  end function create_trig_samples

  function compute_dlml(l,m,sincos,cos2,normalized) result(dlml)
    !! Computes the Wigner little d_lmn(beta) for n=l
    !! using their exact expression
    !! $$d^l_{m,l}(\beta) = \Sqrt{\frac{2l+1}{2}} \sqrt{\frac{2l!}{(l+m)!(l-m)!}} \cos(\frac{\beta}{2})^{l+m} \sin(\frac{\beta}{2})^{l-m} $$
    !!
    !! as well as the following recursion relationships
    !!
    !!  $$d^{l+1}_{m,l+1}(\beta) = d^l_{m,l}*\sqrt{\frac{(2l+2)(2l+3)}{(l+m+1)(l-m+1)}} \cos(\frac{\beta}{2}) \sin(\frac{\beta}{2})$$
    !!  $$d^{l+1}_{m,l+1}(\beta) = d^l_{m,l}*\sqrt{\frac{(2l+2)(2l+3)}{(l+ m+1)(l+m+2)}} \cos(\frac{\beta}{2})^2 $$
    !!
    !! Note the normalization factor $\Sqrt{\frac{2l+1}{2}}$ in the above equations, for unormalized dlml the equations wold change slightly
    !!
    !! Using the recurrence relations remains stable to higher l than using the direct relation.
    !! This is because while d^l_ml(beta) remains relatively bounded, the values for the SQRT factor explode at large l which is compensated
    !! by the strongly decreasing sin,cos part. In contrast the multiplication coefficients listed above are bounded between (0,4) for all l,m,beta.
    !!
    !! It seems that iterating between the two recurrences is the most stable computation approach. No overflow/underflow occurs at least till bw=5000.

    !! Indices $l,m$ at which to compute d^l_{m,l}
    integer(kind = dp),intent(in) :: l,m
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(size(cos2,1))
    integer(kind = dp) :: i,j,ll,mm,nc1,nc2

    dlml = MERGE(SQRT(0.5_dp),1._dp,normalized) ! In agreement with our normalization $\sqrt{\frac{2l+1}{2}}$ at l=0
    !Normalization defining constants
    nc1 = MERGE(2_dp,1_dp,normalized)
    nc2 = nc1+1_dp
    
    ll=0_dp
    mm=0_dp
    do i=0,l-1       
       if ((ll<l) .AND. (l-m+mm-ll>0)) then
          ! l -> l+1 at m=0 starting from l = 0
          dlml = dlml * SQRT(real((2_dp*ll+nc1)*(2_dp*ll+nc2),kind=dp)/real((ll+mm+1_dp)*(ll-mm+1_dp),kind=dp))*sincos
          ll = ll + 1_dp
       end if
       if ((mm<m) .AND. (ll<l)) then
          ! l,m -> l+1 m+1
          dlml = dlml * SQRT(real((2_dp*ll+nc1)*(2_dp*ll+nc2),kind=dp)/real((ll+mm+1_dp)*(ll+mm+2_dp),kind=dp))*cos2
          ll = ll + 1_dp
          mm = mm + 1_dp 
       end if
       if ((ll==l) .AND. (mm==m)) exit
    end do
  end function compute_dlml
  subroutine dlml_recursion_l_contiguous(dlml,m,bw,sincos,cos2,normalized)
    !! dlml for all l at fixed m. Has to be of size [Size(sincos,1),bw]
    real(kind = dp),intent(inout) :: dlml(:,:)
    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: m,bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    integer(kind = dp) :: i,l,swap_id,nc1,nc2

    !Normalization defining constants
    nc1 = MERGE(2_dp,1_dp,normalized)
    nc2 = nc1+1_dp
    
    swap_id = m+1
    do i=1,min(bw-m,swap_id)
       !(l,m-1 -> l+1,m)
       l=m+(i-1_dp)-1_dp
       ! Don't substitute the following by dlml(:,i)*1._dp for m==0 
       ! This will cause tiny rounding errors that can accumulate for very high iterative recursion calls.
       dlml(:,i) = MERGE(dlml(:,i), dlml(:,i) * SQRT(real((2_dp*l+nc1)*(2_dp*l+nc2),kind=dp)/real((l+m)*(l+m+1_dp),kind=dp))*cos2,m==0)
    end do
    do i=swap_id,bw-m-1
       !(l-> l+1)
       l=m+(i-1_dp)
       dlml(:,i+1) = dlml(:,i) * SQRT(real((2_dp*l+nc1)*(2_dp*l+nc2),kind=dp)/real((l+m+1_dp)*(l-m+1_dp),kind=dp))*sincos
    end do
  end subroutine dlml_recursion_l_contiguous
  function compute_all_dlml_l_contiguous(bw,sincos,cos2,normalized) result(dlml)
    !! Computes all Wigner little d coefficients of the form d_lml(beta) with l<bw and stores the values in an l contiguous way.
    !! That is dlml is stored at triangular_to_flat_index(m,l,bw)
    !! 
    !! $$d^l_{m,l}(\beta) = \Sqrt{\frac{2l+1}{2}} \sqrt{\frac{2l!}{(l+m)!(l-m)!}} \cos(\frac{\beta}{2})^{l+m} \sin(\frac{\beta}{2})^{l-m} $$
    !!
    !! as well as the following recursion relationships
    !!
    !!  $$d^{l+1}_{m,l+1}(\beta) = d^l_{m,l}*\sqrt{\frac{(2l+2)(2l+3)}{(l+m+1)(l-m+1)}} \cos(\frac{\beta}{2}) \sin(\frac{\beta}{2})$$
    !!  $$d^{l+1}_{m,l+1}(\beta) = d^l_{m,l}*\sqrt{\frac{(2l+2)(2l+3)}{(l+m+1)(l+m+2)}} \cos(\frac{\beta}{2})^2 $$
    !!
    !! Note the normalization factor $\Sqrt{\frac{2l+1}{2}}$ in the above equations, for unormalized dlml the equations wold change slightly
    !!
    !! It seems that iterating between the two recurrences is the most stable computation approach. No overflow/underflow occurs at least till bw=5000.

    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(SIZE(cos2,1),(bw*(bw+1))/2_dp),dlml_tmp(Size(cos2,1),bw)
    integer(kind=dp) :: l,m,swap_id,lm_id,next_lm_id
    
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
  subroutine dlml_recursion_m_contiguous(dlml,l,bw,sincos,cos2,normalized)
    !! dlml for all l at fixed m. Has to be of size [Size(sincos,1),bw]
    real(kind = dp),intent(inout) :: dlml(:,:)
    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: l,bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    integer(kind = dp) :: i,m,swap_id,nc1,nc2

    !Normalization defining constants
    nc1 = MERGE(2_dp,1_dp,normalized)
    nc2 = nc1+1_dp
    
    swap_id = l/2
    do m=l,swap_id,-1
       !(l,m -> l+1,m+1 )
       dlml(:,m+2) = dlml(:,m+1) * SQRT(real((2_dp*l+nc1)*(2_dp*l+nc2),kind=dp)/real((l+m+1_dp)*(l+m+2_dp),kind=dp))*cos2
    end do
    do m=0,swap_id
       !(l-> l+1)
       dlml(:,m+1) = dlml(:,m+1) * SQRT(real((2_dp*l+nc1)*(2_dp*l+nc2),kind=dp)/real((l+m+1_dp)*(l-m+1_dp),kind=dp))*sincos
    end do
    
  end subroutine dlml_recursion_m_contiguous
  function compute_all_dlml_m_contiguous(bw,sincos,cos2,normalized) result(dlml)
    !! Computes all Wigner little d coefficients of the form d_lml(beta) with l<bw and stores the values in an m contiguous way.
    !! That is dlml is stored at triangular_to_flat_index_reversed(l,m,bw)
    !! 
    !! Computation is done via the exact expression
    !! $$d^l_{m,l}(\beta) = \Sqrt{\frac{2l+1}{2}} \sqrt{\frac{2l!}{(l+m)!(l-m)!}} \cos(\frac{\beta}{2})^{l+m} \sin(\frac{\beta}{2})^{l-m} $$
    !! using the following recursion relationships
    !!
    !!  $$d^{l+1}_{m,l+1}(\beta) = d^l_{m,l}*\sqrt{\frac{(2l+2)(2l+3)}{(l+m+1)(l-m+1)}} \cos(\frac{\beta}{2}) \sin(\frac{\beta}{2})$$
    !!  $$d^{l+1}_{m,l+1}(\beta) = d^l_{m,l}*\sqrt{\frac{(2l+2)(2l+3)}{(l+m+1)(l+m+2)}} \cos(\frac{\beta}{2})^2 $$
    !!
    !! Note the normalization factor $\Sqrt{\frac{2l+1}{2}}$ in the above equations, for unormalized dlml the equations wold change slightly
    !!
    !! It seems that iterating between the two recurrences is the most stable computation approach. No overflow/underflow occurs at least till bw=5000.

    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    !! Whether to compute normalized or unnormalized d matrices 
    logical, intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(SIZE(cos2,1),(bw*(bw+1))/2_dp),dlml_tmp(SIZE(cos2,1),bw)
    integer(kind=dp) :: l,m,swap_id,lm_id,next_lm_id
    
    ! Compute dlml values for m=0
    dlml(:,1) = MERGE(SQRT(0.5_dp),1._dp,normalized)     ! In agreement with our normalization $\sqrt{\frac{2l+1}{2}}$ at l=0
    ! Setup workspace
    dlml_tmp(:,1) = dlml(:,1)
    
    ! Do stair recursion alternating between (l -> l+1) and (l,m -> l+1,m+1 ) steps.
    do l=0,bw-2
       call dlml_recursion_m_contiguous(dlml_tmp, l, bw, sincos, cos2,normalized)
       lm_id = triangular_to_flat_index_reversed(l+1,0_dp)
       next_lm_id = triangular_to_flat_index_reversed(l+1,l+1)
       dlml(:,lm_id:next_lm_id) = dlml_tmp(:,:l+2)
    end do
  end function compute_all_dlml_m_contiguous

  subroutine wig_l_recurrence(workspace,cos,l,m1,m2,normalized)
    !
    integer(kind = dp), intent(in) :: l,m1,m2
    real(kind = dp),intent(in) :: cos(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    real(kind = dp),intent(inout) :: workspace(:,:)
    integer(kind = dp) :: dl_id,dl_1id,o,i,j
    real(kind = dp) :: c_common,c_l,c_l_1

    i = l-m2+1_dp
    dl_1id = MODULO(i-1_dp,3)+1_dp
    dl_id = MODULO(i,3)+1_dp
    o = MODULO(i+1_dp,3)+1_dp
    
       
    c_common = real(l+1_dp,kind = dp)*SQRT(1._dp/real( ((l+1_dp)**2-m1**2)*((l+1_dp)**2-m2**2),kind = dp))
    c_common = MERGE(c_common*SQRT( real(2_dp*l+3_dp,kind=dp)),c_common,normalized)
    
    c_l = real(2_dp*l+1_dp,kind = dp)
    c_l = MERGE(SQRT(c_l),c_l,normalized)

    if (l==0) then
       workspace(:,o) = c_common*c_l*cos*workspace(:,dl_id)
    else if (l==m2) then
       workspace(:,o) = c_common*c_l*( cos-(real(m1*m2,kind=dp)/real(l*(l+1_dp),kind = dp)) )*workspace(:,dl_id)
    else
       workspace(:,o) = c_common*c_l*( cos-(real(m1*m2,kind=dp)/real(l*(l+1_dp),kind = dp)) )*workspace(:,dl_id)          
       c_l_1 = SQRT(real((l**2-m1**2)*(l**2-m2**2),kind = dp))/real(l,kind = dp)
       c_l_1 = MERGE(c_l_1/SQRT(real(2_dp*l-1_dp,kind = dp)),c_l_1,normalized)
       workspace(:,o) = workspace(:,o) - c_common*c_l_1*workspace(:,dl_1id)
       !do j=1,Size(workspace,1)
       !   if (workspace(j,o)>1._dp) then
       !      print *, 'l',l
       !      print *, 'm1',m1
       !      print *, 'm2',m2
       !      print *, c_common*c_l,c_common*c_l_1
       !      print *, workspace(j,o)
       !      print *, 'beta',acos(cos(j))
       !      print *
       !   end if
       !end do
    end if
    
  end subroutine wig_l_recurrence
  function genWig_L2(m1,m2,bw,cos,dlml) result(wigners_m1m2)
    !  Given orders 0<=m1<=m2<=bw, and a bandwidth bw, this function will
    !  generate all the Wigner little d functions whose orders
    !  are (m1, m2) and degrees are j = max(|m1|, |m2|) = m2 through j = bw - 1
    !  using the 3-term recurrence. 
    !  
    !  All of these Wigners will have L2 norm = 1
    !  
    !  
    !  let j = max(|m1|, |m2|)
    !  
    !  The functions generated will be
    !  
    !  d_{m1,m2}^j, d_{m1,m2}^{j+1}, ..., d_{m1,m2}^{bw-1}
    !  
    !  Each of these functions will be evaluated at the n = 2*bw-many
    !  points
    !  
    !  pi*(2 * [0..n-1] + 1) / ( 2 * n )
    !  
    !  If beta(k) = pi*(2*k+1)/(2*n), then what's returned will be the
    !  array
    !  
    !  d_{m1,m2}^j(beta(0)) ... d_{m1,m2}^{bw-1}(beta(0))
    !  d_{m1,m2}^j(beta(1)) ... d_{m1,m2}^{bw-1}(beta(1))
    !  d_{m1,m2}^j(beta(2)) ... d_{m1,m2}^{bw-1}(beta(2)) ...
    !  d_{m1,m2}^j(beta(n-1)) ... d_{m1,m2}^{bw-1}(beta(n-1))
    !  
    !  arguments: m1, m2 = orders of the functions
    !             bw = bandwidth
    !             trig_samples = array containing cos(beta),cos(beta/2) and sine(beta) values
    !               
    !             result = array to store result, length (bw-m)*n   (as (n,bw-m) matrix)
    !                      where m = max( |m1|, |m2| ); 
    !  
    !  The routine won't be efficient, but it'll work.
    !  NOTE: compared to the python/c version the arrays are transposed:
    !        As mentioned before, here we compute the small wigner d matrices as
    !        (n,bw-m) where as in the C/python version they have the shape (bw-m,n)
    integer(kind = dp) :: m1,m2,bw
    real(kind = dp) :: cos(:),dlml(:)
    real(kind = dp) :: wigners_m1m2(Size(cos,1),(bw-max(abs(m1),abs(m2))))
    
    real(kind = dp) :: workspace(Size(cos,1),3)
    integer(kind = dp) :: i,l_start,l,n_samples,o

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
       l=l_start+(i-1_dp)
       call wig_l_recurrence(workspace,cos,l,m1,m2,.TRUE.)
       o = MODULO(i+1_dp,3)+1_dp
       wigners_m1m2(:,i+1) = workspace(:,o)
    end do
  end function genWig_L2
  function genWig_L2_trsp(m1,m2,bw,cos,dlml) result(wigners_m1m2)
    !! Same as genWig_L2 but computes the transposed wigners_m1m2 array
    !! This is faster than Transpose(genWig_L2)
    integer(kind = dp) :: m1,m2,bw
    real(kind = dp) :: cos(:),dlml(:)
    real(kind = dp) :: wigners_m1m2((bw-max(abs(m1),abs(m2))),Size(cos,1))
    
    real(kind = dp) :: workspace(Size(cos,1),3)
    integer(kind = dp) :: i,j,l_start,l,n_samples,o

    
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
       l=l_start+(i-1_dp)
       call wig_l_recurrence(workspace,cos,l,m1,m2,.TRUE.)
       o = MODULO(i+1_dp,3)+1_dp
       wigners_m1m2(i+1,:) = workspace(:,o)
    end do
  end function genWig_L2_trsp

  function compute_dmn(m1,m2,bw,betas,normalized) result(wigners_m1m2)
    integer(kind = dp) :: m1,m2,bw
    real(kind = dp) :: betas(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    real(kind = dp) :: wigners_m1m2(Size(betas,1),(bw-max(abs(m1),abs(m2)))),cos_vals(Size(betas,1)),sincos(Size(betas,1)),cos2(Size(betas,1)),dlml(Size(betas,1))
    real(kind = dp) :: workspace(Size(betas,1),3)
    integer(kind = dp) :: i,l_start,l,n_samples,o

    if (.NOT.(0<=m1 .AND. m1<=m2 .AND. m2<=bw )) then
       print *, 'Invalid arguments:  0<=m1<=m2<=bw is not sattisfied.'
    end if

    cos_vals = COS(betas)
    sincos = SIN(betas/2._dp)*COS(betas/2._dp)
    cos2 = COS(betas/2._dp)**2
    
    l_start = m2 !max(abs(m1),abs(m2))
    n_samples = size(betas,1)
    workspace(:,1) = 0
    workspace(:,2) = compute_dlml(l_start, m1, sincos, cos2, normalized)
    do i=1,Size(workspace,1)
       if (workspace(i,2)<1.0D-300) then
          workspace(i,2)=0._dp
       end if
    end do
    wigners_m1m2(:,1) = workspace(:,2)
    !print * , 'lstart', l_start,'m1m2',m1,m2
    do i=1, bw-l_start-1
       l=l_start+(i-1_dp)
       call wig_l_recurrence(workspace,cos_vals,l,m1,m2,normalized)
       o = MODULO(i+1_dp,3)+1_dp
       wigners_m1m2(:,i+1) = workspace(:,o)
    end do
  end function compute_dmn
  function wigner_l(l,betas,normalized) result(dl)
    !! compute the small wigner-d matrix $d^l_{m_1,m_2}(\beta)$ for fixed l
    integer(kind = dp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = TRUE
    real(kind = dp) :: betas2(2*SIZE(betas,1))
    real(kind = dp) :: dl(2*SIZE(betas,1),2*l+1,2*l+1)
    real(kind = dp) :: trig(2*SIZE(betas,1),3)
    real(kind = dp) :: workspace(2*SIZE(betas,1),3)
    real(kind=dp) :: dlml_workspace(2*SIZE(betas,1),l+1),sym_const_m,sym_const_n,sym_const_l,norm
    integer(kind = dp) :: m,n,i,l_start,mid,nid,neg_mid,neg_nid,o,bw,nbeta2

    dl = 0
    sym_const_l = (-1._dp)**l
    betas2(:Size(betas,1)) = betas
    betas2(Size(betas,1)+1:) = Pi-betas(Size(betas,1):1:-1)
    nbeta2 = Size(betas2,1)
    trig(:,1) = cos(betas2)                      
    trig(:,2) = cos(betas2/2_dp)*sin(betas2/2_dp) 
    trig(:,3) = cos(betas2/2_dp)**2
    
    bw = l+1_dp

    dlml_workspace(:,1) = MERGE(SQRT(0.5_dp),1._dp,normalized)
    do m=0,l
       mid = l+1_dp+m
       neg_mid = l+1_dp-m
       sym_const_m = (-1._dp)**m
       call dlml_recursion_l_contiguous(dlml_workspace, m, bw, trig(:,2), trig(:,3),normalized)
       do n=m,l
          nid = l+1_dp+n
          neg_nid = l+1_dp-n
          sym_const_n = (-1._dp)**n
          !norm = SQRT(2._dp/real(2_dp*l+1_dp,kind=dp))
          
          workspace(:,1) = 0
          workspace(:,2) = dlml_workspace(:,1+n-m)
          workspace(:,3) = workspace(:,2)
          
          do i=n,l-1
             call wig_l_recurrence(workspace,trig(:,1),i,m,n,normalized)
          end do
          
          o = MODULO(l-n+1_dp,3)+1_dp
          !remove norm from wigner d values
          workspace(:,o) = workspace(:,o)!*norm

          ! populate the small wigner d-matrix using all of the available symmetries
          dl(:,nid,mid)=workspace(:,o)          
          if ((m==0 .AND. n==0)) cycle
          dl(:,neg_mid,neg_nid) = workspace(:,o)
          if (.NOT. m==n) then
             dl(:,mid,nid) = (sym_const_m*sym_const_n)*workspace(:,o)
             dl(:,neg_nid,neg_mid) = (sym_const_m*sym_const_n)*workspace(:,o)
          end if
          if (m==0 .or. n==0) cycle
          dl(:,neg_nid,mid)=workspace(nbeta2:1:-1,o)*(sym_const_l*sym_const_m)
          dl(:,nid,neg_mid)=workspace(nbeta2:1:-1,o)*(sym_const_l*sym_const_n)
          if (m==n) cycle
          dl(:,neg_mid,nid)=workspace(nbeta2:1:-1,o)*(sym_const_l*sym_const_m)
          dl(:,mid,neg_nid)=workspace(nbeta2:1:-1,o)*(sym_const_l*sym_const_n)
          
       end do
    end do
  end function wigner_l

  function genWigAll(bw) result(wigners)
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
    
    integer(kind = dp) :: bw
    integer (kind = dp) :: n_samples, wsize
    real(kind = dp) :: betas(2*bw),trig_samples(2*bw,3)
    ! wigner should have size total_n_wigners(bw) but seems like f2py does not support array size asignments trhough pure functions ...
    real(kind=dp) :: wigners(bw*2_dp,(bw * (2_dp + 3_dp*bw + bw*bw))/6_dp)
    real(kind=dp) :: dlml_workspace(2_dp*bw,bw)
    integer(kind = dp) :: m1,m2,slize(2)
    
    wsize = size_wigner_d(bw)
    if (wsize <= 0) then
       print*, "Error: calculated size is non-positive."
       stop
    end if
    
    wigners = 0.0_dp
    dlml_workspace(:,1) = SQRT(0.5_dp)
    trig_samples = create_trig_samples(bw)
    do m1=0,bw-1
       call dlml_recursion_l_contiguous(dlml_workspace, m1, bw, trig_samples(:,2), trig_samples(:,3),.TRUE.)
       do m2=m1, bw-1
          slize = mnl_to_flat_l_slice(m1,m2,bw)
          !print*, m1,m2
          !print*, genWig_L2(m1,m2,bw,trig_samples)
          wigners(:,slize(1):slize(2)) = genWig_L2(m1,m2,bw,trig_samples(:,1),dlml_workspace(:,1+m2-m1))
       end do
    end do
  end function genWigAll
  subroutine genWigAll_preallocated(bw,wigners)
    ! Same as genWigAll but writes the wigner d matricies into the provided wigners array
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
    
    integer(kind = dp), intent(in) :: bw
    real(kind=dp), intent(inout) :: wigners(:,:)
    real(kind = dp) :: trig_samples(2*bw,3)
    real(kind=dp) :: dlml_workspace(2_dp*bw,bw)
    ! wigner should have size total_n_wigners(bw) but seems like f2py does not support array size asignments trhough pure functions ...
    integer(kind = dp) :: m1,m2,slize(2)
    
    if (.NOT. size(wigners)==size_wigner_d(bw)) then
       print '(A, I3)','Fortran Error: Provided wigners array has wrong length, it shoud be of shape:', size_wigner_d(bw) 
       stop
    end if
    
    wigners = 0.0_dp
    dlml_workspace(:,1) = SQRT(0.5_dp)
    trig_samples = create_trig_samples(bw)    
    do m1=0,bw-1
       call dlml_recursion_l_contiguous(dlml_workspace, m1, bw, trig_samples(:,2), trig_samples(:,3),.TRUE.)
       do m2=m1, bw-1
          slize = mnl_to_flat_l_slice(m1,m2,bw)
          !print*, m1,m2
          !print*, genWig_L2(m1,m2,bw,trig_samples)
          wigners(:,slize(1):slize(2)) = genWig_L2(m1,m2,bw,trig_samples(:,1),dlml_workspace(:,m2-m1+1))
       end do
    end do
  end subroutine genWigAll_preallocated

  function wigner_d_recurrence(d_l_1,l,trig_vals,sqrts) result(d_l)
    ! given a wigner little-d matrix of degree L-1 evaluated
    ! at some angle beta, construct the wigner little-d
    ! matrix of degree L, EVALUATED AT THAT SAME BETA.
    
    real(kind = dp), intent(in) :: trig_vals(:),sqrts(:),d_l_1(:)
    integer(kind = dp), intent(in)  :: l
    real(kind = dp) :: d_l((2*l+1)*(2*l+1)),temp((2*l)*(2*l)),cos_beta,sin_beta,inv_deg
    integer(kind = dp) :: deg,tmpdim,i,j
    
    cos_beta = trig_vals(1)
    sin_beta = trig_vals(2)

    if (l==0) then
       d_l = 1
    else if (l==1) then
       d_l(1) = cos_beta * sin_beta
       d_l(2) = sqrts(3) * cos_beta * sin_beta
       d_l(3) = sin_beta**2

       d_l(4) = -d_l(2)
       d_l(5) = d_l(1)-d_l(3)
       d_l(6) = d_l(2)

       d_l(7) = d_l(3)
       d_l(8) = -d_l(2)
       d_l(9) = d_l(1)
      
    else
       d_l = 0
       temp(:(2*l-1)**2) = d_l_1
       temp((2*l-1)**2+1:) = 0
       do deg=(2*l-1), 2*l
          inv_deg = (1._dp/real(deg,kind=dp))
          tmpdim = deg + 1_dp
          d_l(:tmpdim**2) = 0
          do i=1,deg
             do j=1,deg
                d_l((i*tmpdim)+j) = d_l((i*tmpdim)+j) &
                     & + inv_deg*sqrts(deg-i)*sqrts(deg-j)   * temp((i*deg)+j)*cos_beta
                d_l((i*tmpdim)+j+1_dp) = d_l((i*tmpdim)+j+1_dp) &
                     & + inv_deg*sqrts(deg-i)*sqrts(j+1_dp)  * temp((i*deg)+j)*sin_beta
                d_l(((i+1_dp)*tmpdim)+j) = d_l(((i+1_dp)*tmpdim)+j) &
                     & - inv_deg*sqrts(i+1_dp)*sqrts(deg-j)  * temp((i*deg)+j)*sin_beta
                d_l(((i+1_dp)*tmpdim)+j+1_dp) = d_l((i*tmpdim)+j) &
                     & + inv_deg*sqrts(i+1_dp)*sqrts(j+1_dp) * temp((i*deg)+j)*cos_beta
             end do
          end do
          if (deg == 2*l-1) then
             temp = d_l(:tmpdim**2)
          end if
       end do
    end if
  end function wigner_d_recurrence
end module make_wigner

module softclass
  use precision
  use utils
  use make_wigner
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'
  
  ! f2py has issues if there is no variable declaration before the type declaration
  ! for unknown reasons so dont exclude the following integer definition
  integer :: f2py_bug  = 1
  type :: so3ft
     integer(kind = dp) :: wigner_size
     integer(kind = dp) :: bw,Lmax
     real(kind = dp), allocatable :: wigner_d(:,:)
     real(kind = dp), allocatable :: wigner_d_trsp(:)
     real(kind = dp), allocatable :: wigner_dlml(:,:)
     real(kind = dp), allocatable :: trig_samples(:,:)
     real(kind = dp), allocatable :: legendre_weights(:)
     real(kind = dp), allocatable :: fft_r2c_out(:,:,:)
     complex(kind = dp), allocatable :: fft_c2c_in(:,:,:)
     complex(kind = dp), allocatable :: fft_c2r_in(:,:,:)
     complex(kind = dp), allocatable :: fft_c2c_out(:,:,:)
     integer(kind=dp) :: plan_c2c_forward,plan_c2c_backward,plan_r2c_forward,plan_c2r_backward
     integer(kind = sp) :: fftw_flags   ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
     logical :: plans_allocated_c  = .FALSE.
     logical :: plans_allocated_r = .FALSE.
   contains
     procedure :: destroy
     procedure :: destroy_fft
     procedure :: alloc_fft_arrays
     procedure :: dealloc_fft_arrays
     procedure :: init
     procedure :: init_fft
     procedure :: init_wigners

     ! transforms
     procedure :: inverse_wigner_trf_cmplx
     procedure :: inverse_wigner_loop_body_cmplx_alloc => inverse_wigner_loop_body_cmplx_alloc_
     procedure :: inverse_wigner_loop_body_cmplx => inverse_wigner_loop_body_cmplx_
     procedure :: forward_wigner_trf_cmplx
     procedure :: forward_wigner_loop_body_cmplx_alloc => forward_wigner_loop_body_cmplx_alloc_
     procedure :: forward_wigner_loop_body_cmplx => forward_wigner_loop_body_cmplx_
     procedure :: inverse_wigner_trf_real
     procedure :: inverse_wigner_loop_body_real_alloc => inverse_wigner_loop_body_real_alloc_
     procedure :: inverse_wigner_loop_body_real => inverse_wigner_loop_body_real_
     procedure :: forward_wigner_trf_real
     procedure :: forward_wigner_loop_body_real_alloc => forward_wigner_loop_body_real_alloc_
     procedure :: forward_wigner_loop_body_real => forward_wigner_loop_body_real_
     procedure :: isoft
     procedure :: soft
     procedure :: rsoft
     procedure :: irsoft
     procedure :: isoft_many
     procedure :: soft_many
     procedure :: rsoft_many
     procedure :: irsoft_many

     ! derived functions
     procedure :: integrate_over_so3_cmplx
     procedure :: integrate_over_so3_real
     procedure :: inverse_wigner_loop_body_corr_cmplx_alloc => inverse_wigner_loop_body_corr_cmplx_alloc_
     procedure :: inverse_wigner_loop_body_corr_cmplx => inverse_wigner_loop_body_corr_cmplx_
     procedure :: cross_correlation_ylm_cmplx
     procedure :: corss_correlation_ylm_cmplx_3d
     procedure :: inverse_wigner_loop_body_corr_real_alloc => inverse_wigner_loop_body_corr_real_alloc_
     procedure :: inverse_wigner_loop_body_corr_real => inverse_wigner_loop_body_corr_real_
     procedure :: cross_correlation_ylm_real
     procedure :: corss_correlation_ylm_real_3d

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
     subroutine inverse_wigner_interface(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
       import :: dp  ! Otherwise dp is undefined in _abstract interface blocks
       import :: so3ft
       class(so3ft),intent(in),target :: self
       complex(kind = dp), intent(in) :: coeff(:)
       complex(kind = dp), intent(inout) :: so3func(:,:,:)
       integer(kind=dp), intent(in) :: m1,m2
       real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
     end subroutine inverse_wigner_interface
     subroutine forward_wigner_interface(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
       import :: dp  ! Otherwise dp is undefined in _abstract interface blocks
       import :: so3ft
       class(so3ft),intent(in),target :: self
       complex(kind = dp), intent(inout) :: coeff(:)
       complex(kind = dp), intent(in) :: so3func(:,:,:)
       integer(kind=dp), intent(in) :: m1,m2
       real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
     end subroutine forward_wigner_interface
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
  subroutine destroy_fft(self,apply_to_real_fft)
    class(so3ft),intent(inout) :: self
    logical, intent(in) :: apply_to_real_fft
    
    
    if (self%plans_allocated_c .AND. (.NOT. apply_to_real_fft)) then
       call self%dealloc_fft_arrays(.TRUE.)
       call dfftw_destroy_plan(self%plan_c2c_forward)
       call dfftw_destroy_plan(self%plan_c2c_backward)
       self%plans_allocated_c = .FALSE.
    end if

    if (self%plans_allocated_r .AND. apply_to_real_fft) then
       call self%dealloc_fft_arrays(.FALSE.)
       call dfftw_destroy_plan(self%plan_r2c_forward)
       call dfftw_destroy_plan(self%plan_c2r_backward)
       self%plans_allocated_r = .FALSE.
    end if
    
  end subroutine destroy_fft
  subroutine init_wigners(self)
    class(so3ft),intent(inout) :: self
    integer(kind = dp) :: m1, m2, d_shape(2),l_slice(2),trsp_slice(2),wig_len,lid,bw
    bw = self%bw
    d_shape = wigner_d_shape(bw)
    allocate(self%wigner_d(d_shape(1),d_shape(2)))
    call genWigAll_preallocated(bw,self%wigner_d)
    
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
    integer(kind = sp) :: fft_rank,fft_n(2),fft_howmany,fft_inembed(2),fft_istride,fft_idist,fft_onembed(2),fft_ostride,fft_odist,bw2
    integer(kind = dp) :: elshape(3)

    ! Destroy FFTs if they already exist.
    call self%destroy_fft(use_real_fft)
    
    elshape = euler_shape(self%bw)
    bw2 = int(2_dp*self%bw,kind = sp)
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
       !fft_odist = product(fft_onembed)
       call self%alloc_fft_arrays(.TRUE.)
       
       call dfftw_plan_many_dft_r2c(self%plan_r2c_forward,&
            & fft_rank,fft_n,fft_howmany,&
            & self%fft_r2c_out,fft_inembed,fft_istride,fft_idist,&
            & self%fft_c2r_in,fft_onembed,fft_ostride,fft_odist,&
            & self%fftw_flags)       
       
       call dfftw_plan_many_dft_c2r(self%plan_c2r_backward,&
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
       call self%alloc_fft_arrays(.FALSE.)

       call dfftw_plan_many_dft(self%plan_c2c_forward,&
            & fft_rank,fft_n,fft_howmany,&
            & self%fft_c2c_in,fft_inembed,fft_istride,fft_idist,&
            & self%fft_c2c_out,fft_onembed,fft_ostride,fft_odist,&
            & FFTW_FORWARD,self%fftw_flags)
       call dfftw_plan_many_dft(self%plan_c2c_backward,&
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
    integer(kind = dp) :: bw2
    bw2=2_dp*self%bw
    if (use_real_fft) then
       allocate(self%fft_r2c_out(bw2,bw2,bw2))       
       allocate(self%fft_c2r_in(bw2,bw2/2+1,bw2))
    else
       allocate(self%fft_c2c_in(bw2,bw2,bw2))
       allocate(self%fft_c2c_out(bw2,bw2,bw2))       
    end if
  end subroutine alloc_fft_arrays
  subroutine dealloc_fft_arrays(self,use_real_fft)
    logical, intent(in) :: use_real_fft
    class(so3ft), intent(inout) :: self
    
    if (use_real_fft) then
       if (allocated(self%fft_c2r_in))then
          deallocate(self%fft_c2r_in)
       end if
       if (allocated(self%fft_r2c_out))then
          deallocate(self%fft_r2c_out)
       end if
    else
       if (allocated(self%fft_c2c_in))then
          deallocate(self%fft_c2c_in)
       end if
       if (allocated(self%fft_c2c_out))then
          deallocate(self%fft_c2c_out)
       end if
    end if
  end subroutine dealloc_fft_arrays
  subroutine init(self,bw,lmax,precompute_wigners,init_ffts,fftw_flags)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    logical, intent(in) :: init_ffts,precompute_wigners
    integer(kind = sp), intent(in) :: fftw_flags
    class(so3ft),intent(inout) :: self
    if ( (bw/=self%bw) .OR. (.NOT. precompute_wigners)) then
       call self%destroy()
       self%bw = bw
       allocate(self%legendre_weights(2*bw))
       self%legendre_weights = legendre_quadrature_weights(bw)
       allocate(self%trig_samples(2*bw,3))
       self%trig_samples = create_trig_samples(bw)
       self%wigner_size = size_wigner_d(bw)
       if (precompute_wigners) then
          call self%init_wigners()
       else
          allocate(self%wigner_dlml(2*bw,triangular_size(bw)))
          self%wigner_dlml = compute_all_dlml_l_contiguous(bw, self%trig_samples(:,2), self%trig_samples(:,3),.TRUE.)
       end if
    else
       if (precompute_wigners .AND. (.NOT. allocated(self%wigner_d)) ) then
          call self%init_wigners()
       end if
    end if
    
    self%lmax = lmax
    self%fftw_flags = fftw_flags
    if (init_ffts) then
       call self%init_fft(.FALSE.)
       call self%init_fft(.TRUE.)
    end if
  end subroutine init
  function init_soft(bw,lmax,precompute_wigners,init_ffts,fftw_flags) result(self)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    logical, intent(in) :: init_ffts
    logical, intent(in) :: precompute_wigners
    integer(kind = sp), intent(in) :: fftw_flags
    type(so3ft) :: self
    call self%init(bw,lmax,precompute_wigners,init_ffts,fftw_flags)
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
    call self%destroy_fft(.True.)
    call self%destroy_fft(.False.)
    self%bw=0
    self%lmax=0
  end subroutine destroy

  function get_so3func_part_halfcomplex(bw,m1,m2,so3func) result(so3func_part)
    complex(kind = dp),intent(in) :: so3func(:,:,:)
    integer(kind = dp),intent(in) :: m1,m2,bw
    complex(kind = dp) :: so3func_part(size(so3func,1))
    integer(kind = dp) :: s_ids(2)

    if (m1<0) then
       s_ids = order_to_ids(-m1,-m2,bw)
       so3func_part = CONJG(so3func(:,s_ids(1),s_ids(2)))
    else
       s_ids = order_to_ids(m1,m2,bw)
       so3func_part = so3func(:,s_ids(1),s_ids(2))
    end if
       
  end function get_so3func_part_halfcomplex

  subroutine inverse_wigner_loop_body_cmplx_alloc_(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
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
    bw2 = 2_dp*bw
    
    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))

    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    c_slice = coeff_slice(m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))
    so3func(:,m1+1,m2+1) = matmul(coeff_part,wig_mat)
    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    c_slice = coeff_slice(-m2,-m1,bw)
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
       c_slice = coeff_slice(m2,m1,bw)
       coeff_part = coeff(c_slice(1):c_slice(2))*(sym_const_m1*sym_const_m2)
       so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice(-m1,-m2,bw)
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
    c_slice = coeff_slice(m1,-m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice(-m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice(m2,-m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)


    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice(-m2,m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
  end subroutine inverse_wigner_loop_body_cmplx_alloc_
  subroutine inverse_wigner_loop_body_cmplx_(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
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
    bw2 = 2_dp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genWig_L2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    c_slice = coeff_slice(m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))
    so3func(:,m1+1,m2+1) = matmul(coeff_part,wig_mat)
    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    c_slice = coeff_slice(-m2,-m1,bw)
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
       c_slice = coeff_slice(m2,m1,bw)
       coeff_part = coeff(c_slice(1):c_slice(2))*(sym_const_m1*sym_const_m2)
       so3func(:,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice(-m1,-m2,bw)
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
    c_slice = coeff_slice(m1,-m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice(-m1,m2,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
        
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice(m2,-m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)


    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice(-m2,m1,bw)
    coeff_part = coeff(c_slice(1):c_slice(2))*sym_array(m2+1:)*sym_const_m2
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
  end subroutine inverse_wigner_loop_body_cmplx_
  subroutine inverse_wigner_trf_cmplx(self,coeff,so3func,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    integer(kind = dp) :: i,m1,m2,L
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(inverse_wigner_interface),pointer :: loop_body => Null()

    ! Select loop_body
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_cmplx_alloc_
    else
       loop_body => inverse_wigner_loop_body_cmplx_
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
       do i=1,(L+1)*(L+2)/2
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
  end subroutine inverse_wigner_trf_cmplx
  subroutine forward_wigner_loop_body_cmplx_alloc_(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
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
    bw2 = 2_dp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    l_slice = mnl_to_flat_l_slice(m1, m2, bw)
    wig_mat(1:bw2,1:bw-m2) => self%wigner_d(:,l_slice(1):l_slice(2))
    !call get_wigner_matrix(m1,m2,wig_mat,wig_mat_arr)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    c_slice = coeff_slice(m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(so3func(:,m1+1,m2+1),wig_mat)
    !return
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!

    s_ids = order_to_ids(-m2,-m1,bw)
    !print * , -m2,-m1
    !write(*,'(F16.5, F16.5)', advance='yes') so3func(:,s_ids(1),s_ids(2))
    c_slice = coeff_slice(-m2,-m1,bw)
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
       c_slice = coeff_slice(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(so3func(:,s_ids(1),s_ids(2)),wig_mat)

       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice(-m1,-m2,bw)
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
    c_slice = coeff_slice(m1,-m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice(-m1,m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m2
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice(m2,-m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice(-m2,m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m2
    coeff(c_slice(1):c_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
  end subroutine forward_wigner_loop_body_cmplx_alloc_
  subroutine forward_wigner_loop_body_cmplx_(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
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
    bw2 = 2_dp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genWig_L2(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id))

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    c_slice = coeff_slice(m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(:,m1+1,m2+1),wig_mat)
    !return
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!

    s_ids = order_to_ids(-m2,-m1,bw)
    !print * , -m2,-m1
    !write(*,'(F16.5, F16.5)', advance='yes') so3func(:,s_ids(1),s_ids(2))
    c_slice = coeff_slice(-m2,-m1,bw)
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
       c_slice = coeff_slice(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(self%legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)

       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice(-m1,-m2,bw)
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
    c_slice = coeff_slice(m1,-m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m2*sym_array(m2+1:))
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice(m2,-m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(self%legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m2*sym_array(m2+1:))
  end subroutine forward_wigner_loop_body_cmplx_
  subroutine forward_wigner_trf_cmplx(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: i,m1,m2,L
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(forward_wigner_interface),pointer :: loop_body => Null()

    ! Select loop_body
    if (allocated(self%wigner_d)) then
       loop_body => forward_wigner_loop_body_cmplx_alloc_
    else
       loop_body => forward_wigner_loop_body_cmplx_
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
  end subroutine forward_wigner_trf_cmplx

  subroutine inverse_wigner_loop_body_real_alloc_(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2
    real(kind = dp),target :: wig_mat_arr(self%bw-m2,2*self%bw)
    integer(kind=dp) :: s_ids(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw,bw2,t_slice(2)

    sym_const_m2 = (-1.0_dp)**m2
    bw = self%bw
    bw2 = 2_dp*bw

    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice(m1,m2,bw)
    so3func(:,m1+1,m2+1) = matmul(coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)),wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    c_nm2nm1_slice = coeff_slice(-m2,-m1,bw)
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
    c_pm1nm2_slice = coeff_slice(m1,-m2,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)),wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice(m2,-m1,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)),wig_mat)    
  end subroutine inverse_wigner_loop_body_real_alloc_
  subroutine inverse_wigner_loop_body_real_(self,coeff,so3func,m1,m2,sym_array,sym_const_m1)
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

    sym_const_m2 = (-1.0_dp)**m2
    bw = self%bw
    bw2 = 2_dp*bw

    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genWig_L2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice(m1,m2,bw)
    so3func(:,m1+1,m2+1) = matmul(coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)),wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    c_nm2nm1_slice = coeff_slice(-m2,-m1,bw)
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
    c_pm1nm2_slice = coeff_slice(m1,-m2,bw)
    coeff_part = coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice(m2,-m1,bw)
    coeff_part = coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2))*sym_array(m2+1:)*sym_const_m1
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = matmul(coeff_part,wig_mat)    
  end subroutine inverse_wigner_loop_body_real_
  subroutine inverse_wigner_trf_real(self,coeff,so3func,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(in) :: coeff(:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: i,m1,m2,L,s_ids(2),s_ids_sym(2),bw
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(inverse_wigner_interface),pointer :: loop_body

    bw  = self%bw    
    ! Select loop_body
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_real_alloc_
    else
       loop_body => inverse_wigner_loop_body_real_
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
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
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
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
          so3func(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(:,s_ids(1),s_ids(2)))
       end do
    end if
  end subroutine inverse_wigner_trf_real
  subroutine forward_wigner_loop_body_real_alloc_(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    complex(kind = dp), intent(inout) :: coeff(:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2
    real(kind = dp),target :: wig_mat_arr(2*self%bw,self%bw-m2)
    complex(kind = dp) :: so3func_part(2*self%bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw,bw2,l_slice(2)

    sym_const_m2 = (-1.0)**m2
    bw = self%bw
    bw2 = 2_dp*bw

    ! load wigner d matrices
    l_slice = mnl_to_flat_l_slice(m1, m2, bw)
    wig_mat(1:bw2,1:bw-m2) => self%wigner_d(:,l_slice(1):l_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice(m1,m2,bw)
    coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)) = matmul(so3func(:,m1+1,m2+1),wig_mat)

    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    so3func_part = get_so3func_part_halfcomplex(bw,-m2,-m1,so3func)
    c_nm2nm1_slice = coeff_slice(-m2,-m1,bw)
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
       c_slice = coeff_slice(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)))
       
       !! -m1,-m2 !!
       c_slice = coeff_slice(-m1,-m2,bw)
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
    c_pm1nm2_slice = coeff_slice(m1,-m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m1,m2 !!
    !! Use real symmetry
    c_slice = coeff_slice(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice(m2,-m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)) = matmul(so3func_part,wig_mat)*sym_array(m2+1:)
    
    !! -m2,m1 !!
    c_slice = coeff_slice(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)))
  end subroutine forward_wigner_loop_body_real_alloc_
  subroutine forward_wigner_loop_body_real_(self,so3func,coeff,m1,m2,sym_array,sym_const_m1)
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
    bw2 = 2_dp*bw
    
    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genWig_L2(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id))

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice(m1,m2,bw)
    coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)) = matmul(self%legendre_weights*so3func(:,m1+1,m2+1),wig_mat)
    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    so3func_part = get_so3func_part_halfcomplex(bw,-m2,-m1,so3func)
    c_nm2nm1_slice = coeff_slice(-m2,-m1,bw)
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
       c_slice = coeff_slice(m2,m1,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)))
       
       !! -m1,-m2 !!
       c_slice = coeff_slice(-m1,-m2,bw)
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
    c_pm1nm2_slice = coeff_slice(m1,-m2,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)) = matmul(self%legendre_weights*so3func_part,wig_mat)*sym_array(m2+1:)
   
    !! -m1,m2 !!
    !! Use real symmetry
    c_slice = coeff_slice(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice(m2,-m1,bw)
    so3func_part = so3func(bw2:1:-1,s_ids(1),s_ids(2))*sym_const_m1
    coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)) = matmul(self%legendre_weights*so3func_part,wig_mat)*sym_array(m2+1:)
    
    
    !! -m2,m1 !!
    c_slice = coeff_slice(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)))
  end subroutine forward_wigner_loop_body_real_
  subroutine forward_wigner_trf_real(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: i,m1,m2,L
    real(kind=dp) :: sym_const_m1,sym_array(self%bw)
    procedure(forward_wigner_interface),pointer :: loop_body

    ! Select loop_body
    if (allocated(self%wigner_d)) then
       loop_body => forward_wigner_loop_body_real_alloc_
    else
       loop_body => forward_wigner_loop_body_real_
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
  end subroutine forward_wigner_trf_real
  
  subroutine isoft(self,coeff,so3func,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    self%fft_c2c_in=0.0_dp
    call self%inverse_wigner_trf_cmplx(coeff,self%fft_c2c_in,use_mp)

    call dfftw_execute_dft(self%plan_c2c_forward,self%fft_c2c_in,so3func)
    so3func = so3func * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine isoft
  subroutine soft(self,so3func,coeff,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp

    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    self%fft_c2c_out=0.0_dp

    call dfftw_execute_dft(self%plan_c2c_backward,so3func,self%fft_c2c_out)

    self%fft_c2c_out = self%fft_c2c_out * (2.0_dp*pi/real(2_dp*self%bw,kind=dp)**2) ! * 1/(2*bw) * 2*pi/(2*bw)
    
    call self%forward_wigner_trf_cmplx(self%fft_c2c_out,coeff,use_mp)

  end subroutine soft
  subroutine irsoft(self,coeff,so3func,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(in) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    self%fft_c2r_in=0.0_dp
    call self%inverse_wigner_trf_real(coeff,self%fft_c2r_in,use_mp)

    self%fft_c2r_in = CONJG(self%fft_c2r_in) ! to correct for the fact that we have to compute the forward not the backward fft.
    call dfftw_execute_dft_c2r(self%plan_c2r_backward,self%fft_c2r_in,so3func)
    so3func = so3func * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine irsoft
  subroutine rsoft(self,so3func,coeff,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    real(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    self%fft_c2r_in=0.0_dp
    
    call dfftw_execute_dft_r2c(self%plan_r2c_forward,so3func,self%fft_c2r_in)
    self%fft_c2r_in = self%fft_c2r_in * (2.0_dp*pi/real(2_dp*self%bw,kind=dp)**2) ! * 1/(2*bw) * 2*pi/(2*bw)
    self%fft_c2r_in = CONJG(self%fft_c2r_in) ! to correct for the fact that we have to compute the backward not the forward fft.
    call self%forward_wigner_trf_real(self%fft_c2r_in,coeff,use_mp)    
  end subroutine rsoft
  subroutine soft_many(self,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    complex(kind=dp),intent(in) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(inout) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    integer(kind=dp) :: n,i
    n = size(so3funcs,1)

    if (use_mp) then       
       call self%dealloc_fft_arrays(.False.)
       !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
       call self%alloc_fft_arrays(.False.)
       !$OMP DO
       do i=1,n
          call self%soft(so3funcs(:,:,:,i),coeffs(:,i),.FALSE.)       
       end do
       !$OMP END DO
       call self%dealloc_fft_arrays(.False.)
       !$OMP END PARALLEL
       call self%alloc_fft_arrays(.False.)
    else
       do i=1,n
          call self%soft(so3funcs(:,:,:,i),coeffs(:,i),.FALSE.)       
       end do
    end if
       
  end subroutine soft_many
  subroutine isoft_many(self,coeffs,so3funcs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    complex(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(in) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    integer(kind=dp) :: n,i
    n = size(so3funcs,1)

    if (use_mp) then
       call self%dealloc_fft_arrays(.False.)
       !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
       call self%alloc_fft_arrays(.False.)
       !$OMP DO
       do i=1,n
          call self%isoft(coeffs(:,i),so3funcs(:,:,:,i),.FALSE.)       
       end do
       !$OMP END DO
       call self%dealloc_fft_arrays(.False.)
       !$OMP END PARALLEL
       call self%alloc_fft_arrays(.False.)
    else
       do i=1,n
          call self%isoft(coeffs(:,i),so3funcs(:,:,:,i),.FALSE.)       
       end do
    end if
  end subroutine isoft_many
  subroutine rsoft_many(self,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    real(kind=dp),intent(in) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(inout) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    integer(kind=dp) :: n,i
    n = size(so3funcs,1)

    if (use_mp) then
       call self%dealloc_fft_arrays(.False.)
       !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
       call self%alloc_fft_arrays(.False.)
       !$OMP DO
       do i=1,n
          call self%rsoft(so3funcs(:,:,:,i),coeffs(:,i),.FALSE.)       
       end do
       !$OMP END DO
       call self%dealloc_fft_arrays(.False.)
       !$OMP END PARALLEL
       call self%alloc_fft_arrays(.False.)
    else
       do i=1,n
          call self%rsoft(so3funcs(:,:,:,i),coeffs(:,i),.FALSE.)       
       end do
    end if       
  end subroutine rsoft_many
  subroutine irsoft_many(self,coeffs,so3funcs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    real(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(in) :: coeffs(:,:)
    logical, intent(in) :: use_mp
    integer(kind=dp) :: n,i
    !f2py integer(kind = dp) :: nthreads = 1_dp
    n = size(so3funcs,1)

    if (use_mp) then
       call self%dealloc_fft_arrays(.False.)
       !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
       call self%alloc_fft_arrays(.False.)
       !$OMP DO
       do i=1,n
          call self%irsoft(coeffs(:,i),so3funcs(:,:,:,i),.FALSE.)       
       end do
       !$OMP END DO
       call self%dealloc_fft_arrays(.False.)
       !$OMP END PARALLEL
       call self%alloc_fft_arrays(.False.)
    else
       do i=1,n
          call self%irsoft(coeffs(:,i),so3funcs(:,:,:,i),.FALSE.)       
       end do
    end if
  end subroutine irsoft_many

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

  subroutine inverse_wigner_loop_body_corr_cmplx_alloc_(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:),nm1_slice(:)
    real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn(self%bw-m2)
    integer(kind=dp) :: s_ids(2),bw,bw2,nm2_slice(2),pm2_slice(2),l_start,t_slice(2)

    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
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
  end subroutine inverse_wigner_loop_body_corr_cmplx_alloc_
  subroutine inverse_wigner_loop_body_corr_cmplx_(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2,pm1_slice(:),nm1_slice(:)
    real(kind=dp),intent(in) :: sym_array(:),wig_norm(:),sym_const_m1
    real(kind = dp) :: wig_mat(self%bw-m2,2*self%bw)
    real(kind = dp) :: sym_const_m2,sym_const_m1m2
    complex(kind = dp) :: cc_lmn(self%bw-m2)
    integer(kind=dp) :: s_ids(2),bw,bw2,nm2_slice(2),pm2_slice(2),l_start,dlml_id

    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
    pm2_slice = MLc_slice(m2,bw)
    nm2_slice = MLc_slice(-m2,bw) 
    
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genWig_L2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id))
    
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
  end subroutine inverse_wigner_loop_body_corr_cmplx_
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
    complex(kind = dp) ::  f_ml(size(f_lm)),g_ml(size(f_lm))
    integer(kind = dp) :: i,m1,m2,m,L,bw,pm1_slice(2),nm1_slice(2)    
    real(kind=dp) :: sym_const_m1,sym_array(self%bw),wig_norm(self%bw)
    procedure(wigner_corr_interface),pointer :: loop_body

    bw = self%bw
    
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_corr_cmplx_alloc_
    else
       loop_body => inverse_wigner_loop_body_corr_cmplx_
    end if
    
    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if

    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    self%fft_c2c_in = 0._dp
    
    ! initiallizing some constants
    L = self%lmax
    do i=0,L
       sym_array(i+1) = (-1.0)**i
       wig_norm = 2._dp*pi*SQRT(2._dp/real(2_dp*i+1_dp,kind = dp))
    end do
    
    if (use_mp) then
       
       !$OMP PARALLEL PRIVATE(i,l,m,m1,m2,sym_const_m1,pm1_slice,nm1_slice) SHARED(f_ml,f_lm,g_ml,g_lm,self,sym_array,wig_norm)

       ! "Transpose" input coefficient layout to be adapted to the wigner memory layout (l contiguous)   
       !$OMP DO 
       do i=1,bw**2
          call flat_to_pyramid_index(l,m,i)
          f_ml(MLc(m,l,bw)) = f_lm(LMc(l,m))
          g_ml(MLc(m,l,bw)) = g_lm(LMc(l,m))
       end do
       !$OMP END DO
       
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2
       !$OMP DO
       do i=1,(L+1)*(L+2)/2
          call flat_to_triangular_index(m1,m2,i,L)
          sym_const_m1 = (-1.0)**m1
          pm1_slice= MLc_slice(m1,bw)
          nm1_slice= MLc_slice(-m1,bw)
          call loop_body(self,f_ml,g_ml,self%fft_c2c_in,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       ! "Transpose" input coefficient layout to be adapted to the wigner memory layout (l contiguous)
       do l=0,bw-1
          do m=-l,l
             f_ml(MLc(m,l,bw)) = f_lm(LMc(l,m))
             g_ml(MLc(m,l,bw)) = g_lm(LMc(l,m))
          end do
       end do

       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2       
       do m1=0, L
          sym_const_m1 = (-1.0)**m1
          pm1_slice= MLc_slice(m1,bw)
          nm1_slice= MLc_slice(-m1,bw)
          do m2=m1, L          
             call loop_body(self,f_ml,g_ml,self%fft_c2c_in,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice,nm1_slice)
          end do
       end do
    end if

    ! Compute inverse fft
    call dfftw_execute_dft(self%plan_c2c_forward,self%fft_c2c_in,cc)
    cc = cc * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)   
  end subroutine cross_correlation_ylm_cmplx
  subroutine corss_correlation_ylm_cmplx_3d(self,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
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
    integer(kind = dp), intent(in) :: radial_limits(:)
    logical, intent(in) :: use_mp
    real(kind = dp) :: inv_radial_range,radial_step
    integer(kind=dp) :: rid

    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_c) then
       call self%init_fft(.FALSE.)
    end if
    
    cc = 0._dp
    radial_step = radial_sampling_points(2)-radial_sampling_points(1)
    inv_radial_range = 1._dp/(radial_sampling_points(radial_limits(2))-radial_sampling_points(radial_limits(1)) + radial_step)
    if (use_mp) then
       call self%dealloc_fft_arrays(.False.)
       !$OMP PARALLEL PRIVATE(rid) SHARED(f_lms,g_lms,radial_step,radial_sampling_points)
       call self%alloc_fft_arrays(.False.)
       !$omp DO reduction(+:cc)
       do rid=radial_limits(1),radial_limits(2)
          call self%cross_correlation_ylm_cmplx(f_lms(:,rid),g_lms(:,rid),self%fft_c2c_out,.FALSE.)
          cc = cc + self%fft_c2c_out*radial_sampling_points(rid)**2
       end do
       !$OMP END DO
       call self%dealloc_fft_arrays(.False.)
       !$OMP END PARALLEL
       call self%alloc_fft_arrays(.False.)       
    else
       do rid=radial_limits(1),radial_limits(2)
          call self%cross_correlation_ylm_cmplx(f_lms(:,rid),g_lms(:,rid),self%fft_c2c_out,.FALSE.)
          cc = cc+self%fft_c2c_out*radial_sampling_points(rid)**2
       end do
    end if
    cc = cc*inv_radial_range
  end subroutine corss_correlation_ylm_cmplx_3d
  subroutine inverse_wigner_loop_body_corr_real_alloc_(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
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

    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
    pm2_slice = MLr_slice(m2,bw)

    ! get wigner matrix
    t_slice = mnl_to_flat_l_slice_padded(m1,m2,bw,2*bw)
    wig_mat(1:bw-m2,1:2*bw) => self%wigner_d_trsp(t_slice(1):t_slice(2))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1m2
    so3func(:,m1+1,m2+1) = matmul(cc_lmn,wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    ! conjugation to go from m2,m1 to -m2,-m1 => this also multiplies another sym_const_m1m2 thus cancelling the already present sym_const_m1m2
    cc_lmn = wig_norm(l_start:) * CONJG(f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)))
    so3func(:,s_ids(1),s_ids(2))=CONJG(matmul(cc_lmn,wig_mat))

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
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * CONJG(g_ml(pm2_slice(1):pm2_slice(2))) * sym_const_m1 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * CONJG(g_ml(pm1_slice(1):pm1_slice(2))) * sym_const_m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)    
  end subroutine inverse_wigner_loop_body_corr_real_alloc_
  subroutine inverse_wigner_loop_body_corr_real_(self,f_ml,g_ml,so3func,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
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

    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
    pm2_slice = MLr_slice(m2,bw)

    ! get wigner matrix
    dlml_id = triangular_to_flat_index(m1,m2,bw)
    wig_mat = genWig_L2_trsp(m1,m2,bw,self%trig_samples(:,1),self%wigner_dlml(:,dlml_id))
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * g_ml(pm2_slice(1):pm2_slice(2)) * sym_const_m1m2
    so3func(:,m1+1,m2+1) = matmul(cc_lmn,wig_mat)    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    ! conjugation to go from m2,m1 to -m2,-m1 => this also multiplies another sym_const_m1m2 thus cancelling the already present sym_const_m1m2
    cc_lmn = wig_norm(l_start:) * CONJG(f_ml(pm2_slice(1):pm2_slice(2)) * g_ml(pm1_slice(1):pm1_slice(2)))
    so3func(:,s_ids(1),s_ids(2))=CONJG(matmul(cc_lmn,wig_mat))

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
    cc_lmn = wig_norm(l_start:) * f_ml(pm1_slice(1):pm1_slice(2)) * CONJG(g_ml(pm2_slice(1):pm2_slice(2))) * sym_const_m1 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm(l_start:) * f_ml(pm2_slice(1):pm2_slice(2)) * CONJG(g_ml(pm1_slice(1):pm1_slice(2))) * sym_const_m2 &
         & * sym_array(l_start:)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(cc_lmn,wig_mat)    
  end subroutine inverse_wigner_loop_body_corr_real_
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
    integer(kind = dp) :: i,m,m1,m2,L,bw,pm1_slice(2),s_ids(2),s_ids_sym(2)
    real(kind=dp) :: sym_const_m1,sym_array(self%bw),wig_norm(self%bw)
    procedure(wigner_corr_real_interface),pointer :: loop_body

    bw = self%bw
    if (allocated(self%wigner_d_trsp)) then
       loop_body => inverse_wigner_loop_body_corr_real_alloc_
    else
       loop_body => inverse_wigner_loop_body_corr_real_
    end if

    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    
    ! initiallizing some constants
    L = self%lmax
    do i=0,L
       sym_array(i+1) = (-1.0)**i
       wig_norm = 2._dp*pi*SQRT(2._dp/real(2_dp*i+1_dp,kind = dp))
    end do

    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    self%fft_c2r_in = 0._dp

    if (use_mp) then
       
       !$OMP PARALLEL PRIVATE(i,l,m,m1,m2,sym_const_m1,pm1_slice) SHARED(f_ml,g_ml,self,sym_array,wig_norm)
       
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2
       !$OMP DO
       do i=1,(L+1)*(L+2)/2
          call flat_to_triangular_index(m1,m2,i,L)
          sym_const_m1 = (-1.0)**m1
          pm1_slice= MLc_slice(m1,bw)
          call loop_body(self,f_ml,g_ml,self%fft_c2c_in,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
       end do
       !$OMP END DO

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       !$OMP DO
       do m2=1,L
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
          self%fft_c2r_in(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(self%fft_c2r_in(:,s_ids(1),s_ids(2)))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else

       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2       
       do m1=0, L
          sym_const_m1 = (-1.0)**m1
          pm1_slice= MLc_slice(m1,bw)
          do m2=m1, L          
             call loop_body(self,f_ml,g_ml,self%fft_c2r_in,m1,m2,sym_array,wig_norm,sym_const_m1,pm1_slice)
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,L
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
          self%fft_c2r_in(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(self%fft_c2r_in(:,s_ids(1),s_ids(2)))
       end do
    end if

    ! Compute inverse fft
    self%fft_c2r_in = CONJG(self%fft_c2r_in)
    call dfftw_execute_dft_c2r(self%plan_c2r_backward,self%fft_c2r_in,cc)
    cc = cc * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)   
  end subroutine cross_correlation_ylm_real
  subroutine corss_correlation_ylm_real_3d(self,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
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
    real(kind = dp), intent(inout) :: cc(:,:,:)
    real(kind = dp), intent(in) :: radial_sampling_points(:)
    integer(kind = dp), intent(in) :: radial_limits(:)
    logical, intent(in) :: use_mp
    real(kind = dp) :: inv_radial_range,radial_step
    integer(kind=dp) :: rid

    ! Make sure fft plans and matrices are allocated
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    
    cc = 0._dp
    radial_step = radial_sampling_points(2)-radial_sampling_points(1)
    inv_radial_range = 1._dp/(radial_sampling_points(radial_limits(2))-radial_sampling_points(radial_limits(1)) + radial_step)
    if (use_mp) then
       call self%dealloc_fft_arrays(.False.)
       !$OMP PARALLEL PRIVATE(rid) SHARED(f_lms,g_lms,radial_step,radial_sampling_points)
       call self%alloc_fft_arrays(.False.)
       !$omp DO reduction(+:cc)
       do rid=radial_limits(1),radial_limits(2)
          call self%cross_correlation_ylm_real(f_lms(:,rid),g_lms(:,rid),self%fft_r2c_out,.FALSE.)
          cc = cc + self%fft_r2c_out*radial_sampling_points(rid)**2
       end do
       !$OMP END DO
       call self%dealloc_fft_arrays(.False.)
       !$OMP END PARALLEL
       call self%alloc_fft_arrays(.False.)       
    else
       do rid=radial_limits(1),radial_limits(2)
          call self%cross_correlation_ylm_real(f_lms(:,rid),g_lms(:,rid),self%fft_r2c_out,.FALSE.)
          cc = cc+self%fft_r2c_out*radial_sampling_points(rid)**2
       end do
    end if
    cc = cc*inv_radial_range
  end subroutine corss_correlation_ylm_real_3d
  
  function import_fftw_wisdom(file_path) result(error_code)
    ! error_code == 0 means something went wrong in fftw
    character(len=*),intent(in) :: file_path
    integer(kind = dp) :: i
    integer(C_INT) :: error_code
    character(kind=c_char), dimension(LEN(file_path)+1) :: c_path
    ! convert fortran to C string.
    do i = 1, LEN(file_path)
       c_path(i) = file_path(i:i)
    end do
    c_path(LEN(file_path)+1) = C_NULL_CHAR
    error_code = fftw_import_wisdom_from_filename(c_path)
  end function import_fftw_wisdom
  function export_fftw_wisdom(file_path) result(error_code)
    ! error_code == 0 means something went wrong in fftw
    character(len=*),intent(in) :: file_path
    integer(kind = dp) :: i
    integer(C_INT) :: error_code
    character(kind=c_char), dimension(LEN(file_path)+1) :: c_path
    ! convert fortran to C string.
    do i = 1, LEN(file_path)
       c_path(i) = file_path(i:i)
    end do
    c_path(LEN(file_path)+1) = C_NULL_CHAR
    error_code = fftw_export_wisdom_to_filename(c_path)
  end function export_fftw_wisdom

  subroutine ylm_rotation_cmplx_(f_lm,rot_f_lm,bw,trigs,sqrts,exp_a,exp_g)
    ! Computes all wigner D matrices needed for rotating all
    ! spherical harmonic coefficients f^l_m with l<bw
    ! Euler angles follow the Z-Y-Z convention for \alpha,\beta,\gamma
    integer(kind = dp),intent(in) :: bw
    complex(kind = dp),intent(in) :: f_lm(:)
    complex(kind = dp),intent(inout) :: rot_f_lm(:)
    real(kind = dp),intent(in) :: trigs(:),sqrts(:)
    complex(kind = dp),intent(in) :: exp_a(:),exp_g(:)
    real(kind = dp) :: d_l_1((2*bw-1)**2)
    complex(kind = dp),target :: D((2*bw+1)**2)
    complex(kind = dp),pointer :: Dl(:,:)
    integer(kind = dp) :: m,n,l,nm,m_start,lm_slice(2)

    ! rotate coeff while computing wigners on the fly.
    rot_f_lm(1) = f_lm(1)
    do l=1,bw-1
       ! Compute Wigner d by recurrence
       nm = 2_dp*l+1_dp
       m_start = bw-l
       d_l_1 = wigner_d_recurrence(d_l_1,l,trigs,sqrts)
       do m=1,nm
          do n=1,nm
             D(m*nm+n) = exp_a(m_start+m)*d_l_1(m*nm+n)*exp_g(m_start+n)
          end do
       end do
       Dl(1:2*l+1,1:2*l+1) => D(:(2*l+1)**2)

       ! Do the rotation
       lm_slice = LMc(l,m)
       rot_f_lm(lm_slice(1):lm_slice(2)) = MATMUL(Dl,f_lm(lm_slice(1):lm_slice(2)))
    end do
  end subroutine ylm_rotation_cmplx_
  function ylm_rotation_cmplx(f_lm,bw,euler_angles) result(rot_f_lm)
    integer(kind = dp),intent(in) :: bw
    complex(kind = dp),intent(in) :: f_lm(:)
    real(kind = dp),intent(in) :: euler_angles(:)
    complex(kind = dp) :: rot_f_lm(size(f_lm))
    real(kind = dp) :: trigs(2),sqrts(2*bw)
    complex(kind = dp) :: exp_a(2*bw),exp_g(2*bw),za,zg
    integer(kind = dp) :: i,bw2

    ! prepare variables
    bw2 = 2_dp*bw
    trigs(1) = COS(euler_angles(2))
    trigs(2) = SIN(euler_angles(2))
    za = (0._dp,-1._dp)*euler_angles(1)
    zg = (0._dp,-1._dp)*euler_angles(3)
    do i=1,bw2
       sqrts(i) = SQRT(real(i,kind = dp))
       exp_a(i) = EXP(za*real(-bw+i,kind = dp))
       exp_g(i) = EXP(zg*real(-bw+i,kind = dp))
    end do

    call ylm_rotation_cmplx_(f_lm, rot_f_lm, bw, trigs, sqrts, exp_a, exp_g)
    
  end function ylm_rotation_cmplx
  function ylm_rotation_cmplx_many(f_lms,bw,euler_angles) result(rot_f_lms)
    integer(kind = dp),intent(in) :: bw
    complex(kind = dp),intent(in) :: f_lms(:,:)
    real(kind = dp),intent(in) :: euler_angles(:)
    complex(kind = dp) :: rot_f_lms(size(f_lms,1),size(f_lms,2))
    real(kind = dp) :: trigs(2),sqrts(2*bw)
    complex(kind = dp) :: exp_a(2*bw),exp_g(2*bw),za,zg
    integer(kind = dp) :: i,bw2

    ! prepare variables
    bw2 = 2_dp*bw
    trigs(1) = COS(euler_angles(2))
    trigs(2) = SIN(euler_angles(2))
    za = (0._dp,-1._dp)*euler_angles(1)
    zg = (0._dp,-1._dp)*euler_angles(3)
    do i=1,bw2
       sqrts(i) = SQRT(real(i,kind = dp))
       exp_a(i) = EXP(za*real(-bw+i,kind = dp))
       exp_g(i) = EXP(zg*real(-bw+i,kind = dp))
    end do

    do i=1,size(f_lms,2)
       call ylm_rotation_cmplx_(f_lms(:,i), rot_f_lms(:,i), bw, trigs, sqrts, exp_a, exp_g)
    end do
  end function ylm_rotation_cmplx_many
  
end module softclass


module py
  !! Contains versions of the type bound procedures of so3ft that can be wrapped with f2py.
  use precision
  use softclass, only: so3ft,so3ft_ptr
  implicit none
contains
  subroutine set_nthreads(nthreads)
    integer(kind = dp), intent(in) :: nthreads
    call OMP_set_num_threads(nthreads)
  end subroutine set_nthreads
  
  function py_init_soft(bw,lmax,precompute_wigners,init_ffts,fftw_flags) result(self)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    !f2py integer :: lmax = bw - 1
    logical, intent(in) :: init_ffts
    !f2py logical :: init_ffts = FALSE
    logical, intent(in) :: precompute_wigners
    !f2py logical :: precompute_wigners = FALSE
    integer(kind = sp), intent(in) :: fftw_flags
    !f2py integer :: fftw_flags = 0
    ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
    type(so3ft_ptr) :: self_ptr
    integer(kind = dp) :: self
    ALLOCATE(self_ptr%p)
    self_ptr%p = so3ft(bw,lmax,precompute_wigners,init_ffts,fftw_flags)
    self = TRANSFER(self_ptr,self)
  end function py_init_soft

  subroutine int_to_soft_pointer(self_int,self_ptr,self)
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr),intent(inout) :: self_ptr
    type(so3ft),intent(inout),pointer :: self
    self_ptr = TRANSFER(self_int,self_ptr)
    self => self_ptr%p
  end subroutine int_to_soft_pointer
  function py_get_bw(self_int) result(bw)
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = dp) :: bw
    call int_to_soft_pointer(self_int,self_ptr,self)
    bw = self%bw
  end function py_get_bw
  function py_get_lmax(self_int) result(lmax)
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = dp) :: lmax
    call int_to_soft_pointer(self_int,self_ptr,self)
    lmax = self%lmax
  end function py_get_lmax
  subroutine py_set_lmax(self_int,lmax) 
    integer(kind = dp),intent(in) :: self_int
    integer(kind = dp),intent(in) :: lmax
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    self%lmax = lmax
  end subroutine py_set_lmax
  function py_get_fftw_flags(self_int) result(fftw_flags)
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = dp) :: fftw_flags
    call int_to_soft_pointer(self_int,self_ptr,self)
    fftw_flags = self%fftw_flags
  end function py_get_fftw_flags
  subroutine py_reset(self_int,bw,lmax,precompute_wigners,init_ffts,fftw_flags)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    !f2py integer :: lmax = bw - 1
    logical, intent(in) :: init_ffts
    !f2py logical :: init_ffts = FALSE
    logical, intent(in) :: precompute_wigners
    !f2py logical :: precompute_wigners = FALSE
    integer(kind = sp), intent(in) :: fftw_flags
    !f2py integer :: fftw_flags = 0
    ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call self%init(bw,lmax,precompute_wigners,init_ffts,fftw_flags)
  end subroutine py_reset
  subroutine py_destroy(self_int)
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%destroy()

    DEALLOCATE(self)
  end subroutine py_destroy
  subroutine py_init_fft(self_int,use_real_fft)
    integer(kind = dp),intent(in) :: self_int
    logical, intent(in) :: use_real_fft
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%init_fft(use_real_fft)
  end subroutine py_init_fft

  subroutine py_inverse_wigner_trf_cmplx(self_int,coeff,so3func,use_mp)
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
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%soft(so3func,coeff,use_mp)
  end subroutine py_soft
  subroutine py_irsoft(self_int,coeff,so3func,use_mp)
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
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeff(:)
    real(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%rsoft(so3func,coeff,use_mp)
  end subroutine py_rsoft
  subroutine py_isoft_many(self_int,coeffs,so3funcs,use_mp)
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
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeffs(:,:)
    complex(kind = dp), intent(in) :: so3funcs(:,:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%soft_many(so3funcs,coeffs,use_mp)
  end subroutine py_soft_many
  subroutine py_irsoft_many(self_int,coeffs,so3funcs,use_mp)
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
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp), intent(inout) :: coeffs(:,:)
    real(kind = dp), intent(in) :: so3funcs(:,:,:,:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%rsoft_many(so3funcs,coeffs,use_mp)
  end subroutine py_rsoft_many

  function py_integrate_over_so3_cmplx(self_int,f) result(integral)
    integer(kind=dp),intent(in) :: self_int
    complex(kind = dp),intent(in) :: f(:,:,:)
    complex(kind = dp) :: integral
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    integral =  self%integrate_over_so3_cmplx(f)
  end function py_integrate_over_so3_cmplx
  function py_integrate_over_so3_real(self_int,f) result(integral)
    integer(kind=dp),intent(in) :: self_int
    real(kind = dp),intent(in) :: f(:,:,:)
    real(kind = dp) :: integral
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    integral =  self%integrate_over_so3_real(f)
  end function py_integrate_over_so3_real
  
end module py

