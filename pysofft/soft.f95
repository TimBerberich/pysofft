! Define constants
module precision
implicit none
integer(kind=8), parameter :: dp = 8
integer(kind=8), parameter :: sp = 4
end module precision

module math_constants
  use precision
  implicit none
  real(kind=dp), parameter :: pi = 4._dp*atan(1._dp)
end module math_constants

! Start of routine code
module make_wigner
  use precision
  use math_constants
  implicit none
contains
  function size_wigner_d(bw) result(wsize)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: wsize
    wsize = bw*bw * int(2_dp + 3_dp*bw + bw*bw,dp)/3_dp
  end function size_wigner_d
  
  function create_beta_samples(n) result(betas)
    ! returns n uniformly sampled angles in (0,pi)
    ! Note: these are the angles used to create Chebyshev nodes.
    integer(kind = dp) :: n,i
    real(kind = dp) :: betas(n),factor
    factor = pi / real(2_dp * n, dp)
    betas = [(real(2_dp * i - 1_dp, dp) * factor, i = 1, n)]
  end function create_beta_samples

  function create_trig_samples(bw) result(trig_samples)
    integer(kind=dp) :: bw
    real(kind = dp) :: trig_samples(2*bw,3),betas(2*bw)
    betas = create_beta_samples(2*bw)
    trig_samples(:,1) = cos(betas) ! <= Chebyshev nodes
    trig_samples(:,2) = cos(betas/2_dp)
    trig_samples(:,3) = sin(betas/2_dp)
  end function create_trig_samples
  
  function wigSpec_L2(m1,m2,sin_samples,cos_samples) result(wigner_m1m2_min_l)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! L2 normed wigner little d, WHERE THE DEGREE OF THE FUNCTION IS EQUAL
    ! TO ONE OF ITS ORDERS. This is the function to use when starting the
    ! three-term recurrence at orders (m1,m2)
    ! 
    !   arguments: m1, m2 - order of the function
    !   sinEval - sine values of evaluation pts
    !   cosEval - cosine values of the evaluation pts
    !   n - how many points there are, i.e. length of sinEval, cosEval
    !   result - where to place the result, length n
    !   
    !   
    !   Note that, by definition, since I am starting the recurrence with this
    !   function, that the degree j of the function is equal to max(abs(m1), abs(m2) ).
    !   
    !   
    !   This has been tested and stably computes Pmm functions thru bw=512
    integer(kind = dp) :: m1,m2
    real(kind = dp) :: sin_samples(:)
    real(kind = dp) :: cos_samples(size(sin_samples,1))
    real(kind = dp) :: wigner_m1m2_min_l(size(cos_samples,1))
    integer(kind = dp) ::n_samples,l,l_op,n_ls,abs_m1,abs_m2,cos_power,sin_power,sin_sign
    real(kind = dp) :: norm_factor
    integer(kind = dp) :: i
    
    
    n_samples = size(cos_samples,1)
    abs_m1 = abs(m1)
    abs_m2 = abs(m2)
    l = max(abs_m1,abs_m2)
    l_op = min(abs_m1,abs_m2)
    n_ls = l-l_op
    
    norm_factor = 1
    do i=0, n_ls-1
       norm_factor = norm_factor*sqrt(real(2_dp*l-i,dp)/real(i+1_dp,dp))
    end do
    ! need to adjust to make the L2-norm equal to 1 
    norm_factor = norm_factor * sqrt(real(2_dp*l+1_dp)/2._dp)
    
    cos_power = l + sign(1_dp,m1*m2)*l_op
    sin_power = l - sign(1_dp,m1*m2)*l_op
    !the sin_sign expression is a bit tricky
    ! -1 if (m1+m2 odd) and ( (l = |m1| & m1>=0) or (l=|m2| & m2<0) ) 
    !  1 otherwise
    sin_sign = sign(1_dp,(l-m1)*(l+m2)-Modulo((abs_m1+abs_m2),2_dp)) 
    
    do i = 1, n_samples
       wigner_m1m2_min_l(i) = norm_factor*sin_sign &
            &* sin_samples(i)**sin_power * cos_samples(i)**cos_power
    end do
    
  end function wigSpec_L2
  
  function L2_aN_so3(l,m1,m2) result(out)
    integer(kind = dp) :: l,m1,m2
    real(kind = dp) :: rl,rm1,rm2
    real(kind = dp) :: out
    rl = real(l)
    rm1= real(m1)
    rm2= real(m2)
    if (l>0) then
       out = -sqrt( (2_dp*rl+3_dp) / (2_dp*rl-1_dp) ) *&
            & (rl+1_dp) /sqrt( ((rl+1_dp)**2-rm1**2) * ((rl+1_dp)**2-rm2**2) )*&
            & sqrt((rl**2-rm1**2)*(rl**2-rm2**2))/rl
    else
       out =0
    end if
    !print*,l,m1,m2,out
  end function L2_aN_so3
  
  function L2_bN_so3(l,m1,m2) result(out)
    integer(kind = dp) :: l,m1,m2
    real(kind = dp) :: rl,rm1,rm2
    real(kind = dp) :: out
    rl = real(l)
    rm1= real(m1)
    rm2= real(m2)
    out = sqrt((2_dp*rl + 3_dp)/(2_dp*rl+1_dp)) *&
         &(rl+1_dp)*(2_dp*rl+1_dp)/sqrt(((rl+1_dp)**2-rm1**2)*((rl+1_dp)**2-m2**2))
  end function L2_bN_so3
  
  function L2_cN_so3(l,m1,m2) result(out)
    integer(kind = dp) :: l,m1,m2
    real(kind = dp) :: rl,rm1,rm2
    real(kind = dp) :: out
    rl = real(l)
    rm1= real(m1)
    rm2= real(m2)
    if (l>0) then
       out = -L2_bN_so3(l,m1,m2)*rm1*rm2/(rl*(rl+1_dp))
    else
       out = 0
    end if
  end function L2_cN_so3
  
  function genWig_L2(m1,m2,bw,trig_samples) result(wigners_m1m2)
    !  Given orders m1, m2, and a bandwidth bw, this function will
    !  generate all the Wigner little d functions whose orders
    !  are (m1, m2) and degrees are j = max(|m1|, |m2|) through j = bw - 1
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
    !        As mentioned before here we compute the small wigner d matrices as
    !        (n,bw-m) where as in the C/python version they have the shape (bw-m,n)
    integer(kind = dp) :: m1,m2,bw
    real(kind = dp) :: trig_samples(:,:)
    real(kind = dp) :: wigners_m1m2((bw-max(abs(m1),abs(m2)))*2*bw)
    
    real(kind = dp) :: workspace(2*bw,7)
    integer(kind = dp) :: i,l_start,l,n_samples
    
    l_start = max(abs(m1),abs(m2))
    n_samples = size(trig_samples,1)
    workspace(:,1) = 0
    workspace(:,2) = wigSpec_L2(m1,m2,trig_samples(:,3),trig_samples(:,2))
    wigners_m1m2(:n_samples) = workspace(:,2)
    
    do i=1, bw-l_start-1
       l=l_start+(i-1_dp)
       !print*, L2_aN_so3(l+i,m1,m2)
       workspace(:,3) = workspace(:,1) * L2_aN_so3(l,m1,m2)
       workspace(:,4) = workspace(:,2) * L2_bN_so3(l,m1,m2)
       workspace(:,5) = workspace(:,4) * trig_samples(:,1)
       workspace(:,6) = workspace(:,2) * L2_cN_so3(l,m1,m2)
       workspace(:,7) = workspace(:,3) + workspace(:,5) + workspace(:,6)
       ! workspace[:,7] now contains d_{m1,m2}^{l}
       ! storing its values to wigners_m1m2
       !print*,workspace(:,7)
       wigners_m1m2(i*n_samples+1:(i+1)*n_samples) = workspace(:,7)
       !print*,workspace(:,7)
       ! update workspace for the next iteration of the recurrence
       workspace(:,1) = workspace(:,2)
       workspace(:,2) = workspace(:,7)         
    end do
  end function genWig_L2
  
  function wigner_slice(m1,m2,bw) result(slice)
    ! For given m1,m2 computes the start and stop indices
    ! of the corresponging segment int he genWigAll output array.
    ! starting points for different m1=M1+1 m2=M1+1 is given by:
    integer(kind = dp) m1,m2,bw
    integer(kind = dp) L,m,size,slice(2)
    if (.NOT. (m1>=0 .AND. m1<=m2 .AND. m1<bw .AND. m2<bw )) then
       print *, "Invalid arguments: m1,m2 have to satisfy 0<m1<=m2<bw"
    end if
    L  = bw-1_dp
    m=max(abs(m1),abs(m2))
    size  = 2_dp*bw*(bw-m)
    !start = 1
    !do n1=0,m1-1
    !    do n2=n1,L
    !        start=start+(bw-n2)*2*bw
    !    end do
    !end do
    !do n2=m1,m2
    !    start=start+(bw-n2)*2*bw
    !end do
    !slice(2) = start -1_dp
    slice(2) = -(bw*m1*(bw*(-6_dp*L+3_dp*m1-9_dp)+3_dp*L*L+3_dp*L-m1*m1+3_dp*m1-2_dp))/3_dp
    slice(2) = slice(2) - bw*(m1-m2-1_dp)*(2_dp*bw-m1-m2)
    slice(1)=slice(2)-size+1_dp
  end function wigner_slice
  
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
    real(kind=dp) :: wigners(bw*bw * (2_dp + 3_dp*bw + bw*bw)/3)
    integer(kind = dp) :: m1,m2,slize(2)
    
    wsize = size_wigner_d(bw)
    if (wsize <= 0) then
       print*, "Error: calculated size is non-positive."
       stop
    end if
    
    wigners = 0.0_dp 
    
    n_samples = 2*bw
    betas = create_beta_samples(n_samples)
    trig_samples(:,1) = cos(betas) ! <= Chebyshev nodes
    trig_samples(:,2) = cos(betas/2_dp)
    trig_samples(:,3) = sin(betas/2_dp)
    
    do m1=0,bw-1
       do m2=m1, bw-1
          slize = wigner_slice(m1,m2,bw)
          !print*, m1,m2
          !print*, genWig_L2(m1,m2,bw,trig_samples)
          wigners(slize(1):slize(2)) = genWig_L2(m1,m2,bw,trig_samples)
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
    real(kind=dp), intent(inout) :: wigners(:)
    real(kind = dp) :: betas(2*bw),trig_samples(2*bw,3)
    ! wigner should have size total_n_wigners(bw) but seems like f2py does not support array size asignments trhough pure functions ...
    integer(kind = dp) :: m1,m2,slize(2)
    
    if (.NOT. size(wigners)==size_wigner_d(bw)) then
       print '(A, I3)','Fortran Error: Provided wigners array has wrong length, it shoud be of shape:', size_wigner_d(bw) 
       stop
    end if
    
    wigners = 0.0_dp 
    trig_samples = create_trig_samples(bw)    
    do m1=0,bw-1
       do m2=m1, bw-1
          slize = wigner_slice(m1,m2,bw)
          !print*, m1,m2
          !print*, genWig_L2(m1,m2,bw,trig_samples)
          wigners(slize(1):slize(2)) = genWig_L2(m1,m2,bw,trig_samples)
       end do
    end do
  end subroutine genWigAll_preallocated

  subroutine trsp_wigAll(bw,wigners,wigners_trsp)
    integer(kind=dp), intent(in) :: bw
    real(kind = dp), intent(in) :: wigners(:)
    real(kind=dp), intent(inout) :: wigners_trsp(:)
    integer(kind=dp) :: slice(2),wig_shape(2),m1,m2
    !(bw-max(abs(m1),abs(m2)))*2*bw
    do m1=0,bw-1
       do m2=m1, bw-1
          slice = wigner_slice(m1,m2,bw)
          wig_shape = [2*bw,bw-max(abs(m1),abs(m2))]
          !print*, m1,m2
          !print*, genWig_L2(m1,m2,bw,trig_samples)
          wigners_trsp(slice(1):slice(2)) = reshape(transpose(reshape(wigners(slice(1):slice(2)),wig_shape)),[product(wig_shape)])
       end do
    end do
  end subroutine trsp_wigAll
end module make_wigner

module so3ft_utils
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
    ! It is chosen such that in loops of the follwing type memory acces for
    ! for the coefficients happens continuousely without jumps.
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
end module so3ft_utils

module so3ft
  use precision
  use make_wigner
  use so3ft_utils
  use, intrinsic :: iso_c_binding
  implicit none
  include 'fftw3.f03'

  integer(kind = dp) :: wigner_size
  integer(kind = dp) :: bw
  real(kind = dp), allocatable, target :: wigner_d(:)
  real(kind = dp), allocatable, target :: wigner_d_trsp(:)
  real(kind = dp), allocatable, target :: trig_samples(:,:)
  real(kind = dp), allocatable, target :: legendre_weights(:)
  real(kind = dp), allocatable,target :: fft_r2c_out(:,:,:)
  complex(kind = dp), allocatable, target :: fft_c2c_in(:,:,:)
  complex(kind = dp), allocatable, target :: fft_c2r_in(:,:,:)
  complex(kind = dp), allocatable, target :: fft_c2c_out(:,:,:)
  integer(kind=dp) :: plan_c2c_forward,plan_c2c_backward,plan_r2c_forward,plan_c2r_backward
  integer(kind = sp) :: fftw_flags   ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
  logical :: plans_allocated_c,plans_allocated_r = .FALSE.
  !$OMP THREADPRIVATE(fft_r2c_out)
  !$OMP THREADPRIVATE(fft_c2c_in)
  !$OMP THREADPRIVATE(fft_c2r_in)
  !$OMP THREADPRIVATE(fft_c2c_out)

contains
  subroutine destroy_fft()
    if (allocated(fft_r2c_out)) then
       deallocate(fft_r2c_out)
    end if
    if (allocated(fft_c2c_in)) then
       deallocate(fft_c2c_in)
    end if
    if (allocated(fft_c2r_in)) then
       deallocate(fft_c2r_in)
    end if
    if (allocated(fft_c2c_out)) then
       deallocate(fft_c2c_out)
    end if
    
    if (plans_allocated_c) then
       call dfftw_destroy_plan(plan_c2c_forward)
       call dfftw_destroy_plan(plan_c2c_backward)
       plans_allocated_c = .FALSE.
    end if

    if (plans_allocated_r) then
       call dfftw_destroy_plan(plan_r2c_forward)
       call dfftw_destroy_plan(plan_c2r_backward)
       plans_allocated_r = .FALSE.
    end if
  end subroutine destroy_fft
  subroutine init_wigners()
    allocate(wigner_d(wigner_size))
    allocate(wigner_d_trsp(wigner_size))
    call genWigAll_preallocated(bw,wigner_d)
    call trsp_wigAll(bw,wigner_d,wigner_d_trsp)
  end subroutine init_wigners
  subroutine init_fft(use_real_fft)
    logical, intent(in) :: use_real_fft
    integer(kind = sp) :: fft_rank,fft_n(2),fft_howmany,fft_inembed(2),fft_istride,fft_idist,fft_onembed(2),fft_ostride,fft_odist,bw2
    integer(kind = dp) :: elshape(3),output_shape_real(3)

    elshape = euler_shape(bw)

    bw2 = int(2_dp*bw,kind = sp)
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

       allocate(fft_r2c_out(bw2,bw2,bw2))       
       allocate(fft_c2r_in(bw2,bw2/2+1,bw2))
       !wrong result but fast
       !allocate(fft_c2r_in(bw2/2+1,bw2,bw2))

       !wrong result but fast
       !fft_odist = product(fft_onembed)
       
       call dfftw_plan_many_dft_r2c(plan_r2c_forward,&
            & fft_rank,fft_n,fft_howmany,&
            & fft_r2c_out,fft_inembed,fft_istride,fft_idist,&
            & fft_c2r_in,fft_onembed,fft_ostride,fft_odist,&
            & fftw_flags)
       
       
       call dfftw_plan_many_dft_c2r(plan_c2r_backward,&
            & fft_rank,fft_n,fft_howmany,&
            & fft_c2r_in,fft_onembed,fft_ostride,fft_odist,&
            & fft_r2c_out,fft_inembed,fft_istride,fft_idist,&
            & fftw_flags)
       
       plans_allocated_r = .TRUE.
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
       allocate(fft_c2c_in(elshape(1),elshape(2),elshape(2)))
       allocate(fft_c2c_out(elshape(1),elshape(2),elshape(2)))       

       call dfftw_plan_many_dft(plan_c2c_forward,&
            & fft_rank,fft_n,fft_howmany,&
            & fft_c2c_in,fft_inembed,fft_istride,fft_idist,&
            & fft_c2c_out,fft_onembed,fft_ostride,fft_odist,&
            & FFTW_FORWARD,fftw_flags)
       call dfftw_plan_many_dft(plan_c2c_backward,&
            & fft_rank,fft_n,fft_howmany,&
            & fft_c2c_out,fft_onembed,fft_ostride,fft_odist,&
            & fft_c2c_in,fft_inembed,fft_istride,fft_idist,&
            & FFTW_BACKWARD,fftw_flags)
       
       plans_allocated_c = .TRUE.
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
  subroutine alloc_fft_arrays(use_real_fft)
    logical, intent(in) :: use_real_fft
    integer(kind = dp) :: bw2
    bw2=2_dp*bw
    if (use_real_fft) then
       allocate(fft_r2c_out(bw2,bw2,bw2))       
       allocate(fft_c2r_in(bw2,bw2/2+1,bw2))
    else
       allocate(fft_c2c_in(bw2,bw2,bw2))
       allocate(fft_c2c_out(bw2,bw2,bw2))       
    end if
  end subroutine alloc_fft_arrays
  subroutine dealloc_fft_arrays(use_real_fft)
    logical, intent(in) :: use_real_fft
    if (use_real_fft) then
       if (allocated(fft_c2r_in))then
          deallocate(fft_c2r_in)
       end if
       if (allocated(fft_r2c_out))then
          deallocate(fft_r2c_out)
       end if
    else
       if (allocated(fft_c2c_in))then
          deallocate(fft_c2c_in)
       end if
       if (allocated(fft_c2c_out))then
          deallocate(fft_c2c_out)
       end if
        
    end if       
  end subroutine dealloc_fft_arrays
  subroutine init(bandwidth,precompute_wigners,init_ffts,fftw_flags_)
    integer(kind = dp), intent(in) :: bandwidth
    logical, intent(in) :: init_ffts,precompute_wigners
    !f2py logical :: init_ffts = False
    integer(kind = sp), intent(in) :: fftw_flags_
    !f2py integer :: fftw_flags_ = 64
    ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
    integer(kind = dp) :: n_coeff
    
    if ( (bandwidth/=bw) .OR. (.NOT. precompute_wigners)) then
       call destroy()
       bw = bandwidth
       allocate(legendre_weights(2*bw))
       legendre_weights = legendre_quadrature_weights(bw)
       allocate(trig_samples(2*bw,3))
       trig_samples = create_trig_samples(bw)
       wigner_size = size_wigner_d(bandwidth)
       if (precompute_wigners) then
          call init_wigners()
       end if
    else
       if (precompute_wigners .AND. (.NOT. allocated(wigner_d)) ) then
          call init_wigners()
       end if
       call destroy_fft()
    end if
    
    fftw_flags = fftw_flags_
    if (init_ffts) then
       call init_fft(.FALSE.)
       call init_fft(.TRUE.)
    end if
  end subroutine init
  subroutine destroy()
    ! Reset soft module variables
    if (allocated(wigner_d)) then
       deallocate(wigner_d)
    end if
    if (allocated(wigner_d_trsp)) then
       deallocate(wigner_d_trsp)
    end if
    if (allocated(legendre_weights)) then
       deallocate(legendre_weights)
    end if
    if (allocated(trig_samples)) then
       deallocate(trig_samples)
    end if
    call destroy_fft()
    bw=0
  end subroutine destroy

  subroutine get_wigner_matrix(m1,m2,wig_mat_p,wig_mat_arr)
    ! This function returns the part of the wigners array that corresponds
    ! to d_m1,m2^l(beta) for all possible l and beta
    integer(kind = dp), intent(in):: m1,m2
    real(kind = dp), pointer,intent(inout) :: wig_mat_p(:,:)
    real(kind = dp),contiguous,target,intent(inout) :: wig_mat_arr(:,:)
    integer(kind = dp) slice(2)

    if (.NOT. (m1>=0 .AND. m1<=m2 .AND. m1<bw .AND. m2<bw )) then
       !Note: because of this if statement the following is true
       ! m = max(abs(m1),abs(m2)) = m2)
       print *, "Invalid arguments: m1,m2 have to satisfy 0<m1<=m2<bw"
       ERROR STOP
    end if
    
    slice = wigner_slice(m1,m2,bw)

    ! Return pointer to precoputed wigner_d matrices or
    ! if no wigner_ds are allocated compute the relevant matrix.
    if (allocated(wigner_d)) then
       ! be carefull fortran ordering has to be used
       slice = wigner_slice(m1,m2,bw)
       wig_mat_p(1:2*bw,1:bw-m2) => wigner_d(slice(1):slice(2))
    else
       wig_mat_arr = reshape(genWig_L2(m1,m2,bw,trig_samples),[2*bw,bw-m2])
       wig_mat_p(1:2*bw,1:bw-m2) => wig_mat_arr
    end if
    !call print_2d_real_pointer(wig_mat_p)
  end subroutine get_wigner_matrix
  subroutine get_wigner_matrix_trsp(m1,m2,wig_mat_p,wig_mat_arr) 
    ! This function returns the part of the wigners array that corresponds
    ! to d_m1,m2^l(beta) for all possible l and beta
    integer(kind = dp), intent(in):: m1,m2
    real(kind = dp), pointer,intent(inout) :: wig_mat_p(:,:)
    real(kind = dp),contiguous,target,intent(inout) :: wig_mat_arr(:,:)
    integer(kind = dp) slice(2)

    
    if (.NOT. (m1>=0 .AND. m1<=m2 .AND. m1<bw .AND. m2<bw )) then
       !Note: because of this if statement the following is true
       ! m = max(abs(m1),abs(m2)) = m2)
       print *, "Invalid arguments: m1,m2 have to satisfy 0<m1<=m2<bw"
       ERROR STOP
    end if

    ! Return pointer to precoputed wigner_d matrices or
    ! if no wigner_ds are allocated compute the relevant matrix.
    if (allocated(wigner_d)) then
       ! be carefull fortran ordering has to be used
       slice = wigner_slice(m1,m2,bw)
       wig_mat_p(1:bw-m2,1:2*bw) => wigner_d_trsp(slice(1):slice(2))
    else
       wig_mat_arr = transpose(reshape(genWig_L2(m1,m2,bw,trig_samples),[2*bw,bw-m2]))
       wig_mat_p(1:bw-m2,1:2*bw) => wig_mat_arr
    end if
    !call print_2d_real_pointer(wig_mat_p)
  end subroutine get_wigner_matrix_trsp
  function get_wigner_matrix_copy(m1,m2,bw) result(wig_mat)
    ! This function returns the part of the wigners array that corresponds
    ! to d_m1,m2^l(beta) for all possible l and beta
    integer(kind = dp), intent(in):: m1,m2,bw
    integer(kind = dp) L,m,slice(2)
    real(kind = dp) :: wig_mat(2*bw,bw-max(abs(m1),abs(m2)))
    
    if (.NOT. (m1>=0 .AND. m1<=m2 .AND. m1<bw .AND. m2<bw )) then
       print *, "Invalid arguments: m1,m2 have to satisfy 0<m1<=m2<bw"
    end if
    
    L  = bw-1_dp
    m=max(abs(m1),abs(m2))

    slice = wigner_slice(m1,m2,bw)

    ! be carefull fortran ordering has to be used
    wig_mat = reshape(wigner_d(slice(1):slice(2)),[2*bw,bw-m])
  end function get_wigner_matrix_copy
  
  function get_so3func_part_halfcomplex(m1,m2,so3func) result(so3func_part)
    complex(kind = dp),intent(in) :: so3func(:,:,:)
    integer(kind = dp),intent(in) :: m1,m2
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
  
  subroutine inverse_wigner_loop_body_cmplx(coeff,so3func,m1,m2,sym_array,sym_const_m1,sym_const_m2)
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1,sym_const_m2
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp),target :: wig_mat_arr(bw-m2,2*bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),m,bw2

    bw2 = 2_dp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    
    call get_wigner_matrix_trsp(m1,m2,wig_mat,wig_mat_arr)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    c_slice = coeff_slice(m1,m2,bw)
    so3func(:,m1+1,m2+1) = matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    c_slice = coeff_slice(-m2,-m1,bw)
    so3func(:,s_ids(1),s_ids(2)) = matmul(coeff(c_slice(1):c_slice(2)),wig_mat)

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       c_slice = coeff_slice(m2,m1,bw)
       so3func(:,s_ids(1),s_ids(2)) = (sym_const_m1*sym_const_m2)*matmul(coeff(c_slice(1):c_slice(2)),wig_mat)       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice(-m1,-m2,bw)
       so3func(:,s_ids(1),s_ids(2)) = (sym_const_m1*sym_const_m2)*matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
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
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
    
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice(-m1,m2,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice(m2,-m1,bw)
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice(-m2,m1,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
  end subroutine inverse_wigner_loop_body_cmplx
  subroutine inverse_wigner_trf_cmplx(coeff,so3func)
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind = dp) :: i,m1,m2,m,L,s_slice(2),c_slice(2)
    real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)
    
    ! initiallizing some constants
    L = bw-1
    do i=0,L
       sym_array(i+1) = (-1.0)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform        
    do m1=0, L
       sym_const_m1 = (-1.0)**m1
       do m2=m1, L          
          sym_const_m2 = (-1.0)**m2
          call inverse_wigner_loop_body_cmplx(coeff,so3func,m1,m2,sym_array,sym_const_m1,sym_const_m2)
       end do
    end do
  end subroutine inverse_wigner_trf_cmplx
  subroutine forward_wigner_loop_body_cmplx(so3func,coeff,m1,m2,sym_array,sym_const_m1,sym_const_m2)
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp),contiguous, intent(in) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1,sym_const_m2
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp),target :: wig_mat_arr(2*bw,bw-m2) ! possible source of optimization can make it allocatable 
    integer(kind=dp) :: s_ids(2),c_slice(2),bw2
    
    bw2 = 2_dp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    
    call get_wigner_matrix(m1,m2,wig_mat,wig_mat_arr)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    c_slice = coeff_slice(m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(legendre_weights*so3func(:,m1+1,m2+1),wig_mat)
    !return
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!

    s_ids = order_to_ids(-m2,-m1,bw)
    !print * , -m2,-m1
    !write(*,'(F16.5, F16.5)', advance='yes') so3func(:,s_ids(1),s_ids(2))
    c_slice = coeff_slice(-m2,-m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)
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
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)

       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice(-m1,-m2,bw)
       coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*matmul(legendre_weights*so3func(:,s_ids(1),s_ids(2)),wig_mat)
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
    coeff(c_slice(1):c_slice(2)) = matmul(legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m2*sym_array(m2+1:))
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice(m2,-m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = matmul(legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m2*sym_array(m2+1:))
  end subroutine forward_wigner_loop_body_cmplx
  subroutine forward_wigner_trf_cmplx(so3func,coeff)
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp),contiguous, intent(in) :: so3func(:,:,:)
    integer(kind = dp) :: i,m1,m2,m,L,s_slice(2),c_slice(2)
    real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)

    ! initiallizing some constants
    L = bw-1
    do i=0,L
       sym_array(i+1) = (-1)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform        
    do m1=0, L
       sym_const_m1 = (-1.0)**m1
       do m2=m1, L
          sym_const_m2 = (-1.0)**m2
          call forward_wigner_loop_body_cmplx(so3func,coeff,m1,m2,sym_array,sym_const_m1,sym_const_m2)
       end do
    end do    
  end subroutine forward_wigner_trf_cmplx
  subroutine forward_wigner_loop_body_real(so3func,coeff,m1,m2,sym_array,sym_const_m1,sym_const_m2)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    complex(kind = dp), intent(inout) :: coeff(:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1,sym_const_m2
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp),target :: wig_mat_arr(2*bw,bw-m2)
    complex(kind = dp) :: so3func_part(2*bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw2

    bw2 = 2_dp*bw
  
    call get_wigner_matrix(m1,m2,wig_mat,wig_mat_arr)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    c_pm1pm2_slice = coeff_slice(m1,m2,bw)
    coeff(c_pm1pm2_slice(1):c_pm1pm2_slice(2)) = matmul(legendre_weights*so3func(:,m1+1,m2+1),wig_mat)

    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    so3func_part = get_so3func_part_halfcomplex(-m2,-m1,so3func)
    c_nm2nm1_slice = coeff_slice(-m2,-m1,bw)
    coeff(c_nm2nm1_slice(1):c_nm2nm1_slice(2)) = matmul(legendre_weights*so3func_part,wig_mat)

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
    coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2)) = matmul(legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m1,m2 !!
    !! Use real symmetry
    c_slice = coeff_slice(-m1,m2,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm1nm2_slice(1):c_pm1nm2_slice(2))) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    c_pm2nm1_slice = coeff_slice(m2,-m1,bw)
    coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)) = matmul(legendre_weights*so3func(bw2:1:-1,s_ids(1),s_ids(2)),wig_mat)*(sym_const_m1*sym_array(m2+1:))
    
    !! -m2,m1 !!
    c_slice = coeff_slice(-m2,m1,bw)
    coeff(c_slice(1):c_slice(2)) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_pm2nm1_slice(1):c_pm2nm1_slice(2)))
  end subroutine forward_wigner_loop_body_real
  subroutine forward_wigner_trf_real(so3func,coeff)
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    integer(kind = dp) :: i,m1,m2,m,L,s_slice(2),c_slice(2)
    real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)
    
    ! initiallizing some constants
    L = bw-1
    do i=0,L
       sym_array(i+1) = (-1)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform        
    do m1=0, L
       sym_const_m1 = (-1.0)**m1
       do m2=m1, L
          sym_const_m2 = (-1.0)**m2          
          call forward_wigner_loop_body_real(so3func,coeff,m1,m2,sym_array,sym_const_m1,sym_const_m2)
       end do
    end do    
  end subroutine forward_wigner_trf_real
  subroutine inverse_wigner_loop_body_real(coeff,so3func,m1,m2,sym_array,sym_const_m1,sym_const_m2)
    ! This subroutine assumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1,sym_const_m2
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp),target :: wig_mat_arr(bw-m2,2*bw)
    complex(kind = dp) :: so3func_part(2*bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),c_pm1pm2_slice(2),c_nm2nm1_slice(2),c_pm1nm2_slice(2),c_pm2nm1_slice(2),bw2

    bw2 = 2_dp*bw
  
    call get_wigner_matrix_trsp(m1,m2,wig_mat,wig_mat_arr)
    
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
  end subroutine inverse_wigner_loop_body_real
  subroutine inverse_wigner_trf_real(coeff,so3func)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(in) :: coeff(:)
    integer(kind = dp) :: i,m1,m2,m,L,s_ids(2),s_ids_sym(2),c_slice(2)
    real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)
    
    ! initiallizing some constants
    L = bw-1
    do i=0,L
       sym_array(i+1) = (-1)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform        
    do m1=0, L
       sym_const_m1 = (-1.0)**m1
       do m2=m1, L
          sym_const_m2 = (-1.0)**m2
          call inverse_wigner_loop_body_real(coeff,so3func,m1,m2,sym_array,sym_const_m1,sym_const_m2)
       end do
    end do

    ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
    do m2=1,L
       s_ids = order_to_ids(0_dp,m2,bw)
       s_ids_sym = order_to_ids(0_dp,-m2,bw)
       so3func(:,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(:,s_ids(1),s_ids(2)))
    end do
  end subroutine inverse_wigner_trf_real
  
  subroutine isoft(coeff,so3func)
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)

    if (.NOT. plans_allocated_c) then
       call init_fft(.FALSE.)
    end if
    fft_c2c_in=0.0_dp
    !print * , fft_c2c_in
    call inverse_wigner_trf_cmplx(coeff,fft_c2c_in)
    !print * , fft_c2c_in
    call dfftw_execute_dft(plan_c2c_forward,fft_c2c_in,so3func)
    so3func = so3func * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine isoft
  subroutine soft(so3func,coeff)
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)

    if (.NOT. plans_allocated_c) then
       call init_fft(.FALSE.)
    end if


    call dfftw_execute_dft(plan_c2c_backward,so3func,fft_c2c_out)

    fft_c2c_out = fft_c2c_out * (2.0_dp*pi/real(2_dp*bw,kind=dp)**2) ! * 1/(2*bw) * 2*pi/(2*bw)
    !write(*,'(F16.5, F16.5)', advance='yes') reshape(fft_c2c_out,[(2*bw)**3])
    !print *,'fft done'
    call forward_wigner_trf_cmplx(fft_c2c_out,coeff)

  end subroutine soft
  subroutine irsoft(coeff,so3func)
    complex(kind = dp), intent(in) :: coeff(:)
    real(kind = dp), intent(inout) :: so3func(:,:,:)

    
    if (.NOT. plans_allocated_r) then
       call init_fft(.TRUE.)
    end if
    fft_c2r_in=0.0_dp
    call inverse_wigner_trf_real(coeff,fft_c2r_in)

    fft_c2r_in = CONJG(fft_c2r_in) ! to correct for the fact that we have to compute the forward not the backward fft.
    call dfftw_execute_dft_c2r(plan_c2r_backward,fft_c2r_in,so3func)
    so3func = so3func * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine irsoft
  subroutine rsoft(so3func,coeff)
    complex(kind = dp), intent(inout) :: coeff(:)
    real(kind = dp), intent(in) :: so3func(:,:,:)

    if (.NOT. plans_allocated_r) then
       call init_fft(.TRUE.)
    end if
    fft_c2r_in=0.0_dp
    
    call dfftw_execute_dft_r2c(plan_r2c_forward,so3func,fft_c2r_in)
    fft_c2r_in = fft_c2r_in * (2.0_dp*pi/real(2_dp*bw,kind=dp)**2) ! * 1/(2*bw) * 2*pi/(2*bw)
    fft_c2r_in = CONJG(fft_c2r_in) ! to correct for the fact that we have to compute the backward not the forward fft.
    !write(*,'(F16.5, F16.5)', advance='yes') reshape(fft_c2r_in,[(2*bw)**2*(bw+1)])
    !print * , 'rfft done'
    call forward_wigner_trf_real(fft_c2r_in,coeff)    
  end subroutine rsoft

  subroutine soft_many(so3funcs,coeffs,nthreads)
    !f2py threadsafe
    complex(kind=dp),intent(in) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(inout) :: coeffs(:,:)
    integer(kind = dp), intent(in) :: nthreads
    integer(kind=dp) :: n,i
    !f2py integer(kind = dp) :: nthreads = 1_dp
    n = size(so3funcs,1)

    !set number of used threads
    call OMP_set_num_threads(nthreads)

    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
    if (nthreads>1) call alloc_fft_arrays(.False.)
    !$OMP DO
    do i=1,n
       call soft(so3funcs(:,:,:,i),coeffs(:,i))       
    end do
    !$OMP END DO
    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP END PARALLEL
    if (nthreads>1) call alloc_fft_arrays(.False.)
  end subroutine soft_many
  subroutine isoft_many(coeffs,so3funcs,nthreads)
    !f2py threadsafe
    complex(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(in) :: coeffs(:,:)
    integer(kind = dp), intent(in) :: nthreads
    integer(kind=dp) :: n,i
    !f2py integer(kind = dp) :: nthreads = 1_dp
    n = size(so3funcs,1)

    !set number of used threads
    call OMP_set_num_threads(nthreads)

    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
    if (nthreads>1) call alloc_fft_arrays(.False.)
    !$OMP DO
    do i=1,n
       call isoft(coeffs(:,i),so3funcs(:,:,:,i))       
    end do
    !$OMP END DO
    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP END PARALLEL
    if (nthreads>1) call alloc_fft_arrays(.False.)
  end subroutine isoft_many
  subroutine rsoft_many(so3funcs,coeffs,nthreads)
    !f2py threadsafe
    real(kind=dp),intent(in) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(inout) :: coeffs(:,:)
    integer(kind = dp), intent(in) :: nthreads
    integer(kind=dp) :: n,i
    !f2py integer(kind = dp) :: nthreads = 1_dp
    n = size(so3funcs,1)

    !set number of used threads
    call OMP_set_num_threads(nthreads)

    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
    if (nthreads>1) call alloc_fft_arrays(.False.)
    !$OMP DO
    do i=1,n
       call rsoft(so3funcs(:,:,:,i),coeffs(:,i))       
    end do
    !$OMP END DO
    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP END PARALLEL
    if (nthreads>1) call alloc_fft_arrays(.False.)
  end subroutine rsoft_many
  subroutine irsoft_many(coeffs,so3funcs,nthreads)
    !f2py threadsafe
    real(kind=dp),intent(inout) :: so3funcs(:,:,:,:)
    complex(kind=dp),intent(in) :: coeffs(:,:)
    integer(kind = dp), intent(in) :: nthreads
    integer(kind=dp) :: n,i
    !f2py integer(kind = dp) :: nthreads = 1_dp
    n = size(so3funcs,1)

    !set number of used threads
    call OMP_set_num_threads(nthreads)

    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP PARALLEL PRIVATE(i) SHARED(so3funcs,coeffs,n)
    if (nthreads>1) call alloc_fft_arrays(.False.)
    !$OMP DO
    do i=1,n
       call irsoft(coeffs(:,i),so3funcs(:,:,:,i))       
    end do
    !$OMP END DO
    if (nthreads>1) call dealloc_fft_arrays(.False.)
    !$OMP END PARALLEL
    if (nthreads>1) call alloc_fft_arrays(.False.)
  end subroutine irsoft_many
  
  subroutine fft(f1,f2)
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft(plan_c2c_forward,f1,f2)
    f2 = f2 * (1/(2.0_dp*real(bw,kind=dp))) !* 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine fft
  subroutine ifft(f1,f2)
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft(plan_c2c_backward,f1,f2)
    f2 = f2 * (1/(2.0_dp*real(bw,kind=dp))) !* 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine ifft
  subroutine rfft(f1,f2)
    real(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft_r2c(plan_r2c_forward,f1,f2)
    f2 = f2 * (1/(2.0_dp*real(bw,kind=dp))) !* 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine rfft
  subroutine irfft(f1,f2)
    complex(kind = dp), intent(in) :: f1(:,:,:)
    real(kind = dp), intent(inout) :: f2(:,:,:)    
    call dfftw_execute_dft_c2r(plan_c2r_backward,f1,f2)
    f2 = f2 * (1/(2.0_dp*real(bw,kind=dp))) !* 1/(2*bw) * (2*bw)/(2*pi)
  end subroutine irfft

end module so3ft

! Debug tools
module debug
  use precision
  use make_wigner
  use so3ft_utils
  use, intrinsic :: iso_c_binding
  implicit none

contains
  function my_matmul(a,b,n,m) result(c)
    complex(kind=dp) :: a(n),c(m)
    real(kind=dp) :: b(n,m)
    integer(kind=dp) :: n,m
    c=matmul(a,b)
  end function my_matmul
  subroutine print_ms(bw)
    integer(kind = dp),intent(in) :: bw
    integer(kind=dp) :: L,m1,m2,c
    c=0
    L = bw-1
    do m1=0, L
       do m2=m1, L
          print *, m1,m2
          c=c+1
          if (m1 ==0 .AND. m2 ==0) cycle

          print * ,-m2,-m1
          c=c+1
          if (.NOT. m1==m2) then
             print * , m2,m1
             c=c+1
             print * , -m1,-m2
             c=c+1
          end if

          if (m1==0 .or. m2==0) cycle       ! prevents sign swaps on 0 ids which are already covered

          print *, m1,-m2
          c=c+1
          print *, -m1,m2
          c=c+1
          
          if (m1 == m2) cycle               ! prevents duplicates due to swapping equal numbers
          print *, -m2,m1
          c=c+1
          print *, m2,-m1
          c=c+1
       end do
    end do

    print *
    print *,c
    print *

    do m1=-L, L
       do m2=-L, L
          print*, m1,m2
       end do
    end do
    print *
    print *,(2*L+1)**2
    
  end subroutine print_ms
  subroutine print_coeff_ids(bw)
      integer(kind = dp),intent(in) :: bw
      integer(kind = dp) :: i,m1,m2,m,L,s_slice(2),c_slice(2),sum,sum_tot
      real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)
      sum = 0
      ! initiallizing some constants
      L = bw-1
      do i=0,L
         sym_array(i+1) = (-1)**i 
      end do
      sum_tot = 0
      sum = 0
      
      ! non-fft part of the SO(3) fourier transform        
      do m1=0, L
         sym_const_m1 = (-1.0)**m1
         do m2=m1, L
            print *,''
            sym_const_m2 = (-1.0)**m2
            sum_tot = sum_tot + sum
            sum = 0
            c_slice = coeff_slice(m1,m2,bw)           
            print * ,c_slice
            sum = sum+1
            
            if (m1 ==0 .AND. m2 ==0) cycle    ! prevents m1=m2=0 from beeing evaluated twice
            
            c_slice = coeff_slice(-m2,-m1,bw)
            print * ,c_slice
            sum = sum+1
            
            if (.NOT. m1==m2) then
               !!  m2,m1  !!       
               c_slice = coeff_slice(m2,m1,bw)
               print * ,c_slice
               sum = sum+1
               
               !! -m1,-m2 !!
               c_slice = coeff_slice(-m1,-m2,bw)
               print * ,c_slice
               sum = sum+1
            end if
            
            if (m1==0 .or. m2==0) cycle       ! prevents sign swaps on 0 ids which are already covered
            
            !! m1,-m2 !!
            c_slice = coeff_slice(m1,-m2,bw)
            print * ,c_slice
            sum = sum+1
            !! -m1,m2 !!
            c_slice = coeff_slice(-m1,m2,bw)
            print * ,c_slice
            sum = sum+1
            
            if (m1 == m2) cycle               ! prevents duplicates due to swapping equal numbers
            
            !! m2,-m1 !!
            c_slice = coeff_slice(m2,-m1,bw)
            print * ,c_slice
            sum = sum+1
            !! -m2,m1 !!
            c_slice = coeff_slice(-m2,m1,bw)
            print * ,c_slice
            sum = sum+1
            
         end do
      end do
      sum_tot = sum_tot + sum
      print * , sum,bw-m2
      print * , 'total' , sum_tot
    end subroutine print_coeff_ids
end module debug

module so3ft_debug
  use precision
  use make_wigner
  use so3ft_utils
  use so3ft
  use, intrinsic :: iso_c_binding
  implicit none
  
contains
  !subroutine plot_matmul_wig(m1,m2)
  !  integer(kind = dp), intent(in):: m1,m2
  !  integer(kind=dp) :: i,j,sz
  !  real(kind = dp), pointer :: wig_mat(:,:),wig_mat_trsp(:,:),res(:,:)
  !  real(kind = dp),allocatable,target :: res_arr(:,:),weighted_wigmat(:,:),coeff(:),coeff2(:),so3fct(:)
  !  wig_mat_trsp => get_wigner_matrix_trsp(m1,m2)
  !  wig_mat => get_wigner_matrix(m1,m2)
  !  sz = size(matmul(wig_mat_trsp,wig_mat),1)
  !  print *, shape(wig_mat)
  !  print *, shape(wig_mat_trsp)
  !  allocate(res_arr(sz,sz))
  !  allocate(weighted_wigmat(size(wig_mat,1),size(wig_mat,2)))
  !  allocate(coeff(size(wig_mat,2)))
  !  allocate(coeff2(size(wig_mat,2)))
  !  allocate(so3fct(size(wig_mat,1)))
  !  
  !  coeff = 1.0
  !  so3fct = matmul(coeff,wig_mat_trsp)
  !  
  !  do j=1,size(wig_mat,2)
  !     do i=1, size(wig_mat,1)
  !        weighted_wigmat(i,j)= wig_mat(i,j)*legendre_weights(i)
  !     end do
  !  end do
  !  
  !  coeff2 = matmul(so3fct,weighted_wigmat)
  !  res_arr = matmul(wig_mat_trsp,weighted_wigmat)
  !  res => res_arr
  !  call print_2d_real_pointer(res)
  !  print * , coeff2
  !  deallocate(res_arr)
  !  deallocate(weighted_wigmat)
  !  deallocate(coeff)
  !  deallocate(coeff2)
  !  deallocate(so3fct)
  !end subroutine plot_matmul_wig
  !subroutine plot_wig_matrix_pointer(m1,m2,trsp)
  !  integer(kind = dp), intent(in):: m1,m2
  !  logical, intent(in) :: trsp
  !  integer(kind=dp) :: i,j
  !  real(kind = dp), pointer :: wig_mat(:,:)
  !  if (trsp) then
  !     wig_mat => get_wigner_matrix_trsp(m1,m2)
  !  else
  !     wig_mat => get_wigner_matrix(m1,m2)
  !  end if
  !  call print_2d_real_pointer(wig_mat)    
  !end subroutine plot_wig_matrix_pointer
  
  subroutine inverse_wigner_loop_body_cmplx2(so3func,coeff,m1,m2,sym_array,sym_const_m1,sym_const_m2)
      complex(kind = dp), intent(inout) :: so3func(:,:,:)
      complex(kind = dp), intent(in) :: coeff(:)
      integer(kind=dp), intent(in) :: m1,m2
      real(kind=dp),intent(in) :: sym_array(:),sym_const_m1,sym_const_m2
      real(kind = dp),pointer :: wig_mat(:,:)
      real(kind = dp),target :: wig_mat_arr(bw-m2,2*bw)
      integer(kind=dp) :: s_ids(2),c_slice(2),m,bw2
      
      bw2 = 2_dp*bw
      
      ! This method assiumes 0<=m1<=m2<=bw
      ! which also means m = max(abs(m1),abs(m2)) = m2
      
      call get_wigner_matrix_trsp(m1,m2,wig_mat,wig_mat_arr)
      
      !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
      !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
      !! m1,m2 !!
      !s_slice = sample_slice(m1,m2,bw)
      c_slice = coeff_slice(m1,m2,bw)
      so3func(m1+1,m2+1,:) = matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
      
      if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
      
      !! -m2,-m1 !!
      !s_slice = sample_slice(-m2,-m1,bw)
      s_ids=order_to_ids(-m2,-m1,bw)
      c_slice = coeff_slice(-m2,-m1,bw)
      so3func(s_ids(1),s_ids(2),:) = matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
      
      !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
      !! uses the symmetries:                                          !!
      !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
      !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
      !! They cause a constant sign swap by (-1)^(m1-m2)               !!
      if (.NOT. m1==m2) then
         !!  m2,m1  !!
         s_ids=order_to_ids(m2,m1,bw)
         c_slice = coeff_slice(m2,m1,bw)
         so3func(s_ids(1),s_ids(2),:) = (sym_const_m1*sym_const_m2)*matmul(coeff(c_slice(1):c_slice(2)),wig_mat)       
         !! -m1,-m2 !!
         s_ids=order_to_ids(-m1,-m2,bw)
         c_slice = coeff_slice(-m1,-m2,bw)
         so3func(s_ids(1),s_ids(2),:) = (sym_const_m1*sym_const_m2)*matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
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
      so3func(s_ids(1),s_ids(2),bw2:1:-1) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
      
      !! -m1,m2 !!
      s_ids=order_to_ids(-m1,m2,bw)
      c_slice = coeff_slice(-m1,m2,bw)
      so3func(s_ids(1),s_ids(2),bw2:1:-1) = sym_const_m2*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
      
      if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
      
      !! m2,-m1 !!
      s_ids=order_to_ids(m2,-m1,bw)
      c_slice = coeff_slice(m2,-m1,bw)
      so3func(s_ids(1),s_ids(2),bw2:1:-1) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
      
      !! -m2,m1 !!
      s_ids=order_to_ids(-m2,m1,bw)
      c_slice = coeff_slice(-m2,m1,bw)
      so3func(s_ids(1),s_ids(2),bw2:1:-1) = sym_const_m2*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
    end subroutine inverse_wigner_loop_body_cmplx2
  subroutine inverse_wigner_trf_cmplx2(coeff,so3func)
      complex(kind = dp), intent(in) :: coeff(:)
      complex(kind = dp),contiguous,intent(inout) :: so3func(:,:,:)
      real(kind = dp), pointer :: wig_mat(:,:)
      integer(kind = dp) :: i,m1,m2,m,L,s_slice(2),c_slice(2)
      real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)
      
      ! initiallizing some constants
      L = bw-1
      do i=0,L
         sym_array(i+1) = (-1.0)**i 
      end do
      
      ! non-fft part of the SO(3) fourier transform        
      do m1=0, L
         sym_const_m1 = (-1.0)**m1
         do m2=m1, L          
            sym_const_m2 = (-1.0)**m2
            call inverse_wigner_loop_body_cmplx2(so3func,coeff,m1,m2,sym_array,sym_const_m1,sym_const_m2)
         end do
      end do
    end subroutine inverse_wigner_trf_cmplx2

  subroutine inverse_wigner_loop_body_cmplx3(coeff,so3func,m1,m2,sym_array,sym_const_m1,sym_const_m2)
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind=dp), intent(in) :: m1,m2
    real(kind=dp),intent(in) :: sym_array(:),sym_const_m1,sym_const_m2
    real(kind = dp),pointer :: wig_mat(:,:)
    real(kind = dp),target :: wig_mat_arr(bw-m2,2*bw)
    integer(kind=dp) :: s_ids(2),c_slice(2),m,bw2

    bw2 = 2_dp*bw
    
    ! This method assiumes 0<=m1<=m2<=bw
    ! which also means m = max(abs(m1),abs(m2)) = m2
    
    call get_wigner_matrix_trsp(m1,m2,wig_mat,wig_mat_arr)
    
    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
    !! m1,m2 !!
    !s_slice = sample_slice(m1,m2,bw)
    c_slice = coeff_slice_legacy(m1,m2,bw)
    so3func(:,m1+1,m2+1) = matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice
    
    !! -m2,-m1 !!
    !s_slice = sample_slice(-m2,-m1,bw)
    s_ids=order_to_ids(-m2,-m1,bw)
    c_slice = coeff_slice_legacy(-m2,-m1,bw)
    so3func(:,s_ids(1),s_ids(2)) = matmul(coeff(c_slice(1):c_slice(2)),wig_mat)

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids=order_to_ids(m2,m1,bw)
       c_slice = coeff_slice_legacy(m2,m1,bw)
       so3func(:,s_ids(1),s_ids(2)) = (sym_const_m1*sym_const_m2)*matmul(coeff(c_slice(1):c_slice(2)),wig_mat)       
       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       c_slice = coeff_slice_legacy(-m1,-m2,bw)
       so3func(:,s_ids(1),s_ids(2)) = (sym_const_m1*sym_const_m2)*matmul(coeff(c_slice(1):c_slice(2)),wig_mat)
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
    c_slice = coeff_slice_legacy(m1,-m2,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
    
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    c_slice = coeff_slice_legacy(-m1,m2,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    c_slice = coeff_slice_legacy(m2,-m1,bw)
    so3func(bw2:1:-1,s_ids(1),s_ids(2)) = sym_const_m1*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    c_slice = coeff_slice_legacy(-m2,m1,bw)
    so3func(bw2:1:-1, s_ids(1),s_ids(2)) = sym_const_m2*matmul(sym_array(m2+1:)*coeff(c_slice(1):c_slice(2)),wig_mat)
  end subroutine inverse_wigner_loop_body_cmplx3
  subroutine inverse_wigner_trf_cmplx3(coeff,so3func)
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    integer(kind = dp) :: i,m1,m2,m,L,s_slice(2),c_slice(2)
    real(kind=dp) :: sym_const_m1,sym_const_m2,sym_array(bw)
    
    ! initiallizing some constants
    L = bw-1
    do i=0,L
       sym_array(i+1) = (-1.0)**i 
    end do
    
    ! non-fft part of the SO(3) fourier transform        
    do m1=0, L
       sym_const_m1 = (-1.0)**m1
       do m2=m1, L          
          sym_const_m2 = (-1.0)**m2
          call inverse_wigner_loop_body_cmplx3(coeff,so3func,m1,m2,sym_array,sym_const_m1,sym_const_m2)
       end do
    end do
  end subroutine inverse_wigner_trf_cmplx3
end module so3ft_debug
