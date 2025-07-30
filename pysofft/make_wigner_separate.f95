module precision
    implicit none
    integer(kind=8), parameter :: dp = 8
end module precision

module math_constants
    use precision
    real(kind=dp), parameter :: pi = 4._dp*atan(1._dp)
end module math_constants

module make_wigner
    use precision
    use math_constants
    implicit none
contains
    pure function total_n_wigners(bw) result(size)
        integer(kind = dp), intent(in) :: bw
        integer(kind = dp) :: size
        size = bw*bw * (2_dp + 3_dp*bw + bw*bw)/3
    end function total_n_wigners

    function create_beta_samples(n) result(betas)
        ! returns n uniformly sampled angles in (0,pi)
        ! Note: these are the angles used to create Chebyshev nodes.
        integer(kind = dp) :: n,i
        real(kind = dp) :: betas(n),factor
        factor = pi / real(2_dp * n, dp)
        betas = [(real(2_dp * i - 1_dp, dp) * factor, i = 1, n)]
    end function create_beta_samples

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
        integer(kind = dp) :: ifcond
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

    function L2_bN_so3(l,m1,m2) result(out)       !
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
        integer(kind = dp) m1,m2,bw,n1,n2
        integer(kind = dp) L,m,N,size,slice(2)
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
      integer(kind = dp) :: n_samples,size
      real(kind = dp) :: betas(2*bw),trig_samples(2*bw,3)
      ! wigner should have size total_n_wigners(bw) but seems like f2py does not support array size asignments trhough pure functions ...
      real(kind=dp) :: wigners(bw*bw * (2_dp + 3_dp*bw + bw*bw)/3)
      integer(kind = dp) :: m1,m2,m,slize(2)
      
      size = total_n_wigners(bw)
      if (size <= 0) then
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
end module make_wigner
