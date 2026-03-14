!> ------
!! @brief Defines kind aliases sp,dp,qp
module precision
implicit none
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)
end module precision

!> ------
!! @brief Defines pi
module math_constants
  use precision
  implicit none
  real(kind=dp), parameter :: pi = 4.0_dp*atan2(1.0_dp,1.0_dp)
end module math_constants

!> ---------
!! @brief Module containing mostly indexing helper function.
module utils
  use precision
  use math_constants
  implicit none
contains
  !> ----------
  !! @brief $\mathrm{bw}-\max(|m\\_1|,|m\\_2|)$
  !!
  !! Returns the number of small wigner-d coefficients for fixed m1, m2 and bandwidth bw,i.e.
  !!
  !! $$bw-\max(|m\\_1|,|m\\_2|)$$
  !!
  function num_coeffs(m1,m2,bw) result(num_coeff)
    ! not used 
    integer(kind=dp) :: m1,m2,bw,num_coeff
    if (abs(m1)>=bw .OR. abs(m2)>=bw) then
       print *, "Out of bounds error (n_coeffs): |m1|>bw or |m2|>bw"
       stop
    end if
    num_coeff = bw - max(abs(m1),abs(m2))
  end function num_coeffs
  !> -----------
  !! @brief $(4\mathrm{bw}^3-\mathrm{bw})/3$
  !!
  !! Returns the total number of Wigner coefficients $f^l\\_{m,n}$
  !! computed in the SO(3) fourier transform of bandwidth `bw`, i.e.
  !!
  !! $$\sum\\_{l=0}^{\mathrm{bw}-1} (2l+1)^2 =  \frac{(4\mathrm{bw}^3-\mathrm{bw})}{3} $$
  function total_num_coeffs(bw) result(num_coeff)
    ! not used 
    integer(kind=8) bw,num_coeff
    num_coeff = (4_dp*(bw*bw*bw)-bw)/3_dp
  end function total_num_coeffs

  !> -----------
  !! @brief Returns slice into $f(\alpha,\beta,\gamma)$ needed to compute $f\\_{m\\_1,m\\_2}^l$ .
  !!
  !! For order $m\\_1$, $m\\_2$ and bandwidth $\mathrm{bw}$, returns the slice in the 3-D,
  !! fft'd array of the SO(3) function $f(\alpha,\beta,\gamma)$ that is needed  to compute
  !! $f\\_{m\\_1,m\\_2}^l$ for all possible $l$.
  !!
  function order_to_ids(m1,m2,bw) result(ids)
    ! Note: That the array for $f$ is of size (2*bw,2*bw,2*bw}), so I
    ! need to multiply by that bw2 (I'm sampling at twice the bandwidth)
    integer(kind=dp) :: m1,m2,bw,bw2
    integer(kind=dp) :: ids(2)
    if (abs(m1)>=bw .OR. abs(m2)>=bw) then
       print *, "Out of bounds error (orders_to_ids): |m1|>bw or |m2|>bw"
       stop
    end if
    bw2=2*bw
    ids(1) = (1_dp-sign(1_dp,m1))/2_dp*bw2 + m1+1_dp
    ids(2) = (1_dp-sign(1_dp,m2))/2_dp*bw2 + m2+1_dp
  end function order_to_ids

  !> ---------
  !! @brief $n\\_1,n\\_2,\mathrm{bw} \rightarrow$ slice into mnl ordered $f^l\\_{n\\_1,n\\_2}$ array.
  !!
  !! Returns a slice into the coefficient array $f^l\\_{n\\_1,n\\_2}$ for all possible l valaues.
  !! It is chosen such that in loops of the follwing type memory access
  !! for the coefficients is contigous without jumps.
  !!
  !! /// info | Code
  !! ```Fortran
  !! do m1=0,bw-1
  !!   do m2=m1,bw-1
  !!     print *, coeff_slice_mnl(m1,m2,bw)
  !!     if (m1==0 .and. m2==0), cycle
  !!     print *, coeff_slice_mnl(-m2,-m1,bw)
  !!     if (m1/=m2) then
  !!        print *, coeff_slice_mnl(m2,m1,bw)
  !!        print *, coeff_slice_mnl(-m1,-m2,bw)
  !!     end if
  !!     if (m1==0), cycle
  !!     print *, coeff_slice_mnl(m1,-m2,bw)
  !!     print *, coeff_slice_mnl(-m1,m2,bw)
  !!     if (m1==m2), cycle
  !!     print *, coeff_slice_mnl(m2,-m1,bw)
  !!     print *, coeff_slice_mnl(-m2,m1,bw)
  !!   end do
  !! end do
  !! ```
  !! ///
  function coeff_slice_mnl(n1,n2,bw) result(slice)
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
  end function coeff_slice_mnl
  function coeff_location_mnl(m1,m2,l,bw) result(id)
    ! This function returns the index of the coefficient array corresponding to f_{m1,m2}^l.
    ! Note that this index has to lie within the slice given by coefLoc_so3.
    integer(kind = dp) :: m1,m2,bw,l,m,slice(2),id
    slice = coeff_slice_mnl(m1,m2,bw)
    m =  max(abs(m1),abs(m2))
    if (l<m .OR. l>bw) then
       print *, "Out of bounds error (coeff_location): l does not satisfy max(|m1|,|m2|) <= l < bw "
       stop
    end if
    id = slice(1)+l-m
  end function coeff_location_mnl
  function get_coeff_degrees(bw) result(lmn)
    !! Returns an array containing all valid l,m,n coefficient combinations
    !! in the native order they are stored in.
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: lmn((4_dp*(bw*bw*bw)-bw)/3_dp,3)
    integer(kind = dp) :: m1,m2,cslice(2),ls(bw),l
    do l=0,bw-1
       ls(l+1)=l
    end do
    
    do m1 = 0, bw-1
       do m2 = m1, bw-1
          cslice = coeff_slice_mnl(m1,m2,bw)
          lmn(cslice(1):cslice(2),1) = ls(m2+1:)
          lmn(cslice(1):cslice(2),2) = m1
          lmn(cslice(1):cslice(2),3) = m2

          if ((m1==0) .and. (m2==0)) cycle

          cslice = coeff_slice_mnl(-m2,-m1,bw)
          lmn(cslice(1):cslice(2),1) = ls(m2+1:)
          lmn(cslice(1):cslice(2),2) = -m2
          lmn(cslice(1):cslice(2),3) = -m1

          if (m1/=m2) then
             cslice = coeff_slice_mnl(m2,m1,bw)
             lmn(cslice(1):cslice(2),1) = ls(m2+1:)
             lmn(cslice(1):cslice(2),2) = m2
             lmn(cslice(1):cslice(2),3) = m1

             cslice = coeff_slice_mnl(-m1,-m2,bw)
             lmn(cslice(1):cslice(2),1) = ls(m2+1:)
             lmn(cslice(1):cslice(2),2) = -m1
             lmn(cslice(1):cslice(2),3) = -m2
          end if

          if ((m1==0) .or. (m2==0)) cycle
          
          cslice = coeff_slice_mnl(m1,-m2,bw)
          lmn(cslice(1):cslice(2),1) = ls(m2+1:)
          lmn(cslice(1):cslice(2),2) = m1
          lmn(cslice(1):cslice(2),3) = -m2
          
          cslice = coeff_slice_mnl(-m1,m2,bw)
          lmn(cslice(1):cslice(2),1) = ls(m2+1:)
          lmn(cslice(1):cslice(2),2) = -m1
          lmn(cslice(1):cslice(2),3) = m2

          if (m1==m2) cycle

          cslice = coeff_slice_mnl(m2,-m1,bw)
          lmn(cslice(1):cslice(2),1) = ls(m2+1:)
          lmn(cslice(1):cslice(2),2) = m2
          lmn(cslice(1):cslice(2),3) = -m1

          cslice = coeff_slice_mnl(-m2,m1,bw)
          lmn(cslice(1):cslice(2),1) = ls(m2+1:)
          lmn(cslice(1):cslice(2),2) = -m2
          lmn(cslice(1):cslice(2),3) = m1
       end do
    end do
  end function get_coeff_degrees
 
  function coeff_slice_lmn(l,n1,n2) result(slice)
    ! Returns the slice of coefficients f_(m,n)^l for all possible symmetry variants of m,n.
    ! It is chosen such that in loops of the follwing type, memory access
    ! for the coefficients is continuouse without jumps.
    !
    ! do l=0,bw-1
    !    do m=0,l
    !       do n=m,l
    !          slice = coeff_slice_lmn(l,m,n)
    !          f(slice(1))      !f_{m,n} 
    !          if m==0 and n==0, cycle  
    !          f(slice(1)+1)    !f_{-n,-m}
    !          if m/=n
    !             f(slice(1)+2) !f_{n,m}
    !             f(slice(1)+3) !f_{-m,-n}
    !          if m==0 or n==0, cycle
    !          o=MERGE(2,0,m/=n)
    !          f(slice(1)+2+o)  !f_{m,-n}
    !          f(slice(1)+3+o)  !f_{-m,n}
    !          if m==n, cycle
    !          f(slice(1)+6)    !f_{n,-m}
    !          f(slice(1)+7)    !f_{-n,m}
    !
    ! Insight is 
    ! 1 4 4 4 4    top to bottom is m, left to right is n.
    !   4 8 8 8    Numbers are number of function occurances in the above loop
    !     4 8 8
    !       4 8
    !         4

    integer(kind=dp), intent(in) :: n1,n2,l
    integer(kind=dp) :: slice(2)
    integer(kind=dp) :: previouse_l,previouse_m,previouse_n,length ,m,n
    logical :: m_not_0, m_g_1,n_g_m
    m = MIN(ABS(n1),ABS(n2))
    n = MAX(ABS(n1),ABS(n2))
    
    m_not_0 = (m/=0)
    m_g_1 = (m>1)
    n_g_m = (n > m)
    
    ! sum_0^(l-1) (2*l+1)^2
    previouse_l = int((l*(4_dp*l**2-1_dp))/3,dp) 
    
    previouse_m = MERGE(1_dp+l*4_dp + MERGE(4*(m-1)*(2*l+1-m),0_dp,m_g_1),0_dp,m_not_0)
    previouse_n = MERGE(1_dp+4_dp*(n-m-1_dp) + MERGE(3_dp+4_dp*(n-m-1_dp),0_dp,m_not_0),0_dp,n_g_m)
    length = MERGE(1_dp,MERGE(4_dp,8_dp,(m==n) .OR. (m==0)),n==0)
    
    slice(1) = previouse_l+previouse_m+previouse_n
    slice(2) = slice(1)+length
    slice(1) = slice(1) + 1_dp ! 1 indexing
  end function coeff_slice_lmn
  function coeff_location_lmn(l,n1,n2) result(id)
        integer(kind=dp), intent(in) :: n1,n2,l
        integer(kind=dp) :: id,slice(2),m,n
        
        n=MAX(ABS(n1),ABS(n2))
        m=MIN(ABS(n1),ABS(n2))
        slice = coeff_slice_lmn(l,m,n)
        if (n1==m .AND. n2==n) then
           id = slice(1) 
        else if (n1==-n .AND. n2==-m) then
           id = slice(1) + 1_dp
        else if (n1==n .AND. n2==m) then
           id = slice(1) + 2_dp
        else if (n1==-m .AND. n2 == -n) then
           id = slice(1) + 3_dp
        else if (n1==m .AND. n2 == -n) then
           if (m/=n) then
              id = slice(1) + 4_dp
           else
              id = slice(1) + 2_dp
           end if
        else if (n1==-m .AND. n2 == n) then
           if (m/=n) then
              id = slice(1) + 5_dp
           else
              id = slice(1) + 3_dp
           end if
        else if (n1==n .AND. n2 == -m) then
           id = slice(1) + 6_dp
        else ! (n1==-n .AND. n2 == m) then
           id = slice(1) + 7_dp
        end if
  end function coeff_location_lmn
  function get_coeff_degrees_risbo(bw) result(lmn)
    !! Returns an array containing all valid l,m,n coefficient combinations
    !! in the native order they are stored in.
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: lmn((4_dp*(bw*bw*bw)-bw)/3_dp,3)
    integer(kind = dp) :: m1,m2,start,cslice(2),l,id,o
    do l = 0, bw-1
       do m1 = 0, l
          do m2 = m1, l
             cslice = coeff_slice_lmn(l,m1,m2)
             start = cslice(1)
             lmn(start,1) = l
             lmn(start,2) = m1
             lmn(start,3) = m2

             if ((m1==0) .and. (m2==0)) cycle
             id = start+1_dp
             lmn(id,1) = l
             lmn(id,2) = -m2
             lmn(id,3) = -m1

             if (m1/=m2) then
                id = start+2_dp
                lmn(id,1) = l
                lmn(id,2) = m2
                lmn(id,3) = m1

                id = start+3_dp
                lmn(id,1) = l
                lmn(id,2) = -m1
                lmn(id,3) = -m2
             end if

             if ((m1==0) .or. (m2==0)) cycle
             
             o=MERGE(2_dp,0_dp,m1/=m2)
             id = start+2_dp+o
             lmn(id,1) = l
             lmn(id,2) = m1
             lmn(id,3) = -m2
             
             id = start+3_dp+o
             lmn(id,1) = l
             lmn(id,2) = -m1
             lmn(id,3) = m2
             
             if (m1==m2) cycle
             
             id = start+6_dp
             lmn(id,1) = l
             lmn(id,2) = m2
             lmn(id,3) = -m1
             
             id = start+7_dp
             lmn(id,1) = l
             lmn(id,2) = -m2
             lmn(id,3) = m1
          end do
       end do
    end do
  end function get_coeff_degrees_risbo
  
  function euler_shape(bw) result(eshape)
    integer(kind=dp), intent(in) :: bw
    integer(kind=dp) :: eshape(3)
    eshape = 2*bw
  end function euler_shape

  !> ----------
  !! @brief Computes even legendre weights
  !!
  !! This function computes the even legendre weights
  !! used in the quadrature of the SOFT algorithm
  !! For a fixed bandwith (bw) the weights are the
  !! unique solutions to the problem:
  !!
  !! $$\sum\\_{k=0}^{2\text{bw}-1} w\\_{\text{bw}}(k) P\\_m(\cos(\beta\\_k)) = \delta\\_{0,m} \quad \forall 0 \leq m \leq \text{bw}$$ 
  !!
  !! where the sampling angles are given by $\beta\\_k = \frac{\pi(2k+1)}{4\text{bw}}$
  !! Their closed form expression is given by:
  !!
  !! $$ w\\_{\text{bw}}(k) = \frac{2}{\text{bw}}\sin\left(\frac{\pi(2k+1)}{4\text{bw}}\right) \sum\\_{j=0}^{\text{bw}-1}\frac{1}{2j+1}\sin\left((2k+1)(2j+1)\frac{\pi}{4\text{bw}}\right)$$
  !!
  !! both of these Formulas can be found in equations (2.13) and (2.14) of  
  !! :cite: P.J. Kostelec and D.N. Rockmore, J Fourier Anal Appl (2008) 14: 145–179  
  !! The proof of the closed form is contained in  
  !! :cite: Driscoll, J.R., Healy, D.:  Proc. 34th IEEE FOCS (1989), pp. 344–349. Adv. in Appl. Math., vol. 15, pp. 202–250 (1994)  
  function legendre_quadrature_weights(bw) result(weights)
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
  function get_empty_coeff_many(bw,n) result(coeff)
    integer(kind = dp) :: bw,n
    complex(kind = dp) :: coeff((4_dp*(bw*bw*bw)-bw)/3_dp,n)
    coeff = 0.0
  end function get_empty_coeff_many
  function get_empty_so3func_cmplx(bw) result(so3func)
    integer(kind = dp) :: bw
    complex(kind = dp) :: so3func(2*bw,2*bw,2*bw)
    so3func = 0.0
  end function get_empty_so3func_cmplx
  function get_empty_so3func_cmplx_many(bw,n) result(so3func)
    integer(kind = dp) :: bw,n
    complex(kind = dp) :: so3func(2*bw,2*bw,2*bw,n)
    so3func = 0.0
  end function get_empty_so3func_cmplx_many
  function get_empty_so3func_halfcmplx(bw) result(so3func)
    integer(kind = dp) :: bw
    complex(kind = dp) :: so3func(2*bw,bw+1,2*bw)
    so3func = 0.0
  end function get_empty_so3func_halfcmplx
  function get_empty_so3func_halfcmplx_many(bw,n) result(so3func)
    integer(kind = dp) :: bw,n
    complex(kind = dp) :: so3func(2*bw,bw+1,2*bw,n)
    so3func = 0.0
  end function get_empty_so3func_halfcmplx_many
  function get_empty_so3func_real(bw) result(so3func)
    integer(kind = dp) ::  bw
    real(kind = dp) :: so3func(2*bw,2*bw,2*bw)
    so3func = 0.0
  end function get_empty_so3func_real
  function get_empty_so3func_real_many(bw,n) result(so3func)
    integer(kind = dp) ::  bw,n
    real(kind = dp) :: so3func(2*bw,2*bw,2*bw,n)
    so3func = 0.0
  end function get_empty_so3func_real_many
  !> --------
  !! @brief Enforces real symmetry in mnl ordered coefficients
  !!
  !! Enforces the symmetry $f^l\\_{m,n}^* = f^l\\_{-m,-n} (-1)^{m+n}$ in
  !! m,n,l ordered coefficients
  !!
  subroutine enforce_real_sym_mnl(coeff,bw)
    complex(kind=dp) ,intent(inout) :: coeff(:)
    integer(kind=dp) ,intent(in) :: bw
    real(kind=dp) :: sym_const_m1,sym_const_m2
    integer(kind=dp) :: m1,m2,c_slice(2),c_slice_sym(2)

    c_slice = coeff_slice_mnl(0_dp,0_dp,bw)
    coeff(c_slice(1):c_slice(2)) = coeff(c_slice(1):c_slice(2))%re
    do m1=0,bw-1
       sym_const_m1 = (-1._dp)**m1
       do m2=-bw+1,bw-1
          sym_const_m2 = (-1._dp)**m2
          c_slice_sym = coeff_slice_mnl(m1,m2,bw)
          c_slice = coeff_slice_mnl(-m1,-m2,bw)
          coeff(c_slice(1):c_slice(2)) = sym_const_m1*sym_const_m2*CONJG(coeff(c_slice_sym(1):c_slice_sym(2)))       
       end do
    end do
  end subroutine enforce_real_sym_mnl

  !> --------
  !! @brief many version of enforce_real_sym
  !!
  !! Enforces the symmetry $f^l\\_{m,n}^* = f^l\\_{-m,-n} (-1)^{m+n}$ in
  !! m,n,l ordered coefficients.
  !! 
  !! coeff_many is now a 2D array where the second index labels differenct coefficient arrays.
  subroutine enforce_real_sym_mnl_many(coeff_many,bw)
    complex(kind=dp) ,intent(inout) :: coeff_many(:,:)
    integer(kind=dp) ,intent(in) :: bw
    integer(kind=dp) :: n

    do n=1,SIZE(coeff_many,2)
       call enforce_real_sym_mnl(coeff_many(:,n),bw)
    end do
  end subroutine enforce_real_sym_mnl_many
  !> --------
  !! @brief Enforces real symmetry in lmn ordered coefficients
  !!
  !! Enforces the symmetry $f^l\\_{m,n}^* = f^l\\_{-m,-n} (-1)^{m-n}$ in
  !! l,m,n ordered coefficients
  !!
  subroutine enforce_real_sym_lmn(coeff,bw)
    complex(kind=dp) ,intent(inout) :: coeff(:)
    integer(kind=dp) ,intent(in) :: bw
    real(kind=dp) :: sym_const_m1,sym_const_m2
    integer(kind=dp) :: m1,l,m2,c_slice(2),c_loc,o

    do l=0,bw-1_dp
       do m1=0,l
          sym_const_m1 = (-1.0)**m1
          do m2=m1,l
             sym_const_m2 = (-1.0)**m2
             c_slice = coeff_slice_lmn(l,m1,m2)
             c_loc = c_slice(1)
             if (m2==0) then
                coeff(c_loc) = coeff(c_loc)%re
             end if
             if (m1==0 .AND. m2==0) cycle
             
             if (m1/=m2) then
                !! m2,m1 !!
                !coeff(c_loc+2_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 1_dp))
                coeff(c_loc+1_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 2_dp))
                !! -m1,-m2 !!
                coeff(c_loc+3_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc)) 
             else
                coeff(c_loc+1_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc))
             end if
             if (m1==0) cycle
             o=MERGE(2_dp,0_dp,(m1/=m2))

             !! m1,-m2 !!
             !! -m1,m2 !!
             coeff(c_loc + o +3_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc+o+2_dp))
             if (m1==m2) cycle
             !! m2,-m1 !!
             !! -m2,m1 !!
             coeff(c_loc + 7_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc+6_dp))
          end do
       end do
    end do
  end subroutine enforce_real_sym_lmn

  !> --------
  !! @brief many version of enforce_real_sym_lmn
  !!
  !! Enforces the symmetry $f^l\\_{m,n}^* = f^l\\_{-m,-n} (-1)^{m+n}$ in
  !! l,m,n ordered coefficients.
  !! 
  !! coeff_many is now a 2D array where the second index labels differenct coefficient arrays.
  subroutine enforce_real_sym_lmn_many(coeff_many,bw)
    complex(kind=dp) ,intent(inout) :: coeff_many(:,:)
    integer(kind=dp) ,intent(in) :: bw
    integer(kind=dp) :: n
    do n=1,SIZE(coeff_many,2)
       call enforce_real_sym_mnl(coeff_many(:,n),bw)
    end do
  end subroutine enforce_real_sym_lmn_many
  
  
  ! Indexing tricks section
  !< --------------
  !! 
  function triangular_size(bw) result(tri_size)
    !! Consider the triangular index 0<=i<=j<bw
    !! This function returns the total number of possible pairs i,j
    integer(kind = dp),intent(in) :: bw
    integer(kind = dp) :: tri_size
    tri_size = (bw*(bw+1))/2_dp
  end function triangular_size
  function pyramid_size(bw) result(pyr_size)
    !! consider the pyramid index 0<=i<bw and -i<=j<=i
    !!This function returns the total number of possible pairs i,j
    integer(kind = dp),intent(in) ::bw
    integer(kind = dp) :: pyr_size
    pyr_size = bw**2
  end function pyramid_size
  subroutine flat_to_pyramid_index(i,j,k)
    ! Converts a running index k=0 to bw*(bw+1)/2-1 into
    ! a triangular double index i,j with  0<=i<bw and i<=j<=bw
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
  !> ---------------------
  !! @brief Converts a 1d index k to a triagonal double index (i,j).
  !!
  !! Converts a running index k=0 to (N+1)*(N+2)/2-1 into
  !! a triangular double index i,j with  0<=i<kw and i<=j<=bw
  !! This allows to reformulate Triangular loops as simple loops via
  !!\code{.fortran}
  !! do i=0,N
  !!   do j=i,N
  !!     print*,i,j   
  !!   end do
  !! end do
  !!\endcode
  !! with 
  !!\code{.fortran}
  !! do k=1,((N+1)*(N+2))/2_dp
  !!   triangular_to_flat(i,j,k,N)
  !!   print*,i,j    
  !! end do
  !!\endcode
  subroutine flat_to_triangular_index(i,j,k,N)
    integer(kind = dp), intent(inout) :: i,j
    integer(kind = dp), intent(in) :: k,N
    i = int(-SQRT(real((2_dp*N+3_dp)**2-8_dp*(k-1_dp),kind = dp))+real(2_dp*N+3_dp,kind=dp),kind=dp)/2_dp
    j =  (k-1_dp)-i*(2_dp*N+1_dp-i)/2_dp
    !i = (-int(SQRT(real((2_dp*N+3_dp)**2-8_dp*(k-1_dp),kind = dp)),kind = dp)+2_dp*N+3_dp)/2_dp
    !j =  (k-1_dp)-i*(2_dp*N+1_dp-i)/2_dp
  end subroutine flat_to_triangular_index
  !> ------------------
  !! @brief Converts triangular index (i,j) to 1d index.
  !!
  !! Consider the triangular index $0<=i<=j<\text{bw}$.
  !! 
  !! This function converts an index (i,j) to a one dimensional index.
  !! For fixed $i$ the differenct possible values of $j$ are contiguous in the one dimesional representation. 
  function triangular_to_flat_index(i,j,bw) result(id)
    integer(kind = dp), intent(in) :: i,j,bw 
    integer(kind = dp) :: id 
    id = (i*(2_dp*bw-i+1_dp))/2_dp + j-i + 1_dp
  end function triangular_to_flat_index
  !> Consider the triangular index 0<=i<=j<bw 
  !! This function returns the slice of of a flattened array that corresponds
  !! to all valid i for a fixed j.
  !! This function allows to store a triangular array in an i contiguous way
  function triangular_to_flat_index_reversed(j,i) result(id)
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

  !> ----------
  !! Computes number of complex spherical harmonic coefficients for
  !! a given order bandwidth bw.
  function n_LMc(bw) result(n)
    integer(kind=dp), intent(in) :: bw
    integer(kind=dp) :: n
    n=bw*bw
  end function n_LMc

  !> ----------
  !! Computes number of real spherical harmonic coefficients for
  !! a given order bandwidth bw.
  function n_LMr(bw) result(n)
    integer(kind=dp), intent(in) :: bw
    integer(kind=dp) :: n
    n= (bw*(bw+1_dp))/2_dp
  end function n_LMr
  
  !> ----------
  !! Index for the spherical harmonic coefficient Y_lm for complex data.
  !! This is the same convention used by the software SHTNS
  !! https://nschaeff.bitbucket.io/shtns/index.html
  !! Complex data => Coefficients are ordered by l and m contiguous
  !! Loops of the following kind are contiguous in memory with this indexing
  !! do l in [0,1,..., bw-1]
  !!   do m in [-(bw-1),...,bw-1]
  function LMc(l,m) result(ind)
    !! Result ind was previousely called index this caused errors in f2py
    !! it automatically canged the kind of the output to sp=4 instead of 8.
    integer(kind = dp),intent(in) :: l,m
    integer(kind = dp) :: ind
    ind = l*(l+1_dp)+m+1_dp
  end function LMc

  function LMc_slice(l) result(slice)
    integer(kind = dp),intent(in) :: l
    integer(kind = dp) :: slice(2)
    slice(1) = LMc(l,-l)
    slice(2) = LMc(l,l)
  end function LMc_slice
  !> ----------
  !! Index for the spherical harmonic coefficient Y_lm for real data.
  !! This is the same convention used by the software SHTNS
  !! https://nschaeff.bitbucket.io/shtns/index.html
  !! Real data => Coefficients are ordered by m and l contiguous
  !! Loops of the following kind are contiguous in memory with this indexing
  !! do m in [0,1,..., bw-1]
  !!   do l in [m,...,bw-1]
  !! ( ((((m*(2*lmax + 2 - m))>>1) + (l) )
  function MLr(m,l,bw) result(ind)
    !! Result ind was previousely called index this caused errors in f2py
    !! it automatically canged the kind of the output to sp=4 instead of 8.
    integer(kind = dp),intent(in) :: l,m,bw
    integer(kind = dp) :: ind
    ind = (m*(2_dp * bw - 1_dp - m))/2_dp + l + 1_dp
  end function MLr

  function MLr_slice(m,bw) result(slice)
    integer(kind = dp), intent(in) :: m,bw
    integer(kind = dp) :: slice(2)

    slice(1) = MLr(m,m,bw)
    slice(2) = MLr(m,bw-1_dp,bw)
  end function MLr_slice

  !> ----------
  !! Index for the spherical harmonic coefficient Y_lm for real data.
  !! Custom vesion that is ordered by m and l contiguous,
  !! while positive and negative m are interleaved.
  !! Loops of the following kind are contiguous in memory with this indexing
  !! do m in [0,1,-1,2,-2,..., bw-1,-(bw-1)]
  !!   do l in [|m|,...,bw-1]
  function MLc(m,l,bw) result(ind)
    !! Result ind was previousely called index this caused errors in f2py
    !! it automatically canged the kind of the output to sp=4 instead of 8.
    integer(kind = dp),intent(in) :: l,m,bw
    integer(kind = dp) :: ind
    ind = MAX(0, ABS(m)*(2_dp*bw - ABS(m)) - bw) + MERGE(0_dp,bw-abs(m),m>=0) + l + 1_dp
  end function MLc

  !> ------------
  !! Index for the spherical harmonic coefficient Y_lm for complex data.
  !! Custom vesion that is ordered by m and l contiguous,
  !! while positive and negative m are interleaved.
  !! Loops of the following kind are contiguous in memory with this indexing
  !! do m in [0,1,-1,2,-2,..., bw-1,-(bw-1)]
  !!   do l in [|m|,...,bw-1]
  function MLc_slice(m,bw) result(slice)
    integer(kind = dp),intent(in) :: m,bw
    integer(kind = dp) :: slice(2)
    
    slice(1) = MLc(m,0_dp,bw)+abs(m)
    slice(2) = slice(1) + bw-abs(m)-1_dp
  end function MLc_slice

end module utils

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
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: wsize
    wsize = bw*bw * int(2_dp + 3_dp*bw + bw*bw,dp)/3_dp
  end function size_wigner_d
  
  !> ------
  !! @brief Shape of all non-symmetric wigner-d matrix elements, $l<\\mathrm{bw}$.
  !!
  !! $\\left[2\\,\\mathrm{bw},\\frac{\\mathrm{bw}(2+3\\mathrm{bw}+\\mathrm{bw}^2)}{6}\\right]$
  function wigner_d_shape(bw) result(d_shape)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp) :: d_shape(2)
    d_shape = [2_dp*bw,(bw*(2_dp + 3_dp*bw + bw*bw))/6_dp]
  end function wigner_d_shape

  !> ------
  !! @brief Returns n uniformly sampled angles in $(0,\\pi)$
  !!
  !! $\frac{(2 i + 1) \pi}{2 n}$ for $i=0,\\ldots,n-1$
  !! Note: These are angles used in Chebyshev nodes of the first kind.
  function create_beta_samples(n) result(betas)
    integer(kind = dp) :: n,i
    real(kind = dp) :: betas(n),factor
    factor = pi / real(2_dp * n, dp)
    betas = [(real(2_dp * i - 1_dp, dp) * factor, i = 1, n)]
  end function create_beta_samples

  !> -------
  !! @brief Returns trig samples for Wigner-d computations via kostelec.
  !! 
  !! Does the following:
  !! ```fortran
  !!  betas = create_beta_samples(2*bw)
  !!  trig_samples(:,1) = cos(betas)                      ! <= Chebyshev nodes
  !!  trig_samples(:,2) = cos(betas/2_dp)*sin(betas/2_dp) ! <= needed for d^l_ml computation
  !!  trig_samples(:,3) = cos(betas/2_dp)**2              ! <= needed for d^l_ml computation
  !! ```
  function create_trig_samples(bw) result(trig_samples)
    integer(kind=dp) :: bw
    real(kind = dp) :: trig_samples(2*bw,3),betas(2*bw)
    betas = create_beta_samples(2*bw)
    trig_samples(:,1) = cos(betas)                      ! <= Chebyshev nodes
    trig_samples(:,2) = cos(betas/2_dp)*sin(betas/2_dp) ! <= needed for d^l_ml computation
    trig_samples(:,3) = cos(betas/2_dp)**2              ! <= needed for d^l_ml computation
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
  !! @brief Returns trig samples for Wigner-d computations via risbo.
  !! 
  !! Does the following:
  !! ```fortran
  !!  betas = create_beta_samples(2*bw)
  !!  trig_samples(:,1) = cos(betas/2)
  !!  trig_samples(:,2) = sin(betas/2)
  !! ```
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
    integer(kind = dp) :: n,i
    real(kind = dp) :: alpha(n),factor
    factor = 2._dp*pi / real(n, dp)
    alpha = [(real(i-1_dp, dp) * factor, i = 1, n)]
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
  !! The recursion can be started at $m\\_1=0$ by providing teh input array `dlml` with `dlml(1) = SQRT(0.5_dp)` (i.e.: $d^0\\_{0,0}(\beta) = \\sqrt{\frac{1}{2}}$ ). 
  !!
  !! 
  !! /// info|Recursion scheme\n
  !! ![dlml\\_l]{../images/dlml.svg#only-light}(width=400 align=left)\n
  !! ![dlml\\_l]{../images/dlml_white.svg#only-dark}(width=400 align=left)\n
  !! In the lefthand diagram shows the recursion scheme for `bw=6`. Diagonal arrows correspond to usage of the first recursion relation, i.e. $d^l\\_{m_1-1,l}\rightarrow d^{l+1}\\_{m_1,l}$ , and horizontal steps imply usage of the second recursion relation, i.e. $d^l\\_{m_1,l}\rightarrow d^{l+1}\\_{m_1,l}$.\n
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
    integer(kind = dp) :: i,l,swap_id,swap_id2,nc1,nc2

    !Normalization defining constants
    nc1 = MERGE(2_dp,1_dp,normalized)
    nc2 = nc1+1_dp
    
    swap_id = m+1_dp
    swap_id2 = Merge(m,swap_id,(m==0))
    do i=1,min(bw-m,swap_id2)
       !(l,m-1 -> l+1,m)
       l=m+(i-1_dp)-1_dp
       dlml(:,i) = dlml(:,i) * SQRT(real((2_dp*l+nc1)*(2_dp*l+nc2),kind=dp)/real((l+m)*(l+m+1_dp),kind=dp))*cos2
    end do
    do i=swap_id,bw-m-1_dp
       !(l -> l+1)
       l=m+(i-1_dp)
       dlml(:,i+1) = dlml(:,i) * SQRT(real((2_dp*l+nc1)*(2_dp*l+nc2),kind=dp)/real((l+m+1_dp)*(l-m+1_dp),kind=dp))*sincos
    end do
  end subroutine dlml_recursion_l_contiguous
  
  !> ----------
  !! @brief Computes all $d^l\\_{m\\_1,m\\_2}(\beta)$ for $m\\_2=l<\\mathrm{bw}$.
  !!
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
  function compute_all_dlml_l_contiguous(bw,sincos,cos2,normalized) result(dlml)
    !! bandwidth 0<=l<bw
    integer(kind=dp),intent(in) :: bw
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(SIZE(cos2,1),(bw*(bw+1))/2_dp),dlml_tmp(Size(cos2,1),bw)
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
  !! DO NOT use for perfomance critical computations!  
  !! This is a convenience function that computes the small wigner-d matrix elements $d^l_{m,l}$ for specific $l,m$.
  !! It does so using the l-contiguous recurrence relation `dlml_recursion_l_contuguous`.
  function compute_dlml(l,m,sincos,cos2,normalized) result(dlml)
    !! Indices $l,m$ at which to compute d^l_{m,l}
    integer(kind = dp),intent(in) :: l,m
    !! Arrays containing $sin(\frac{\beta}{2})*cos(\frac{\beta}{2})$ and $cos(\frac{\beta}{2})**2$
    real(kind = dp),intent(in) :: sincos(:),cos2(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    !! Array of output Wigner small d values 
    real(kind = dp) :: dlml(size(cos2,1)), dlml_tmp(size(cos2,1),l+1)
    integer(kind = dp) :: mm

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
  !!
  subroutine wig_l_recurrence_kostelec(workspace,cos,l,m1,m2,normalized)
    !
    integer(kind = dp), intent(in) :: l,m1,m2
    real(kind = dp),intent(in) :: cos(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp),intent(inout) :: workspace(:,:)
    integer(kind = dp) :: dl_id,dl_1id,o,i
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
  !!  All of these Wigners will have L2 norm = 1
  !!  
  !!  
  !!  let j = max(|m1|, |m2|)
  !!  
  !!  The functions generated will be
  !!  
  !!  d_{m1,m2}^j, d_{m1,m2}^{j+1}, ..., d_{m1,m2}^{bw-1}
  !!  
  !!  Each of these functions will be evaluated at the n = 2*bw-many
  !!  points
  !!  
  !!  pi*(2 * [0..n-1] + 1) / ( 2 * n )
  !!  
  !!  If beta(k) = pi*(2*k+1)/(2*n), then what's returned will be the
  !!  array
  !!  
  !!  d_{m1,m2}^j(beta(0)) ... d_{m1,m2}^{bw-1}(beta(0))
  !!  d_{m1,m2}^j(beta(1)) ... d_{m1,m2}^{bw-1}(beta(1))
  !!  d_{m1,m2}^j(beta(2)) ... d_{m1,m2}^{bw-1}(beta(2)) ...
  !!  d_{m1,m2}^j(beta(n-1)) ... d_{m1,m2}^{bw-1}(beta(n-1))
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
    integer(kind = dp),intent(in) :: m1,m2,bw
    real(kind = dp),intent(in) :: cos(:),dlml(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
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
       call wig_l_recurrence_kostelec(workspace,cos,l,m1,m2,normalized)
       o = MODULO(i+1_dp,3)+1_dp
       wigners_m1m2(:,i+1) = workspace(:,o)
    end do
  end function genwig_l2

  !> ---------
  !! @brief Same as genwig_l2 but transposed 
  !!
  !! This routine is faster than calling genwig_l2 and then transposing.
  function genwig_l2_trsp(m1,m2,bw,cos,dlml,normalized) result(wigners_m1m2)
    integer(kind = dp), intent(in) :: m1,m2,bw
    real(kind = dp), intent(in) :: cos(:),dlml(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp) :: wigners_m1m2((bw-max(abs(m1),abs(m2))),Size(cos,1))
    
    real(kind = dp) :: workspace(Size(cos,1),3)
    integer(kind = dp) :: i,l_start,l,n_samples,o

    
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
       call wig_l_recurrence_kostelec(workspace,cos,l,m1,m2,normalized)
       o = MODULO(i+1_dp,3)+1_dp
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
    integer(kind = dp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind = dp) :: dl(SIZE(betas,1),2*l+1,2*l+1)
    real(kind = dp) :: trig(SIZE(betas,1),3)
    real(kind = dp) :: workspace(SIZE(betas,1),3)
    real(kind=dp) :: dlml_workspace(SIZE(betas,1),l+1),sym_const_m,sym_const_n,sym_const_l
    integer(kind = dp) :: m,n,i,mid,nid,neg_mid,neg_nid,o,bw,nbeta

    dl = 0
    sym_const_l = (-1._dp)**l
    trig(:,1) = cos(betas)                      
    trig(:,2) = cos(betas/2_dp)*sin(betas/2_dp) 
    trig(:,3) = cos(betas/2_dp)**2
    nbeta = SIZE(betas,1)
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
             call wig_l_recurrence_kostelec(workspace,trig(:,1),i,m,n,normalized)
          end do
          
          o = MODULO(l-n+1_dp,3)+1_dp

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
    
    integer(kind = dp), intent(in) :: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
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
    
    integer(kind = dp),intent(in):: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp) :: wigners(bw*2_dp,(bw * (2_dp + 3_dp*bw + bw*bw))/6_dp)
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
    integer(kind = dp), intent(in)  :: l
    real(kind = dp) :: d_l(Size(d_l_1,1),(2*l+1)*(2*l+1)),temp(Size(d_l_1,1),(2*l)*(2*l)),inv_deg
    integer(kind = dp) :: deg,i,j,ij,ippj,ijpp,ippjpp,pdeg
    
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
    integer(kind = dp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical, intent(in) :: normalized
    real(kind = dp) :: dl_tmp(Size(betas,1),(2*l+1)*(2*l+1)), dl(SIZE(betas,1),2*l+1,2*l+1)
    real(kind=dp) :: cos_b(Size(betas,1)),sin_b(Size(betas,1)),ls_sqrt(2*l+1)
    integer(kind = dp) :: i,j
    ! Setup
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    do i=0_dp,2_dp*l
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl_tmp = 0._dp
    do i=0_dp,l
       dl_tmp(:,1:(2*i+1)*(2*i+1)) = wigner_mn_recurrence_risbo(dl_tmp(:,1:(2*i-1)*(2*i-1)), i, cos_b, sin_b, ls_sqrt)
    end do

    ! Save transposed array to output
    do i=1,2*l+1
       do j=1,2*l+1
          dl(:,j,i) = dl_tmp(:,(i-1_dp)*(2*l+1_dp)+j)
       end do
    end do

    if (normalized) then
       dl = dl * SQRT(real(2*l+1_dp,kind=dp)/2._dp)
    end if
  end function wigner_dl_risbo


  !>--------
  !! @brief Symmetry reduced risbo recurrence for $d^l\\_{m,n}$
  !!
  !! Risbo recurrence for symmetry reduce small wigner-d matrices.
  function wigner_recurrence_risbo_reduced(dl_1,l,cos_beta,sin_beta,sqrts,normalized) result(dl)
    real(kind = dp), intent(in) :: cos_beta(:),sin_beta(:),sqrts(:),dl_1(:,:)
    integer(kind = dp), intent(in)  :: l
    logical, intent(in) :: normalized
    real(kind = dp) :: dl(Size(cos_beta,1),((l+1)*(l+2))/2),dl_halve(Size(cos_beta,1),(l*(l+1))/2),inv_deg
    real(kind = dp) :: l_sym,lpp_sym,n_sym,norm_const,norm_const2
    integer(kind = dp) :: ll,m,n,mn,nbeta
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

       l_sym   = real((-1_dp)**l,dp)
       lpp_sym = real((-1_dp)**(l+1_dp),dp)
       nbeta = SIZE(cos_beta,1)
       !!!!!!!!!!!!!!!!!!
       ! Step l -> l+1/2!
       !!!!!!!!!!!!!!!!!!
       inv_deg = 1._dp/real(2_dp*l-1_dp,dp)
       ll=l-1_dp
       norm_const = MERGE(inv_deg*SQRT(real(2_dp*l,dp)/real(2_dp*l-1,dp)),inv_deg,normalized)
       !! 0 <= m <= n <= l-2
       do m=0,l-2_dp
          !!! m=n
          mn = tri(m,m,l)
          !print *,m,m,mn
          dl_halve(:,mn)= ( dl_1(:,mn)                  *real(l+m,dp)                *cos_beta &
                        &  -dl_1(:,tri(m,m+1_dp,l))     *SQRT(real((l+m)*(ll-m),dp)) *sin_beta &
                        &  -dl_1(:,tri(m,m+1_dp,l))     *SQRT(real((ll-m)*(l+m),dp)) *sin_beta &
                        &  +dl_1(:,tri(m+1_dp,m+1_dp,l))*real(ll-m,dp)               *cos_beta &
                        & )*norm_const
          !print *, dl_1(1,:)
          !print *, tri(m,m+1_dp,l),tri(m+1_dp,m+1_dp,l),dl_halve(1,mn)
          
          do n=m+1_dp,l-2_dp
             mn = tri(m,n,l)
             !print *,m,n,mn
             dl_halve(:,mn)= ( dl_1(:,mn)                  *SQRT(real((l+m)*(l+n),dp))  *cos_beta &
                           &  -dl_1(:,tri(m,n+1_dp,l))     *SQRT(real((l+m)*(ll-n),dp)) *sin_beta &
                           &  +dl_1(:,tri(m+1_dp,n,l))     *SQRT(real((ll-m)*(l+n),dp)) *sin_beta &
                           &  +dl_1(:,tri(m+1_dp,n+1_dp,l))*SQRT(real((ll-m)*(ll-n),dp))*cos_beta &
                           & )*norm_const
          end do
          !!! n=l-1
          mn = tri(m,l-1_dp,l)
          !print *,m,l-1_dp,mn
          dl_halve(:,mn)=( dl_1(:,mn)                   *SQRT(real((l+m)*(l+l-1_dp),dp))  *cos_beta &
                        &  +dl_1(:,tri(m+1_dp,l-1_dp,l))*SQRT(real((ll-m)*(l+l-1_dp),dp)) *sin_beta &
                        & )*norm_const

       end do
       !! m=n=l-1       
       mn = tri(l-1_dp,l-1_dp,l)
       !print *,l-1_dp,l-1_dp,mn
       dl_halve(:,mn)= dl_1(:,mn) * real((l+l-1_dp),dp)*cos_beta * norm_const


       
       !!!!!!!!!!!!!!!!!!!!
       ! Step l+1/2 -> l+1!
       !!!!!!!!!!!!!!!!!!!!
       inv_deg = 1._dp/real((2_dp*l),dp)
       norm_const = MERGE(inv_deg*SQRT(real(2_dp*l+1,dp)/real(2_dp*l,dp)),inv_deg,normalized)
       norm_const2 = MERGE(0.5_dp*SQRT(real(2_dp*l+1,dp)/real(2_dp*l,dp)),0.5_dp,normalized)
       !! m=n=0
       mn = tri(0_dp,0_dp,l+1)
       dl(:,mn)  = (  2._dp  *dl_halve(:,tri(0_dp,0_dp,l))         *cos_beta &
                 &  - lpp_sym*dl_halve(nbeta:1:-1,tri(0_dp,0_dp,l))*sin_beta &
                 &  + l_sym  *dl_halve(nbeta:1:-1,tri(0_dp,0_dp,l))*sin_beta &
                 & )*norm_const2
       !! m=0 1<=n<=l-1
       do n=1,l-1
          mn = tri(0_dp,n,l+1_dp)
          n_sym = real((-1_dp)**(-n),dp)

          dl(:,mn)  = ( l_sym*n_sym  *dl_halve(nbeta:1:-1,tri(0_dp,n-1_dp,l))*SQRT(real(l*(l+n),dp))*cos_beta &
                    &  -lpp_sym*n_sym*dl_halve(nbeta:1:-1,tri(0_dp,n,l))     *SQRT(real(l*(l-n),dp))*sin_beta &
                    &  +              dl_halve(:,tri(0_dp,n-1_dp,l))         *SQRT(real(l*(l+n),dp))*sin_beta &
                    &  +              dl_halve(:,tri(0_dp,n,l))              *SQRT(real(l*(l-n),dp))*cos_beta &
                    & )*norm_const
       end do
       !! m=0 n=l
       mn = tri(0_dp,l,l+1_dp)
       dl(:,mn)  = ( dl_halve(nbeta:1:-1,tri(0_dp,l-1_dp,l))*cos_beta &
                 &  +dl_halve(:,tri(0_dp,l-1_dp,l))         *sin_beta &
                 & )*norm_const2*SQRT(2._dp)
       !! 1<=m<=n<=l-1
       do m=1,l-1
          !!! m=n
          mn = tri(m,m,l+1_dp)
          dl(:,mn)= ( dl_halve(:,tri(m-1_dp,m-1_dp,l))*real(l+m,dp)               *cos_beta &
                  &  -dl_halve(:,tri(m-1_dp,m,l))     *SQRT(real((l+m)*(l-m),dp)) *sin_beta &
                  &  -dl_halve(:,tri(m-1_dp,m,l))     *SQRT(real((l-m)*(l+m),dp)) *sin_beta &
                  &  +dl_halve(:,tri(m,m,l))          *real(l-m,dp)               *cos_beta &
                  & )*norm_const
          do n=m+1,l-1
             mn = tri(m,n,l+1_dp)
             dl(:,mn)= ( dl_halve(:,tri(m-1_dp,n-1_dp,l))*SQRT(real((l+m)*(l+n),dp))*cos_beta &
                     &  -dl_halve(:,tri(m-1_dp,n,l))     *SQRT(real((l+m)*(l-n),dp))*sin_beta &
                     &  +dl_halve(:,tri(m,n-1_dp,l))     *SQRT(real((l-m)*(l+n),dp))*sin_beta &
                     &  +dl_halve(:,tri(m,n,l))          *SQRT(real((l-m)*(l-n),dp))*cos_beta &
                     & )*norm_const
          end do
          !!! n=l
          mn = tri(m,l,l+1_dp)
          dl(:,mn)= ( dl_halve(:,tri(m-1_dp,l-1_dp,l))*SQRT(real((l+m)*(l+l),dp))*cos_beta &
                  &  +dl_halve(:,tri(m,l-1_dp,l))     *SQRT(real((l-m)*(l+l),dp))*sin_beta &
                  & )*norm_const
       end do
       !! m=n=l
       mn = tri(l,l,l+1_dp)
       dl(:,mn)= dl_halve(:,tri(l-1_dp,l-1_dp,l)) * real((l+l),dp)*cos_beta * norm_const
    end if
  end function wigner_recurrence_risbo_reduced

  !> ---------
  !! @brief Computes symmetry reduced part of the small Wigner d matrix $d^l\\_{m,n}(\beta)$ for fixed $l$.
  !!
  !! This function uses the recurrence by Risbo et. al. from function `wigner_mn_recurrence_risbo_reduced`.
  !! It only computes the non symmetry equivalent part that is its values for  0<=m<=n<=l.
  function wigner_dl_risbo_reduced(l,betas,normalized) result(dl)
    integer(kind = dp),intent(in) :: l
    real(kind = dp),intent(in) :: betas(:)
    logical, intent(in) :: normalized
    real(kind = dp) :: dl(Size(betas,1),((l+1)*(l+2))/2)
    real(kind=dp) :: cos_b(Size(betas,1)),sin_b(Size(betas,1)),ls_sqrt(2*l+1)
    integer(kind = dp) :: i,ndl,ndlmm
    ! Setup
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    do i=0_dp,2_dp*l
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl=0
    do i=0_dp,l
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
    integer(kind = dp),intent(in) :: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp), intent(inout) :: wigners(:,:)
    real(kind = dp) :: dl(2_dp*bw,(bw*(bw+1_dp))/2)
    real(kind = dp) :: betas(2_dp*bw),cos_b(2_dp*bw),sin_b(2_dp*bw),ls_sqrt(2_dp*bw)
    integer(kind=dp) :: ndl,ndlmm,i,l,mn,mnl,m,n

    betas = create_beta_samples(2_dp*bw)
    cos_b = COS(betas/2._dp)
    sin_b = SIN(betas/2._dp)
    do i=0_dp,2_dp*bw
       ls_sqrt(i+1) = SQRT(real(i,kind=dp))
    end do

    dl=0
    do l=0_dp,bw-1_dp
       ndl   = ((l+1)*(l+2))/2
       ndlmm = (l*(l+1))/2
       dl(:,1:ndl) = wigner_recurrence_risbo_reduced(dl(:,1:ndlmm), l, cos_b, sin_b, ls_sqrt,normalized)
       do m=0_dp,l
          do n=m,l
             mn = triangular_to_flat_index(m,n,l+1_dp)
             mnl = mnl_to_flat_index(m,n,l,bw)
             !$OMP SIMD
             do i=1,2_dp*bw
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
    
    integer(kind = dp),intent(in):: bw
    logical,intent(in) :: normalized
    !f2py logical :: normalized = 1
    real(kind=dp) :: wigners(bw*2_dp,(bw * (2_dp + 3_dp*bw + bw*bw))/6_dp)

    call genwig_all_risbo_preallocated(bw,wigners,normalized)
  end function genwig_all_risbo


  !> ------
  !! @brief Convert symmetry reduced Wigner-d to full Wigner-d
  !!
  !! Converts symmetry reduced $d^l\\_{m,n}$, i.e. with $0\\leq m \\leq n \leq l$, to
  !! the full matrix with $-l \leq m,n \leq l$.
  function sym_reduced_to_full_wigner(sym_dl,l) result(dl)
        real(kind = dp), intent(in) :: sym_dl(:,:)
        integer(kind = dp), intent(in) :: l
        real(kind = dp) :: dl(SIZE(sym_dl,1),2*l+1,2*l+1)
        !real(kind = dp), intent(out) :: dl(SIZE(sym_dl,1),int(SQRT(9_dp/4_dp-2_dp*real(SIZE(sym_dl,2),dp))-3_dp/2_dp) ,int(SQRT(9_dp/4_dp - 2_dp*real(SIZE(sym_dl,2),dp))-3_dp/2_dp)) ! would get rid of the degree l as input but f2py cant handle this ...
        integer(kind = dp) :: m,n,mn,nbeta,mid,nid,neg_mid,neg_nid
        real(kind=dp) :: sym_const_m,sym_const_n,sym_const_l
        
        sym_const_l=(-1._dp)**l
        nbeta=SIZE(sym_dl,1)
        dl=0
        do m=0,l
           mid = l+1_dp+m
           neg_mid = l+1_dp-m
           sym_const_m = (-1._dp)**m

           do n=m,l
              nid = l+1_dp+n
              neg_nid = l+1_dp-n
              sym_const_n = (-1._dp)**n

              mn = triangular_to_flat_index(m,n,l+1_dp)
              
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
     integer(kind = dp) :: wigner_size
     integer(kind = dp) :: bw = 0
     integer(kind = dp) :: Lmax = 0
     integer(kind = dp) :: recurrence_type = 0
     real(kind = dp), allocatable :: wigner_d(:,:)
     real(kind = dp), allocatable :: wigner_d_trsp(:)
     real(kind = dp), allocatable :: wigner_dlml(:,:)
     real(kind = dp), allocatable :: trig_samples(:,:)
     real(kind = dp), allocatable :: trig_samples_risbo(:,:)
     real(kind = dp), allocatable :: sqrts_risbo(:)
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
  subroutine check_recurrence_type(self)
    class(so3ft),intent(inout) :: self
    if ((self%recurrence_type /= kostelec_recurrence) .AND. (self%recurrence_type /= risbo_recurrence)) then
        print *, "Error: unknown recurrence type", self%recurrence_type
        stop
     end if
  end subroutine check_recurrence_type
  subroutine init(self,bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    integer(kind = dp), intent(in) :: recurrence_type
    logical, intent(in) :: init_ffts,precompute_wigners
    integer(kind = sp), intent(in) :: fftw_flags
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
    allocate(self%wigner_dlml(2_dp*bw,triangular_size(bw)))
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
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    integer(kind = dp), intent(in) :: recurrence_type
    logical, intent(in) :: init_ffts
    logical, intent(in) :: precompute_wigners
    integer(kind = sp), intent(in) :: fftw_flags
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
    call self%destroy_fft(.True.)
    call self%destroy_fft(.False.)
    self%bw=0
    self%lmax=0
  end subroutine destroy

  ! inverse wigner complex
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
    bw2 = 2_dp*bw
    
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
    bw2 = 2_dp*bw
    
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
    integer(kind = dp) :: i,m1,m2,L
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
    bw2 = 2_dp*bw
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
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + (sym_const_m1*sym_const_m2*coeff(c_loc+2_dp))*dlmn(i)
       end do
       !$OMP END SIMD

       !! -m1,-m2 !!
       s_ids=order_to_ids(-m1,-m2,bw)
       !$OMP SIMD
       do i=1,bw2
          so3func(i,s_ids(1),s_ids(2)) = so3func(i,s_ids(1),s_ids(2)) + (sym_const_m1*sym_const_m2*coeff(c_loc+3_dp))*dlmn(i)
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
    o = MERGE(2_dp,0_dp,(m1/=m2))
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+o+2_dp))*dlmn(i)
    end do
    !$OMP END SIMD

    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m2*coeff(c_loc+o+3_dp))*dlmn(i)
    end do
    !$OMP END SIMD
        
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i,s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i,s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+6_dp))*dlmn(i)
    end do
    !$OMP END SIMD

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i,s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i,s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m2*coeff(c_loc + 7_dp))*dlmn(i)
    end do
    !$OMP END SIMD
  end subroutine inverse_wigner_loop_body_cmplx_risbo
  subroutine inverse_wigner_trf_cmplx_risbo(self,coeff,so3func,use_mp)
    !f2py threadsafe
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(in) :: coeff(:)
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    integer(kind = dp) :: m1,m2,l,mnid
    real(kind=dp) :: sym_const_m1,sym_const_l,dl(2*self%bw,(self%bw*(self%bw+1))/2)
       
    ! zeroing so3func needed since it will be populated by +=
    so3func = 0

    if (use_mp) then       
       ! non-fft part of the SO(3) fourier transform
       !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,dl,sym_const_l,l) SHARED(so3func,coeff)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_dp,((l+1_dp)*(l+2_dp))/2_dp
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
                mnid = triangular_to_flat_index(m1,m2,l+1_dp)
                call inverse_wigner_loop_body_cmplx_risbo(self,coeff,so3func,dl(:,mnid),l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
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
  ! forward wigner complex
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
    bw2 = 2_dp*bw
    
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
    bw2 = 2_dp*bw
    
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
    integer(kind = dp) :: i,m1,m2,L
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
    bw2 = 2_dp*bw
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
    coeff(c_loc+1_dp) = DOT_PRODUCT(dlmn,so3func(:,s_ids(1),s_ids(2)))
   

    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    if (.NOT. m1==m2) then
       !!  m2,m1  !!
       s_ids = order_to_ids(m2,m1,bw)
       coeff(c_loc + 2_dp) = DOT_PRODUCT(dlmn,so3func(:,s_ids(1),s_ids(2)))*sym_const_m1*sym_const_m2
       
       !! -m1,-m2 !!
       s_ids = order_to_ids(-m1,-m2,bw)       
       coeff(c_loc + 3_dp) = DOT_PRODUCT(dlmn,so3func(:,s_ids(1),s_ids(2)))*sym_const_m1*sym_const_m2
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
    o=MERGE(2_dp,0_dp,(m1/=m2))
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    coeff(c_loc + o + 2_dp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_m1*sym_const_l
    
    !! -m1,m2 !!
    s_ids = order_to_ids(-m1,m2,bw)
    coeff(c_loc + o + 3_dp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_m2*sym_const_l

    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    coeff(c_loc + 6_dp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2))) *sym_const_m1*sym_const_l
    
    !! -m2,m1 !!
    s_ids = order_to_ids(-m2,m1,bw)
    coeff(c_loc + 7_dp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_m2*sym_const_l

  end subroutine forward_wigner_loop_body_cmplx_risbo
  subroutine forward_wigner_trf_cmplx_risbo(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: m1,m2,l,mnid
    real(kind=dp) :: sym_const_m1,sym_const_l,dl(2*self%bw,(self%bw*(self%bw+1))/2),dl_tmp(2_dp*self%bw)

      
    if (use_mp) then
       !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,dl_tmp,dl,sym_const_l,l) SHARED(so3func,coeff)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_dp,((l+1_dp)*(l+2_dp))/2_dp
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
                mnid = triangular_to_flat_index(m1,m2,l+1_dp)
                dl_tmp = dl(:,mnid)*self%legendre_weights
                call forward_wigner_loop_body_cmplx_risbo(self,so3func,coeff,dl_tmp,l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
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
  ! inverse wigner real
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

    sym_const_m2 = (-1.0_dp)**m2
    bw = self%bw
    bw2 = 2_dp*bw

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

    sym_const_m2 = (-1.0_dp)**m2
    bw = self%bw
    bw2 = 2_dp*bw

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
    integer(kind = dp) :: i,m1,m2,L,s_ids(2),s_ids_sym(2),bw
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

    sym_const_m2 = (-1.0_dp)**m2
    bw = self%bw
    bw2 = 2_dp*bw
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
    o = MERGE(2_dp,0_dp,(m1/=m2))
    
    !! m1,-m2 !!
    s_ids=order_to_ids(m1,-m2,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+o+2_dp))*dlmn(i)
    end do
    !$OMP END SIMD
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i,s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i,s_ids(1),s_ids(2)) + (sym_const_l*sym_const_m1*coeff(c_loc+6_dp))*dlmn(i)
    end do
    !$OMP END SIMD       
  end subroutine inverse_wigner_loop_body_real_risbo
  subroutine inverse_wigner_trf_real_risbo(self,coeff,so3func,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: so3func(:,:,:)
    complex(kind = dp), intent(in) :: coeff(:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: j,i,m1,m2,l,mnid,s_ids(2),s_ids_sym(2),bw,bw2
    real(kind=dp) :: sym_const_m1,sym_const_l,dl(2*self%bw,(self%bw*(self%bw+1))/2)

    bw  = self%bw
    bw2  = 2_dp*bw

    if (use_mp) then
       ! non-fft part of the SO(3) fourier transform
       !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,s_ids,s_ids_sym,i,j,dl,sym_const_l,l) SHARED(so3func,coeff,bw,bw2)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_dp,((l+1_dp)*(l+2_dp))/2_dp
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
       do j=1,(self%lmax+1_dp)*bw2
          m2 = ((j-1_dp)/bw2)
          i = MOD(j,bw2)+1_dp
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
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
                mnid = triangular_to_flat_index(m1,m2,l+1_dp)
                call inverse_wigner_loop_body_real_risbo(self,coeff,so3func,dl(:,mnid),l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,self%lmax
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
          do i=1,bw2
             so3func(i,s_ids_sym(1),s_ids_sym(2)) = CONJG(so3func(i,s_ids(1),s_ids(2)))
          end do
       end do
    end if
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
  ! forward wigner real
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
    bw2 = 2_dp*bw

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
    bw2 = 2_dp*bw
    
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
    integer(kind = dp) :: i,m1,m2,L
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
    bw2 = 2_dp*bw
    c_slice = coeff_slice_lmn(l,m1,m2)
    c_loc = c_slice(1)

    !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
    !! m1,m2 !!
    coeff(c_loc) = DOT_PRODUCT(dlmn,so3func(:,m1+1,m2+1))    
    
    if (m1 ==0 .AND. m2 ==0) return    ! prevents m1=m2=0 from beeing evaluated twice

    !! -m2,-m1 !!
    s_ids = order_to_ids(m2,m1,bw)
    coeff(c_loc+1_dp) = DOT_PRODUCT(dlmn,CONJG(so3func(:,s_ids(1),s_ids(2))))



    !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
    !! uses the symmetries:                                          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
    !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
    !! They cause a constant sign swap by (-1)^(m1-m2)               !!
    !! Due to real input both transforms have already been computed  !!
    !! Use symmetry $$D^l_{m1,m2}(\alpha,\beta,\gamma) =(-1)^{m2-m1} (D^l_{-m1,-m2}(\alpha,\beta,\gamma))^*$$
    if (m1/=m2) then
       !!  m2,m1  !!
       coeff(c_loc + 2_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 1_dp))
       
       !! -m1,-m2 !!
       coeff(c_loc + 3_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc))
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
    o=MERGE(2_dp,0_dp,(m1/=m2))
    
    !! m1,-m2 !!
    s_ids = order_to_ids(m1,-m2,bw)
    coeff(c_loc + o + 2_dp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_l*sym_const_m1
   
    !! -m1,m2 !!
    !! Use real symmetry
    coeff(c_loc + o + 3_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + o + 2_dp)) 
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    coeff(c_loc + 6_dp) = DOT_PRODUCT(dlmn,so3func(bw2:1:-1,s_ids(1),s_ids(2)))*sym_const_l*sym_const_m1    
    
    !! -m2,m1 !!
    coeff(c_loc + 7_dp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 6_dp))
  end subroutine forward_wigner_loop_body_real_risbo
  subroutine forward_wigner_trf_real_risbo(self,so3func,coeff,use_mp)
    class(so3ft),intent(in),target :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: m1,m2,l,mnid
    real(kind=dp) :: sym_const_m1,sym_const_l,dl(2*self%bw,(self%bw*(self%bw+1))/2),dl_tmp(2_dp*self%bw)

   if (use_mp) then
      !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1,dl_tmp,dl,sym_const_l,l) SHARED(so3func,coeff)
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          !$OMP DO
          do mnid=1_dp,((l+1_dp)*(l+2_dp))/2_dp
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
                mnid = triangular_to_flat_index(m1,m2,l+1_dp)
                dl_tmp = dl(:,mnid)*self%legendre_weights
                call forward_wigner_loop_body_real_risbo(self,so3func,coeff,dl_tmp,l,m1,m2,sym_const_l,sym_const_m1)
             end do
          end do
       end do
    end if
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

    fft_array=0.0_dp
    call self%inverse_wigner_trf_cmplx(coeff,fft_array,use_mp)
    call dfftw_execute_dft(self%plan_c2c_forward,fft_array,so3func)
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
    complex(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    fft_array = 0.0_dp
    call dfftw_execute_dft(self%plan_c2c_backward,so3func,fft_array)
    fft_array = fft_array * (1._dp/real(2_dp*self%bw,kind=dp)**2) ! * 1/(2*bw) * 1/(2*bw)    
    call self%forward_wigner_trf_cmplx(fft_array,coeff,use_mp)
  end subroutine soft_
  subroutine soft(self,so3func,coeff,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    complex(kind = dp), intent(in) :: so3func(:,:,:)
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

    fft_array=0.0_dp
    call self%inverse_wigner_trf_real(coeff,fft_array,use_mp)
    fft_array = CONJG(fft_array) ! to correct for the fact that we have to compute the forward not the backward fft.
    call dfftw_execute_dft_c2r(self%plan_c2r_backward,fft_array,so3func)
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
    real(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    fft_array=0.0_dp
    call dfftw_execute_dft_r2c(self%plan_r2c_forward,so3func,fft_array)
    fft_array = fft_array * (1._dp/real(2_dp*self%bw,kind=dp)**2) ! * 1/(2*bw) * 1/(2*bw)
    fft_array = CONJG(fft_array) ! to correct for the fact that we have to compute the backward not the forward fft.
    call self%forward_wigner_trf_real(fft_array,coeff,use_mp)    
  end subroutine rsoft_
  subroutine rsoft(self,so3func,coeff,use_mp)
    class(so3ft),intent(inout) :: self
    complex(kind = dp), intent(inout) :: coeff(:)
    real(kind = dp), intent(in) :: so3func(:,:,:)
    logical,intent(in) :: use_mp
    
    if (.NOT. self%plans_allocated_r) then
       call self%init_fft(.TRUE.)
    end if
    call self%rsoft_(so3func,coeff,self%fft_c2r_in,use_mp)
  end subroutine rsoft
  subroutine soft_many(self,so3funcs,coeffs,use_mp)
    !f2py threadsafe
    class(so3ft),intent(inout) :: self
    complex(kind=dp),intent(in) :: so3funcs(:,:,:,:)
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
    real(kind=dp),intent(in) :: so3funcs(:,:,:,:)
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
    integer(kind = dp) :: i,n,m1,m2,m,l,lmax,bw,pm1_slice(2),nm1_slice(2),pm1_tmp(2),nm1_tmp(2)   
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
       wig_norm(i+1_dp) = 2._dp*pi*SQRT(2._dp/real(2_dp*i+1_dp,kind = dp))
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

    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
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
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
    
    
    !! -m1,m2 !!
    s_ids=order_to_ids(-m1,m2,bw)
    cc_lmn = wig_norm * f_lm(nm1) * CONJG(g_lm(pm2)) * sym_const_m1 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD

    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers
    
    !! m2,-m1 !!
    s_ids=order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm * f_lm(pm2) * CONJG(g_lm(nm1)) * sym_const_m2 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD

    !! -m2,m1 !!
    s_ids=order_to_ids(-m2,m1,bw)
    cc_lmn = wig_norm * f_lm(nm2) * CONJG(g_lm(pm1)) * sym_const_m1 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
  end subroutine inverse_wigner_loop_body_corr_cmplx_risbo
  subroutine inverse_wigner_trf_corr_cmplx_risbo(self,f_lm,g_lm,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_lm(:),g_lm(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: mnid,m1,m2,l,bw
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
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_dp*l+1_dp,kind = dp))
          !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1) SHARED(fft_array,f_lm,g_lm,dl,sym_const_l,l)
          !$OMP DO
          do mnid=1_dp,((l+1_dp)*(l+2_dp))/2_dp
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
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_dp*l+1_dp,kind = dp))
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_dp)
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
    call dfftw_execute_dft(self%plan_c2c_forward,fft_array,cc)
    cc = cc * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)   
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
    integer(kind = dp), intent(in) :: radial_limits(:)
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

    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
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
    integer(kind = dp) :: i,n,m,m1,m2,lmax,bw,pm1_slice(2),pm1_tmp(2),s_ids(2),s_ids_sym(2)
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
       wig_norm(i+1) = 2._dp*pi*SQRT(2._dp/real(2_dp*i+1_dp,kind = dp))
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
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
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
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
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


    sym_const_m2 = (-1.0_dp)**m2
    sym_const_m1m2 = sym_const_m1*sym_const_m2
    l_start = m2+1_dp
    bw = self%bw
    bw2 = 2_dp*bw
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
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
    
    if (m1 == m2) return               ! prevents duplicates due to swapping equal numbers

    !! m2,-m1 !!
    s_ids = order_to_ids(m2,-m1,bw)
    cc_lmn = wig_norm * f_ml(pm2) * g_ml(pm1) * sym_const_m1m2 * sym_const_l
    !$OMP SIMD
    do i=1,bw2
       so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) = so3func(bw2+1_dp-i, s_ids(1),s_ids(2)) + cc_lmn*dlmn(i)
    end do
    !$OMP END SIMD
  end subroutine inverse_wigner_loop_body_corr_real_risbo
  subroutine inverse_wigner_trf_corr_real_risbo(self,f_ml,g_ml,fft_array,use_mp)
    class(so3ft),intent(inout),target :: self
    complex(kind = dp),target,intent(in) :: f_ml(:),g_ml(:)
    complex(kind = dp), intent(inout) :: fft_array(:,:,:)
    logical, intent(in) :: use_mp
    integer(kind = dp) :: mnid,m1,m2,l,bw,bw2,i,j,s_ids(2),s_ids_sym(2)
    real(kind=dp) :: sym_const_m1,sym_const_l,wig_norm,dl(2*self%bw,(self%bw*(self%bw+1))/2)

    bw = self%bw
    bw2 = 2_dp*bw
    
    ! zero fft array
    ! Important since not all elements will be written to before doing the fft
    fft_array = 0._dp

    if (use_mp) then       
       ! non-fft part of the SO(3) fourier transform
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_dp*l+1_dp,kind = dp))
          !$OMP PARALLEL PRIVATE(mnid,m1,m2,sym_const_m1) SHARED(fft_array,f_ml,g_ml,dl,sym_const_l,l,wig_norm)
          !$OMP DO
          do mnid=1_dp,((l+1_dp)*(l+2_dp))/2_dp
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
       do j=1,(self%lmax+1_dp)*bw2
          m2 = ((j-1_dp)/bw2)
          i = MOD(j,bw2)+1_dp
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
          fft_array(i,s_ids_sym(1),s_ids_sym(2)) = CONJG(fft_array(i,s_ids(1),s_ids(2)))
       end do
       !$OMP END DO
       !$OMP END PARALLEL
    else
       ! non-fft part of the SO(3) fourier transform + assembly of cc_lmn = wig_norm * f_ml_part * g_ml_part * sym_const_m1 * sym_const_m2       
       do l=0,self%lmax
          dl(:,1:((l+1)*(l+2))/2) = wigner_recurrence_risbo_reduced(dl(:,1:(l*(l+1))/2),l,self%trig_samples_risbo(:,1),self%trig_samples_risbo(:,2),self%sqrts_risbo,.True.)
          sym_const_l = (-1._dp)**l
          wig_norm = 2._dp*pi*SQRT(2._dp/real(2_dp*l+1_dp,kind = dp))
          do m1=0, l
             sym_const_m1 = (-1.0)**m1
             do m2=m1, l
                mnid = triangular_to_flat_index(m1,m2,l+1_dp)
                call inverse_wigner_loop_body_corr_real_risbo(self,f_ml,g_ml,fft_array,dl(:,mnid),l,m1,m2,wig_norm,sym_const_l,sym_const_m1)
             end do
          end do
       end do

       ! fill remaining 2d real fft symmetry values using f_{0,m1}=f_{0,-m1}^*
       do m2=1,self%lmax
          s_ids = order_to_ids(0_dp,m2,bw)
          s_ids_sym = order_to_ids(0_dp,-m2,bw)
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
    call dfftw_execute_dft_c2r(self%plan_c2r_backward,fft_array,cc)
    cc = cc * (1/(2.0_dp*pi)) ! * 1/(2*bw) * (2*bw)/(2*pi)   
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
    integer(kind = dp), intent(in) :: radial_limits(:)
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
  function load_fftw_wisdom(file_path) result(error_code)
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
  end function load_fftw_wisdom
  function save_fftw_wisdom(file_path) result(error_code)
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
  end function save_fftw_wisdom
  
end module softclass

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
  use softclass, only: so3ft,so3ft_ptr
  use omp_lib
  implicit none
contains
  subroutine OMP_set_num_threads_(nthreads)
    integer(kind = sp), intent(in) :: nthreads
    call OMP_set_num_threads(nthreads)
  end subroutine OMP_set_num_threads_
  function OMP_get_max_threads_(f2py_bug) result(nthreads)
    !! This function does not need an input ...
    !! but f2py generates an error if it does not have one ...
    logical, intent(in) :: f2py_bug
    logical :: f2py_bug_dummy
    !f2py logical :: f2py_bug = 0
    integer(kind = sp) :: nthreads
    f2py_bug_dummy = f2py_bug ! only there to suppress warning during f2py compilation.
    nthreads = OMP_get_max_threads()
  end function OMP_get_max_threads_
  
  function py_init_soft(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags) result(self)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    integer(kind = dp), intent(in) :: recurrence_type
    logical, intent(in) :: init_ffts
    logical, intent(in) :: precompute_wigners
    integer(kind = sp), intent(in) :: fftw_flags
    type(so3ft_ptr) :: self_ptr
    integer(kind = dp) :: self
    ALLOCATE(self_ptr%p)
    self_ptr%p = so3ft(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
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
    !f2py threadsafe
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = dp) :: bw
    call int_to_soft_pointer(self_int,self_ptr,self)
    bw = self%bw
  end function py_get_bw
  function py_get_lmax(self_int) result(lmax)
    !f2py threadsafe
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = dp) :: lmax
    call int_to_soft_pointer(self_int,self_ptr,self)
    lmax = self%lmax
  end function py_get_lmax
  function py_wigners_are_precomputed(self_int) result(precomputed)
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    logical :: precomputed
    call int_to_soft_pointer(self_int,self_ptr,self)
    precomputed = (allocated(self%wigner_d) .AND. allocated(self%wigner_d_trsp))
  end function py_wigners_are_precomputed
  function py_get_recurrence_type(self_int) result(recurrence_type)
    !f2py threadsafe
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    integer(kind = dp) :: recurrence_type
    call int_to_soft_pointer(self_int,self_ptr,self)
    recurrence_type = self%recurrence_type
  end function py_get_recurrence_type
  subroutine py_set_lmax(self_int,lmax)
    !f2py threadsafe
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
  subroutine py_reset(self_int,bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
    integer(kind = dp), intent(in) :: bw
    integer(kind = dp), intent(in) :: lmax
    !f2py integer :: lmax = bw - 1
    logical, intent(in) :: init_ffts
    !f2py logical :: init_ffts = 0
    logical, intent(in) :: precompute_wigners
    !f2py logical :: precompute_wigners = 0
    integer(kind = dp), intent(in) :: recurrence_type
    !f2py integer :: recurrence_type = 0
    ! kistelec_recurrence=0, risbo_recurrence=1
    integer(kind = sp), intent(in) :: fftw_flags
    !f2py integer :: fftw_flags = 0
    ! FFTW_ESTIMATE=64, FFTW_MEASURE=0
    integer(kind = dp),intent(in) :: self_int
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%init(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
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
    complex(kind = dp), intent(in) :: so3func(:,:,:)
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
    real(kind = dp), intent(in) :: so3func(:,:,:)
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
    complex(kind = dp), intent(in) :: so3funcs(:,:,:,:)
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
    real(kind = dp), intent(in) :: so3funcs(:,:,:,:)
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
    integer(kind = dp),intent(in) :: self_int
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
    integer(kind = dp),intent(in) :: self_int
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
    integer(kind = dp), intent(in) :: radial_limits(:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%cross_correlation_ylm_cmplx_3d(f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
  end subroutine py_cross_correlation_ylm_cmplx_3d
  subroutine py_cross_correlation_ylm_real(self_int,f_lm,g_lm,cc,use_mp)
    !f2py threadsafe
    integer(kind = dp),intent(in) :: self_int
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
    integer(kind = dp), intent(in) :: radial_limits(:)
    logical,intent(in) :: use_mp
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%cross_correlation_ylm_real_3d(f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
  end subroutine py_cross_correlation_ylm_real_3d

  subroutine py_fft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = dp), intent(in) :: self_int
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%fft(f1,f2)
  end subroutine py_fft
  subroutine py_ifft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = dp), intent(in) :: self_int
    complex(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%ifft(f1,f2)
  end subroutine py_ifft
  subroutine py_rfft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = dp), intent(in) :: self_int
    real(kind = dp), intent(in) :: f1(:,:,:)
    complex(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%rfft(f1,f2)
  end subroutine py_rfft
  subroutine py_irfft(self_int,f1,f2)
    !f2py threadsafe
    integer(kind = dp), intent(in) :: self_int
    complex(kind = dp), intent(in) :: f1(:,:,:)
    real(kind = dp), intent(inout) :: f2(:,:,:)
    type(so3ft_ptr) :: self_ptr
    type(so3ft),pointer :: self
    call int_to_soft_pointer(self_int,self_ptr,self)
    call self%irfft(f1,f2)
  end subroutine py_irfft

end module py
