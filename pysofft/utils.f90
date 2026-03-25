!> ---------
!! @brief Module containing mostly indexing helper function.
module utils
  use, intrinsic :: iso_c_binding, only : c_int, c_int64_t, c_double, c_double_complex, c_bool,c_int32_t
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
    integer(kind=idp) :: m1,m2,bw,num_coeff
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
    integer(kind=idp) bw,num_coeff
    num_coeff = (4_idp*(bw*bw*bw)-bw)/3_idp
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
    integer(kind=idp) :: m1,m2,bw,bw2
    integer(kind=idp) :: ids(2)
    if (abs(m1)>=bw .OR. abs(m2)>=bw) then
       print *, "Out of bounds error (orders_to_ids): |m1|>bw or |m2|>bw"
       stop
    end if
    bw2=2*bw
    ids(1) = (1_idp-sign(1_idp,m1))/2_idp*bw2 + m1+1_idp
    ids(2) = (1_idp-sign(1_idp,m2))/2_idp*bw2 + m2+1_idp
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
    integer(kind=idp), intent(in) :: n1,n2,bw
    integer(kind=idp) :: m(2),slice(2),an1,an2,temp1,temp2,temp3
    logical :: swapped,n1_neg,n2_neg
    an1 = ABS(n1)
    an2 = ABS(n2)
    m = MERGE([an1,an2],[an2,an1],an1<an2)

    temp1 = 4*(bw*m(1)*(bw-m(1)+1))-2*(bw**2+m(1)**2)-bw + m(1)*(4*m(1)**2+2)/3
    slice(1) = MERGE(temp1,0_idp,m(1)>0)
    
    temp1 = bw+2*(2*bw-m(2))*(m(2)-1)
    temp2 = 4*((bw-m(1))+(m(2)-1 - m(1))*(2*bw - m(1) - m(2)))
    slice(1) = slice(1)+MERGE(MERGE(temp1,temp2,m(1)==0),0_idp,m(1)/=m(2))

    swapped = (m(1)/=ABS(n1))
    n1_neg = (n1<0)
    n2_neg = (n2<0)
  
    temp1 = MERGE(1_idp,MERGE(2_idp,3_idp,(n2_neg .AND.  n1/=0_idp ) .OR. n2==0_idp )&
         &, n1_neg .AND. (n2_neg .OR. n2==0_idp  ) )
    
    temp2 = MERGE(MERGE(MERGE(1_idp,3_idp,swapped),MERGE(7_idp,5_idp,swapped),n2_neg) &
         &,MERGE(MERGE(6_idp,4_idp,swapped),2_idp,n2_neg),n1_neg)
  
    temp3 = MERGE(temp1,temp2, (m(1)==0) .OR. (m(2)==0) .OR. (m(1)==m(2)) )
    slice(1) = slice(1) + MERGE(temp3*(bw-m(2)),0_idp,(m(1)/=n1) .OR. (m(2)/=n2))
    slice(2) = slice(1) + bw-m(2)
    slice(1) = slice(1) + 1_idp ! 1 indexing
  end function coeff_slice_mnl

  !> ----------
  !! @brief Returns the location of $f^l\\_{n\\_1,n\\_2}$ for mnl ordered Wigner coefficients.
  !!
  !! Returns: `coeff_slice_mnl(m1,m2,bw)(1)+l-max(abs(m1),abs(m2))`
  function coeff_location_mnl(m1,m2,l,bw) result(id)
    ! This function returns the index of the coefficient array corresponding to f_{m1,m2}^l.
    ! Note that this index has to lie within the slice given by coefLoc_so3.
    integer(kind = idp) :: m1,m2,bw,l,m,slice(2),id
    slice = coeff_slice_mnl(m1,m2,bw)
    m =  max(abs(m1),abs(m2))
    if (l<m .OR. l>bw) then
       print *, "Out of bounds error (coeff_location): l does not satisfy max(|m1|,|m2|) <= l < bw "
       stop
    end if
    id = slice(1)+l-m
  end function coeff_location_mnl
  
  !> ------
  !! @brief Return all valid (l,m,n) ids in mnl order.
  !!
  !! Returns an array containing all valid l,m,n coefficient combinations
  !! in the native order they are stored in for mnl ordered Wigner coefficients.
  function get_coeff_degrees(bw) result(lmn)
    !! Returns an array containing all valid l,m,n coefficient combinations
    !! in the native order they are stored in.
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp) :: lmn((4_idp*(bw*bw*bw)-bw)/3_idp,3)
    integer(kind = idp) :: m1,m2,cslice(2),ls(bw),l
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

  !> --------
  !! @brief  $l,n\\_1,n\\_2 \rightarrow$ slice into lmn ordered $f^l\\_{n\\_1,n\\_2}$ array.
  !!
  !! Returns the slice of coefficients $f\\_(m,n)^l$ for all possible symmetry variants of m,n.
  !! It is chosen such that in loops of the follwing type, memory access
  !! for the coefficients is continuouse without jumps.
  !! ```
  !! do l=0,bw-1
  !!    do m=0,l
  !!       do n=m,l
  !!          slice = coeff_slice_lmn(l,m,n)
  !!          f(slice(1))      !f_{m,n} 
  !!          if m==0 and n==0, cycle  
  !!          f(slice(1)+1)    !f_{-n,-m}
  !!          if m/=n
  !!             f(slice(1)+2) !f_{n,m}
  !!             f(slice(1)+3) !f_{-m,-n}
  !!          if m==0 or n==0, cycle
  !!          o=MERGE(2,0,m/=n)
  !!          f(slice(1)+2+o)  !f_{m,-n}
  !!          f(slice(1)+3+o)  !f_{-m,n}
  !!          if m==n, cycle
  !!          f(slice(1)+6)    !f_{n,-m}
  !!          f(slice(1)+7)    !f_{-n,m}
  !! ```
  function coeff_slice_lmn(l,n1,n2) result(slice)
    ! Insight is 
    ! 1 4 4 4 4    top to bottom is m, left to right is n.
    !   4 8 8 8    Numbers are number of function occurances in the above loop
    !     4 8 8
    !       4 8
    !         4

    integer(kind=idp), intent(in) :: n1,n2,l
    integer(kind=idp) :: slice(2)
    integer(kind=idp) :: previouse_l,previouse_m,previouse_n,length ,m,n
    logical :: m_not_0, m_g_1,n_g_m
    m = MIN(ABS(n1),ABS(n2))
    n = MAX(ABS(n1),ABS(n2))
    
    m_not_0 = (m/=0)
    m_g_1 = (m>1)
    n_g_m = (n > m)
    
    ! sum_0^(l-1) (2*l+1)^2
    previouse_l = int((l*(4_idp*l**2-1_idp))/3,idp) 
    
    previouse_m = MERGE(1_idp+l*4_idp + MERGE(4*(m-1)*(2*l+1-m),0_idp,m_g_1),0_idp,m_not_0)
    previouse_n = MERGE(1_idp+4_idp*(n-m-1_idp) + MERGE(3_idp+4_idp*(n-m-1_idp),0_idp,m_not_0),0_idp,n_g_m)
    length = MERGE(1_idp,MERGE(4_idp,8_idp,(m==n) .OR. (m==0)),n==0)
    
    slice(1) = previouse_l+previouse_m+previouse_n
    slice(2) = slice(1)+length
    slice(1) = slice(1) + 1_idp ! 1 indexing
  end function coeff_slice_lmn
  
  !> --------
  !! @brief  Location of $f^l\\_{n\\_1,n\\_2}$ in an lmn ordered coefficient array.
  function coeff_location_lmn(l,n1,n2) result(id)
        integer(kind=idp), intent(in) :: n1,n2,l
        integer(kind=idp) :: id,slice(2),m,n
        
        n=MAX(ABS(n1),ABS(n2))
        m=MIN(ABS(n1),ABS(n2))
        slice = coeff_slice_lmn(l,m,n)
        if (n1==m .AND. n2==n) then
           id = slice(1) 
        else if (n1==-n .AND. n2==-m) then
           id = slice(1) + 1_idp
        else if (n1==n .AND. n2==m) then
           id = slice(1) + 2_idp
        else if (n1==-m .AND. n2 == -n) then
           id = slice(1) + 3_idp
        else if (n1==m .AND. n2 == -n) then
           if (m/=n) then
              id = slice(1) + 4_idp
           else
              id = slice(1) + 2_idp
           end if
        else if (n1==-m .AND. n2 == n) then
           if (m/=n) then
              id = slice(1) + 5_idp
           else
              id = slice(1) + 3_idp
           end if
        else if (n1==n .AND. n2 == -m) then
           id = slice(1) + 6_idp
        else ! (n1==-n .AND. n2 == m) then
           id = slice(1) + 7_idp
        end if
  end function coeff_location_lmn

  !> ------
  !! @brief Return all valid (l,m,n) ids in lmn order.
  !!
  !! Returns an array containing all valid l,m,n coefficient combinations
  !! in the native order they are stored in for lmn ordered Wigner coefficients.
  function get_coeff_degrees_risbo(bw) result(lmn)
    !! Returns an array containing all valid l,m,n coefficient combinations
    !! in the native order they are stored in.
    integer(kind = idp), intent(in) :: bw
    integer(kind = idp) :: lmn((4_idp*(bw*bw*bw)-bw)/3_idp,3)
    integer(kind = idp) :: m1,m2,start,cslice(2),l,id,o
    do l = 0, bw-1
       do m1 = 0, l
          do m2 = m1, l
             cslice = coeff_slice_lmn(l,m1,m2)
             start = cslice(1)
             lmn(start,1) = l
             lmn(start,2) = m1
             lmn(start,3) = m2

             if ((m1==0) .and. (m2==0)) cycle
             id = start+1_idp
             lmn(id,1) = l
             lmn(id,2) = -m2
             lmn(id,3) = -m1

             if (m1/=m2) then
                id = start+2_idp
                lmn(id,1) = l
                lmn(id,2) = m2
                lmn(id,3) = m1

                id = start+3_idp
                lmn(id,1) = l
                lmn(id,2) = -m1
                lmn(id,3) = -m2
             end if

             if ((m1==0) .or. (m2==0)) cycle
             
             o=MERGE(2_idp,0_idp,m1/=m2)
             id = start+2_idp+o
             lmn(id,1) = l
             lmn(id,2) = m1
             lmn(id,3) = -m2
             
             id = start+3_idp+o
             lmn(id,1) = l
             lmn(id,2) = -m1
             lmn(id,3) = m2
             
             if (m1==m2) cycle
             
             id = start+6_idp
             lmn(id,1) = l
             lmn(id,2) = m2
             lmn(id,3) = -m1
             
             id = start+7_idp
             lmn(id,1) = l
             lmn(id,2) = -m2
             lmn(id,3) = m1
          end do
       end do
    end do
  end function get_coeff_degrees_risbo

  !> ------
  !! @brief `(2*bw,2*bw,2*bw)`
  function euler_shape(bw) result(eshape)
    integer(kind=idp), intent(in) :: bw
    integer(kind=idp) :: eshape(3)
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
    integer(kind = idp) :: bw,j,k
    real(kind = dp) :: weights(2*bw),tempsum,k_odd,j_odd,xi
    xi = pi/real(4_idp*bw,kind=dp)
    
    do k=0,2*bw-1
       k_odd = real(2_idp*k+1_idp,kind=dp)
       tempsum = 0
       do j=0,bw-1
          j_odd = real(2_idp*j+1_idp,kind=dp)
          tempsum = tempsum + (1._dp/j_odd)*SIN(k_odd*j_odd*xi)
       end do
       tempsum = tempsum*( (2._dp/real(bw,kind=dp)) * SIN(k_odd*xi) )
       weights(k+1_idp)=tempsum
    end do
  end function legendre_quadrature_weights

  !> -------
  !! @brief Zeroed Wigner coefficients for a given bandwidth `bw`.
  !!
  !! Returns a complex array of shape $\left(\frac{(4\mathrm{bw}^3-\mathrm{bw})}{3},\right)$.
  function get_empty_coeff(bw) result(coeff)
    integer(kind = idp) :: bw
    complex(kind = dp) :: coeff((4_idp*(bw*bw*bw)-bw)/3_idp)
    coeff = 0.0
  end function get_empty_coeff

  !> -------
  !! @brief n zeroed Wigner coefficients for a given bandwidth `bw`.
  !!
  !! Returns a complex array of shape $\left(\frac{(4\mathrm{bw}^3-\mathrm{bw})}{3},n\right)$.
  function get_empty_coeff_many(bw,n) result(coeff)
    integer(kind = idp) :: bw,n
    complex(kind = dp) :: coeff((4_idp*(bw*bw*bw)-bw)/3_idp,n)
    coeff = 0.0
  end function get_empty_coeff_many

  !> --------
  !! @brief Zeroed complex function over SO(3) for a given bandwidth `bw`.
  !!
  !! Returns a complex array of shape `(2*bw,2*bw,2*bw)`
  function get_empty_so3func_cmplx(bw) result(so3func)
    integer(kind = idp) :: bw
    complex(kind = dp) :: so3func(2*bw,2*bw,2*bw)
    so3func = 0.0
  end function get_empty_so3func_cmplx

  !> --------
  !! @brief n zeroed complex functions over SO(3) for a given bandwidth `bw`.
  !!
  !! Returns a complex array of shape `(2*bw,2*bw,2*bw,n)`
  function get_empty_so3func_cmplx_many(bw,n) result(so3func)
    integer(kind = idp) :: bw,n
    complex(kind = dp) :: so3func(2*bw,2*bw,2*bw,n)
    so3func = 0.0
  end function get_empty_so3func_cmplx_many

  !> --------
  !! @brief Zeroed real function over SO(3) for a given bandwidth `bw`.
  !!
  !! Returns a real array of shape `(2*bw,2*bw,2*bw)`
  function get_empty_so3func_real(bw) result(so3func)
    integer(kind = idp) ::  bw
    real(kind = dp) :: so3func(2*bw,2*bw,2*bw)
    so3func = 0.0
  end function get_empty_so3func_real
  
  !> --------
  !! @brief n zeroed real functions over SO(3) for a given bandwidth `bw`.
  !!
  !! Returns a real array of shape `(2*bw,2*bw,2*bw,n)`
  function get_empty_so3func_real_many(bw,n) result(so3func)
    integer(kind = idp) ::  bw,n
    real(kind = dp) :: so3func(2*bw,2*bw,2*bw,n)
    so3func = 0.0
  end function get_empty_so3func_real_many
  
  !> --------
  !! @brief Enforces real symmetry in mnl ordered coefficients
  !!
  !! Enforces the symmetry ${f^l\\_{m,n}}^* = f^l\\_{-m,-n} (-1)^{m+n}$ in
  !! m,n,l ordered coefficients
  subroutine enforce_real_sym_mnl(coeff,bw)
    complex(kind=dp) ,intent(inout) :: coeff(:)
    integer(kind=idp) ,intent(in) :: bw
    real(kind=dp) :: sym_const_m1,sym_const_m2
    integer(kind=idp) :: m1,m2,c_slice(2),c_slice_sym(2)

    c_slice = coeff_slice_mnl(0_idp,0_idp,bw)
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
  !! Enforces the symmetry ${f^l\\_{m,n}}^* = f^l\\_{-m,-n} (-1)^{m+n}$ in
  !! m,n,l ordered coefficients.
  !! 
  !! coeff_many is now a 2D array where the second index labels differenct coefficient arrays.
  subroutine enforce_real_sym_mnl_many(coeff_many,bw)
    complex(kind=dp) ,intent(inout) :: coeff_many(:,:)
    integer(kind=idp) ,intent(in) :: bw
    integer(kind=idp) :: n

    do n=1,SIZE(coeff_many,2)
       call enforce_real_sym_mnl(coeff_many(:,n),bw)
    end do
  end subroutine enforce_real_sym_mnl_many
  
  !> --------
  !! @brief Enforces real symmetry in lmn ordered coefficients
  !!
  !! Enforces the symmetry ${f^l\\_{m,n}}^{\*} = f^l\\_{-m,-n} (-1)^{m-n}$ in
  !! l,m,n ordered coefficients
  !!
  subroutine enforce_real_sym_lmn(coeff,bw)
    complex(kind=dp) ,intent(inout) :: coeff(:)
    integer(kind=idp) ,intent(in) :: bw
    real(kind=dp) :: sym_const_m1,sym_const_m2
    integer(kind=idp) :: m1,l,m2,c_slice(2),c_loc,o

    do l=0,bw-1_idp
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
                !coeff(c_loc+2_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 1_idp))
                coeff(c_loc+1_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc + 2_idp))
                !! -m1,-m2 !!
                coeff(c_loc+3_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc)) 
             else
                coeff(c_loc+1_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc))
             end if
             if (m1==0) cycle
             o=MERGE(2_idp,0_idp,(m1/=m2))

             !! m1,-m2 !!
             !! -m1,m2 !!
             coeff(c_loc + o +3_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc+o+2_idp))
             if (m1==m2) cycle
             !! m2,-m1 !!
             !! -m2,m1 !!
             coeff(c_loc + 7_idp) = (sym_const_m1*sym_const_m2)*CONJG(coeff(c_loc+6_idp))
          end do
       end do
    end do
  end subroutine enforce_real_sym_lmn

  !> --------
  !! @brief many version of enforce_real_sym_lmn
  !!
  !! Enforces the symmetry ${f^l\\_{m,n}}^{\*} = f^l\\_{-m,-n} (-1)^{m+n}$ in
  !! l,m,n ordered coefficients.
  !! 
  !! coeff_many is now a 2D array where the second index labels differenct coefficient arrays.
  subroutine enforce_real_sym_lmn_many(coeff_many,bw)
    complex(kind=dp) ,intent(inout) :: coeff_many(:,:)
    integer(kind=idp) ,intent(in) :: bw
    integer(kind=idp) :: n
    do n=1,SIZE(coeff_many,2)
       call enforce_real_sym_mnl(coeff_many(:,n),bw)
    end do
  end subroutine enforce_real_sym_lmn_many
  
  
  ! Indexing tricks section
  !< --------------
  !! @brief Total number of elements of a triangular index of bandwidth `bw`
  !!
  !! Consider the triangular index $0<=i<=j<bw$
  !! This function returns the total number of possible pairs $i,j$
  function triangular_size(bw) result(tri_size)
    integer(kind = idp),intent(in) :: bw
    integer(kind = idp) :: tri_size
    tri_size = (bw*(bw+1))/2_idp
  end function triangular_size

  !> -----------
  !! @brief Total number of elements of a triangular index of bandwidth `bw`
  !!
  !! consider the pyramid index $0<=i<bw$ and $-i<=j<=i$.
  !! This function returns the total number of possible pairs $i,j$.
  function pyramid_size(bw) result(pyr_size)
    integer(kind = idp),intent(in) ::bw
    integer(kind = idp) :: pyr_size
    pyr_size = bw**2
  end function pyramid_size

  !> ---------
  !! @brief Running index to pyramid index.
  !!
  !! Converts a running index k=0 to bw*(bw+1)/2-1 into
  !! a pyramid double index i,j with  0<=i<bw and -i<=j<=i
  !! This allows to reformulate pyramid loops as simple loops via
  !! ```
  !! do k=1, (N+1)**2
  !!    call flat_to_pyramid_index(i,j,k)
  !! ```
  !! is the same as
  !! ```
  !! do i=0,N
  !!   do j=-i,i
  !! ```
  subroutine flat_to_pyramid_index(i,j,k)
    integer(kind = idp), intent(inout) :: i,j
    integer(kind = idp), intent(in) :: k
    i = int(SQRT(real(k-1,dp)),kind = idp)
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
  !! do k=1,((N+1)*(N+2))/2_idp
  !!   triangular_to_flat(i,j,k,N)
  !!   print*,i,j    
  !! end do
  !!\endcode
  subroutine flat_to_triangular_index(i,j,k,N)
    integer(kind = idp), intent(inout) :: i,j
    integer(kind = idp), intent(in) :: k,N
    i = int(-SQRT(real((2_idp*N+3_idp)**2-8_idp*(k-1_idp),kind = dp))+real(2_idp*N+3_idp,kind=dp),kind= idp)/2_idp
    j =  (k-1_idp)-i*(2_idp*N+1_idp-i)/2_idp
    !i = (-int(SQRT(real((2_idp*N+3_idp)**2-8_idp*(k-1_idp),kind = dp)),kind = dp)+2_idp*N+3_idp)/2_idp
    !j =  (k-1_idp)-i*(2_idp*N+1_idp-i)/2_idp
  end subroutine flat_to_triangular_index
  !> ------------------
  !! @brief Converts triangular index (i,j) to 1d index.
  !!
  !! Consider the triangular index $0<=i<=j<\text{bw}$.
  !! 
  !! This function converts an index (i,j) to a one dimensional index.
  !! For fixed $i$ the differenct possible values of $j$ are contiguous in the one dimesional representation. 
  function triangular_to_flat_index(i,j,bw) result(id)
    integer(kind = idp), intent(in) :: i,j,bw 
    integer(kind = idp) :: id 
    id = (i*(2_idp*bw-i+1_idp))/2_idp + j-i + 1_idp
  end function triangular_to_flat_index
  !> ---------
  !! @brief d
  !!
  !! Consider the triangular index 0<=i<=j<bw 
  !! This function returns the slice of of a flattened array that corresponds
  !! to all valid i for a fixed j.
  !! This function allows to store a triangular array in an i contiguous way
  function triangular_to_flat_index_reversed(j,i) result(id)
    integer(kind = idp), intent(in) :: j,i
    integer(kind = idp) :: id
    id = (j*(j+1_idp))/2_idp + i +1_idp
  end function triangular_to_flat_index_reversed

  function lmn_to_flat_index(l,m,n) result(id)
    !! Considers all possible indexes 0<=m<=n<=l
    !! such that the following loop is cintiguous in memory
    !! do i=0,l
    !!    do j=0,m
    !!       do k=j,n
    integer(kind = idp),intent(in) :: l,m,n
    integer(kind = idp) :: id
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
    integer(kind = idp),intent(in) :: l,m,n,bw
    integer(kind = idp) :: id
    if (.NOT. (0<=m .AND. m<=n .AND. n<=l .AND. l<bw )) then
       print *, "Invalid arguments: l,m,n have to satisfy 0<=m<=n<=l<bw"
    end if
    id = (m*(3*bw**2+m**2+2-3*bw*(m-2)-3*m))/6+((n-m)*(2*bw+1-m-n))/2 + l-n+1
  end function mnl_to_flat_index
  function mnl_to_flat_l_slice(m,n,bw) result(slice)
    integer(kind = idp),intent(in) :: m,n,bw
    integer(kind = idp) :: slice(2)
    slice(1) = mnl_to_flat_index(m,n,n,bw)
    slice(2) = slice(1) + bw - 1_idp - n 
  end function mnl_to_flat_l_slice
  function mnl_to_flat_l_slice_padded(m,n,bw,pad_size) result(slice)
    integer(kind = idp),intent(in) :: m,n,bw,pad_size
    integer(kind = idp) :: slice(2)
    slice(1) = (mnl_to_flat_index(m,n,n,bw)-1_idp)*pad_size+1_idp
    slice(2) = slice(1) + pad_size*(bw - n) - 1_idp 
  end function mnl_to_flat_l_slice_padded

  !> ----------
  !! Computes number of complex spherical harmonic coefficients for
  !! a given order bandwidth bw.
  function n_LMc(bw) result(n)
    integer(kind=idp), intent(in) :: bw
    integer(kind=idp) :: n
    n=bw*bw
  end function n_LMc

  !> ----------
  !! Computes number of real spherical harmonic coefficients for
  !! a given order bandwidth bw.
  function n_LMr(bw) result(n)
    integer(kind=idp), intent(in) :: bw
    integer(kind=idp) :: n
    n= (bw*(bw+1_idp))/2_idp
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
    integer(kind = idp),intent(in) :: l,m
    integer(kind = idp) :: ind
    ind = l*(l+1_idp)+m+1_idp
  end function LMc

  function LMc_slice(l) result(slice)
    integer(kind = idp),intent(in) :: l
    integer(kind = idp) :: slice(2)
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
    integer(kind = idp),intent(in) :: l,m,bw
    integer(kind = idp) :: ind
    ind = (m*(2_idp * bw - 1_idp - m))/2_idp + l + 1_idp
  end function MLr

  function MLr_slice(m,bw) result(slice)
    integer(kind = idp), intent(in) :: m,bw
    integer(kind = idp) :: slice(2)

    slice(1) = MLr(m,m,bw)
    slice(2) = MLr(m,bw-1_idp,bw)
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
    integer(kind = idp),intent(in) :: l,m,bw
    integer(kind = idp) :: ind
    ind = MAX(0, ABS(m)*(2_idp*bw - ABS(m)) - bw) + MERGE(0_idp,bw-abs(m),m>=0) + l + 1_idp
  end function MLc

  !> ------------
  !! Index for the spherical harmonic coefficient Y_lm for complex data.
  !! Custom vesion that is ordered by m and l contiguous,
  !! while positive and negative m are interleaved.
  !! Loops of the following kind are contiguous in memory with this indexing
  !! do m in [0,1,-1,2,-2,..., bw-1,-(bw-1)]
  !!   do l in [|m|,...,bw-1]
  function MLc_slice(m,bw) result(slice)
    integer(kind = idp),intent(in) :: m,bw
    integer(kind = idp) :: slice(2)
    
    slice(1) = MLc(m,0_idp,bw)+abs(m)
    slice(2) = slice(1) + bw-abs(m)-1_idp
  end function MLc_slice

end module utils
