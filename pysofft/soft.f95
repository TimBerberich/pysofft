module precission
implicit none
integer(kind=8), parameter :: dp = 8
end module precission

module soft
use precission
use, intrinsic :: iso_c_binding
implicit none
include 'fftw3.f03'

contains
    function n_coeffs(m1,m2,bw) result(n_coeff)
        !  For orders m1, m2, and bandwidth bw, returns how many coefficients
        !  will be computed.
        !  
        !  let m = Max( |m1|, |m2| )
        !  
        !  The number of coefficients is = bw - m
        integer(kind=8) :: m1,m2,bw,n_coeff
        if (abs(m1)>=bw .OR. abs(m2)>=bw) then
            print *, "Out of bounds error (n_coeffs): |m1|>bw or |m2|>bw"
            stop
        end if
        n_coeff = bw - max(abs(m1),abs(m2))
    end function n_coeffs
  
    function total_n_coeffs(bw) result(n_coeff)
        !  For bandwidth bw, returns the TOTAL NUMBER of coefficients
        !  that will be computed when doing a full SO3 forward transform.
        !
        !  total number = 1/3 * bw * (4 * bw^2 - 1 ) 
        !
        !  note that for integer bw, 4 * bw^3 - bw ,is always a multiple of 3,
        !  so integer division won't be messy.   
        integer(kind=8) bw,n_coeff
        n_coeff = (4_dp*(bw*bw*bw)-bw)/3_dp
    end function total_n_coeffs

    function sample_slice(m1,m2,bw) result(slice)
        ! For order m1, m2 and bandwidth bw, returns the slice in the 3-D,
        ! fft'd array of the data necessary for the Wigner-(m1,m2)
        ! transform. I.e. The slice of the fft'd data
        ! needed to compute f_{m1,m2}^l for all legal l. 
        !
        ! Note that the array is of size 2*bw x 2*bw x 2*bw, so I
        ! need to multiply by that 2*bw (I'm sampling at twice
        ! the bandwidth)
        integer(kind=dp) :: m1,m2,bw,bw2
        integer(kind=dp) :: slice(2)
        if (abs(m1)>=bw .OR. abs(m2)>=bw) then
            print *, "Out of bounds error (sample_slice): |m1|>bw or |m2|>bw"
            stop
        end if
        bw2=2*bw
        slice(1) = bw2*(bw2*( (1_dp-sign(1_dp,m1))/2*bw2 + m1 + (1_dp-sign(1_dp,m2))/2_dp ) + m2) 
        slice(1) = slice(1) + 1_dp ! Fortran indexing starting at 1 
        slice(2)=slice(1)+bw2
    end function sample_slice
    
    function coeff_slice(m1,m2,bw) result(slice)
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
        start = total_n_coeffs(bw)*(1_dp-sign(1_dp,m1))/2_dp &
                &+ sqrbw*(m1+(1_dp-sign(1_dp,m2))/2_dp) &
                &- (m1-sign(1_dp,m2))*m1*(2_dp*m1 - sign(1_dp,m2))/6_dp

        start = start + lsum
        start = start + 1_dp ! +1 due to fortran indexing starting at 1   
        slice(1) = start
        slice(2) = start + size
    end function coeff_slice

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

    function wigLen_so3(m1,m2,bw) result(wigLen)
        ! returns the number of small wigner d values d_{m1,m2}^l(\beta)
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
    
    function get_wigner_matrix(m1,m2,bw,wigners) result(wig_mat) !! still a mock
        ! returns the number of small wigner d values d_{m1,m2}^l(\beta)
        ! for a given choice of m1,m2 and bw.
        ! Since beta is sampled at twice the bandwith the number of values is
        !
        ! number of values = (number of possible l) * 2*bw
        integer(kind = dp) :: m1,m2,m,bw,nl,nbeta, wstart,wstop
        real(kind = dp), intent(in), target :: wigners(:)
        real(kind = dp), pointer :: wig_mat(:,:)
        
        nbeta = bw*2
        m = max(abs(m1),abs(m2))
        nl = bw-m
        if (m>=bw) then
            print *, "Out of bounds error (wigLen_so3): |m1|>bw or |m2|>bw"
            stop
        end if

        wstart = 1
        wstop = wstart + nbeta*nl
        wig_mat(1:nbeta,1:nl) => wigners(wstart:wstop)
    end function get_wigner_matrix


    function get_coeff_part(m1,m2,bw,coeff) result(coeff_part)
        integer(kind = dp) :: m1,m2,bw,c_slice(2)
        complex(kind = dp), pointer :: coeff_part(:)
        complex(kind = dp), target :: coeff(:)
        c_slice = coeff_slice(m1,m2,bw)
        coeff_part => coeff(c_slice(1):c_slice(2))
    end function get_coeff_part

    subroutine inverse_so3_naive_fft_pc(bw,coeff,wigners,so3func,fft_plans)
        integer(kind = dp), intent(in) :: bw
        real(kind = dp), intent(in), target :: wigners(:)
        complex(kind = dp), intent(in) :: coeff(:)
        complex(kind = dp), intent(inout), target :: so3func(2*bw,2*bw,2*bw)

        real(kind = dp), pointer :: wig_mat(:,:)
        complex(kind = dp), pointer :: so3func_flat(:)
        integer(kind = dp) :: i,m1,m2,m,L,sym_const,s_slice(2),c_slice(2),sym_array(bw),nls
        integer(kind = dp) :: wig_slice_shape(4)
        integer(kind=8) , intent (in) :: fft_plans(:)

        ! 1d pointer int o 3d array representing a function f(\alpha,\beta,\gamma) on SO(3)
        so3func_flat(1:(2*bw)**3) => so3func

        ! initiallizing some constants
        L = bw-1
        do i=0,L
            sym_array(i+1) = (-1)**i 
        end do

        ! non-fft part of the SO(3) fourier transform        
        do m1=0, L
            do m2=m1, L
                wig_mat = get_wigner_matrix(m1,m2,bw,wigners)

                !! normal branch for m1<=m2 and sgn(m1)==sgn(m2)                  !!
                !! use of symmetries does not cause a change in d_{m1,m2}^l(beta) !!
                !! m1,m2 !!
                s_slice = sample_slice(m1,m2,bw) 
                so3func_flat(s_slice(1):s_slice(2)) = matmul(wig_mat,get_coeff_part(m1,m2,bw,coeff))
                
                if (m1 ==0 .AND. m2 ==0) cycle    ! prevents m1=m2=0 from beeing evaluated twice
                !! -m2,-m1 !!
                s_slice = sample_slice(-m2,-m1,bw) 
                so3func_flat(s_slice(1):s_slice(2)) = matmul(wig_mat,get_coeff_part(-m2,-m1,bw,coeff))
                
                !! branch for m1>m2 and sgn(m1)==sgn(m2)                         !!
                !! uses the symmetries:                                          !!
                !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{m2,m1}^l(\beta)          !!
                !! d_{m1,m2}^l(\beta) = (-1)^(m1-m2)*d_{-m1,-m2}^l(\beta)        !!
                !! They cause a constant sign swap by (-1)^(m1-m2)               !!
                if (.NOT. m1==m2) then
                    sym_const = (-1)**(m1-m2)
                    !!  m2,m1  !!
                    s_slice = sample_slice(m2,m1,bw) 
                    so3func_flat(s_slice(1):s_slice(2)) = sym_const*matmul(wig_mat,get_coeff_part(m2,m1,bw,coeff))
                    !! -m1,-m2 !!
                    s_slice = sample_slice(-m1,-m2,bw) 
                    so3func_flat(s_slice(1):s_slice(2)) = sym_const*matmul(wig_mat,get_coeff_part(-m1,-m2,bw,coeff))
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
                if (m1==0 .or. m2==0) cycle       ! prevents sign swaps on 0 ids which are already covered
                sym_const = (-1)**m1
                m = max(abs(m1),abs(m2))
                !! m1,-m2 !!
                s_slice = sample_slice(m1,-m2,bw)
                so3func_flat(s_slice(2):s_slice(1):-1) = sym_const*matmul(wig_mat,sym_array(m:)*get_coeff_part(m1,-m2,bw,coeff))
                 
                if (m1 == m2) cycle               ! prevents duplicates due to swapping equal numbers
                sym_const = (-1)**m2
                !! -m1,m2 !!
                s_slice = sample_slice(-m1,m2,bw) 
                so3func_flat(s_slice(2):s_slice(1):-1) = sym_const*matmul(wig_mat,sym_array(m:)*get_coeff_part(-m1,m2,bw,coeff))
                !! -m2,m1 !!
                s_slice = sample_slice(-m2,m1,bw) 
                so3func_flat(s_slice(2):s_slice(1):-1) = sym_const*matmul(wig_mat,sym_array(m:)*get_coeff_part(-m2,m1,bw,coeff))   
                !! m2,-m1 !!
                s_slice = sample_slice(m2,-m1,bw)
                so3func_flat(s_slice(2):s_slice(1):-1) = sym_const*matmul(wig_mat,sym_array(m:)*get_coeff_part(m2,-m1,bw,coeff))      
            end do
        end do


    end subroutine inverse_so3_naive_fft_pc
    

    !subroutine inverse_so3_naive_fft_pc_real(bw,coeff,wigners)
    !end subroutine inverse_so3_naive_fft_pc_real
end module soft 
