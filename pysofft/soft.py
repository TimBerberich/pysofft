import numpy as np
import time
from numba import njit, objmode, types,int64,float64,complex128
from pysofft.make_wiegner import SinEvalPts,SinEvalPts2,CosEvalPts,CosEvalPts2,n_wigners
from pysofft.make_wiegner import genWigAll,genWigAllTrans,get_euler_angles,wig_little_d_grid,wig_little_d_lookup
from pysofft.wignerTransform import wigNaiveAnalysis_fftw,wigNaiveSynthesis_fftw
from pysofft.wignerTransform import wigNaiveAnalysis_fftwX,wigNaiveSynthesis_fftwX
from pysofft.wignerTransform import wigNaiveAnalysis_fftwY,wigNaiveSynthesis_fftwY
from pysofft.wignerWeights import makeweights2
from pysofft.rotate import rotate_coeff_multi,rotate_coeff


complex_3d = types.complex128[:,:,:]
pi = np.pi

##################################################################
#  At bandwidth bw, what are the TOTAL NUMBER of coefficients
#  that will be computed when doing a full SO3 forward transform ?
#
#  total number = 1/3 * bw * (4 * bw^2 - 1 )

@njit(int64(int64),cache=True)
def totalCoeffs_so3(bw):
    # note that 4 * bw^3 - bw is always a multiple of 3,
    # so integer division won't be messy */    
    return (4*(bw**3) - bw)/3


###########################################################
#  At orders m1, m2, and bandwidth bw, how many coefficients
#  will be computed ?
#  
#  let m = Max( |m1|, |m2| )
#  
#  The number of coefficients is = bw - m
@njit(int64(int64,int64,int64),cache=True)
def howMany_so3(m1,m2,bw):
    return bw - max( abs(m1),abs(m2) )


###############################################################
# At order m1, m2, bandwidth bw, where in the 3-D,
# fft'd array does the data necessary for the Wigner-(m1,m2)
# transform start? I.e. I need the location of the fft'd data
# in order to compute f_{m1,m2}^l for all legal l. With respect
# to the beginning of that block of memory, where does this
# occur?
#
# Note that the array is of size 2bw x 2bw x 2bw, so I
# need to multiply by that 2*bw (I'm sampling at twice
# the bandwidth)
@njit(int64[:](int64,int64,int64),cache=True)
def sampLoc_so3(m1,m2,bw):    
    if m1 >= 0:
        if m2 >= 0:
            start = (2*bw)*(2*bw*m1 + m2)
        else:
            start = (2*bw)*(2*bw*m1 + 2*bw + m2)
    else:
        if m2 >= 0:
            start = (2*bw)*(2*bw*(2*bw + m1) + m2)
        else:
            start = (2*bw)*(2*bw*(2*bw + m1) + 2*bw + m2)
    stop = start + 2*bw
    return np.array((start,stop))


#############################################################
# So the above function tells me where in the array
# I should start taking the samples from, in order
# to do the order (m1,m2)-Wigner-transform on the
# just fft'd data.
#
# This function, coefLoc, tells where in the coefficient
# array I should start placing this evaluated inner-products
# of this bandwidth bw transform.
@njit(int64[:](int64,int64,int64),cache=True)
def coefLoc_so3(m1,m2,bw):
    m = max( abs(m1),abs(m2) )
    size = bw-m
    if m1 >= 0:
        if m2 >= 0:
            start = (bw**2)*m1
            start -= (m1 - 1)*m1*(2*m1 - 1)/6
            for k in range(m2):
                start += howMany_so3(m1,k,bw)
        else:
            start = (bw**2)*(m1+1) 
            start -= (m1 + 1)*m1*(2*m1 + 1)/6 
            for k in range(m2,0):
                start -= howMany_so3(m1,k,bw)
    else:
        if m2 >= 0:
            start = totalCoeffs_so3(bw) 
            start -= (bw**2)*(-m1) 
            start += (-m1 + 1)*(-m1)*(2*(-m1) + 1)/6 
            for k in range(m2):
                start += howMany_so3(m1,k,bw)
        else:
            start = totalCoeffs_so3(bw) 
            start -= (bw**2)*(-m1 - 1)
            start += (-m1 - 1)*(-m1)*(2*(-m1) - 1)/6            
            for k in range(m2,0):
                start -= howMany_so3(m1,k,bw)
    stop = start + size
    #print([start,stop])
    return np.array((start,stop),dtype = np.int64)


#################################################################
#/*
#  so3CoefLoc: Where, in the coefficient array of the
#              SO(3) signal of bandwidth bw does the
#	      coefficient f_{MM'}^l reside?
#
#  Unfortunately, the indexing is a little hairy.
#
#                 bw: bandwidth
#		 l : degree
#		 m, mp : orders
#
#  NOTE: Here's the convention I'm observing:
#
#        D_{MM'}^l(alpha, beta, gamma) =
#
#	exp(-I M alpha) * d_{MM'}^l(beta) * exp(-I M gamma)
#
#	where m == M, mp == M'
#
#*/
@njit(int64(int64,int64,int64,int64))
def so3CoefLoc( m, mp, l, bw):
    k=0
    tmpA = 0
    tmpB =0
    
    if m>=0:
        if mp>=0:
            tmpA = ( 6*bw*bw + 3*m - 2*m*m - 1)
            tmpA *= m 
            tmpA /= 6 	  
            tmpB = 0 
            for k  in range(mp):
                tmpB += bw - max(m,k) #HPP(m, k, bw) #define HPP(M, K, BW) ((BW) - MAX((M), (K))) 
            tmpA += tmpB 
            tmpA += (l - max(m,mp))	
        else: # /* mp < 0 */	
            tmpA = -(m*(1+m)*(1+2*m))/6 
            tmpA += bw*bw*(1+m)
            tmpB = 0;
            for k in range(mp,0):# ( k = mp; k < 0 ; k ++ )
                tmpB += bw - max(m,-k) #HPM(m, k, bw) #define HPM(M, K, BW) ((BW) - MAX((M), (-K)))
            tmpA -= tmpB 
            tmpA += (l - max(m,-mp))
    else: # /* so m < 0 */
        if mp>=0 :
            tmpA = (4*bw*bw*bw-bw)/3 
            tmpB = (-m)*(1-m)*(1-2*m) 
            tmpB /= -6 
            tmpB += bw*bw*(-m) 
            tmpA -= tmpB 	  
            tmpB = 0 
            for k in range(mp): #( k = 0 ; k < mp ; k ++ )
                tmpB += bw - max(-m,k) #HMP(m, k, bw) #define HMP(M, K, BW) ((BW) - MAX((-M), (K)))
            tmpA += tmpB 
            tmpA += (l - max(-m,mp))    
        else: # /* mp < 0 */	
            tmpA = (4*bw*bw*bw-bw)/3 
            tmpB = ((-m)-1)*(6*bw*bw-m-2*(m*m))
            tmpB /= 6
            tmpA -= tmpB 
            tmpB = 0
            for k in range(mp,0): #( k = mp; k < 0 ; k ++ )
                tmpB += bw - max(-m,-k) #HMM(m, k, bw) #define HMM(M, K, BW) ((BW) - MAX((-M), (-K)))
            tmpA -= tmpB 
            tmpA += (l - max(-m,-mp))
    return tmpA


def coef_id_to_lnk(bw):
    _nk_range = np.roll(np.arange(-bw+1,bw),bw)
    id_to_lnk = np.zeros((totalCoeffs_so3(bw),3),dtype=np.int32)
    pos = 0
    for n in _nk_range:
        for k in _nk_range:
            for l in range(max(abs(n),abs(k)),bw):
                id_to_lnk[pos] = (l,n,k)
                pos += 1
    return id_to_lnk
    
@njit(int64(int64,int64,int64),cache=True)
def wigLen_so3(m1,m2,bw):
    m = max( abs(m1),abs(m2) )
    size = bw-m
    return (bw-m)*2*bw

def d_to_coeff_info(bw):
    '''
    helper function that returns index connections between the compuded small wigner d matrices for 0<=l,n,k<bw and the full small wigner d matrices for  0<=l<bw and  -bw+1<=l,n,k<bw in the native format that is ( n,k in 0,...,L,-L,...,1 and l in 0,...,L) where L = bw-1.

    Returns:
        np.ndarray: np.int32,((4*(bw**3) - bw)/3,3) Array of indices that converts the small wigner d array generated by genWigAll into the coefficient format for l,n,k.
    (n,0) start or stop index in genWigAll output of the n'th coefficient value, (n,1) +1 or -1 beta step of the n's coefficient, (n,2) sign factor +1 or -1. There are always 2*bw beta values for each coefficent.
        dict: Keys are of the form (l,n,k) and the corresponding values are entries of the above np.ndarray
    '''
    ids = np.zeros((totalCoeffs_so3(bw),3),dtype = np.int32)
    lnk_to_ids = wig_little_d_lookup(bw)
    lookup = {}
    start = time.time()
    _nk_range = np.roll(np.arange(-bw+1,bw),bw)
    pos = 0
    for n in _nk_range:
        for k in _nk_range:
            for l in range(max(abs(n),abs(k)),bw):                                
                sign_arg = 0
                absn=abs(n)
                absk=abs(k)
                swap_nk = k!=0 and (n==0 or absn>absk)
                if swap_nk:
                    # use d^l_nk(beta) = (-1)^{n-k} d^l_{kn}(beta)
                    nn = k
                    kk = n
                    absnn=absk
                    abskk=absn
                    sign_arg = n-k
                else:
                    nn=n
                    kk=k
                    absnn=absn
                    abskk=absk
                nneg = nn<0
                kneg = kk<0
                reverse_angles = (nneg ^ kneg) 
                if reverse_angles:
                    beta_step = -1
                else:
                    beta_step = 1
                
                sl = lnk_to_ids[(l,absnn,abskk)]
                
                if nneg:
                    # use d^l_nk(beta) = (-1)^{l-k} d^l_{-n,k}(beta)
                    sign_arg += (l-kk)
                if kneg:
                    # use d^l_nk(beta) = (-1)^{l+n} d^l_{n,-k}(beta)
                    sign_arg += (l+nn)                    
                if sign_arg%2:
                    sign = -1
                else:
                    sign = 1
                    
                if beta_step==1:
                    ids[pos] = (sl.start,1,sign)
                    #lookup[(l,nn,kk)] = (sl,sign)
                else:
                    ids[pos] = (sl.stop-1,-1,sign)
                    #lookup[(l,nn,kk)] = (slice(sl.stop,sl.start,-1),sign)
                lookup[(l,n,k)] = ids[pos]
                pos += 1
    return ids,lookup


def wigner_normalization_factor(l):
    '''
    Computes normalization constant for a fixed order l that is present in the wigners computed by make_wigner.py
    wigners = d^l_nk * \sqrt(\frac{2l+1}{2})
    '''
    return np.sqrt((2*l+1)/2)
def wigner_normalization_factors(bw):
    '''
    Creates an array of normalization factors used by the wigner d matrices calculated in make_wigner.py
    wigners = d^l_nk * normalization_factors
    '''
    lnks = coef_id_to_lnk(bw)
    return wigner_normalization_factor(lnks[...,0])

@njit(int64(int64))
def calc_m1_0_fudge(m1):
    if m1 % 2 == 0:
        fudge = 1 
    else:
        fudge = -1
    return fudge

@njit(int64(int64,int64))
def calc_m1_m2_fudge(m1,m2):
    if ((m2-m1) % 2) == 0:
        fudge = 1 
    else:
        fudge = -1
    return fudge


############################################################################
#  ok, the inverse transform
#
#  Function arguments are as follows:
#
#  bw = bandwidth of transform
#  coeffs: plain COMPLEX array of size (4*bw^3-bw)/3, will contain the
#          coefficients of the signal
#  data: FFTW_COMPLEX array of size (2*bw)^3 containing the (output)
#        signal ->
#	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
#	FFTW_COMPLEX!!! (although I will sometimes treat it
#	as the plain COMPLEX array it is)
#        
#  workspace_cx: scratch space FFTW_COMPLEX array of size (2*bw)^3
#	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
#	FFTW_COMPLEX!!! (although I will sometimes treat it
#	as the plain COMPLEX array it is)
#        
#  workspace_cx2: scratch space FFTW_COMPLEX array of size (2*bw)
#	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
#	FFTW_COMPLEX!!! (although I will sometimes treat it
#	as the plain COMPLEX array it is)
#        
#  workspace_re: REAL scratch space of size 12*n + n*bw
#		where n = 2*bw
#                
#  NOTE: I don't need as much tmp space as Forward_SO3_Naive because I'll
#        be able to use the output DATA array as tmp storage
#	within this function. I couldn't use the COEFFS array
#	in Forward_SO3_Naive() because they weren't large enough ... I
#	needed (2 bw)^3 space and they're only (1/3 * bw * ( 4 bw^2 - 1 ))
#        
#  p1: pointer to FFTW plan for correctly ffting WORKSPACE_CX array and
#      placing the result in DATA; this is a FORWARD
#      FFT you're doing!!!
#        
#  wigners: pointer to array of precomputed TRANSPOSED Wigner little-d's
#  
#  flag: = 0 : data is COMPLEX
#        = 1 : data is REAL
@njit(complex128[:](int64,complex128[:],float64[:],int64),cache=True)
def Inverse_SO3_Naive_fft_pc(bw,coeffs,wigners,data_is_complex):
    n = 2*bw
    
    data = np.zeros(n**3,np.complex128)    
    
    # Stage 1: Do the Inverse Wigner transform. The rcoeffs, icoeffs    
    # arrays are assumed to be in the same "arrangement" as that produced
    # by Forward_SO3_Naive().
    # 
    # Since I'm working with two order indeces, m1 and m2, the
    # for-loops will be more intricate than in the case of the
    # "ordinary" spherical transform on S^2.
    # 
    # Also, I will be taking advantage of the symmetries of the
    # Wigner-d functions. As long as I keep my signs and flips
    # right, the Wigner-d's I precompute for an order (m1, m2)
    # transform can generally  be used in seven more transforms:
    # (m1,-m2), (m2,m1), (m2,-m1), (-m2,m1), (-m2,-m1), (-m1,m2)
    # and (-m1,-m2).
    # 
    # 
    # The for-loops will be "designed" as follows. They will be
    # divided into cases according to the orders:
    # 
    # 0) {f_{0,0}} inverse transform
    # 
    # 1) for 0 <= m1 <= bw-1
    # compute inverse transform of
    # i)   {f_{ m1, m1}}
    # ii)  {f_{-m1,-m1}}
    # iii) {f_{-m1, m1}}
    # iv)  {f_{ m1,-m1}}
    # 
    # 2) for 1 <= m1 <= bw-1
    # compute inverse transform of
    # i)   {f_{ m1,  0}}
    # ii)  {f_{-m1,  0}}
    # iii) {f_{  0, m1}}
    # iv)  {f_{  0,-m1}}
    # 
    # 3) for 1 <= m1 <= bw-1
    # for m1+1 <= m2 <= bw-1
    # compute inverse transform 
    # i)    {f_{ m1, m2}}
    # ii)   {f_{-m1,-m2}}
    # iii)  {f_{ m1,-m2}}
    # iv)   {f_{-m1, m2}}
    # v)    {f_{ m2, m1}}
    # vi)   {f_{-m2,-m1}}
    # vii)  {f_{ m2,-m1}}
    # viii) {f_{-m2, m1}}
    # 
    # If assumptions are made regarding the original input signal,
    # e.g. it's strictly real, then one may take advantage of
    # symmetries of the big D wigners (i.e. function of all 3
    # parameters alpha, beta, gamma) and so simplify the for-loops
    # some and hence increase the speed of the program. However,
    # the for-loops to follow will make no such assumptions.
    # Whether the signal is real or complex, these for-loops will
    # handle it.
    # 
    # 
    # Fasten your seatbelt, folks. It's going to be a bumpy ride.


    ########################
    # f_{0,0} coefficients #
    ########################

    #  now, get the locations of where the
    #  samples I have to transform are, and
    #  where the coefficients have to

    sampHere = sampLoc_so3( 0, 0, bw ) 
    coefHere = coefLoc_so3( 0, 0, bw )
    coef = coeffs[coefHere[0]:coefHere[1]]
    wig = wigners[:wigLen_so3(0,0,bw)]
    
    data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(0, 0, bw, coef, wig)
    wignerPos = len(wig)
    
    # 0 <= m1 <= bw-1
    for m1 in range(1,bw):
        wig = wigners[wignerPos:wignerPos+wigLen_so3(m1,m1,bw)] 
    
        ###########################
        # {f_{m1,m1}} coefficient #
        ###########################
        sampHere = sampLoc_so3( m1, m1, bw ) 
        coefHere = coefLoc_so3( m1, m1, bw )    
        coef = coeffs[coefHere[0]:coefHere[1]]
        data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(m1, m1, bw, coef, wig)

        #############################
        # {f_{-m1,-m1}} coefficient #
        #############################
        if data_is_complex:
            sampHere = sampLoc_so3( -m1, -m1, bw ) 
            coefHere = coefLoc_so3( -m1, -m1, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(-m1, -m1, bw, coef, wig)
        else: # use symmetry
            sampHere = sampLoc_so3( m1, m1, bw )
            sampHere2 = sampLoc_so3( -m1, -m1, bw )
            data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()

        ############################
        # {f_{-m1,m1}} coefficient #
        ############################
        sampHere = sampLoc_so3( -m1, m1, bw ) 
        coefHere = coefLoc_so3( -m1, m1, bw )    
        coef = coeffs[coefHere[0]:coefHere[1]]
        data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwY(-m1, m1, bw, coef, wig)
        
        ############################
        # {f_{m1,-m1}} coefficient #
        ############################
        if data_is_complex:
            sampHere = sampLoc_so3( m1, -m1, bw ) 
            coefHere = coefLoc_so3( m1, -m1, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwY(m1, -m1, bw, coef, wig)
        else: # use symmetry
            sampHere = sampLoc_so3( -m1, m1, bw )
            sampHere2 = sampLoc_so3( m1, -m1, bw )
            data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()
            
        wignerPos += len(wig)
    #print(data[-100:])
    # for 1 <= m1 <= bw-1
    for m1 in range(1,bw):
        wig = wigners[wignerPos:wignerPos+wigLen_so3(m1,0,bw)]
        
        ##########################
        # {f_{m1,0}} coefficient #
        ##########################
        sampHere = sampLoc_so3( m1, 0, bw ) 
        coefHere = coefLoc_so3( m1, 0, bw )    
        coef = coeffs[coefHere[0]:coefHere[1]]
        data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(m1, 0, bw, coef, wig)

        ###########################
        # {f_{-m1,0}} coefficient #
        ###########################
        if data_is_complex:
            sampHere = sampLoc_so3( -m1, 0, bw ) 
            coefHere = coefLoc_so3( -m1, 0, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwX(-m1, 0, bw, coef, wig)
        else: # use symmetry
            sampHere = sampLoc_so3( m1, 0, bw )
            sampHere2 = sampLoc_so3( -m1, 0, bw )
            data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()
            
        ##########################
        # {f_{0,m1}} coefficient #
        ##########################
        sampHere = sampLoc_so3( 0, m1, bw ) 
        coefHere = coefLoc_so3( 0, m1, bw )    
        coef = coeffs[coefHere[0]:coefHere[1]]
        data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwX(0, m1, bw, coef, wig)

        ###########################
        # {f_{0,-m1}} coefficient #
        ###########################
        if data_is_complex:
            sampHere = sampLoc_so3( 0, -m1, bw ) 
            coefHere = coefLoc_so3( 0, -m1, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(0, -m1, bw, coef, wig)

        else: # use symmetry
            sampHere = sampLoc_so3( 0, m1, bw )
            sampHere2 = sampLoc_so3( 0, -m1, bw )
            data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()            
            
        wignerPos += len(wig)
    #print(data[-100:])
    # 1 <= m1 <= bw-1
    # m1+1 <= m2 <= bw-1
    for m1 in range(1,bw):
        for m2 in range(m1+1,bw):
            wig = wigners[wignerPos:wignerPos+wigLen_so3(m1,m2,bw)] 

            ###########################
            # {f_{m1,m2}} coefficient #
            ###########################
            sampHere = sampLoc_so3( m1, m2, bw ) 
            coefHere = coefLoc_so3( m1, m2, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(m1, m2, bw, coef, wig)
            
            #############################
            # {f_{-m1,-m2}} coefficient #
            #############################            
            if data_is_complex:
                sampHere = sampLoc_so3( -m1, -m2, bw ) 
                coefHere = coefLoc_so3( -m1, -m2, bw )    
                coef = coeffs[coefHere[0]:coefHere[1]]
                data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwX(-m1, -m2, bw, coef, wig)            
            else: # use symmetry
                sampHere = sampLoc_so3( m1, m2, bw )
                sampHere2 = sampLoc_so3( -m1, -m2, bw )
                data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()

            ###########################
            # {f_{m1,-m2}} coefficient #
            ###########################
            sampHere = sampLoc_so3( m1, -m2, bw ) 
            coefHere = coefLoc_so3( m1, -m2, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwY(m1, -m2, bw, coef, wig)
            
            #############################
            # {f_{-m1,m2}} coefficient #
            #############################
            if data_is_complex:
                sampHere = sampLoc_so3( -m1, m2, bw ) 
                coefHere = coefLoc_so3( -m1, m2, bw )    
                coef = coeffs[coefHere[0]:coefHere[1]]
                data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwY(-m1, m2, bw, coef, wig)     
            else: # use symmetry
                sampHere = sampLoc_so3( m1, -m2, bw )
                sampHere2 = sampLoc_so3( -m1, m2, bw )
                data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()
           
            ################## same but for switched m1, m2 ##################
            
            ###########################
            # {f_{m2,m1}} coefficient #
            ###########################
            sampHere = sampLoc_so3( m2, m1, bw ) 
            coefHere = coefLoc_so3( m2, m1, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwX(m2, m1, bw, coef, wig)
            
            #############################
            # {f_{-m2,-m1}} coefficient #
            #############################
            if data_is_complex:
                sampHere = sampLoc_so3( -m2, -m1, bw ) 
                coefHere = coefLoc_so3( -m2, -m1, bw )    
                coef = coeffs[coefHere[0]:coefHere[1]]
                data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftw(-m2, -m1, bw, coef, wig)            
            else: # use symmetry
                sampHere = sampLoc_so3( m2, m1, bw )
                sampHere2 = sampLoc_so3( -m2, -m1, bw )
                data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()

            ###########################
            # {f_{m2,-m1}} coefficient #
            ###########################
            sampHere = sampLoc_so3( m2, -m1, bw ) 
            coefHere = coefLoc_so3( m2, -m1, bw )    
            coef = coeffs[coefHere[0]:coefHere[1]]
            data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwY(m1, -m2, bw, coef, wig)
            
            #############################
            # {f_{-m2,m1}} coefficient #
            #############################
            if data_is_complex:
                sampHere = sampLoc_so3( -m2, m1, bw ) 
                coefHere = coefLoc_so3( -m2, m1, bw )    
                coef = coeffs[coefHere[0]:coefHere[1]]
                data[sampHere[0]:sampHere[1]] = wigNaiveSynthesis_fftwY(-m1, m2, bw, coef, wig)            
            else: # use symmetry
                sampHere = sampLoc_so3( m2, -m1, bw )
                sampHere2 = sampLoc_so3( -m2, m1, bw )
                data[sampHere2[0]:sampHere2[1]] = data[sampHere[0]:sampHere[1]].conj()
                                            
            wignerPos += len(wig)
    #print(data)


    # matrix n2,n3,n1 ---> n1,n2,n3     
    #data = data.reshape(n**2,n).T.reshape(n,n,n)
    with objmode(data=complex_3d):
        #data = np.fft.fft2(data,norm='ortho')
        data = np.fft.fft2(data.reshape(n,n,n),axes=(0,1),norm='ortho')
    # matrix n1,n3,n2 ---> n1,n2,n3     
    data = np.transpose(data,(0,2,1))    
    # normalization (1/n is contained in numpy fft)
    data = data.flatten() * bw/pi
    #print(data.shape)
    return data

    n = 2*bw
    coeffs = np.zeros(int((4*bw**3-bw)/3 + 0.5),np.complex128)    

####################################################################
#  ok, the forward transform
#
#  Function arguments are as follows:
#
#  bw = bandwidth of transform
#  data: FFTW_COMPLEX array of size (2*bw)^3 containing the input
#        signal ->
#	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
#	FFTW_COMPLEX!!! (although I will sometimes treat it
#	as the plain COMPLEX array it is)
#  coeffs: plain COMPLEX array of size (4*bw^3-bw)/3, will contain the
#          coefficients of the signal
#
#  workspace_cx: scratch space FFTW_COMPLEX array of size (2*bw)^3
#	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
#	FFTW_COMPLEX!!! (although I will sometimes treat it
#	as the plain COMPLEX array it is)
#
#  workspace_cx2: scratch space FFTW_COMPLEX array of size (2*bw)^3
#	MUST BE ALLOCATED BY CALLING FFTW_MALLOC WITH SIZOEF
#	FFTW_COMPLEX!!! (although I will sometimes treat it
#	as the plain COMPLEX array it is)
#                
#  workspace_re: double scratch space of size 12*n + n*bw
#		where n = 2*bw
#
#  weights: ptr to double array of size 2*bw - this array holds
#	   the PRECOMPUTED quadrature weights
#
#  p1: pointer to FFTW plan for correctly ffting WORKSPACE_CX2
#      array and placing the result in WORKSPACE_CX; this is an
#      INVERSE FFT you're doing!!!
#
#  wigners: pointer to array of precomputed Wigner little-d's
#
#  flag: = 0 : data is COMPLEX
#        = 1 : data is double

#@njit(complex128[:](int64,complex128[:],float64[:],float64[:],int64),cache=True)
@njit(cache=True)
def Forward_SO3_Naive_fft_pc(bw,data,weights,wigners,data_is_complex):    
    n = 2*bw
    coeffs = np.zeros(int((4*bw**3-bw)/3 + 0.5),np.complex128)    

    # n1,n2,n3 ------> n1,n3,n2
    data=np.transpose(data.reshape(n,n,n),(0,2,1))
    # n1,n2,n3 ------> n2,n3,n1
    #data = data.reshape(n,n**2).T.reshape(n,n,n)
    with objmode(data=complex_3d):
        data = np.fft.ifft2(data,axes=(0,1),norm='ortho')
    data = data.flatten() * pi/bw
    # normalize data (numpy ifft already contains contains 1/n factor) 

    # Stage 5: Do the Wigner transforms. This is the tricky bit.
    # 
    # Since I'm working with two order indeces, m1 and m2, the
    # for-loops will be more intricate than in the case of the
    # "ordinary" spherical transform on S^2.
    # 
    # Also, I will be taking advantage of the symmetries of the
    # Wigner-d functions. As long as I keep my signs and flips
    # right, the Wigner-d's I precompute for an order (m1, m2)
    # transform can generally  be used in seven more transforms:
    # (m1,-m2), (m2,m1), (m2,-m1), (-m2,m1), (-m2,-m1), (-m1,m2)
    # and (-m1,-m2).
    # 
    # I say "general" because, of course, I'll be transforming
    # at orders (m1,m1), (m1,0) and (0,m1), so I won't get such
    # a huge reduction. Still, this should save time.
    # 
    # If assumptions are made regarding the original input signal,
    # e.g. it's strictly real, then one may take advantage of
    # symmetries of the big D wigners (i.e. function of all 3
    # parameters alpha, beta, gamma) and so simplify the for-loops
    # some and hence increase the speed of the program. However,
    # the for-loops to follow will make no such assumptions.
    # Whether the signal is real or complex, these for-loops will
    # handle it.
    # 
    # The for-loops will be "designed" as follows. They will be
    # divided into cases according to the orders:
    # 
    # 0) {f_{0,0}}
    # 
    # 1) for 0 <= m1 <= bw-1
    # compute the coefficients
    # i)   {f_{ m1, m1}}
    # ii)  {f_{-m1,-m1}}
    # iii) {f_{-m1, m1}}
    # iv)  {f_{ m1,-m1}}
    # 
    # 2) for 1 <= m1 <= bw-1
    # compute the coefficients
    # i)   {f_{ m1,  0}}
    # ii)  {f_{-m1,  0}}
    # iii) {f_{  0, m1}}
    # iv)  {f_{  0,-m1}}
    # 
    # 3) for 1 <= m1 <= bw-1
    # for m1+1 <= m2 <= bw-1
    # compute the coefficients
    # i)    {f_{ m1, m2}}
    # ii)   {f_{-m1,-m2}}
    # iii)  {f_{ m1,-m2}}
    # iv)   {f_{-m1, m2}}
    # v)    {f_{ m2, m1}}
    # vi)   {f_{-m2,-m1}}
    # vii)  {f_{ m2,-m1}}
    # viii) {f_{-m2, m1}}
    # 
    # 
    # Fasten your seatbelt, folks. It's going to be a bumpy ride.

    
    ########################
    # f_{0,0} coefficients #
    ########################

    #  now, get the locations of where the
    #  samples I have to transform are, and
    #  where the coefficients have to

    sampHere = sampLoc_so3( 0, 0, bw ) 
    coefHere = coefLoc_so3( 0, 0, bw )
    dat = data[sampHere[0]:sampHere[1]]
    wig = wigners[:wigLen_so3(0,0,bw)]

    coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(0, 0, bw, dat, wig, weights)
    wignerPos = len(wig)

    # 1 <= m1 <= bw-1
    for m1 in range(1,bw):
        wig = wigners[wignerPos:wignerPos+wigLen_so3(m1,m1,bw)] 
    
        ###########################
        # {f_{m1,m1}} coefficient #
        ###########################
        sampHere = sampLoc_so3( m1, m1, bw ) 
        coefHere = coefLoc_so3( m1, m1, bw )    
        dat = data[sampHere[0]:sampHere[1]]

        coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(m1, m1, bw, dat, wig, weights)

        #############################
        # {f_{-m1,-m1}} coefficient #
        #############################
        if data_is_complex:
            sampHere = sampLoc_so3( -m1, -m1, bw ) 
            coefHere = coefLoc_so3( -m1, -m1, bw )

            dat = data[sampHere[0]:sampHere[1]]

            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(-m1, -m1, bw, dat, wig, weights)
        else: # use symmetry
            coefHere = coefLoc_so3( m1, m1, bw )
            coefHere2 = coefLoc_so3( -m1, -m1, bw )
            coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()

        ############################
        # {f_{-m1,m1}} coefficient #
        ############################
        sampHere = sampLoc_so3( -m1, m1, bw ) 
        coefHere = coefLoc_so3( -m1, m1, bw )
        dat = data[sampHere[0]:sampHere[1]]

        coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwY(-m1, m1, bw, dat, wig, weights)

        ############################
        # {f_{m1,-m1}} coefficient #
        ############################
        if data_is_complex:
            sampHere = sampLoc_so3( m1, -m1, bw ) 
            coefHere = coefLoc_so3( m1, -m1, bw )
            dat = data[sampHere[0]:sampHere[1]]
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwY(m1, -m1, bw, dat, wig, weights)
        else: # use symmetry
            coefHere = coefLoc_so3( -m1, m1, bw )
            coefHere2 = coefLoc_so3( m1, -m1, bw )
            coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()
            
        wignerPos += len(wig)

    #print(data[-100:])
    # for 1 <= m1 <= bw-1
    for m1 in range(1,bw):
        wig = wigners[wignerPos:wignerPos+wigLen_so3(m1,0,bw)]
        
        ##########################
        # {f_{m1,0}} coefficient #
        ##########################
        sampHere = sampLoc_so3( m1, 0, bw ) 
        coefHere = coefLoc_so3( m1, 0, bw )        
        dat = data[sampHere[0]:sampHere[1]]
            
        coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(m1, 0, bw, dat, wig, weights)

        ###########################
        # {f_{-m1,0}} coefficient #
        ###########################
        if data_is_complex:
            sampHere = sampLoc_so3( -m1, 0, bw ) 
            coefHere = coefLoc_so3( -m1, 0, bw )
            dat = data[sampHere[0]:sampHere[1]]
            
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwX(-m1, 0, bw, dat, wig, weights)
            
        else: # use symmetry
            coefHere = coefLoc_so3( m1, 0, bw )
            coefHere2 = coefLoc_so3( -m1, 0, bw )
            fudge = calc_m1_0_fudge(m1)
            coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()*fudge
            
        ##########################
        # {f_{0,m1}} coefficient #
        ##########################
        sampHere = sampLoc_so3( 0, m1, bw ) 
        coefHere = coefLoc_so3( 0, m1, bw )
        dat = data[sampHere[0]:sampHere[1]]
        coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwX(0, m1, bw, dat, wig, weights)

        ###########################
        # {f_{0,-m1}} coefficient #
        ###########################
        if data_is_complex:
            sampHere = sampLoc_so3( 0, -m1, bw ) 
            coefHere = coefLoc_so3( 0, -m1, bw )
            dat = data[sampHere[0]:sampHere[1]]            
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(0, -m1, bw, dat, wig, weights)

        else: # use symmetry
            coefHere = coefLoc_so3( 0, m1, bw )
            coefHere2 = coefLoc_so3( 0, -m1, bw )
            fudge = calc_m1_0_fudge(m1)
            coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()*fudge
            
        wignerPos += len(wig)


    #print(data[-100:])
    # 1 <= m1 <= bw-1
    # m1+1 <= m2 <= bw-1
    for m1 in range(1,bw):
        for m2 in range(m1+1,bw):
            wig = wigners[wignerPos:wignerPos+wigLen_so3(m1,m2,bw)] 

            ###########################
            # {f_{m1,m2}} coefficient #
            ###########################
            sampHere = sampLoc_so3( m1, m2, bw ) 
            coefHere = coefLoc_so3( m1, m2, bw )
            dat = data[sampHere[0]:sampHere[1]]            
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(m1, m2, bw, dat, wig, weights)
            
            #############################
            # {f_{-m1,-m2}} coefficient #
            #############################            
            if data_is_complex:
                sampHere = sampLoc_so3( -m1, -m2, bw ) 
                coefHere = coefLoc_so3( -m1, -m2, bw )
                dat = data[sampHere[0]:sampHere[1]]                
                coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwX(-m1, -m2, bw, dat, wig, weights)            
            else: # use symmetry
                coefHere = coefLoc_so3( m1, m2, bw )
                coefHere2 = coefLoc_so3( -m1, -m2, bw )
                fudge = calc_m1_m2_fudge(m1,m2)
                coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()*fudge

            ############################
            # {f_{m1,-m2}} coefficient #
            ############################
            sampHere = sampLoc_so3( m1, -m2, bw ) 
            coefHere = coefLoc_so3( m1, -m2, bw )
            dat = data[sampHere[0]:sampHere[1]]
            
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwY(m1, -m2, bw, dat, wig, weights)
            
            #############################
            # {f_{-m1,m2}} coefficient #
            #############################
            if data_is_complex:
                sampHere = sampLoc_so3( -m1, m2, bw ) 
                coefHere = coefLoc_so3( -m1, m2, bw )
                dat = data[sampHere[0]:sampHere[1]]                
                coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwY(-m1, m2, bw, dat, wig, weights)            
            else: # use symmetry
                coefHere = coefLoc_so3( m1, -m2, bw )
                coefHere2 = coefLoc_so3( -m1, m2, bw )
                fudge = calc_m1_m2_fudge(m1,m2)
                coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()*fudge
           
            ################## same but for switched m1, m2 ##################
            
            ###########################
            # {f_{m2,m1}} coefficient #
            ###########################
            sampHere = sampLoc_so3( m2, m1, bw ) 
            coefHere = coefLoc_so3( m2, m1, bw )
            dat = data[sampHere[0]:sampHere[1]]            
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwX(m2, m1, bw, dat, wig, weights)
            
            #############################
            # {f_{-m2,-m1}} coefficient #
            #############################
            if data_is_complex:
                sampHere = sampLoc_so3( -m2, -m1, bw ) 
                coefHere = coefLoc_so3( -m2, -m1, bw )    
                dat = data[sampHere[0]:sampHere[1]]                
                coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftw(-m2, -m1, bw, dat, wig, weights)            
            else: # use symmetry
                coefHere = coefLoc_so3( m2, m1, bw )
                coefHere2 = coefLoc_so3( -m2, -m1, bw )
                fudge = calc_m1_m2_fudge(m1,m2)
                coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()*fudge

            ############################
            # {f_{m2,-m1}} coefficient #
            ############################
            sampHere = sampLoc_so3( m2, -m1, bw ) 
            coefHere = coefLoc_so3( m2, -m1, bw )
            dat = data[sampHere[0]:sampHere[1]]            
            coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwY(m1, -m2, bw, dat, wig, weights)
            
            ############################
            # {f_{-m2,m1}} coefficient #
            ############################
            if data_is_complex:
                sampHere = sampLoc_so3( -m2, m1, bw ) 
                coefHere = coefLoc_so3( -m2, m1, bw )
                dat = data[sampHere[0]:sampHere[1]]                
                coeffs[coefHere[0]:coefHere[1]] = wigNaiveAnalysis_fftwY(-m1, m2, bw, dat, wig, weights)            
            else: # use symmetry
                coefHere = coefLoc_so3( m2, -m1, bw )
                coefHere2 = coefLoc_so3( -m2, m1, bw )
                fudge = calc_m1_m2_fudge(m1,m2)
                coeffs[coefHere2[0]:coefHere2[1]] = coeffs[coefHere[0]:coefHere[1]].conj()*fudge
                                            
            wignerPos += len(wig)

    return coeffs


#######################################
# let f,g be two square integrable functions on the 2 sphere 
# Define C: SO(3) ---> R,   R |---> <f,g \circ R> = \int_{S^2} dx f(x)*\overline{g(Rx)}
# This function calculates C(R) for all R defined in make_SO3_grid
#
# arguments :
#   f_coeff: f_{l,m} spherical harmonic coefficients of f 
#   g_coeff: g_{l,m} spherical harmonic coefficients of g 
#            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
#   split _ids: ids that split coefficients in 2*bw+1 sub arrays indexed by m
#
# output :
#   C_values: array of shape (2*bw,2*bw,2*bw)
@njit()
def calc_C_ml(bw,f_coeff,g_coeff,split_ids,wigners_transposed,data_is_complex):
    n=2*bw
    coeff = np.zeros(totalCoeffs_so3(bw),np.complex128)
    f_coeff = np.split(f_coeff,split_ids)
    g_coeff = np.split(g_coeff,split_ids)
    for l in range(bw):
        wigNorm = 2.*pi*np.sqrt(2./(2.*l+1.))
        for m1 in range(-l,l+1):
            lm1_id = l-abs(m1)
            f_lm1 = f_coeff[m1][lm1_id]#[-m1][l]
            #/* need to reset the -1 fudge factor */
            if ( (m1 + l) % 2 ):
                fudge = -1 
            else:
                fudge = 1
            
            for m2 in range(-l,l+1):
                lm2_id = l-abs(m2)
                #first take the CONJUGATE of the pattern g_coef
                #and multiply it by the fudge factor
                g_lm2 = fudge * f_coeff[-m2][lm2_id].conjugate()
                
                #/* now multiply the signal f_coef by the pattern g_coef,
                # and save it in the so3 coefficient array */
                
                index = so3CoefLoc(m1,m2,l,bw)
                coeff[index] = wigNorm*f_lm1*g_lm2
                #/* multiply fudge factor by -1 */
                fudge *= -1
                #coeff_pos += 1
    #print('nans in coeffs = ',np.isnan(coeff).any())
    #print('max abs real coeff = ',np.max(np.abs(coeff.real)))
    #print('max abs imag coeff = ',np.max(np.abs(coeff.imag)))
    C = Inverse_SO3_Naive_fft_pc(bw,coeff,wigners_transposed,data_is_complex)
    #print('max abs = ',np.max(np.abs(C)))
    #print('arg max abs = ',np.argmax(np.abs(C)))
    return C


#######################################
# let f,g be two square integrable functions on the 2 sphere 
# Define C: SO(3) ---> R,   R |---> <f,g \circ R> = \int_{S^2} dx f(x)*\overline{g(Rx)}
# This function calculates C(R) for all R defined in make_SO3_grid
#
# arguments :
#   f_coeff: f_{l,m} spherical harmonic coefficients of f 
#   g_coeff: g_{l,m} spherical harmonic coefficients of g 
#            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
#   split _ids: ids that split coefficients in 2*bw+1 sub arrays indexed by m
#
# output :
#   C_values: array of shape (2*bw,2*bw,2*bw)
#   The maximum in C corresponds to the Rotation that best maps f to g
#
#   The Idis of the maximum in C i_max,j_max,k_max correspond to the euler anglees beta,alpha,gamma in this order.
#   Somehow there is still a bug. The resulting euler angles need to be modified as follows to yield correct results
#   alpha,beta,gamma -----> 2*pi - alpha, beta , 2*pi-gamma
@njit()
def calc_C_lm(bw,f_coeff,g_coeff,split_ids,wigners_transposed,data_is_complex):    
    coeff = np.zeros(totalCoeffs_so3(bw),np.complex128)
    f_coeff = np.split(f_coeff,split_ids)
    g_coeff = np.split(g_coeff,split_ids)
    #    coeff_pos = 0
    for l in range(bw):
        wigNorm = 2.*pi*np.sqrt(2./(2.*l+1.))
        f_l = f_coeff[l]
        g_l = g_coeff[l]
        m0_pos=l # position of the m = 0 component in arrays with components -l,...,l
        for m1 in range(-l,l+1):
            f_lm1 = f_l[m0_pos - m1] #should give f_{l-m1}
            for m2 in range(-l,l+1):
                #first take the CONJUGATE of the pattern g_coef              
                g_lm2 =  g_l[m0_pos - m2].conjugate() #fudge* g_l[m0_pos - m2].conjugate() # f_l[m0_pos - m2] should give f_{l-m2}                
                #/* now multiply the signal f_coef by the pattern g_coef,
                # and save it in the so3 coefficient array */                
                index = so3CoefLoc(m1,m2,l,bw)
                coeff[index] = wigNorm * f_lm1 * g_lm2 * (-1)**(m1-m2)
    # Calculate C by inverse SOFT
    C = Inverse_SO3_Naive_fft_pc(bw,coeff,wigners_transposed,data_is_complex)

    return C


#######################################
# Calculates the mean of above C over 0- axis of f_coeff,g_coeff.
# let f,g be two square integrable functions on the 2 sphere 
# Define C: SO(3) ---> R,   R |---> <f,g \circ R> = \int_{S^2} dx f(x)*\overline{g(Rx)}
# This function calculates C(R) for all R defined in make_SO3_grid
#
# arguments :
#   f_coeff: f_{l,m} spherical harmonic coefficients of f 
#   g_coeff: g_{l,m} spherical harmonic coefficients of g 
#            f_coeff, g_coeff are complex numpy arrays of shape (N,bw*(bw+1)+bw+1)
#            for an asbitrary integer N.
#   split _ids: ids that split coefficients in 2*bw+1 sub arrays indexed by m
#   r_cutoff_id: int64, 0<=r_cutoff_id<=N only the first r_cutoff_id radial coordinates are averaged over
#
@njit()
def calc_mean_C_array(bw,f_coeff,g_coeff,r_lower_limit,r_upper_limit,split_ids,wigners_transposed,data_is_complex):
    n= 2*bw
    C = np.zeros(n**3,np.complex128)
    n_steps = r_upper_limit - r_lower_limit
    for i in range(r_lower_limit,r_upper_limit):
        #print('find_rotation step',i+1,' of ',n_steps)
        tmp = calc_C_lm(bw,f_coeff[i],g_coeff[i],split_ids,wigners_transposed,data_is_complex)
        #print_angles(bw,tmp)
        C += tmp        
    return C/n_steps


#########################################
# Calculate Rotated spherical harmonic coefficients
# F_{lm} ----> \sum_n D^l_{m,n}(\alpha,\beta,\gamma) F_{ln}
# 
#
#@njit()
#def rotate_harmonic_coefficients(bw,f_coeff,split_ids,euler_angle_ids,wigners,exponentials):
#    rotated_coeff = np.split(np.zeros_like(f_coeff),split_ids)
#    f_coeff = np.split(f_coeff,split_ids)
#    print(f_coeff)
#    wigner_pos = 0
#    n = 2*bw
#    for m2 in range(-bw+1,bw):        
#        for m1 in range(-bw+1,bw):
#            lm1_start = max(0,abs(m1)-abs(m2))
#            lm2_start = max(0,abs(m2)-abs(m1))
#            a_id,b_id,g_id = euler_angle_ids
#            wigner_D = exponentials[m1][a_id]*wigners[:,b_id] *exponentials[m2][g_id]
#            #print("m: ",m,' n: ',n)
#            #print('len _wigner: ',stop-start, ' lm_start: ',lm_start,' ln_start: ',ln_start)
#            #print(wigner_D*f_coeff[n][ln_start:])
#            rotated_coeff[m2][lm2_start:] += wigner_D*f_coeff[m1][lm1_start:]
#    return rotated_coeff

####################################################################
# The following functions are just for testing purposes
@njit()
def print_angles(bw,C_vals):
    n = 2*bw
    argmax = np.argmax(np.abs(C_vals))
    ii = argmax//(n**2)
    jj = (argmax - ii*n**2)//n
    kk = (argmax-ii*n**2-jj*n)
    beta =(2*ii+1)/(4*bw)
    alpha = 2 - jj/bw
    gamma = 2 - kk/bw
    #print("angles ids = ",(jj,ii,kk))
    #print("angles = ",(alpha,beta,gamma),' * pi')

    
# inverse fft is part of forward SOFT
@njit(cache=True)
def inverse_fft(bw,data):
    n = 2*bw
    # n1,n2,n3 ------> n2,n3,n1
    data = data.reshape(n,n**2).T.reshape(n,n,n)
    with objmode(data=complex_3d):
        data = np.fft.ifft2(data,axes=(0,1),norm='ortho')
    
    data = data.flatten()*pi/bw
    return data

# forward fft is part of inverse SOFT
@njit(cache=True)
def forward_fft(bw,data):
    n = 2*bw
    # matrix n2,n3,n1 ---> n1,n2,n3     
    data = data.reshape(n**2,n).T.reshape(n,n,n)
    with objmode(data=complex_3d):
        data = np.fft.fft2(data,norm='ortho')
        
    # normalization (1/n is contained in numpy fft)
    data = data.flatten() * bw/pi #ffted_data * bw/pi
    return data



def zeros_order_forward_so3(signal,wigner000,weights):
    N=len(signal)
    temp = np.sum(signal,axis = (0,2))*1/N*np.pi/(N//2) # Why use bw = N//2  and not N ??? maybe error in soft rewrite from above ? NO is fine should be 2*np.pi/(2*bw)=np.pi/bw = np.pi/(N//2). The additional factor of 1/N is to comply with the numpy inverse fourier transform used in Forward_SO3_Naive_fft_pc
    return wigner000 @ (temp*weights)
def integrate_over_so3_normalized(signal,wigner000,weights):
    N=len(signal)
    temp = np.sum(signal,axis = (0,2))/N**2 # 
    return np.sum(temp*weights/2) #
def integrate_over_so3(signal,wigner000,weights):
    N=len(signal)
    temp = np.sum(signal,axis = (0,2))*(2*np.pi/N)**2 # 
    return np.sum(temp*weights) #

def get_so3_probabily_weights(legendre_weights):
    N = len(legendre_weights)
    alpha_gamma_weights = np.full(N,1/N**2)
    beta_weights = legendre_weights
    return alpha_gamma_weights[:,None,None]*beta_weights[None,:,None]*alpha_gamma_weights[:,None,None]

##########################################

def list_of_lnk(bw):
    L=bw-1
    ls = []
    ns = []
    ks = []
    for l in range(bw):
        for n in range(-l,l+1):
            for k in range(-l,l+1):
                ls.append(l)
                ns.append(n)
                ks.append(k)
    return np.stack([ls,ns,ks]).T

##########################################


class Soft():
    def __init__(self,bw):
        self.bw = bw
        tmp = self.generate_data()
        self.grid = self.make_SO3_grid()
        self.coeff_grid = self.make_coeff_grid()
        self.wigners = tmp[0]
        self.wigners_transposed = tmp[1] #small d matrices with m1,m2 swapped
        self.legendre_weights = tmp[2]
        self.n_coeffs = totalCoeffs_so3(bw)
        self.n_points = (2*bw)**3
        self.so_shape=(2*bw,)*3
        self.wigners000=self.wigners[:2*bw]
        
    def generate_data(self):
        band_width = self.bw
        wigners = genWigAll(band_width)
        wigners_transposed = genWigAllTrans(band_width)
        legendre_weights = makeweights2(band_width)
        return wigners, wigners_transposed, legendre_weights
    
    def forward_cmplx(self,data):        
        if data.ndim>1:
            data=data.flatten()
        if (data.dtype != np.dtype(np.complex128)):
            print('Warning: input dtype is {} but {} is required. Trying to change dtype.'.format(data.dtype,np.dtype(np.complex128)))
            data = data.astype(np.complex128)
        if (len(data) != self.n_points):
            raise AssertionError('input data needs to be of length {}. Given length = {}'.format(self.n_points,data.dtype,len(data)))        
        coeffs = Forward_SO3_Naive_fft_pc(self.bw,data,self.legendre_weights,self.wigners,True)        
        return coeffs
    
    def inverse_cmplx(self,coeff):
        if (coeff.dtype != np.dtype(np.complex128)):
            print('Warning: input dtype is {} but {} is required. Trying to change dtype.'.format(coeff.dtype,np.dtype(np.complex128)))
            coeff = coeff.astype(np.complex128)
        if (len(coeff) != self.n_coeffs) or (coeff.dtype != np.dtype(np.complex128)):
            raise AssertionError('input coeffs need to be of type complex128 and of length {}. Given data dtype = {} length = {}'.format(self.n_coeffs,coeff.dtype,len(coeff)))        
        data = Inverse_SO3_Naive_fft_pc(self.bw,coeff,self.wigners_transposed,True)
        return data.reshape(self.so_shape)

    def make_SO3_grid(self):
        alpha,beta,gamma = get_euler_angles(self.bw)

        grid = np.stack(np.meshgrid(alpha,beta,gamma,indexing='ij'),3)
        #grid = GridFactory.construct_grid('uniform',[alpha,beta,gamma])[:]
        return grid

    def make_coeff_grid(self):
        bw= self.bw
        lnks = list_of_lnk(bw)
        #print(lnks.shape)
        coeff_grid = np.zeros((totalCoeffs_so3(bw),3))
        for lnk in lnks:
            l,n,k = lnk
            _id = so3CoefLoc(n,k,l,bw)
            #print(_id)
            coeff_grid[_id]=lnk
        return coeff_grid
    
    def get_empty_coeff(self):
        return np.zeros(self.n_coeffs,dtype = complex)
        
    def integrate_over_so3(self,signal):
        return integrate_over_so3(signal,self.wigners000,self.legendre_weights)

    #################################################
    # computes the rotated harmonic coefficients
    # f_{l,m} ---> \sum_n D^l_{n,m} f_{l,n}
    # assumes f_{l,m}=coeff to be a complex array of length bw**2
    # split_ids such that np.split(coeff,split_ids) gives a list of coefficients indexed by l.
    # euler_angles is an array of length 3 containing the euler angles (alpha,beta,gamma) in ZYZ format
    def rotate_coeff(self,coeff,split_ids,euler_angles):
        #print('coeff.shape before rotation = {}'.format(coeff[0].shape))
        rotated_coeff = rotate_coeff_multi(self.bw, coeff,split_ids, euler_angles)
        return rotated_coeff

    def rotate_coeff_single(self,coeff,split_ids,euler_angles):
        #print('coeff.shape before rotation = {}'.format(coeff.shape))
        rotated_coeff = rotate_coeff(self.bw, coeff,split_ids, euler_angles)
        return rotated_coeff

    
    #######################################
    # let f,g be two square integrable functions on the 2 sphere 
    # Define C: SO(3) ---> R,   R |---> <f,g \circ R> = \int_{S^2} dx f(x)*\overline{g(Rx)}
    # This function calculates C(R) for all R defined in make_SO3_grid
    #
    # arguments :
    #   f_coeff: f_{l,m} spherical harmonic coefficients of f 
    #   g_coeff: g_{l,m} spherical harmonic coefficients of g 
    #            f_coeff, g_coeff are complex numpy arrays of shape bw*(bw+1)+bw+1
    #   split _ids: ids that split coefficients in 2*bw+1 sub arrays indexed by m
    #
    def calc_mean_C(self,f_coeff,g_coeff,r_split_ids,ml_split_ids):
        r_split_lower=r_split_ids[0]
        r_split_upper=r_split_ids[1]
        mean_C = calc_mean_C_array(self.bw,f_coeff,g_coeff,r_split_lower,r_split_upper,ml_split_ids,self.wigners_transposed,True)
        return mean_C


    #testing 
    def combine_coeffs(self,f_coeff,g_coeff,split_ids):
        return combine_harmonic_coeffs(self.bw,f_coeff,g_coeff,split_ids)
