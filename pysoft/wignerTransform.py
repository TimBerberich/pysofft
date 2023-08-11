import numpy as np
from numba import njit, objmode, types
from numba.types import float64,int64,complex128


#   wigNaiveAnalysis
#
#   Given a bandwidth bw, and orders m1, m2, this function will analyse
#   the signal of length n = 2*bw by projecting it onto the L2-normed
#   little d wigners of degrees l = Max(|m1|, |m2|) ... bw - 1, i.e. I
#   am doing the integrals
#
#   <signal, d_{m1,m2}^l> for l = Max(|m1|, |m2|) ... bw - 1
#
#   where the little d is normalized such that
#
#   \int_0^pi (d_{m1,m2}^l) (d_{m1,m2}^{lx}) \sin\theta d\theta = \delta_{l,lx}
#
#   NOTE: to use this routine as part of the full transform on SO(3), I still
#   need to multiply these coefficients by (PI/bw) ... the PI is because
#   my basis functions are the big D's now, the bw because of the adjustment
#   in scaling doing the Fourier transform introduces. This multiplying HAS
#   TO HAPPEN outside of this routine.
#
#
#   NOTE that this function was written to be part of a full
#   SO(3) harmonic transform, so a certain amount of precomputation
#   has been assumed.
#   
#   arguments: m1, m2 = orders of the transform
#              bw = bandwidth
#              signal = array of sample values, length n = 2*bw
#              wigners = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
#                        this array holds the wigner little d's that I
#                        need - PRECOMPUTED by genWig_L2() 
#              coeffs = array of length (bw - max(|m1|,|m2|)), to store the
#                        final results
#              workspace = scratch area, of length n
#              weights = ptr to length 2*bw array containing the
#                        quadrature weights - PRECOMPUTED by makeweights()


@njit(complex128[:](int64,int64,int64,complex128[:],float64[:],float64[:]))
def wigNaiveAnalysis_fftw(m1,m2,bw,signal,wigners,weights):
    #  m is the degree of the "first" wigner function at
    # this order */
    m = max( abs(m1),abs(m2) )
    n = 2*bw
    size = bw-m
    # weight the signal with the appropriate quadrature weights
    weighted_signal = weights*signal

    # and now I start analysing, i.e. multiply the matrix "wigners" by
    # the vector weighted_signal
    coeffs = np.zeros(size,dtype = np.complex128)
    for i in range(size):
        wigners_part = wigners[i*n:(i+1)*n] 
        tmpA = 0
        for j in range(n):
            tmpA += wigners_part[j] * weighted_signal[j]
        
        coeffs[i] = tmpA
    # and that should be all
    return coeffs


#   wigNaiveAnalysis_fftwX
#
#   Given a bandwidth bw, and orders m1, m2, this function will analyse
#   the signal of length n = 2*bw by projecting it onto the L2-normed
#   little d wigners of degrees l = Max(|m1|, |m2|) ... bw - 1, i.e. I
#   am doing the integrals
#
#   <signal, d_{m1,m2}^l> for l = Max(|m1|, |m2|) ... bw - 1
#
#   where the little d is normalized such that
#
#   \int_0^pi (d_{m1,m2}^l) (d_{m1,m2}^{lx}) \sin\theta d\theta = \delta_{l,lx}
#
#   NOTE: to use this routine as part of the full transform on SO(3), I still
#   need to multiply these coefficients by (PI/bw) ... the PI is because
#   my basis functions are the big D's now, the bw because of the adjustment
#   in scaling doing the Fourier transform introduces. This multiplying HAS
#   TO HAPPEN outside of this routine.
#
#
#   NOTE that this function was written to be part of a full
#   SO(3) harmonic transform, so a certain amount of precomputation
#   has been assumed.
#
#   NOTE: This routine assumes that the SIGNS of m1 and m2 are
#         IDENTICAL. This routine is written under the assumption
#	 that I precomputed the Wigner-ds for order (|m1|,|m2|),
#	 i.e. both positive. Since I could conceivably use this
#	 routine for (|m2|,|m1|), (-|m1|,-|m2|) or (-|m2|,-|m1|), 
#	 I need to compute a sign "fudge factor" -> I need to get
#	 the correct power of -1!!!
#
#
#   arguments: m1, m2 = orders of the transform
#              bw = bandwidth
#	      signal = COMPLEX array of sample values, length n = 2*bw
#	      wigners = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
#                        this array holds the wigner little d's that I
#			need - PRECOMPUTED by genWig_L2() 
#	      weights -> ptr to double array of size 2*bw - this array holds
#	         the PRECOMPUTED quadrature weights
#	      coeffs = COMPLEX array of length (bw - max(|m1|,|m2|)), to
#	               store the final results
#	      workspace = COMPLEX scratch area, of length n
#


@njit(complex128[:](int64,int64,int64,complex128[:],float64[:],float64[:]))
def wigNaiveAnalysis_fftwX(m1,m2,bw,signal,wigners,weights):
    #  m is the degree of the "first" wigner function at
    # this order */    
    m = max( abs(m1),abs(m2) )
    n = 2*bw
    size = bw-m

    if abs( m1 - m2 ) % 2:
        fudge = -1 
    else:
        fudge = 1 
    
    # weight the signal with the appropriate quadrature weights
    weighted_signal = weights*signal

    # and now I start analysing, i.e. multiply the matrix "wigners" by
    # the vector weighted_signal
    coeffs =np.zeros(size,dtype = np.complex128)
    for i in range(size):
        wigners_part = wigners[i*n:(i+1)*n] 
        tmpA = 0
        for j in range(n):
            tmpA += wigners_part[j] * weighted_signal[j]
        coeffs[i] = tmpA*fudge
    # and that should be all
    return coeffs



#  wigNaiveAnalysis_fftwY
#
#   Given a bandwidth bw, and orders m1, m2, this function will analyse
#   the signal of length n = 2*bw by projecting it onto the L2-normed
#   little d wigners of degrees l = Max(|m1|, |m2|) ... bw - 1, i.e. I
#   am doing the integrals
#
#   <signal, d_{m1,m2}^l> for l = Max(|m1|, |m2|) ... bw - 1
#
#   where the little d is normalized such that
#
#   \int_0^pi (d_{m1,m2}^l) (d_{m1,m2}^{lx}) \sin\theta d\theta = \delta_{l,lx}
#
#   NOTE: to use this routine as part of the full transform on SO(3), I still
#   need to multiply these coefficients by (PI/bw) ... the PI is because
#   my basis functions are the big D's now, the bw because of the adjustment
#   in scaling doing the Fourier transform introduces. This multiplying HAS
#   TO HAPPEN outside of this routine.
#
#
#   NOTE that this function was written to be part of a full
#   SO(3) harmonic transform, so a certain amount of precomputation
#   has been assumed.
#
#   NOTE: This routine assumes that the SIGNS of m1 and m2 are
#         DIFFERENT. This routine is written under the assumption
#	 that I precomputed the Wigner-ds for order (|m1|,|m2|),
#	 i.e. both positive. Since I could conceivably use this
#	 routine for (-|m1|,|m2|), (-|m2|,|m1|), (-|m1|,|m2|) or
#	 (-|m2|,|m1|), 
#	 I need to compute a sign "fudge factor" -> I need to get
#	 the correct power of -1!!! Note that, given the identities,
#	 I need to FLIP either the signal or wigner-d's.
#
#
#   arguments: m1, m2 = orders of the transform
#              bw = bandwidth
#	      signal = COMPLEX array of sample values, length n = 2*bw
#	      wigners = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
#                        this array holds the wigner little d's that I
#			need - PRECOMPUTED by genWig_L2() 
#	      weights -> ptr to double array of size 2*bw - this array holds
#	         the PRECOMPUTED quadrature weights
#	      coeffs = COMPLEX array of length (bw - max(|m1|,|m2|)), to
#	               store the final results
#	      workspace = COMPLEX scratch area, of length n


@njit(complex128[:](int64,int64,int64,complex128[:],float64[:],float64[:]))
def wigNaiveAnalysis_fftwY(m1,m2,bw,signal,wigners,weights):
    #  m is the degree of the "first" wigner function at
    # this order */    
    m = max( abs(m1),abs(m2) )
    n = 2*bw
    size = bw-m

    if  m1 < 0:
        if (m - m2) % 2:
            fudge = -1 
        else:
            fudge = 1 
    else:
        if  (m + m1) % 2:
            fudge = -1 
        else:
            fudge = 1 
    
    # weight the signal with the appropriate quadrature weights
    weighted_signal = weights*signal

    # and now I start analysing, i.e. multiply the matrix "wigners" by
    # the vector weighted_signal
    coeffs =np.zeros(size,dtype=np.complex128)
    for i in range(size):
        wigners_part = wigners[i*n:(i+1)*n] 
        tmpA = 0
        for j in range(n):
            tmpA += wigners_part[j] * weighted_signal[n-1-j]
        coeffs[i] = tmpA*fudge
        fudge *= -1
    # and that should be all
    return coeffs


#   wigNaiveSynthesis
#
#   Given a bandwidth bw, and orders m1, m2, this function will synthesize
#   the signal of length n = 2*bw by "summing up the coefficients". More
#   plainly, this is the inverse transform of wigNaiveAnalysis.
#
#   Let l = Max(|m1|, |m2|). In matrix-lingo, wigNaiveAnalysis may be
#   written as:
#
#   c = P W f
#
#   where f is the data vector, W is the quadrature matrix (i.e. weights),
#   P is the (bw-l) x n matrix of sample values of the L2-normed wigners
#   d_{m1,m2}^l d_{m1,m2}^{l+1} ... d_{m1,m2}^{bw-1}, and c is the
#   wigner series representation (i.e. coefficients) of f (c is
#   a vector of length bw-l).
#
#   So wigNaiveSynthesis can be written as
#
#   f = Transpose(P) c
#
#   No quadrature matrix is necessary.
#
#   NOTE that this function was written to be part of a full
#   SO(3) harmonic transform, so a certain amount of precomputation
#   has been assumed.
#
#   NOTE: to use this routine as part of the full transform on SO(3), I still
#   need to multiply these coefficients by (PI/bw) ... the PI is because
#   my basis functions are the big D's now, the bw because of the adjustment
#   in scaling doing the Fourier transform introduces. This multiplying HAS
#   TO HAPPEN outside of this routine.
#
#   arguments: m1, m2 = orders of the transform
#              bw = bandwidth
#              coeffs = array of coefficients, length (bw - max(|m1|,|m2|))
#              wignersTrans = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
#                             this array holds the wigner little d's that I
#                             need at this order - PRECOMPUTED genWigTrans_L2()
#              signal = array of length n = 2*bw, to store the final results,
#                       the reconstructed sample values
#              workspace = scratch area, of length 0 * n
#                          (that's right, 0 * n ... I'm keeping this
#                          argument here just so that this function
#                          call looks identical to wigNaiveAnalysis)

@njit(complex128[:](int64,int64,int64,complex128[:],float64[:]))
def wigNaiveSynthesis_fftw(m1,m2,bw,coeff,wignersTrans):
    #  m is the degree of the "first" wigner function at
    # this order */    
    m = max( abs(m1),abs(m2) )
    n = 2*bw
    size = bw-m
    signal = np.zeros(2*bw,dtype=np.complex128)
    for i in range(n):
        wignersTrans_part = wignersTrans[i*size:(i+1)*size] 
        tmpA = 0
        for j in range(size):
            tmpA += wignersTrans_part[j] * coeff[j]
        signal[i] = tmpA
    # and that should be all
    return signal

#   wigNaiveSynthesis_fftwX
#
#   Given a bandwidth bw, and orders m1, m2, this function will synthesize
#   the signal of length n = 2*bw by "summing up the coefficients". More
#   plainly, this is the inverse transform of wigNaiveAnalysis.
#
#   Let l = Max(|m1|, |m2|). In matrix-lingo, wigNaiveAnalysis may be
#   written as:
#
#   c = P W f
#
#   where f is the data vector, W is the quadrature matrix (i.e. weights),
#   P is the (bw-l) x n matrix of sample values of the L2-normed wigners
#   d_{m1,m2}^l d_{m1,m2}^{l+1} ... d_{m1,m2}^{bw-1}, and c is the
#   wigner series representation (i.e. coefficients) of f (c is
#   a vector of length bw-l).
#
#   So wigNaiveSynthesis can be written as
#
#   f = Transpose(P) c
#
#   No quadrature matrix is necessary.
#
#   NOTE that this function was written to be part of a full
#   SO(3) harmonic transform, so a certain amount of precomputation
#   has been assumed.
#
#   NOTE: to use this routine as part of the full transform on SO(3), I still
#   need to multiply these coefficients by (PI/bw) ... the PI is because
#   my basis functions are the big D's now, the bw because of the adjustment
#   in scaling doing the Fourier transform introduces. This multiplying HAS
#   TO HAPPEN outside of this routine.
#
#
#   NOTE: This routine assumes that the SIGNS of m1 and m2 are
#         IDENTICAL. This routine is written under the assumption
#	 that I precomputed the Wigner-ds for order (|m1|,|m2|),
#	 i.e. both positive. Since I could conceivably use this
#	 routine for (|m2|,|m1|), (-|m1|,-|m2|) or (-|m2|,-|m1|), 
#	 I need to compute a sign "fudge factor" -> I need to get
#	 the correct power of -1!!!
#
#
#   arguments: m1, m2 = orders of the transform
#              bw = bandwidth
#	      coeffs = COMPLEX array of coefficients,
#	               length (bw - max(|m1|,|m2|))
#	      wignersTrans = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
#                             this array holds the wigner little d's that I
#			     need at this order - PRECOMPUTED genWigTrans_L2()
#	      signal = COMPLEX array of length n = 2*bw, to store the
#	               final results, the reconstructed sample values
#	      workspace = scratch area, of length 0 * n
#                          (that's right, 0 * n ... I'm keeping this
#			   argument here just so that this function
#			   call looks identical to wigNaiveAnalysis)

@njit(complex128[:](int64,int64,int64,complex128[:],float64[:]))
def wigNaiveSynthesis_fftwX(m1,m2,bw,coeff,wignersTrans):
    #  m is the degree of the "first" wigner function at
    # this order */    
    m = max( abs(m1),abs(m2) )
    n = 2*bw
    size = bw-m

    if abs( m1 - m2 ) % 2:
        fudge = -1 
    else:
        fudge = 1
        
    signal = np.zeros(2*bw,dtype = np.complex128)
    for i in range(n):
        wignersTrans_part = wignersTrans[i*size:(i+1)*size] 
        tmpA = 0
        for j in range(size):
            tmpA += wignersTrans_part[j] * coeff[j]
        signal[i] = tmpA*fudge
    # and that should be all
    return signal


#  wigNaiveSynthesis_fftwY
#
#   Given a bandwidth bw, and orders m1, m2, this function will synthesize
#   the signal of length n = 2*bw by "summing up the coefficients". More
#   plainly, this is the inverse transform of wigNaiveAnalysis.
#
#   Let l = Max(|m1|, |m2|). In matrix-lingo, wigNaiveAnalysis may be
#   written as:
#
#   c = P W f
#
#   where f is the data vector, W is the quadrature matrix (i.e. weights),
#   P is the (bw-l) x n matrix of sample values of the L2-normed wigners
#   d_{m1,m2}^l d_{m1,m2}^{l+1} ... d_{m1,m2}^{bw-1}, and c is the
#   wigner series representation (i.e. coefficients) of f (c is
#   a vector of length bw-l).
#
#   So wigNaiveSynthesis can be written as
#
#   f = Transpose(P) c
#
#   No quadrature matrix is necessary.
#
#   NOTE that this function was written to be part of a full
#   SO(3) harmonic transform, so a certain amount of precomputation
#   has been assumed.
#
#   NOTE: to use this routine as part of the full transform on SO(3), I still
#   need to multiply these coefficients by (PI/bw) ... the PI is because
#   my basis functions are the big D's now, the bw because of the adjustment
#   in scaling doing the Fourier transform introduces. This multiplying HAS
#   TO HAPPEN outside of this routine.
#
#   NOTE: This routine assumes that the SIGNS of m1 and m2 are
#         DIFFERENT. This routine is written under the assumption
#	 that I precomputed the Wigner-ds for order (|m1|,|m2|),
#	 i.e. both positive. Since I could conceivably use this
#	 routine for (-|m1|,|m2|), (-|m2|,|m1|), (-|m1|,|m2|) or
#	 (-|m2|,|m1|), 
#	 I need to compute a sign "fudge factor" -> I need to get
#	 the correct power of -1!!! Note that, given the identities,
#	 I need to FLIP either the signal or wigner-d's.
#
#   arguments: m1, m2 = orders of the transform
#              bw = bandwidth
#	      coeffs = COMPLEX array of coefficients,
#	               length (bw - max(|m1|,|m2|))
#	      wignersTrans = array of length (bw - max(|m1|,|m2|))*( 2 * bw );
#                             this array holds the wigner little d's that I
#			     need at this order - PRECOMPUTED genWigTrans_L2()
#	      signal = COMPLEX array of length n = 2*bw, to store the
#	               final results, the reconstructed sample values
#	      workspace = COMPLEX scratch area, of length 1 * n

@njit(complex128[:](int64,int64,int64,complex128[:],float64[:]))
def wigNaiveSynthesis_fftwY(m1,m2,bw,in_coeff,wignersTrans):
    #  m is the degree of the "first" wigner function at
    # this order */    
    m = max( abs(m1),abs(m2) )
    n = 2*bw
    size = bw-m

    if  m1 < 0:
        if (m - m2) % 2:
            fudge = -1 
        else:
            fudge = 1 
    else:
        if  (m + m1) % 2:
            fudge = -1 
        else:
            fudge = 1
            
    coeff = in_coeff.copy()
    for i in range(size):
      coeff[i] *= fudge
      fudge *= -1 

    signal = np.zeros(2*bw,dtype = np.complex128)
    for i in range(n):
        #wignersTrans_part = wignersTrans[i*size:(i+1)*size]
        wignersTrans_part = wignersTrans[(n - i)*size - size: (n - i)*size]
        tmpA = 0
        for j in range(size):
            tmpA += wignersTrans_part[j] * coeff[j]
        signal[i] = tmpA
    # and that should be all
    return signal.copy()
