import numpy as np
from numba import njit, objmode, types

pi = np.pi

###################################################
## returns an array of n Chebyshev nodes, i.e.     
##   cos( EvalPts )
@njit()
def CosEvalPts(N):
    two_N = 2*N
    eval_pts = (2*np.arange(N)+1)*pi/two_N
    return np.cos(eval_pts)


###################################################
## returns an array of n Chebyshev nodes, i.e.     
##     sin( EvalPts )
@njit()
def SinEvalPts(N):
    two_N = 2*N
    eval_pts = (2*np.arange(N)+1)*pi/two_N
    return np.sin(eval_pts)

#################################################################
##  returns an array of n slightly changed Chebyshev nodes, i.e. 
##   cos( EvalPts/2 )
@njit()
def CosEvalPts2(N):
    two_N = 2*N
    eval_pts = (2*np.arange(N)+1)*pi/two_N
    return np.cos(eval_pts/2)

#################################################################
##  returns an array of n slightly changed Chebyshev nodes, i.e. 
##   sin( EvalPts/2 )
@njit()
def SinEvalPts2(N):
    two_N = 2*N
    eval_pts = (2*np.arange(N)+1)*pi/two_N
    return np.sin(eval_pts/2)


#########################################################################
## Recurrence coefficients 
#########################################################################
#
#  Recurrence coefficients for the L2-normed Wigner d functions.
#  When using these coefficients, make sure that the initial
#  d^l_{m1,m2} function is also normed.
#  
#  NOTE: These functions are normed in the sense that, as functions
#  only of theta,
#  
#  \int_0^pi (d_{m1,m2}^l) (d_{m1,m2}^lx) \sin\theta d\theta = delta_{l,lx}
#  
#  SO I'm treating these functions "on their own" and not as part of
#  the big D function (a function of all three Euler angles):
#  
#  bigD_{m1,m2}^{j}(alpha, theta, gamma) =
#  exp(-i * m1 * alpha ) * d_{m1,m2}^l(theta) * exp(-i * m2 * gamma )
#  
#  NOTE: when I say "wigner", I mean the little d function (of one angle)
#  and not the bigD function (three angles)
#  
#  If
#  
#  normWig[l, m1, m2] = Sqrt(2/(2 l + 1))
#  ( so \int_0^pi (d_{m1,m2}^l / normWig[l,m1,m2])^2 \sin\theta d\theta = 1 )
# 
#  then
#  
#  normD[l, m1, m2] = 2 * PI * normWig(l, m1, m2)
#  
#  
#  
#  Notation: j - degree
#            m1, m2, orders
#
#  
#  d_{m1,m2}^{j+1}(theta) =
#          L2_aN_so3(j, m1, m2) * d_{m1,m2}^{j-1}(theta) +
#          L2_bN_so3(j, m1, m2) * cos(theta) * d_{m1,m2}^{j}(theta) +
#          L2_cN_so3(j, m1, m2) * d_{m1,m2}^{j}(theta)
@njit()
def L2_aN_so3(j,m1,m2):
    dj = j
    dm1 = m1
    dm2 = m2

    if j > 0:
        return -np.sqrt( (2*dj + 3)/(2*dj - 1) ) * (dj + 1)/np.sqrt(((dj+1)*(dj+1) - dm1*dm1) * ((dj+1)*(dj+1) - dm2*dm2)) * np.sqrt( (dj*dj - dm1*dm1)*(dj*dj - dm2*dm2) ) / dj
    else:
        return 0

@njit()
def L2_bN_so3(j,m1,m2):
    dj = j
    dm1 = m1
    dm2 = m2
  
    return  np.sqrt( (2*dj + 3)/(2*dj + 1) ) * (dj + 1)*(2*dj + 1)/np.sqrt(((dj+1)*(dj+1) - dm1*dm1) * ((dj+1)*(dj+1) - dm2*dm2))
    
@njit()
def L2_cN_so3(j,m1,m2):
    dj = j
    dm1 = m1
    dm2 = m2

    if j != 0:
        return -L2_bN_so3(j, m1, m2) * dm1 * dm2 / ( dj * ( dj + 1 ) )
    else:
        return 0
    

#################################################################
# L2 normed wigner little d, WHERE THE DEGREE OF THE FUNCTION IS EQUAL
# TO ONE OF ITS ORDERS. This is the function to use when starting the
# three-term recurrence at orders (m1,m2)
# 
#   arguments: m1, m2 - order of the function
#   sinEval - sine values of evaluation pts
#   cosEval - cosine values of the evaluation pts
#   n - how many points there are, i.e. length of sinEval, cosEval
#   result - where to place the result, length n
#   
#   
#   Note that, by definition, since I am starting the recurrence with this
#   function, that the degree j of the function is equal to max(abs(m1), abs(m2) ).
#   
#   
#   This has been tested and stably computes Pmm functions thru bw=512
@njit()
def wigSpec_L2(m1,m2,sinEval,cosEval,n):    
    absM1 = abs(m1)
    absM2 = abs(m2)
    l = max(absM1,absM2)
    delta = l - min(absM1,absM2)
    
    dl = l
    dm1 = m1
    dm2 = m2
    sinSign = 1
    normFactor = 1

    for i in range(delta):
        normFactor *= np.sqrt( (2.*dl - i)/( i + 1.) )

    # need to adjust to make the L2-norm equal to 1 
    normFactor *= np.sqrt((2*dl+1)/2)

    if l == absM1:
        if m1 >= 0:
            cosPower = l + m2
            sinPower = l - m2
            if (l-m2)%2:
                sinSign = -1
        else:
            cosPower = l - m2
            sinPower = l + m2
    elif m2 >= 0:
        cosPower = l + m1
        sinPower = l - m1
    else:
        cosPower = l - m1
        sinPower = l + m1
        if (l+m1)%2:
            sinSign = -1

    dCP = cosPower
    dSP = sinPower

    result = np.zeros(n)
    for i in range(n):
        result[i] = normFactor * sinSign * sinEval[i]**dSP * cosEval[i]**dCP
        
    return result


#################################################################
#  Given orders m1, m2, and a bandwidth bw, this function will
#  generate the all the Wigner little d functions whose orders
#  are (m1, m2) and degrees are j = max(|m1|, |m2|) through j = bw - 1
#  using the 3-term recurrence. 
#  
#  All of these Wigners will have L2 norm = 1
#  
#  
#  let j = max(|m1|, |m2|)
#  
#  The functions generated will be
#  
#  d_{m1,m2}^j, d_{m1,m2}^{j+1}, ..., d_{m1,m2}^{bw-1}
#  
#  Each of these functions will be evaluated at the n = 2*bw-many
#  points
#  
#  pi*(2 * [0..n-1] + 1) / ( 2 * n )
#  
#  If theta(k) = pi*(2*k+1)/(2*n), then what's returned will be the
#  array
#  
#  d_{m1,m2}^j(theta(0)) ... d_{m1,m2}^{bw-1}(theta(0))
#  d_{m1,m2}^j(theta(1)) ... d_{m1,m2}^{bw-1}(theta(1))
#  d_{m1,m2}^j(theta(2)) ... d_{m1,m2}^{bw-1}(theta(2)) ...
#  d_{m1,m2}^j(theta(n-1)) ... d_{m1,m2}^{bw-1}(theta(n-1))
#  
#  arguments: m1, m2 = orders of the functions
#             bw = bandwidth
#             sinEval = sine values of evaluation pts, length n = 2*bw
#             cosEval = cosine values of the evaluation pts, length n = 2*bw
#               
#             I need the two arrays for wigSpec purposes
#              
#             sinEval2 = sine values of evaluation pts, length n = 2*bw
#             cosEval2 = cosine values of the evaluation pts, length n = 2*bw
#               
#             result = array to store result, length (bw-m)*n   (as (bw-m,n) matrix)
#                      where m = max( |m1|, |m2| ); 
#             workspace = scratch area, length 8 * n = 16 * bw ;
#  
#  The routine won't be efficient, but it'll work.
@njit()
def genWig_L2(m1,m2,bw,cosEval,sinEval2,cosEval2):
    tmpM1 = abs(m1)
    tmpM2 = abs(m2)
    l = max(tmpM1,tmpM2)
    m = l
    n = 2*bw

    result= np.zeros((bw-m)*n)
    workspace = np.zeros((7,n))
    
    # generate and save first Wigner fct by call to wigSpec_L2     
    workspace[1] = wigSpec_L2( m1, m2,sinEval2, cosEval2, n)
    result[:n] = workspace[1]
    
    # now start doing the recurrence
    for i in range(bw-m-1):
        workspace[2] = workspace[0]*L2_aN_so3(m+i,m1,m2)
        workspace[3] = workspace[1]*L2_bN_so3(m+i,m1,m2)
        workspace[4] = workspace[3]*cosEval
        workspace[5] = workspace[1]*L2_cN_so3(m+i,m1,m2)
        workspace[6] = workspace[2] + workspace[4] + workspace[5] # workspace[6] contains d_{m1,m2}^{m+i+1}
        # store the function values
        result[(i+1)*n:(i+2)*n] = workspace[6]
        
        # update for the next iteration of the recurrence
        workspace[0] = workspace[1]
        workspace[1] = workspace[6]
                    
    return result

@njit()
def genWigTrans_L2(m1,m2,bw,cosEval,sinEval2,cosEval2):
    tmpM1 = abs(m1)
    tmpM2 = abs(m2)
    m = max(tmpM1,tmpM2)
    n = 2*bw
    W = genWig_L2(m1,m2,bw,cosEval,sinEval2,cosEval2).reshape(bw-m,n)
    return W.T.flatten()


#######################################################################
# create_SO3_grid: C
#
# arguments:
#    bw: bandwidth of the transform
# output: the corresponding euler angles for j in [0,2*bw[
#    alpha = 2*\pi*j/(2*bw)
#    beta  = \pi(2*j+1)/(4*bw)
#    gamma = 2*\pi*j/(2*bw)
#
@njit()
def get_euler_angles(bw):
    j=np.arange(int(2*bw))
    alpha =  2*np.pi*j/(2*bw)
    gamma = alpha.copy()
    beta = np.pi*(2*j+1)/(4*bw)
    return alpha,beta,gamma

########################################################################
# generates the exponentials that are needed to calculate the big wigner D
# from the small d matrices generated in above methods via the formula
#    D^l_{m,n} = e^{-i m \alpha} d^l_{m,n}(\beta) e^{-i n \gamma}
# 
#
@njit()
def get_exponentials(bw):
    alpha,beta,gamma = get_euler_angles(bw)
    ls = np.arange(bw)
    ms = np.concatenate((ls,-1*ls[::-1]))
    return np.exp(-1.j*np.expand_dims(ms,-1)*np.expand_dims(alpha,0))

#######################################################################
# genAllWig: make ALL the Wigner little-d's necessary to do a full
#            FORWARD SOFT (i.e. SO(3)) transform. Designed to be used in
#            the SOFT routines which rely on the symmetries of the
#            Wigner little-d's.
# 
#  arguments:
# 
#   bw: bandwidth of transform
# 
#   wigners: array to store ALL the wigners, of size (gulp)
# 
#      1/3 * bw^2 * (2 + 3*bw + bw^2)
# 
#   workspace: scratch space, of size 12 * n, where n = 2*bw
def genWigAll(bw):
    n = 2 * bw
    wigners = np.zeros( int(0.5+1/3 * bw**2 * (2 + 3*bw + bw**2)) )
    
    # precompute appropriate sine and cosine values 
    sinPts = SinEvalPts(n)
    cosPts = CosEvalPts(n)
    sinPts2 = SinEvalPts2(n)
    cosPts2 = CosEvalPts2(n)
    
    # precompute Wigner little-d's for m1 = m2 = 0
    pt = 0
    wigners[:bw*n] = genWig_L2( 0, 0, bw, cosPts, sinPts2, cosPts2)
    pt += bw*n
    
    # precompute Wigner little-d's for abs(m1)=abs(m2)
    for m1 in range(1,bw):
        tmp = genWig_L2( m1, m1,bw, cosPts, sinPts2, cosPts2)
        size = len(tmp)
        wigners[pt:pt+size] = tmp
        pt += size

    # precompute Wigner little-d's for one order being 0,
    for m1 in range(1,bw):
        tmp = genWig_L2( m1, 0,bw, cosPts, sinPts2, cosPts2)
        size = len(tmp)
        wigners[pt:pt+size] = tmp
        pt += size        

    # precompute Wigner little-d's for m1, m2
    for m1 in range(1,bw):
        for m2 in range(m1+1,bw):
            tmp = genWig_L2( m1, m2,bw, cosPts, sinPts2, cosPts2)
            size = len(tmp)
            wigners[pt:pt+size] = tmp
            pt += size
    return wigners


#######################################################################
# genAllWigTrans: make ALL the Wigner little-d's necessary to do a full
#          INVERSE SOFT (i.e. SO(3)) transform. Designed to be used in
#	   the SOFT routines which rely on the symmetries of the
#	   Wigner little-d's.
# 
#  arguments:
# 
#   bw: bandwidth of transform
# 
#   wigners: array to store ALL the wigners, of size (gulp)
# 
#      1/3 * bw^2 * (2 + 3*bw + bw^2)
# 
#   workspace: scratch space, of size 12 * n, where n = 2*bw
def genWigAllTrans(bw):
    n = 2 * bw
    wigners = np.zeros( int(0.5+1/3 * bw**2 * (2 + 3*bw + bw**2)) )

    # precompute appropriate sine and cosine values 
    sinPts = SinEvalPts(n)
    cosPts = CosEvalPts(n)
    sinPts2 = SinEvalPts2(n)
    cosPts2 = CosEvalPts2(n)

    # precompute Wigner little-d's for m1 = m2 = 0
    pt = 0
    wigners[:bw*n] = genWigTrans_L2( 0, 0, bw, cosPts, sinPts2, cosPts2)
    pt += bw*n

    # precompute Wigner little-d's for abs(m1)=abs(m2)
    for m1 in range(1,bw):
        tmp = genWigTrans_L2( m1, m1,bw, cosPts, sinPts2, cosPts2)
        size = len(tmp)
        wigners[pt:pt+size] = tmp
        pt += size

    # precompute Wigner little-d's for one order being 0,
    for m1 in range(1,bw):
        tmp = genWigTrans_L2( m1, 0,bw, cosPts, sinPts2, cosPts2)
        size = len(tmp)
        wigners[pt:pt+size] = tmp
        pt += size        

    # precompute Wigner little-d's for m1, m2
    for m1 in range(1,bw):
        for m2 in range(m1+1,bw):
            tmp = genWigTrans_L2( m1, m2,bw, cosPts, sinPts2, cosPts2)
            size = len(tmp)
            wigners[pt:pt+size] = tmp
            pt += size
    return wigners


if __name__ == '__main__':

    bw = 10
    w = genWigAllTrans(bw)
   # bw=10
   # m1=0
   # m2=0
    
   # c = CosEvalPts(2*bw)
   # s2 = SinEvalPts2(2*bw)
   # c2 = CosEvalPts2(2*bw)
   # r = genWig_L2(m1,m2,bw,c,s2,c2)
   # rt = genWigTrans_L2(m1,m2,bw,c,s2,c2)
   # print(r)
