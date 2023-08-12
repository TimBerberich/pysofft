import numpy as np
from numba import njit, objmode, types,int64,float64,complex128


complex_3d = types.complex128[:,:,:]
pi = np.pi

@njit()
def wignerdmat(L,matIn,trigs,sqrts):
    cosVal = trigs[0]
    sinVal = trigs[1]
    matOut = np.zeros((2*L+1)**2)
    if L>1:
        tmp = np.zeros((2*L)**2)
        tmp[:(2*L-1)**2] = matIn
        for deg in range(2*L-1,2*L+1):
            rdeg = 1/deg
            tmpdim = deg + 1
            matOut[:tmpdim**2]=0
            for i in range(deg):
                for j in range(deg):
                    matOut[(i*tmpdim)+j] += rdeg*sqrts[deg-i]*sqrts[deg-j]* tmp[(i*deg)+j]*cosVal 
                    matOut[((i+1)*tmpdim)+j] -= rdeg*sqrts[i+1]*sqrts[deg-j]*  tmp[(i*deg)+j]*sinVal 
                    matOut[(i*tmpdim)+(j+1)] +=  rdeg*sqrts[deg-i]*sqrts[j+1]* tmp[(i*deg)+j]*sinVal
                    matOut[((i+1)*tmpdim)+(j+1)] +=  rdeg*sqrts[i+1]*sqrts[j+1]* tmp[(i*deg)+j]*cosVal
                    #print("i,j = ",i,j, 'c: ',tmp[(i*deg)+j])            
            if deg == 2*L-1:
                tmp[:] = matOut[:(2*L)**2]
    elif L == 0:
        matOut[0] = 1
    else:
        matOut[0] = cosVal*cosVal 
        matOut[1] = sqrts[2]*cosVal*sinVal 
        matOut[2] = sinVal*sinVal 
        
        matOut[3] = -matOut[1] 
        matOut[4] = matOut[0]-matOut[2]
        matOut[5] = matOut[1] 
        
        matOut[6] = matOut[2] 
        matOut[7] = -matOut[1] 
        matOut[8] = matOut[0]
    #print('matOut', matOut)
    return matOut


#Assumes f_{L,m} is given as coeff[m] for m in -L<=m<=L
#Assumes expA,expG is given as expA[m] for m in -Lmax<=m<=Lmax (where Lmax>=L)
@njit()
def rotateCoefDegree(L,coeff,wigner,expA,expG):
    dim = 2*L+1
    Lp1 = L+1
    #print('alive')
    D = np.zeros((2*L+1,2*L+1)).astype(np.complex128)
    #assemble D matrix D_m1,m2 = D[m1,m2] -L<=m1,m2<=L
    #texpG = np.array([expG[m] for m in range(-L,L+1)])
        
    for i in range(dim):
        #m=i-L
        #print(expA[m],np.roll(wigner[i*dim:(i+1)*dim],Lp1).shape,texpG.shape)
        D[i,:]= expA[i]*wigner[i*dim:(i+1)*dim]*expG
    rotated_coeff = np.dot(D.T,coeff)
    #print(rotated_coeff.shape)
    return rotated_coeff
    
    


############################################
#
#
#
#
#
@njit()
def rotate_coeff(bw, coeff,split_ids, euler_angles):
    n = 2*bw
    alpha,beta,gamma = euler_angles
    sqrts = np.sqrt(np.arange(n))
    trigs = np.array([np.cos(beta/2),np.sin(beta/2)])
    ms = np.arange(-bw+1,bw)
    #print('len ms = ',len(ms))
    expA = np.exp(-1j*ms*alpha)
    expG = np.exp(-1j*ms*gamma)

    rotated_coeff=np.zeros_like(coeff)
    coeff = np.split(coeff,split_ids)
    matIn = np.zeros((n-1)**2)
    matOut= np.zeros_like(matIn)
    #    matOut= matIn + (n-1)**2
    wig = np.zeros(1)

    pos=0
    for deg in range(bw):
        dim = 2*deg + 1
        if deg>0:
            matIn = wig.copy()
        #grab the spherical coefficients at this degree
        l_coeff = coeff[deg]
        #l_coeff = np.roll(l_coeff,bw)
        wig = wignerdmat(deg,matIn,trigs,sqrts)
        exp_start =  bw-deg-1
        exp_stop= exp_start + 2*deg+1
        rotated_coeff[pos: pos+dim] = rotateCoefDegree(deg,l_coeff,wig,expA[exp_start:exp_stop],expG[exp_start:exp_stop])
        pos+=dim
    return rotated_coeff

@njit()
def rotate_coeff_multi(bw, coeffs,split_ids, euler_angles):
    rotated_coeff = np.zeros_like(coeffs)
    for i in range(len(coeffs)):
        rotated_coeff[i] = rotate_coeff(bw, coeffs[i],split_ids, euler_angles)
    return rotated_coeff
