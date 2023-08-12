import numpy as np
from numba import njit, objmode, types


pi = np.pi

#  makeweights: given a bandwidth bw, make weights for
#  both even *and* odd order Legendre transforms.
#  
#  bw -> bandwidth of transform
#  weights -> pointer to double array of length 4*bw; this
#             array will contain the even and odd weights;
#             even weights start at weight[0], and odd weights
#             start at weights[2*bw]
@njit()
def makeweights(bw):
    weights = np.zeros(4*bw)
    fudge = pi/(4*bw)
    for j in range(2*bw):
        tmpsum = 0
        for k in range(bw):
            tmpsum += 1/(2*k+1) * np.sin((2*j+1) * (2*k+1) * fudge)
        tmpsum *= np.sin((2*j+1)*fudge)
        tmpsum *= 2/bw

        weights[j] = tmpsum
        weights[j + 2*bw] = tmpsum * np.sin((2*j+1)*fudge)
    return weights

@njit()
def makeweights2(bw):
    weights = np.zeros(2*bw)
    fudge = pi/(4*bw)
    for j in range(2*bw):
        tmpsum = 0
        for k in range(bw):
            tmpsum += 1/(2*k+1) * np.sin((2*j+1) * (2*k+1) * fudge)
        tmpsum *= np.sin((2*j+1)*fudge)
        tmpsum *= 2/bw

        weights[j] = tmpsum
    return weights
