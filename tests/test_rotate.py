import time
import sys
import numpy as np

from pysofft.make_wiegner import CosEvalPts,CosEvalPts2,SinEvalPts,SinEvalPts2,genWigTrans_L2
from pysofft.make_wiegner import genWigAll,genWigAllTrans,get_exponentials
from pysofft.wignerTransform import wigNaiveSynthesis_fftw
from pysofft.wignerWeights import makeweights2
from pysofft.soft import Inverse_SO3_Naive_fft_pc,Forward_SO3_Naive_fft_pc,coefLoc_so3,sampLoc_so3,forward_fft,inverse_fft,forward_wigner,inverse_wigner, rotate_harmonic_coefficients
from c_arrays import tcoeff_forward, tcoeff_forward_2, tdata_inverse, tdata_inverse_ffted, tdata_inverse_ffted_2 


def rotate(bw =2):
    data = np.zeros((2*bw)**3)+0.j
    f_coeff = np.zeros(bw**2)+0.j
    ls = np.arange(bw+1)
    index = (ls*(ls+1)/2).astype(int)
    split_ids = np.concatenate((index[-1] - index[-2::-1],index[-1]+index[1:-2]))
    
    wigners = genWigAll(bw)
    wignersT = genWigAllTrans(bw)
    exponentials = get_exponentials(bw)
    weights = makeweights2(bw)    
    x = rotate_harmonic_coefficients(bw,f_coeff,split_ids,(0,0,0),wigners,exponentials)
    return locals()

if __name__ == "__main__" :
    locals().update(rotate(bw = 10))
