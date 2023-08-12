#from pySOFT import soft as soft
#from pySOFT import make_wiegner as make_wigner
#from pySOFT import wignerWeights as wignerWeights
#from pySOFT import wignerTransform as wignerTransform

from pysofft import soft 
from pysofft import make_wiegner 
from pysofft import wignerWeights
from pysofft import wignerTransform


#The following is a bit of a hack
#In calling numba functions they will get compiled when the package is loaded and not when they are called the first time.
#import numpy as np
#bw = 2
#data = np.zeros(4**3)+0.j
#f_coeff = np.zeros((2,2+2))+0.j
#g_coeff = np.zeros((2,2+2))+0.j
#split_ids=np.array([2,3],dtype=int)
#r_limit_ids = [0,2]
#
#wigners = make_wigner.genWigAll(bw)
#wignersT = make_wigner.genWigAllTrans(bw)
#exponentials = make_wigner.get_exponentials(bw)
#weights = wignerWeights.makeweights2(bw)
#coef = soft.Forward_SO3_Naive_fft_pc(bw,data,weights,wigners,True)
#soft.Inverse_SO3_Naive_fft_pc(bw,data,wignersT,True)
#soft.calc_mean_C_array(bw,f_coeff,g_coeff,r_limit_ids,split_ids,wignersT,True)
##x = soft.rotate_harmonic_coefficients(bw,f_coeff,split_ids,(0,0,0),wigners,exponentials)
