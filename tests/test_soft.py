import time
import sys
import numpy as np

from pysofft.make_wiegner import CosEvalPts,CosEvalPts2,SinEvalPts,SinEvalPts2,genWigTrans_L2
from pysofft.make_wiegner import genWigAll,genWigAllTrans
from pysofft.wignerTransform import wigNaiveSynthesis_fftw
from pysofft.wignerWeights import makeweights2
from pysofft.soft import Inverse_SO3_Naive_fft_pc,Forward_SO3_Naive_fft_pc,coefLoc_so3,sampLoc_so3,forward_fft,inverse_fft,forward_wigner,inverse_wigner 
from c_arrays import tcoeff_forward, tcoeff_forward_2, tdata_inverse, tdata_inverse_ffted, tdata_inverse_ffted_2 



def test_forward():
    bw = 5
    n=2*bw
    cs = CosEvalPts(2*bw)
    si2 = SinEvalPts2(2*bw)
    cs2 = CosEvalPts2(2*bw)
    
    wigners = genWigAll(bw)
    wignersT = genWigAllTrans(bw)
    weights = makeweights2(bw)
    
    n_total_coeffs = int((4 * bw * bw * bw - bw ) / 3 +0.5)

    #in_data = np.full(n**3,1+1j)#np.random.rand(n**3)+np.random.rand(n**3)*1j
    in_data = np.arange(n**3)+1j #np.random.rand(n**3)+np.random.rand(n**3)*1j
    in_data[0]=0
    in_data = tdata_inverse_ffted_2
    tcoeff = np.arange(n_total_coeffs) +1j
    
    data_is_complex = True
    outcoeff = Forward_SO3_Naive_fft_pc(bw,in_data,weights,wigners,data_is_complex)
    print(outcoeff)

    is_close_mask = np.isclose(outcoeff,tcoeff,atol = 1e-15)
    if is_close_mask.all():
        print("test passed!")
    else:
        print("TEST FAILED!")
        print("failing indices are \n {}".format(~is_close_mask))
        print('max error ={}'.format(np.max(np.abs(tcoeff_forward_2-outcoeff))))
    return locals()

def test_inverse():
    bw = 5
    n=2*bw
    
    wignersT = genWigAllTrans(bw)
    
    n_total_coeffs = int((4 * bw * bw * bw - bw ) / 3 +0.5)
    
    in_coeff = np.full(n_total_coeffs,1+1j).astype(np.complex128) #np.random.rand(n**3)+np.random.rand(n**3)*1j
    in_coeff = np.arange(n_total_coeffs) +1j
    #in_coeff[:bw*n]=1+1j
    
    data_is_complex = True
    data = Inverse_SO3_Naive_fft_pc(bw,in_coeff,wignersT,data_is_complex)
    
    is_close_mask = np.isclose(tdata_inverse_ffted_2,data,atol = 1e-15)
    if is_close_mask.all():
        print("test passed!")
    else:
        print("TEST FAILED!")
        print("failing indices are \n {}".format(~is_close_mask))
        print('max error ={}'.format(np.max(np.abs(tdata_inverse_ffted_2-data))))
    return locals()

def test_IF_loop(iterations = 1):
    bw = 5
    n=2*bw
    start = time.time()

    wigners = genWigAll(bw)
    wignersT = genWigAllTrans(bw)
    weights = makeweights2(bw)    
    n_total_coeffs = int((4 * bw * bw * bw - bw ) / 3 +0.5)

    #in_coeff = np.full(n_total_coeffs,1+1j)#np.random.rand(n**3)+np.random.rand(n**3)*1j
    #in_coeff = np.arange(n_total_coeffs) +1.j
    in_coeff = np.random.rand(n_total_coeffs)+0j#+np.random.rand(n_total_coeffs)*1j
    data_is_complex = True
    coeff = in_coeff.copy()
    print('data generation took {}s'.format(time.time()-start))

    inverse_time = 0
    forward_time = 0
    for i in range(iterations):
        print('current iteration = {} of {}'.format(i+1,iterations))        
        start_inverse = time.time()
        data = Inverse_SO3_Naive_fft_pc(bw,coeff,wignersT,data_is_complex)
        inverse_time += time.time()-start_inverse
        start_forward = time.time()
        coeff = Forward_SO3_Naive_fft_pc(bw,data,weights,wigners,data_is_complex)
        forward_time += time.time()-start_forward
        
    print('forward SOFT took on average {}s using {} iterations'.format(forward_time/iterations,iterations))
    print('inverse SOFT took on average {}s using {} iterations'.format(inverse_time/iterations,iterations))    
    is_close_mask = np.isclose(coeff,in_coeff,atol = 1e-15)
    if is_close_mask.all():
        print("test passed!")
    else:
        print("TEST FAILED!")
        print("failed indices mask \n {}".format(~is_close_mask))
        print('max error ={}'.format(np.max(np.abs(in_coeff-coeff))))
    return locals()

# incoeff = np.random.rand(n_total_coeffs)+np.random.rand(n_total_coeffs)*1j
# incoeff = np.full(n_total_coeffs,1+1j)
# in_data = np.full(n**3,1+1j)#np.random.rand(n**3)+np.random.rand(n**3)*1j
# in_data[:1]=0
# data_is_complex = True
# data = Inverse_SO3_Naive_fft_pc(bw,incoeff,wignersT,data_is_complex)
# outcoeff = Forward_SO3_Naive_fft_pc(bw,data,weights,wigners,data_is_complex)

#data = Inverse_00(bw,incoeff,wignersT,data_is_complex)
#outcoeff = Forward_00(bw,in_data,weights,wigners,data_is_complex)

#out_data= fft_test(bw,in_data)

if __name__ == "__main__" :
    #locals().update(test_forward())
    #locals().update(test_inverse())
    locals().update(test_IF_loop(iterations = 2))
