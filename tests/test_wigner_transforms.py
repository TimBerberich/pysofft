import numpy as np
from pysofft.make_wiegner import CosEvalPts,CosEvalPts2,SinEvalPts,SinEvalPts2,genWigTrans_L2,genWig_L2
from pysofft.wignerTransform import wigNaiveSynthesis_fftw,wigNaiveAnalysis_fftw
from pysofft.wignerTransform import wigNaiveSynthesis_fftwX,wigNaiveAnalysis_fftwX
from pysofft.wignerTransform import wigNaiveSynthesis_fftwY,wigNaiveAnalysis_fftwY
from pysofft.wignerWeights import makeweights2
bw = 5
cs = CosEvalPts(2*bw)
si2 = SinEvalPts2(2*bw)
cs2 = CosEvalPts2(2*bw)

m1 = 1
m2 = -3 
wignersT = genWigTrans_L2(m1,m2,bw,cs,si2,cs2)
wigners = genWig_L2(m1,m2,bw,cs,si2,cs2)
weights = makeweights2(bw)


# test coeffs (are only real in the C tests)
coeffs = np.array([
    0.000272879644053, #0.000272879644053,
    0.392601207495749,
]).astype(np.complex128)

data = wigNaiveSynthesis_fftw(m1,m2,bw,coeffs,wignersT)

coeffs_2 = wigNaiveAnalysis_fftw(m1,m2,bw,data,wigners,weights)
diff = coeffs-coeffs_2

#Resuts from C routine test_genWig for m1=1 m2=-3 bw=5
c_data = np.array([
    0.000410911594818,
    0.028245186685174,
    0.154460845028305,
    0.336313302784190,
    0.368693739070161,
    0.116571851969348,
    -0.259138736025564,
    -0.429647297137376,
    -0.275152781041740,
    -0.039521272387210,
]).astype(np.complex128)


is_close_mask = np.isclose(data,c_data,atol=1e-15)
if is_close_mask.all():
    print('test passed!')
else:
    print('TEST FAILED!')
    failing_indices = np.nonzero(~is_close_mask) # the '~' means logical not
    print('Differences detected in array indices \n {}'.format(failing_indices))
    
## the following does not make much sense for given m1,m2 since here fftw=fftwX=fftwY routines
dataX = wigNaiveSynthesis_fftwX(m1,m2,bw,coeffs,wignersT)
dataY = wigNaiveSynthesis_fftwY(m1,m2,bw,coeffs,wignersT)

coeffs_2X = wigNaiveAnalysis_fftwX(m1,m2,bw,data,wigners,weights)
coeffs_2Y = wigNaiveAnalysis_fftwY(m1,m2,bw,data,wigners,weights)


diff_X = coeffs-coeffs_2X
diff_Y = coeffs-coeffs_2Y
