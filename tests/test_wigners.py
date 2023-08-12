import numpy as np
from pysofft.make_wiegner import CosEvalPts,CosEvalPts2,SinEvalPts,SinEvalPts2,genWig_L2,genWigAll,genWigAllTrans

bw = 5
cs = CosEvalPts(2*bw)
si2 = SinEvalPts2(2*bw)
cs2 = CosEvalPts2(2*bw)

m1 = 1
m2 = -3 
wigners = genWig_L2(m1,m2,bw,cs,si2,cs2)

#Resuts from C routine test_genWig for m1=1 m2=-3 bw=5
c_wigners = np.array([
    0.000272879644053,
    0.020346218162353,
    0.132638311882183,
    0.392601208495749,
    0.745329561274668,
    1.021763878146642,
    1.045473434192605,
    0.773072734754658,
    0.353001232422972,
    0.044055774208316,
    0.001046448973227,
    0.071929566199260,
    0.393337176226859,
    0.856355412644639,
    0.938586906588108,
    0.296211603248880,
    -0.660782541750537,
    -1.094897936246974,
    -0.701090831959001,
    -0.100695803924461,
])

if np.isclose(wigners,c_wigners,atol=1e-15).all():
    print('test passed!')
else:
    print('TEST FAILED!')
    failing_indices = np.nonzero(wigners-c_wigners)
    print('Differences detected in array indices \n {}'.format(c_wigners))
    
