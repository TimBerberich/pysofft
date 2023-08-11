import rotate
import numpy as np
from xframe2.externalLibraries.shtns_plugin import sh

bw=4
shtns = sh(3)

split_ids = np.array([l**2+2*l+1 for l in  range(bw-1)])
coeff = np.random.rand(bw**2) + np.random.rand(bw**2)*1j

r_coeff = rotate.rotate_coeff(bw,coeff,split_ids,np.array([0.,0.,0.]))
