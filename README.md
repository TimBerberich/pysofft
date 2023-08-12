# PySOFFT
This package allows to compute Fourier transforms on the rotation group SO(3).\
It is a partial python port of soft-2.0 released by Peter Kostelec and Daniel Rockmore.

For details see their paper:\
FFTs on the Rotation Group\
J Fourier Anal Appl (2008) 14: 145–179\
DOI 10.1007/s00041-008-9013-5

The present numba code is a more or less straight copy/adaptation of the original C implementation.\
PySOFT is made available with conset of the original authors and under the same GPL3 license.


## Installation
```
pip install pysofft
```
Dependencies are: numpy and numba
## Usage
Creating a Transform class for a fixed bandwidth (maximal considered spherical harmonic order +1)
```
from pysofft.soft import Soft

bandwidth = 32
soft_obj = Soft(bandwidth)
```
SO(3) grid points (gives as (2 bandwidth,2 bandwidth,2 bandwidth,3) of euler angles $\alpha,\beta,\gamma$):
```
soft_obj.grid
```
FFT Coefficient grid points labeled by $n,m,l$ corresponding to the Wigner-D matrix $D^l_{n,m}$
```
soft_obj.coeff_grid
```

Creating random coefficients and computing the FFT on SO(3)
```
import numpy as np
Flnm=np.random.rand(*soft_obj.coeff_grid.shape[:-1]).astype('complex')

f = soft_obj.inverse_cmplx(Flnm) # inverse Transform   (Dlnm Coefficients -> SO(3))
Flnm2 = soft_obj.forward_cmplx(f) # forward Transform  (SO(3) -> Dlnm Coefficients)
```

## Known Issues
There are some bugs when spawning child processes after instanciation of the class Soft.
So, as of now consider Soft to be not fork save.
If you want to use Soft in a multiprocess environment initiallize Soft in each child process individually.

