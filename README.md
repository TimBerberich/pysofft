[![Lates Release](https://img.shields.io/github/v/release/European-XFEL/xframe)](https://github.com/TimBerberich/pysofft/releases)
![License](https://img.shields.io/github/license/TimBerberich/pysofft)
![Language](https://img.shields.io/badge/language-python-blue)
![Language](https://img.shields.io/badge/language-fortran-blue)

# PySOFFT
## Fast Fourier transforms on the 3D rotation group $\mathrm{SO}(3)$
PySOFFT provides high performance routines for harmonic analysis on the 3D rotation group written in Fortran and wrapped in Python.  

$$
f:\mathrm{SO}(3) \rightarrow \mathbb{C} \quad \overset{\mathrm{PySOFFT}}{\longleftrightarrow} \quad f^l_{m,n} \text{ with } |m|,|n|\leq l<\infty
$$

## Some applications:
* Rotational alignment of datasets given by their spherical harmonic coefficients.
* Statistical analysis over SO(3).
* X-ray scattering simulations of randomly oriented particles.

## Main features:
* __Fast__ (see docs for more details)
* __OpenMP__ routines to speed up single transforms or compute many in parallel. 
* Dedicated faster transforms for real data.
* Compute __rotational cross-correlations__.
* Built in Python wrapper.
* On-the-fly computation of Wigner matrices: saving memory.

## Origin
PySOFFT started as a partial python port of soft-2.0 released by Peter Kostelec and Daniel Rockmore.
Now it is an entire rewrite in Fortran featuring several improvements, e.g higher precission computation of Wigner-D matrices, OpenMP routines and convenient python wrappers.

For details on soft-2.0 see:  
FFTs on the Rotation Group  
J Fourier Anal Appl (2008) 14: 145–179  
DOI [10.1007/s00041-008-9013-5](https://doi.org/10.1007/s00041-008-9013-5) ([pdf](https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Kostelec08.pdf))

PySOFT is made available with consent of the original soft-2.0 authors and under the same GPL3 license.

## Installation
The easiest installation option is via pip.
	
	pip install pysofft
	
The only python dependency is __numpy__.
Non-python dependencies are __fftw__, __openmp__, __meson__, __gcc__ and __gfortran__.

### Pixi
 If you use [pixi](https://pixi.prefix.dev/latest/) you can use the following pixi.toml for installation.
 ``` toml
  [workspace]
 channels = ["conda-forge"]
 description = "Workspace for pysofft"
 name = "temp"
 platforms = ["linux-64"]
 version = "0.1.0"
 
 [system-requirements]
 libc = { family = "glibc", version = "2.34" }
 
 [dependencies]
 gfortran = ">=15.2,<15.5.0"
 gcc = ">=15.2,<15.5.0"
 cpython = ">=3.13.11,<4"
 python = ">=3.13.11,<3.14"
 pip = ">=26.0.1,<27"
 numpy = ">=2.2.6,<2.2.7"
 meson = ">=1.8.1,<1.9.0"
 meson-python = ">=0.19.0,<0.20"
 fftw = ">=3.3.10,<4"
 openmp = ">=8.0.1,<9"
 ```
 
## Basic Usage Python
	
Forward and inverse transforms

```Python
from pysofft import Soft
bw = 64
s = Soft(bw)

# complex case: inverse then forward 
f_lmn = s.get_coeff(random=True)
f = s.isoft(f_lmn)
f_lmn2 = s.soft(f)

# real case: inverse then forward
g_lmn = s.get_coeff(real=True,random=True)
g = s.irsoft(g_lmn)
g_lmn2 = s.rsoft(g)
```
Accessing individual harmonic coefficients

```Python
from pysofft import Soft
bw = 64
s = Soft(bw)

f_lmn = s.get_coeff(random=True)

l=4
m=-1
n=2

# access single coefficient
print(f_lmn.lmn[l,m,n])

# slicing is also possible
print(f_lmn.lmn[l,m,:])

# as well as value asignment
f_lmn.lmn[l,m,:] = 1 + 1.j
```
