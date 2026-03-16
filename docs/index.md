# Fast Fourier transforms on the 3D rotation group $\mathrm{SO}(3)$
PySOFFT provides high performance routines for harmonic analysis on the 3D rotation group written in Fortran and wrapped in Python.  

$$
f:\mathrm{SO}(3) \rightarrow \mathbb{C} \quad \overset{\mathrm{PySOFFT}}{\longleftrightarrow} \quad f^l_{m,n} \text{ with } |m|,|n|\leq l<\infty
$$

## Some applications:
* Rotational alignment of datasets given by their spherical harmonic coefficients.
* Statistical analysis over SO(3).
* X-ray scattering simulations of randomly oriented particles.

## Main features:
* [__Fast__](speed.md)
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
DOI [10.1007/s00041-008-9013-5](https://doi.org/10.1007/s00041-008-9013-5){target="_blank"} ([pdf](https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Kostelec08.pdf){target="_blank"})

PySOFT is made available with consent of the original soft-2.0 authors and under the same GPL3 license.

## Installation
PySOFFT has the following non-python dependencies
 [__fftw__](https://fftw.org/){target="_blank"},[__openmp__](https://www.openmp.org/){target="_blank"}, [__meson__](https://mesonbuild.com/){target="_blank"} and [__gcc__](https://gcc.gnu.org/){target="_blank"}.  
 The only python dependency is [__numpy__](https://numpy.org/){target="_blank"}.
### [Pixi](https://pixi.prefix.dev/latest/){target="_blank"}
 The advantage of pixi is that all non-python dependencies are installed automatically.
 You can use the following pixi.toml for installation.
 /// info | pixi.toml
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
 
 [pypi-dependencies]
 pysofft = ">=0.9.0, <2.0.0"
 ```
 ///
 after placing the above `pixi.toml` in an empty directory you can install the environment using
 
 ```bash
 pixi update
 ```
 and finally access a pixi shell using 
 
 ```bash
 pixi shell
 ```
 
 
### [PyPi](https://pypi.org/project/pysofft/){target="_blank"}
PySOFFT can be installed via 
```bash
pip install pysofft
```
NOTE: This install will fail if some of the non-python dependencies are not installed on your system.

### Git Clone + pip
You can also git clone and pip install.

```bash
git clone https://github.com/TimBerberich/pysofft
cd pysofft
pip install .
```
NOTE: This install will fail if some of the non-python dependencies are not installed on your system.  

Editable installs are possible by substituting the last line with
```bash
pip install --no-build-isolation -e .
```

/// Info | CPU specific optimization
This can significally speed up computations.  
You can use them by adding `-march=native` to the lines in 
`pysofft/pysofft/meson.build` starting with `c_args` and `fortran_args`, before calling `pip install .`.  
The relevant lines should then be
```
  fortran_args: ['-fopenmp','-lfftw3', '-O3', '-fno-math-errno', '-fno-trapping-math', '-march=native'],	
```

and 

```
  c_args: ['-fopenmp','-lfftw3', '-O3', '-fno-math-errno', '-fno-trapping-math', '-march=native'],	
```
///

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
	
	
	
