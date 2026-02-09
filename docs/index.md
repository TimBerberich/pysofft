# Fast Fourier transforms on the 3D rotation group $\mathrm{SO}(3)$
PySOFFT provides high performance routines for harmonic analysis on the 3D rotation group written in Fortran wrapped in Python.  

$$
f:\mathrm{SO}(3) \rightarrow \mathbb{C} \quad \overset{\mathrm{PySOFFT}}{\longleftrightarrow} \quad f^l_{m,n} \text{ with } |m|,|n|\leq l<\infty
$$

## Some applications:
* Statistical analysis over SO(3).
* Rotational alignment of datasets given by their spherical harmonic coefficients.
* X-ray scattering simulations of randomly oriented particles.

## Main features:
* __Fast__
* __OpenMP__ routines to speed up single transforms or compute many in parallel. 
* Dedicated faster transforms for real data.
* Compute __rotational cross-correlations__.
* Built in Python wrapper.
* On-the-fly computation of Wigner matrices: saving memory.

## Origin
PySOFFT started as a partial python port of soft-2.0 released by Peter Kostelec and Daniel Rockmore.
Now it is entire rewrite from scratch in Fortran with several improvements, e.g better precission in Wigner-D matrix computations of high orders as well as a new data layout for harmonic coefficients.

For details see their paper:
FFTs on the Rotation Group
J Fourier Anal Appl (2008) 14: 145–179
DOI 10.1007/s00041-008-9013-5

PySOFT is made available with conset of the original soft-2.0 authors and under the same GPL3 license.

## Installation
The easiest installation option is via pip.
	
	pip install pysofft
	
Non-python dependencies are __fftw__ and __gfortran__.
The only python dependency is __numpy__.

Compiling the fortran code only can be done via

## Basic Usage Python
	
Forward and inverse transforms

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
	
Accessing individual harmonic coefficients

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
	
	
	
	
