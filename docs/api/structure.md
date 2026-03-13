# Package structure

Upon calling 
```py
import pysofft
```
the following modules are available

| Name                                               | Description                                                                                                |
|----------------------------------------------------|------------------------------------------------------------------------------------------------------------|
| pysofft._soft                                      | Module containing all f2py wraped fortran modules                                                          |
| [`pysofft.soft`][pysofft.soft]                     | Module containing all SO(3) Fourier transform related custom coded python wrappers and extra routines.     |
| [`pysofft.Soft`][pysofft.soft.Soft]                | Main transform class. Alias for pysofft.soft.Soft                                                          |
| [`pysofft.utils`](/pysofft/namespaceutils)         | Containing fortran indixing tools for Wigner coefficients. Alias for pysofft._soft.utils                   |
| [`pysofft.wigner`](/pysofft/namespacemake__wigner) | Contains all routines to generate small wigner d matrices. Alias for pysofft._soft.wigner                  |
| [`pysofft.fftw`][pysofft.fftw]                     | Contains fucntions to manage FFTW wisdom related settings. (Access specific functions from pysofft._soft)  |
| [`pysofft.omp`][pysofft.omp]                       | Contains functions to set and get number of openMP threads. (Access specific functions from pysofft._soft) |
| [`pysofft.stats`][pysofft.stats]                   | (NOT FULLY IMPLEMENTED) Python module containing Classes targeting statistical analysis over SO(3).        |



