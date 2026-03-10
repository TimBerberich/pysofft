# Package structure

Upon calling 
```py
import pysofft
```
the following modules are available

| Name           | Description                                                                                                                      |
|----------------|----------------------------------------------------------------------------------------------------------------------------------|
| pysofft._soft  | Module containing all f2py wraped fortran modules                                                                                |
| pysofft.soft   | Module containing all SO(3) Fourier transform related custom coded python wrappers and extra routines.                           |
| pysofft.Soft   | Main transform class. Alias for pysofft.soft.Soft                                                                                |
| pysofft.utils  | Containing fortran indixing tools for Wigner coefficients. Alias for pysofft._soft.utils                                         |
| pysofft.wigner | Contains all routines to generate small wigner d matrices. Alias for pysofft._soft.wigner                                        |
| pysofft.fftw   | Contains fucntions to manage FFTW wisdom related settings. (Convenience module to access specific functions from pysofft._soft)  |
| pysofft.omp    | Contains functions to set and get number of openMP threads. (Convenience module to access specific functions from pysofft._soft) |
| pysofft.stats  | (NOT FULLY IMPLEMENTED) Python module containing Classes targeting statistical analysis over SO(3).                              |



