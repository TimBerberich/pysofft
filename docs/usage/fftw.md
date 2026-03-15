# FFTW
PySOFFT uses [FFTW](https://www.fftw.org){target="_blank"} for the computation of normal fast Fourier transforms.
## Flags
In the creation of a transform instance you can specify [FFTW flags](https://www.fftw.org/fftw3_doc/Planner-Flags.html#Planner-Flags){target="_blank"} that are used in the creation of fft plans.
Many of these flags are acessible at
```py
from pysofft import fftw
fftw.flags
```
There main use is to tell the planner how thouroughly it should test different evaluation strategies for the given fftw plan to find the fastest strategy for your hardware, e.g.:
```py
from pysofft import Soft,fftw
flags = fftw.flags
s = Soft(32,init_ffts = True,fftw_flags=flags.fftw_measure)
```
`fftw_measure` can bring significant speedups compared to the default `fftw_estimate`.

## Wisdom
[FFTW wisdom](https://www.fftw.org/fftw3_doc/Words-of-Wisdom_002dSaving-Plans.html#Words-of-Wisdom_002dSaving-Plans){target="_blank"} is the idea to save the optimal fft execution strategies to a file and reload them at a later time, since plan creation using `fftw_measure` or even `fftw_exhaustive` can take a lot of time.
You can make use of this feature as follows:
```py
from pysofft import Soft,fftw

s = Soft(32,
	use_fftw_wisdom = True,
	init_ffts=True,
	fftw_flags = fftw.flags.fftw_measure)
```
By default this stores fftw wisdom at `~/.config/pysofft/fftw_wisdom.dat`.  
You can change the location of the used wisdom file:
```py
from pysofft import fftw
fftw.set_wisdom_file_path('~/some/other/path/your_wisdom.dat')
```
This change is not percistent and only remains for the current pysofft import. 

/// note | multiprocessing
To avoid race conditions and possible file corruptions, saving wisdom is only permitted from the master process in python.
The program wont crash but wisdom is simply not saved if computed from a child process.
Loading wisdom is always possible.
///
