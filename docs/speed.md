# Performance metrics
The following performence metric where generated on a compute node with two AMD EPYC 7543, a total of 64 physical CPU cores and 512GB of RAM. 
Results are shown for both [Kostelec](defs/wigner.md#kostelec-recurrence) and [Risbo](defs/wigner.md#risbo-recurrence) recurrence shemes. 

/// Info | on-the-fly computation VS. precomputed Wigners
All results on this site where generated using the on-the-fly computation of Wigner matrices.  
The precomputed variants give a slight speed increase and eliminate all performance differences between Kostelec and Risbo reccurences.
Memory consumption becomes a nightmare though. The number of Wigner coefficients scales with $O(\mathrm{bw}^4)$ and their required memory reaches $~180\,$GB at `bw = 512`.
///

<p float="left" align="middle">
  <img src="single_core_speed_kostelec.png" width="490" />
  <img src="single_core_speed_risbo.png" width="490" /> 
</p>

__Kostelec recurrence__

| method\\bw | 16                    | 32                      | 64                     | 128                     | 256                    | 512                 |
|------------|-----------------------|-------------------------|------------------------|-------------------------|------------------------|---------------------|
| soft       | $450\mu$s $\pm 1\mu$s | $4.91$ms $\pm 0.05$ms   | $61.6$ms $\pm 0.3$ms   | $0.8466$s $\pm 0.0006$s | $12.0$s $\pm 0.4$s     | $172$s $\pm 4$s     |
| isoft      | $501\mu$s $\pm 4\mu$s | $5.28$ms $\pm 0.02$ms   | $64.6$ms $\pm 0.2$ms   | $0.8962$s $\pm 0.0006$s | $13.328$s $\pm 0.004$s | $186$s $\pm 1$s     |
| rsoft      | $273\mu$s $\pm 2\mu$s | $2.917$ms $\pm 0.004$ms | $35.5$ms $\pm 0.2$ms   | $0.5517$s $\pm 0.0001$s | $7.8$s $\pm 0.2$s      | $103.7$s $\pm 0.7$s |
| irsoft     | $273\mu$s $\pm 2\mu$s | $3.01$ms $\pm 0.02$ms   | $38.1.6$ms $\pm 0.2$ms | $0.5783$s $\pm 0.0002$s | $9.00$s $\pm 0.07$s    | $120$s $\pm 1$s     |

__Risbo recurrence__

| method\\bw | 16                    | 32                    | 64                   | 128                   | 256                  | 512                 |
|------------|-----------------------|-----------------------|----------------------|-----------------------|----------------------|---------------------|
| soft       | $401\mu$s $\pm 2\mu$s | $5.35$ms $\pm 0.04$ms | $81.0$ms $\pm 0.3$ms | $1.457$s $\pm 0.009$s | $20.46$s $\pm 0.02$s | $172$s $\pm 4$s     |
| isoft      | $353\mu$s $\pm 1\mu$s | $4.00$ms $\pm 0.03$ms | $60$ms $\pm 1$ms     | $1.197$s $\pm 0.009$s | $17.11$s $\pm 0.02$s | $186$s $\pm 1$s     |
| rsoft      | $238\mu$s $\pm 1\mu$s | $3.21$ms $\pm 0.03$ms | $44.7$ms $\pm 0.2$ms | $0.878$s $\pm 0.001$s | $14.18$s $\pm 0.02$s | $103.7$s $\pm 0.7$s |
| irsoft     | $133\mu$s $\pm 1\mu$s | $2.34$ms $\pm 0.03$ms | $31.0$ms $\pm 0.2$ms | $0.720$s $\pm 0.001$s | $12.19$s $\pm 0.06$s | $120$s $\pm 1$s     |
	
The given errors are simply the standard deviations of the datasets computed by the code below.
/// note | Code: Single-core speed test 
WARNING: The bw=512 computation uses $\sim 70$GB of RAM and bw=256 uses $\sim 8$GB. Depending on your PC you might need to skip these.
``` py
from pysofft import Soft
import pysofft
import timeit
import numpy as np

soft_execution_times=[[],[]]
rsoft_execution_times=[[],[]]
isoft_execution_times=[[],[]]
irsoft_execution_times=[[],[]]
for bw in [16,32,64,128,256,512]:
    for recurrence_type in [0,1]:
        # instanciating Soft class + generating in/out arrays
        s = Soft(bw,init_ffts=True,
                 use_fftw_wisdom=True,
                 fftw_flags=pysofft._soft.softclass.fftw_measure,
                 precompute_wigners=False,
                 recurrence_type=recurrence_type)
        coeff = s.get_coeff(random=True)
        coeffr = s.get_coeff(random=True,real=True)
        func = s.get_so3func()
        funcr = s.get_so3func(real=True)
        if recurrence_type==0:
            recurrence = 'Kostelec'
        else:
            recurrence = 'Risbo'
        # speed tests
        tmpi = np.array(timeit.repeat('s.isoft(coeff,out = func)',number=10,repeat=7,globals=globals()))/10
        isoft_execution_times[recurrence_type].append(tmpi)
        tmpir = np.array(timeit.repeat('s.irsoft(coeffr, out = funcr)',number=10,repeat=7,globals=globals()))/10
        irsoft_execution_times[recurrence_type].append(tmpir)
        tmp = np.array(timeit.repeat('s.soft(func,out = coeff)',number=10,repeat=7,globals=globals()))/10
        soft_execution_times[recurrence_type].append(tmp)
        tmpr = np.array(timeit.repeat('s.rsoft(funcr, out = coeffr)',number=10,repeat=7,globals=globals()))/10
        rsoft_execution_times[recurrence_type].append(tmpr)
        np.save(f"soft_{recurrence}.npy",soft_execution_times[recurrence_type])
        np.save(f"rsoft_{recurrence}.npy",rsoft_execution_times[recurrence_type])
        np.save(f"isoft_{recurrence}.npy",isoft_execution_times[recurrence_type])
        np.save(f"irsoft_{recurrence}.npy",irsoft_execution_times[recurrence_type])
```
///

	
	
## Memory consumption
The memory consumption (in Byte) of individual transforms roughly scales with 

$$
  64 (2\ \mathrm{bw})^3
$$

Note that $(2\ \mathrm{bw})^3$ is the number of SO(3) grid points for a given bandwidth.
In paticular this is much lower than the number of Wigner-D matrix elements involvend in the trasform, which scales with
$O(\mathrm{bw}^4)$.


## Multi-core performance
PySOFFT allows parallelization in two different ways:  
 1. Compute transforms for many datasets in parallel.  
 2. Speed up individual transforms (Currently no parallel executions of the FFT part)  


<p float="left" align="middle">
  <img src="mp_many_kostelec.png" width="490" />
  <img src="mp_many_risbo.png" width="490" /> 
</p>

/// note | Code: Multi-core multi-transform speed test
``` py
import pysofft
from pysofft import Soft
import timeit
import numpy as np

soft_execution_times=[[],[]]
rsoft_execution_times=[[],[]]
isoft_execution_times=[[],[]]
irsoft_execution_times=[[],[]]
bw = 64
N = 256
for nthreads in [1,2,4,8,16,32,64]:
    for recurrence_type in [0,1]:
        # instanciating Soft class + generating in/out arrays
        s = Soft(bw,init_ffts=True,
                 use_fftw_wisdom=True,
                 fftw_flags=pysofft._soft.softclass.fftw_measure,
                 precompute_wigners=False,
                 recurrence_type=recurrence_type)
        
        pysofft.omp.set_num_threads(nthreads)
        coeff = s.get_coeff(random=True,howmany=N)
        coeffr = s.get_coeff(random=True,real=True,howmany=N)
        func = s.get_so3func(howmany=N)
        funcr = s.get_so3func(real=True,howmany=N)
        print(nthreads)
        if recurrence_type==0:
            recurrence = 'Kostelec'
        else:
            recurrence = 'Risbo'
        print('alive')
        tmp = np.array(timeit.repeat('s.isoft_many(coeff,out = func,use_mp=True)',number=10,repeat=1,globals=globals()))/10
        isoft_execution_times[recurrence_type].append(tmp)
        tmpr = np.array(timeit.repeat('s.irsoft_many(coeffr, out = funcr,use_mp=True)',number=10,repeat=1,globals=globals()))/10
        irsoft_execution_times[recurrence_type].append(tmpr)
        tmp = np.array(timeit.repeat('s.soft_many(func,out = coeff,use_mp=True)',number=10,repeat=1,globals=globals()))/10
        soft_execution_times[recurrence_type].append(tmp)
        tmpr = np.array(timeit.repeat('s.rsoft_many(funcr, out = coeffr,use_mp=True)',number=10,repeat=1,globals=globals()))/10
        rsoft_execution_times[recurrence_type].append(tmpr)
        
        np.save(f"soft_{recurrence}_mp_many.npy",soft_execution_times[recurrence_type])
        np.save(f"rsoft_{recurrence}_mp_many.npy",rsoft_execution_times[recurrence_type])
        np.save(f"isoft_{recurrence}_mp_many.npy",isoft_execution_times[recurrence_type])
        np.save(f"irsoft_{recurrence}_mp_many.npy",irsoft_execution_times[recurrence_type])
```
///

<p float="left" align="middle">
  <img src="mp_single_kostelec.png" width="490" />
  <img src="mp_single_risbo.png" width="490" /> 
</p>

/// note |Code: Multi-core single-transform speed test
``` py
import pysofft
from pysofft import Soft
import timeit
import numpy as np

soft_execution_times=[[],[]]
rsoft_execution_times=[[],[]]
isoft_execution_times=[[],[]]
irsoft_execution_times=[[],[]]
bw = 64
for nthreads in [1,2,4,8,16,32,64]:
    for recurrence_type in [0,1]:
        # instanciating Soft class + generating in/out arrays
        s = Soft(bw,init_ffts=True,
                 use_fftw_wisdom=True,
                 fftw_flags=pysofft._soft.softclass.fftw_measure,
                 precompute_wigners=False,
                 recurrence_type=recurrence_type)

        pysofft.omp.set_num_threads(nthreads)
        coeff = s.get_coeff(random=True)
        coeffr = s.get_coeff(random=True,real=True)
        func = s.get_so3func()
        funcr = s.get_so3func(real=True)
        print(nthreads)
        if recurrence_type==0:
            recurrence = 'Kostelec'
        else:
            recurrence = 'Risbo'
        print('alive')
        tmp = np.array(timeit.repeat('s.isoft(coeff,out = func,use_mp=True)',number=10,repeat=7,globals=globals()))/10
        isoft_execution_times[recurrence_type].append(tmp)
        tmpr = np.array(timeit.repeat('s.irsoft(coeffr, out = funcr,use_mp=True)',number=10,repeat=7,globals=globals()))/10
        irsoft_execution_times[recurrence_type].append(tmpr)
        tmp = np.array(timeit.repeat('s.soft(func,out = coeff,use_mp=True)',number=10,repeat=7,globals=globals()))/10
        soft_execution_times[recurrence_type].append(tmp)
        tmpr = np.array(timeit.repeat('s.rsoft(funcr, out = coeffr,use_mp=True)',number=10,repeat=7,globals=globals()))/10
        rsoft_execution_times[recurrence_type].append(tmpr)

	    np.save(f"soft_{recurrence}_mp_single.npy",soft_execution_times[recurrence_type])
        np.save(f"rsoft_{recurrence}_mp_single.npy",rsoft_execution_times[recurrence_type])
        np.save(f"isoft_{recurrence}_mp_single.npy",isoft_execution_times[recurrence_type])
        np.save(f"irsoft_{recurrence}_mp_single.npy",irsoft_execution_times[recurrence_type])
```
///
