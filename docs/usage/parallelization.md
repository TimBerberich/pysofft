# Parallelization

## Native
PySOFT uses openMP in fortran to implement multiprocessing.  
You can set and get the number of used threads via:
```py
from pysofft import omp

print(f'Current max number of trheads = {omp.get_max_threads()}')
omp.set_num_threads(2)
print(f'New max number of trheads = {omp.get_max_threads()}')
```

All transforms have an optional boolean agument `use_mp` by which you can enable parallelization.  
For multiple transform routines, i.e.`s.soft_many,s.isoft_many,s.rsoft_many,s.irsoft_many`, parallelizes over the full individual transforms and the speedup scales almost linear with the number of threads as long as no hardware limits are hit, see [performace_metrics](../speed.md#multi-core-performance).
```py
from pysofft import Soft,omp
s = Soft(32)
many_flmn = s.get_coeff(random=True,howmany=10)
many_f = s.isoft_many(many_flmn,use_mp=True)
```
For single transforms, i.e.`s.soft,s.isoft,s.rsoft,s.irsoft`, parallel computation speeds up the Wigner transform part but not the computation of ordinary ffts. The achivable speedup is limited to about 2-3 times compared to single threaded transforms, see [performance metrics](../speed.md#multi-core-performance).
```py
from pysofft import Soft,omp
s = Soft(32)
flmn = s.get_coeff(random=True)
f = s.isoft(flmn,use_mp=True)
```

## Python
All routines in PySOFFT are fork save. 
Although the native routines are faster, there is no problem with writing code like the following.
```py
from multiprocessing import Pool
import numpy as np
from pysofft import Soft

def forked_soft(args):
    s,flmn = args
    f = s.isoft(flmn,use_mp=False)
    return f

s = Soft(16)
flmn = s.get_coeff(howmany=100,random=True)

with Pool(8) as p:
    f2 = p.map(forked_soft,[(s,c) for c in flmn])
    f2 = np.array(f2)

f = s.isoft_many(flmn,use_mp=True)
print(f'Pool computation gives same result as s.isoft_many = {np.allclose(f2,f)}')
```
