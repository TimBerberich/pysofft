# Many transforms
PySOFFT allows you to perform many independent transforms in one go.
This can be usefull especially in combination with [parallelization](parallelization.md).
The general api is similar to the single transfoms shown in [simple transforms](simple.md)

```py
from pysofft import Soft
s = Soft(32)

# create 100 coefficient arrays and fill them with random numbers.
many_flmn = s.get_coeff(howmany=100,random=True)

# Perform 100 inverse SO(3) Fourier transforms.
many_f = s.isoft_many(many_flmn)
print(many_f.shape)
```
