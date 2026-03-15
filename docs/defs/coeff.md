# Wigner Coefficients

By Wigner coefficients of a function $f$  
 
$$ f:\mathrm{SO}(3) \rightarrow \mathbb{K} \quad \text{with } \mathbb{K}=\mathbb{C} \text{ or } \mathbb{R},$$

we mean the result of the following integral for all $l \in \mathbb{N}$ and $m,n \in \mathbb{Z}$ with $|m|,|n|\leq l$

$$ f^l_{m,n} = \int_{\mathrm{SO}(3)} dg\ f(g) D^l_{m,n}(g)^*, $$

where $D^l_{m,n}(g)$ are the [Wigner-D matrices](wigner.md).
PySOFFT computes a discrete version of the above integral which uses a finite sampling of $SO(3)$ (see [euler angles](so3_grid.md)) as well as an upper limit of $l$ given by the transform bandwidth \mathrm{bw}, i.e. $l<\mathrm{bw}$.

The total number of coefficients of a transform with bandwith $\mathrm{bw}$ is

$$ \sum_{l=0}^{\mathrm{bw}-1} (2l+1)^2 = (4\mathrm{bw}^3-\mathrm{bw})/3 $$

## Coefficient access
The easiest way to access Wigner coefficients is via the provided ndarray subclass [CoeffSO3](:) in this way you dont have to deal with the internal data layout at all.
```py
from pysofft import Soft
s = Soft(32)
flmn = s.get_coeff()
l=10
m=-2
n=8
c = flmn.lmn[l,m,n]
print(f"f^{l}_"+"{"+f"{m},{n}"+"}"+f" = {c}")

# You can also assign data
flmn.lmn[l,m,n] = 42
print(f"f^{l}_"+"{"+f"{m},{n}"+"}"+f" = {c}")

# As well as use slicing
cs = flmn.lmn[8:16,m,n]
print("f^{8...16}_"+"{"+f"{m},{n}"+"}"+f" = {cs}")
```

## Data layout
Internally PySOFT uses two different layouts depending on the [recurrence type](wigner.md#computation-via-recurrence) that is used to generate the Wigner small-d matrices. If [Kostelec recurrence](wigner.md#kostelec-recurrence) is selected (default) or the wigner small-d matrices have been precomputed, the data layout is such that the following loops iterate contiguousely over the 1D coefficient array.

/// Example | "mnl" layout (Kostelec recurrence)
```py
from pysofft import Soft
import numpy as np
s = Soft(3)
flmn = s.get_coeff()
print(f"The Wigner coeff order is: {flmn.coeff_order}")
flmn[:] = np.arange(len(flmn))

for m in range(s.bw):
    for n in range(m,s.bw):
        print(flmn.lmn[:,m,n])

        if m==0 and n==0: continue

        print(flmn.lmn[:,-n,-m])

        if m!=n:
            print(flmn.lmn[:,n,m])
            print(flmn.lmn[:,-m,-n])

        if m==0: continue
        print(flmn.lmn[:,m,-n])
        print(flmn.lmn[:,-m,n])

        if m==n: continue
        print(flmn.lmn[:,n,-m])
        print(flmn.lmn[:,-n,m])
```
///

If [Risbo recurrence](wigner.md#risbo-recurrence) is used and the wigner small-d matrices have __not__ been precomputedthe, the following loops are contigous instead.
/// Example | "lmn" layout (Risbo recurrence)
```py
from pysofft import Soft
import numpy as np
s = Soft(3,recurrence_type = s.recurrence_types.risbo)
flmn = s.get_coeff()
print(f"The Wigner coeff order is: {flmn.coeff_order}")
flmn[:] = np.arange(len(flmn))

_id = 0
for l in range(s.bw):
    for m in range(l+1):
        for n in range(m,l+1):
            
            print(flmn.lmn[l,m,n])
  
            if m==0 and n==0: continue
            print(flmn.lmn[l,-n,-m])
            
            if m!=n:
                print(flmn.lmn[l,n,m])
                print(flmn.lmn[l,-m,-n])

            if m==0: continue
            print(flmn.lmn[l,m,-n])
            print(flmn.lmn[l,-m,n])

            if m==n: continue
            print(flmn.lmn[l,n,-m])
            print(flmn.lmn[l,-n,m])
```
///

## Wigner coefficients of real functions
The Wigner coefficents of a real valued function $f$ satisfy the following symmetry relation

$$ {f^l_{m,n}}^* = f^l_{-m,-n} (-1)^{m-n} $$

Currenctly PySOFFT does not have a special Wigner coefficient type that reduces the number of stored coefficients using the above symmetry.
Instead the real transforms ([rsoft][pysofft.soft.Soft.rsoft] and [irsoft][pysofft.soft.Soft.irsoft]) will only access half of the values in the coefficient array.

If you want to manually populate wigner coefficients for real functions, write only to indices $l,m,n$ where $m$ is positive and afterwards call [enforce_real_symmetry][pysofft.soft.Soft.enforce_real_symmetry], i.e.
```py
from pysofft import Soft
import numpy as np
s = Soft(16)
flmn = s.get_coeff()

for m in range(s.bw):
    for n in range(-s.bw+1,s.bw):
        flmn.lmn[:,m,n] = np.random.rand(len(flmn.lmn[:,m,n]))
s.enforce_real_symmetry(flmn)

f2lmn = s.rsoft(s.irsoft(flmn))
np.allclose(f2lmn,flmn)
```









