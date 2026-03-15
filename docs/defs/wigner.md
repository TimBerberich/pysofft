# Wigner-D matrix

## Definition
In this package we use the following definition for Wigner D-matrices $D^l_{m,n}(\alpha,\beta,\gamma)$ as function of euler angles $\alpha,\beta,\gamma$, using the $ZYZ$ convention.

$$
D^l_{m,n}(\alpha,\beta,\gamma) = e^{-im\alpha} d^l_{m,n}(\beta) e^{-in\gamma}
$$

$$
d^l_{m,n}(\beta) = A_l(-1)^{(m-n+\mu)/2}\sqrt{\frac{s!(s+\mu+\nu)!}{(s+\mu)!(s+\nu)!}} \left(\sin\frac{\beta}{2}\right)^\mu \left(\cos\frac{\beta}{2}\right)^\nu P^{\mu,\nu}_s(\cos\beta)
$$

where $\mu=|m-n|,\ \nu=|m+n|,\ s=l-\frac{\mu+\nu}{2},\quad  P^{\mu,\nu}_s$ are [Jacobi polynomials](https://en.wikipedia.org/wiki/Jacobi_polynomials){target="_blank"} and $A_l=\sqrt{\frac{2l+1}{2}}$ is an optional normalization term.

## Symmetries
Not all of the values of $d^l_{m,n}(\beta)$ are independent. PySOFFT uses the following symmeries to reduce the number of elements to compute by a factor of 8. 

$$
\begin{aligned}
d^l_{m,n}(\beta)&= d^l_{-n,-m}(\beta)\\
d^l_{m,n}(\beta)&=(-1)^{m-n}d^l_{n,m}(\beta)\\
d^l_{m,n}(\beta)&=(-1)^{l+2m-n}d^l_{m,-n}(\pi-\beta)
\end{aligned}
$$

In order to take advantage of the last symmetry we have to use a $\pi$ symmetric sampling of $\beta$. In that case we only need to compute the values of $d^l_{m,n}(\beta)$ for $0\leq m \leq n \leq l < \mathrm{bw}$ for a given bandwith $\mathrm{bw}$.

## Computation via Recurrence
PySOFFT implements two different recursion based computation routines for Wigner small-d matrices.
### Kostelec recurrence
!!! note 
	The basic idea given in [Kostelecs paper](https://doi.org/10.1007/s00041-008-9013-5){target="_blank"}([pdf](https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Kostelec08.pdf){target="_blank"}) is to use known values for $d^n_{m,n}(\beta)$ to start a three term recurrence in $l$ that allows to compute Wigner $d^l_{m,n}$ at $l+1$ from known values at $l-1$ and $l$, i.e.:

	$$
	\{d^{l-1}_{m,n}(\beta),d^{l}_{m,n}(\beta)\} \longrightarrow d^{l+1}_{m,n}(\beta)
	$$

The recursion start is given by the equation

$$
d^n_{m,n}(\beta) = \sqrt{\frac{2n+1}{2}} \sqrt{\frac{2n!}{(n+m)!(n-m)!}} \cos\left(\frac{\beta}{2}\right)^{n+m} \sin\left(\frac{\beta}{2}\right)^{n-m} 
$$

To compute this relation PySOFFT employs the following recursion formulas

$$
\begin{aligned}
n,m-1 &\rightarrow n+1,m & {d^{n+1}_{m,n+1}}(\beta)&={d^n_{m-1,n}}\sqrt{\frac{(2n+2)(2n+3)}{(n+m)(n+m+1)}} \cos(\frac{\beta}{2})^2 \\
n,m &\rightarrow n+1,m & d^{n+1}_{m,n+1}(\beta)&=d^l_{m,n}\sqrt{\frac{(2n+2)(2n+3)}{(n+m+1)(n-m+2)}} \cos(\frac{\beta}{2}) \sin(\frac{\beta}{2})
\end{aligned}
$$

Starting from $d^{0}_{0,0}(\beta)=\sqrt{\frac{1}{2}}$ (i.e. the normed definition of $d^l_{m,n}$ these two recursion formulas are applied according to the following scheme 

<figure markdown="span">
	![dnmn](/images/dlml_simple.svg#only-light){width="300"}
	![dnmn](/images/dlml_simple_white.svg#only-dark){width="300"}
	<figcaption style="width:1000px">Example recursion scheme for $0\leq m\leq n<6$. Diagonal arrows correspond to usage of the first recursion relation and horizontal arrows imply the application of the second recursion relation. The shown pattern continues for higher $n,m$. </figcaption>
</figure>

This recursion is implemented in [`compute_all_dlml_l_contiguous`](/pysofft/namespacemake__wigner/#function-compute_all_dlml_l_contiguous).

Based on these vales we can use a three-term recurrence formulat to sucessively increas the order and compute all remaining values for 

for $0\leq m \leq n$, the following recurrence is used to compute all $d^l_{mn}$ for $m \leq n \leq l<bw$.

$$\begin{aligned}
 d^{l+1}_{mn}(\beta) = \sqrt{\frac{2l+3}{2l+1}} \frac{(l+1)(2l+1)}{\sqrt{\left[(l+1)^2-m^2\right]\left[(l+1)^2-n^2\right]}} \left(\cos\beta -\frac{mn}{l(l+1)}\right)\ &d^l_{mn}(\beta) \\
	                 -\quad \sqrt{\frac{2l+3}{2l-1}} \frac{\sqrt{[l^2-m^2][l^2-n^2]}}{\sqrt{\left[(l+1)^2-m^2\right]\left[(l+1)^2-n^2\right]}}\frac{l+1}{l}\ &d^{l-1}_{mn}(\beta)					 
\end{aligned} $$

All other Wigner-d matrix elements are connected to these values by [symmetries](wigner.md#symmetries).  

Kostelec recurrence is currently the default option and is used when you instanciate a transform object as follwos
```py
from pysofft import Soft
bw = 32
s = Soft(bw)
# or equivalently
s = Soft(bw,recurrence_type=Soft.recurrence_types.kostelec)
```

/// Info | Advantages of Kostelec recurrence in PySOFT
__Speed__  
This recurrence does not mix different $m,n$ values it this allows to compute the wigner transforms conceptually as follows
```Fortran
for m=0,bw-1
	for n=m, bw-1
		!compute all relevant wigners
		d_mn = ...
		!Perform wigner transform
		flmn = matmul(f,d_mn)
```
The important thing here is that one can use highly efficient matrix multiplication to compute the inner sum over $l$. 
For the same reason the Kostelec recurrence behaves much better when using OpenMP to accellerate computations.
For details take a look at the [performance page](/speed).
///

/// Example | compute $d^l_{m,n}$ for fixed $l$. 
If you want to compute a full Wigner small-d matrix at fixed $l$ for some $\beta$ you can do so via
```py
from pysofft import wigner as wig
import numpy as np
l=32
betas = np.array([np.pi*4/11,np.pi-np.pi*4/11])
normalization = False
assert np.allclose(np.pi-betas,betas[::-1]), 'np.pi -betas musst be the same as betas[::1] because of the used symmetries.'

dlmn = wig.wigner_dl_kostelec(l,betas,normalization)

print(f'{dlmn.shape} (betas,m,n)')
# Check that dlmn is Orthogonal as it should be
np.allclose(dlmn[0]@dlmn[0].T,np.eye(2*l+1))
```
Note that to compute $d^l_{m,n}$ this routinte has to also compute $d^i_{m,n}$ for all $i<l$, so it is quite inefficient if you want to compute multiple orders at once. If you want to do this, thake a look at [`genwig_all_kostelec`](/pysofft/namespacemake__wigner/#function-genwig_all_kostelec) or the implementation of [`wigner_dl_kostelec`](/pysofft/namespacemake__wigner/#function-wigner_dl_kostelec) together with [`wig_l_recurrence_kostelec`](/pysofft/namespacemake__wigner/#function-wig_l_recurrence_kostelec).

///

/// Warning | Stability at l < 600
See [accuracy](/accuracy) for details.
For certain $\beta$ values and $l > 600$ the recurrence seems to brake down.
I currently dont fully understand why, my best guess is that the recursion start $d^n_{m,n}(\beta)$
has true values lower than double precission ($< \sim 10^{-308}$).
///


### Risbo recurrence
The recurrence givnen in [Risbos paper](https://doi.org/10.1007/bf01090814){target="_blank"} starts from $d^1_{m,n}$ for all $m,n$:

$$\begin{bmatrix}
	\cos(\frac{\beta}{2}) & -\sqrt{3}\cos(\frac{\beta}{2})\sin(\frac{\beta}{2}) & \sin(\frac{\beta}{2})^2\\
	\sqrt{3}\cos(\frac{\beta}{2})\sin(\frac{\beta}{2}) & \cos(\frac{\beta}{2})-\sin(\frac{\beta}{2})^2 & -\sqrt{3}\cos(\frac{\beta}{2})\sin(\frac{\beta}{2})\\
	\sin(\frac{\beta}{2})^2 & \sqrt{3}\cos(\frac{\beta}{2})\sin(\frac{\beta}{2}) & \cos(\frac{\beta}{2})
\end{bmatrix} 	 $$

and then increases the degree $l-\frac{1}{2} \rightarrow l$ via

$$
\begin{aligned}
	d^{l+\frac{1}{2}}_{m,n}(\beta) &=  d^{l-\frac{1}{2}}_{m-1,n-1}(\beta) \frac{\sqrt{(l+m)(l+n))}}{2l}\cos(\frac{\beta}{2}) \\
	                               &  -d^{l-\frac{1}{2}}_{m-1,n}(\beta) \frac{\sqrt{(l+m)(l-n)}}{2l}\sin(\frac{\beta}{2}) \\ 
								   &  +d^{l-\frac{1}{2}}_{m,n-1}(\beta) \frac{\sqrt{(l-m)(l+n)}}{2l}\sin(\frac{\beta}{2}) \\ 
								   &  +d^{l-\frac{1}{2}}_{m,n}(\beta) \frac{\sqrt{(l-m)(l-n)}}{2l}\cos(\frac{\beta}{2})
\end{aligned}
$$

graphically this can be represented as follows:

![dnmn](/images/risbo_example_sym.svg#only-light){width="420"}
![dnmn](/images/risbo_example_sym_white.svg#only-dark){width="420"}
![dnmn](/images/risbo_kernel.svg#only-light){width="365"}
![dnmn](/images/risbo_kernel_white.svg#only-dark){width="365"}

On the left you can see the how the recurrence looks for the example of $d^{1}_{m,n} \rightarrow d^{\frac{3}{2}}_{m,n}$.
On the right is a graphical representation of the recurrence kernel that has been given in the equation from above.

Risbo recurrence can be used when you instanciate a transform object as follwos
```py
from pysofft import Soft
bw = 32
s = Soft(bw,recurrence_type=Soft.recurrence_types.risbo)
```
/// Info | Advantages of Risbo recurrence in PySOFT
__Precission__  
The Risbo recurrence has higher acurracy at high $l$ and does not brake down for $l>600$, for details see [accuracy](/accuracy).
///

///  Example | compute $d^l_{m,n}$ for fixed $l$. 
If you want to compute a full Wigner small-d matrix at fixed $l$ for some $\beta$ you can do so via
```py
from pysofft import wigner as wig
import numpy as np
l=32
betas = np.array([np.pi*4/11,np.pi-np.pi*4/11])
normalization = False
assert np.allclose(np.pi-betas,betas[::-1]), 'np.pi -betas musst be the same as betas[::1] because of the used symmetries.'

dlmn_red = wig.wigner_dl_risbo_reduced(l,betas,normalization)
dlmn = wig.sym_reduced_to_full_wigner(dlmn,l)

print(f'{dlmn.shape} (betas,m,n)')
# Check that dlmn is Orthogonal as it should be
np.allclose(dlmn[0]@dlmn[0].T,np.eye(2*l+1))
```
Note that to compute $d^l_{m,n}$ this routinte has to also compute $d^i_{m,n}$ for all $i<l$, so it is quite inefficient if you want to compute multiple orders at once. If you want to do this, thake a look at [`genwig_all_risbo`](/pysofft/namespacemake__wigner/#function-genwig_all_risbo) the implementation of [`wigner_dl_risbo_reduced`](/pysofft/namespacemake__wigner/#function-wigner_dl_risbo_reduced) and [`wigner_recurrence_risbo_reduced`](/pysofft/namespacemake__wigner/#function-wigner_recurrence_risbo_reduced).
///
