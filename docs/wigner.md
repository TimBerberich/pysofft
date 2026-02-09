# Wigner-D matrix

## Definition
In this package we use the following definition for Wigner D-matrices $D^l_{m,n}(\alpha,\beta,\gamma)$ as function of euler angles $\alpha,\beta,\gamma$, using the $ZYZ$ convention.

??? warning "Euler angles are disgusting!"
	To the true connoisseurs among you, especially those coming over from [spherical.readthedocs.io](https://spherical.readthedocs.io/){target="_blank"}.  
	I know, I share your pain ... euler angles are utterly disgusting and indeed quaternions are way better.  
	Buuuut, they managed to bribe me with FFTs, for now that is. 

$$
D^l_{m,n}(\alpha,\beta,\gamma) = e^{-im\alpha} d^l_{m,n}(\beta) e^{-in\gamma}
$$

$$
d^l_{m,n}(\beta) = A_l(-1)^{(m-n+\mu)/2}\sqrt{\frac{s!(s+\mu+\nu)!}{(s+\mu)!(s+\nu)!}} \left(\sin\frac{\beta}{2}\right)^\mu \left(\cos\frac{\beta}{2}\right)^\nu P^{\mu,\nu}_s(\cos\beta)
$$

where $\mu=|m-n|,\ \nu=|m+n|,\ s=l-\frac{\mu+\nu}{2},\quad  P^{\mu,\nu}_s$ are [Jacobi polynomials](https://en.wikipedia.org/wiki/Jacobi_polynomials){target="_blank"} and $A_l=\sqrt{\frac{2l+1}{2}}$ is an optional normalization term.

## Computation via Recurrence
PySOFFT implements two different recursion based computation routines for Wigner small-d matrices.
### Kostelec recurrence
See []()
Starting from 

$$
d^n_{m,n}(\beta) = \sqrt{\frac{2n+1}{2}} \sqrt{\frac{2n!}{(n+m)!(n-m)!}} \cos(\frac{\beta}{2})^{n+m} \sin(\frac{\beta}{2})^{n-m} $$

for $0\leq m \leq n$, the following recurrence is used to compute all $d^l_{mn}$ for $m \leq n \leq l<bw$.

$$\begin{aligned}
 d^{l+1}_{mn}(\beta) = \sqrt{\frac{2l+3}{2l+1}} \frac{(l+1)(2l+1)}{\sqrt{\left[(l+1)^2-m^2\right]\left[(l+1)^2-n^2\right]}} \left(\cos\beta -\frac{mn}{l(l+1)}\right)\ &d^l_{mn}(\beta) \\
	                 -\quad \sqrt{\frac{2l+3}{2l-1}} \frac{\sqrt{[l^2-m^2][l^2-n^2]}}{\sqrt{\left[(l+1)^2-m^2\right]\left[(l+1)^2-n^2\right]}}\frac{l+1}{l}\ &d^{l-1}_{mn}(\beta)					 
\end{aligned} $$

All other Wigner-d matrix elements are connected to these values by variouse symmetry 



### Risbo recurrence
## Accuracy
