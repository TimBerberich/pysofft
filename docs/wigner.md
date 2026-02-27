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
	The basic idea given in [Kostelecs paper](https://doi.org/10.1007/s00041-008-9013-5){target="_blank"}([pdf](https://www.cs.jhu.edu/~misha/ReadingSeminar/Papers/Kostelec08.pdf){target="_blank"}) is to use known values for $d^n_{m,n}(\beta)$ to start a three term recurrence in $l$ that allows to compute Wigner-d at $l+1$ from known vales at $l-1$ and $l$, i.e.:

	$$
	\{d^{l-1}_{mn}(\beta),d^{l}_{mn}(\beta)\} \longrightarrow d^{l+1}_{mn}(\beta)
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
	![dnmn](images/dlml_simple.svg#only-light){width="300"}
	![dnmn](images/dlml_simple_white.svg#only-dark){width="300"}
	<figcaption style="width:1000px">Example recursion scheme for $0\leq m\leq n<6$. Diagonal arrows correspond to usage of the first recursion relation and horizontal arrows imply the application of the second recursion relation. The shown pattern continues for higher $n,m$. </figcaption>
</figure>

This recursion is implemented in [`compute_all_dlml_l_contiguous`](../pysofft/namespacemake__wigner/#function-compute_all_dlml_l_contiguous).

/// info | Numerical stability for m\leq n<5000
As the following graph shows there are now underflows acc
///

Based on these vales we can use a three-term recurrence formulat to sucessively increas the order and compute all remaining values for 

for $0\leq m \leq n$, the following recurrence is used to compute all $d^l_{mn}$ for $m \leq n \leq l<bw$.

$$\begin{aligned}
 d^{l+1}_{mn}(\beta) = \sqrt{\frac{2l+3}{2l+1}} \frac{(l+1)(2l+1)}{\sqrt{\left[(l+1)^2-m^2\right]\left[(l+1)^2-n^2\right]}} \left(\cos\beta -\frac{mn}{l(l+1)}\right)\ &d^l_{mn}(\beta) \\
	                 -\quad \sqrt{\frac{2l+3}{2l-1}} \frac{\sqrt{[l^2-m^2][l^2-n^2]}}{\sqrt{\left[(l+1)^2-m^2\right]\left[(l+1)^2-n^2\right]}}\frac{l+1}{l}\ &d^{l-1}_{mn}(\beta)					 
\end{aligned} $$

All other Wigner-d matrix elements are connected to these values by variouse symmetry 



### Risbo recurrence
## Accuracy
