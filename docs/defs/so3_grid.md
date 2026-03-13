# Euler angles 

PySOFFT uses Euler-angles in the ZYZ-convention as coordinates for rotation group $\mathrm{SO}(3)$.

/// warning | Euler angles are disgusting!
To the true connoisseurs among you, especially those coming over from [spherical.readthedocs.io](https://spherical.readthedocs.io/){target="_blank"}.  
I know, I share your pain ... euler angles are utterly disgusting and indeed quaternions,matrices ,or the axis-angle form are way better. Buuuut Euler angles bribed me ... with FFTs ... and I could not resist, for now at least. 
///

![Euler angles](/images/eulerangles.png){ width="400" align=left } Using the corotation reference picture a rotation $(\alpha,\beta,\gamma)\in\mathrm{SO}(3)$ acts via:

1. Rotate by $\alpha$ around $\boldsymbol{z}$, which results in a new coordinate frame $(\boldsymbol{x'} ,\boldsymbol{y'} , \boldsymbol{z})$.
2. Rotate by $\beta$ around $\boldsymbol{y'}$, which results in a new coordinate frame $(\boldsymbol{x''} ,\boldsymbol{y'} , \boldsymbol{z'})$.
3. Rotate by $\gamma$ around $\boldsymbol{z'}$.

Resulting in the final displayed coordinate axes $(\boldsymbol{x_{\mathrm{rot}}} ,\boldsymbol{y_{\mathrm{rot}}} , \boldsymbol{z_{\mathrm{rot}}})$.
<br style="clear: both;" />

PySOFFT uses the following euler angle sampling for a given bandwith `bw`.

$$
\begin{aligned}
	\alpha_i,\gamma_i &= \frac{i \pi}{\mathrm{bw}} \quad \text{for } i\in\{0,\ldots,2\mathrm{bw}-1\}   & 
	\beta_i &= \frac{(i+\frac{1}{2})\pi}{2\mathrm{bw}} \quad \text{for } i\in\{0,\ldots,2\mathrm{bw}-1\} 
\end{aligned}
$$

These angles can be acessed form a transform object via
```Python
from pysofft import Soft
s = Soft(8)

a =  s.euler_angles['alpha']
g =  s.euler_angles['gamma']
b =  s.euler_angles['beta']
```
