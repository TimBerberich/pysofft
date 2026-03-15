# Functions over SO(3)
PySOFFT samples functions over SO(3) at the Euler angles defined [here](so3_grid.md).
For a given bandwidth `bw` a function $f$ over SO(3) is given by a real or complex [numpy array](https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html){target="_blank"} of shape

$$ (2\mathrm{bw},2\mathrm{bw},2\mathrm{bw}) \quad \text{with axes} \quad (\gamma,\alpha,\beta)$$ 

/// info | Why $(\gamma,\alpha,\beta)$ and not  $(\alpha,\beta,\gamma)$ ?
In short, for performance reasons.  
The implemented SO(3) Fourier transform consists of 2 normal FFTs over the axes belonging to $\alpha$ and $\gamma$ followed by a Wigner transform over the $\beta$ axes.
For perfomance reason it is beneficial to compute the FFTs over adjecient axes.
Finially the C contiguouse array $(\gamma,\alpha,\beta)$ has the Fortran contiguous coordinates $(\beta,\alpha,\gamma)$ with 
$\beta$ beeing the innermost, i.e. fast, axis, over which the final Wigner transform is computed.
///
