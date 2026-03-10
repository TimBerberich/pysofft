# Simple transforms

```Python
from pysofft import Soft
import numpy as np

# You can create a transform instance via.
bw = 32
s = Soft(bw)

# Get coefficient arrays and SO(3) function arrays 
f = s.get_so3func()
flmn = s.get_coeff()

# The euler angles on which f is sampled are given in
euler_anges = s.euler_angles

# Access to individual Wigner coefficients works as follows
# Note that indexes have to satisfy 0<=|m|,|n|<=l<bw
l=10
m=-10
n=2
flmn.lmn[l,m,n] = 1+1.j

# Slicing is also possible
print(f'flmn.lmn[10:12,m,n]     = {flmn.lmn[10:12,m,n]}') 
# There are some limits though e.g. indexing by 11:9:-1 does not invert the result order.
print(f'flmn.lmn[11:9:-1,m,n]   = {flmn.lmn[11:9:-1,m,n]}') 

#############################
# inverse complex transform #

# The following creates a new output array
g = s.isoft(flmn)
# This uses the already defined f array as output
s.isoft(flmn,out=f)
print(f'f and g are the same    = {np.allclose(f,g)}')

#############################
# forward complex transform #
glmn = s.soft(g)
print(f'Transform of g at l,m,n = {glmn.lmn[l,m,n]}')


# Real coefficient arrays, initiallized with random numbers
flmn_real = s.get_coeff(real=True,random=True)
# Be aware that there are pitfalls when populating real coefficient arrays by 
# explicitely assigning values to l,m,n etries.
# They have to satisfy f^l_mn=f^l_-m-n * (-1)^(m+n).

f_real = s.get_so3func(real=True)

##########################
# inverse real transform #
g_real = s.irsoft(flmn_real)
s.irsoft(flmn_real,out=f_real)
print(f'f_real and g_real are the same        = {np.allclose(f_real,g_real)}')

##########################
# forward real transform #
glmn_real = s.rsoft(g_real)
print(f'flmn_real and glmn_real are the same  = {np.allclose(glmn_real,flmn_real)}')
```



