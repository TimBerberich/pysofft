# Harmonic Coefficient 
## Accessing hamronic coefficients
In __python__ this is best done via the CoeffSO3 class, e.g.

	from pysofft import Soft
	bw = 64
	s = Soft(bw)
	
	f_lmn = s.get_coeff(random=True)
	
	l=4
	m=-1
	n=2
	
	# access single coefficient
	print(f_lmn.lmn[l,m,n])
	
	# slicing is also possible
	print(f_lmn.lmn[l,m,:])
	
	# as well as value asignment
	f_lmn.lmn[l,m,:] = 1 + 1.j

In Fortran the utils module contains the functions __coeff_slice__ and __coeff_location__
that return the fortran index of a specific coefficient.

## Internal memmory layout

The internal memory layout of the harmonic coefficients $f_{lmn}$ is optimized for transform speed.
They are orderd such that access via loops of the following type are contiguous in memory.

    do m1=0,bw-1
       do m2=m1,bw-1
         f_{m1,m2}

	     if m1==0 and m2==0, cycle  
         f_{-m2,-m1}
         
		 if m1/=m2
            f_{m2,m1}
            f_{-m1,-m2}
         
		 if m1==0 or m2==0
         f_{m1,-m2}
         f_{-m1,m2}
         
		 if m1==m2, cycle
         f_{m2,-m1}
         f_{-m2,m1}

where f_{m1,m2} stands for the indices $f_{l,m1,m2}$ for all $l=\mathrm{max}(|m1|,|m2|),\ldots,bw-1$. 

