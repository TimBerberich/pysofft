import numpy as np
from pysofft import _soft
from pysofft._soft import py
from pysofft._fftw_aligned_alloc import create_float64,create_complex128
from pysofft import fftw
import os
from collections import namedtuple
import multiprocessing
utils = _soft.utils

#Temporary untill propper logging is implemented
def log(txt):
    print(txt)

class Soft:
    r"""
    Python wrapper arround fortran class, that handles all transforms.

    Attributes
    ----------
    recurrence_types : namedtuple
        Contains all available Wigner recurrence types (currently kostelec and risbo)
    bw:int
        Bandwidth of the SO(3) Fourier transforn defining the SO(3) sampling.
    """
    _fortran_pointer = None
    _wisdom_path = os.path.expanduser('~/.config/pysofft/fftw_wisdom.dat')
    recurrence_types = namedtuple('RecurrenceTypes',['kostelec','risbo'])(int(py.get_kostelec_recurrence_val()),int(py.get_risbo_recurrence_val()))

    def __init__(self,
                 bw:int,
                 lmax:int|None = None,
                 precompute_wigners:bool = False,
                 recurrence_type:int|None = None,
                 init_ffts:bool = False,
                 fftw_flags:int = 0,
                 use_fftw_wisdom:bool = False):
        r"""
        Create a transform instance.
        
        Parameters
        ----------
        bw : int64
           Bandwidth of the SO(3) Fourier transform defining the SO(3) sampling.
        lmax: int64
           Maximum considered Wigner degree $l$. lmax has to satisfy 0<=lmax<bw.
        precompute_wigners: bool            
           Whether or not to precompute and store the $O(\mathrm{bw}^4)$ small Wigner-d matrix values $d^l_{m,n}(\\beta)$.
        recurrence_type: int64
           Selects the recurrence that is used to compute the $d^l_{m,n}(\\beta)$.
        init_ffts: bool
           Whether to allocate memory and create fft plans during initialization rather than lacy on the first call to a transform.
        """
        
        self._init_process_name = multiprocessing.current_process().name
        
        if recurrence_type not in self.recurrence_types:
            recurrence_type = self.recurrence_types.kostelec
        if lmax is None:
            lmax = bw-1
        self.use_fftw_wisdom = use_fftw_wisdom
        if use_fftw_wisdom:
            fftw.load_wisdom_from_file()
        self._fortran_pointer = py.py_init_soft(bw,lmax,precompute_wigners,init_ffts,recurrence_type,fftw_flags)
        if use_fftw_wisdom:
            fftw.save_wisdom_to_file()
        self._lmns = None
            
    def __del__(self):
        if multiprocessing.current_process().name == self._init_process_name:
            # free up any used fortran memory
            py.py_destroy(self._fortran_pointer)
    
    @property
    def bw(self):
        return py.py_get_bw(self._fortran_pointer)
    @property
    def lmax(self):
        return py.py_get_lmax(self._fortran_pointer)
    @lmax.setter
    def lmax(self,value):
        if isinstance(value,int):
            value = np.array(value)
        if isinstance(value,np.ndarray):            
            py.py_set_lmax(self._fortran_pointer,value)
        else:
            log(f'Warning: lmax ({lmax}) not an integer, abbort changing lmax.')
    @property
    def recurrence_type(self):
        return py.py_get_recurrence_type(self._fortran_pointer)
    @property
    def wigners_are_precomputed(self):
        return bool(py.py_wigners_are_precomputed(self._fortran_pointer))
    @property
    def fftw_flags(self):
        return py.py_get_fftw_flags(self._fortran_pointer)

    @property
    def coeff_indices(self):
        if self._lmns is None:
            if (self.recurrence_type == self.recurrence_types.kostelec) or self.wigners_are_precomputed:
                self._lmns=utils.get_coeff_degrees(self.bw)
            else:
                self._lmns=utils.get_coeff_degrees_risbo(self.bw)
            self._lmns.flags.writeable=False
        return self._lmns
    @property
    def euler_angles(self):
        beta = _soft.make_wigner.create_beta_samples(self.bw*2)
        alpha = _soft.make_wigner.create_alpha_gamma_samples(self.bw*2)
        gamma = alpha.copy()
        return {'gamma':gamma,'alpha':alpha,'beta':beta}
        
    @staticmethod
    def _fill_random(arr,seed=12345):
        rng = np.random.default_rng(seed)
        arr[...] = rng.random(arr.shape)
        if np.isdtype(arr.dtype,np.dtype(complex)):
            arr += 1.j*rng.random(arr.shape)
        return arr

    def get_coeff(self, real:bool=False, random:bool = False, seed:int=12345, raw:bool=False, howmany:int=0):
        r'''
        Create SO(3) coefficients array.

        Parameters
        ----------
        real:bool
           Whether to generate a coefficient array for real or complex functions over SO(3).
        random:bool
           If `True` fill with random numbers, if `False` fill with zeros.
        seed:int
           Seed for the random number generator
        raw:bool
           If `True` returns the pure numpy array without wrapping it in the CoeffSO3 class.
        howmany:int
           If grater than 0, create array of shape (howmany,n_coefficients) that is compatible with the *_many transforms.
        
        Returns
        -------
        coeff:ndarray
            `((4*bw^3-bw)/3,)|(howmany,(4*bw^3-bw)/3,) complex`: Wigner coefficient array.
        
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        s.get_coeff(random=True)
        ```
        '''
        n_coeff = utils.total_num_coeffs(self.bw)
        if howmany>0:
            coeff = create_complex128((howmany,n_coeff))
        else:
            coeff = create_complex128(n_coeff)
        if random:
            coeff = self._fill_random(coeff,seed=seed)
            if real:
                if (self.recurrence_type == self.recurrence_types.kostelec) or self.wigners_are_precomputed:
                    if howmany>0:
                        utils.enforce_real_sym_mnl_many(coeff.T,self.bw)
                    else:
                        utils.enforce_real_sym_mnl(coeff,self.bw)
                else:
                    if howmany>0:
                        utils.enforce_real_sym_lmn_many(coeff.T,self.bw)
                    else:
                        utils.enforce_real_sym_lmn(coeff,self.bw)
        else:
            coeff[:]=0
            
        if not raw:
            if (self.recurrence_type == self.recurrence_types.kostelec) or self.wigners_are_precomputed:
                coeff = CoeffSO3(coeff,self.coeff_indices,coeff_order = 'mnl')
            else:
                coeff = CoeffSO3(coeff,self.coeff_indices,coeff_order = 'lmn')        
        return coeff

    def get_so3func(self,real=False,random=False,seed=12345,howmany=0):
        r'''
        Create SO(3) function array.

        Parameters
        ----------
        real:bool
           Whether to generate a function array for real or complex functions over SO(3).
        random:bool
           If `True` fill with random numbers, if `False` fill with zeros.
        seed:int
           Seed for the random number generator
        howmany:int
           If grater than 0, create array of shape (howmany,2*bw,2*bw,2*bw) that is compatible with the *_many transforms.
        
        Returns
        -------
        so3func:ndarray
            `(2*bw,2*bw,2*bw)|(howmany,2*bw,2*bw,2*bw)   complex|real`: Array representing a function over SO(3) sampled on self.euler_angles.
        
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        s.get_so3func(random=True)
        ```
        '''
        euler_shape = utils.euler_shape(self.bw)
        if real:
            if howmany>0:
                func = create_float64((howmany,)+euler_shape)
            else:
                func = create_float64(euler_shape)
        else:
            if howmany>0:
                func = create_complex128((howmany,)+euler_shape)
            else:
                func = create_complex128(euler_shape)
        if random:
            func = self._fill_random(func,seed=seed)
        else:
            func[:]=0
        return func

    def enforce_real_symmetry(self,coeff):
        r'''
        Enforces the Wigner symmetry $f^l_{m,n} = f^l_{-m,-n} (-1)^{m+n}$ that occurs if they belong to a real
        valued function over SO(3).

        Parameters
        ----------
        coeff:ndarray
            `((4*bw**3-bw)/3,) complex`: Wigner coefficient array.
        
        Examples
        --------
        ```py
        from pysofft import Soft
        import numpy as np
        s = Soft(16)
        flmn = s.get_coeff(random=True)
        f = s.isoft(flmn)
        print(f'Max complex magnitude of f is = {np.max(np.abs(f.imag))}')
        flmn_real = flmn.copy()
        s.enforce_real_symmetry(flmn_real)
        f_real = s.isoft(flmn_real)
        print(f'Max complex magnitude of f_real is = {np.max(np.abs(f_real.imag))}')
        ```
        '''
        if (self.recurrence_type == self.recurrence_types.kostelec) or self.wigners_are_precomputed:
            utils.enforce_real_sym_mnl(coeff,self.bw)
        else:
            utils.enforce_real_sym_lmn(coeff,self.bw)
    # Transforms
    def _inverse_wigner_trf_cmplx(self,coeff,so3func,use_mp = False):
        py.py_inverse_wigner_trf_cmplx(self._fortran_pointer,coeff,(so3func.T),use_mp)
    def _forward_wigner_trf_cmplx(self,so3func,coeff,use_mp = False):
        py.py_forward_wigner_trf_cmplx(self._fortran_pointer,(so3func.T),coeff,use_mp)
    def _inverse_wigner_trf_real(self,coeff,so3func,use_mp = False):
        py.py_inverse_wigner_trf_real(self._fortran_pointer,coeff,(so3func.T),use_mp)
    def _forward_wigner_trf_real(self,so3func,coeff,use_mp = False):
        py.py_forward_wigner_trf_real(self._fortran_pointer,(so3func.T),coeff,use_mp)

    
    def soft(self,so3func:np.ndarray,out:np.ndarray|None=None,use_mp:bool=False):
        r'''
        Computes forward complex valued Fourier transform over SO(3)
        
        Parameters
        ----------
        so3func:ndarray
            `(2*bw,2*bw,2*bw) complex`: SO(3) function sampled at the eulerangles in `self.euler_angles`.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
           `((4*bw^3-bw)/3,) complex`: Computed Wigner coefficients.
                   
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        so3func = s.get_so3func(random=True)
        s.soft(so3func)
        ```
        '''
        if out is None:
            out=self.get_coeff()
        py.py_soft(self._fortran_pointer,(so3func.T),out,use_mp)
        return out
    def isoft(self,coeff:np.ndarray,out:np.ndarray|None=None,use_mp:bool=False):
        r'''
        Computes inverse complex valued Fourier transform over SO(3)

        Parameters
        ----------
        coeff:ndarray
            `((4*bw^3-bw)/3,) complex`: Wigner coefficients.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
           `(2*bw,2*bw,2*bw) complex`: Computed SO(3) function. 
        
        
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        coeff = s.get_coeff(random=True)
        s.isoft(coeff)
        ```
        '''
        if out is None:
            out=self.get_so3func()
        py.py_isoft(self._fortran_pointer,coeff,out.T,use_mp)
        return out
    def rsoft(self,so3func:np.ndarray,out:np.ndarray|None=None,use_mp:bool=False):
        r'''
        Computes forward real valued Fourier transform over SO(3)

        Parameters
        ----------
        so3func:ndarray
            `(2*bw,2*bw,2*bw) float`: SO(3) function sampled at the eulerangles in `self.euler_angles`.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `((4*bw^3-bw)/3) complex`: Computed Wigner coefficients.
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        so3func = s.get_so3func(random=True,real=True)
        s.rsoft(so3func)
        ```        
        '''
        if out is None:
            out=self.get_coeff(real=True)
        py.py_rsoft(self._fortran_pointer,so3func.T,out,use_mp)
        return out
    def irsoft(self,coeff:np.ndarray,out:np.ndarray|None=None,use_mp:bool=False):
        r'''
        Computes inverse real valued Fourier transform over SO(3)

        Parameters
        ----------
        coeff:ndarray
            `((4*bw^3-bw)/3,) complex`: Wigner coefficients.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
           `(2*bw,2*bw,2*bw) float`: Computed SO(3) function.
        
        Examples
        --------
        ```py
        from pysofft import Soft
        s = Soft(16)
        coeff = s.get_coeff(random=True,real=True)
        s.irsoft(coeff)
        ```
        '''
        if out is None:
            out=self.get_so3func(real=True)
        py.py_irsoft(self._fortran_pointer,coeff,out.T,use_mp)
        return out
    
    def soft_many(self,so3funcs:np.ndarray,out:np.ndarray|None=None,use_mp:bool=False):                
        r'''
        Computes forward complex valued Fourier transform over SO(3) for many functions at once.  
        Parallel processing is done over the individual function arrays if enabled.
        
        Parameters
        ----------
        so3func:ndarray
            `(n,2*bw,2*bw,2*bw) complex`: containing `n` separate functions over SO(3).
        out:None|ndarray
            If provided the transform result is written to this array.        
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `(n,(4*bw**3-bw)/3) complex`, Wigner coefficient arrays.
        
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        so3func = s.get_so3func(random=True,howmany=100)
        s.soft_many(so3func)
        ```
        '''
        if out is None:
            out=self.get_coeff(howmany=so3funcs.shape[0])
        py.py_soft_many(self._fortran_pointer,so3funcs.T,out.T,use_mp)
        return out
    def isoft_many(self,coeffs,out=None,use_mp=False):
        r'''
        Computes inverse complex valued Fourier transform over SO(3) for many coefficient arrays at once.  
        Parallel processing is done over the individual coefficient arrays if enabled.
        
        Parameters
        ----------
        so3func:ndarray
            `(n,(4*bw^3-bw)/3) complex`: containing `n` separate Wiegner coefficient arrays.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.

        Returns
        -------
        out:ndarray
            `(n,2*bw,2*bw,2*bw) complex`: Computed SO(3) functions. 
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        coeff = s.get_coeff(random=True,howmany=100)
        s.isoft_many(coeff)
        ```
        '''
        if out is None:
            out=self.get_so3func(howmany=coeffs.shape[0])
        py.py_isoft_many(self._fortran_pointer,coeffs.T,out.T,use_mp)
        return out
    
    def rsoft_many(self,so3funcs:np.ndarray,out:None|np.ndarray=None,use_mp:bool=False):
        r'''
        Computes forward real valued Fourier transform over SO(3) for many functions at once.  
        Parallel processing is done over the individual function arrays if enabled.
        
        Parameters
        ----------
        so3func: ndarray
            `(n,2*bw,2*bw,2*bw) float` storing `n` separate SO(3) functions.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `(n,(4*bw^3-bw)/3) complex`: Wigner coefficient arrays.
        
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        so3func = s.get_so3func(random=True,howmany=100,real=True)
        s.rsoft_many(so3func)
        ```
        '''
        if out is None:
            out=self.get_coeff(real=True,howmany=so3funcs.shape[0])
        py.py_rsoft_many(self._fortran_pointer,so3funcs.T,out.T,use_mp)
        return out
    
    def irsoft_many(self,coeffs:np.ndarray,out:None|np.ndarray=None,use_mp:bool=False):
        r'''
        Computes inverse real valued Fourier transform over SO(3) for many coefficient arrays at once.  
        Parallel processing is done over the individual coefficient arrays if enabled.
        
        Parameters
        ----------
        so3func: ndarray
            `(n,(4*bw^3-bw)/3) complex`: storing `n` separate Wigner coefficient arrays.
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
            `(n,2*bw,2*bw,2*bw) float`: computed SO(3) functions.        
        
        Examples
        -------
        ```py
        from pysofft import Soft
        s = Soft(16)
        coeff = s.get_coeff(random=True,howmany=100,real=True)
        s.irsoft_many(coeff)
        ```
        '''
        if out is None:
            out=self.get_so3func(real=True,howmany=coeffs.shape[0])
        py.py_irsoft_many(self._fortran_pointer,coeffs.T,out.T,use_mp)
        return out
    
    def integrate_over_so3(self,f:np.ndarray):
        r'''
        Nummerically integrate a function over SO(3)
        
        Parameters
        ----------
        f:ndarray
           ´(2*bw,2*bw,2*bw) complex|float`: SO(3) function to be integrated
        
        Returns
        -------
        complex|float
            Numerical íntegral of f
        
        Examples
        --------
        ```py
        from pysofft import Soft
        import numpy as np
        s = Soft(16)
        func = s.get_so3func()
        func[...] = 1
        vol_SO3 = s.integrate_over_so3(func)
        np.isclose(vol_SO3,8*np.pi**2)
        ```
        '''
        if f.dtype == complex:
            return py.py_integrate_over_so3_cmplx(self._fortran_pointer,f.T)
        else:
            return py.py_integrate_over_so3_real(self._fortran_pointer,f.T)
        
    def _inverse_wigner_trf_corr_cmplx(self,f_lm,g_lm,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(real=False)
        py.py_inverse_wigner_trf_corr_cmplx(self._fortran_pointer,f_lm,g_lm,out.T,use_mp)
        return out
    def _inverse_wigner_trf_corr_real(self,f_ml,g_ml,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(real=False)
        py.py_inverse_wigner_trf_corr_real(self._fortran_pointer,f_ml,g_ml,out.T,use_mp)
        return out
    
    def cross_correlation_ylm_cmplx(self,f_lm:np.ndarray,g_lm:np.ndarray,out:None|np.ndarray=None,use_mp:bool=False):
        r'''
        Computes the rotational cross-correlation between two complex functions $f,g: S^2 \rightarrow \mathbb{C}$ over the 2-sphere given by their spherical harmonic coefficients.
        That is, it computes
        
        $$C(\omega) = \int_{S^2} dx f(x) g(R_\omega x)^* $$
        
        where $\omega$ are elements of SO(3), $x$ are points on the 2-sphere and $R_\omega x$ describes the action of $\omega$ on a point $x$.
        
        Parameters
        ----------
        f_lm:ndarray
            `(bw**2,) complex`: Spherical harmonic coefficients of $f$. ($|m|\leq l$ and $f^l_m$ is stored at index `l*(l+1)+m`)
        g_lm:ndarray
            `(bw**2,) complex`: Spherical harmonic coefficients of $g$. ($|m|\leq l$ and $g^l_m$ is stored at index `l*(l+1)+m`)
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `(2*bw,2*bw,2*bw) complex` Cross-correlation between.
        
        Examples
        --------
        For an example see [rotational cross-correlation](/usage/corr).
        '''
        if out is None:
            out=self.get_so3func(real=False)
        py.py_cross_correlation_ylm_cmplx(self._fortran_pointer,f_lm,g_lm,out.T,use_mp)
        return out
    def cross_correlation_ylm_real(self,f_lm,g_lm,out=None,use_mp=False):
        r'''
        Computes the rotational cross-correlation between two real functions $f,g: S^2 \rightarrow \mathbb{R}$ over the 2-sphere given by their spherical harmonic coefficients.
        That is, it computes
        
        $$C(\omega) = \int_{S^2} dx f(x) g(R_\omega x)^* $$
        
        where $\omega$ are elements of SO(3), $x$ are points on the 2-sphere and $R_\omega x$ describes the action of $\omega$ on a point $x$.
        
        Parameters
        ----------
        f_lm:ndarray
            `(bw**2,) complex`: Spherical harmonic coefficients of $f$. ($0\leq m\leq l$ and $f^l_m$ is stored at index `(m*(2*bw-1-m))/2+l`)
        g_lm:ndarray
            `(bw**2,) complex`: Spherical harmonic coefficients of $g$. ($0\leq m\leq l$ and $g^l_m$ is stored at index `(m*(2*bw-1-m))/2+l`)
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `(2*bw,2*bw,2*bw) real` Cross-correlation between.
        
        Examples
        --------
        For an example see [rotational cross-correlation](/usage/corr).
        '''
        if out is None:
            out=self.get_so3func(real=True)            
        py.py_cross_correlation_ylm_real(self._fortran_pointer,f_lm,g_lm,out.T,use_mp)
        return out
    def cross_correlation_ylm_cmplx_3d(self,f_lms,g_lms,out = None,radial_sampling_points=None,radial_limits=None,use_mp=False):
        r'''
        Computes the rotational cross-correlation between two complex functions $f,g: \mathbb{R}^3 \rightarrow \mathbb{C}$ defined in spherical coordinates given by `n` sets of spherical harmonic coefficients, where `n` is the number or radial grid points
        That is, it computes
        
        $$C(\omega) = \int_{r} dr r^2 \int_{S^2} dx f(r,x) g(r,R_\omega x)^* $$
        
        where $\omega$ are elements of SO(3), $x$ are points on the 2-sphere, $r$ is the radial coordinate and $R_\omega x$ describes the action of $\omega$ on a point $x$. (Written in azimutal and polar angels we have $x=(\theta,\phi)$ and $dx=d\theta d\phi sin(\theta)$.)
        
        Parameters
        ----------
        f_lms:ndarray
            `(n,bw**2,) complex`: Spherical harmonic coefficients of $f$. ($|m|\leq l$ and $f^l_m$ is stored at index `l*(l+1)+m`)
        g_lms:ndarray
            `(n,bw**2,) complex`: Spherical harmonic coefficients of $g$. ($|m|\leq l$ and $g^l_m$ is stored at index `l*(l+1)+m`)
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `(2*bw,2*bw,2*bw) complex` Cross-correlation between.
        
        Examples
        --------
        For an example see [rotational cross-correlation](/usage/corr).
        '''
        if out is None:
            out=self.get_so3func(howmany=len(f_lms))
        if radial_sampling_points is None:
            radial_sampling_points = np.linspace(0,1,len(f_lmns))
        if radial_limits == None:
            radial_limits = [0,len(f_lms)]
        fortran_rlims = [radial_limits[0]+1,radial_limits[1]]
            
        py.py_cross_correlation_ylm_cmplx_3d(self._fortran_pointer,f_lms.T,g_lms.T,out.T,radial_sampling_points,fortran_rlims,use_mp)
        return out
    def cross_correlation_ylm_real_3d(self,f_lms,g_lms,out=None,radial_sampling_points=None,radial_limits=None,use_mp=False):
        r'''
        Computes the rotational cross-correlation between two real functions $f,g: \mathbb{R}^3 \rightarrow \mathbb{R}$ defined in spherical coordinates given by `n` sets of spherical harmonic coefficients, where `n` is the number or radial grid points
        That is, it computes
        
        $$C(\omega) = \int_{r} dr r^2 \int_{S^2} dx f(r,x) g(r,R_\omega x)^* $$
        
        where $\omega$ are elements of SO(3), $x$ are points on the 2-sphere, $r$ is the radial coordinate and $R_\omega x$ describes the action of $\omega$ on a point $x$. (Written in azimutal and polar angels we have $x=(\theta,\phi)$ and $dx=d\theta d\phi sin(\theta)$.)
        
        Parameters
        ----------
        f_lm:ndarray
            `(n,bw**2,) complex`: Spherical harmonic coefficients of $f$. ($0\leq m\leq l$ and $f^l_m$ is stored at index `(m*(2*bw-1-m))/2+l`)
        g_lm:ndarray
            `(n,bw**2,) complex`: Spherical harmonic coefficients of $g$. ($0\leq m\leq l$ and $g^l_m$ is stored at index `(m*(2*bw-1-m))/2+l`)
        out:None|ndarray
            If provided the transform result is written to this array.
        use_mp:bool
            Whether to allow parallel execution using OpenMP.
        
        Returns
        -------
        out:ndarray
            `(n,2*bw,2*bw,2*bw) real` Cross-correlation between.
        
        Examples
        --------
        For an example see [rotational cross-correlation](/usage/corr).
        '''
        if out is None:
            out=self.get_so3func(real=True,howmany=len(f_lms))
        if radial_sampling_points is None:
            radial_sampling_points = np.linspace(0,1,len(f_lms))
        if radial_limits == None:
            radial_limits = [0,len(f_lms)]
        fortran_rlims = [radial_limits[0]+1,radial_limits[1]]
        py.py_cross_correlation_ylm_real_3d(self._fortran_pointer,f_lms.T,g_lms.T,out.T,radial_sampling_points,fortran_rlims,use_mp)
        return out
    def _fft(self,f1,f2):
        py.py_fft(self._fortran_pointer,f1,f2)
    def _ifft(self,f1,f2):
        py.py_ifft(self._fortran_pointer,f1,f2)
    def _rfft(self,f1,f2):
        py.py_rfft(self._fortran_pointer,f1,f2)
    def _irfft(self,f1,f2):
        py.py_irfft(self._fortran_pointer,f1,f2)

    @staticmethod
    def rotate_ylm_cmplx(ylms,euler_angles):
        '''
        Rotates an array of spherical harmonic coefficients of a complex function
        by an array of rotations given as ZXZ euler angles.
    
        Input
        ylms         : (n_ylmns,spherical_harmonic_coefficients)  : (n,(2*bw-1)**2)   : complex
        euler_angles : (m_rotations,euler_angles)                 : (m,3)             : int
    
        Output
        ylns_rot     : (m_rotations,n_ylms,harmonic_coefficients) : (m,n,(2*bw-1)**2) : complex   
        '''
        euler_angles = np.atleast_2d(euler_angles)
        ylms = np.atleast_2d(ylms)
        #setup working arrays
        cos_b = np.cos(np.array(euler_angles[:,1])/2)
        sin_b = np.sin(np.array(-euler_angles[:,1])/2)
        # I need to invert the beta angles because my fortran code produces transposed d-matrices in python/c code
        # and we have d^l_mn(\beta) = d^l_(nm)(-beta) so substituting beta with -beta is the same as transposing.
        # because of this there is then minus in  np.sin(np.array(-euler_angles[:,1])/2)
        
        bw = int(np.sqrt(ylms.shape[-1]))
        sqrts = np.sqrt(np.arange(2*bw-1))
    
        #exp_a = np.exp(-1.j*euler_angles[:,0,None]*(np.arange(2*bw-1)[None,:]-bw+1))
        #exp_g = np.exp(-1.j*euler_angles[:,2,None]*(np.arange(2*bw-1)[None,:]-bw+1))
        exp_a = np.exp(-1.j*euler_angles[:,0,None]*np.arange(-bw+1,bw)[None,:])
        exp_g = np.exp(-1.j*euler_angles[:,2,None]*np.arange(-bw+1,bw)[None,:])
        
        dl_tmp = np.zeros((len(euler_angles),(2*bw-1)*(2*bw-1)),float,order="C")
        D_conj = np.zeros((len(euler_angles),(2*bw-1),(2*bw-1)),complex,order="C")
        ylms_rot = np.zeros((len(euler_angles),)+ylms.shape,dtype=complex,order="C") 
        for l in range(bw):
            # use recursion to compute the Wigner small-d matrices for the current oder l
            dl_tmp[:,0:(2*l+1)*(2*l+1)] = _soft.make_wigner.wigner_mn_recurrence_risbo(dl_tmp[:,0:(2*l-1)*(2*l-1)],l,cos_b,sin_b,sqrts)
            # slice of relevant orders in exp_a, exp_b and D
            l_slice_sym = slice(bw-1-l,bw+l)
            # populate Wiegner D^l_mn(alpha,beta,gamma)
            D_conj[:,l_slice_sym,l_slice_sym] =  dl_tmp[:,0:(2*l+1)*(2*l+1)].reshape(len(euler_angles),(2*l+1),(2*l+1))*exp_g[:,None,l_slice_sym]*exp_a[:,l_slice_sym,None]
            # slice to order l part of the input ylmns
            tmp = _soft.utils.lmc_slice(l)
            lmc_slice = slice(tmp[0]-1,tmp[1])
            lmc_slice_reverse = slice(tmp[1]-1,tmp[0])
            # Perform the rotation in order l
            np.matmul(ylms[:,lmc_slice],D_conj[:,l_slice_sym,l_slice_sym],out=ylms_rot[...,lmc_slice])
    
        # If only one rotation or one set of harmonic coefficients is provided
        # reduce the unecesary dimension in the output array
        if euler_angles.shape[0]==1 or ylms.shape[0]==1:
            ylms_rot = np.squeeze(ylms_rot)
        return ylms_rot
    @staticmethod
    def rotate_ylm_real(ymls,euler_angles):
        '''
        Terribly inefficient conveniece function.
        Rotates an array of spherical harmonic coefficients of a complex function
        by an array of rotations given as ZXZ euler angles.
        Works by converting to complex representation and apply rotate_ylm_cmplx.
        
        Parameters
        ----------
        ylms:ndarray
            `(n,(bw*(bw+1))/2) complex`: set of spherical harmonic coefficients. (`bw` here is not the `bw` of of `Soft` object but of the input coefficient array)
        euler_angles:ndarray 
            `(m,3) float` Array of euler angels.
        
        Returns
        -------
        ylns_rot:ndarray
            `(m,n,(bw*(bw+1))/2) complex` Rotated spherical harmonic coefficients.
        
        '''
        bw = int( (np.sqrt(1+8*ymls.shape[-1])-1)/2 )
        ymls = np.atleast_2d(ymls)
        ylms = np.zeros((ymls.shape[:-1]+(bw**2,)),dtype = complex)
    
        
        for l in range(bw):
            for m in range(l+1):
                ml_id = _soft.utils.mlr(m,l,bw)-1
                lm_id = _soft.utils.lmc(l,m)-1
                lm_neg_id = _soft.utils.lmc(l,-m)-1
                ylms[...,lm_id] = ymls[...,ml_id]
                ylms[...,lm_neg_id] = ylms[...,lm_id].conj()*(-1)**m
        
        ylms_rot = Soft.rotate_ylm_cmplx(ylms,euler_angles)
        ymls_rot = np.zeros_like(ymls)
    
        for m in range(bw):
            for l in range(m,bw):
                ml_id = _soft.utils.mlr(m,l,bw)-1
                lm_id = _soft.utils.lmc(l,m)-1
                ymls_rot[...,ml_id] = ylms_rot[...,lm_id]
        return ymls_rot


# Syntactic sugar section
class CoeffSO3(np.ndarray):
    r"""
    Numpy ndarray subclass adding easy access to individual Wigner coefficients $f^l_{n,k}$ via `self.lmn[]`.
    
    Attributes
    ----------
    bw:int
        Bandwidth of the distribution.
    lmn:CoeffView
        allows to acces $f^l_{m,n}$ by its indices $l,n,k$.
    """
    _coeff_orders = ['mnl','lmn']
    def __new__(cls, coeff,lmn_indices,coeff_order = 'mnl'):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(coeff).view(cls)
        
         # add the new attribute to the created instance
        lmn_indices.flags.writeable = False
        obj._lmns = lmn_indices
        obj.lmn = CoeffView(obj,lmn_indices,coeff_order)
        obj.bw = obj.lmn.bw
        obj.coeff_order = obj.lmn.coeff_order
        
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None:
            return
        
        self._lmns = getattr(obj, '_lmns', None)
        if self._lmns is None:
            self.lmn = None
            self.bw = None
            self.coeff_order = None
        else:
            self.lmn = CoeffView(self,self._lmns)
            self.bw = obj.lmn.bw
            self.coeff_order = obj.lmn.coeff_order

    def __array_ufunc__(self, ufunc, method, *inputs, out=None, **kwargs):
        # Enables ufuncs like +,-,*,/ to preserve the additional attributes _lmns,lmn,bw
        
        # convert args and output to ndarray
        args = []
        for i, input_ in enumerate(inputs):
            if isinstance(input_, CoeffSO3):
                args.append(input_.view(np.ndarray))
            else:
                args.append(input_)
        outputs = out
        if outputs:
            out_args = []
            for j, output in enumerate(outputs):
                if isinstance(output, CoeffSO3):
                    out_args.append(output.view(np.ndarray))
                else:
                    out_args.append(output)
            kwargs['out'] = tuple(out_args)
        else:
            outputs = (None,) * ufunc.nout

        # compute ufuncs
        results = super().__array_ufunc__(ufunc, method, *args, **kwargs)
        if results is NotImplemented:
            return NotImplemented

        if ufunc.nout == 1:
            results = (results,)
            
        if method == 'at':
            if isinstance(inputs[0], c):
                inputs[0].info = info
            return
        
        # Convert result back to CoeffSO3
        results = list((np.asarray(result).view(CoeffSO3)
                         if output is None else output)
                        for result, output in zip(results, outputs))
        
        #If the shape did not change add back the propper attributes
        if len(inputs)>0:
            coeff_inputs = tuple(i for i in inputs if isinstance(i,CoeffSO3))
            for res_id,res in enumerate(results):
                res_modified = False
                for inp in coeff_inputs:
                    if res.shape == inp.shape:
                        res._lmns = inp._lmns
                        res.lmn = CoeffView(res,res._lmns)
                        res.bw = res.lmn.bw
                        res_modified = True
                        break
                if not res_modified:
                    results[res_id]=res.view(np.ndarray)
                
        return results[0] if len(results) == 1 else results
    
class CoeffView:
    _coeff_orders = ['mnl','lmn']
    """
    Helper class that allows to access a Wigner coefficients $f^l_{n,k}$ by its degrees l,n,k,
    i.e. `coeffview[l,n,k]` returns  $f^l_{n,k}$. For example, it is also possible to use slice notation like
    charview[:,0,:], which returns a view into $f^l_{n,k}$ corresponding to all values with n=0. 

    Attributes
    ----------
    bw:int
        Bandwidth
    array : ndarray
        `(4*bw^3-bw)/3,) complex`: Wigner coefficient array
    ls:ndarray
        `(bw,), int64`: l indices 0,...,bw
    ns:ndarray,
        `(2*bw+1,), int64`: n indices  0,...,bw,-bw+1,...,-1
    ks:ndarray
        `(2*bw+1,), int64`: k indices  0,...,bw,-bw+1,...,-1
    _lmns:ndarray
        `(3,(4*bw**4-bw)/3), int64`: l,n,k values associated to each entry in array.
    """
    def __init__(self,coeff_array,lmn_indices,coeff_order='mnl'):
        self.array = coeff_array
        self._lmns = lmn_indices
        self.ls = np.sort(np.unique(lmn_indices[:,0]))
        bw = len(self.ls)
        self.bw = bw
        self.ns = np.roll(np.sort(np.unique(lmn_indices[:,1])),bw)
        self.ks = np.roll(np.sort(np.unique(lmn_indices[:,2])),bw)
        if coeff_order in self._coeff_orders:
            self.coeff_order = coeff_order
        else:
            AssertionError(f'provided coeff_order {coeff_order} not contained in {self._coeff_orders}.')
    def _check_index(self,l,m,n):
        try:
            assert (0<=l and l<self.bw), f'l={l} does not satisfy 0<=l<bw={self.bw}'
            if self.coeff_order == 'lmn':
                assert (abs(m)<=l and abs(n)<=l), f'm={m} or n={n} do not satisfy |m|<=l={l} |n|<=l={l}'
            else:
                mn = max(abs(m),abs(n))
                assert (l>=mn), f'l={l} does not satisfy  l>=Max(|m|,|n|)={mn}'
        except AssertionError as e:
            raise e
        
    def _get_mask(self,items):
        if not isinstance(items,tuple):
            items = (items,)
        assert len(items)<=3
        if (len(items)==3) and np.prod(tuple(isinstance(i,int) for i in items)):
            self._check_index(items[0],items[1],items[2])
            if self.coeff_order == 'lmn':
                id_ = utils.coeff_location_lmn(items[0],items[1],items[2])-1
            else:
                id_ = utils.coeff_location_mnl(items[1],items[2],items[0],self.bw)-1
            out = slice(id_,id_+1)
        else:
            selection_ids=[ids[sel] for sel,ids in zip(items,[self.ls,self.ns,self.ks])]
            masks = [np.isin(lmn,sids) for sids,lmn in zip(selection_ids,self._lmns.T)]
            out = np.prod(masks,axis = 0,dtype=bool)
        return out
    def __getitem__(self, items):        
        mask = self._get_mask(items)
        return self.array[mask]
    def __setitem__(self,items,value):
        m = self._get_mask(items)
        self.array[m] = value
