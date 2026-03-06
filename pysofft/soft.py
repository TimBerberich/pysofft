import numpy as np
from pysofft import _soft
from pysofft._soft import py
from pysofft._soft import softclass
from pysofft import fftw
import os
from collections import namedtuple
import multiprocessing
utils = _soft.utils

#Temporary untill propper logging is implemented
def log(txt):
    print(txt)

#> ---
## @brief fubar
##
class Soft:
    r"""
    Python wrapper arround fortran class, that handles all transforms.

    Attributes
    ----------
    bw : int64
        Bandwidth of the SO(3) fourier Transform.
    lmax: int64
        Maximum considered Wigner degree $l$. lmax has to satisfy 0<=lmax<bw.
    precompute_wigners: bool
        Whether or not to precompute and store the $O(\mathrm{bw}^4)$ small Wigner-d matrix values $d^l_{m,n}(\\beta)$.
    recurrence_type: int64
        Selects the recurrence that is used to compute the $d^l_{m,n}(\\beta)$.
    init_ffts: bool
        Whether or not to allocate memory and create fft plans during 
    """
    _fortran_pointer = None
    _wisdom_path = os.path.expanduser('~/.config/pysofft/fftw_wisdom.dat')
    enable_fftw_wisdom = False
    recurrence_types = namedtuple('RecurrenceTypes',['kostelec','risbo'])(int(softclass.kostelec_recurrence),int(softclass.risbo_recurrence))

    
    
    def __init__(self,
                 bw,
                 lmax=None,
                 precompute_wigners = False,
                 recurrence_type = None,
                 init_ffts=False,
                 fftw_flags = 0,
                 use_fftw_wisdom=False):
        r"""
        Parameters
        ----------
        bw : int64
           Bandwidth of the SO(3) fourier Transform.
        lmax: int64
           Maximum considered Wigner degree $l$. lmax has to satisfy 0<=lmax<bw.
        precompute_wigners: bool
           Whether or not to precompute and store the $O(\mathrm{bw}^4)$ small Wigner-d matrix values $d^l_{m,n}(\\beta)$.
        recurrence_type: int64
           Selects the recurrence that is used to compute the $d^l_{m,n}(\\beta)$.
        init_ffts: bool
           Whether or not to allocate memory and create fft plans during 
        """
        
        self._init_proces_name = multiprocessing.current_process().name
        
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
        if multiprocessing.current_process().name == self._init_proces_name:
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
    def fftw_flags(self):
        return py.py_get_fftw_flags(self._fortran_pointer)

    @property
    def coeff_indices(self):
        if self._lmns is None:
            self._lmns=utils.get_coeff_degrees(self.bw)
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

    def get_coeff(self, real=False, random = False, seed=12345, raw=False, howmany=0):
        if howmany>0:
            coeff = utils.get_empty_coeff_many(self.bw,howmany).T #Transpose to go from F_CONTIGUOUS to C_CONTIGUOUS array
        else:
            coeff = utils.get_empty_coeff(self.bw)
        if random:
            coeff = self._fill_random(coeff,seed=seed)
            if real:
                if self.recurrence_type == self.recurrence_types.kostelec:
                    utils.enforce_real_sym(coeff,self.bw)
                else:
                    utils.enforce_real_sym_lmn(coeff,self.bw)                    
        if not raw:
            coeff = CoeffSO3(coeff,self.coeff_indices)
        return coeff

    def get_so3func(self,real=False,random=False,seed=12345,howmany=0):
        if real:
            if howmany>0:
                func = utils.get_empty_so3func_real_many(self.bw,howmany)
            else:
                func = utils.get_empty_so3func_real(self.bw)
        else:
            if howmany>0:
                func = utils.get_empty_so3func_cmplx_many(self.bw,howmany)
            else:
                func = utils.get_empty_so3func_cmplx(self.bw)
        if random:
            func = self._fill_random(func,seed=seed)
        return func.T #Transpose to go from F_CONTIGUOUS to C_CONTIGUOUS array
    
    # Transforms
    def _inverse_wigner_trf_cmplx(self,coeff,so3func,use_mp = False):
        py.py_inverse_wigner_trf_cmplx(self._fortran_pointer,coeff,(so3func.T),use_mp)
    def _forward_wigner_trf_cmplx(self,so3func,coeff,use_mp = False):
        py.py_forward_wigner_trf_cmplx(self._fortran_pointer,(so3func.T),coeff,use_mp)
    def _inverse_wigner_trf_real(self,coeff,so3func,use_mp = False):
        py.py_inverse_wigner_trf_real(self._fortran_pointer,coeff,(so3func.T),use_mp)
    def _forward_wigner_trf_real(self,so3func,coeff,use_mp = False):
        py.py_forward_wigner_trf_real(self._fortran_pointer,(so3func.T),coeff,use_mp)

    
    def soft(self,so3func,out=None,use_mp=False):
        if out is None:
            out=self.get_coeff()
        py.py_soft(self._fortran_pointer,(so3func.T),out,use_mp)
        return out
    def isoft(self,coeff,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func()
        py.py_isoft(self._fortran_pointer,coeff,out.T,use_mp)
        return out
    def rsoft(self,so3func,out=None,use_mp=False):
        if out is None:
            out=self.get_coeff(real=True)
        py.py_rsoft(self._fortran_pointer,so3func.T,out,use_mp)
        return out
    def irsoft(self,coeff,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(real=True)
        py.py_irsoft(self._fortran_pointer,coeff,out.T,use_mp)
        return out
    def soft_many(self,so3funcs,out=None,use_mp=False):
        if out is None:
            out=self.get_coeff(howmany=so3funcs.shape[0])
        py.py_soft_many(self._fortran_pointer,so3funcs.T,out.T,use_mp)
        return out
    def isoft_many(self,coeffs,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(howmany=coeffs.shape[0])
        py.py_isoft_many(self._fortran_pointer,coeffs.T,out.T,use_mp)
        return out
    def rsoft_many(self,so3funcs,out=None,use_mp=False):
        if out is None:
            out=self.get_coeff(real=True,howmany=so3funcs.shape[0])
        py.py_rsoft_many(self._fortran_pointer,so3funcs.T,out.T,use_mp)
        return out
    def irsoft_many(self,coeffs,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(real=True,howmany=coeffs.shape[0])
        py.py_irsoft_many(self._fortran_pointer,coeffs.T,out.T,use_mp)
        return out
    def integrate_over_so3_cmplx(self,f):
        return py.py_integrate_over_so3_cmplx(self._fortran_pointer,f.T)
    def integrate_over_so3_real(self,f):
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
    def cross_correlation_ylm_cmplx(self,f_lm,g_lm,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(real=False)
        py.py_cross_correlation_ylm_cmplx(self._fortran_pointer,f_lm,g_lm,out.T,use_mp)
        return out
    def cross_correlation_ylm_cmplx_3d(self,f_lms,g_lms,out = None,radial_sampling_points=None,radial_limits=None,use_mp=False):
        if out is None:
            out=self.get_so3func(howmany=len(f_lms))
        if radial_sampling_points is None:
            radial_sampling_points = np.linspace(0,1,len(f_lmns))
        if radial_limits == None:
            radial_limits = [0,len(f_lms)]
        fortran_rlims = [radial_limits[0]+1,radial_limits[1]]
            
        py.py_cross_correlation_ylm_cmplx_3d(self._fortran_pointer,f_lms.T,g_lms.T,out.T,radial_sampling_points,fortran_rlims,use_mp)
        return out
    def cross_correlation_ylm_real(self,f_lm,g_lm,out=None,use_mp=False):
        if out is None:
            out=self.get_so3func(real=True)            
        py.py_cross_correlation_ylm_real(self._fortran_pointer,f_lm,g_lm,out.T,use_mp)
        return out
    def cross_correlation_ylm_real_3d(self,f_lms,g_lms,out=None,radial_sampling_points=None,radial_limits=None,use_mp=False):
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
def rotate_ylm_real(ymls,euler_angles):
    '''
    Terribly inefficient conveniece function.
    Rotates an array of spherical harmonic coefficients of a complex function
    by an array of rotations given as ZXZ euler angles.
    Works by converting to complex representation and apply rotate_ylm_cmplx.

    Input
    ylms         : (n_ylmns,spherical_harmonic_coefficients)  : (n,(2*bw-1)**2)   : complex
    euler_angles : (m_rotations,euler_angles)                 : (m,3)             : int

    Output
    ylns_rot     : (m_rotations,n_ylms,harmonic_coefficients) : (m,n,(2*bw-1)**2) : complex   
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
    
    ylms_rot = rotate_ylm_cmplx(ylms,euler_angles)
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
    A class to wrap SO3 fourier coefficients that provides easy access to individual coefficents
    
    Attributes
    ----------
    bw : int
    Bandwidth of the distribution.
    lnk : CharView
    allows to acces $M_{lnk}$ by its indices l,n,k
    """
    def __new__(cls, coeff,lmns):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(coeff).view(cls)
        
         # add the new attribute to the created instance
        lmns.flags.writeable = False
        obj._lmns = lmns
        obj.lmn = CoeffView(obj,lmns)
        obj.bw = obj.lmn.bw
        
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
        else:
            self.lmn = CoeffView(self,self._lmns)
            self.bw = obj.lmn.bw

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
    """
    Helper class that allows to access a characteristic function $M_{l,n,k}$ by its indices l,n,k.
    i.e. charview[l,n,k] returns  $M_{l,n,k}$. It is also possible to use slice notation i.e.
    charview[:,0,:] returns a view into $M_{l,n,k}$ corresponding to all values with n=0. 

    Attributes
    ----------
    bw : int
    Bandwidth
    array : np.ndarray, ((4*bw**4-bw)/3,), complex
    characteristic function array
    ls : np.ndarray, (bw,), int64
    l indices 0,...,bw
    ns : np.ndarray, (2*bw+1), int64
    n indices  0,...,bw,-bw+1,...,-1
    ks : np.ndarray, (2*bw+1), int64
    k indices  0,...,bw,-bw+1,...,-1
    _lnks: np.ndarray, (3,(4*bw**4-bw)/3), int64
    l,n,k value grid
    _lookups : np.ndarray, ((4*bw**4-bw)/3,), int64
    array ids associated to each entry of _lnks  
    """
    def __init__(self,coeff_array,lnks):
        self.array = coeff_array
        self._lnks = lnks
        self.ls = np.sort(np.unique(lnks[:,0]))
        bw = len(self.ls)
        self.bw = bw
        self.ns = np.roll(np.sort(np.unique(lnks[:,1])),bw)
        self.ks = np.roll(np.sort(np.unique(lnks[:,2])),bw)
    def _get_mask(self,items):
        if not isinstance(items,tuple):
            items = (items,)
        assert len(items)<=3
        if (len(items)==3) and np.prod(tuple(isinstance(i,int) for i in items)):
            id_ = utils.coeff_location(items[1],items[2],items[0],self.bw)-1
            out = slice(id_,id_+1)
        elif (len(items)==3) and np.prod(tuple(isinstance(i,int) for i in items[1:])) and isinstance(items[0],slice):
            id_ = utils.coeff_slice(items[1],items[2],self.bw)
            id_[0] -= 1
            if (items[0].start is None) and (items[0].stop is None):
                step = items[0].step
                if step is None:
                    out = slice(id_[0],id_[1])
                elif step>=0:
                    out = slice(id_[0],id_[1],step)
                else:
                    out = slice(id_[1]-1,id_[0]-1,step)
            else:
                selection_ids=[ids[sel] for sel,ids in zip(items,[self.ls,self.ns,self.ks])]
                masks = [np.isin(lnk,sids) for sids,lnk in zip(selection_ids,self._lnks.T)]
                out = np.prod(masks,axis = 0,dtype=bool)
        else:
            selection_ids=[ids[sel] for sel,ids in zip(items,[self.ls,self.ns,self.ks])]
            masks = [np.isin(lnk,sids) for sids,lnk in zip(selection_ids,self._lnks.T)]
            out = np.prod(masks,axis = 0,dtype=bool)
        return out
    def __getitem__(self, items):        
        mask = self._get_mask(items)
        return self.array[mask]
    def __setitem__(self,items,value):
        m = self._get_mask(items)
        self.array[m] = value
