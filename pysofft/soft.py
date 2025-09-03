import numpy as np
from pysofft import _soft2
from pysofft._soft2 import py
from pysofft._soft2 import softclass
import os

#Temporary untill propper logging is implemented
def log(txt):
    print(txt)
    
class Soft:
    _fortran_pointer = None
    _wisdom_path = os.path.expanduser('~/.config/pysofft/fftw_wisdom.dat')
    enable_fftw_wisdom = False
    
    def __init__(self,
                 bw,
                 lmax=None,
                 precompute_wigners = False,
                 init_ffts=False,
                 fftw_flags = 0,
                 enable_fftw_wisdom=False,
                 fftw_wisdom_path=None):
        
        self.enable_fftw_wisdom = enable_fftw_wisdom
        if enable_fftw_wisdom:
            if isinstance(fftw_wisdom_path,str):
                self._wisdom_path = fftw_wisdom_path
            wisdom_dir = os.path.dirname(self._wisdom_path)
            if not os.path.exists(wisdom_dir):
                os.makedirs(wisdom_dir)
            softclass.import_fftw_wisdom(self._wisdom_path)
        self._fortran_pointer = py.py_init_soft(bw,lmax,precompute_wigners,init_ffts,fftw_flags)
        if enable_fftw_wisdom:
            softclass.export_fftw_wisdom(self._wisdom_path)
            
    def __del__(self):
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
    def fftw_flags(self):
        return py.py_get_fftw_flags(self._fortran_pointer)
    def reset(self,
              bw,
              lmax=None,
              precompute_wigners = False,
              init_ffts=False,
              fftw_flags = 0,
              enable_fftw_wisdom=False):
        self.enable_fftw_wisdom = enable_fftw_wisdöm
        if enable_fftw_wisdom:
            softclass.import_fftw_wisdom(self._wisdom_path)
        py.py_reset(self._fortran_pointer,bw,lmax,precompute_wigners,init_ffts,fftw_flags)
        if enable_fftw_wisdom:
            softclass.export_fftw_wisdom(self._wisdom_path)
    def init_ffts(self,real_fft=False):
        py.py_init_fft(self._fortran_pointer,real_fft)
        
    # Transforms
    def _inverse_wigner_trf_cmplx(self,coeff,so3func,use_mp = False):
        py.py_inverse_wigner_trf_cmplx(self._fortran_pointer,coeff,so3func,use_mp)
    def _forward_wigner_trf_cmplx(self,so3func,coeff,use_mp = False):
        py.py_forward_wigner_trf_cmplx(self._fortran_pointer,so3func,coeff,use_mp)
    def _inverse_wigner_trf_real(self,coeff,so3func,use_mp = False):
        py.py_inverse_wigner_trf_real(self._fortran_pointer,coeff,so3func,use_mp)
    def _forward_wigner_trf_real(self,so3func,coeff,use_mp = False):
        py.py_forward_wigner_trf_real(self._fortran_pointer,so3func,coeff,use_mp)
    def soft(self,so3func,coeff,use_mp=False):
        py.py_soft(self._fortran_pointer,so3func,coeff,use_mp)
    def isoft(self,coeff,so3func,use_mp=False):
        py.py_isoft(self._fortran_pointer,coeff,so3func,use_mp)
    def rsoft(self,so3func,coeff,use_mp=False):
        py.py_rsoft(self._fortran_pointer,so3func,coeff,use_mp)
    def irsoft(self,coeff,so3func,use_mp=False):
        py.py_irsoft(self._fortran_pointer,coeff,so3func,use_mp)
    def soft_many(self,so3funcs,coeffs,use_mp=False):
        py.py_soft_many(self._fortran_pointer,so3funcs,coeffs,use_mp)
    def isoft_many(self,coeffs,so3funcs,use_mp=False):
        py.py_isoft_many(self._fortran_pointer,coeffs,so3funcs,use_mp)
    def rsoft_many(self,so3funcs,coeffs,use_mp=False):
        py.py_rsoft_many(self._fortran_pointer,so3funcs,coeffs,use_mp)
    def irsoft_many(self,coeffs,so3funcs,use_mp=False):
        py.py_irsoft_many(self._fortran_pointer,coeffs,so3funcs,use_mp)
    def cross_correlation_ylm_cmplx(self,f_lm,g_lm,cc,use_mp=False):
        py.py_cross_correlation_ylm_cmplx(self._fortran_pointer,f_lm,g_lm,cc,use_mp)
    def corss_correlation_ylm_cmplx_3d(self,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp=False):
        py.py_cross_correlation_ylm_cmplx_3d(self._fortran_pointer,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
    def cross_correlation_ylm_real(self_int,f_lm,g_lm,cc,use_mp=False):
        py.py_cross_correlation_ylm_real(self._fortran_pointer,f_lm,g_lm,cc,use_mp)
    def corss_correlation_ylm_real_3d(self_int,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp=False):
        py.py_cross_correlation_ylm_real_3d(self._fortran_pointer,f_lms,g_lms,cc,radial_sampling_points,radial_limits,use_mp)
    def fft(self,f1,f2):
        py.py_fft(self._fortran_pointer,f1,f2)
    def ifft(self,f1,f2):
        py.py_ifft(self._fortran_pointer,f1,f2)
    def rfft(self,f1,f2):
        py.py_rfft(self._fortran_pointer,f1,f2)
    def irfft(self,f1,f2):
        py.py_irfft(self._fortran_pointer,f1,f2)
