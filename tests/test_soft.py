import numpy as np
import pysofft
from math import factorial
from pysofft import _soft
from pysofft import soft,Soft
from multiprocessing import Pool
import pytest

try:
    import shtns
except Exception:
    shtns=None

class TestMultiprocessingCompatibility:
    def test_fork_safety(self):
        s = Soft(32,use_fftw_wisdom=True,init_ffts=True)
        coeff = s.get_coeff(howmany=100,random=True)
        
        def forked_soft(i):
            f = s.isoft(coeff[i],use_mp=False)
            return f
        
        def reinit_soft(i):
            s = Soft(32,use_fftw_wisdom=True,init_ffts=True)
            f = s.isoft(coeff[i])
            return f
        
        print(locals())
        with Pool(8) as p:
            d2 = p.map(forked_soft,np.arange(len(coeff)))
        d2 = np.array(d)

        with Pool(8) as p:
            d3 = p.map(reinit_soft,np.arange(len(coeff)))
        d3 = np.array(d3)

        d = s.isoft(coeff,use_mp=True)
        
        assert np.allclose(d2,d),'Pool using forked Soft instance is not the same as native computation.'
        assert np.allclose(d3,d),'Pool using separate Soft instances is not the same as native computation.'
        
    
class TestRotate:
    @pytest.mark.skipif(shtns is None, reason="python pakage 'shtns' is not installed but required for this test.")
    def test_rotate(self):
        #setup
        lmax = 64
        sh = shtns.sht(lmax)
        n_coeff = (lmax+1)**2
        n_phi = 256
        n_theta = 128
        sh.set_grid(polar_opt=0,flags=shtns.sht_gauss)
        sh.set_grid(nlat=n_theta,nphi=n_phi,polar_opt=0,flags=shtns.sht_gauss)
        phis=2*np.pi*np.arange(n_phi)/(n_phi*sh.mres)
        thetas=np.arccos(sh.cos_theta)
        
        
        rs = np.arange(10)
        y_dens = np.zeros((len(thetas),len(phis)),dtype = complex)
        delta = np.pi/64
        y_dens[:,(phis<delta) + (phis>(2*np.pi-delta)) ]=1
        y_dens_lm = sh.analys_cplx(y_dens)
        
        euler_angles = np.array([[0,0,0],[np.pi/2,0,0],[np.pi,0,0],[np.pi*3/2,0,0]])
        res = np.squeeze(soft.rotate_ylm_cmplx(y_dens_lm[None,...],euler_angles))
        
        y_out = [sh.synth_cplx(r) for r in res]
        y_argmax = tuple(np.argmax(o.real.mean(axis = 0)) for o in y_out)

        assert y_argmax == (0,64,128,192), 'Density maxima do not coincide with used righthanded rotations around z-axis.'
        
        z_dens = np.zeros((len(thetas),len(phis)),dtype = complex)
        delta = np.pi/64
        z_dens[(thetas<(np.pi/4+delta)) * (thetas>(np.pi/4-delta)) ,:]=1
        z_dens *= y_dens
        
        z_dens_lm = sh.analys_cplx(z_dens)
        euler_angles = np.array([[0,0,0],[0,np.pi/2,0],[0,np.pi,0],[0,np.pi*3/2,0]])
        res = np.squeeze(soft.rotate_ylm_cmplx(z_dens_lm[None,...],euler_angles))
        
        z_out = [sh.synth_cplx(r) for r in res]
        phi_ids = tuple(np.argmax(o.real.mean(axis = 0)) for o in z_out)
        theta_ids = tuple(np.argmax(o.real.mean(axis = 1)) for o in z_out)
        assert (phi_ids == (0,0,128,128)) and (theta_ids == (31,96,96,31)), 'Density maxima do not coincide with right handed rotations around y-axis'

        #print(f"phi_loc = {[np.argmax(o.real.mean(axis=0))for o in z_out]}")
        #print(f"theta_loc = {[np.argmax(o.real.mean(axis=1))for o in z_out]}")
        
