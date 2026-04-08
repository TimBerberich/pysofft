import numpy as np
import pysofft
from math import factorial
from pysofft import _soft
from pysofft import soft,Soft
from pysofft._fftw_aligned_alloc import create_float64,create_complex128
from multiprocessing import Pool
import pytest
import weakref
import gc

try:
    import shtns
except Exception:
    shtns=None

def forked_soft(args):
    s,c = args
    f = s.isoft(c,use_mp=False)
    return f

class TestMultiprocessingCompatibility:
    def test_fork_safety(self):
        s = Soft(16,use_fftw_wisdom=True,init_ffts=True)
        coeff = s.get_coeff(howmany=100,random=True)
                
        with Pool(8) as p:
            d2 = p.map(forked_soft,[(s,c) for c in coeff])
        d2 = np.array(d2)
        
        d = s.isoft_many(coeff,use_mp=True)
        
        assert np.allclose(d2,d),'Pool using forked Soft instance is not the same as native computation.'
        
    
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
        res = np.squeeze(Soft.rotate_ylm_cmplx(y_dens_lm[None,...],euler_angles))
        
        y_out = [sh.synth_cplx(r) for r in res]
        y_argmax = tuple(np.argmax(o.real.mean(axis = 0)) for o in y_out)

        assert y_argmax == (0,64,128,192), 'Density maxima do not coincide with used righthanded rotations around z-axis.'
        
        z_dens = np.zeros((len(thetas),len(phis)),dtype = complex)
        delta = np.pi/64
        z_dens[(thetas<(np.pi/4+delta)) * (thetas>(np.pi/4-delta)) ,:]=1
        z_dens *= y_dens
        
        z_dens_lm = sh.analys_cplx(z_dens)
        euler_angles = np.array([[0,0,0],[0,np.pi/2,0],[0,np.pi,0],[0,np.pi*3/2,0]])
        res = np.squeeze(Soft.rotate_ylm_cmplx(z_dens_lm[None,...],euler_angles))
        
        z_out = [sh.synth_cplx(r) for r in res]
        phi_ids = tuple(np.argmax(o.real.mean(axis = 0)) for o in z_out)
        theta_ids = tuple(np.argmax(o.real.mean(axis = 1)) for o in z_out)
        assert (phi_ids == (0,0,128,128)) and (theta_ids == (31,96,96,31)), 'Density maxima do not coincide with right handed rotations around y-axis'

class TestCoeffSO3:
    def test_mem_leak(self):
        bw =32
        n_coeff = _soft.utils.total_num_coeffs(bw)
        lmns = _soft.utils.get_coeff_degrees(bw)
        c1 = create_complex128(n_coeff)
        c2 = soft.CoeffSO3(c1,lmns,'mnl')
        c1w = weakref.ref(c1)
        c2w = weakref.ref(c2)
    
        del c1
        del c2
        assert (c1w() is None) and (c2w() is None), 'CoeffSO3 memory leak detected, for "mnl" order '
        
        lmns = _soft.utils.get_coeff_degrees_risbo(bw)
        c1 = create_complex128(n_coeff)
        c2 = soft.CoeffSO3(c1,lmns,'mnl')
        c1w = weakref.ref(c1)
        c2w = weakref.ref(c2)
        del c1
        del c2
        assert (c1w() is None) and (c2w() is None), 'CoeffSO3 memory leak detected, for "lmn" order '

class TestIntegration:
    def test_CoeffSO3_same_as_manual_indexer(self):
        bw = 32
        for recurrence_type in Soft.recurrence_types:
            s = Soft(bw,recurrence_type = recurrence_type)
            c = s.get_coeff(random=True)
            assert np.allclose(c.lmn[4:10,4,:],c[s.coeff_indexer[4:10,4,:]]) , f'Rucurrence type {rucurrence_type}: Indexing mismatch between CoeffSO3.lmn[] and Soft.coeff_indexer[].'
        

        
