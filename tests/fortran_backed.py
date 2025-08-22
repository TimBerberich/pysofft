import numpy as np
import pysofft
from pysofft import _soft


class TestMakeWigner:
    bw = 15
    def test_size_wigner_d(self):
        bw = self.bw
        n_wigners = (bw**2 * (2 + 3*bw + bw**2))//3
        assert n_wigners == _soft.make_wigner.size_wigner_d(bw), 'Computed size mismatch of wigner_d'
        
    def test_create_beta_samples(self):
        N=10
        ks = np.arange(N)
        chebyshev_nodes = (2*ks+1)*np.pi/(2*N)
        chebyshev_nodes = np.array(chebyshev_nodes)
        print(chebyshev_nodes)
        assert np.allclose(chebyshev_nodes, _soft.make_wigner.create_beta_samples(N)), 'Chebyshev nodes are not correct'
        
    def test_genWig_L2(self):
        # See https://en.wikipedia.org/wiki/Wigner_D-matrix
        # First test uses that the small wigner_d matrices for m1=m2=0 are given by
        # Legendre polynomials ,i.e.:      d^l_{0,0}(\beta) = P_l(\cos(\beta))  
        
        Lmax = 128
        n_beta = 2*(Lmax+1)
        ks = np.arange(n_beta)
        betas = (2*ks+1)*np.pi/(2*n_beta)
        trig_samples = np.zeros((n_beta,3),order='F')
        trig_samples[:,0] = np.cos(betas)
        trig_samples[:,1] = np.cos(betas/2)
        trig_samples[:,2] = np.sin(betas/2)
        # Normalization parameter for wigner small d matrices used in SOFT 
        wig_norm = np.sqrt((2*np.arange(Lmax+1)+1)/2)
        
        P_l = np.array([np.polynomial.legendre.legval(trig_samples[:,0],[0]*l+[1]) for l in range(0,Lmax+1)]).T
        P_l *= wig_norm
        
        trig = _soft.make_wigner.create_trig_samples(Lmax+1)

        dlml = _soft.make_wigner.compute_dlml(0,0,trig[:,1],trig[:,2])
        d_l00 = _soft.make_wigner.genwig_l2(0,0,Lmax+1,trig[:,0],dlml)
        assert np.allclose(P_l,d_l00), 'Wigner small d differ from Legendre polynomial for m1=m2=0.'
        
        
        
        
class TestUtils:
    bw = 16
    def test_flat_to_triangular_index(self):
        N=self.bw-1
        k=0
        fi = np.array(1)
        fj = np.array(1)
        for i in range(N+1):
            for j in range(i,N+1):
                k+=1
                _soft.utils.flat_to_triangular_index(fi,fj,k,N)
                assert i==int(fi) and j==int(fj), f'Index mismatch between i,j={(i,j)} and k={k}, fi,fj={(fi,fj)}.'
                
    def test_flat_to_pyramid_index(self):
        N=self.bw-1
        k=0
        fi = np.array(1)
        fj = np.array(1)
        for i in range(N+1):
            for j in range(-i,i+1):
                k+=1
                _soft.utils.flat_to_pyramid_index(fi,fj,k)
                assert i==int(fi) and j==int(fj), f'Index mismatch between i,j={(i,j)} and k={k}, fi,fj={(fi,fj)}.'
