import numpy as np
import pysofft
from math import factorial
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

    def compute_dlml_naiv(self,l,m,betas):
        # This exact formula is only good for small l !
        f = factorial
        dlml = np.sqrt((2*l+1)/2)*np.sqrt(f(2*l)/(f(l+m)*f(l-m)))*np.cos(betas/2)**(l+m)*np.sin(betas/2)**(l-m)
        return dlml
    
    def test_compute_dlml(self):
        Lmax = 4
        n_beta = 2*(Lmax+1)
        ks = np.arange(n_beta)
        betas = (2*ks+1)*np.pi/(2*n_beta)
        trigs = np.zeros((n_beta,3),order='F')
        trigs[:,0] = np.cos(betas)
        trigs[:,1] = np.sin(betas/2)*np.cos(betas/2)
        trigs[:,2] = np.cos(betas/2)**2
        for m in range(0,Lmax):
            for l in range(m,Lmax):
                dlml = self.compute_dlml_naiv(l,m,betas)
                assert np.allclose(_soft.make_wigner.compute_dlml(l,m,trigs[:,1],trigs[:,2]), dlml), f'd^{l}_({m},{l}) does not coincide with its exact formula!'
    def test_compute_all_dlml_l_contiguous(self):
        # compare to naiv computation for small l
        pass
    def test_compute_all_dlml_l_contiguous(self):
        # compare to naiv computation for small l
        pass
        
    def test_genwig_l2(self):
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
        
    def test_genwigall(self):
        # test Orthogonality with legendre weights
        pass
    def test_genwigall_preallocated(self):
        # compare to test_genwigall
        pass
    def test_wigner_l(self):
        # test Orthogonality of dlmn as m,n matrix
        pass
        
        
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

class TestSo3ft:    
    # Maybe test wigner transform separately
    
    # test forward and inverse for small l with known solutions soft and isoft
    
    # make sure real,alloc,nthreads>1 and _many versions yield the same results as the simple complex routines soft and isoft
    # dont forget lmax setting

    # test correlation with simple example
    # make sure the real and nthread versions

    # test ylm rotation for simple exact example

    # test Fork safety
    
    def test_init(self):
        # exhaustive test init for segfaults or similar by trying out all obtions for bandwidth below a threshold.
        # Test is passed if the program does not crash.
        fft_flags=[0,64]
        init_ffts = [True,False]
        precompute_wigners = [True,False]
        n_threads = [1,2]
        bandwidth = np.arange(30)

        grid = np.meshgrid(*[bandwidth,n_threads,precompute_wigners,init_ffts,fft_flags])
        grid = [g.flatten() for g in grid]

        for i in np.arange(len(grid[0])):
            _soft.so3ft.init(grid[0][i],
                             grid[0][i],
                             grid[1][i],
                             grid[2][i],
                             grid[3][i],
                             grid[4][i]
                             )


