import numpy as np
import pysofft
from math import factorial
from pysofft import _soft

def compute_dlml_naiv(l,m,betas):
    # This exact formula is only good for small l !
    f = factorial
    dlml = np.sqrt((2*l+1)/2)*np.sqrt(f(2*l)/(f(l+m)*f(l-m)))*np.cos(betas/2)**(l+m)*np.sin(betas/2)**(l-m)
    return dlml
def wigner_Dlmn_limited_l(l,m,n,alpha,beta,gamma):
    if l==0:
        d = 1
    elif l==1:
        normalization = np.sqrt((2*1+1)/2)
        if m==-1 and n==-1:
            d = np.cos(beta/2)**2
        elif m==-1 and n==0:
            d = -np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==-1 and n==1:
            d[0,0]= np.sin(beta/2)**2
        elif m==0 and n==-1:
            d = np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==0 and n==0:
            d = np.cos(beta)
        elif m==0 and n==1:
            d = -np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==1 and n==-1:
            d = np.sin(beta/2)**2
        elif m==1 and n==0:
            d = np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==1 and n==1:
            d = np.cos(beta/2)**2
        else:
            raise Exception(f'm={m} and n={n} not possible for l=1')
    elif l==2 :
        normalization = np.sqrt((2*2+1)/2)
        if m==-2 and n==-2:
            d=np.cos(beta/2)**2
        elif m==-2 and n==-1:
            d=-2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==-2 and n==0:
            d= np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==-2 and n==1:
            d= -2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==-2 and n==2:
            d=np.sin(beta/2)**4
        elif m==-1 and n==-2:
            d=2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==-1 and n==-1:
            d=1/2*np.cos(beta/2)**2*(4*np.cos(beta)-2)
        elif m==-1 and n==0:
            d=np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==-1 and n==1:
            d=1/2*np.sin(beta/2)**2*(4*np.cos(beta)+2)
        elif m==-1 and n==2:
            d=-2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==0 and n==-2:
            d=np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==0 and n==-1:
            d=np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==0 and n==0:
            d=1/2*(3*np.cos(beta)-1)
        elif m==0 and n==1:
            d=-np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==0 and n==2:
            d=np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==1 and n==-2:
            d=2*np.sin(beta/2)**2*np.cos(beta/2)
        elif m==1 and n==-1:
            d=1/2*np.sin(beta/2)**2*(4*np.cos(beta)+2)
        elif m==1 and n==0:
            d=np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==1 and n==1:
            d=1/2*np.cos(beta/2)**2*(4*np.cos(beta)-2)
        elif m==1 and n==2:
            d=-2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==2 and n==-2:
            d=np.sin(beta/2)**4
        elif m==2 and n==-1:
            d=2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==2 and n==0:
            d=np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==2 and n==1:
            d=2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==2 and n==2:
            d=np.cos(beta/2)**4
        else:
            raise Exception(f'm={m} and n={n} not possible for l=2')                
    elif l>2:
        normalization = np.sqrt((2*l+1)/2)
        if n==l:
            d = compute_dlml_naiv(l,m,beta)
        else:
            raise NotImplementedError()
        d *=normalization
        D = np.exp(-1.j*m*alpha)*d*np.exp(-1.j*n*gamma)
    return D
    
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

    def test_compute_dlml(self):
        # Makes sure that the computed dlml agree with its exact formula.
        Lmax = 10
        n_beta = 2*(Lmax+1)
        ks = np.arange(n_beta)
        betas = (2*ks+1)*np.pi/(2*n_beta)
        trigs = np.zeros((n_beta,3),order='F')
        trigs[:,0] = np.cos(betas)
        trigs[:,1] = np.sin(betas/2)*np.cos(betas/2)
        trigs[:,2] = np.cos(betas/2)**2
        for m in range(0,Lmax):
            for l in range(m,Lmax):
                dlml = compute_dlml_naiv(l,m,betas)
                assert np.allclose(_soft.make_wigner.compute_dlml(l,m,trigs[:,1],trigs[:,2]), dlml), f'd^{l}_({m},{l}) does not coincide with its exact formula!'
                
    def test_compute_all_dlml_l_contiguous(self):
        # Makes sure that the computed dlml agree with its exact formula.
        Lmax = 10
        n_beta = 2*(Lmax+1)
        ks = np.arange(n_beta)
        betas = (2*ks+1)*np.pi/(2*n_beta)
        trigs = np.zeros((n_beta,3),order='F')
        trigs[:,0] = np.cos(betas)
        trigs[:,1] = np.sin(betas/2)*np.cos(betas/2)
        trigs[:,2] = np.cos(betas/2)**2
        dlmls = _soft.make_wigner.compute_all_dlml_l_contiguous(Lmax+1,trigs[:,1],trigs[:,2],True)
        for m in range(Lmax+1):
            start = _soft.utils.triangular_to_flat_index(m,m,Lmax+1)-1
            stop = _soft.utils.triangular_to_flat_index(m,Lmax,Lmax+1)
            dlml = dlmls[:,start:stop]
            dlml_naiv = np.array([compute_dlml_naiv(l,m,betas) for l in range(m,Lmax+1)]).T
            assert np.allclose(dlml,dlml_naiv), f'dlml are not correct at m={m}'            
    def test_compute_all_dlml_m_contiguous(self):
        # Makes sure that the computed dlml agree with its exact formula.
        Lmax = 10
        n_beta = 2*(Lmax+1)
        ks = np.arange(n_beta)
        betas = (2*ks+1)*np.pi/(2*n_beta)
        trigs = np.zeros((n_beta,3),order='F')
        trigs[:,0] = np.cos(betas)
        trigs[:,1] = np.sin(betas/2)*np.cos(betas/2)
        trigs[:,2] = np.cos(betas/2)**2

        dlmls = _soft.make_wigner.compute_all_dlml_m_contiguous(Lmax+1,trigs[:,1],trigs[:,2],True)
        for l in range(Lmax+1):
            start = _soft.utils.triangular_to_flat_index_reversed(l,0)-1
            stop = _soft.utils.triangular_to_flat_index_reversed(l,l)
            dlml = dlmls[:,start:stop]
            dlml_naiv = np.array([compute_dlml_naiv(l,m,betas) for m in range(0,l+1)]).T
            assert np.allclose(dlml,dlml_naiv), f'dlml are not correct at m={l}'


    def _get_dl0n(self,l,n,betas):
        from scipy.special import assoc_legendre_p
        dl0n = assoc_legendre_p(l,n,np.cos(betas),norm=True)[0]
        # The (-1)**n is due to P^l_m = d^l_{m,0} = (-1)**(m-0)*d^l_{0,m}
        dl0n *= (-1)**(n)
        return dl0n
    def test_genwig_l2(self):
        # See https://en.wikipedia.org/wiki/Wigner_D-matrix
        # First test uses that the small wigner_d matrices for m1=m2=0 are given by
        # Legendre polynomials ,i.e.:      d^l_{0,0}(\beta) = P_l(\cos(\beta))
        # Second test uses that the small wigner_d matrices for m1=0 are given by
        # Associated Legendre polynomials ,i.e.:      d^l_{0,n}(\beta) = P^{|n|}_l(\cos(\beta))*(-1)^m  
        
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
        
        for n in range(1,Lmax+1):
            dlml = _soft.make_wigner.compute_dlml(n,0,trig[:,1],trig[:,2])
            d_l0n = _soft.make_wigner.genwig_l2(0,n,Lmax+1,trig[:,0],dlml)
            Plm = np.array([self._get_dl0n(l,n,betas) for l in range(n,Lmax+1)]).T
            assert np.allclose(d_l0n,Plm), f'Wigner small d differ from asssociated Legendre polynomial for m1=0 m2={n}.'            
    def test_genwig_l2_trsp(self):
        # genwig_l2_trsp is supposed to compute the transpose of genwig_l2
        Lmax = 128
        m1=32
        m2=56
        trig = _soft.make_wigner.create_trig_samples(Lmax+1)
        dlml = _soft.make_wigner.compute_dlml(Lmax,m1,trig[:,1],trig[:,2])
        d_lmn = _soft.make_wigner.genwig_l2(m1,m2,Lmax+1,trig[:,0],dlml)
        d_lmn_trsp = _soft.make_wigner.genwig_l2_trsp(m1,m2,Lmax+1,trig[:,0],dlml)
        assert np.allclose(d_lmn,d_lmn_trsp.T), f'genwig_l2_trsp is not the transpose of genwig_l2'
        
    def test_wigner_dl_kostelec(self):
        # makes sure that the computed small wigner d matrices are elements of  O(2l+1)
        Lmax = 128
        betas = np.array([np.pi*5/11])
        for l in range(2,Lmax+1):
            dl = _soft.make_wigner.wigner_dl_kostelec(l,betas,False)[0]
            assert np.allclose(dl@dl.T,np.eye(2*l+1)), f'wigner matrix of degree {l} is not orthogonal'
    def test_wigner_dl_risbo(self):
        # makes sure that the computed small wigner d matrices are elements of  O(2l+1)
        Lmax = 128
        betas = np.array([np.pi*5/11])
        for l in range(2,Lmax+1):
            dl = _soft.make_wigner.wigner_dl_risbo(l,betas,False)[0]
            assert np.allclose(dl@dl.T,np.eye(2*l+1)), f'wigner matrix of degree {l} is not orthogonal'    
    def test_wigner_dl_comparative(self):
        Lmax = 128
        betas = np.array([np.pi*5/11])
        for l in range(2,Lmax+1):
            dl_risbo = _soft.make_wigner.wigner_dl_risbo(l,betas,False)[0]
            dl_kostelec = _soft.make_wigner.wigner_dl_kostelec(l,betas,False)[0]
            assert np.allclose(dl_risbo,dl_kostelec), f'Missmatch between risbo an kosteleic in wigner matrix of degree {l}.'
        
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
        bandwidth = np.arange(1,30)
        
        grid = np.meshgrid(*[bandwidth,precompute_wigners,init_ffts,fft_flags])
        grid = [g.flatten() for g in grid]
        
        for i in np.arange(len(grid[0])):
            s_int=_soft.py.py_init_soft(grid[0][i],
                                        grid[0][i],
                                        grid[1][i],
                                        grid[2][i],
                                        grid[3][i])
            _soft.py.py_destroy(s_int)
            
    def test_isoft(self):
        bw = 64
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw,precompute_wigners,True,0)

        coeff = _soft.utils.get_empty_coeff(bw)
        
