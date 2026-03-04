import numpy as np
import pysofft
from math import factorial
from pysofft import _soft
from pysofft.soft import rotate_ylm_cmplx,rotate_ylm_real

def compute_dlml_naiv(l,m,betas):
    # This exact formula is only good for small l !
    f = factorial
    dlml = np.sqrt((2*l+1)/2)*np.sqrt(f(2*l)/(f(l+m)*f(l-m)))*np.cos(betas/2)**(l+m)*np.sin(betas/2)**(l-m)
    return dlml

def wigner_Dlmn_limited_l(l,m,n,alpha,beta,gamma):
    if l==0:
        normalization = np.sqrt(1/2)
        d = np.ones(len(beta))
    elif l==1:
        normalization = np.sqrt((2*1+1)/2)
        if m==1 and n==1:
            d = np.cos(beta/2)**2
        elif m==1 and n==0:
            d = -np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==1 and n==-1:
            d = np.sin(beta/2)**2
        elif m==0 and n==1:
            d = np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==0 and n==0:
            d = np.cos(beta)
        elif m==0 and n==-1:
            d = -np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==-1 and n==1:
            d = np.sin(beta/2)**2
        elif m==-1 and n==0:
            d = np.sqrt(2)*np.sin(beta/2)*np.cos(beta/2)
        elif m==-1 and n==-1:
            d = np.cos(beta/2)**2
        else:
            raise Exception(f'm={m} and n={n} not possible for l=1')
    elif l==2 :
        normalization = np.sqrt((2*2+1)/2)
        if m==2 and n==2:
            d=np.cos(beta/2)**4
        elif m==2 and n==1:
            d=-2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==2 and n==0:
            d= np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==2 and n==-1:
            d= -2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==2 and n==-2:
            d=np.sin(beta/2)**4
        elif m==1 and n==2:
            d=2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==1 and n==1:
            d=1/2*np.cos(beta/2)**2*(4*np.cos(beta)-2)
        elif m==1 and n==0:
            d=-np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==1 and n==-1:
            d=1/2*np.sin(beta/2)**2*(4*np.cos(beta)+2)
        elif m==1 and n==-2:
            d=-2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==0 and n==2:
            d=np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==0 and n==1:
            d=np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==0 and n==0:
            d=1/2*(3*np.cos(beta)**2-1)
        elif m==0 and n==-1:
            d=-np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==0 and n==-2:
            d=np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==-1 and n==2:
            d=2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==-1 and n==1:
            d=1/2*np.sin(beta/2)**2*(4*np.cos(beta)+2)
        elif m==-1 and n==0:
            d=np.sqrt(6)*np.sin(beta/2)*np.cos(beta/2)*np.cos(beta)
        elif m==-1 and n==-1:
            d=1/2*np.cos(beta/2)**2*(4*np.cos(beta)-2)
        elif m==-1 and n==-2:
            d=-2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==-2 and n==2:
            d=np.sin(beta/2)**4
        elif m==-2 and n==1:
            d=2*np.sin(beta/2)**3*np.cos(beta/2)
        elif m==-2 and n==0:
            d=np.sqrt(6)*np.sin(beta/2)**2*np.cos(beta/2)**2
        elif m==-2 and n==-1:
            d=2*np.sin(beta/2)*np.cos(beta/2)**3
        elif m==-2 and n==-2:
            d=np.cos(beta/2)**4
        else:
            raise Exception(f'm={m} and n={n} not possible for l=2')                
    elif l>2:
        normalization = 1#np.sqrt((2*l+1)/2)
        if n==l:
            d = compute_dlml_naiv(l,m,beta)
        else:
            raise NotImplementedError()
    d *=normalization
    D = np.exp(-1.j*m*alpha[None,:,None])*d[:,None,None]*np.exp(-1.j*n*gamma[None,None,:])
    return D

def cos_2eta_coeff(l,eta):
    from math import lgamma
    if l<=(2*eta) and l%2==0:
        kcos = (2*eta+1)*np.sqrt((2*l+1)/2)
        ll=l//2
        return kcos*np.exp(lgamma(2*eta+1)+lgamma(eta+ll+1)+(l+1)*np.log(2)-lgamma(eta-l//2+1)-lgamma(2*eta+l+2))
    else:
        return 0
def cos_2eta_func(betas,eta):
    kcos = 2*eta+1
    return kcos*np.cos(betas)**(eta*2)

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
    def test_wigner_dl_risbo_reduced(self):
        Lmax = 128
        betas = np.array([np.pi*5/11,np.pi-np.pi*5/11])
        for l in range(2,Lmax+1):
            dl_red = _soft.make_wigner.wigner_dl_risbo_reduced(l,betas,False)
            dl = _soft.make_wigner.sym_reduced_to_full_wigner(dl_red,l)[0]
            assert np.allclose(dl@dl.T,np.eye(2*l+1)), f'wigner matrix of degree {l} is not orthogonal'
    def test_wigner_dl_comparative(self):
        Lmax = 128
        betas = np.array([np.pi*5/11,np.pi-np.pi*5/11])
        for l in range(2,Lmax+1):
            dl_risbo = _soft.make_wigner.wigner_dl_risbo(l,betas,False)[0]
            dl_tmp = _soft.make_wigner.wigner_dl_risbo_reduced(l,betas,False)
            dl_risbo_reduced = _soft.make_wigner.sym_reduced_to_full_wigner(dl_tmp,l)[0]
            dl_kostelec = _soft.make_wigner.wigner_dl_kostelec(l,betas,False)[0]
            assert np.allclose(dl_risbo,dl_kostelec), f'Missmatch between risbo and kostelec in wigner matrix of degree {l}.'
            assert np.allclose(dl_risbo_reduced,dl_kostelec), f'Missmatch between risbo reduced and kostelec in wigner matrix of degree {l}.' 
        
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
    def test_init(self):
        # exhaustive test init for segfaults or similar by trying out all obtions for bandwidth below a threshold.
        # Test is passed if the program does not crash.
        fft_flags=[0,64]
        init_ffts = [True,False]
        precompute_wigners = [True,False]
        bandwidth = np.arange(1,30)
        recurrence_types = [0,1]
        
        grid = np.meshgrid(*[bandwidth,precompute_wigners,init_ffts,recurrence_types,fft_flags])
        grid = [g.flatten() for g in grid]
        
        for i in np.arange(len(grid[0])):
            s_int=_soft.py.py_init_soft(grid[0][i],
                                        grid[0][i],
                                        grid[1][i],
                                        grid[2][i],
                                        grid[3][i],
                                        grid[4][i])
            _soft.py.py_destroy(s_int)
            
    # test_soft and test_isoft are the important tests all other tests
    # check transform properties and self consistency but these two check against analytic truth.
    def test_isoft(self):
        '''
        1. Test whether inverse soft of a single nonzero coefficient gives Wigner D matrices as result.
        2. Test using that the SO3 coefficients of (2*n+1)*cos(beta)**{2*n} are given by
            flmn nonzero only for m=n=0 and even l such that l<=2*n with coefficients
                (2*n+1)*(2**(l+1)*(2*n)!*(n+l/2+1)!)/((n-l/2)!*(2*n+l+2)!)
            see equation 199 in my thesis
            Tim Berberich
            Fluctuation X-ray Scattering for Systems of Particles Obeying Arbitrary Orientational Distributions: From Theory to Applications
            https://www.proquest.com/openview/ef531770a6c7bf2df0989399196adfa5/1?pq-origsite=gscholar&cbl=2026366&diss=y
        '''
        bw = 32
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,0,0)
        
        coeff = _soft.utils.get_empty_coeff(bw)
        beta = _soft.make_wigner.create_beta_samples(2*bw)
        alpha = _soft.make_wigner.create_alpha_gamma_samples(2*bw)
        gamma = alpha.copy()
        
        for l in range(10):
            for m in range(-l,l+1):
                if l<3:
                    for n in range(-l,l+1):
                        #print(f'l,n,m = {l,n,m}')
                        d_test = wigner_Dlmn_limited_l(l,m,n,alpha,beta,gamma)
                        coeff = _soft.utils.get_empty_coeff(bw)
                        d =  _soft.utils.get_empty_so3func_cmplx(bw)
                        d2 =  _soft.utils.get_empty_so3func_cmplx(bw)
                        cid = _soft.utils.coeff_location(m,n,l,bw)-1
                        coeff[cid]=1
                        _soft.py.py_isoft(s_int,coeff,d,False)
                        assert np.allclose(d_test,d), f'wigner mismatch for l,m,n,= {l,m,n}'
                        _soft.py.omp_set_num_threads_(4)
                        _soft.py.py_isoft(s_int,coeff,d2,True)
                        assert np.allclose(d_test,d2), f'wigner mismatch using OMP for l,m,n,= {l,m,n}'
                else:
                    n=l
                    #print(f'l,n,m = {l,n,m}')
                    d_test = wigner_Dlmn_limited_l(l,m,n,alpha,beta,gamma)
                    coeff = _soft.utils.get_empty_coeff(bw)
                    d =  _soft.utils.get_empty_so3func_cmplx(bw)
                    d2 =  _soft.utils.get_empty_so3func_cmplx(bw)
                    cid = _soft.utils.coeff_location(m,n,l,bw)-1
                    coeff[cid]=1
                    _soft.py.py_isoft(s_int,coeff,d,False)
                    assert np.allclose(d_test,d), f'wigner mismatch for l,m,n,= {l,m,n}'
                    _soft.py.omp_set_num_threads_(4)
                    _soft.py.py_isoft(s_int,coeff,d2,True)
                    assert np.allclose(d_test,d2), f'wigner mismatch using OMP for l,m,n,= {l,m,n}'
                    
        for eta in range(1,bw//2):
            d =  _soft.utils.get_empty_so3func_cmplx(bw)
            d2 =  _soft.utils.get_empty_so3func_cmplx(bw)
            coeff = _soft.utils.get_empty_coeff(bw)
            for l in range(bw):
                cid = _soft.utils.coeff_location(0,0,l,bw)-1
                coeff[cid]=cos_2eta_coeff(l,eta)
            _soft.py.py_isoft(s_int,coeff,d,False)
            d_test = cos_2eta_func(beta,eta)
            assert np.allclose(d_test,d[:,0,0]),f'isoft & cos(beta)^(2*eta) mismatch at eta={eta}'
            _soft.py.omp_set_num_threads_(4)
            _soft.py.py_isoft(s_int,coeff,d2,True)
            assert np.allclose(d_test,d2[:,0,0]),f'isoft OMP & cos(beta)^(2*eta) mismatch at eta={eta}'
        _soft.py.py_destroy(s_int)
        _soft.py.omp_set_num_threads_(1)        
    def test_soft(self):
        '''
        1. Test whether soft of a Wigner D matrics has a singe nonzero coefficient value as result.
        2. Test using that the SO3 coefficients of (2*n+1)*cos(beta)**{2*n} are given by
            flmn nonzero only for m=n=0 and even l such that l<=2*n with coefficients
                (2*n+1)*(2**(l+1)*(2*n)!*(n+l/2+1)!)/((n-l/2)!*(2*n+l+2)!)
            see equation 199 in my thesis
            Tim Berberich
            Fluctuation X-ray Scattering for Systems of Particles Obeying Arbitrary Orientational Distributions: From Theory to Applications
            https://www.proquest.com/openview/ef531770a6c7bf2df0989399196adfa5/1?pq-origsite=gscholar&cbl=2026366&diss=y
        '''
        bw = 32
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,0,0)
        
        coeff = _soft.utils.get_empty_coeff(bw)
        beta = _soft.make_wigner.create_beta_samples(2*bw)
        alpha = _soft.make_wigner.create_alpha_gamma_samples(2*bw)
        gamma = alpha.copy()
        
        for l in range(10):
            for m in range(-l,l+1):
                if l<3:
                    for n in range(-l,l+1):
                        #print(f'l,n,m = {l,n,m}')
                        coeff_test = _soft.utils.get_empty_coeff(bw)
                        cid = _soft.utils.coeff_location(m,n,l,bw)-1
                        coeff_test[cid]=1
                        coeff = _soft.utils.get_empty_coeff(bw)
                        coeff2 = _soft.utils.get_empty_coeff(bw)
                        d = wigner_Dlmn_limited_l(l,m,n,alpha,beta,gamma)
                        _soft.py.py_soft(s_int,d,coeff,False)
                        assert np.allclose(coeff_test,coeff), f'wigner mismatch for l,m,n,= {l,m,n}'
                        _soft.py.omp_set_num_threads_(4)
                        _soft.py.py_soft(s_int,d,coeff2,True)
                        assert np.allclose(coeff_test,coeff2), f'wigner mismatch OMP for l,m,n,= {l,m,n}'
                else:
                    n=l
                    #print(f'l,n,m = {l,n,m}')
                    coeff_test = _soft.utils.get_empty_coeff(bw)
                    cid = _soft.utils.coeff_location(m,n,l,bw)-1
                    coeff_test[cid]=1
                    coeff = _soft.utils.get_empty_coeff(bw)
                    d = wigner_Dlmn_limited_l(l,m,n,alpha,beta,gamma)
                    _soft.py.py_soft(s_int,d,coeff,False)
                    assert np.allclose(coeff_test,coeff), f'wigner mismatch for l,m,n,= {l,m,n}'
                    _soft.py.omp_set_num_threads_(4)
                    _soft.py.py_soft(s_int,d,coeff2,True)
                    assert np.allclose(coeff_test,coeff2), f'wigner mismatch OMP for l,m,n,= {l,m,n}'
                    
        for eta in range(1,bw//2):
            d =  _soft.utils.get_empty_so3func_cmplx(bw)
            d[:,:,:] =  cos_2eta_func(beta,eta)[:,None,None]
            coeff = _soft.utils.get_empty_coeff(bw)
            coeff_test = _soft.utils.get_empty_coeff(bw)
            for l in range(bw):
                cid = _soft.utils.coeff_location(0,0,l,bw)-1
                coeff_test[cid]=cos_2eta_coeff(l,eta)
            _soft.py.py_soft(s_int,d,coeff,False)
            assert np.allclose(coeff_test,coeff),f'isoft & cos(beta)^(2*eta) mismatch at eta={eta}'
            _soft.py.omp_set_num_threads_(4)
            _soft.py.py_soft(s_int,d,coeff2,True)
            assert np.allclose(coeff_test,coeff2),f'isoft OMP & cos(beta)^(2*eta) mismatch at eta={eta}'
        _soft.py.py_destroy(s_int)
        _soft.py.omp_set_num_threads_(1)
        
    # test invertability and equality between soft,isoft,rsoft,irsoft
    def test_soft_isoft_invertible(self):
        '''
        Tests invertability of the complex soft and isoft routines
        '''
        bw = 32
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,0,0)
        for recurrence_type in [0,1]:
            for bw in range(1,32):
                precompute_wigners = False
                s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
                coeff = _soft.utils.get_empty_coeff(bw)
                coeff2 = _soft.utils.get_empty_coeff(bw)
                d = _soft.utils.get_empty_so3func_cmplx(bw)
                d2 = _soft.utils.get_empty_so3func_cmplx(bw)
                
                coeff[...] = np.random.rand(*coeff.shape) + 1.j*np.random.rand(*coeff.shape)
                _soft.py.py_isoft(s_int,coeff,d,False)
                _soft.py.py_soft(s_int,d,coeff2,False)
                assert np.allclose(coeff,coeff2), f'isoft soft not identity for bw = {bw}'
                
                _soft.py.py_isoft(s_int,coeff2,d2,False)
                assert np.allclose(d,d2), f'soft isoft not identity for bw = {bw}'
                
                _soft.py.py_destroy(s_int)
        
    def test_rsoft_same_as_soft(self):
        '''
        Tests that the complex transform soft, restricted to real inputs, gives the same result as the real version rsoft.
        '''
        bw = 32
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,0,0)
        for recurrence_type in [0,1]:
            for bw in range(1,32):        
                s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
                d = _soft.utils.get_empty_so3func_cmplx(bw)            
                coeff = _soft.utils.get_empty_coeff(bw)
                coeff2 = _soft.utils.get_empty_coeff(bw)
                coeff3 = _soft.utils.get_empty_coeff(bw)
                d[...] = np.random.rand(*d.shape)
                _soft.py.py_soft(s_int,d,coeff,False)
                _soft.py.py_rsoft(s_int,d.real,coeff2,False)
                assert np.allclose(coeff,coeff2),f'rsoft,soft mismatch for bw={bw}'
                _soft.py.omp_set_num_threads_(4)
                _soft.py.py_rsoft(s_int,d.real,coeff3,True)
                assert np.allclose(coeff,coeff3),f'rsoft OMP,soft mismatch for bw={bw}'
                _soft.py.omp_set_num_threads_(1)
                _soft.py.py_destroy(s_int)            
    def test_irsoft_same_as_isoft(self):
        '''
        Tests that the complex transform soft, restricted to real inputs, gives the same result as the real version rsoft.
        '''
        bw = 32
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,0,0)
        for recurrence_type in [0,1]:
            for bw in range(1,32):        
                s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
                coeff = _soft.utils.get_empty_coeff(bw)
                d = _soft.utils.get_empty_so3func_cmplx(bw)
                d2 = _soft.utils.get_empty_so3func_real(bw)
                d3 = _soft.utils.get_empty_so3func_real(bw)            
                coeff[...] = np.random.rand(*coeff.shape) + 1.j*np.random.rand(*coeff.shape)
                if recurrence_type==0:
                    _soft.utils.enforce_real_sym(coeff,bw)
                else:
                    _soft.utils.enforce_real_sym_lmn(coeff,bw)
                
                _soft.py.py_isoft(s_int,coeff,d,False)
                _soft.py.py_irsoft(s_int,coeff,d2,False)
                assert np.allclose(d.real,d2),f'irsoft,isoft mismatch for bw={bw} and recurrence_type = {recurrence_type}'
                _soft.py.omp_set_num_threads_(4)
                _soft.py.py_irsoft(s_int,coeff,d3,True)
                assert np.allclose(d.real,d3),f'irsoft OMP,isoft  mismatch for bw={bw} and recurrence_type = {recurrence_type}'
                _soft.py.omp_set_num_threads_(1)
                _soft.py.py_destroy(s_int)         
    def test_rsoft_irsoft_invertible(self):
        '''
        Tests invertability of the real rsoft and irsoft routines
        '''
        bw = 32
        precompute_wigners = False
        for recurrence_type in [0,1]:
            for bw in range(1,32):        
                s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
                coeff = _soft.utils.get_empty_coeff(bw)
                coeff2 = _soft.utils.get_empty_coeff(bw)
                d = _soft.utils.get_empty_so3func_real(bw)
                d2 = _soft.utils.get_empty_so3func_real(bw)
                
                coeff[...] = np.random.rand(*coeff.shape) + 1.j*np.random.rand(*coeff.shape)
                if recurrence_type==0:
                    _soft.utils.enforce_real_sym(coeff,bw)
                else:
                    _soft.utils.enforce_real_sym_lmn(coeff,bw)
                _soft.py.py_irsoft(s_int,coeff,d,False)
                _soft.py.py_rsoft(s_int,d,coeff2,False)
                assert np.allclose(coeff,coeff2), f'isoft soft not identity for bw = {bw} and recurrence_type={recurrence_type}'
                
                _soft.py.py_irsoft(s_int,coeff2,d2,False)
                assert np.allclose(d,d2), f'soft isoft not identity for bw = {bw} and recurrence_type={recurrence_type}'
                
                _soft.py.py_destroy(s_int)                        

    # test _many versions of all transforms
    def test_soft_many_same_as_soft(self):
        '''
        Tests invertability of the real rsoft and irsoft routines
        '''
        bw = 32
        howmany=10
        precompute_wigners = False
        for recurrence_type in [0,1]:
            s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
            d = _soft.utils.get_empty_so3func_cmplx_many(bw,howmany)
            coeff = _soft.utils.get_empty_coeff_many(bw,howmany)
            coeff2 = _soft.utils.get_empty_coeff_many(bw,howmany)
            coeff3 = _soft.utils.get_empty_coeff_many(bw,howmany)
            d[:]=np.random.rand(*d.shape)
            for i in range(howmany):
                _soft.py.py_soft(s_int,d[...,i],coeff[:,i],False)
            _soft.py.py_soft_many(s_int,d,coeff2,False)
            assert np.allclose(coeff,coeff2),'Mismatch between soft and soft_many'
            _soft.py.omp_set_num_threads_(4)
            _soft.py.py_soft_many(s_int,d,coeff3,True)
            assert np.allclose(coeff,coeff3),'Mismatch between soft and soft_many using OMP'        
            _soft.py.py_destroy(s_int)
            _soft.py.omp_set_num_threads_(1)
    def test_rsoft_many_same_as_rsoft(self):
        '''
        Tests invertability of the real rsoft and irsoft routines
        '''
        bw = 31
        howmany=10
        precompute_wigners = False
        for recurrence_type in [0,1]:
            s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
            d = _soft.utils.get_empty_so3func_real_many(bw,howmany)
            d2 = _soft.utils.get_empty_so3func_real_many(bw,howmany)
            d3 = _soft.utils.get_empty_so3func_real_many(bw,howmany)
            coeff = _soft.utils.get_empty_coeff_many(bw,howmany)
            coeff2 = _soft.utils.get_empty_coeff_many(bw,howmany)
            coeff3 = _soft.utils.get_empty_coeff_many(bw,howmany)
            d[:]=np.random.rand(*d.shape)
            d2[:] = d
            d3[:] = d
            for i in range(howmany):
                _soft.py.py_rsoft(s_int,d[...,i],coeff[:,i],False)
            _soft.py.py_rsoft_many(s_int,d2,coeff2,False)
            assert np.allclose(coeff,coeff2),'Mismatch between rsoft and rsoft_many'
            _soft.py.omp_set_num_threads_(4)
            _soft.py.py_rsoft_many(s_int,d2,coeff3,True)
            assert np.allclose(coeff,coeff3),'Mismatch between rsoft and rsoft_many using OMP'        
            _soft.py.py_destroy(s_int)
            _soft.py.omp_set_num_threads_(1)
            
    def test_isoft_many_same_as_isoft(self):
        '''
        Tests invertability of the real rsoft and irsoft routines
        '''
        bw = 33
        howmany=9
        precompute_wigners = False
        for recurrence_type in [0,1]:
            s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
            coeff = _soft.utils.get_empty_coeff_many(bw,howmany)
            d = _soft.utils.get_empty_so3func_cmplx_many(bw,howmany)
            d2 = _soft.utils.get_empty_so3func_cmplx_many(bw,howmany)
            d3 = _soft.utils.get_empty_so3func_cmplx_many(bw,howmany)
            coeff[:]=np.random.rand(*coeff.shape)+1.j*np.random.rand(*coeff.shape)
            for i in range(howmany):
                _soft.py.py_isoft(s_int,coeff[:,i],d[...,i],False)
            _soft.py.py_isoft_many(s_int,coeff,d2,False)
            assert np.allclose(d,d2),'Mismatch between rsoft and rsoft_many'
            _soft.py.omp_set_num_threads_(4)
            _soft.py.py_isoft_many(s_int,coeff,d3,True)
            assert np.allclose(d,d3),'Mismatch between isoft and isoft_many using OMP'        
            _soft.py.py_destroy(s_int)
            _soft.py.omp_set_num_threads_(1)
    def test_irsoft_many_same_as_irsoft(self):
        '''
        Tests invertability of the real rsoft and irsoft routines
        '''
        bw = 27
        howmany=11
        precompute_wigners = False
        for recurrence_type in [0,1]:
            s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,True,recurrence_type,0)
            coeff = _soft.utils.get_empty_coeff_many(bw,howmany)
            coeff2 = _soft.utils.get_empty_coeff_many(bw,howmany)
            d = _soft.utils.get_empty_so3func_real_many(bw,howmany)        
            d2 = _soft.utils.get_empty_so3func_real_many(bw,howmany)
            d3 = _soft.utils.get_empty_so3func_real_many(bw,howmany)
            coeff[:]=np.random.rand(*coeff.shape)+1.j*np.random.rand(*coeff.shape)
            _soft.utils.enforce_real_sym(coeff,bw)
            coeff2[:] = coeff
            for i in range(howmany):
                _soft.py.py_irsoft(s_int,coeff[:,i],d[...,i],False)
            _soft.py.py_irsoft_many(s_int,coeff2,d2,False)
            #assert np.allclose(d,d2),'Mismatch between rsoft and rsoft_many'
            #assert((d3==0).all())
            _soft.py.omp_set_num_threads_(4)
            _soft.py.py_irsoft_many(s_int,coeff2,d3,True)
            assert np.allclose(d,d3),'Mismatch between irsoft and irsoft_many using OMP'        
            _soft.py.py_destroy(s_int)
            _soft.py.omp_set_num_threads_(1)
        
    def test_integrate_over_so3(self):
        r'''Tests that integral of a constant 1 functions returns the volume of SO(3),
        which is $$8 \pi^2$$ in the Haar measure.
        '''
        bw = 64
        precompute_wigners = False
        s_int = _soft.py.py_init_soft(bw,bw-1,precompute_wigners,False,0,0)
        dr = _soft.utils.get_empty_so3func_real(bw)
        dc = _soft.utils.get_empty_so3func_cmplx(bw)
        dr[:]=1
        dc[:]=1
        
        vr = _soft.py.py_integrate_over_so3_real(s_int,dr)
        vc = _soft.py.py_integrate_over_so3_cmplx(s_int,dc)
        assert np.isclose(vr,8*np.pi**2), f'Incorrect real integral, value is {vr}'
        assert np.isclose(vc.real,8*np.pi**2), f'Incorrect real integral, value is {vc}'

    # test correlation coputations 
    def test_cross_correlate_cmplx(self):
        '''
        Complex case
        Makes sure that cross_correlating two, with respect to each other, rotated spherical harmonic coefficients
        gives a maximal value at precisely their relative rotation.
        '''
        bw = 32
        coeff = np.zeros(_soft.utils.n_lmc(bw),dtype=complex)
        coeff[...] = np.random.rand(*coeff.shape)+1.j*np.random.rand(*coeff.shape)
        betas = _soft.make_wigner.create_beta_samples(2*bw)
        albe = _soft.make_wigner.create_alpha_gamma_samples(2*bw)
        
        for precomputed_wigners in [False,True]:
            for recurrence_type in [0,1]:
                s_int = _soft.py.py_init_soft(bw,bw-1,precomputed_wigners,True,recurrence_type,0)
                for multiprocessing in [False,True]:                
                    for i in range(5):
                        rotation_index = (np.random.rand(3)*len(albe)).astype(int)
                        eulers = (albe[rotation_index[0]],betas[rotation_index[1]],albe[rotation_index[2]])
                        rot_coeff = rotate_ylm_cmplx(coeff[None,...],eulers)
                        so_coeff= _soft.utils.get_empty_coeff(bw).copy()
                        corr = _soft.utils.get_empty_so3func_cmplx(bw)
                        _soft.py.py_cross_correlation_ylm_cmplx(s_int,coeff,rot_coeff,corr,multiprocessing)
                        found_rotation_index = np.array(np.unravel_index(np.argmax(np.abs(corr)),corr.shape))
                        tmp = found_rotation_index[0]
                        found_rotation_index[0] = found_rotation_index[1]
                        found_rotation_index[1] = tmp
                        assert np.all(rotation_index == found_rotation_index), f'Multiprocessing={multiprocessing},precomputed_wigners = {precomputed_wigners},recurrence_type{recurrence_type}: True rotation_index = {rotation_index} but found was {found_rotation_index}.'
    def test_cross_correlate_real(self):
        '''
        Real case
        Makes sure that cross_correlating two, with respect to each other, rotated spherical harmonic coefficients
        gives a maximal value at precisely their relative rotation.
        '''
        bw = 32
        coeff = np.zeros(_soft.utils.n_lmr(bw),dtype=complex)
        coeff[...] = np.random.rand(*coeff.shape)
        betas = _soft.make_wigner.create_beta_samples(2*bw)
        albe = _soft.make_wigner.create_alpha_gamma_samples(2*bw)
        
        for precomputed_wigners in [False,True]:
            for recurrence_type in [0,1]:
                s_int = _soft.py.py_init_soft(bw,bw-1,precomputed_wigners,True,recurrence_type,0)
                for multiprocessing in [False,True]:                
                    for i in range(5):
                        rotation_index = (np.random.rand(3)*len(albe)).astype(int)
                        eulers = (albe[rotation_index[0]],betas[rotation_index[1]],albe[rotation_index[2]])
                        rot_coeff = rotate_ylm_real(coeff[None,...],eulers)
                        so_coeff= _soft.utils.get_empty_coeff(bw).copy()
                        corr = _soft.utils.get_empty_so3func_real(bw)
                        _soft.py.py_cross_correlation_ylm_real(s_int,coeff,rot_coeff,corr,multiprocessing)
                        found_rotation_index = np.array(np.unravel_index(np.argmax(np.abs(corr)),corr.shape))
                        tmp = found_rotation_index[0]
                        found_rotation_index[0] = found_rotation_index[1]
                        found_rotation_index[1] = tmp
                        assert np.all(rotation_index == found_rotation_index), f'Multiprocessing={multiprocessing},precomputed_wigners = {precomputed_wigners},recurrence_type{recurrence_type}: True rotation_index = {rotation_index} but found was {found_rotation_index}.'
    def test_cross_correlate_cmplx_3d(self):
        '''
        Complex case
        Makes sure that cross_correlating two, with respect to each other, rotated spherical harmonic coefficients
        gives a maximal value at precisely their relative rotation.
        '''
        bw = 32
        coeff = np.zeros((33,_soft.utils.n_lmc(bw)),dtype=complex)
        coeff[...] = np.random.rand(*coeff.shape)
        betas = _soft.make_wigner.create_beta_samples(2*bw)
        albe = _soft.make_wigner.create_alpha_gamma_samples(2*bw)
        
        for precomputed_wigners in [False,True]:
            for recurrence_type in [0,1]:
                s_int = _soft.py.py_init_soft(bw,bw-1,precomputed_wigners,True,recurrence_type,0)
                for multiprocessing in [False,True]:                
                    for i in range(5):
                        rotation_index = (np.random.rand(3)*len(albe)).astype(int)
                        eulers = (albe[rotation_index[0]],betas[rotation_index[1]],albe[rotation_index[2]])
                        rot_coeff = rotate_ylm_cmplx(coeff,eulers).T
                        coefft = coeff.T
                        rad_points = np.arange(coefft.shape[1]).astype(float)
                        rad_lim = np.array((1,coefft.shape[1]))
                    
                        corr = _soft.utils.get_empty_so3func_cmplx(bw)
                        _soft.py.py_cross_correlation_ylm_cmplx_3d(s_int,coefft,rot_coeff,corr,rad_points,rad_lim,multiprocessing)
                        found_rotation_index = np.array(np.unravel_index(np.argmax(np.abs(corr)),corr.shape))
                        tmp = found_rotation_index[0]
                        found_rotation_index[0] = found_rotation_index[1]
                        found_rotation_index[1] = tmp
                        assert np.all(rotation_index == found_rotation_index), f'Multiprocessing={multiprocessing},precomputed_wigners = {precomputed_wigners},recurrence_type{recurrence_type}: True rotation_index = {rotation_index} but found was {found_rotation_index}.'
    def test_cross_correlate_real_3d(self):
        '''
        Complex case
        Makes sure that cross_correlating two, with respect to each other, rotated spherical harmonic coefficients
        gives a maximal value at precisely their relative rotation.
        '''
        bw = 32
        coeff = np.zeros((33,_soft.utils.n_lmr(bw)),dtype=complex)
        coeff[...] = np.random.rand(*coeff.shape)
        betas = _soft.make_wigner.create_beta_samples(2*bw)
        albe = _soft.make_wigner.create_alpha_gamma_samples(2*bw)
        
        for precomputed_wigners in [False,True]:
            for recurrence_type in [0,1]:
                s_int = _soft.py.py_init_soft(bw,bw-1,precomputed_wigners,True,recurrence_type,0)
                for multiprocessing in [False,True]:                
                    for i in range(5):
                        rotation_index = (np.random.rand(3)*len(albe)).astype(int)
                        eulers = (albe[rotation_index[0]],betas[rotation_index[1]],albe[rotation_index[2]])
                        rot_coeff = rotate_ylm_real(coeff,eulers).T
                        coefft = coeff.T
                        rad_points = np.arange(coefft.shape[1]).astype(float)
                        rad_lim = np.array((1,coefft.shape[1]))
                    
                        corr = _soft.utils.get_empty_so3func_real(bw)
                        _soft.py.py_cross_correlation_ylm_real_3d(s_int,coefft,rot_coeff,corr,rad_points,rad_lim,multiprocessing)
                        found_rotation_index = np.array(np.unravel_index(np.argmax(np.abs(corr)),corr.shape))
                        tmp = found_rotation_index[0]
                        found_rotation_index[0] = found_rotation_index[1]
                        found_rotation_index[1] = tmp
                        assert np.all(rotation_index == found_rotation_index), f'Multiprocessing={multiprocessing},precomputed_wigners = {precomputed_wigners},recurrence_type{recurrence_type}: True rotation_index = {rotation_index} but found was {found_rotation_index}.'
        
