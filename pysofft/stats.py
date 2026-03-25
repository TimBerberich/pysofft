import numpy as np
from pysofft import _soft,Soft,wigner,utils
try:
    from scipy.stats import rv_discrete
except ImportError as e:
    raise ImportError(
        "pysofft.stats requires the optional 'stats' dependencies. "
        "Install them with: pip install 'pysofft[stats]'"
    ) from e

########################################
## Probability distributions on SO(3) ##
class CharView:
    """
    Helper class that allows to acces a characteristic function $M_{l,n,k}$ by its indices l,n,k.
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
    def __init__(self,char_array,lnks):
        self.array = char_array
        self._lnks = lnks
        self.ls = np.sort(np.unique(lnks[:,0]))
        bw = len(self.ls)
        self.ns = np.roll(np.sort(np.unique(lnks[:,1])),bw)
        self.ks = np.roll(np.sort(np.unique(lnks[:,2])),bw)
    def _get_mask(self,items):
        if not isinstance(items,tuple):
            items = (items,)
        assert len(items)<=3
        selection_ids=[ids[sel] for sel,ids in zip(items,[self.ls,self.ns,self.ks])]
        masks = [np.isin(lnk,sids) for sids,lnk in zip(selection_ids,self._lnks.T)]
        mask = np.prod(masks,axis = 0,dtype=bool)
        return mask
    def __getitem__(self, items):        
        mask = self._get_mask(items)
        return self.array[mask]
    def __setitem__(self,items,value):
        m = self._get_mask(items)
        self.array[m] = value

class CharFuncSO3(np.ndarray):
    r"""
    A class to represent Charactersitic functions $M_{lnk}$ of distributions $\rho(\alpha,\beta,\gamma)$ on SO(3).
    $$M_{lnk} = \int_{SO(3)} \rho(\alpha,\beta,\gamma) D^l_{nk}(\alpha,\beta,\gamma) sin(\beta)\,d\alpha d\beta d\gamma$$
    where $D^l_{nk}(\alpha,\beta,\gamma)$ are Wigner-D matrices.

    Acts as a normal numpy array with additional attributes
    
    Attributes
    ----------
    bw : int
    Bandwidth of the distribution.
    soft : Soft
    Fourier transform instance on SO(3)
    Methods
    lnk : CharView
    allows to acces $M_{lnk}$ by its indices l,n,k
    distrib : np.ndarray, (2*bw,)*3, complex
    property that calculates the SO3 distribution by inverse fourier transform of the current characteristic function.
    """
    def __new__(cls,bw : int,char_func=None,density=None):
        s = Soft(bw)
        _density_shape = (2*bw)*3
        char = s.get_coeff(raw = True)
        char[0]=1
        if char_func is not None:
            assert char_func.shape == char.shape, f'specified characteristic function has wrong shape, expected {char.shape} given shape {char_func.shape}'
            char = char_func
        elif density is not None:
            assert density.shape == _density_shape, f'specified density has wrong shape, expected {_density_shape} given shape {density.shape}'
            char[:]= s.soft(density)

        char = char.view(cls)

        char._soft = s
        char._density_shape=_density_shape
        
        lnks = s.coeff_indices
        char.lnk = CharView(char,lnks)
        char.so3_const = np.sqrt(8*np.pi**2) # sqrt of SO(3) volume.
        char.so3_grid = np.stack(np.meshgrid(*s.euler_angles.values(),indexing='ij'),3) #beta,gamma,alpha
        char.legendre_weights = _soft.utils.legendre_quadrature_weights(s.bw) 
        return char
    @property
    def distrib(self):
        return self._soft.isoft(self)/self.so3_const
    @property
    def pmf(self):
        density = self.distrib
        N = len(self.so3_grid)
        pmf = (density.real*(2*np.pi/N)**2*self.legendre_weights[None,None,:])
        min_pmf = pmf.min()
        if min_pmf<0:
            pmf[pmf<np.abs(min_pmf)]=0
            pmf/=np.sum(pmf)
        return pmf
        
    @property
    def integrate_distrib(self):
        distrib = self._soft.inverse_cmplx(self)/self.so3_const
        return self._soft.integrate_over_so3(distrib)
    
    def sample(self,size = 100):
        pmf = self.pmf.copy()
        flat_pmf = pmf.reshape(-1)
        vals = np.arange(len(flat_pmf))
        flat_pos_pmf = flat_pmf[flat_pmf>0]
        pos_vals = vals[flat_pmf>0]
        rv = rv_discrete(name='so3_delta',values = (pos_vals,flat_pos_pmf))

        samples = rv.rvs(size = size)
        abg_samples = self.so3_grid.reshape(-1,3)[samples,:]
        return abg_samples

    def get_sampling_func(self):
        pmf = self.pmf.copy()
        flat_pmf = pmf.reshape(-1)
        vals = np.arange(len(flat_pmf))
        flat_pos_pmf = flat_pmf[flat_pmf>0]
        pos_vals = vals[flat_pmf>0]
        rv = rv_discrete(name='so3_delta',values = (pos_vals,flat_pos_pmf))
        flat_so3_grid = self.so3_grid.reshape(-1,3)

        def sample_func(size = 100):
            samples = rv.rvs(size = size)
            abg_samples = flat_so3_grid[samples,:].copy()
            return abg_samples
        return sample_func
        
        
    def mean(self,so3_function):
        '''
        Computes the average of a scalar/ndarray valued function on SO3 using the current probability distribution.
        Expects the function values to be given in the shape shape (M1,...,Mn,2*bw,2*bw,2*bw) and (2*bw,2*bw,2*bw) for scalar functions.
        Where bw is the bandwidth of the this probability distribution i.e. self._soft.bw and M1,...,Mn are the value array dimensions)
        '''
        bw=self._soft.bw
        shape = so3_function.shape        
        so3_function = so3_function.reshape(-1,2*bw,2*bw,2*bw)
        
        mean = np.zeros(so3_function.shape[0],dtype=complex)
        for _id in range(len(mean)):
            dconj_lnk = self._soft.forward_cmplx(so3_function[_id].conj())
            mean[_id] = np.sum(dconj_lnk*self)
        if len(shape)>3:
            mean=mean.reshape(shape[:-3])
        else:
            mean=mean[0]            
        return mean/self.so3_const

    def mean_density(self,density):
        '''Expects the density to be of shape (...,Nr,Ntheta,Nphi) where the last three dimensions correspond to the spherical coordinates and ... can be any shape'''
        bw = self._soft.bw
        shape = density.shape
        density = density.reshape(-1,shape[-3],shape[-2],shape[-1])
        Nr,Ntheta,Nphi = shape[-3:]
        sh = get_spherical_harmonic_transform_obj(bw-1,'complex',n_phi=Nphi,n_theta=Ntheta)
        sh_coeff = np.zeros((density.shape[0],Nr,bw**2),dtype=complex)
        for _id in range(len(density)):
            sh_coeff[_id] = sh._forward_direct(density[_id])
        mean_sh_coeff = self.mean_spherical_coeff(sh_coeff)
        mean_density=np.zeros(density.shape,dtype=complex)
        for _id in range(len(density)):
            mean_density[_id] = sh._inverse_direct(mean_sh_coeff[_id])

        if len(shape)>3:
            mean_density = mean_density.reshape(shape)
        else:
            mean_density = mean_density[0]
        return mean_density
        

    def mean_spherical_coeff(self,spherical_harmonic_coeff):
        '''
        Expects input coeff to be stored in a 1d array with the l,m's coefficient beeing at location $ l(l+1)+m $.
        '''
        bw = self._soft.bw        
        shape = spherical_harmonic_coeff.shape
        spherical_harmonic_coeff = spherical_harmonic_coeff.reshape(-1,bw**2)

        n_coeffs = spherical_harmonic_coeff.shape[0]
        mean_coeff = np.zeros((n_coeffs,bw**2),dtype=complex)
        for _id in range(n_coeffs):
            coeff = spherical_harmonic_coeff[_id]            
            for l in range(bw):
                l_start = l*(l+1)-l
                l_stop = l*(l+1)+l+1
                Ml_mat = self.lnk[l].reshape(2*l+1,2*l+1).T
                #np.roll is needed to convert -l,...,l format of spherical coeff to 0,..l,-l,...-1 in order to match the convention in Ml_mat that comes from ffts.
                mean_coeff_part = Ml_mat@np.roll(coeff[l_start:l_stop],l+1)
                #mean_coeff_part = Ml_mat@coeff[l_start:l_stop]
                mean_coeff[_id,l_start:l_stop] = np.roll(mean_coeff_part,-(l+1))
            
        if len(shape)>1:
            mean_coeff = mean_coeff.reshape(shape)
        else:
            mean_coeff = mean_coeff[0]
        return mean_coeff
        
        
    @property
    def normalize(self):
        so3_coeff =  self._soft.inverse_cmplx(self).reshape(self.so3_grid.shape[:-1])
        new_coeff =  self._soft.forward_cmplx(so3_coeff.flatten())
        scale = new_coeff[0].real*np.sqrt(2)
        so3_coeff*=scale
        self[:] = self._soft.forward_cmplx(so3_coeff.flatten())
        
class CharFuncFactory:
    """
    Factory that can create various characteristic functions of probability distribitions on SO(3)
    All of the distributions here use orthonormalized versions of the Wigner D matrices
    ...

    Methods
    -------
    uniform(bw):
        creates char function of uniform distribution
    """
    #General
    @staticmethod
    def uniform(bw = 16):
        r'''
        Creates uniform characteristif function
        $$ M_{lnk} = \delta_{l,0}\delta_{n,0}\delta_{k,0} $$
        '''
        u = CharFuncSO3(bw)
        u[0]= 1#/np.sqrt(2) #1/(np.sqrt(2)*2*np.pi)
        u[1:]=0
        return u

    @staticmethod
    def delta(bw = 16,angles=None):
        r'''
        Creates characteristic function for the deltadistribution 
        $$ M_{lnk} = D^l_{nk}(\alpha,\beta,\gamma) $$        
        '''        
        delta = CharFuncSO3(bw)
        if angles is not None:
            a,b,g = angles
            delta[:]=wigner.full_big_d_risbo(bw,a,b,g,True).conj()
        else:
            lnks = delta.lnk._lnks
            delta[:]=0
            for _id,(l,n,k) in enumerate(lnks):
                if n==k:
                    delta[_id] = np.sqrt((2*l+1)/(8*np.pi**2))
        return delta
    

    #1D distribution to SO(3) distribution via Bernstein functions
    @staticmethod
    def gaussian(bw = 16, sigma=1,center_rotation=None):
        if sigma<np.sqrt(2*np.pi)/bw:
            log.warning('sigma < sqrt(2*pi)/bw : distribution will contain negative values due to finite sampling. Incriase bw or decrease sigma to avoid this situation.')
        gauss = CharFuncSO3(bw)
        coeffs = gauss.lnk._lnks
        if center_rotation is None:
            for _id,(l,n,k) in enumerate(coeffs):
                if n == k:
                    gauss[_id] = np.sqrt((2*l+1)/(8*np.pi**2)) * np.exp(-sigma**2*l*(l+1)/2)
                else:
                    gauss[_id]=0
        else:
            a,b,g = center_rotation
            gauss[:]=wigner.full_big_d_risbo(bw,a,b,g,True)
            for _id,(l,n,k) in enumerate(coeffs):
                gauss[_id]*=np.exp(-sigma**2*l*(l+1)/2)
        return gauss

    @staticmethod
    def cauchy(bw = 16, gamma=1,center_rotation=None):
        cauchy = CharFuncSO3(bw)
        coeffs = cauchy._soft.lnks
        if center_rotation is None:
            for _id,(l,n,k) in enumerate(coeffs):
                if n == k:
                    cauchy[_id]=np.exp(-gamma*np.sqrt(l*(l+1)))
                else:
                    cauchy[_id]=0
        else:
            D = cauchy._soft.wigners
            cauchy[:]=D.abg[center_rotation[0],center_rotation[1],center_rotation[2]]
            for _id,(l,n,k) in enumerate(coeffs):
                cauchy[_id]*=np.exp(-gamma*np.sqrt(l*(l+1)))
        return cauchy

    @staticmethod
    def laplace(bw = 16, b=1,center_rotation=None):
        laplace = CharFuncSO3(bw)
        coeffs = laplace._soft.lnks
        if center_rotation is None:
            for _id,(l,n,k) in enumerate(coeffs):
                if n == k:
                    laplace[_id]=1/(1+b**2*l*(l+1))
                else:
                    laplace[_id]=0
        else:
            D = laplace._soft.wigners
            laplace[:]=D.abg[center_rotation[0],center_rotation[1],center_rotation[2]]
            for _id,(l,n,k) in enumerate(coeffs):
                laplace[_id]/=(1+b**2*l*(l+1))
        return laplace

    @staticmethod
    def voigt(bw = 16, sigma=1,gamma=1,center_rotation=None):
        voigt = CharFuncSO3(bw)
        coeffs = voigt._soft.lnks
        if center_rotation is None:
            for _id,(l,n,k) in enumerate(coeffs):
                if n == k:
                    voigt[_id]=np.exp(-gamma*np.sqrt(l*(l+1))-sigma**2*l*(l+1)/2)
                else:
                    voigt[_id]=0
        else:
            D = voigt._soft.wigners
            voigt[:]=D.abg[center_rotation[0],center_rotation[1],center_rotation[2]]
            for _id,(l,n,k) in enumerate(coeffs):
                voigt[_id]*=np.exp(-gamma*np.sqrt(l*(l+1))-sigma**2*l*(l+1)/2)
        return voigt    
            

    #Spherical Distributions
    @staticmethod
    def cos2n(bw=16,n=1):
        nn=n
        cos2n = CharFuncSO3(bw)
        Kcos = (2*nn+1)
        coeffs = cos2n._soft.lnks
        for _id,(l,n,k) in enumerate(coeffs):
            nonzero_value =  (l%2==0) and (n==0) and (k==0) and (l<=2*nn)
            if nonzero_value:
                #The following two equations are equivalent
                #cos2n[_id]=Kcos*np.sqrt(np.pi)*np.exp(loggamma(2*nn+1)-(2*nn+1)*np.log(2)-loggamma(nn-l//2+1)-loggamma(nn+l//2+3/2))
                cos2n[_id]=Kcos*np.exp(loggamma(2*nn+1)+loggamma(nn+l//2+1)+(l+1)*np.log(2)-loggamma(nn-l//2+1)-loggamma(2*nn+l+2))
            else:
                cos2n[_id]=0
        return cos2n

    @staticmethod
    def mises_fisher(bw=16,kappa=1):
        mf = CharFuncSO3(bw)
        Kmf = kappa/(2*np.sinh(kappa))
        coeffs = mf._soft.lnks
        for _id,(l,n,k) in enumerate(coeffs):
            nonzero_value = (n==0) and (k==0)
            if nonzero_value:
                mf[_id]= Kmf*np.sqrt(2*np.pi/kappa)*gsl.bessel_Inu(l+1/2,kappa)
            else:
                mf[_id]=0
        return mf
    @staticmethod
    def watson(bw=16,kappa=1):
        watson = CharFuncSO3(bw)
        so3_grid = watson.so3_grid
        Kw=2/(np.sqrt(np.pi/kappa)*error_function_imag(np.sqrt(kappa))) #/(1.j*np.sqrt(kappa)*error_function(1.j*np.sqrt(kappa)))
        distribution = Kw*np.exp(kappa*np.cos(so3_grid[...,1])**2).astype(complex)
        watson[:] = watson._soft.forward_cmplx(distribution/watson.so3_const)
        return watson
