'''
Wrapper for fftw options
'''
from pysofft import _soft
from multiprocessing import cpu_count

def set_num_threads(nthreads:int=None):
    r'''
    Set the maximal number of OpenMP threads.
    
    Parameters
    ----------
    nthreads:int64
        Maximal number of threads. By default nthreads is set to half the available cpu threads.
    '''
    if nthreads is None:
        nthreads = cpu_count()//2
    
    _soft.py.omp_set_num_threads_(nthreads)
    
def get_max_threads()->int:
    r'''
    Get the current maximal number of OpenMP threads.
    
    Returns
    -------
    nthreads:int64
        Current maximal number of threads.
    '''
    return _soft.py.omp_get_max_threads_()
#set_num_threads = _soft.py.omp_set_num_threads_
#get_max_threads = _soft.py.omp_get_max_threads_
