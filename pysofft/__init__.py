from pysofft import _soft
from pysofft import fftw
from pysofft.soft import Soft
from pysofft import stats
from multiprocessing import cpu_count

OMP_set_num_threads = _soft.py.omp_set_num_threads_
OMP_get_max_threads = _soft.py.omp_get_max_threads_
# set default num threads
OMP_set_num_threads(cpu_count()//2)

utils = _soft.utils
wigner = _soft.make_wigner

