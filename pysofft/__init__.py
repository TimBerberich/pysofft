from pysofft import _soft
from pysofft.soft import Soft
from pysofft import fftw
from pysofft import omp
from multiprocessing import cpu_count


# set default num threads
omp.set_num_threads(cpu_count()//2)

utils = _soft.utils
wigner = _soft.make_wigner
