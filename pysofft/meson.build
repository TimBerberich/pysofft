add_languages('fortran','c')

omp = dependency('openmp')
fftw = dependency('fftw3')
fftw_include = fftw.get_variable('includedir')
fs = import('fs')
py_mod = import('python')

py = py_mod.find_installation(required: true, pure: false)
py_dep = py.dependency()



incdir_numpy = run_command(py,
  ['-c', 'import os; import numpy; print(os.path.relpath(numpy.get_include()))'],
  check : true
).stdout().strip()

incdir_f2py = run_command(py,
    ['-c', 'import os; import numpy.f2py; print(os.path.relpath(numpy.f2py.get_include()))'],
    check : true
).stdout().strip()

inc_np = include_directories(incdir_numpy, incdir_f2py,fftw_include)

f2py_skip = fs.read('f2py_skip.txt',encoding: 'utf-8')

# Order matters for f2py scanning: provider module first
fortran_sources = files(
  'precision.f90',
  'math_constants.f90',
  'utils.f90',
  'wigner.f90',
  'soft.f08',
  'py_wrap.f90'	
)

pysofft_source = custom_target('_softmodule.c',
  input : fortran_sources,  # .f so no F90 wrappers
  output : ['_softmodule.c', '_soft-f2pywrappers2.f90'],
  command : [py, '-m', 'numpy.f2py','-I','/usr/include','-m', '_soft',
             '--lower','skip:',f2py_skip,':','@INPUT@', '--build-dir', 'pysofft'],
)


run_command(py,
  ['-c', 'import os; import numpy; print("f2py is done")'],
  check : true
).stdout().strip()

py.extension_module('_soft',
  [
    fortran_sources,
	pysofft_source
  ],
  incdir_f2py / 'fortranobject.c',
  # for debugging: -g -fbacktrace in fortran_args and -g in c_args
  # fortran_args: ['-fopenmp','-lfftw3', '-O3', '-fno-math-errno', '-fno-trapping-math','-g','-fbacktrace'],
  fortran_args: ['-fopenmp','-lfftw3', '-O3', '-fno-math-errno', '-fno-trapping-math'],	
  # c_args: ['-fopenmp','-lfftw3', '-O3', '-fno-math-errno', '-fno-trapping-math','-g'],
  c_args: ['-fopenmp','-lfftw3', '-O3', '-fno-math-errno', '-fno-trapping-math'],	
  link_args : ['-fopenmp','-lfftw3','-lm'],
  include_directories: inc_np,
  dependencies : [py_dep, omp, fftw],
  subdir: 'pysofft',
  install : true,
)
#add_global_arguments('--opt="-O3 -fno-math-errno -fno-trapping-math"', language : 'fortran')

python_sources = [
  'fftw.py',
   'omp.py',		
  'soft.py',
  'stats.py',
  '__init__.py'
]

py.install_sources(
  python_sources,
  pure: false,
  subdir: 'pysofft'
)

