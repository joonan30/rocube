from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy 

ext_modules = [Extension("rocube", ["rocube.pyx"])]

setup(
  name = 'rocube',
  cmdclass = {'build_ext': build_ext},
  ext_modules = cythonize("rocube.pyx", include_path = [numpy.get_include()])
)

# python setup.py build_ext --inplace