from distutils.core import setup
from distutils.extension import Extension
import numpy

USE_CYTHON = False

ext = '.pyx' if USE_CYTHON else '.cpp'

extensions = [Extension("algorithms", ["algorithms"+ext], language="c++", extra_compile_args = ["-std=c++11"], include_dirs=[numpy.get_include()])]

if USE_CYTHON:
    from Cython.Build import cythonize
    extensions = cythonize(extensions)

setup( ext_modules = extensions )

