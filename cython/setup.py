
from distutils.core import setup
from Cython.Build import cythonize

setup(
    ext_modules = cythonize("seq_to_vector.pyx")
)
