from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

ext_modules = [
    Extension('xcorspc',
              ['xcorspc.pyx'])
]

setup(
  name = 'XCOR_C app',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)
