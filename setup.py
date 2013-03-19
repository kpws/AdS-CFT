from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


setup(cmdclass = {'build_ext': build_ext},
      ext_modules = [Extension("solveBulk",
          sources=["solveBulk.pyx","solveODE.c","bulkODE.c"],
          include_dirs=[],
          libraries=['gsl','gslcblas','m'],
          language="c")])
