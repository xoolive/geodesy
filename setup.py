from setuptools import setup
from setuptools.extension import Extension
# from distutils.command.build_ext import build_ext
from Cython.Build import cythonize

INCLUDE = []

try:
    import numpy
    INCLUDE.append(numpy.get_include())  # __path__[0] + "/core/include")
except ImportError:
    print "numpy is required in order to use geodesy toolkit"
    raise

# copt = {'msvc': ['/arch:SSE2'],
#         }
# lopt = {'msvc': [''],
#         }

# class build_ext_subclass(build_ext):
#     def build_extensions(self):
#         c = self.compiler.compiler_type
#         if c in copt:
#             for e in self.extensions:
#                 e.extra_compile_args = copt[c]
#         if c in lopt:
#             for e in self.extensions:
#                 e.extra_link_args = lopt[c]
#         build_ext.build_extensions(self)

extensions = [
    Extension("geodesy.sphere",
              ["geodesy/sphere.pyx"],
              include_dirs=INCLUDE, language="c++"),
    Extension("geodesy.wgs84",
              ["geodesy/wgs84.pyx"],
              include_dirs=INCLUDE, language="c++"),
]

setup(name="geodesy",
      version="1.0",
      author="Xavier Olive",
      author_email="xavier.olive@onera.fr",
      description="Toolkit for geodesy",
      ext_modules=cythonize(extensions),
      packages=['geodesy', ],
      )
#       cmdclass={'build_ext': build_ext_subclass},
