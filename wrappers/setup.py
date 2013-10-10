# on windows:
# set VS90COMNTOOLS=%VS100COMNTOOLS%
# python setup.py build_ext --compiler=msvc

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import glob

sources = glob.glob('../source/*.cpp')
sources.append('terran.pyx')

ext_modules = [Extension('terran',
                        include_dirs = ['../include'],
                        language = 'c++',
                        sources = sources,
                        define_macros = [('TERRAN_BUILDING_SHARED_LIBRARY', '1')],
                        )]
                        
setup(
    name = 'TERRAN',
    cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules
)