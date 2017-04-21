from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import os, sys, subprocess, numpy

# get paths to CRTM stuff.
CRTM_incdir = os.environ.get('CRTM_INCDIR')
if CRTM_incdir is None:
    CRTM_incdir='CRTM_REL-2.2.3/libsrc'
CRTM_libdir = os.environ.get('CRTM_LIBDIR')
if CRTM_libdir is None:
    CRTM_libdir='CRTM_REL-2.2.3/libsrc'

# build iso_c_binding fortran wrapper.
compiler = 'gfortran'
fname = 'pycrtm_interface'
fopts = '-c -O2 -fPIC -I %s' % CRTM_incdir

fsource = '%s.f90' % fname
strg = '%s %s %s' % (compiler, fopts, fsource)
sys.stdout.write('executing "%s"\n' % strg)
subprocess.call(strg,shell=True)

if os.path.exists('pycrtm.c'): os.remove('pycrtm.c') # trigger a rebuild

# build c extension that calls fortran.
libs = ['crtm','gfortran']
inc_dirs = [CRTM_incdir]
inc_dirs.append(numpy.get_include())
lib_dirs = [CRTM_libdir,'/opt/local/lib/gcc6']
objs = ['%s.o' % fname] # fortran object to link
ext_modules = [Extension('pycrtm',                       # module name
                        ['pycrtm.pyx'],                  # cython source file
                        libraries     = libs,
                        extra_objects = objs,
                        include_dirs  = inc_dirs,
                        library_dirs  = lib_dirs)]

setup(name = 'pycrtm',
      cmdclass = {'build_ext': build_ext},
      ext_modules  = ext_modules)
