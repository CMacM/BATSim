import os, glob, sys
import pybind11
import distutils.sysconfig
from setuptools import setup, find_packages, Extension

# Get the path to the GalSim directory
galsim_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                          'extern/GalSim'))

# Find Eigen, below code is copied from GalSims setup.py
try_dirs = []

# Start with a user-specified directory.
if 'EIGEN_DIR' in os.environ:
    try_dirs.append(os.environ['EIGEN_DIR'])
    try_dirs.append(os.path.join(os.environ['EIGEN_DIR'], 'include'))

# Add the python system include directory.
try_dirs.append(distutils.sysconfig.get_config_var('INCLUDEDIR'))

# If using Anaconda, add their lib dir in case fftw is installed there.
# (With envs, this might be different than the sysconfig LIBDIR.)
if 'CONDA_PREFIX' in os.environ:
    try_dirs.append(os.path.join(os.environ['CONDA_PREFIX'],'lib'))

# Some standard install locations:
try_dirs.extend(['/usr/local/include', '/usr/include'])
if sys.platform == "darwin":
    try_dirs.extend(['/sw/include', '/opt/local/include'])

# Also if there is a C_INCLUDE_PATH, check those dirs.
for path in ['C_INCLUDE_PATH']:
    if path in os.environ:
        for dir in os.environ[path].split(':'):
            try_dirs.append(dir)

# Finally, (last resort) check our own download of eigen.
if os.path.isdir('extern/GalSim/downloaded_eigen'):
    try_dirs.extend(glob.glob(os.path.join('extern/GalSim/downloaded_eigen','*')))

for dir in try_dirs:
    if dir is None: continue
    if not os.path.isdir(dir): continue
    if os.path.isfile(os.path.join(dir, 'Eigen/Core')):
        eigen_dir = dir
    if os.path.isfile(os.path.join(dir, 'eigen3', 'Eigen/Core')):
        dir = os.path.join(dir, 'eigen3')
        # Only print this if the eigen3 addition was key to finding it.
        print("\n  ", dir, "  (yes)")
        eigen_dir = dir

# Now get include dirs
include_dirs = [
    pybind11.get_include(),
    os.path.join(galsim_dir, 'include'),  # Use the include directory under GALSIM_DIR
    os.path.join(galsim_dir, 'include/galsim'),
    eigen_dir
]

galsim_lib_path = os.path.join(galsim_dir, 'build/shared_clib')

gsinterface = Extension(
    'batsim._gsinterface',
    sources=['src/batsim/_gsinterface.cpp'],
    include_dirs=include_dirs,
    library_dirs=[galsim_lib_path],
    runtime_library_dirs=[galsim_lib_path],
    libraries=['galsim', 'galsim.2.5'],
    language='c++',
    extra_compile_args=['-std=c++11', '-fopenacc', '-fopenmp', '-O3'],
    extra_link_args=['-flto', '-fopenacc', '-fopenmp','-Wl,-rpath,'+galsim_lib_path]
)

scripts = ['src/batsim']
scripts = [os.path.join('bin',f) for f in scripts ]

setup(
    name='batsim',
    author='Charlie MacMahon, Andy Park',
    author_email='c.macmahon@ncl.ac.uk, chanhyup@andrew.cmu.edu',
    license='MIT',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Operating System :: Unix',
        'Programming Language :: Python',
    ],
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'ipykernel',
        'numpy',
        'matplotlib',
        'scipy',
        'pybind11>=2.2'
    ],
    extras_require={
        'dev': [
            'pytest',
            'pytest-cov',
            'pre-commit',
            'sphinx',
            'sphinx-rtd-theme',
            'sphinx-autoapi',
            'pylint',
            'nbconvert',
            'nbsphinx',
            'ipython',
            'asv==0.5.1',
            # "fpfs"
        ],
    },
    ext_modules=[gsinterface],
)