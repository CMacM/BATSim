import os
import pybind11
from setuptools import setup, find_packages, Extension

# Get the path to the GalSim directory
galsim_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), 
                                          'extern/GalSim'))

# Now get include dirs
include_dirs = [
    pybind11.get_include(),
    os.path.join(galsim_dir, 'include'),  # Use the include directory under GALSIM_DIR
    os.path.join(galsim_dir, 'include/galsim'),
    os.path.join(galsim_dir, 'downloaded_eigen/eigen-3.4.0')
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
    extra_compile_args=['-std=c++11', '-fopenmp', '-O3'],
    extra_link_args=['-flto','-fopenmp','-Wl,-rpath,'+galsim_lib_path]
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