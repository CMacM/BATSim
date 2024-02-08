import os
from setuptools import setup, Extension
import pybind11

galsim_dir = os.getenv('GALSIM_DIR', '/default/path/if/any')
include_dirs = [
    pybind11.get_include(),
    os.path.join(galsim_dir, 'include'),  # Use the include directory under GALSIM_DIR
    os.path.join(galsim_dir, 'downloaded_eigen/eigen-3.4.0/Eigen')
]

gsinterface_module = Extension(
    '_gsinterface',
    sources=['src/_gsinterface.cpp'],
    include_dirs=include_dirs,
    language='c++',
    extra_compile_args=['-std=c++11'],
)

setup(
    name='YourPackageName',
    version='0.1',
    description='Your package description',
    ext_modules=[gsinterface_module],
    install_requires=['pybind11>=2.2'],
)
