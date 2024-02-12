import os, sys, glob, re
from setuptools import setup, Extension
import pybind11
import ctypes
import ctypes.util
import urllib.request as urllib2
import tarfile
import shutil
from setuptools import setup, find_packages, Extension
from setuptools.command.build_ext import build_ext

# Turn this on for more verbose debugging output about compile attempts.
debug = False

# Check for the fftw3 library in some likely places
def find_fftw_lib(output=False):
    import distutils.sysconfig

    if debug: output = True
    try_libdirs = []

    # Start with the explicit FFTW_DIR, if present.
    if 'FFTW_DIR' in os.environ:
        try_libdirs.append(os.environ['FFTW_DIR'])
        try_libdirs.append(os.path.join(os.environ['FFTW_DIR'],'lib'))

    # Add the python system library directory.
    try_libdirs.append(distutils.sysconfig.get_config_var('LIBDIR'))

    # If using Anaconda, add their lib dir in case fftw is installed there.
    # (With envs, this might be different than the sysconfig LIBDIR.)
    if 'CONDA_PREFIX' in os.environ:
        try_libdirs.append(os.path.join(os.environ['CONDA_PREFIX'],'lib'))

    # Try some standard locations where things get installed
    try_libdirs.extend(['/usr/local/lib', '/usr/lib'])
    if sys.platform == "darwin":
        try_libdirs.extend(['/sw/lib', '/opt/local/lib'])

    # Check the directories in LD_LIBRARY_PATH.  This doesn't work on OSX >= 10.11
    for path in ['LIBRARY_PATH', 'LD_LIBRARY_PATH', 'DYLD_LIBRARY_PATH']:
        if path in os.environ:
            for dir in os.environ[path].split(':'):
                try_libdirs.append(dir)

    # The user's home directory is often a good place to check.
    try_libdirs.append(os.path.join(os.path.expanduser("~"),"lib"))

    # If the above don't work, the fftw3 module may have the right directory.
    try:
        import fftw3
        try_libdirs.append(fftw3.lib.libdir)
    except ImportError:
        pass

    if sys.platform == "darwin":
        lib_ext = '.dylib'
    else:
        lib_ext = '.so'
    name = 'libfftw3' + lib_ext
    if output: print("Looking for ",name)
    tried_dirs = set()  # Keep track, so we don't try the same thing twice.
    for dir in try_libdirs:
        if dir == '': continue  # This messes things up if it's in there.
        if dir in tried_dirs: continue
        else: tried_dirs.add(dir)
        if not os.path.isdir(dir): continue
        libpath = os.path.join(dir, name)
        if not os.path.isfile(libpath): continue
        if output: print("  ", dir, end='')
        try:
            lib = ctypes.cdll.LoadLibrary(libpath)
            if output: print("  (yes)")
            return libpath
        except OSError as e:
            if output: print("  (no)")
            # Some places use lib64 rather than/in addition to lib.  Try that as well.
            if dir.endswith('lib') and os.path.isdir(dir + '64'):
                dir += '64'
                try:
                    libpath = os.path.join(dir, name)
                    if not os.path.isfile(libpath): continue
                    lib = ctypes.cdll.LoadLibrary(libpath)
                    if output: print("  ", dir, "  (yes)")
                    return libpath
                except OSError:
                    pass

    # If we didn't find it anywhere, but the user has set FFTW_DIR, trust it.
    if 'FFTW_DIR' in os.environ:
        libpath = os.path.join(os.environ['FFTW_DIR'], name)
        print("WARNING:")
        print("Could not find an installed fftw3 library named %s"%(name))
        print("Trusting the provided FFTW_DIR=%s for the library location."%(libpath))
        print("If this is incorrect, you may have errors later when linking.")
        return libpath

    # Last ditch attempt.  Use ctypes.util.find_library, which sometimes manages to find it
    # when the above attempts fail.
    try:
        libpath = ctypes.util.find_library('fftw3')
        if libpath == None:
            raise OSError
        if os.path.split(libpath)[0] == '':
            # If the above doesn't return a real path, try this instead.
            libpath = ctypes.util._findLib_gcc('fftw3')
            if libpath == None:
                raise OSError
        libpath = os.path.realpath(libpath)
        lib = ctypes.cdll.LoadLibrary(libpath)
    except Exception as e:
        print("Could not find fftw3 library.  Make sure it is installed either in a standard ")
        print("location such as /usr/local/lib, or the installation directory is either in ")
        print("your LIBRARY_PATH or FFTW_DIR environment variable.")
        raise
    else:
        dir, name = os.path.split(libpath)
        if output:
            if dir == '': dir = '[none]'
            print("  ", dir, "  (yes)")
        return libpath


# Check for Eigen in some likely places
def find_eigen_dir(output=False):
    if debug: output = True
    import distutils.sysconfig

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
    if os.path.isdir('downloaded_eigen'):
        try_dirs.extend(glob.glob(os.path.join('downloaded_eigen','*')))

    if output: print("Looking for Eigen:")
    for dir in try_dirs:
        if dir is None: continue
        if not os.path.isdir(dir): continue
        if output: print("  ", dir, end='')
        if os.path.isfile(os.path.join(dir, 'Eigen/Core')):
            if output: print("  (yes)")
            return dir
        if os.path.isfile(os.path.join(dir, 'eigen3', 'Eigen/Core')):
            dir = os.path.join(dir, 'eigen3')
            if output:
                # Only print this if the eigen3 addition was key to finding it.
                print("\n  ", dir, "  (yes)")
            return dir
        if output: print("  (no)")

    if output:
        print("Could not find Eigen in any of the standard locations.")
        print("Will now try to download it from gitlab.com. This requires an internet")
        print("connection, so it will fail if you are currently offline.")
        print("If Eigen is installed in a non-standard location, and you want to use that")
        print("instead, you should make sure the right directory is either in your")
        print("C_INCLUDE_PATH or specified in an EIGEN_DIR environment variable.")

    try:
        dir = 'downloaded_eigen'
        if os.path.isdir(dir):
            # If this exists, it was tried above and failed.  Something must be wrong with it.
            print("Previous attempt to download eigen found. Deleting and trying again.")
            shutil.rmtree(dir)
        os.mkdir(dir)
        url = 'https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.bz2'
        if output:
            print("Downloading eigen from ",url)
        # Unfortunately, gitlab doesn't allow direct downloads. We need to spoof the request
        # so it thinks we're a web browser.
        # cf. https://stackoverflow.com/questions/42863240/how-to-get-round-the-http-error-403-forbidden-with-urllib-request-using-python
        page=urllib2.Request(url,headers={'User-Agent': 'Mozilla/5.0'})
        data=urllib2.urlopen(page).read()
        fname = 'eigen.tar.bz2'
        with open(fname, 'wb') as f:
            f.write(data)
        if output:
            print("Downloaded %s.  Unpacking tarball."%fname)
        with tarfile.open(fname) as tar:
            
            def is_within_directory(directory, target):
                
                abs_directory = os.path.abspath(directory)
                abs_target = os.path.abspath(target)
            
                prefix = os.path.commonprefix([abs_directory, abs_target])
                
                return prefix == abs_directory
            
            def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
                # Avoid security vulnerability in tar.extractall function.
                # This bit of code was added by the Advanced Research Center at Trellix in PR #1188.
                # For more information about the security vulnerability, see
                # https://github.com/advisories/GHSA-gw9q-c7gh-j9vm
            
                for member in tar.getmembers():
                    member_path = os.path.join(path, member.name)
                    if not is_within_directory(path, member_path):
                        raise Exception("Attempted Path Traversal in Tar File")
            
                tar.extractall(path, members, numeric_owner=numeric_owner) 
            
            safe_extract(tar, dir)
        os.remove(fname)
        # This actually extracts into a subdirectory with a name eigen-eigen-5a0156e40feb/
        # I'm not sure if that name is reliable, so use glob to get it.
        dir = glob.glob(os.path.join(dir,'*'))[0]
        if os.path.isfile(os.path.join(dir, 'Eigen/Core')):
            return dir
        elif output:
            print("Downloaded eigen, but it didn't have the expected Eigen/Core file.")
    except Exception as e:
        if output:
            print("Error encountered while downloading Eigen from the internet")
            print(e)

    raise OSError("Could not find Eigen")

galsim_dir = os.getenv('GALSIM_DIR', '')
eigen_dir = find_eigen_dir()
include_dirs = [
    pybind11.get_include(),
    os.path.join(galsim_dir, 'include'),  # Use the include directory under GALSIM_DIR
    os.path.join(galsim_dir, 'include/galsim'),
    os.path.join(eigen_dir)
]

gsinterface = Extension(
    '_gsinterface',
    sources=['src/_gsinterface.cpp'],
    include_dirs=include_dirs,
    libraries=[os.getenv('GALSIM_DIR')],
    language='c++',
    extra_compile_args=['-std=c++11'],
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
        'galsim',
        'numpy',
        'matplotlib',
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