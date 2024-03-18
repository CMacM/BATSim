import os
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


class CustomBuildExt(build_ext):
    def run(self):
        # Dynamically set include_dirs and library_dirs before building
        # extensions
        include_dirs, lib_dirs = self.find_galsim_paths()

        for ext in self.extensions:
            for _ in include_dirs:
                ext.include_dirs.append(_)
            for _ in lib_dirs:
                ext.library_dirs.append(_)
                ext.runtime_library_dirs.append(_)  # For Linux

        super().run()

    def find_galsim_paths(self):
        # Implement logic to locate GalSim's include and library directories
        # This is a placeholder implementation; adjust based on your setup

        # Example: Assuming GalSim is installed in a conda environment
        conda_prefix = os.environ.get("CONDA_PREFIX")
        include_dirs = []
        lib_dirs = []
        if conda_prefix:
            include_dirs.append(os.path.join(conda_prefix, "include"))
            include_dirs.append(os.path.join(conda_prefix, "include/galsim"))
            include_dirs.append(os.path.join(conda_prefix, "include/Eigen"))
            lib_dirs.append(os.path.join(conda_prefix, "lib"))
        else:
            # Fallback or other logic to locate GalSim
            include_dirs.append("/usr/local/include")
            include_dirs.append("/usr/local/include/galsim")
            include_dirs.append("/usr/local/include/Eigen")
            lib_dirs.append("/usr/local/lib")

        return include_dirs, lib_dirs


# Define your extension module
gsinterface = Extension(
    "batsim._gsinterface",
    sources=["src/batsim/_gsinterface.cpp"],
    libraries=["galsim"],
    language="c++",
    extra_compile_args=["-std=c++11", "-fopenmp", "-O3"],
    extra_link_args=["-flto", "-fopenmp"],
)

setup(
    name="batsim",
    author="Charlie MacMahon, Andy Park",
    author_email="c.macmahon@ncl.ac.uk, chanhyup@andrew.cmu.edu",
    license="MIT",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Programming Language :: Python",
    ],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=["numpy", "pybind11>=2.2"],
    ext_modules=[gsinterface],
    cmdclass={"build_ext": CustomBuildExt},
)
