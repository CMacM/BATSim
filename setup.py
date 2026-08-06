import os
import sys
import sysconfig

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

version_ns = {}
with open(os.path.join("src", "batsim", "_version.py")) as version_file:
    exec(version_file.read(), version_ns)


class BuildExt(build_ext):
    def build_extensions(self):
        pybind11_include_dirs = self._find_pybind11_paths()
        include_dirs, library_dirs = self._find_galsim_paths()

        if not self._has_library(library_dirs, "galsim"):
            search_paths = ", ".join(library_dirs) or "<none>"
            raise RuntimeError(
                "Could not find GalSim's C++ shared library libgalsim. "
                "Install GalSim from conda-forge, or build GalSim's shared C++ "
                "library from the GalSim source tree with "
                "`python setup.py build_shared_clib` and set GALSIM_LIB_DIR "
                "and GALSIM_INCLUDE_DIR before installing BATSim. "
                f"Searched library paths: {search_paths}"
            )

        for ext in self.extensions:
            ext.include_dirs.extend(
                [
                    *pybind11_include_dirs,
                    *include_dirs,
                ]
            )
            ext.library_dirs.extend(library_dirs)
            ext.runtime_library_dirs.extend(library_dirs)
            ext.libraries.append("galsim")

        super().build_extensions()

    def _find_pybind11_paths(self):
        import pybind11

        candidates = [
            pybind11.get_include(),
            os.path.join(sys.prefix, "include"),
        ]
        include_dirs = []
        for directory in self._dedupe(candidates):
            if self._has_pybind11_conduit(directory):
                include_dirs.append(directory)

        if include_dirs:
            return include_dirs

        search_paths = ", ".join(self._dedupe(candidates)) or "<none>"
        raise RuntimeError(
            "BATSim requires pybind11 >= 3 with the conduit headers. "
            f"Searched pybind11 include paths: {search_paths}"
        )

    def _has_pybind11_conduit(self, include_dir):
        return os.path.exists(os.path.join(include_dir, "pybind11", "conduit", "pybind11_conduit_v1.h"))

    def _find_galsim_paths(self):
        include_dirs = []
        library_dirs = []

        include_override = os.environ.get("GALSIM_INCLUDE_DIR")
        if include_override:
            include_dirs.extend([include_override, os.path.join(include_override, "galsim")])

        lib_override = os.environ.get("GALSIM_LIB_DIR")
        if lib_override:
            library_dirs.append(lib_override)

        prefixes = [
            sys.prefix,
            os.environ.get("CONDA_PREFIX"),
            os.environ.get("PREFIX"),
            sysconfig.get_config_var("prefix"),
        ]

        for prefix in prefixes:
            if not prefix:
                continue
            include_dirs.extend(
                [
                    os.path.join(prefix, "include"),
                    os.path.join(prefix, "include", "galsim"),
                    os.path.join(prefix, "include", "eigen3"),
                ]
            )
            lib_dir = os.path.join(prefix, "lib")
            if os.path.isdir(lib_dir):
                library_dirs.append(lib_dir)

        if not self._has_header(include_dirs, "GalSim.h"):
            try:
                import galsim

                galsim_include = getattr(galsim, "include_dir", None)
                if galsim_include:
                    include_dirs.extend([galsim_include, os.path.join(galsim_include, "galsim")])
            except ImportError:
                pass

        return self._dedupe(include_dirs), self._dedupe(library_dirs)

    def _has_header(self, include_dirs, name):
        return any(os.path.exists(os.path.join(directory, name)) for directory in include_dirs)

    def _has_library(self, library_dirs, name):
        patterns = [f"lib{name}.so", f"lib{name}.dylib", f"{name}.lib"]
        return any(
            os.path.exists(os.path.join(directory, pattern))
            for directory in library_dirs
            for pattern in patterns
        )

    def _dedupe(self, values):
        out = []
        for value in values:
            if value and value not in out:
                out.append(value)
        return out


gsinterface = Extension(
    "batsim._gsinterface",
    sources=["src/batsim/_gsinterface.cpp"],
    include_dirs=[],
    libraries=[],
    language="c++",
    extra_compile_args=["-std=c++17", "-fopenmp", "-O3"],
    extra_link_args=["-fopenmp"],
)

setup(
    name="batsim",
    version=version_ns["__version__"],
    author="Charlie MacMahon, Andy Park",
    author_email="c.macmahon@ncl.ac.uk, chanhyup@andrew.cmu.edu",
    license="MIT",
    python_requires=">=3.10,<3.13",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: Unix",
        "Programming Language :: Python",
    ],
    packages=find_packages(where="src"),
    py_modules=["batsim_install_gpu"],
    package_dir={"": "src"},
    ext_modules=[gsinterface],
    cmdclass={"build_ext": BuildExt},
    install_requires=[
        "numpy>=1.26,<2.0",
        "galsim",
        "fitsio",
        "matplotlib>=3.8,<3.9",
        "astropy>=6.0,<6.1",
    ],
    extras_require={
        "benchmark": ["asv>=0.6"],
        "cuda12": ["cupy-cuda12x>=13.6,<14"],
        "cuda13": ["cupy-cuda13x>=13.6,<14"],
    },
    entry_points={
        "console_scripts": [
            "batsim-install-gpu=batsim_install_gpu:main",
        ],
    },
)
