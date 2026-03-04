import os
from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext
import pybind11


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
        include_dirs = []
        lib_dirs = []

        # Prefer GalSim's own include directory if import works
        try:
            import galsim
            inc = galsim.include_dir              # .../site-packages/galsim/include
            include_dirs.append(inc)
            include_dirs.append(os.path.join(inc, "galsim"))   # <-- add this
        except Exception:
            print("Error: Could not import GalSim to find include directory.")

        # Conda-build uses PREFIX for the host env (headers live here)
        prefixes = [
            os.environ.get("PREFIX"),        # conda-build host env
            os.environ.get("CONDA_PREFIX"),  # active env fallback
        ]

        for p in prefixes:
            if not p:
                continue
            include_dirs.append(os.path.join(p, "include"))
            include_dirs.append(os.path.join(p, "include", "galsim"))
            include_dirs.append(os.path.join(p, "include", "eigen3"))
            lib_dirs.append(os.path.join(p, "lib"))

        # 3) Last-resort system paths
        include_dirs += ["/usr/local/include", "/usr/include", "/usr/include/eigen3"]
        lib_dirs += ["/usr/local/lib", "/usr/lib"]

        # De-duplicate preserving order
        def uniq(xs):
            out = []
            for x in xs:
                if x and x not in out:
                    out.append(x)
            return out

        return uniq(include_dirs), uniq(lib_dirs)


# Define your extension module
gsinterface = Extension(
    "batsim._gsinterface",
    sources=["src/batsim/_gsinterface.cpp"],
    include_dirs=[pybind11.get_include()],
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
    install_requires=[
        "numpy", "pybind11>=2.2", "fitsio", "matplotlib", "astropy",
    ],
    ext_modules=[gsinterface],
    cmdclass={"build_ext": CustomBuildExt},
)
