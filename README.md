# BATSim

BATSim is designed to enable the application of non-affine shear
transformations to GalSim images. It works by trasnforming the stamp pixel
grid, distorting the location of the pixel centers. When a GalSim image is
sampled onto this distorted grid, the resulting image will be sheared due to
the transformation of the pixel grid.

Non-affine transforms can simualate complex shear effects such as intrinsic
alignment, flexion, and optical field distrotion maps. Custom transform
functions can also be passed to Stamp objects.

## Building and installing BATSim from source (Conda + GalSim C++)

The package is not currently available on PyPI and only installable on Linux while in early development. We suggest using conda/mamba to build the package and install its dependencies. This is because the primary dependency, Galsim, does not ship a pre-built C++ library via pip.

The below instructions will work with pure conda, but we recommend using mamba to make the installation of dependencies much quicker.

First, make sure to clone and switch to the repository root:

```bash
git clone https://github.com/CMacM/BATSim.git
cd BATSim
```

## Regular Install (not editable)

### 1. Create build environment and activate
```bash
mamba create -n batsim -c conda-forge -c defaults python=3.10 conda-build boa
mamba activate batsim
```

### 2. Build the package
From the repository root:
```bash
# mambabuild is sometimes only recognised as a conda command
conda mambabuild --override-channels -c conda-forge -c defaults conda/recipe
```

This will:
- Create an isolated environment to build and run batsim
- Install all required dependencies
- Compile the BATSim C++ extension and link it to Galsim's C++ API
- Produce a conda package in:
```bash
$CONDA_PREFIX/conda-bld/linux-64/
```

### 3. Install the built package
Install the locally built package into the new environment:
```bash
mamba install -c local batsim
```

If this fails, install directly from the build artifact:
```bash
mamba install $CONDA_PREFIX/conda-bld/linux-64/batsim-*.tar.bz2
```

## Development Installation (Editable)

For development you can install BATSim in editable mode so that Python changes
take effect immediately without reinstalling.

### 1. Create a development environment and activate

Note: You may encounter issues if some of these packages are already installed locally and have different builds. To fix, remove them and ensure they are installed through conda-forge.

```bash
mamba create -n batsim-dev -c conda-forge -c defaults \
    python=3.10 \
    galsim \
    eigen \
    pybind11 \
    numpy \
    fitsio \
    astropy \
    matplotlib
mamba activate batsim-dev
```

### 2. Install BATSim in editable mode

```bash
pip install -e .
```

After following either of the above installation routes, with the conda environment active, verify the installation:

```bash
python -c "import batsim; import batsim._gsinterface; print('BATSim installed successfully')"
```

## Pip only installation

We currently do not provide a pip only installation route as Galsim requires the shared C++ library to be built and linked manually when installed via pip. If you wish to build via pip only, you will need to install all dependencies and build and link the Galsim C++ headers. You can find details on how to do this [here][https://galsim-developers.github.io/GalSim/_build/html/install_pip.html]. You may then need to update your LD\_LIBRARY\_PATH, LIBRARY\_PATH, and CPLUS\_INCLUDE\_PATH to point to build and include folders for the GalSim C++ shared library.

We plan to add a pure pip installation route in future.

![BATSim Logo](./image/batsim_logo.png)

(repo icon created using images by Frepik and brgfx on Freepik)
