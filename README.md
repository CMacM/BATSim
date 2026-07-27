# BATSim

BATSim renders GalSim surface-brightness profiles on coordinate grids that can
be transformed before sampling. This makes it possible to apply non-affine shear
fields, such as intrinsic alignment and flexion, while keeping the rendering
pipeline close to GalSim's image conventions.

The current renderer samples profiles on a supersampled grid, optionally applies
sub-pixel block integration, and performs PSF and pixel convolution in Fourier
space. The default render path uses:

- `psf_mode="kvalue"` for analytic PSF Fourier sampling through BATSim's C++ layer
- `force_input_flux=True` to preserve the input galaxy flux after convolution
- `use_true_center=True` to match GalSim's true-image-center convention
- `integration_order=2` for Gauss-Legendre block integration
- `compensate_integration="quadrature"` to remove the matching Gauss-Legendre
  transfer function; pass `"exact_sinc"` for ideal top-hat compensation or
  `None` to disable compensation
- `backend="np"` by default, with optional CuPy support for FFT-heavy steps

## Quickstart

```python
import galsim
import batsim

galaxy = galsim.Sersic(n=1.0, half_light_radius=0.7, flux=1.0)
psf = galsim.Gaussian(fwhm=0.7)
transform = batsim.IaTransform(scale=0.2, hlr=0.7, A=0.02)

image = batsim.simulate_galaxy(
    galaxy,
    scale=0.2,
    ngrid=64,
    transform_obj=transform,
    psf_obj=psf,
)
```

`simulate_galaxy` returns a NumPy array. NumPy is used by default; pass
`backend="cp"` to request CuPy for FFT-heavy steps.

## Public API

The stable top-level API is:

- `batsim.simulate_galaxy`
- `batsim.clear_backend_memory`
- `batsim.Stamp`
- `batsim.Transform` base class for custom coordinate transforms
- `batsim.LensTransform` and `batsim.AffineLensingTransform`
- `batsim.IaTransform` and `batsim.IATransform`
- `batsim.FlexionTransform`
- `batsim.experimental` for non-stable helpers such as SIP WCS parsing

`batsim.pltutil` remains available for legacy plotting utilities used by older
examples and validation notebooks.

## Building and Installing BATSim from Source

BATSim currently expects GalSim's C++ shared library to be available at build
time. The recommended route is to use conda or mamba with dependencies from
conda-forge.

First, clone the repository and switch to the repository root:

```bash
git clone https://github.com/CMacM/BATSim.git
cd BATSim
```

### Regular Install

Create and activate a build environment:

```bash
mamba create -n batsim -c conda-forge -c defaults python=3.10 conda-build boa
mamba activate batsim
```

Build the package:

```bash
conda mambabuild --override-channels -c conda-forge -c defaults conda/recipe
```

This creates an isolated build environment, installs dependencies, compiles the
BATSim C++ extension, links it to GalSim's C++ API, and produces a conda package
under:

```bash
$CONDA_PREFIX/conda-bld/linux-64/
```

Install the built package:

```bash
mamba install -c local batsim
```

If needed, install directly from the build artifact:

```bash
mamba install $CONDA_PREFIX/conda-bld/linux-64/batsim-*.tar.bz2
```

### Development Install

For development, install BATSim in editable mode so Python changes take effect
without reinstalling:

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
pip install -e .
```

Verify the installation:

```bash
python -c "import batsim; import batsim._gsinterface; print('BATSim installed successfully')"
```

## Running Benchmarks

BATSim uses ASV for performance benchmarks. Install the optional benchmark extra
into a working BATSim development environment:

```bash
pip install -e ".[benchmark]"
```

Then run discovery and a quick smoke benchmark:

```bash
cd benchmarks
asv check
asv run --quick --show-stderr
```

The default benchmark suite uses deterministic analytic GalSim profiles and
CPU-backed NumPy renders. Optional CuPy benchmarks skip automatically when a
working CUDA runtime is not available.

## Pip-Only Installation

A pure pip installation route is not currently provided because BATSim links
against GalSim's shared C++ library. If you build GalSim's C++ library manually,
set `GALSIM_LIB_DIR` and `GALSIM_INCLUDE_DIR` before installing BATSim. GalSim's
pip build notes are available in the
[GalSim documentation](https://galsim-developers.github.io/GalSim/_build/html/install_pip.html).

![BATSim Logo](./image/batsim_logo.png)

Repository icon created using images by Frepik and brgfx on Freepik.
