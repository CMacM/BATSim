# BATSim

BATSim is designed to enable the application of non-affine shear
transformations to GalSim images. It works by trasnforming the stamp pixel
grid, distorting the location of the pixel centers. When a GalSim image is
sampled onto this distorted grid, the resulting image will be sheared due to
the transformation of the pixel grid.

Non-affine transforms can simualate complex shear effects such as intrinsic
alignment, flexion, and optical field distrotion maps. Custom transform
functions can also be passed to Stamp objects.

## Installation

The package is not currently available on PyPI while in early development.

First clone the repository:
```shell
git clone https://github.com/CMacM/BATSim.git
```

Then install the package
```shell
cd BATSim
conda install --file requirements.txt -c conda-forge
pip install . --user
```

![BATSim Logo](./image/batsim_logo.png)

(repo icon created using images by Frepik and brgfx on Freepik)
