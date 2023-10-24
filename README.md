# BATSim

BATSim is designed to enable the application of non-affine shear transformations to GalSim images. It works by trasnforming the stamp pixel grid, distroting the location of the pixel centers. When a GalSim image is sampled on this distorted grid, the resulting image will be sheared due to the transformation of the pixel grid. 

Non-affine transforms can simualate complex shear effects such as intrinsic alignment, flexion, and field distrotion maps. Custom transform functions can also be passed to Stamp objects.

## Installation

The package is not currently available on PyPI while in early development.

To install, first clone the repository:
```shell
git clone https://github.com/CMacM/BATSim.git
```
Then change to the cloned directiory and run pip install to build from the setup.toml file.
```shell
cd BATSim
pip install -e . --user
```

![BATSim Logo](batsim_logo.png)

(repo icon created using images by Frepik and brgfx on Freepik)
