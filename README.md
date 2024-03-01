# BATSim

BATSim is designed to enable the application of non-affine shear transformations to GalSim images. It works by trasnforming the stamp pixel grid, distorting the location of the pixel centers. When a GalSim image is sampled onto this distorted grid, the resulting image will be sheared due to the transformation of the pixel grid. 

Non-affine transforms can simualate complex shear effects such as intrinsic alignment, flexion, and field distrotion maps. Custom transform functions can also be passed to Stamp objects.

## Installation

The package is not currently available on PyPI while in early development.

To install, first clone the repository:
```shell
git clone https://github.com/CMacM/BATSim.git
```

Next we need to build GalSim as a submodule and build its shared c++ library for BATSim to interface with. IMPORTANT: GalSim should be built from the submodule directory as specified, this is to ensure the GalSim python bindings are compiled against the same libraries as the BATSim libraries.
```shell
cd BATSim/extern/GalSim
pip install . --user
python setup.py build_shared_clib
```

Finally, we can change back to the BATSim parent directory and pip install.
```shell
cd ../..
pip install . --user
```

![BATSim Logo](batsim_logo.png)

(repo icon created using images by Frepik and brgfx on Freepik)
