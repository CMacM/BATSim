BATSim
======

BATSim renders GalSim profiles on coordinate grids that can be transformed
before sampling. This supports non-affine shear fields, including
intrinsic-alignment-style radial shear and flexion, while preserving the current
GalSim image-centering and flux conventions used by the renderer.

Current Rendering Defaults
--------------------------

The public render entry point is :func:`batsim.simulate_galaxy`. Its current
default pipeline uses:

* ``psf_mode="kvalue"`` to sample the analytic PSF Fourier profile through the
  compiled BATSim/GalSim bridge.
* ``force_input_flux=True`` to preserve the input galaxy flux after convolution.
* ``use_true_center=True`` to align the fine grid with GalSim's true image
  center convention.
* ``integration_order=2`` to use Gauss-Legendre block integration before the
  FFT convolution path.
* ``compensate_integration="quadrature"`` to remove the matching
  Gauss-Legendre transfer function. Pass ``"exact_sinc"`` for ideal top-hat
  compensation, or None to disable compensation.
* ``backend="np"`` by default. Passing ``backend=None`` also selects NumPy;
  pass ``backend="cp"`` to request CuPy.

Quickstart
----------

.. code-block:: python

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

``image`` is returned as a NumPy array, even when CuPy is used internally for
FFT-heavy work.

Installation Notes
------------------

BATSim links against GalSim's shared C++ library. The recommended installation
route is currently a conda or mamba environment using conda-forge packages for
GalSim, Eigen, pybind11, and NumPy, followed by ``pip install -e .`` for local
development. See the repository README for full build instructions.

Public API
----------

The stable top-level API includes:

* :func:`batsim.simulate_galaxy`
* :func:`batsim.clear_backend_memory`
* :class:`batsim.Stamp`
* :class:`batsim.LensTransform` and :class:`batsim.AffineLensingTransform`
* :class:`batsim.IaTransform` and :class:`batsim.IATransform`
* :class:`batsim.FlexionTransform`

Experimental helpers live under :mod:`batsim.experimental`. They are available
for development and validation work, but are not part of the stable top-level
API.


.. toctree::
   :hidden:

   Home page <self>
   API Reference <autoapi/index>
   Notebooks <notebooks>
