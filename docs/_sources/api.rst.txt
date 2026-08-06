API Reference
=============

This page groups the stable BATSim API by how users normally look for it. The
generated module pages are still available at the bottom for lower-level
implementation details.

Rendering
---------

``simulate_galaxy`` is the main public entry point. It renders an analytic
GalSim profile on a BATSim coordinate grid, optionally applying a coordinate
transform and PSF convolution.

* :func:`batsim.simulate_galaxy`
* :func:`batsim.clear_backend_memory`

Image Stamps And Grids
----------------------

``Stamp`` owns the render grid, backend arrays, and sampling coordinates used by
the renderer. Most users will not need to construct stamps directly, but it is
the main lower-level object for debugging or building custom render workflows.

* :class:`batsim.Stamp`

Coordinate Transforms
---------------------

Transforms map output coordinates onto source-plane coordinates before profile
sampling. Built-in transforms cover affine lensing, radial intrinsic alignment,
and flexion. Subclass :class:`batsim.Transform` for custom non-affine fields.

* :class:`batsim.Transform`
* :class:`batsim.LensTransform`
* :class:`batsim.AffineLensingTransform`
* :class:`batsim.IaTransform`
* :class:`batsim.IATransform`
* :class:`batsim.FlexionTransform`

Supporting Modules
------------------

These modules are public enough to inspect and use for advanced work, but the
stable v1.0 surface is the top-level API above.

* :mod:`batsim.backend`
* :mod:`batsim.grid`
* :mod:`batsim.fft`
* :mod:`batsim.sampling`
* :mod:`batsim.pltutil`

Experimental Area
-----------------

Experimental helpers are available for development and validation work, but they
are not part of the stable top-level API.

* :mod:`batsim.experimental`
* :mod:`batsim.experimental.wcs`

Generated Module Pages
----------------------

.. toctree::
   :maxdepth: 1

   Raw generated package index <autoapi/batsim/index>
