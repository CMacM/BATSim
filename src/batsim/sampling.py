"""Real-space GalSim profile sampling helpers."""

import time

import numpy as np

from .backend import _to_numpy, sync_if_gpu
from .grid import _integration_offsets

try:
    from . import _gsinterface
except ImportError:
    _gsinterface = None


def _require_gsinterface():
    """Return the compiled C++ sampling backend, or raise a clear error."""
    if _gsinterface is None:
        raise ImportError(
            "batsim._gsinterface is not compiled. Reinstall BATSim with "
            "GalSim's C++ shared library built and linked. Instructions can be found here: "
            "https://galsim-developers.github.io/GalSim/_build/html/install_pip.html#installing-the-c-shared-library"
        )
    return _gsinterface


def _prepare_transforms(transform_obj, xp, dtype=None):
    """
    Move transform objects to the selected backend once per sampling call.
    """
    if transform_obj is None:
        return None

    if not isinstance(transform_obj, (list, tuple)):
        transform_obj = [transform_obj]

    prepared = []
    for trf in transform_obj:
        if hasattr(trf, "to_backend"):
            try:
                trf = trf.to_backend(backend=xp, dtype=dtype)
            except TypeError:
                trf = trf.to_backend(backend=xp)
        prepared.append(trf)

    return prepared


def _apply_prepared_transforms(coords, transform_obj, xp):
    """
    Apply backend-prepared transforms to a coordinate array.
    """
    coords = xp.asarray(coords)

    if transform_obj is None:
        return coords

    for trf in transform_obj:
        coords = trf.transform(coords)

    return coords


def _sample_galaxy_profile(gal_obj, coords, fine_scale, real_dtype):
    """
    Sample the galaxy surface brightness using the compiled C++ backend.
    """
    gsinterface = _require_gsinterface()
    use_single = np.dtype(real_dtype) == np.dtype(np.float32)

    sampler = gsinterface.getFluxVec32 if use_single else gsinterface.getFluxVec64

    return sampler(
        scale=fine_scale,
        gsobj=gal_obj._sbp,
        xy_coords=coords,
    )


def _sample_galaxy_profile_on_stamp(
    gal_obj,
    stamp,
    transform_obj,
    fine_scale,
    real_dtype,
    xp,
    integration_order=2,
    profile=False,
):
    """
    Sample the locally integrated galaxy profile on the fine simulation grid.
    """
    transform_obj = _prepare_transforms(
        transform_obj,
        xp=xp,
        dtype=real_dtype,
    )

    if profile:
        sync_if_gpu(xp)
    start = time.time()

    offsets, weights = _integration_offsets(
        xp=xp,
        fine_scale=fine_scale,
        integration_order=integration_order,
        dtype=real_dtype,
    )

    if offsets is None:
        coords_backend = stamp.coords
    else:
        n_offsets = offsets.shape[0]
        coords_backend = stamp.coords[None, :, :] + offsets[:, :, None]
        coords_backend = xp.transpose(coords_backend, (1, 0, 2))
        coords_backend = coords_backend.reshape(2, n_offsets * stamp.nn * stamp.nn)

    coords_backend = _apply_prepared_transforms(
        coords=coords_backend,
        transform_obj=transform_obj,
        xp=xp,
    )

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        print(f"Build coordinates and transforms took {end - start:.3f} seconds")

    if profile:
        sync_if_gpu(xp)
    start = time.time()

    coords = _to_numpy(xp, coords_backend)
    del coords_backend

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        print(f"Transfer coords to CPU took {end - start:.3f} seconds")

    start = time.time()

    gal_prof = _sample_galaxy_profile(
        gal_obj=gal_obj,
        coords=coords,
        fine_scale=fine_scale,
        real_dtype=real_dtype,
    )

    if integration_order > 1:
        gal_prof = np.asarray(gal_prof).reshape(
            integration_order * integration_order,
            stamp.nn,
            stamp.nn,
        )
        weights = _to_numpy(xp, weights).astype(np.dtype(real_dtype), copy=False)
        gal_prof = np.tensordot(
            weights,
            gal_prof,
            axes=(0, 0),
        ).astype(np.dtype(real_dtype), copy=False)

    end = time.time()
    if profile:
        print(f"Galaxy sampling took {end - start:.3f} seconds")

    del coords

    if profile:
        sync_if_gpu(xp)
    start = time.time()
    gal_prof = xp.asarray(gal_prof, dtype=real_dtype)

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        target = "GPU" if xp is not np else "NumPy"
        print(f"Transfer galaxy profile to {target} took {end - start:.3f} seconds")

    return gal_prof


def _sample_psf_spectrum(psf_obj, n, fine_scale, complex_dtype):
    """
    Sample the analytic PSF Fourier transform on an rFFT frequency grid.
    """
    gsinterface = _require_gsinterface()
    sampler = (
        gsinterface.getPsfKValue32
        if np.dtype(complex_dtype) == np.dtype(np.complex64)
        else gsinterface.getPsfKValue64
    )

    return sampler(
        scale=fine_scale,
        gsobj=psf_obj._sbp,
        n=int(n),
    )
