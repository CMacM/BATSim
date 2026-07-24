import numpy as np
import time

from .backend import (
    _precision_dtypes,
    _release_backend_memory,
    _resolve_array_backend,
    _to_numpy,
    clear_backend_memory,
    sync_if_gpu,
)
from . import sampling as _sampling
from .fft import (
    _apply_centering_phase,
    _apply_pixel_response,
    _center_crop_or_pad,
    _compensate_block_integration,
    _convolve_psf_fft,
    _draw_centered_psf,
    _extract_centered_coarse_image,
    _extract_centered_real_space_coarse_image,
    _rfft_centered_image,
)
from .grid import (
    FineGrid,
    _determine_supersampling,
    _integration_offsets,
    _make_fine_grid,
    _resolve_integration_sampling,
    _resolve_simulation_ngrid,
)
from .sampling import (
    _apply_prepared_transforms,
    _prepare_transforms,
    _require_gsinterface,
    _sample_psf_spectrum,
)
from .stamp import Stamp

_sample_galaxy_profile_impl = _sampling._sample_galaxy_profile


def _resolve_integration_compensation(compensate_integration):
    """Return the requested block-integration compensation mode."""
    if compensate_integration is None:
        return None

    if compensate_integration in ("quadrature", "exact_sinc"):
        return compensate_integration

    raise ValueError(
        "Invalid compensate_integration value "
        f"'{compensate_integration}'; use 'quadrature', 'exact_sinc', or None."
    )


def _needs_fourier_render(psf_obj, draw_method):
    """Return True when PSF or pixel convolution requires the FFT path."""
    return psf_obj is not None or draw_method == "auto"


def _sample_galaxy_profile(gal_obj, coords, fine_scale, real_dtype):
    """Compatibility wrapper for the sampling backend."""
    return _sample_galaxy_profile_impl(gal_obj, coords, fine_scale, real_dtype)


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

    This wrapper preserves legacy monkeypatching of
    ``batsim.sim._sample_galaxy_profile`` in internal tests and notebooks.
    """
    original_sampler = _sampling._sample_galaxy_profile
    _sampling._sample_galaxy_profile = _sample_galaxy_profile
    try:
        return _sampling._sample_galaxy_profile_on_stamp(
            gal_obj=gal_obj,
            stamp=stamp,
            transform_obj=transform_obj,
            fine_scale=fine_scale,
            real_dtype=real_dtype,
            xp=xp,
            integration_order=integration_order,
            profile=profile,
        )
    finally:
        _sampling._sample_galaxy_profile = original_sampler


def simulate_galaxy(
    gal_obj,
    scale=None,
    ngrid=None,
    transform_obj=None,
    psf_obj=None,
    draw_method="auto",
    safety=2.0,
    max_supersample=64,
    min_supersample=4,
    integration_order=2,
    max_fine_grid=4096,
    pad=16,
    precision="single",
    backend="np",
    profile=False,
    pix_scale=None,
    psf_mode="kvalue",
    force_input_flux=False,
    compensate_integration="quadrature",
    use_true_center=True,
):
    """
    Render a GalSim object through BATSim's non-affine sampling pipeline.

    The galaxy profile is sampled on a supersampled coordinate grid, optional
    coordinate transforms are applied before profile evaluation, and optional
    PSF/pixel convolution is performed in Fourier space. The returned image is
    a NumPy array regardless of the backend used internally.

    Parameters
    ----------
    gal_obj : galsim.GSObject
        Galaxy surface-brightness profile to sample.
    scale : float, optional
        Output pixel scale in arcsec.  Either ``scale`` or ``pix_scale`` must
        be supplied.  If both are supplied, they must agree.
    ngrid : int, optional
        Final square output size in pixels.  If omitted, BATSim chooses the
        output size from the GalSim object support at ``scale``.
    transform_obj : object or sequence of objects, optional
        Coordinate transform(s) applied before sampling the galaxy profile.
        Each transform must provide ``transform(coords)`` where ``coords`` has
        shape ``(2, npoints)``.  Transforms with a ``to_backend`` method
        are moved to the selected array backend before use.
    psf_obj : galsim.GSObject, optional
        PSF profile to convolve with the sampled galaxy.
    draw_method : {"auto", "no_pixel"}, optional
        Pixel-response mode.  ``"auto"`` applies the square-pixel response;
        ``"no_pixel"`` omits it.
    safety : float, optional
        Multiplicative safety factor used when estimating the required
        supersampling from the galaxy Fourier bandwidth and Nyquist frequency.
    max_supersample : int, optional
        Maximum total supersampling factor considered by the automatic grid
        selection.
    min_supersample : int, optional
        Minimum FFT supersampling factor used by the renderer. Reducing this
        below 4 is not recommended, as it can result in very corase representations
        of non-affine shear gradients and aliasing.
    integration_order : int, optional
        Gauss-Legendre quadrature order for sub-pixel block integration. 1 is
        equivalent to no sub-pixel integration.  Higher orders reduce
        high-frequency leakage at the cost of more backend evaluations.
    max_fine_grid : int or None, optional
        Maximum fine-grid side length.  Use ``None`` to disable the cap.
    pad : int, optional
        Padding, in coarse pixels, added around the internal simulation grid
        before supersampling.
    precision : {"single", "double"}, optional
        Floating-point precision for FFT-side arrays. Generally ``"single"`` is
        sufficient for most applications.
    backend : {"np", "numpy", "cp", "cupy"} or module or None, optional
        Array backend used for coordinate construction and FFT operations.
        ``None`` selects NumPy. Use ``"cp"`` or ``"cupy"`` to request CuPy.
    profile : bool, optional
        If True, print timing diagnostics for the major rendering stages.
    pix_scale : float, optional
        Deprecated spelling for ``scale`` retained for compatibility.
    psf_mode : {"kvalue", "real"}, optional
        PSF Fourier sampling mode.  ``"kvalue"`` evaluates the analytic PSF
        Fourier profile through the compiled backend; ``"real"`` draws the PSF
        in real space and FFTs it.
    force_input_flux : bool, optional
        If True, normalise the final convolved image to ``gal_obj.flux``.
    compensate_integration : {"quadrature", "exact_sinc"} or None, optional
        Fourier-space compensation applied when ``integration_order > 1``.
        ``"quadrature"`` removes the discrete Gauss-Legendre transfer function
        introduced by block integration. ``"exact_sinc"`` removes the ideal
        top-hat response instead. Using None leaves the block-integration
        smoothing in the rendered image and is not recommended for regular usage.
    use_true_center : bool, optional
        If True, align the fine grid to GalSim's true-image-center convention
        for the eventual coarse output grid. This is a functional mirror of the same
        argument in ``galsim.drawImage``.

    Returns
    -------
    ndarray
        Rendered square image with shape ``(ngrid, ngrid)`` when ``ngrid`` is
        provided, otherwise the automatically selected output shape.

    Raises
    ------
    RuntimeError
        Raised with a compact message if rendering fails after clearing cached
        backend memory.
    """
    clean_error = None
    cleanup_backend = None
    result = None

    try:
        result = _simulate_galaxy_impl(
            gal_obj=gal_obj,
            scale=scale,
            ngrid=ngrid,
            transform_obj=transform_obj,
            psf_obj=psf_obj,
            draw_method=draw_method,
            safety=safety,
            max_supersample=max_supersample,
            min_supersample=min_supersample,
            integration_order=integration_order,
            max_fine_grid=max_fine_grid,
            pad=pad,
            precision=precision,
            backend=backend,
            profile=profile,
            pix_scale=pix_scale,
            psf_mode=psf_mode,
            force_input_flux=force_input_flux,
            compensate_integration=compensate_integration,
            use_true_center=use_true_center,
        )
    except Exception as exc:
        err_type = type(exc).__name__
        err_msg = str(exc)

        try:
            cleanup_backend = _resolve_array_backend(backend)
        except Exception:
            pass

        clean_error = (
            f"simulate_galaxy failed cleanly after releasing backend memory " f"({err_type}: {err_msg})"
        )

    if cleanup_backend is not None:
        _release_backend_memory(cleanup_backend)

    if clean_error is not None:
        raise RuntimeError(clean_error) from None

    return result


def _simulate_galaxy_impl(
    gal_obj,
    scale=None,
    ngrid=None,
    transform_obj=None,
    psf_obj=None,
    draw_method="auto",
    safety=2.0,
    max_supersample=64,
    min_supersample=4,
    integration_order=2,
    max_fine_grid=4096,
    pad=16,
    precision="single",
    backend="np",
    profile=False,
    pix_scale=None,
    psf_mode="kvalue",
    force_input_flux=False,
    compensate_integration="quadrature",
    use_true_center=True,
):
    """
    Implement the render pipeline after public error handling is stripped away.

    The requested output ngrid controls only the final crop. The internal
    simulation grid is chosen from the galaxy/PSF support so that changing
    ngrid does not change the high-resolution simulation footprint.
    """
    if scale is None:
        if pix_scale is None:
            raise TypeError("simulate_galaxy requires scale or pix_scale.")
        scale = pix_scale
    elif pix_scale is not None and not np.isclose(scale, pix_scale):
        raise ValueError("scale and pix_scale were both supplied with different values.")

    if psf_mode not in ("real", "kvalue"):
        raise ValueError(f"Invalid psf_mode '{psf_mode}'; must be 'real' or 'kvalue'.")

    if draw_method not in ("auto", "no_pixel"):
        raise ValueError(f"Invalid draw_method '{draw_method}'; must be 'auto' or 'no_pixel'.")

    integration_compensation_mode = _resolve_integration_compensation(compensate_integration)

    xp = _resolve_array_backend(backend)

    dtypes = _precision_dtypes(xp, precision)
    real_dtype, _ = dtypes

    min_supersample = int(min_supersample)
    max_supersample = int(max_supersample)

    if min_supersample < 1:
        raise ValueError("min_supersample must be >= 1.")
    if max_supersample < 1:
        raise ValueError("max_supersample must be >= 1.")
    if max_supersample < min_supersample:
        max_supersample = min_supersample

    sim_ngrid = _resolve_simulation_ngrid(gal_obj, psf_obj, scale)
    output_ngrid = sim_ngrid if ngrid is None else int(ngrid)

    target_flux = gal_obj.flux if force_input_flux else None

    requested_integration_order = int(integration_order)
    if requested_integration_order < 1:
        raise ValueError("integration_order must be >= 1.")

    effective_integration_order = (
        requested_integration_order if _needs_fourier_render(psf_obj, draw_method) else 1
    )

    supersample = _determine_supersampling(
        gal_obj,
        scale,
        effective_integration_order,
        sim_ngrid=sim_ngrid,
        pad=pad,
        safety=safety,
        max_supersample=max_supersample,
        min_supersample=min_supersample,
        max_fine_grid=max_fine_grid,
    )

    requested_supersample = supersample
    fft_supersample, integration_order = _resolve_integration_sampling(
        supersample=requested_supersample,
        integration_order=effective_integration_order,
    )

    fft_supersample = max(fft_supersample, min_supersample)

    grid = _make_fine_grid(
        gal_obj=gal_obj,
        scale=scale,
        sim_ngrid=sim_ngrid,
        supersample=fft_supersample,
        pad=pad,
        max_fine_grid=max_fine_grid,
    )

    # Check if supersampling is too aggressive while preserving the requested
    # lower bound on the FFT supersampling factor.
    while (
        max_fine_grid is not None
        and grid.fine_compact / max_fine_grid >= 3
        and fft_supersample > min_supersample
    ):
        next_fft_supersample = max(min_supersample, fft_supersample // 2)
        if next_fft_supersample == fft_supersample:
            break

        fft_supersample = next_fft_supersample
        grid = _make_fine_grid(
            gal_obj=gal_obj,
            scale=scale,
            sim_ngrid=sim_ngrid,
            supersample=fft_supersample,
            pad=pad,
            max_fine_grid=max_fine_grid,
        )

    if profile:
        print(
            "Grid sizing: "
            f"sim_ngrid={sim_ngrid}, supersample={requested_supersample}, "
            f"integration_order={integration_order}, "
            f"fft_supersample={fft_supersample}, "
            f"fine_ngrid={grid.fine_ngrid}, fine_compact={grid.fine_compact}, "
            f"max_fine_grid={max_fine_grid}"
        )

    if profile:
        sync_if_gpu(xp)
    start = time.time()

    stamp = Stamp(
        nn=grid.fine_ngrid,
        scale=grid.fine_scale,
        backend=xp,
        dtype=real_dtype,
        use_true_center=use_true_center,
        downsample_ratio=fft_supersample,
    )
    center_index = stamp.center_index

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        print(f"Stamp construction took {end - start:.3f} seconds")

    gal_prof = _sample_galaxy_profile_on_stamp(
        gal_obj=gal_obj,
        stamp=stamp,
        transform_obj=transform_obj,
        fine_scale=grid.fine_scale,
        real_dtype=real_dtype,
        xp=xp,
        integration_order=integration_order,
        profile=profile,
    )
    del stamp

    if profile:
        sync_if_gpu(xp)
    start = time.time()

    needs_fft = _needs_fourier_render(psf_obj, draw_method)

    if needs_fft:
        sim_image = _convolve_psf_fft(
            xp=xp,
            gal_prof=gal_prof,
            target_flux=target_flux,
            scale=grid.fine_scale,
            pix_scale=scale,
            psf_obj=psf_obj,
            downsample_ratio=fft_supersample,
            draw_method=draw_method,
            dtypes=dtypes,
            psf_mode=psf_mode,
            integration_order=integration_order,
            integration_compensation_mode=integration_compensation_mode,
            center_index=center_index,
            profile=profile,
        )
        del gal_prof
    else:
        # The sampled image is still in centred real-space order because it has
        # not passed through the FFT path.
        sim_image = _extract_centered_real_space_coarse_image(
            xp=xp,
            image=gal_prof,
            downsample_ratio=fft_supersample,
        )
        del gal_prof

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        print(f"FFT and convolution took {end - start:.3f} seconds")

    sim_image = _center_crop_or_pad(sim_image, (output_ngrid, output_ngrid))
    if target_flux is not None and not needs_fft:
        sim_image *= target_flux / sim_image.sum()

    return sim_image
