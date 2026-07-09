import numpy as np
import time

from .backend import (
    _get_array_backend,
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
    force_input_flux=True,
    compensate_integration=True,
    use_true_center=True,
):
    """
    Simulate a galaxy image on a pixel grid.

    On failure, CuPy memory pools are cleared and a short exception is raised
    without retaining the large simulation frame in notebook tracebacks.
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
    compensate_integration=True,
    use_true_center=True,
):
    """
    Simulate a galaxy image on a pixel grid.

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

    supersample = _determine_supersampling(
        gal_obj,
        scale,
        integration_order,
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
        integration_order=integration_order,
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

    needs_fft = (
        psf_obj is not None or draw_method == "auto" or (compensate_integration and integration_order > 1)
    )

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
            compensate_integration=compensate_integration,
            center_index=center_index,
            profile=profile,
        )
        del gal_prof
    else:
        # Match the convention used by the convolution path:
        # extract the centred coarse image from the fine sampled image.
        sim_image = _extract_centered_coarse_image(
            xp=xp,
            image=gal_prof,
            downsample_ratio=fft_supersample,
            center_index=center_index,
        )
        del gal_prof

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        print(f"FFT and convolution took {end - start:.3f} seconds")

    sim_image = _center_crop_or_pad(sim_image, (output_ngrid, output_ngrid))

    return sim_image
