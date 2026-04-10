import galsim
import numpy as np
from functools import lru_cache

from . import _gsinterface
from .stamp import Stamp

from time import perf_counter

_PSF_REGISTRY = {}


def _register_psf(psf_obj):
    if len(_PSF_REGISTRY) >= 128:
        _PSF_REGISTRY.clear()
        _cached_effective_kernel.cache_clear()
    key = id(psf_obj)
    _PSF_REGISTRY[key] = psf_obj
    return key


@lru_cache(maxsize=128)
def _cached_effective_kernel(kernel_key):
    return _PSF_REGISTRY[kernel_key]


def _round_up_multiple(n, m=16):
    if n <= 0:
        return m
    return int(m * np.ceil(n / m))


def _choose_fine_sampling_scale(gobj, pix_scale, draw_method, kernel_obj):
    pix_scale = float(pix_scale)

    candidates = [pix_scale, float(gobj.nyquist_scale)]
    if kernel_obj is not None:
        try:
            candidates.append(float(kernel_obj.nyquist_scale))
        except Exception:
            pass

    target = min(candidates)

    if draw_method == "auto":
        target = min(target, pix_scale / 8.0)
        min_allowed = pix_scale / 16.0
    elif kernel_obj is not None:
        target = min(target, pix_scale / 4.0)
        min_allowed = pix_scale / 8.0
    else:
        min_allowed = pix_scale

    # temporarily force simulation at Pixel Nyquist scale for testing
    return 0.01


def _choose_downsample_ratio(pix_scale, fine_scale, max_ratio=24):
    raw_ratio = float(pix_scale) / float(fine_scale)
    ratio = int(2 ** np.ceil(np.log2(max(raw_ratio, 1.0))))
    return max(1, min(ratio, max_ratio))


def _choose_render_grid_size(
    gobj,
    fine_scale,
    truncate_ratio,
    maximum_num_grids,
    ngrid,
    force_ngrid,
    downsample_ratio,
):
    base_nn = int(np.ceil(gobj.getGoodImageSize(fine_scale) * float(truncate_ratio)))
    base_nn = max(base_nn, 1)

    if force_ngrid and ngrid is not None:
        base_nn = max(base_nn, int(ngrid) * int(downsample_ratio))

    # Small physical-support margin only.
    # This is NOT FFT padding.
    support_margin = max(128, int(np.ceil(0.05 * base_nn)))
    render_nn = base_nn + 2 * support_margin

    render_nn = _round_up_multiple(render_nn, 16)
    render_nn = min(render_nn, int(maximum_num_grids))

    if render_nn < base_nn + 2 * support_margin:
        raise ValueError(
            "Required render size exceeds maximum_num_grids. "
            f"Required at least {base_nn + 2 * support_margin}, "
            f"got maximum_num_grids={maximum_num_grids}."
        )

    return int(render_nn)


def _prepare_kernel(psf_obj, pix_scale, draw_method):
    if draw_method == "auto":
        pixel = galsim.Pixel(scale=pix_scale)
        if psf_obj is None:
            return pixel
        return galsim.Convolve([psf_obj, pixel])

    if draw_method == "no_pixel":
        return psf_obj

    raise ValueError(f"do not support draw_method={draw_method}")


def simulate_galaxy(
    gal_obj,
    pix_scale,
    ngrid=None,
    transform_obj=None,
    psf_obj=None,
    truncate_ratio=1.0,
    maximum_num_grids=4096,
    draw_method="auto",
    force_ngrid=False,
    delta_image_x=0.0,
    delta_image_y=0.0,
    profile=False,
    fft_pad_pixels=16,
):
    """
    Render a galaxy using a compact fine-grid sampling stage in Python and
    internal FFT padding in C++.

    GalSim-like semantics:
    - draw_method="auto"     -> include pixel response via galsim.Pixel(scale=pix_scale)
    - draw_method="no_pixel" -> omit pixel response

    The Python render grid is intentionally unpadded for convolution.
    Any guard padding needed for FFT convolution is added only inside C++.
    """
    def _log(msg):
        if profile:
            print(f"[simulate_galaxy] {msg}")

    if truncate_ratio <= 0:
        raise ValueError("truncate_ratio must be > 0")
    if maximum_num_grids <= 0:
        raise ValueError("maximum_num_grids must be > 0")
    if ngrid is not None and int(ngrid) <= 0:
        raise ValueError("ngrid must be > 0 when provided")
    if fft_pad_pixels < 0:
        raise ValueError("fft_pad_pixels must be >= 0")

    t_total = perf_counter() if profile else None

    # Apply requested image-plane shift in detector-pixel units.
    t = perf_counter() if profile else None
    gobj = gal_obj.shift(
        float(delta_image_x) * float(pix_scale),
        float(delta_image_y) * float(pix_scale),
    )
    if profile:
        _log(f"apply_shift={perf_counter() - t:.4e}s")

    # Build GalSim-like effective kernel.
    t = perf_counter() if profile else None
    kernel_obj = _prepare_kernel(
        psf_obj=psf_obj,
        pix_scale=pix_scale,
        draw_method=draw_method,
    )
    kernel_key = _register_psf(kernel_obj) if kernel_obj is not None else None
    if profile:
        _log(f"prepare_kernel={perf_counter() - t:.4e}s")

    # Choose fine sampling scale.
    t = perf_counter() if profile else None
    fine_scale = _choose_fine_sampling_scale(
        gobj=gobj,
        pix_scale=pix_scale,
        draw_method=draw_method,
        kernel_obj=kernel_obj,
    )
    downsample_ratio = _choose_downsample_ratio(
        pix_scale=pix_scale,
        fine_scale=fine_scale,
        max_ratio=4,
    )
    fine_scale = float(pix_scale) / float(downsample_ratio)
    if profile:
        _log(
            f"choose_sampling={perf_counter() - t:.4e}s "
            f"fine_scale={fine_scale:.4e} "
            f"downsample_ratio={downsample_ratio}"
        )

    # Choose compact unpadded render size.
    t = perf_counter() if profile else None
    render_nn = _choose_render_grid_size(
        gobj=gobj,
        fine_scale=fine_scale,
        truncate_ratio=truncate_ratio,
        maximum_num_grids=maximum_num_grids,
        ngrid=ngrid,
        force_ngrid=force_ngrid,
        downsample_ratio=downsample_ratio,
    )
    if profile:
        _log(f"choose_render_grid={perf_counter() - t:.4e}s render_nn={render_nn}")

    # Build only the physically relevant stamp in Python.
    t = perf_counter() if profile else None
    stamp = Stamp(nn=render_nn, scale=fine_scale)

    if isinstance(transform_obj, list):
        gal_coords = stamp.coords
        for trf in transform_obj:
            gal_coords = trf.transform(gal_coords)
    elif transform_obj is not None:
        gal_coords = transform_obj.transform(stamp.coords)
    else:
        gal_coords = stamp.coords

    if profile:
        _log(f"coords_transform={perf_counter() - t:.4e}s")

    # Sample the galaxy on the compact fine grid.
    t = perf_counter() if profile else None
    gal_prof = _gsinterface.getFluxVec(
        scale=fine_scale,
        gsobj=gobj._sbp,
        xy_coords=gal_coords,
    )
    if profile:
        _log(f"sample_flux={perf_counter() - t:.4e}s")

    # Convolve in C++, which adds its own zero-padding right before the FFT.
    if kernel_obj is not None:
        t = perf_counter() if profile else None
        eff_kernel = _cached_effective_kernel(kernel_key)
        gal_prof = _gsinterface.convolvePsfFine(
            scale=fine_scale,
            psf=eff_kernel._sbp,
            gal_prof=gal_prof,
            pad_pixels=int(fft_pad_pixels),
        )
        if profile:
            _log(f"fine_convolution={perf_counter() - t:.4e}s")

    # Choose final output dimension.
    t = perf_counter() if profile else None
    if ngrid is None:
        out_dim = max(1, int(np.round(render_nn * fine_scale / pix_scale)))
    else:
        out_dim = int(ngrid)
    if profile:
        _log(f"choose_output_grid={perf_counter() - t:.4e}s out_dim={out_dim}")

    # Only resample when scales differ. Otherwise crop/pad cheaply.
    t = perf_counter() if profile else None
    if np.isclose(fine_scale, pix_scale):
        if gal_prof.shape[0] != out_dim:
            gal_prof = _gsinterface.centerCropOrPad(gal_prof, out_dim)
    else:
        gal_prof = _gsinterface.resampleToGrid(
            image=gal_prof,
            in_scale=fine_scale,
            out_scale=pix_scale,
            out_dim=out_dim,
        )
    if profile:
        _log(f"final_output_stage={perf_counter() - t:.4e}s")

    if profile:
        _log(
            "stats "
            f"render_nn={render_nn} "
            f"fine_scale={fine_scale:.4e} "
            f"downsample_ratio={downsample_ratio} "
            f"fft_pad_pixels={fft_pad_pixels} "
            f"draw_method={draw_method} "
            f"out_dim={out_dim}"
        )
        _log(f"total={perf_counter() - t_total:.4e}s")

    return gal_prof