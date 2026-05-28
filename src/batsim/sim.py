import warnings
import gc
from dataclasses import dataclass

import galsim
import numpy as np
import time
import types

from .stamp import Stamp

try:
    from . import _gsinterface
except ImportError:
    _gsinterface = None


_ARRAY_BACKEND = None
_CPU_FALLBACK_WARNED = False


@dataclass(frozen=True)
class FineGrid:
    """
    Description of the high-resolution grid used for real-space sampling
    and Fourier-space convolution.
    """
    fine_scale: float
    fine_ngrid: int
    fine_compact: int


def _require_gsinterface():
    """Return the compiled C++ sampling backend, or raise a clear error."""
    if _gsinterface is None:
        raise ImportError(
            "batsim._gsinterface is not compiled. Reinstall BATSim with "
            "GalSim's C++ shared library built and linked. Instructions can be found here: "
            "https://galsim-developers.github.io/GalSim/_build/html/install_pip.html#installing-the-c-shared-library"
        )
    return _gsinterface


def _get_array_backend():
    """
    Return CuPy if available, otherwise NumPy.

    The selected module is cached after the first call.
    """
    global _ARRAY_BACKEND
    global _CPU_FALLBACK_WARNED

    if _ARRAY_BACKEND is not None:
        return _ARRAY_BACKEND

    try:
        import cupy as cp

        cp.cuda.runtime.getDeviceCount()
        _ARRAY_BACKEND = cp
    except Exception:
        if not _CPU_FALLBACK_WARNED:
            warnings.warn(
                "CuPy unavailable; falling back to NumPy FFTs on CPU.",
                RuntimeWarning,
                stacklevel=2,
            )
            _CPU_FALLBACK_WARNED = True
        _ARRAY_BACKEND = np

    return _ARRAY_BACKEND


def _resolve_array_backend(backend):
    """
    Resolve a user-supplied backend selector to an array module.
    """
    if backend is None:
        return _get_array_backend()

    if backend is np or backend in ("np", "numpy"):
        return np

    if backend in ("cp", "cupy"):
        try:
            import cupy as cp

            cp.cuda.runtime.getDeviceCount()
            return cp
        except Exception as exc:
            raise RuntimeError(
                "backend='cp' was requested, but CuPy is unavailable."
            ) from exc

    if isinstance(backend, types.ModuleType):
        return backend

    raise ValueError(
        f"Invalid backend '{backend}'; must be 'np', 'cp', a module, or None."
    )


def _to_numpy(xp, array):
    """Convert a NumPy/CuPy array to a NumPy array."""
    if xp is np:
        return np.asarray(array)
    return xp.asnumpy(array)


def _precision_dtypes(xp, precision):
    """Return real and complex dtypes for the requested FFT precision."""
    if precision == "single":
        return xp.float32, xp.complex64

    if precision == "double":
        return xp.float64, xp.complex128

    raise ValueError(
        f"Invalid precision '{precision}'; must be 'single' or 'double'."
    )


def _release_backend_memory(xp):
    """
    Release cached CuPy allocations after a simulation.
    """
    if xp is np:
        return

    gc.collect()

    try:
        xp.cuda.Stream.null.synchronize()
    except Exception:
        pass

    try:
        xp.fft.config.get_plan_cache().clear()
    except Exception:
        pass

    gc.collect()

    try:
        xp.get_default_memory_pool().free_all_blocks()
    except Exception:
        pass

    try:
        xp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass

    gc.collect()


def clear_backend_memory(backend=None):
    """
    Clear cached backend allocations, primarily useful after CuPy OOMs.
    """
    xp = _resolve_array_backend(backend)
    _release_backend_memory(xp)


def _center_crop_or_pad(array, target_shape):
    """
    Centrally crop or zero-pad an array to the requested shape.
    """
    out = np.zeros(target_shape, dtype=array.dtype)

    src_y, src_x = array.shape
    dst_y, dst_x = target_shape

    copy_y = min(src_y, dst_y)
    copy_x = min(src_x, dst_x)

    src_y0 = (src_y - copy_y) // 2
    src_x0 = (src_x - copy_x) // 2
    dst_y0 = (dst_y - copy_y) // 2
    dst_x0 = (dst_x - copy_x) // 2

    out[dst_y0:dst_y0 + copy_y, dst_x0:dst_x0 + copy_x] = array[
        src_y0:src_y0 + copy_y,
        src_x0:src_x0 + copy_x,
    ]

    return out


def _draw_centered_psf(xp, psf_obj, n, scale, dtype=None):
    """
    Draw the PSF on the full FFT grid.

    Pixel convolution is not included here; the pixel response is applied
    separately in Fourier space.
    """
    psf_image = psf_obj.drawImage(
        nx=n,
        ny=n,
        scale=scale,
        method="no_pixel",
        use_true_center=False,
    ).array

    return xp.asarray(psf_image, dtype=dtype)


def _determine_supersampling(
    obj,
    scale,
    integration_order,
    safety=2.0,
    max_supersample=16,
    min_supersample=4,
    attenuation_target=0.5
):
    """
    Estimate supersampling factor from the object's Fourier bandwidth.
    """
    maxk = obj.maxk
    nyquist = np.pi / scale

    raw = safety * maxk / nyquist

    q = max(1, int(integration_order))

    if q <= 1:
        effective = raw
    else:
        # attenuation_target lets you tune how much high-k leakage you tolerate.
        integration_reduction = q ** 2 / attenuation_target**0.5

        # Avoid absurd reductions for high q.
        integration_reduction = min(integration_reduction, 8.0 * q)

        effective = raw / integration_reduction

        # Do not allow the requested effective supersampling to fall below q,
        # otherwise fft_supersample may collapse too aggressively.
        effective = max(effective, q)
        
    supersample = 2 ** int(np.ceil(np.log2(max(min_supersample, effective))))

    return int(np.clip(supersample, min_supersample, max_supersample))


def _resolve_simulation_ngrid(gal_obj, psf_obj, scale):
    """
    Choose the internal coarse-grid simulation size.

    This is deliberately independent of the requested final output ngrid.
    """
    if psf_obj is None:
        return int(gal_obj.getGoodImageSize(scale))

    conv = galsim.Convolve([gal_obj, psf_obj])
    return int(conv.getGoodImageSize(scale))


def _make_fine_grid(gal_obj, scale, sim_ngrid, supersample, pad, max_fine_grid=None):
    """
    Construct the fine grid used for sampling and convolution.
    """
    fine_scale = scale / supersample

    requested_fine_ngrid = (sim_ngrid + 2 * pad) * supersample
    fine_compact = int(gal_obj.getGoodImageSize(fine_scale))

    fine_ngrid = max(requested_fine_ngrid, fine_compact)

    remainder = fine_ngrid % supersample
    if remainder:
        fine_ngrid += supersample - remainder

    if max_fine_grid is not None:
        max_fine_grid = int(max_fine_grid)
        if max_fine_grid < 1:
            raise ValueError("max_fine_grid must be positive.")

        capped_fine_ngrid = max_fine_grid - (max_fine_grid % supersample)
        if capped_fine_ngrid < supersample:
            raise ValueError(
                "max_fine_grid must be at least one supersampling block."
            )

        fine_ngrid = min(fine_ngrid, capped_fine_ngrid)

    return FineGrid(
        fine_scale=fine_scale,
        fine_ngrid=fine_ngrid,
        fine_compact=fine_compact,
    )


def _resolve_integration_sampling(supersample, integration_order):
    """
    Split the requested anti-aliasing factor between integration and FFT size.
    """
    supersample = int(supersample)
    integration_order = int(integration_order)

    if supersample < 1:
        raise ValueError("supersample must be >= 1.")
    if integration_order < 1:
        raise ValueError("integration_order must be >= 1.")

    integration_order = min(integration_order, supersample)
    fft_supersample = int(np.ceil(supersample / integration_order))

    return fft_supersample, integration_order


def _integration_offsets(xp, fine_scale, integration_order, dtype):
    """
    Return sub-pixel offsets for midpoint block integration.
    """
    if integration_order <= 1:
        return None

    offsets_1d = (
        (xp.arange(integration_order, dtype=dtype) + 0.5) / integration_order
        - 0.5
    ) * fine_scale
    yy, xx = xp.meshgrid(offsets_1d, offsets_1d, indexing="ij")

    return xp.stack(
        [
            xx.ravel(),
            yy.ravel(),
        ],
        axis=1,
    )


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

    offsets = _integration_offsets(
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
        ).mean(axis=0, dtype=np.dtype(real_dtype))

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


def _rfft_centered_image(xp, image, real_dtype, complex_dtype):
    """
    Compute rfft2(ifftshift(image)) without explicitly shifting the image.

    The input image remains in the centred real-space convention. The phase
    correction converts the spectrum to the convention that would have been
    obtained by first applying np.fft.ifftshift(image).
    """
    n = image.shape[0]

    image = xp.asarray(image, dtype=real_dtype)
    image = xp.ascontiguousarray(image)
    spectrum = xp.fft.rfft2(image).astype(complex_dtype, copy=False)

    del image

    _apply_centering_phase(xp, spectrum, n)

    return spectrum


def _apply_centering_phase(xp, spectrum, n):
    """
    Apply the Fourier phase equivalent to np.fft.ifftshift before an FFT.

    This avoids materialising a shifted real-space image, which is expensive
    for very large supersampled stamps.
    """
    shift = -(n // 2)

    ky = xp.fft.fftfreq(n) * n
    kx = xp.fft.rfftfreq(n) * n

    phase_y = xp.exp(-2j * xp.pi * ky * shift / n)
    phase_x = xp.exp(-2j * xp.pi * kx * shift / n)

    phase_y = phase_y.astype(spectrum.dtype, copy=False)
    phase_x = phase_x.astype(spectrum.dtype, copy=False)

    spectrum *= phase_y[:, None]
    spectrum *= phase_x[None, :]


def _apply_pixel_response(xp, spectrum, n, fine_scale, pix_scale):
    """
    Apply the square-pixel Fourier response in-place on an rfft2 grid.

    The response is separable, so this avoids allocating a full 2D pixel FFT.
    """
    ky = 2.0 * xp.pi * xp.fft.fftfreq(n, d=fine_scale)
    kx = 2.0 * xp.pi * xp.fft.rfftfreq(n, d=fine_scale)

    pix_y = xp.sinc(ky * pix_scale / (2.0 * xp.pi))
    pix_x = xp.sinc(kx * pix_scale / (2.0 * xp.pi))

    real_dtype = spectrum.real.dtype
    pix_y = pix_y.astype(real_dtype, copy=False)
    pix_x = pix_x.astype(real_dtype, copy=False)

    spectrum *= pix_y[:, None]
    spectrum *= pix_x[None, :]


def _extract_centered_coarse_image(xp, image, downsample_ratio):
    """
    Return fftshift(image)[::downsample_ratio, ::downsample_ratio] without
    explicitly materialising fftshift(image).

    This keeps the expensive full-resolution image shift out of both GPU and
    CPU memory paths. Only the already-downsampled coarse image is gathered.
    """
    n = image.shape[0]
    s = int(downsample_ratio)

    if s < 1:
        raise ValueError("downsample_ratio must be >= 1")

    coarse_n = n // s

    rows = (xp.arange(coarse_n) * s - n // 2) % n
    cols = (xp.arange(coarse_n) * s - n // 2) % n

    coarse = image[rows[:, None], cols[None, :]]
    del image, rows, cols

    if s > 1:
        coarse = (s ** 2) * coarse

    return _to_numpy(xp, coarse)


def _convolve_psf_fft(
    xp,
    gal_prof,
    scale,
    pix_scale,
    psf_obj,
    downsample_ratio,
    draw_method,
    dtypes,
    psf_mode,
    target_flux=None,
    profile=False
):
    """
    Convolve the sampled galaxy with PSF and/or pixel response using rFFT.

    The expensive path avoids explicit fftshift/ifftshift operations. Instead:
    - centred real-space images are converted by Fourier phase correction;
    - the pixel response is applied separably in-place;
    - the full fine-grid image is inverse FFT'd, then downsampled in real space.
    """
    n = gal_prof.shape[0]
    real_dtype, complex_dtype = dtypes

    spectrum = _rfft_centered_image(
        xp=xp,
        image=gal_prof,
        real_dtype=real_dtype,
        complex_dtype=complex_dtype,
    )

    if psf_obj is not None and psf_mode == "real":
        start = time.time() if profile else None
        psf_image = _draw_centered_psf(
            xp=xp,
            psf_obj=psf_obj,
            n=n,
            scale=scale,
            dtype=real_dtype,
        )
        if profile:
            sync_if_gpu(xp)
            end = time.time()
            print(f"Drawing and transferring PSF took {end - start:.3f} seconds")

        psf_spectrum = _rfft_centered_image(
            xp=xp,
            image=psf_image,
            real_dtype=real_dtype,
            complex_dtype=complex_dtype,
        )

        spectrum *= psf_spectrum
        del psf_image, psf_spectrum

    elif psf_obj is not None and psf_mode == "kvalue":
        start = time.time() if profile else None
        psf_spectrum = _sample_psf_spectrum(
            psf_obj=psf_obj,
            n=n,
            fine_scale=scale,
            complex_dtype=complex_dtype,
        )

        if profile:
            end = time.time()
            print(f"Sampling PSF spectrum took {end - start:.3f} seconds")

        start = time.time() if profile else None
        psf_spectrum = xp.asarray(psf_spectrum, dtype=complex_dtype)
        if profile:
            sync_if_gpu(xp)
            end = time.time()
            target = "GPU" if xp is not np else "NumPy"
            print(f"Transfer PSF spectrum to {target} took {end - start:.3f} seconds")

        spectrum *= psf_spectrum
        del psf_spectrum

    elif psf_obj is not None:
        raise ValueError(
            f"Invalid psf_mode '{psf_mode}'; must be 'real' or 'kvalue'."
        )

    if draw_method == "auto":
        _apply_pixel_response(
            xp=xp,
            spectrum=spectrum,
            n=n,
            fine_scale=scale,
            pix_scale=pix_scale,
        )

    image = xp.fft.irfft2(spectrum, s=(n, n))
    del spectrum

    # Normalise flux based on target input flux
    if target_flux is not None:
        image *= target_flux / image.sum()

    coarse_image = _extract_centered_coarse_image(
        xp=xp,
        image=image,
        downsample_ratio=downsample_ratio,
    )
    del image

    return coarse_image


def simulate_galaxy(
    gal_obj,
    scale=None,
    ngrid=None,
    transform_obj=None,
    psf_obj=None,
    draw_method="auto",
    safety=2.0,
    max_supersample=32,
    min_supersample=4,
    integration_order=2,
    max_fine_grid=4096,
    pad=16,
    precision="single",
    backend=None,
    profile=False,
    pix_scale=None,
    psf_mode="kvalue",
    force_input_flux=True
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
            force_input_flux=force_input_flux
        )
    except Exception as exc:
        err_type = type(exc).__name__
        err_msg = str(exc)

        try:
            cleanup_backend = _resolve_array_backend(backend)
        except Exception:
            pass

        clean_error = (
            f"simulate_galaxy failed cleanly after releasing backend memory "
            f"({err_type}: {err_msg})"
        )

    try:
        cleanup_backend = _resolve_array_backend(backend)
    except Exception:
        cleanup_backend = None

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
    max_supersample=32,
    min_supersample=4,
    integration_order=2,
    max_fine_grid=4096,
    pad=16,
    precision="single",
    backend=None,
    profile=False,
    pix_scale=None,
    psf_mode="kvalue",
    force_input_flux=False
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
        raise ValueError(
            f"Invalid psf_mode '{psf_mode}'; must be 'real' or 'kvalue'."
        )

    if draw_method not in ("auto", "no_pixel"):
        raise ValueError(
            f"Invalid draw_method '{draw_method}'; must be 'auto' or 'no_pixel'."
        )

    xp = _resolve_array_backend(backend)

    dtypes = _precision_dtypes(xp, precision)
    real_dtype, _ = dtypes

    sim_ngrid = _resolve_simulation_ngrid(gal_obj, psf_obj, scale)
    output_ngrid = sim_ngrid if ngrid is None else int(ngrid)

    target_flux = gal_obj.flux if force_input_flux else None

    supersample = _determine_supersampling(
        gal_obj,
        scale,
        integration_order,
        safety=safety,
        max_supersample=max_supersample,
        min_supersample=min_supersample
    )

    requested_supersample = supersample
    fft_supersample, integration_order = _resolve_integration_sampling(
        supersample=requested_supersample,
        integration_order=integration_order,
    )

    grid = _make_fine_grid(
        gal_obj=gal_obj,
        scale=scale,
        sim_ngrid=sim_ngrid,
        supersample=fft_supersample,
        pad=pad,
        max_fine_grid=max_fine_grid,
    )

    # Check if supersampling is too aggresive
    while grid.fine_compact / max_fine_grid > 4:
        fft_supersample = fft_supersample // 2
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
        dtype=real_dtype
    )

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

    if psf_obj is not None or draw_method == "auto":
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
        )
        del gal_prof

    if profile:
        sync_if_gpu(xp)
        end = time.time()
        print(f"FFT and convolution took {end - start:.3f} seconds")

    sim_image = _center_crop_or_pad(sim_image, (output_ngrid, output_ngrid))

    return sim_image


def sync_if_gpu(xp):
    if xp is not np:
        xp.cuda.Stream.null.synchronize()
