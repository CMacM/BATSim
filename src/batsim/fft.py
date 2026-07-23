"""Fourier-space convolution and image extraction helpers."""

import time

import numpy as np

from .backend import _to_numpy, sync_if_gpu
from .sampling import _sample_psf_spectrum


def _center_crop_or_pad(array, target_shape):
    """
    Centrally crop or zero-pad an array to the requested shape.

    This is used after the internal simulation grid has been rendered so the
    caller's requested output size is a centered view of the same simulation.
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

    out[dst_y0 : dst_y0 + copy_y, dst_x0 : dst_x0 + copy_x] = array[
        src_y0 : src_y0 + copy_y,
        src_x0 : src_x0 + copy_x,
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


def _rfft_centered_image(xp, image, real_dtype, complex_dtype, center_index=None):
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

    _apply_centering_phase(xp, spectrum, n, center_index=center_index)

    return spectrum


def _apply_centering_phase(xp, spectrum, n, center_index=None):
    """
    Apply the Fourier phase equivalent to np.fft.ifftshift before an FFT.

    This avoids materialising a shifted real-space image, which is expensive
    for very large supersampled stamps.
    """
    if center_index is None:
        center_index = n // 2

    shift = -float(center_index)

    # Applying the separable phase is equivalent to ifftshift before rfft2,
    # without moving the large real-space array in memory.
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


def _block_integration_response(xp, k, fine_scale, mode="quadrature", integration_order=None):
    """Return the 1D Fourier response introduced by block integration."""
    if mode == "exact_sinc":
        return xp.sinc(k * fine_scale / (2.0 * xp.pi))

    if mode == "quadrature":
        if integration_order is None or integration_order <= 1:
            raise ValueError("integration_order > 1 is required for quadrature compensation.")

        nodes, weights = np.polynomial.legendre.leggauss(int(integration_order))
        offsets = xp.asarray(0.5 * nodes * fine_scale, dtype=k.dtype)
        weights = xp.asarray(0.5 * weights, dtype=k.dtype)

        phase = k[:, None] * offsets[None, :]
        response = xp.sum(weights[None, :] * xp.cos(phase), axis=1)
        return response

    raise ValueError(
        "Invalid block integration compensation mode " f"'{mode}'; must be 'quadrature' or 'exact_sinc'."
    )


def _compensate_block_integration(
    xp,
    spectrum,
    n,
    fine_scale,
    min_response=1.0e-3,
    mode="quadrature",
    integration_order=None,
):
    """
    Remove the fine-pixel averaging response introduced by block integration.

    Gauss-Legendre block integration approximates the mean surface brightness
    over each fine-grid pixel, introducing a separable pre-smoothing response.
    Divide out this response within the represented Fourier band before
    applying the physical PSF convolution.

    ``mode="quadrature"`` removes the discrete Gauss-Legendre transfer
    function for the configured nodes and weights. ``mode="exact_sinc"`` keeps
    the ideal top-hat response compensation available as an explicit option.
    """
    ky = 2.0 * xp.pi * xp.fft.fftfreq(n, d=fine_scale)
    kx = 2.0 * xp.pi * xp.fft.rfftfreq(n, d=fine_scale)

    response_y = _block_integration_response(
        xp=xp,
        k=ky,
        fine_scale=fine_scale,
        mode=mode,
        integration_order=integration_order,
    )
    response_x = _block_integration_response(
        xp=xp,
        k=kx,
        fine_scale=fine_scale,
        mode=mode,
        integration_order=integration_order,
    )

    real_dtype = spectrum.real.dtype
    min_response = xp.asarray(min_response, dtype=real_dtype)
    response_y = response_y.astype(real_dtype, copy=False)
    response_x = response_x.astype(real_dtype, copy=False)

    response_y = xp.maximum(response_y, min_response)
    response_x = xp.maximum(response_x, min_response)

    spectrum /= response_y[:, None]
    spectrum /= response_x[None, :]


def _extract_centered_coarse_image(xp, image, downsample_ratio, center_index=None):
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
    if center_index is None:
        center_index = n // 2

    # The wrapped center must fall on an integer pixel in FFT output order so
    # strided extraction matches fftshift(image)[::s, ::s].
    wrapped_center = -float(center_index)
    if not np.isclose(wrapped_center, round(wrapped_center)):
        raise ValueError(
            "The fine-grid center is not aligned with integer FFT output pixels. "
            "Use an even downsample_ratio, or use_true_center=False."
        )

    wrapped_center = int(round(wrapped_center))

    rows = (xp.arange(coarse_n) * s + wrapped_center) % n
    cols = (xp.arange(coarse_n) * s + wrapped_center) % n

    coarse = image[rows[:, None], cols[None, :]]
    del image, rows, cols

    if s > 1:
        coarse = (s**2) * coarse

    return _to_numpy(xp, coarse)


def _extract_centered_real_space_coarse_image(xp, image, downsample_ratio):
    """
    Extract the aligned coarse grid from a centred real-space fine image.

    Unlike ``_extract_centered_coarse_image``, this is for images that have not
    passed through an inverse FFT and therefore are already in spatial order.
    """
    n = image.shape[0]
    s = int(downsample_ratio)

    if s < 1:
        raise ValueError("downsample_ratio must be >= 1")
    if n % s != 0:
        raise ValueError("fine image size must be divisible by downsample_ratio")

    coarse = image[::s, ::s]

    if s > 1:
        coarse = (s**2) * coarse

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
    integration_order=1,
    integration_compensation_mode="quadrature",
    center_index=None,
    target_flux=None,
    profile=False,
):
    """
    Convolve the sampled galaxy with PSF and/or pixel response using rFFT.

    The expensive path avoids explicit fftshift/ifftshift operations. Instead:
    - centred real-space images are converted by Fourier phase correction;
    - the pixel response is applied separably in-place;
    - the full fine-grid image is inverse FFT'd, then downsampled in real space.

    ``psf_mode="kvalue"`` samples GalSim's analytic Fourier PSF through the C++
    backend, avoiding a real-space PSF draw and an extra FFT.
    """
    n = gal_prof.shape[0]
    real_dtype, complex_dtype = dtypes

    spectrum = _rfft_centered_image(
        xp=xp,
        image=gal_prof,
        real_dtype=real_dtype,
        complex_dtype=complex_dtype,
        center_index=center_index,
    )

    if integration_compensation_mode is not None and integration_order > 1:
        start = time.time() if profile else None
        _compensate_block_integration(
            xp=xp,
            spectrum=spectrum,
            n=n,
            fine_scale=scale,
            mode=integration_compensation_mode,
            integration_order=integration_order,
        )
        if profile:
            sync_if_gpu(xp)
            end = time.time()
            print(f"Block integration compensation took {end - start:.3f} seconds")

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
            center_index=center_index,
        )

        spectrum *= psf_spectrum
        del psf_image, psf_spectrum

    elif psf_obj is not None and psf_mode == "kvalue":
        start = time.time() if profile else None
        # The compiled backend returns the analytic PSF on BATSim's rFFT grid,
        # matching the spectrum convention used for the sampled galaxy.
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
        raise ValueError(f"Invalid psf_mode '{psf_mode}'; must be 'real' or 'kvalue'.")

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

    # Preserve the requested input flux after PSF/pixel convolution and any
    # integration compensation.
    if target_flux is not None:
        image *= target_flux / image.sum()

    coarse_image = _extract_centered_coarse_image(
        xp=xp,
        image=image,
        downsample_ratio=downsample_ratio,
        center_index=center_index,
    )
    del image

    return coarse_image
