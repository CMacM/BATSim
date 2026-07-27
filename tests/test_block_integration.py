import importlib

import galsim
import numpy as np
import pytest

import batsim
from batsim.stamp import Stamp

sim = importlib.import_module("batsim.sim")


def _quadratic_flux(coords, fine_scale, real_dtype):
    x = coords[0]
    y = coords[1]
    values = x * x + 2.0 * x * y + 3.0 * y * y
    values = values * fine_scale**2

    # Match the C++ sampler shape convention.  The Python code should still
    # recover the offset-major ordering when it reshapes this square array.
    n = int(np.sqrt(coords.shape[1]))
    return values.astype(real_dtype, copy=False).reshape(n, n)


class MatrixTransform:
    def __init__(self, matrix):
        self.matrix = np.asarray(matrix)

    def transform(self, coords):
        return self.matrix @ coords


class BandwidthProbe:
    flux = 1.0
    maxk = 8.0 * np.pi

    def getGoodImageSize(self, scale):
        return 8


def _moments(image, scale):
    n = image.shape[0]
    ind = (np.arange(n) - n // 2) * scale
    yy, xx = np.meshgrid(ind, ind, indexing="ij")

    flux = image.sum()
    x0 = (image * xx).sum() / flux
    y0 = (image * yy).sum() / flux

    dx = xx - x0
    dy = yy - y0
    qxx = (image * dx * dx).sum() / flux
    qyy = (image * dy * dy).sum() / flux
    qxy = (image * dx * dy).sum() / flux

    trace = qxx + qyy
    return np.array(
        [
            flux,
            (qxx - qyy) / trace,
            2.0 * qxy / trace,
            trace,
        ]
    )


def _peak_percent_flux_residual(image, reference):
    return 100.0 * np.max(np.abs(image - reference)) / np.max(reference)


def test_integration_offsets_are_fine_pixel_offsets():
    fine_scale = 0.2
    offsets, weights = sim._integration_offsets(
        xp=np,
        fine_scale=fine_scale,
        integration_order=4,
        dtype=np.float64,
    )

    assert offsets.shape == (16, 2)
    assert weights.shape == (16,)
    assert np.max(np.abs(offsets)) < 0.5 * fine_scale

    np.testing.assert_allclose(weights.sum(), 1.0)
    np.testing.assert_allclose(weights @ offsets[:, 0], 0.0, atol=1.0e-15)
    np.testing.assert_allclose(weights @ offsets[:, 1], 0.0, atol=1.0e-15)
    np.testing.assert_allclose(weights @ offsets[:, 0] ** 2, fine_scale**2 / 12.0)
    np.testing.assert_allclose(weights @ offsets[:, 1] ** 2, fine_scale**2 / 12.0)


def test_integration_compensation_modes_are_explicit():
    assert sim._resolve_integration_compensation("quadrature") == "quadrature"
    assert sim._resolve_integration_compensation("exact_sinc") == "exact_sinc"
    assert sim._resolve_integration_compensation(None) is None

    for mode in (True, False, "sinc", "exact", "gauss-legendre"):
        with pytest.raises(ValueError):
            sim._resolve_integration_compensation(mode)


def test_no_psf_no_pixel_uses_fine_grid_without_block_integration_or_fft(monkeypatch):
    records = {}

    def sample_without_block_integration(
        gal_obj,
        stamp,
        transform_obj,
        fine_scale,
        real_dtype,
        xp,
        integration_order=2,
        profile=False,
    ):
        records["integration_order"] = integration_order
        records["downsample_ratio"] = stamp.downsample_ratio
        records["fine_scale"] = fine_scale
        return xp.ones(stamp.shape, dtype=real_dtype) * fine_scale**2

    def fail_if_fft_is_used(*args, **kwargs):
        raise AssertionError("FFT path should not be used without PSF or pixel response.")

    monkeypatch.setattr(sim, "_sample_galaxy_profile_on_stamp", sample_without_block_integration)
    monkeypatch.setattr(sim, "_convolve_psf_fft", fail_if_fft_is_used)

    image = sim._simulate_galaxy_impl(
        gal_obj=BandwidthProbe(),
        scale=1.0,
        ngrid=8,
        psf_obj=None,
        draw_method="no_pixel",
        safety=1.0,
        max_supersample=16,
        min_supersample=1,
        integration_order=4,
        max_fine_grid=None,
        pad=0,
        backend="np",
        force_input_flux=False,
        compensate_integration="quadrature",
    )

    assert records["integration_order"] == 1
    assert records["downsample_ratio"] == 8
    np.testing.assert_allclose(records["fine_scale"], 0.125)
    np.testing.assert_allclose(image, np.ones((8, 8)))


def test_block_integration_preserves_offset_order_and_flux_scale(monkeypatch):
    fine_scale = 0.2
    stamp = Stamp(nn=5, scale=fine_scale, backend=np, dtype=np.float64)

    monkeypatch.setattr(
        sim,
        "_sample_galaxy_profile",
        lambda gal_obj, coords, fine_scale, real_dtype: _quadratic_flux(
            coords,
            fine_scale,
            real_dtype,
        ),
    )

    image = sim._sample_galaxy_profile_on_stamp(
        gal_obj=object(),
        stamp=stamp,
        transform_obj=None,
        fine_scale=fine_scale,
        real_dtype=np.float64,
        xp=np,
        integration_order=4,
    )

    x = stamp.coords[0].reshape(stamp.nn, stamp.nn)
    y = stamp.coords[1].reshape(stamp.nn, stamp.nn)
    expected_mean = x * x + 2.0 * x * y + 3.0 * y * y + fine_scale**2 / 3.0
    expected = expected_mean * fine_scale**2

    np.testing.assert_allclose(image, expected, rtol=1.0e-13, atol=1.0e-15)


def test_transform_is_applied_after_offsets(monkeypatch):
    fine_scale = 0.2
    stamp = Stamp(nn=4, scale=fine_scale, backend=np, dtype=np.float64)
    transform = MatrixTransform([[1.2, 0.3], [-0.4, 0.8]])

    monkeypatch.setattr(
        sim,
        "_sample_galaxy_profile",
        lambda gal_obj, coords, fine_scale, real_dtype: _quadratic_flux(
            coords,
            fine_scale,
            real_dtype,
        ),
    )

    image = sim._sample_galaxy_profile_on_stamp(
        gal_obj=object(),
        stamp=stamp,
        transform_obj=transform,
        fine_scale=fine_scale,
        real_dtype=np.float64,
        xp=np,
        integration_order=3,
    )

    offsets, weights = sim._integration_offsets(
        xp=np,
        fine_scale=fine_scale,
        integration_order=3,
        dtype=np.float64,
    )

    expected = np.zeros(stamp.coords.shape[1], dtype=np.float64)
    for offset, weight in zip(offsets, weights):
        coords = transform.transform(stamp.coords + offset[:, None])
        x = coords[0]
        y = coords[1]
        expected += weight * (x * x + 2.0 * x * y + 3.0 * y * y)

    expected = (expected * fine_scale**2).reshape(stamp.nn, stamp.nn)

    np.testing.assert_allclose(image, expected, rtol=1.0e-13, atol=1.0e-15)


def test_psf_auto_affine_shape_matches_galsim_percent_level():
    scale = 0.2
    ngrid = 64
    gamma1 = 0.12
    gamma2 = 0.04
    kappa = 0.0

    gal = galsim.Gaussian(sigma=0.55, flux=1.0)
    psf = galsim.Gaussian(fwhm=0.45, flux=1.0)

    reduced_g1 = gamma1 / (1.0 - kappa)
    reduced_g2 = gamma2 / (1.0 - kappa)
    mu = 1.0 / ((1.0 - kappa) ** 2 - gamma1**2 - gamma2**2)

    reference = (
        galsim.Convolve(
            [
                gal.lens(g1=reduced_g1, g2=reduced_g2, mu=mu),
                psf,
            ]
        )
        .drawImage(nx=ngrid, ny=ngrid, scale=scale, method="auto")
        .array
    )

    rendered = batsim.simulate_galaxy(
        gal_obj=gal,
        pix_scale=scale,
        ngrid=ngrid,
        transform_obj=batsim.LensTransform(gamma1, gamma2, kappa),
        psf_obj=psf,
        draw_method="auto",
        integration_order=2,
        backend="np",
        force_input_flux=False,
        compensate_integration="quadrature",
    )

    reference_moments = _moments(reference, scale)
    rendered_moments = _moments(rendered, scale)

    np.testing.assert_allclose(rendered_moments[1:3], reference_moments[1:3], atol=1.0e-3)


def test_block_integration_compensation_modes_reduce_peak_flux_residual():
    scale = 0.2
    ngrid = 64
    gamma1 = 0.12
    gamma2 = 0.04
    kappa = 0.0

    gal = galsim.Gaussian(sigma=0.35, flux=1.0)

    reduced_g1 = gamma1 / (1.0 - kappa)
    reduced_g2 = gamma2 / (1.0 - kappa)
    mu = 1.0 / ((1.0 - kappa) ** 2 - gamma1**2 - gamma2**2)

    reference = (
        gal.lens(g1=reduced_g1, g2=reduced_g2, mu=mu)
        .drawImage(
            nx=ngrid,
            ny=ngrid,
            scale=scale,
            method="auto",
            use_true_center=False,
        )
        .array
    )

    common_kwargs = dict(
        gal_obj=gal,
        pix_scale=scale,
        ngrid=ngrid,
        transform_obj=batsim.LensTransform(gamma1, gamma2, kappa),
        psf_obj=None,
        draw_method="auto",
        integration_order=2,
        min_supersample=1,
        max_supersample=2,
        backend="np",
        force_input_flux=False,
        use_true_center=False,
    )

    residuals = {
        "default": _peak_percent_flux_residual(
            batsim.simulate_galaxy(**common_kwargs),
            reference,
        ),
        "none": _peak_percent_flux_residual(
            batsim.simulate_galaxy(
                **common_kwargs,
                compensate_integration=None,
            ),
            reference,
        ),
        "exact_sinc": _peak_percent_flux_residual(
            batsim.simulate_galaxy(
                **common_kwargs,
                compensate_integration="exact_sinc",
            ),
            reference,
        ),
        "quadrature": _peak_percent_flux_residual(
            batsim.simulate_galaxy(
                **common_kwargs,
                compensate_integration="quadrature",
            ),
            reference,
        ),
    }

    assert residuals["none"] > 1.0
    assert residuals["exact_sinc"] < 0.05
    assert residuals["quadrature"] < 0.01
    assert residuals["quadrature"] < 0.25 * residuals["exact_sinc"]
    np.testing.assert_allclose(residuals["default"], residuals["quadrature"])
