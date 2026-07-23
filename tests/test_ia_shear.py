import galsim
import numpy as np

import batsim


def _moments(image, scale):
    n = image.shape[0]
    ind = (np.arange(n) - 0.5 * (n - 1)) * scale
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
            x0,
            y0,
            (qxx - qyy) / trace,
            2.0 * qxy / trace,
            trace,
        ]
    )


def test_beta_zero_ia_transform_matches_lens_transform_coordinates():
    scale = 0.2
    hlr = 1.4
    amplitude = 0.2
    phi = 0.3

    gamma1 = amplitude * np.cos(2.0 * phi)
    gamma2 = amplitude * np.sin(2.0 * phi)

    ia_transform = batsim.IaTransform(A=amplitude, beta=0.0, phi=phi, scale=scale, hlr=hlr)
    lens_transform = batsim.LensTransform(gamma1=gamma1, gamma2=gamma2, kappa=0.0)

    x = np.linspace(-2.0 * hlr, 2.0 * hlr, 9)
    y = np.linspace(-1.5 * hlr, 1.5 * hlr, 7)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    coords = np.stack([xx.ravel(), yy.ravel()], axis=0)

    np.testing.assert_allclose(
        ia_transform.transform(coords),
        lens_transform.transform(coords),
        rtol=1.0e-14,
        atol=1.0e-14,
    )


def test_beta_zero_ia_render_matches_lens_transform_moments():
    flux = 40
    scale = 0.2
    nn = 64
    hlr = 1.4
    g1 = 0.2

    sersic_gal = galsim.Sersic(n=1.0, half_light_radius=hlr, flux=flux, trunc=0)
    lens_transform = batsim.LensTransform(gamma1=g1, gamma2=0.0, kappa=0.0)
    ia_transform = batsim.IaTransform(A=g1, beta=0.0, phi=0.0, scale=scale, hlr=hlr)

    lens_image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=lens_transform,
    )

    ia_image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=ia_transform,
    )

    np.testing.assert_allclose(ia_image.sum(), lens_image.sum(), rtol=1.0e-12)
    np.testing.assert_allclose(
        _moments(ia_image, scale)[3:5],
        _moments(lens_image, scale)[3:5],
        atol=1.0e-12,
    )


def test_ia_power_law_amplitude_is_hlr_normalized():
    scale = 0.2
    hlr = 1.4
    amplitude = 0.2
    beta = 0.8

    ia_transform = batsim.IaTransform(A=amplitude, beta=beta, phi=0.0, scale=scale, hlr=hlr)

    test_g1, test_g2 = ia_transform.get_g1g2(np.array([hlr]), np.array([0.0]))
    np.testing.assert_allclose(test_g1, amplitude)
    np.testing.assert_allclose(test_g2, 0.0)

    inner_g1, _ = ia_transform.get_g1g2(np.array([0.0]), np.array([0.8 * hlr]))
    outer_g1, _ = ia_transform.get_g1g2(np.array([0.0]), np.array([1.2 * hlr]))

    assert inner_g1 < amplitude
    assert outer_g1 > amplitude


if __name__ == "__main__":
    test_beta_zero_ia_transform_matches_lens_transform_coordinates()
    test_beta_zero_ia_render_matches_lens_transform_moments()
    test_ia_power_law_amplitude_is_hlr_normalized()
