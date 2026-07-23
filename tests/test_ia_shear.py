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


def _reduced_shear_inverse_matrix(g1, g2):
    inv_norm = 1.0 / np.sqrt(1.0 - g1 * g1 - g2 * g2)
    return inv_norm * np.array([[1.0 - g1, -g2], [-g2, 1.0 + g1]])


def test_beta_zero_ia_transform_uses_unit_determinant_reduced_shear_matrix():
    scale = 0.2
    hlr = 1.4
    amplitude = 0.2
    phi = 0.3

    gamma1 = amplitude * np.cos(2.0 * phi)
    gamma2 = amplitude * np.sin(2.0 * phi)

    ia_transform = batsim.IaTransform(A=amplitude, beta=0.0, phi=phi, scale=scale, hlr=hlr)

    x = np.linspace(-2.0 * hlr, 2.0 * hlr, 9)
    y = np.linspace(-1.5 * hlr, 1.5 * hlr, 7)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    coords = np.stack([xx.ravel(), yy.ravel()], axis=0)

    expected_matrix = _reduced_shear_inverse_matrix(gamma1, gamma2)

    np.testing.assert_allclose(np.linalg.det(expected_matrix), 1.0, rtol=1.0e-14)
    np.testing.assert_allclose(
        ia_transform.transform(coords),
        expected_matrix @ coords,
        rtol=1.0e-14,
        atol=1.0e-14,
    )


def test_beta_zero_ia_render_matches_galsim_reduced_shear_moments():
    flux = 40
    scale = 0.15
    nn = 96
    hlr = 1.1
    g1 = 0.18
    g2 = 0.07
    amplitude = np.hypot(g1, g2)
    phi = 0.5 * np.arctan2(g2, g1)

    gal = galsim.Gaussian(half_light_radius=hlr, flux=flux)
    ia_transform = batsim.IaTransform(A=amplitude, beta=0.0, phi=phi, scale=scale, hlr=hlr)

    reference = (
        gal.shear(g1=g1, g2=g2)
        .drawImage(
            nx=nn,
            ny=nn,
            scale=scale,
            method="no_pixel",
            use_true_center=False,
        )
        .array
    )

    ia_image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=gal,
        transform_obj=ia_transform,
        draw_method="no_pixel",
        force_input_flux=False,
        use_true_center=False,
    )

    reference_moments = _moments(reference, scale)
    ia_moments = _moments(ia_image, scale)

    np.testing.assert_allclose(ia_moments[0], reference_moments[0], rtol=5.0e-6)
    np.testing.assert_allclose(
        ia_moments[3:6],
        reference_moments[3:6],
        rtol=2.0e-3,
        atol=2.0e-4,
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
    test_beta_zero_ia_transform_uses_unit_determinant_reduced_shear_matrix()
    test_beta_zero_ia_render_matches_galsim_reduced_shear_moments()
    test_ia_power_law_amplitude_is_hlr_normalized()
