# This function test the implementation of affine transform and compare it with
# Galsim. This test does not include PSF convolution.
import galsim
import numpy as np

import batsim

nn = 128
scale = 0.2
gal = (
    galsim.Sersic(
        n=0.5,
        half_light_radius=1.5,
        flux=20,
    )
    .shear(g1=0.24, g2=-0.13)
    .shift(2 * scale, 3 * scale)
)
psf = galsim.Gaussian(fwhm=0.3, flux=1.0)
gal_conv = galsim.Convolve([gal, psf])


def _center_crop(image, ngrid):
    y0 = (image.shape[0] - ngrid) // 2
    x0 = (image.shape[1] - ngrid) // 2
    return image[y0 : y0 + ngrid, x0 : x0 + ngrid]


def _moments(image):
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

    return np.array([flux, x0, y0, (qxx - qyy) / trace, 2.0 * qxy / trace, trace])


def test_psf_no_pixel_default_render_matches_galsim_shape():
    image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=gal,
        psf_obj=psf,
        draw_method="no_pixel",
    )

    reference = gal_conv.drawImage(
        nx=nn,
        ny=nn,
        scale=scale,
        method="no_pixel",
    ).array

    np.testing.assert_allclose(image.sum(), gal.flux, rtol=2.0e-4)
    np.testing.assert_allclose(_moments(image)[3:5], _moments(reference)[3:5], atol=5.0e-3)


def test_output_ngrid_is_center_crop_of_default_render():
    full_image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=gal,
        psf_obj=psf,
        draw_method="no_pixel",
    )
    cropped_image = batsim.simulate_galaxy(
        ngrid=nn // 4,
        pix_scale=scale,
        gal_obj=gal,
        psf_obj=psf,
        draw_method="no_pixel",
    )

    np.testing.assert_allclose(cropped_image, _center_crop(full_image, nn // 4))


def test_convolved_lensed(gamma1=0.2, gamma2=0.0, kappa=0.0):
    # reduced shear and lensing magnification
    g1 = gamma1 / (1 - kappa)
    g2 = gamma2 / (1 - kappa)
    mu = 1 / ((1 - kappa) ** 2 - gamma1**2 - gamma2**2)

    lensed_gal = gal.lens(g1=g1, g2=g2, mu=mu)
    conv_gal = galsim.Convolve([lensed_gal, psf])
    # # apply lensing shear to galaxy
    lens = batsim.LensTransform(gamma1=gamma1, gamma2=gamma2, kappa=kappa)

    # get galaxy array from stamp object
    gal_array = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=gal,
        transform_obj=lens,
        psf_obj=psf,
    )

    reference = conv_gal.drawImage(nx=nn, ny=nn, scale=scale, method="auto").array

    np.testing.assert_allclose(gal_array.sum(), gal.flux, rtol=2.0e-4)
    np.testing.assert_allclose(_moments(gal_array)[3:5], _moments(reference)[3:5], atol=2.0e-2)


def test_draw_methods():

    # set up galaxy object
    galaxy = galsim.Sersic(n=3, half_light_radius=1.0, flux=1.0).shear(e1=0.1, e2=0.03)

    # drawing parameters
    scale = 0.2
    nn = 64

    # set up psf object
    seeing = 0.6
    psf = galsim.Moffat(beta=3.5, fwhm=seeing, trunc=4 * seeing)

    # Create batsim image with no pixel convolution
    batsim_image_np = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=galaxy,
        transform_obj=None,
        psf_obj=psf,
        draw_method="no_pixel",
    )

    # Create batsim image with pixel convolution
    # Pixel profile is defined in function using provided pix_scale
    batsim_image_auto = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=galaxy,
        transform_obj=None,
        psf_obj=psf,
        draw_method="auto",
    )

    np.testing.assert_allclose(batsim_image_np.sum(), galaxy.flux, rtol=5.0e-2)
    np.testing.assert_allclose(batsim_image_auto.sum(), galaxy.flux, rtol=5.0e-2)
    assert not np.allclose(batsim_image_np, batsim_image_auto)


def test_no_psf():

    galsim_image = gal.drawImage(nx=nn, ny=nn, scale=scale, method="auto").array

    batsim_image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=gal,
        transform_obj=None,
        psf_obj=None,
        draw_method="auto",
    )

    np.testing.assert_allclose(batsim_image.sum(), gal.flux, rtol=2.0e-4)
    np.testing.assert_allclose(_moments(batsim_image)[3:5], _moments(galsim_image)[3:5], atol=1.0e-2)


if __name__ == "__main__":
    test_psf_no_pixel_default_render_matches_galsim_shape()
    test_output_ngrid_is_center_crop_of_default_render()
    test_convolved_lensed()
    test_draw_methods()
    test_no_psf()
