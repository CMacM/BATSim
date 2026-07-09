# This function test the implementation of affine transform and compare it with
# Galsim. This test does not include PSF convolution.
import galsim
import numpy as np

import batsim

flux = 40
scale = 0.2
nn = 64
hlr = 1.3
# create galaxy to be sampled by shear stamp objects
sersic_gal = galsim.Sersic(n=1.0, half_light_radius=hlr, flux=flux, trunc=0)
psf = galsim.Moffat(beta=3.5, fwhm=0.85, flux=1.0)


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


def test_affine(gamma1=0.2, gamma2=0.0, kappa=0.0):
    # reduced shear and lensing magnification
    g1 = gamma1 / (1 - kappa)
    g2 = gamma2 / (1 - kappa)
    mu = 1 / ((1 - kappa) ** 2 - gamma1**2 - gamma2**2)

    lensed_gal = sersic_gal.lens(g1=g1, g2=g2, mu=mu)
    # # apply lensing shear to galaxy
    lens = batsim.LensTransform(gamma1=gamma1, gamma2=gamma2, kappa=kappa)

    # get galaxy array from stamp object
    gal_array = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=lens,
        draw_method="no_pixel",
    )

    # The current BATSim default normalises to the input flux, while GalSim's
    # lens call includes lensing magnification.  Compare shape moments rather
    # than exact pixel values.
    gal_galsim = lensed_gal.drawImage(
        nx=nn,
        ny=nn,
        scale=scale,
        method="no_pixel",
    ).array

    np.testing.assert_allclose(gal_array.sum(), flux, rtol=5.0e-3)
    np.testing.assert_allclose(
        _moments(gal_array)[3:5],
        _moments(gal_galsim)[3:5],
        atol=1.0e-2,
    )


def test_affine_psf(gamma1=0.2, gamma2=0.0, kappa=0.0):
    # reduced shear and lensing magnification
    g1 = gamma1 / (1 - kappa)
    g2 = gamma2 / (1 - kappa)
    mu = 1 / ((1 - kappa) ** 2 - gamma1**2 - gamma2**2)

    lensed_gal = sersic_gal.lens(g1=g1, g2=g2, mu=mu)
    smeared_gal = galsim.Convolve([lensed_gal, psf])
    # # apply lensing shear to galaxy
    lens = batsim.LensTransform(gamma1=gamma1, gamma2=gamma2, kappa=kappa)

    # get galaxy array from stamp object
    gal_array = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=lens,
        psf_obj=psf,
        draw_method="auto",
    )

    # Use GalSim as a shape reference; BATSim defaults preserve input flux.
    gal_galsim = smeared_gal.drawImage(
        nx=nn,
        ny=nn,
        scale=scale,
        method="auto",
    ).array

    np.testing.assert_allclose(gal_array.sum(), flux, rtol=5.0e-3)
    np.testing.assert_allclose(
        _moments(gal_array)[3:5],
        _moments(gal_galsim)[3:5],
        atol=1.0e-2,
    )


if __name__ == "__main__":
    test_affine()
    test_affine_psf()
