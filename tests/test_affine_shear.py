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
        draw_method="no_pixel"
    )

    # apply the distortion with Galsim note that we need to use lens instead of
    # shear which is the nature lensing shear that conserve surface density.
    # the galsim.shear function is flux conservative (Based on the report from
    # Gary Yang)
    gal_galsim = (
        lensed_gal.shift(0.5 * scale, 0.5 * scale).drawImage(
            nx=nn, ny=nn, scale=scale, method="no_pixel"
        )
    ).array
    np.testing.assert_array_almost_equal(gal_array, gal_galsim)
    return


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
        draw_method="auto"
    )

    # apply the distortion with Galsim note that we need to use lens instead of
    # shear which is the nature lensing shear that conserve surface density.
    # the galsim.shear function is flux conservative (Based on the report from
    # Gary Yang)
    gal_galsim = (
        smeared_gal.shift(0.5 * scale, 0.5 * scale).drawImage(
            nx=nn, ny=nn, scale=scale, method="auto"
        )
    ).array
    np.testing.assert_array_almost_equal(gal_array, gal_galsim, decimal=4)
    return


if __name__ == "__main__":
    test_affine()
    test_affine_psf()
