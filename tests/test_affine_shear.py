# This function test the implementation of affine transform and compare it with
# Galsim. This test does not include PSF convolution.
import galsim
import numpy as np

import batsim


def test_affine_galsim(gamma1=0.2, gamma2=0.0, kappa=0.0):
    ## create a galaxy with raidall dependent shear
    flux = 40
    scale = 0.2
    nn = 64
    hlr = 1.4

    # create galaxy to be sampled by shear stamp objects
    sersic_gal = galsim.Sersic(n=1.0, half_light_radius=hlr, flux=flux, trunc=0)

    # reduced shear and lensing magnification
    g1 = gamma1 / (1 - kappa)
    g2 = gamma2 / (1 - kappa)
    mu = 1 / ((1 - kappa) ** 2 - gamma1**2 - gamma2**2)

    # # apply lensing shear to galaxy
    lens = batsim.LensTransform(gamma1=gamma1, gamma2=gamma2, kappa=kappa)

    # get galaxy array from stamp object
    gal_array = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=lens,
    )

    # apply the distortion with Galsim
    # note that we need to use lens instead of shear which is the nature
    # lensing shear that conserve surface density.
    # the galsim.shear function is flux conservative
    # (Based on the report from Gary Yang)
    gal_galsim = (
        sersic_gal.lens(g1=g1, g2=g2, mu=mu)
        .shift(0.5 * scale, 0.5 * scale)
        .drawImage(nx=nn, ny=nn, scale=scale, method="no_pixel")
    )
    # use no_pixel to draw image's flux but do not use pixel top-hat response
    gal_galsim = gal_galsim.array
    np.testing.assert_array_almost_equal(gal_array, gal_galsim)
    return


if __name__ == "__main__":
    test_affine_galsim()
