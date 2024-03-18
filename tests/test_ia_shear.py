# This function test the implementation of the IA transform with power 0
# (constant shear) and compares it with BATSim's affine transfrom
import os

import galsim
import numpy as np

import batsim


def test_ia_shear():
    ## create a galaxy with raidall dependent shear
    flux = 40
    scale = 0.2
    nn = 64
    hlr = 1.4

    # create galaxy to be sampled by shear stamp objects
    sersic_gal = galsim.Sersic(n=1.0, half_light_radius=hlr, flux=flux, trunc=0)

    # apply affine shear with BATSim
    Lens = batsim.LensTransform(gamma1=0.2, gamma2=0, kappa=0)

    # sample galaxy object onto stamp
    Lens_gal = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=Lens,
    )

    # determine correct IA ampltidue to match g1 shear
    A_IA = galsim.Shear(g1=0.2).e1

    # apply IA shear to galaxy
    IA = batsim.IaTransform(A=A_IA, beta=0, phi=0, scale=scale, hlr=hlr)

    # get galaxy array from stamp object
    IA_gal = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=sersic_gal,
        transform_obj=IA,
    )

    np.testing.assert_array_almost_equal(IA_gal, Lens_gal)

    # Now we test that the power law gives us the expected value at
    # the half light radius

    # initialise new transform with non-zero power law
    IA_pow = batsim.IaTransform(A=A_IA, beta=0.8, phi=0, scale=scale, hlr=hlr)

    # set coords for hlr
    x = np.array([hlr])
    y = np.array([0])

    # get the expected g1 at the half light radius
    test_g1, _ = IA_pow.get_g1g2(x, y)
    # shear at hlr should = input amplitude
    np.testing.assert_almost_equal(test_g1, 0.2)

    # now do same but for coord outside hlr
    # set coords for hlr
    x = np.array([0])
    y = np.array([1.2 * hlr])

    # get the expected g1 at the half light radius
    test_g1, _ = IA_pow.get_g1g2(x, y)
    # outside shear should be greater than hlr shear
    np.testing.assert_array_less(0.2, test_g1)

    # now do same but for coord inside hlr
    # set coords for hlr
    x = np.array([0])
    y = np.array([0.8 * hlr])

    # get the expected g1 at the half light radius
    test_g1, _ = IA_pow.get_g1g2(x, y)
    # inside shear should be less than hlr shear
    np.testing.assert_array_less(test_g1, 0.2)
    return


if __name__ == "__main__":
    test_ia_shear()
