# This function test the implementation of affine transform and compare it with
# Galsim
import os
import galsim
import numpy as np

from batsim.stamp import Stamp
from batsim.transforms import LensTransform


def test_affine_galsim():
    ## create a galaxy with raidall dependent shear
    flux = 40
    scale = 0.2
    nn = 64
    hlr = 1.4

    # create galaxy to be sampled by shear stamp objects
    sersic_gal = galsim.Sersic(n=0.5, half_light_radius=hlr, flux=flux, trunc=0)

    # lensing shear and kappa convergence
    gamma1 = 0.2
    gamma2 = 0.0
    kappa= 0.0

    # reduced shear and lensing magnification
    g1 = gamma1/(1-kappa)
    g2 = gamma2/(1-kappa)
    mu = 1/((1-kappa)**2 - gamma1**2 - gamma2**2)

    # # apply lensing shear to galaxy
    stamp = Stamp(nn=nn, scale=scale)
    lens = LensTransform(gamma1=g1, gamma2=g2, kappa=kappa, xref=-0.5*scale, yref=-0.5*scale)
    stamp.transform_grids(lens)

    # get galaxy array from stamp object
    gal_array = stamp.sample_galaxy(sersic_gal)

    # apply the distortion with Galsim
    # note that we need to use lens instead of shear which is the nature
    # lensing shear that conserve surface density.
    # the galsim.shear function is flux conservative
    # (Based on the report from Gary Yang)
    gal_galsim = sersic_gal.lens(
        g1=g1,
        g2=g2,
        mu = mu
    ).drawImage(nx=64, ny=64, scale=scale, method='sb')
    gal_galsim = gal_galsim.array
    np.testing.assert_array_almost_equal(gal_array, gal_galsim)
    return


if __name__ == "__main__":
    test_affine_galsim()
