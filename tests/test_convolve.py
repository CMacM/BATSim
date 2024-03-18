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


def test_base():
    ## create a galaxy with raidall dependent shear
    image = (
        gal.shift(0.5 * scale, 0.5 * scale)
        .drawImage(nx=nn, ny=nn, scale=scale, method="no_pixel")
        .array
    )
    image_conv = (
        gal_conv.shift(scale * 0.5, scale * 0.5)
        .drawImage(nx=nn, ny=nn, scale=scale, method="no_pixel")
        .array
    )
    image_conv2 = batsim._gsinterface.convolvePsf(
        scale,
        psf._sbp,
        image,
        downsample_ratio=1,
        ngrid=nn,
    )

    np.testing.assert_array_almost_equal(image_conv2, image_conv, decimal=5)
    return


def test_downsample():
    ## create a galaxy with raidall dependent shear
    image = (
        gal.shift(0.5 * scale, 0.5 * scale)
        .drawImage(nx=nn, ny=nn, scale=scale, method="no_pixel")
        .array
    )
    image_conv = (
        gal_conv.shift(scale, scale)
        .drawImage(nx=nn, ny=nn, scale=scale * 2, method="no_pixel")
        .array
    )
    image_conv2 = batsim._gsinterface.convolvePsf(
        scale,
        psf._sbp,
        image,
        downsample_ratio=2,
        ngrid=nn,
    )

    np.testing.assert_array_almost_equal(image_conv2, image_conv, decimal=5)
    return


def test_truncate():
    ## create a galaxy with raidall dependent shear
    image = (
        gal.shift(0.5 * scale, 0.5 * scale)
        .drawImage(nx=nn, ny=nn, scale=scale, method="no_pixel")
        .array
    )
    image_conv = (
        gal_conv.shift(scale, scale)
        .drawImage(nx=nn // 4, ny=nn // 4, scale=scale * 2, method="no_pixel")
        .array
    )
    image_conv2 = batsim._gsinterface.convolvePsf(
        scale,
        psf._sbp,
        image,
        downsample_ratio=2,
        ngrid=int(nn / 4),
    )

    np.testing.assert_array_almost_equal(image_conv2, image_conv, decimal=5)
    return


if __name__ == "__main__":
    test_base()
    test_downsample()
    test_truncate()