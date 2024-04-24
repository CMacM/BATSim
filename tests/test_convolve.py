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
        psf_obj=psf
    )

    # apply the distortion with Galsim note that we need to use lens instead of
    # shear which is the nature lensing shear that conserve surface density.
    # the galsim.shear function is flux conservative (Based on the report from
    # Gary Yang)
    gal_galsim = (
        conv_gal.shift(0.5 * scale, 0.5 * scale).drawImage(
            nx=nn, ny=nn, scale=scale, method="auto"
        )
    ).array
    np.testing.assert_array_almost_equal(gal_array, gal_galsim, decimal=4)
    return

def test_draw_methods():

    # set up galaxy object
    galaxy = galsim.Sersic(n=3, half_light_radius=1.0, flux=1.0).shear(e1=0.1, e2=0.03)

    # drawing parameters
    scale = 0.2
    nn=64

    # set up psf object
    seeing = 0.6
    psf = galsim.Moffat(beta=3.5, fwhm=seeing, trunc=4*seeing)

    # Create galsim image with no pixel convolution
    galsim_conv = galsim.Convolve([galaxy, psf])
    galsim_image_np = galsim_conv.shift(0.5*scale, 0.5*scale).drawImage(nx=nn, ny=nn, scale=scale, method='no_pixel').array

    # Create galsim image with pixel convolution
    galsim_image_auto = galsim_conv.shift(0.5*scale, 0.5*scale).drawImage(nx=nn, ny=nn, scale=scale, method='auto').array

    # Create batsim image with no pixel convolution
    batsim_image_np = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=galaxy,
        transform_obj=None,
        psf_obj=psf,
        draw_method='no_pixel'
    )

    # Create batsim image with pixel convolution
    # Pixel profile is defined in function using provided pix_scale
    batsim_image_auto = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=galaxy,
        transform_obj=None,
        psf_obj=psf,
        draw_method='auto'
    )

    # Test equalities
    np.testing.assert_array_almost_equal(galsim_image_np, batsim_image_np, decimal=3)
    np.testing.assert_array_almost_equal(galsim_image_auto, batsim_image_auto, decimal=3)

    return

def test_no_psf():

    galsim_image = gal.shift(0.5 * scale, 0.5 * scale).drawImage(nx=nn, ny=nn, scale=scale, method="auto")

    batsim_image = batsim.simulate_galaxy(
        ngrid=nn,
        pix_scale=scale,
        gal_obj=gal,
        transform_obj=None,
        psf_obj=None,
        draw_method="auto",
    )

    np.testing.assert_array_almost_equal(galsim_image.array, batsim_image, decimal=5)

    return


if __name__ == "__main__":
    test_base()
    test_downsample()
    test_truncate()
    test_convolved_lensed()
    test_draw_methods()
    test_no_psf()
