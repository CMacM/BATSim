import galsim
import numpy as np

from . import _gsinterface
from .stamp import Stamp


def simulate_galaxy(
    ngrid,
    pix_scale,
    gal_obj,
    transform_obj=None,
    psf_obj=None,
    truncation_ratio=20.0,
):
    """The function samples the surface density field of a galaxy at the grids
    This function only conduct sampling; PSF and pixel response are not
    included.

    Args:
    ngrid (int):        number of grids
    pix_scale (float):  pixel scale
    gal_obj (galsim):   Galsim galaxy object to sample on the grids
    transform_obj :     Coordinate transform object
    psf_obj (galsim):   Galsim PSF object to smear the image

    Returns:
    outcome (ndarray):  2D galaxy image on the grids
    """
    # Get scale used to make simulation
    if psf_obj is None:
        scale = pix_scale
        psf_obj = galsim.Gaussian(fwhm=1e-10, flux=1.0)
        pad_arcsec = 0.0
        downsample_ratio = 1
    else:
        scale = min(gal_obj.nyquist_scale, psf_obj.nyquist_scale / 4.0)
        scale = min(scale, pix_scale * 0.5)
        pad_arcsec = psf_obj.calculateFWHM()
        # set the maximum oversample ratio to 32
        downsample_ratio = int(min(np.floor(pix_scale / scale), 32))
        scale = pix_scale / downsample_ratio

    # Get number of grids to generate simulation
    npad = int(pad_arcsec / scale + 0.5) * 2
    nn = int(gal_obj.calculateMomentRadius() / scale * truncation_ratio) + npad * 2
    nn = max(int(2 ** np.ceil(np.log2(nn))), ngrid * downsample_ratio)

    # Initialize and Distort Coordinates
    stamp = Stamp(nn=nn, scale=scale)
    if transform_obj is not None:
        # Distort galaxy
        gal_coords = transform_obj.transform(stamp.coords)
    else:
        gal_coords = stamp.coords
    # Record flux
    gal_prof = _gsinterface.getFluxVec(
        scale=scale,
        gsobj=gal_obj._sbp,
        xy_coords=gal_coords,
    )

    # Convolution in Fourier space
    gal_prof = _gsinterface.convolvePsf(
        scale=scale,
        gsobj=psf_obj._sbp,
        gal_prof=gal_prof,
        downsample_ratio=downsample_ratio,
        ngrid=ngrid,
    )
    return gal_prof
