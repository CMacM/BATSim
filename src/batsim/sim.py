import galsim
import numpy as np
from scipy.fft import fftshift, ifftshift, irfft2, rfft2

from . import _gsinterface
from .stamp import Stamp


def simulate_galaxy(
    nn,
    scale,
    gal_obj,
    transform_obj=None,
    psf_obj=None,
    pixel_response=False,
):
    """Sample the surface density field of a galaxy at the grids
    This function only conduct sampling; PSF and pixel response are
    not included.

    Args:
    nn (int):           number of grids
    scale (float):      pixel scale
    gal_obj (galsim):   Galsim galaxy object to sample on the grids
    transform_obj :     Coordinate transform object
    psf_obj (galsim):   Galsim PSF object to smear the image

    Returns:
    outcome (ndarray):  2D galaxy image on the grids
    """
    if psf_obj is not None:
        npad = int(psf_obj.calculateFWHM() / scale + 0.5) * 4
    else:
        npad = 0
    stamp = Stamp(nn=nn + npad * 2, scale=scale)
    if transform_obj is not None:
        # Distort galaxy
        gal_coords = transform_obj.transform(stamp.coords)
    else:
        gal_coords = stamp.coords
    gal_prof = _gsinterface.getFluxVec(gal_obj._sbp, gal_coords)

    if psf_obj is not None:
        # Convolution in Fourier space
        psf_coords = stamp.coords
        gal_kprof = _gsinterface.mulFourier(
            scale,
            psf_obj._sbp,
            rfft2(gal_prof),
        )
        gal_prof = irfft2(gal_kprof, stamp.shape)[npad:-npad, npad:-npad]

    # normailize to flux
    gal_prof = gal_prof * stamp.pixel_area
    return gal_prof
