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
        npad = int(psf_obj.calculateFWHM() / scale + 0.5) * 2
    else:
        npad = 0

    nn_new = max(int(2 ** np.ceil(np.log2(nn + npad * 2))), 64)
    off = int(nn_new - nn) // 2
    nn = nn_new

    stamp = Stamp(nn=nn, scale=scale)
    if transform_obj is not None:
        # Distort galaxy
        gal_coords = transform_obj.transform(stamp.coords)
    else:
        gal_coords = stamp.coords
    gal_prof = _gsinterface.getFluxVec(gal_obj._sbp, gal_coords) * stamp.pixel_area

    if psf_obj is not None:
        # Convolution in Fourier space
        gal_prof = _gsinterface.convolvePsf(
            scale,
            psf_obj._sbp,
            gal_prof,
        )

    if off > 0:
        gal_prof = gal_prof[off:-off, off:-off]
    return gal_prof
