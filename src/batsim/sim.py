import galsim
import numpy as np
import os
import multiprocessing as mp

from . import _gsinterface
from .stamp import Stamp


def simulate_galaxy(
    ngrid,
    pix_scale,
    gal_obj,
    transform_obj=None,
    psf_obj=None,
    truncate_ratio=1.0,
    maximum_num_grids=4096,
    draw_method="auto",
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
    truncate_ratio (float):
                        truncate at truncate_ratio times good_image_size
    maximum_num_grids (int):
                        maximum number of grids for simulation in real space
    draw_method (str):  method to draw the galaxy image, "auto" will convolve with
                        pixel response, "no_pixel" is as it implies

    Returns:
    outcome (ndarray):  2D galaxy image on the grids
    """

    # Initialize variables based on PSF presence
    if psf_obj is None and draw_method == "no_pixel":
        # In this case we just get the fluxes for the requested stamp size
        scale = pix_scale
        nn = int(ngrid)
    else:
        # Compute the effective scale for simulation
        if psf_obj is None:
            scale = pix_scale / 4.0
            pad_arcsec = 0.0
            downsample_ratio = 1
        else:
            scale = min(gal_obj.nyquist_scale, psf_obj.nyquist_scale / 4.0, pix_scale / 4.0)
            pad_arcsec = psf_obj.calculateMomentRadius(size=32, scale=pix_scale / 2.0)
            downsample_ratio = min(int(2 ** np.ceil(np.log2(pix_scale / scale))), 128)
        
        scale = pix_scale / downsample_ratio

        # Calculate the number of grids considering padding and truncation
        npad = int(pad_arcsec / scale + 0.5) * 4
        nn = npad * 2 + min(gal_obj.getGoodImageSize(pixel_scale=scale) 
                            * truncate_ratio, ngrid * downsample_ratio
                            )
        nn = min(int(2 ** np.ceil(np.log2(nn))), maximum_num_grids)

    # Initialize and Distort Coordinates
    stamp = Stamp(nn=nn, scale=scale)
    if transform_obj is not None:
        gal_coords = transform_obj.transform(stamp.coords)
    else:
        gal_coords = stamp.coords

    # Sample the galaxy flux
    gal_prof = _gsinterface.getFluxVec(
        scale=scale,
        gsobj=gal_obj._sbp,
        xy_coords=gal_coords
        )

    # No convolution necessary in this case so just return the fluxes
    if draw_method == "no_pixel" and psf_obj is None:
        return gal_prof

    # Construct pixel response
    pixel_response = galsim.Pixel(scale=pix_scale)
    if psf_obj is None:
        psf_obj = pixel_response
    else:
        psf_obj = galsim.Convolve([psf_obj, pixel_response])

    # Convolution in Fourier space
    gal_prof = _gsinterface.convolvePsf(
        scale=scale,
        gsobj=psf_obj._sbp,
        gal_prof=gal_prof,
        downsample_ratio=downsample_ratio,
        ngrid=ngrid
        )

    return gal_prof

def simulate_galaxy_batch(
        ngrid,
        pix_scale,
        gal_obj_list,
        transform_obj=None,
        psf_obj=None,
        truncate_ratio=1.0,
        maximum_num_grids=4096,
        draw_method="auto",
        nproc=4
):
    
    """
    The function samples the surface density field of a galaxy at the grids

    Args:

    ngrid (int):        number of grids
    pix_scale (float):  pixel scale
    gal_obj_list (list):   List of Galsim galaxy objects to sample on the grids
    transform_obj :     Coordinate transform object
    psf_obj (galsim):   Galsim PSF object to smear the image
    truncate_ratio (float):    truncate at truncate_ratio times good_image_size
    maximum_num_grids (int):   maximum number of grids for simulation in real space
    draw_method (str):  method to draw the galaxy image, "auto" will convolve with
                        pixel response, "no_pixel" is as it implies
    nproc (int):        Number of processors to use for multiprocessing. Default is 4
    """

    original_omp_num_threads = os.environ.get('OMP_NUM_THREADS', None)
    os.environ['OMP_NUM_THREADS'] = '1'

    mp.set_start_method('spawn', force=True)

    with mp.Pool(nproc) as p:
        
        args_list = [
                        (
                        ngrid, 
                        pix_scale, 
                        gal_obj, 
                        transform_obj, 
                        psf_obj, 
                        truncate_ratio, 
                        maximum_num_grids, 
                        draw_method
                        ) for gal_obj in gal_obj_list
                    ]
        
        outcome = p.starmap(simulate_galaxy, args_list)

    if original_omp_num_threads is None:
        del os.environ['OMP_NUM_THREADS']
    else:
        os.environ['OMP_NUM_THREADS'] = original_omp_num_threads

    return outcome

