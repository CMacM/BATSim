import galsim
import numpy as np
import os
import multiprocessing as mp

from . import _gsinterface
from .stamp import Stamp
import heapq


def _round_up_multiple(n: int, m: int) -> int:
    return n if n % m == 0 else n + (m - n % m)


def next_235_heap(B: int):
    if B < 1:
        return 1  # smallest 5-smooth is 1
    primes = (2, 3, 5)
    seen = {1}
    h = [1]
    while h:
        x = heapq.heappop(h)
        if x >= B:
            return x
        for p in primes:
            y = x * p
            if y not in seen:
                seen.add(y)
                heapq.heappush(h, y)


def simulate_galaxy(
    ngrid,
    pix_scale,
    gal_obj,
    transform_obj=None,
    psf_obj=None,
    truncate_ratio=1.0,
    maximum_num_grids=4096,
    draw_method="auto",
    force_ngrid=False,
    delta_image_x=0.0,
    delta_image_y=0.0,
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
    draw_method (str):  method to draw the galaxy image, "auto" will convolve
                        with pixel response, "no_pixel" is as it implies
    force_ngrid (bool): If True, force the number of grids to be ngrid even if
                        a smaller number of grids is sufficient for the
                        simulation
    Returns:
    outcome (ndarray):  2D galaxy image on the grids
    """

    gobj = gal_obj.shift(
        delta_image_x * pix_scale,
        delta_image_y * pix_scale,
    )
    if draw_method == "auto":
        # Construct pixel response
        pixel_response = galsim.Pixel(scale=pix_scale)
        if psf_obj is None:
            psf_obj = pixel_response
        else:
            psf_obj = galsim.Convolve([psf_obj, pixel_response])

    # Initialize variables based on PSF presence
    if psf_obj is None:
        # In this case we just get the fluxes for the requested stamp size
        scale = pix_scale
        nn = int(ngrid)
        downsample_ratio = 1
    else:
        scale = min(gobj.nyquist_scale, pix_scale / 4.0)
        # downsample_ratio = min(
        #     int(2 ** np.ceil(np.log2(pix_scale / scale))),
        #     128,
        # )
        downsample_ratio = min(
            next_235_heap(int(pix_scale / scale + 0.5)),
            100,
        )
        scale = pix_scale / int(downsample_ratio)
        npad = int(1.5 / scale + 0.5)  # pad 1.5 arcsec
        nn = npad * 2 + min(
            gobj.getGoodImageSize(scale) * truncate_ratio,
            ngrid * int(downsample_ratio),
        )
        # nn = min(int(2 ** np.ceil(np.log2(nn))), maximum_num_grids)
        nn = min(nn, maximum_num_grids)
        nn = int(_round_up_multiple(nn, downsample_ratio * 2))

    if force_ngrid and nn < ngrid:
        nn = ngrid
        scale = pix_scale
        downsample_ratio = 1

    # Initialize and Distort Coordinates
    stamp = Stamp(nn=nn, scale=scale)
    if transform_obj is not None:
        gal_coords = transform_obj.transform(stamp.coords)
    else:
        gal_coords = stamp.coords
    # Sample the galaxy flux
    gal_prof = _gsinterface.getFluxVec(
        scale=scale,
        gsobj=gobj._sbp,
        xy_coords=gal_coords
    )
    if psf_obj is not None:
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
        nproc=4,
        force_ngrid=False
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
                        draw_method,
                        force_ngrid
                        ) for gal_obj in gal_obj_list
                    ]

        outcome = p.starmap(simulate_galaxy, args_list)

    if original_omp_num_threads is None:
        del os.environ['OMP_NUM_THREADS']
    else:
        os.environ['OMP_NUM_THREADS'] = original_omp_num_threads

    return outcome
