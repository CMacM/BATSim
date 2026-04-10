import galsim
import numpy as np
import os
import multiprocessing as mp
from functools import lru_cache

from . import _gsinterface
from .stamp import Stamp

from time import perf_counter

_PSF_REGISTRY = {}

def _register_psf(psf_obj):
    # Keep memory bounded and invalidate stale cached entries.
    if len(_PSF_REGISTRY) >= 128:
        _PSF_REGISTRY.clear()
        _cached_pad_arcsec.cache_clear()
        _cached_effective_psf.cache_clear()
    key = id(psf_obj)
    _PSF_REGISTRY[key] = psf_obj
    return key

@lru_cache(maxsize=128)
def _cached_pad_arcsec(psf_key, pix_scale):
    return _PSF_REGISTRY[psf_key].calculateMomentRadius(size=32, scale=pix_scale / 4.0)

@lru_cache(maxsize=128)
def _cached_effective_psf(psf_key):
    psf_obj = _PSF_REGISTRY[psf_key]
    return psf_obj

def simulate_galaxy(
    gal_obj,
    pix_scale,
    ngrid = None,
    transform_obj=None,
    psf_obj=None,
    truncate_ratio=1.0,
    maximum_num_grids=4096,
    draw_method="auto",
    force_ngrid=False,
    delta_image_x=0.0,
    delta_image_y=0.0,
    profile=False,
):
    """The function samples the surface density field of a galaxy at the grids
    This function only conduct sampling; PSF and pixel response are not
    included.

    Args:
    ngrid (int):        number of grids
    pix_scale (float):  pixel scale
    gal_obj (galsim):   Galsim galaxy object to sample on the grids
    transform_obj :     Coordinate transform object or list of transform in order
                        that they should be applied.
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
    profile (bool):     If True, print per-galaxy timings and stats
    Returns:
    outcome (ndarray):  2D galaxy image on the grids
    """
    def _log(msg):
        if profile:
            print(f"[simulate_galaxy] {msg}")

    t_total = perf_counter() if profile else None
    gobj = gal_obj.shift(
        delta_image_x * pix_scale,
        delta_image_y * pix_scale,
    )
    psf_key = _register_psf(psf_obj) if psf_obj is not None else None

    # Initialize variables based on PSF presence
    downsample_ratio = 1
    pad_arcsec = 0.0
    if psf_obj is None and draw_method == "no_pixel":
        # In this case we just get the fluxes for the requested stamp size
        scale = pix_scale
        nn = int(ngrid)
    else:
        t = perf_counter() if profile else None
        # Compute the effective scale for simulation
        if psf_obj is None:
            scale = pix_scale / 4.0
        else:
            scale = min(gobj.nyquist_scale, pix_scale / 4.0)
            pad_arcsec = _cached_pad_arcsec(psf_key, pix_scale)
            downsample_ratio = min(int(2 ** np.ceil(np.log2(pix_scale / scale))), 128)
        if profile:
            _log(f"effective_scale_padding={perf_counter() - t:.4e}s")
        scale = pix_scale / downsample_ratio

        t = perf_counter() if profile else None
        # Calculate the number of grids considering padding and truncation
        npad = int(pad_arcsec / scale + 0.5) * 4
        nn = npad * 2 + min(
            gobj.getGoodImageSize(scale)
            * truncate_ratio, ngrid * downsample_ratio
        )
        nn = min(int(2 ** np.ceil(np.log2(nn))), maximum_num_grids)
        if profile:
            _log(f"grid_sizing={perf_counter() - t:.4e}s")

    if force_ngrid and nn < ngrid:
        nn = ngrid
        scale = pix_scale

    # Initialize and Distort Coordinates in order
    t_section = perf_counter() if profile else None
    stamp = Stamp(nn=nn, scale=scale)

    # Check if transform_obj is a list of transforms and apply them in order
    if isinstance(transform_obj, list):
        gal_coords = stamp.coords
        for trf in transform_obj:
            gal_coords = trf.transform(gal_coords)

    # If transform_obj is a single transform, apply it directly
    elif transform_obj is not None:
        gal_coords = transform_obj.transform(stamp.coords)

    # If no transform is provided, use the original coordinates
    else:
        gal_coords = stamp.coords

    if profile:
        _log(f"coords_transform={perf_counter() - t_section:.4e}s")

    # Sample the galaxy flux
    t = perf_counter() if profile else None
    gal_prof = _gsinterface.getFluxVec(
        scale=scale,
        gsobj=gobj._sbp,
        xy_coords=gal_coords
    )
    if profile:
        _log(f"sample_flux={perf_counter() - t:.4e}s")

    # No convolution necessary in this case so just return the fluxes
    if draw_method == "no_pixel":
        if psf_obj is None:
            if profile:
                _log(
                    f"stats nn={nn} downsample_ratio={downsample_ratio} scale={scale:.4e} "
                    f"pad_arcsec={pad_arcsec:.4e} draw_method={draw_method}"
                )
                _log(f"total={perf_counter() - t_total:.4e}s")
            return gal_prof
        else:
            pass
    elif draw_method == "auto":
        t = perf_counter() if profile else None
        if psf_obj is None:
            psf_obj = galsim.Pixel(scale=pix_scale)
        else:
            psf_obj = _cached_effective_psf(psf_key)
            psf_obj = galsim.Convolve([psf_obj, galsim.Pixel(scale=pix_scale)])
        if profile:
            _log(f"prepare_psf={perf_counter() - t:.4e}s")
    else:
        raise ValueError("do not support draw_method=%s" %draw_method)
    
    t = perf_counter() if profile else None
    # Convolution in Fourier space
    gal_prof = _gsinterface.convolvePsf(
        scale=scale,
        gsobj=psf_obj._sbp,
        gal_prof=gal_prof,
        downsample_ratio=downsample_ratio,
        ngrid=ngrid
    )
    if profile:
        _log(f"convolution_downsample={perf_counter() - t:.4e}s")
        _log(
            f"stats nn={nn} downsample_ratio={downsample_ratio} scale={scale:.4e} "
            f"pad_arcsec={pad_arcsec:.4e} draw_method={draw_method}"
        )
        _log(f"total={perf_counter() - t_total:.4e}s")
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
        force_ngrid=False,
        profile=False
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
    profile (bool):     If True, enable per-galaxy profiling logs in workers
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
                        force_ngrid,
                        0.0,
                        0.0,
                        profile
                        ) for gal_obj in gal_obj_list
                    ]

        outcome = p.starmap(simulate_galaxy, args_list)

    if original_omp_num_threads is None:
        del os.environ['OMP_NUM_THREADS']
    else:
        os.environ['OMP_NUM_THREADS'] = original_omp_num_threads

    return outcome