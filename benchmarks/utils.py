"""Shared fixtures and helpers for BATSim ASV benchmarks."""

import gc

import galsim
import numpy as np

import batsim

SCALE = 0.2
HLR = 0.7
FLUX = 1.0
DEFAULT_MAX_FINE_GRID = 2048


def gaussian_galaxy():
    """Return a cheap analytic galaxy profile."""
    return galsim.Gaussian(sigma=0.55, flux=FLUX)


def sersic_galaxy(n=1.0, half_light_radius=HLR):
    """Return a representative analytic Sersic galaxy."""
    return galsim.Sersic(n=n, half_light_radius=half_light_radius, flux=FLUX)


def compact_sersic_galaxy():
    """Return a compact high-n profile that stresses supersampling decisions."""
    return galsim.Sersic(n=4.0, half_light_radius=0.35, flux=FLUX)


def elliptical_sersic_galaxy():
    """Return a sheared profile that stresses compact image support sizing."""
    return galsim.Sersic(n=3.0, half_light_radius=0.8, flux=FLUX).shear(e1=0.75)


def gaussian_psf():
    """Return a cheap Gaussian PSF."""
    return galsim.Gaussian(fwhm=0.55, flux=1.0)


def moffat_psf():
    """Return a broader Moffat PSF for convolution benchmarks."""
    return galsim.Moffat(beta=3.5, fwhm=0.75, trunc=3.0, flux=1.0)


def transform_for_name(name):
    """Build one of the public transform objects used in render benchmarks."""
    if name == "none":
        return None
    if name == "lens":
        return batsim.LensTransform(gamma1=0.12, gamma2=0.04, kappa=0.0)
    if name == "ia":
        return batsim.IaTransform(scale=SCALE, hlr=HLR, A=0.04, beta=0.8, phi=0.2)
    if name == "flexion":
        return batsim.FlexionTransform(
            gamma1=0.05,
            gamma2=0.02,
            kappa=0.0,
            F1=0.005,
            F2=-0.003,
            G1=0.002,
            G2=0.001,
        )

    raise ValueError(f"Unknown transform benchmark name: {name}")


def render_kwargs(
    *,
    gal_obj=None,
    ngrid=128,
    transform_obj=None,
    psf_obj=None,
    draw_method="auto",
    integration_order=2,
    psf_mode="kvalue",
    backend="np",
    precision="single",
    force_input_flux=True,
    compensate_integration="quadrature",
    use_true_center=True,
    max_fine_grid=DEFAULT_MAX_FINE_GRID,
):
    """Return current-public-API render kwargs shared by benchmarks."""
    return {
        "gal_obj": sersic_galaxy() if gal_obj is None else gal_obj,
        "scale": SCALE,
        "ngrid": int(ngrid),
        "transform_obj": transform_obj,
        "psf_obj": psf_obj,
        "draw_method": draw_method,
        "integration_order": int(integration_order),
        "psf_mode": psf_mode,
        "backend": backend,
        "precision": precision,
        "force_input_flux": force_input_flux,
        "compensate_integration": compensate_integration,
        "use_true_center": use_true_center,
        "max_fine_grid": max_fine_grid,
    }


def coordinate_grid(nn, scale=SCALE, dtype=np.float64):
    """Return flattened ``(2, nn * nn)`` NumPy coordinates for transform benchmarks."""
    ind = (np.arange(nn, dtype=dtype) - 0.5 * (nn - 1)) * scale
    yy, xx = np.meshgrid(ind, ind, indexing="ij")
    return np.stack([xx.ravel(), yy.ravel()], axis=0)


def grid_decision(gal_obj, integration_order=2, ngrid=128, max_fine_grid=DEFAULT_MAX_FINE_GRID):
    """Return the grid sizing decisions used by the public renderer."""
    import batsim.sim as sim

    sim_ngrid = sim._resolve_simulation_ngrid(gal_obj, None, SCALE)
    supersample = sim._determine_supersampling(
        gal_obj,
        SCALE,
        integration_order,
        sim_ngrid=sim_ngrid,
        pad=16,
        max_supersample=64,
        min_supersample=4,
        max_fine_grid=max_fine_grid,
    )
    fft_supersample, resolved_integration_order = sim._resolve_integration_sampling(
        supersample=supersample,
        integration_order=integration_order,
    )
    fft_supersample = max(fft_supersample, 4)
    grid = sim._make_fine_grid(
        gal_obj=gal_obj,
        scale=SCALE,
        sim_ngrid=sim_ngrid,
        supersample=fft_supersample,
        pad=16,
        max_fine_grid=max_fine_grid,
    )

    return {
        "output_ngrid": int(ngrid),
        "sim_ngrid": int(sim_ngrid),
        "supersample": int(supersample),
        "fft_supersample": int(fft_supersample),
        "integration_order": int(resolved_integration_order),
        "fine_ngrid": int(grid.fine_ngrid),
        "fine_compact": int(grid.fine_compact),
    }


def flux_residual(image, target_flux=FLUX):
    """Return absolute flux residual for a rendered image."""
    return abs(float(np.sum(image)) - float(target_flux))


def get_cupy():
    """Return CuPy if a working CUDA runtime is available, otherwise skip."""
    try:
        import cupy as cp

        cp.cuda.runtime.getDeviceCount()
        return cp
    except Exception as exc:
        raise NotImplementedError("CuPy/CUDA is not available for GPU benchmarks.") from exc


def cleanup_backend(backend="np"):
    """Release transient memory after benchmarks that allocate large arrays."""
    gc.collect()
    try:
        batsim.clear_backend_memory(backend)
    except RuntimeError:
        if backend not in ("cp", "cupy"):
            raise
