"""Grid sizing and quadrature helpers for BATSim rendering."""

from dataclasses import dataclass

import galsim
import numpy as np


@dataclass(frozen=True)
class FineGrid:
    """
    Description of the high-resolution grid used by the renderer.

    The fine grid is where BATSim samples the GalSim profile before any
    Fourier-space PSF or pixel convolution. ``fine_compact`` records the
    GalSim-recommended compact support at ``fine_scale`` so callers can inspect
    whether a memory cap forced the simulation grid smaller than that support.
    """

    fine_scale: float
    fine_ngrid: int
    fine_compact: int


def _determine_supersampling(
    obj,
    scale,
    integration_order,
    sim_ngrid,
    pad,
    safety=2.0,
    max_supersample=16,
    min_supersample=4,
    attenuation_target=0.5,
    max_fine_grid=None,
):
    """
    Estimate supersampling from Fourier bandwidth and compact image extent.

    The returned factor is the requested total anti-aliasing budget before it
    is split between sub-pixel integration and FFT supersampling.
    """
    maxk = obj.maxk
    nyquist = np.pi / scale

    raw = safety * maxk / nyquist

    q = max(1, int(integration_order))

    if q <= 1:
        effective = raw
    else:
        # Block integration suppresses high-frequency leakage before the FFT,
        # so a higher quadrature order can safely reduce the FFT-side
        # supersampling budget.
        # attenuation_target lets you tune how much high-k leakage you tolerate.
        integration_reduction = q**2 / attenuation_target**0.5

        # Avoid absurd reductions for high q.
        integration_reduction = min(integration_reduction, 8.0 * q)

        effective = raw / integration_reduction

        # Do not allow the requested effective supersampling to fall below q,
        # otherwise fft_supersample may collapse too aggressively.
        effective = max(effective, q)

    bandwidth_supersample = 2 ** int(np.ceil(np.log2(max(min_supersample, effective))))
    bandwidth_supersample = int(np.clip(bandwidth_supersample, min_supersample, max_supersample))

    if max_fine_grid is None:
        return bandwidth_supersample

    max_fine_grid = int(max_fine_grid)
    trial = bandwidth_supersample
    extent_supersample = min_supersample
    compact_limit = 3.0 * max_fine_grid

    while trial >= min_supersample:
        trial_integration_order = min(q, trial)
        trial_fft_supersample = int(np.ceil(trial / trial_integration_order))
        trial_fft_supersample = max(trial_fft_supersample, min_supersample)
        fine_scale = scale / trial_fft_supersample
        fine_compact = int(obj.getGoodImageSize(fine_scale))
        # GalSim's compact size is allowed to exceed the final grid, but not by
        # so much that a compact source silently creates a very large FFT.
        if fine_compact < compact_limit:
            extent_supersample = trial
            break
        trial //= 2

    supersample = min(bandwidth_supersample, extent_supersample)
    return int(np.clip(supersample, min_supersample, max_supersample))


def _resolve_simulation_ngrid(gal_obj, psf_obj, scale):
    """
    Choose the internal coarse-grid simulation size.

    This is deliberately independent of the requested final output ngrid.
    """
    if psf_obj is None:
        return int(gal_obj.getGoodImageSize(scale))

    conv = galsim.Convolve([gal_obj, psf_obj])
    return int(conv.getGoodImageSize(scale))


def _make_fine_grid(gal_obj, scale, sim_ngrid, supersample, pad, max_fine_grid=None):
    """
    Construct the fine grid used for sampling and convolution.

    The grid is at least large enough for the padded coarse simulation region
    and is rounded up to a whole number of supersampling blocks so later
    downsampling is exact.
    """
    fine_scale = scale / supersample

    requested_fine_ngrid = (sim_ngrid + 2 * pad) * supersample
    fine_compact = int(gal_obj.getGoodImageSize(fine_scale))

    fine_ngrid = max(requested_fine_ngrid, fine_compact)

    remainder = fine_ngrid % supersample
    if remainder:
        fine_ngrid += supersample - remainder

    if max_fine_grid is not None:
        max_fine_grid = int(max_fine_grid)
        if max_fine_grid < 1:
            raise ValueError("max_fine_grid must be positive.")

        capped_fine_ngrid = max_fine_grid - (max_fine_grid % supersample)
        if capped_fine_ngrid < supersample:
            raise ValueError("max_fine_grid must be at least one supersampling block.")

        fine_ngrid = min(fine_ngrid, capped_fine_ngrid)

    return FineGrid(
        fine_scale=fine_scale,
        fine_ngrid=fine_ngrid,
        fine_compact=fine_compact,
    )


def _resolve_integration_sampling(supersample, integration_order):
    """
    Split the anti-aliasing budget between integration and FFT supersampling.

    The product of the returned FFT supersampling and integration order covers
    the requested total supersampling, with the integration order capped so it
    never exceeds that total budget.
    """
    supersample = int(supersample)
    integration_order = int(integration_order)

    if supersample < 1:
        raise ValueError("supersample must be >= 1.")
    if integration_order < 1:
        raise ValueError("integration_order must be >= 1.")

    integration_order = min(integration_order, supersample)
    fft_supersample = int(np.ceil(supersample / integration_order))

    return fft_supersample, integration_order


def _integration_offsets(xp, fine_scale, integration_order, dtype):
    """
    Return sub-pixel offsets and weights for Gauss-Legendre block integration.

    The weights integrate over one fine pixel in normalized coordinates and
    therefore sum to one.
    """
    if integration_order <= 1:
        return None, None

    nodes, weights = np.polynomial.legendre.leggauss(integration_order)
    offsets_1d = xp.asarray(0.5 * nodes * fine_scale, dtype=dtype)
    weights_1d = xp.asarray(0.5 * weights, dtype=dtype)

    yy, xx = xp.meshgrid(offsets_1d, offsets_1d, indexing="ij")
    wy, wx = xp.meshgrid(weights_1d, weights_1d, indexing="ij")

    offsets = xp.stack(
        [
            xx.ravel(),
            yy.ravel(),
        ],
        axis=1,
    )

    return offsets, (wx * wy).ravel()
