import importlib

import galsim

sim = importlib.import_module("batsim.sim")


class CompactSizeProbe:
    maxk = 1000.0

    def getGoodImageSize(self, scale):
        return int(200.0 / scale)


def _supersample(obj, scale, integration_order, sim_ngrid, pad, max_fine_grid):
    return sim._determine_supersampling(
        obj,
        scale,
        integration_order,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_supersample=32,
        min_supersample=4,
        max_fine_grid=max_fine_grid,
    )


def test_compact_high_n_keeps_bandwidth_supersample_with_large_budget():
    scale = 0.2
    pad = 16
    gal = galsim.Sersic(n=6, half_light_radius=0.12, flux=1.0)
    sim_ngrid = sim._resolve_simulation_ngrid(gal, None, scale)

    bandwidth_supersample = _supersample(
        gal,
        scale,
        integration_order=2,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=None,
    )
    extent_supersample = _supersample(
        gal,
        scale,
        integration_order=2,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=4096,
    )

    assert bandwidth_supersample == 32
    assert extent_supersample == bandwidth_supersample


def test_compact_elliptical_galaxy_keeps_supersample_inside_compact_limit():
    scale = 0.2
    pad = 16
    max_fine_grid = 1500
    gal = galsim.Sersic(n=6, half_light_radius=0.12, flux=1.0).shear(e1=0.85)
    sim_ngrid = sim._resolve_simulation_ngrid(gal, None, scale)

    bandwidth_supersample = _supersample(
        gal,
        scale,
        integration_order=1,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=None,
    )
    extent_supersample = _supersample(
        gal,
        scale,
        integration_order=1,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=max_fine_grid,
    )
    fine_compact = int(gal.getGoodImageSize(scale / extent_supersample))

    assert bandwidth_supersample == 32
    assert extent_supersample == bandwidth_supersample
    assert fine_compact < 3 * max_fine_grid


def test_elliptical_galaxy_reduces_until_inside_compact_limit():
    scale = 0.2
    pad = 16
    max_fine_grid = 1500
    gal = galsim.Sersic(n=3, half_light_radius=1.0, flux=1.0).shear(e1=0.85)
    sim_ngrid = sim._resolve_simulation_ngrid(gal, None, scale)

    bandwidth_supersample = _supersample(
        gal,
        scale,
        integration_order=1,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=None,
    )
    extent_supersample = _supersample(
        gal,
        scale,
        integration_order=1,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=max_fine_grid,
    )

    assert bandwidth_supersample == 32
    assert extent_supersample == 8
    assert extent_supersample < bandwidth_supersample

    accepted_fine_compact = int(gal.getGoodImageSize(scale / extent_supersample))
    rejected_fine_compact = int(gal.getGoodImageSize(scale / (extent_supersample * 2)))

    assert accepted_fine_compact < 3 * max_fine_grid
    assert rejected_fine_compact >= 3 * max_fine_grid


def test_extent_limit_uses_fft_supersample_after_integration_split():
    scale = 0.2
    max_fine_grid = 4096
    obj = CompactSizeProbe()

    supersample = _supersample(
        obj,
        scale,
        integration_order=2,
        sim_ngrid=1340,
        pad=16,
        max_fine_grid=max_fine_grid,
    )
    fft_supersample, _ = sim._resolve_integration_sampling(
        supersample,
        integration_order=2,
    )
    rejected_fft_supersample, _ = sim._resolve_integration_sampling(
        supersample * 2,
        integration_order=2,
    )

    assert supersample == 16
    assert fft_supersample == 8
    assert obj.getGoodImageSize(scale / fft_supersample) < 3 * max_fine_grid
    assert obj.getGoodImageSize(scale / rejected_fft_supersample) >= 3 * max_fine_grid


def test_extent_fallback_returns_min_supersample_and_leaves_safety_clip():
    scale = 0.2
    pad = 16
    max_fine_grid = 300
    gal = galsim.Sersic(n=3, half_light_radius=1.0, flux=1.0).shear(e1=0.85)
    sim_ngrid = sim._resolve_simulation_ngrid(gal, None, scale)

    supersample = _supersample(
        gal,
        scale,
        integration_order=1,
        sim_ngrid=sim_ngrid,
        pad=pad,
        max_fine_grid=max_fine_grid,
    )
    grid = sim._make_fine_grid(
        gal,
        scale,
        sim_ngrid,
        supersample,
        pad,
        max_fine_grid=max_fine_grid,
    )

    requested_fine_ngrid = (sim_ngrid + 2 * pad) * supersample

    assert supersample == 4
    assert grid.fine_compact >= 3 * max_fine_grid
    assert requested_fine_ngrid > max_fine_grid
    assert grid.fine_ngrid == max_fine_grid
    assert grid.fine_ngrid < requested_fine_ngrid


def test_typical_sersic_profiles_are_unchanged_by_default_extent_budget():
    scale = 0.2
    pad = 16
    galaxies = [
        galsim.Sersic(n=4.0, half_light_radius=0.5, flux=1.0),
        galsim.Sersic(n=2.0, half_light_radius=0.8, flux=1.0).shear(e1=0.3),
        galsim.Sersic(n=5.5, half_light_radius=0.2, flux=1.0),
    ]

    for gal in galaxies:
        sim_ngrid = sim._resolve_simulation_ngrid(gal, None, scale)
        assert _supersample(
            gal,
            scale,
            integration_order=2,
            sim_ngrid=sim_ngrid,
            pad=pad,
            max_fine_grid=4096,
        ) == _supersample(
            gal,
            scale,
            integration_order=2,
            sim_ngrid=sim_ngrid,
            pad=pad,
            max_fine_grid=None,
        )
