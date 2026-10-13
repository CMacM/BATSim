# Tests the composition order of BATSim transform lists against GalSim.
# Note: this test script was written with help from Claude Code.
#
# BATSim applies transforms to output-plane coordinates in list order, so
# transform_obj=[t1, t2] samples the source at t2(t1(x)).  For a pull map this
# means the LAST transform in the list is the FIRST physical distortion.  To
# model IA physically first and lensing second, pass [LensTransform, IaTransform].
#
# IaTransform(beta=0) is a constant reduced shear A * exp(2i phi), which GalSim
# can reproduce exactly.  Exaggerated shears are used so that the
# non-commuting O(g_L g_IA) cross terms are large enough to tell the orders apart.

import galsim
import numpy as np
import pytest

import batsim

scale = 0.1
nn = 128
hlr = 0.7
gal = galsim.Gaussian(half_light_radius=hlr, flux=25.0)

g_ia = 0.2
phi = 0.4
gamma1 = 0.18
gamma2 = -0.10
kappa = 0.05

g1 = gamma1 / (1.0 - kappa)
g2 = gamma2 / (1.0 - kappa)
mu = 1.0 / ((1.0 - kappa) ** 2 - gamma1**2 - gamma2**2)
ia_shear = galsim.Shear(g=g_ia, beta=phi * galsim.radians)


def _moments(image):
    n = image.shape[0]
    ind = (np.arange(n) - 0.5 * (n - 1)) * scale
    yy, xx = np.meshgrid(ind, ind, indexing="ij")

    flux = image.sum()
    x0 = (image * xx).sum() / flux
    y0 = (image * yy).sum() / flux
    dx = xx - x0
    dy = yy - y0
    qxx = (image * dx * dx).sum() / flux
    qyy = (image * dy * dy).sum() / flux
    qxy = (image * dx * dy).sum() / flux
    trace = qxx + qyy

    return np.array([(qxx - qyy) / trace, 2.0 * qxy / trace, trace])


def _draw_galsim(obj):
    return obj.drawImage(nx=nn, ny=nn, scale=scale, method="no_pixel", use_true_center=False).array


def _render(obj, transforms):
    return batsim.simulate_galaxy(
        gal_obj=obj,
        scale=scale,
        ngrid=nn,
        transform_obj=transforms,
        draw_method="no_pixel",
        use_true_center=False,
    )


def _lens():
    return batsim.LensTransform(gamma1=gamma1, gamma2=gamma2, kappa=kappa)


def _ia():
    return batsim.IaTransform(scale=scale, hlr=hlr, A=g_ia, beta=0.0, phi=phi)


def _max_rel_residual(image, reference):
    return np.abs(image - reference).max() / reference.max()


@pytest.fixture(scope="module")
def ia_then_lens_reference():
    """GalSim render of the physical chain: IA in the source plane, then lensing."""
    return _draw_galsim(gal.shear(ia_shear).lens(g1=g1, g2=g2, mu=mu))


@pytest.fixture(scope="module")
def lens_then_ia_reference():
    """GalSim render of lensing first, then IA applied in the image plane."""
    return _draw_galsim(gal.lens(g1=g1, g2=g2, mu=mu).shear(ia_shear))


def test_galsim_references_are_distinguishable(ia_then_lens_reference, lens_then_ia_reference):
    # Guard against the test shears being too small to separate the two orders.
    assert _max_rel_residual(lens_then_ia_reference, ia_then_lens_reference) > 1.0e-2
    delta_e = _moments(lens_then_ia_reference)[:2] - _moments(ia_then_lens_reference)[:2]
    assert np.max(np.abs(delta_e)) > 1.0e-2


def test_lens_then_ia_list_order(ia_then_lens_reference):
    rendered = _render(gal, [_lens(), _ia()])

    assert _max_rel_residual(rendered, ia_then_lens_reference) < 1.0e-4
    np.testing.assert_allclose(rendered.sum(), ia_then_lens_reference.sum(), rtol=1.0e-5)
    np.testing.assert_allclose(_moments(rendered), _moments(ia_then_lens_reference), atol=1.0e-5)


def test_ia_then_lens_list_order(
    ia_then_lens_reference,
    lens_then_ia_reference,
):
    rendered = _render(gal, [_ia(), _lens()])

    # Wrong physical order: visibly fails the IA-first reference...
    assert _max_rel_residual(rendered, ia_then_lens_reference) > 1.0e-2
    delta_e = _moments(rendered)[:2] - _moments(ia_then_lens_reference)[:2]
    assert np.max(np.abs(delta_e)) > 1.0e-2

    # ...because it is exactly lensing first, then IA in the image plane.
    assert _max_rel_residual(rendered, lens_then_ia_reference) < 1.0e-4
    np.testing.assert_allclose(_moments(rendered), _moments(lens_then_ia_reference), atol=1.0e-5)


def test_galsim_lensed_input_with_batsim_ia_fails(
    ia_then_lens_reference,
    lens_then_ia_reference,
):
    # A GalSim-side transformation on gal_obj is evaluated after all BATSim
    # transforms in the pull map, so it is always the first physical distortion.
    rendered = _render(gal.lens(g1=g1, g2=g2, mu=mu), [_ia()])

    assert _max_rel_residual(rendered, ia_then_lens_reference) > 1.0e-2
    assert _max_rel_residual(rendered, lens_then_ia_reference) < 1.0e-4
    np.testing.assert_allclose(_moments(rendered), _moments(lens_then_ia_reference), atol=1.0e-5)
