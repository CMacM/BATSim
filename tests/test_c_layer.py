import batsim
import galsim
import numpy as np
import pytest


def _expected_flux_image(gsobj, coords, scale):
    return np.array(
        [gsobj.xValue(x, y) * scale * scale for x, y in coords.T],
        dtype=np.float64,
    ).reshape(2, 2)


def test_get_flux_vec():
    scale = 0.2
    coords = np.array(
        [
            [-0.1, 0.1, -0.1, 0.1],
            [-0.2, -0.2, 0.2, 0.2],
        ],
        dtype=np.float64,
    )

    sb_obj = galsim.Sersic(n=4, half_light_radius=0.5)
    trans_obj = sb_obj.shear(g1=0.1, g2=0.2).shift(0.5, 0.5)

    sb_flux = batsim._gsinterface.getFluxVec(scale=scale, gsobj=sb_obj._sbp, xy_coords=coords)
    trans_flux = batsim._gsinterface.getFluxVec(scale=scale, gsobj=trans_obj._sbp, xy_coords=coords)

    np.testing.assert_allclose(sb_flux, _expected_flux_image(sb_obj, coords, scale))
    np.testing.assert_allclose(trans_flux, _expected_flux_image(trans_obj, coords, scale))


def test_get_flux_vec_rejects_non_square_coordinate_count():
    sb_obj = galsim.Sersic(n=4, half_light_radius=0.5)
    coords = np.array([[0.1, 0.1], [0.2, 0.2]], dtype=np.float64)

    with pytest.raises(RuntimeError, match=r"xy_coords\.shape\[1\] must be a perfect square"):
        batsim._gsinterface.getFluxVec(scale=0.2, gsobj=sb_obj._sbp, xy_coords=coords)


if __name__ == "__main__":
    test_get_flux_vec()
