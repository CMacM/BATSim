import galsim
import numpy as np

from batsim import pltutil


def _constant_image(value, nx=2, ny=2, scale=0.2):
    image = galsim.ImageF(nx, ny, scale=scale)
    image.array[:, :] = value
    return image


def test_stitch_images_square_layout_leaves_unused_cells_blank():
    images = [_constant_image(1), _constant_image(2), _constant_image(3)]

    stitched = pltutil.stitch_images(images, direction="square")

    assert stitched.array.shape == (4, 4)
    assert stitched.scale == images[0].scale

    np.testing.assert_allclose(
        stitched.subImage(galsim.BoundsI(xmin=1, xmax=2, ymin=1, ymax=2)).array,
        1.0,
    )
    np.testing.assert_allclose(
        stitched.subImage(galsim.BoundsI(xmin=3, xmax=4, ymin=1, ymax=2)).array,
        2.0,
    )
    np.testing.assert_allclose(
        stitched.subImage(galsim.BoundsI(xmin=1, xmax=2, ymin=3, ymax=4)).array,
        3.0,
    )
    np.testing.assert_allclose(
        stitched.subImage(galsim.BoundsI(xmin=3, xmax=4, ymin=3, ymax=4)).array,
        0.0,
    )
