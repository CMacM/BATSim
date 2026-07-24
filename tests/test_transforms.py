import galsim
import numpy as np
import pytest

import batsim


class IdentityTransform(batsim.Transform):
    def _transform_relative(self, coords_relative, xp, dtype):
        return coords_relative


class ScaleXTransform(batsim.Transform):
    def __init__(self, factor, **kwargs):
        super().__init__(**kwargs)
        self.factor = factor

    def _transform_relative(self, coords_relative, xp, dtype):
        x, y = coords_relative
        return xp.stack([self.factor * x, y], axis=0)


def test_builtin_transforms_inherit_transform_base():
    assert issubclass(batsim.LensTransform, batsim.Transform)
    assert issubclass(batsim.IaTransform, batsim.Transform)
    assert issubclass(batsim.FlexionTransform, batsim.Transform)


def test_transform_base_requires_subclass_logic():
    transform = batsim.Transform()
    coords = np.zeros((2, 4))

    with pytest.raises(NotImplementedError, match="_transform_relative"):
        transform.transform(coords)


def test_custom_transform_gets_center_and_coordinate_management():
    transform = ScaleXTransform(factor=2.0, center=[1.0, -0.5])
    coords = np.array(
        [
            [0.0, 1.0, 2.0],
            [-0.5, 0.0, 0.5],
        ]
    )

    expected = np.array(
        [
            [-1.0, 1.0, 3.0],
            [-0.5, 0.0, 0.5],
        ]
    )

    np.testing.assert_allclose(transform.transform(coords), expected)


def test_transform_to_backend_preserves_behavior_and_dtype():
    transform = ScaleXTransform(factor=1.4, center=[0.25, -0.75])
    copied = transform.to_backend(dtype=np.float32)

    coords = np.array(
        [
            [-1.0, 0.0, 1.0],
            [-0.5, 0.5, 1.5],
        ],
        dtype=np.float32,
    )

    assert copied is not transform
    assert copied.dtype == np.float32
    assert copied.ref_vec.dtype == np.float32
    np.testing.assert_allclose(copied.transform(coords), transform.transform(coords))


def test_custom_transform_can_be_used_in_simulate_galaxy():
    gal = galsim.Gaussian(half_light_radius=0.8, flux=5.0)
    transform = IdentityTransform()

    reference = batsim.simulate_galaxy(
        ngrid=48,
        pix_scale=0.2,
        gal_obj=gal,
        draw_method="no_pixel",
        use_true_center=False,
    )
    rendered = batsim.simulate_galaxy(
        ngrid=48,
        pix_scale=0.2,
        gal_obj=gal,
        transform_obj=transform,
        draw_method="no_pixel",
        use_true_center=False,
    )

    np.testing.assert_allclose(rendered, reference)
