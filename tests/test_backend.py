import numpy as np

from batsim.backend import _get_array_backend, _resolve_array_backend
from batsim.stamp import Stamp


def test_none_backend_selects_numpy():
    assert _get_array_backend() is np
    assert _resolve_array_backend(None) is np


def test_stamp_defaults_to_numpy_backend():
    stamp = Stamp(nn=4, scale=0.2)

    assert stamp.xp is np
    assert stamp.coords.dtype == np.float64
