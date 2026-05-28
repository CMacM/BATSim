import warnings

import numpy as np


_ARRAY_BACKEND = None
_CPU_FALLBACK_WARNED = False


def _get_array_backend():
    """
    Return CuPy if available, otherwise NumPy.

    This mirrors the backend-selection behaviour used in the simulation code.
    """
    global _ARRAY_BACKEND
    global _CPU_FALLBACK_WARNED

    if _ARRAY_BACKEND is not None:
        return _ARRAY_BACKEND

    try:
        import cupy as cp

        cp.cuda.runtime.getDeviceCount()
        _ARRAY_BACKEND = cp
    except Exception:
        if not _CPU_FALLBACK_WARNED:
            warnings.warn(
                "CuPy unavailable; falling back to NumPy stamp coordinates.",
                RuntimeWarning,
                stacklevel=2,
            )
            _CPU_FALLBACK_WARNED = True
        _ARRAY_BACKEND = np

    return _ARRAY_BACKEND


class Stamp:
    """
    2D coordinate stamp for sampling galaxy profiles.

    If CuPy is available, coordinates are created directly on the GPU.
    Otherwise, NumPy is used.
    """

    def __init__(
        self,
        nn: int = 32,
        scale: float = 0.2,
        backend=None,
        dtype=None,
    ):
        """
        Parameters
        ----------
        nn : int
            Number of grid points in x and y.
        scale : float
            Pixel scale in arcsec.
        backend : module, optional
            Array backend, e.g. numpy or cupy. If None, auto-detect.
        dtype : dtype, optional
            Coordinate dtype. Defaults to float32 for CuPy, float64 for NumPy.
        """
        self.xp = _get_array_backend() if backend is None else backend

        if dtype is None:
            dtype = self.xp.float32 if self.xp is not np else np.float64

        self.dtype = dtype
        self.set_coords(nn, scale)

    def set_coords(self, nn, scale):
        """
        Construct coordinates with shape (2, nn*nn), ordered as [x, y].
        """
        nn = int(nn)
        scale = float(scale)

        xp = self.xp

        ind = xp.arange(
            -(nn // 2),
            (nn + 1) // 2,
            dtype=self.dtype,
        ) * scale

        yy, xx = xp.meshgrid(ind, ind, indexing="ij")

        self.coords = xp.stack(
            [
                xx.ravel(),
                yy.ravel(),
            ],
            axis=0,
        )

        self.scale = scale
        self.pixel_area = scale**2
        self.shape = (nn, nn)
        self.nn = nn

    def to_numpy(self):
        """
        Return coordinates as NumPy array.

        Useful when passing coordinates to CPU-only code such as a C++/GalSim
        sampling backend.
        """
        if self.xp is np:
            return self.coords

        return self.xp.asnumpy(self.coords)
