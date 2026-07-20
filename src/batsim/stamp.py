import numpy as np

from .backend import _get_array_backend


class Stamp:
    """
    Coordinate grid used for sampling GalSim profiles.

    A ``Stamp`` stores flattened ``x``/``y`` coordinates with shape
    ``(2, nn * nn)``.  Coordinates are built on the requested array backend so
    transform operations can run on either NumPy or CuPy arrays.

    Parameters
    ----------
    nn : int, optional
        Number of pixels along each side of the square coordinate grid.
    scale : float, optional
        Pixel scale in arcsec.
    backend : module, optional
        Array backend, usually ``numpy`` or ``cupy``.  If omitted, BATSim
        auto-detects CuPy and falls back to NumPy.
    dtype : dtype, optional
        Coordinate dtype.  Defaults to ``float32`` for CuPy and ``float64`` for
        NumPy.
    use_true_center : bool, optional
        If True, use GalSim's true-image-center convention.  When the fine grid
        will later be downsampled, this aligns it to the true center of the
        coarse output grid.
    downsample_ratio : int, optional
        Ratio between the fine grid and the eventual coarse output grid.

    Attributes
    ----------
    coords : array-like
        Flattened coordinate array ordered as ``[x, y]`` with shape
        ``(2, nn * nn)``.
    pixel_area : float
        Fine-grid pixel area in square arcsec.
    center_index : float
        Pixel index used as the coordinate origin.
    """

    def __init__(
        self,
        nn: int = 32,
        scale: float = 0.2,
        backend=None,
        dtype=None,
        use_true_center=True,
        downsample_ratio=1,
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
        use_true_center : bool, optional
            If True, use GalSim's default true-image-center convention.  If the
            stamp will later be downsampled, this aligns the fine grid to the
            true center of the eventual coarse grid.
        downsample_ratio : int, optional
            Coarse-to-fine pixel ratio used for true-center alignment.
        """
        self.xp = (
            _get_array_backend("CuPy unavailable; falling back to NumPy stamp coordinates.")
            if backend is None
            else backend
        )

        if dtype is None:
            dtype = self.xp.float32 if self.xp is not np else np.float64

        self.dtype = dtype
        self.set_coords(
            nn,
            scale,
            use_true_center=use_true_center,
            downsample_ratio=downsample_ratio,
        )

    def set_coords(self, nn, scale, use_true_center=True, downsample_ratio=1):
        """
        Construct coordinates with shape (2, nn*nn), ordered as [x, y].

        Parameters
        ----------
        nn : int
            Number of grid points along each side.
        scale : float
            Pixel scale in arcsec.
        use_true_center : bool, optional
            Whether to use GalSim's true-image-center convention.
        downsample_ratio : int, optional
            Coarse-to-fine pixel ratio used for true-center alignment.
        """
        nn = int(nn)
        scale = float(scale)
        downsample_ratio = int(downsample_ratio)

        xp = self.xp

        if downsample_ratio < 1:
            raise ValueError("downsample_ratio must be >= 1.")

        if use_true_center:
            center_index = 0.5 * (nn - downsample_ratio)
        else:
            center_index = nn // 2

        ind = (xp.arange(nn, dtype=self.dtype) - center_index) * scale

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
        self.use_true_center = bool(use_true_center)
        self.downsample_ratio = downsample_ratio
        self.center_index = center_index

    def to_numpy(self):
        """
        Return coordinates as NumPy array.

        Useful when passing coordinates to CPU-only code such as a C++/GalSim
        sampling backend.

        Returns
        -------
        ndarray
            Coordinate array with shape ``(2, nn * nn)``.
        """
        if self.xp is np:
            return self.coords

        return self.xp.asnumpy(self.coords)
