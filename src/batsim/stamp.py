import numpy as np

class Stamp(object):
    def __init__(self, nn: int = 32, scale: float = 0.2):
        """Initialize the 2D stamp object. This class enables distorting
        an image by changing the samplinng position with non-affine
        transformation

        Args:
        nn (int):      number of grids on x and y direction
        scale (float): pixel scale in units of arcsec
        """
        # Set up the grids
        self.set_coords(nn, scale)
        return

    def set_coords(self, nn, scale):
        inds_1d = (np.arange(nn, dtype=float) - 0.5 * (nn - 1)) * scale
        yy, xx = np.meshgrid(inds_1d, inds_1d, indexing="ij")
        self.coords = np.stack([xx.ravel(), yy.ravel()], axis=0)
        self.scale = scale
        self.pixel_area = self.scale**2.0
        self.shape = (nn, nn)
        self.nn = nn
        return
