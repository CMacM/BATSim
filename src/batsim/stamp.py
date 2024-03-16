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
        indx = np.arange(-int(nn / 2), int((nn + 1) / 2), 1) * scale
        indy = np.arange(-int(nn / 2), int((nn + 1) / 2), 1) * scale
        inds = np.meshgrid(indy, indx, indexing="ij")
        # coords in shape of (2, npoints), in order of [x, y]
        self.coords = np.vstack([np.ravel(_) for _ in inds[::-1]])
        self.scale = scale
        self.pixel_area = self.scale**2.0
        self.shape = (nn, nn)
        self.nn = nn
        return
