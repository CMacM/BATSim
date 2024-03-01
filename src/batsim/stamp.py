import galsim
import numpy as np
import warnings
from . import _gsinterface

class Stamp(object):
    def __init__(self, coords=None, nn=32, scale=0.2, centering='fpfs'):
        """Initialize the 2D stamp object. This class enables distorting
        an image by changing the samplinng position with non-affine
        transformation

        Args:
            nn (int):      number of grids on x and y direction
            scale (float): pixel scale in units of arcsec
        """
        self.scale = scale
        if centering  == 'galsim':
            self.centering = 0.5*scale
        elif centering == 'fpfs':
            self.centering = 0
        else:
            self.centering = centering
        
        if coords is None:
            indx = (np.arange(-int(nn / 2),
                             int((nn + 1) / 2), 1) * scale)
            indy = (np.arange(-int(nn / 2),
                             int((nn + 1) / 2), 1) * scale)
            inds = np.meshgrid(indy, indx, indexing="ij")
            self.coords = np.vstack([np.ravel(_) for _ in inds[::-1]])
        else:
            self.coords = coords
        self.pixel_values = None
        self.transformed = False
        if self.coords.shape[0] != 2:
            raise ValueError("Something wrong with the coordinate shape")
        self.shape = (nn, nn)
        return

    def sample_galaxy(self, gal_obj):
        """Sample the surface density field of a galaxy at the grids
        This function only conduct sampling; PSF and pixel response are
        not included.

        Args:
            gal_obj (galsim):   Galsim galaxy object to sample on the grids
        Returns:
            outcome (ndarray):  2D galaxy image on the grids
        """
        
        #pixel_values = np.array([gal_obj.xValue(cc + self.centering) for cc in self.coords.T])
        coords = (self.coords.T + self.centering)
        x_coords = coords[:,0].tolist()
        y_coords = coords[:,1].tolist()
        pixel_values = _gsinterface.getFluxVec(gal_obj._sbp, x_coords, y_coords)
        self.flux_array = np.reshape(pixel_values, self.shape)
        return self.flux_array

    def transform_grids(self, transform_obj):
        if not hasattr(transform_obj, "transform"):
            raise TypeError("transform_obj is not in correct data type")
        self.coords = transform_obj.transform(self.coords)
        self.transformed = True
        return
    
def interp_and_convolve(gal_arr, psf_obj, scale):
    """
    Interpolates the given galaxy array as a galsim image, 
    and then convolves it with a given PSF object.

    Parameters:
    - gal_arr (numpy.ndarray): The galaxy flux array.
    - psf_obj (galsim.GSObject): The Point Spread Function (PSF) object.
    - scale (float): The scale of the image.

    Returns:
    - convolved (galsim.GSObject): The convolved object.

    """

    # draw an image from the batsim flux array
    gal_im = galsim.Image(gal_arr, scale=scale)
    
    # intepolate the batsim image as a galsim object
    gal_obj = galsim.InterpolatedImage(gal_im, scale=scale,  normalization='sb', depixelize=False)
    
    # convolve interpolated object with PSF 
    convolved = galsim.Convolve([gal_obj, psf_obj])
    
    return convolved
