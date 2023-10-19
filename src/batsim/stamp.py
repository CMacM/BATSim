import galsim
import numpy as np

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
<<<<<<< HEAD:src/batsim/stamp.py
=======
<<<<<<< HEAD
>>>>>>> main:batsim/stamp.py
            self.centering = (0.5*scale)
        elif centering == 'fpfs':
            self.centering = 0
        else:
            self.centering = centering
        
        if coords is None:
            indx = (np.arange(-int(nn / 2), 
<<<<<<< HEAD:src/batsim/stamp.py
                             int((nn + 1) / 2), 1) * scale)
            indy = (np.arange(-int(nn / 2), 
                             int((nn + 1) / 2), 1) * scale)
=======
                             int((nn + 1) / 2), 1) * scale) + self.centering
            indy = (np.arange(-int(nn / 2), 
                             int((nn + 1) / 2), 1) * scale) + self.centering
=======
            self.centering = (0.5 * scale)
        else:
            print("Centering type not implemented yet",
                "please shift the galsim object directly",
            )
            self.centering = 0

        if coords is None:
            indx = np.arange(-int(nn / 2),
                             int((nn + 1) / 2), 1) * scale
            indy = np.arange(-int(nn / 2),
                             int((nn + 1) / 2), 1) * scale
>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
>>>>>>> main:batsim/stamp.py
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
<<<<<<< HEAD:src/batsim/stamp.py
        
        pixel_values = np.array([gal_obj.xValue(cc + self.centering) for cc in self.coords.T])
=======
<<<<<<< HEAD
        
        pixel_values = np.array([gal_obj.xValue(cc) for cc in self.coords.T])
=======
        pixel_values = np.array([gal_obj.xValue(cc + self.centering)
                                                 for cc in self.coords.T])
>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
>>>>>>> main:batsim/stamp.py

        return np.reshape(pixel_values, self.shape)

    def transform_grids(self, transform_obj):
        if not hasattr(transform_obj, "transform"):
            raise TypeError("transform_obj is not in correct data type")
        self.coords = transform_obj.transform(self.coords)
        self.transformed = True
        return
<<<<<<< HEAD:src/batsim/stamp.py
=======
<<<<<<< HEAD
        
        
        
    

    
        
    


        
=======
>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
>>>>>>> main:batsim/stamp.py
