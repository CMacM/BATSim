import galsim
import numpy as np

class Stamp(object):
    def __init__(self, coords=None, nn=32, scale=0.2):
        """Initialize the 2D stamp object. This class enables distorting
        an image by changing the samplinng position with non-affine
        transformation

        Args:
            nn (int):      number of grids on x and y direction
            scale (float): pixel scale in units of arcsec
        """
        if coords is None:
            indx = np.arange(-int(nn / 2), int((nn + 1) / 2), 1) * scale
            indy = np.arange(-int(nn / 2), int((nn + 1) / 2), 1) * scale
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
        pixel_values = np.array([gal_obj.xValue(cc) for cc in self.coords.T])

        return np.reshape(pixel_values, self.shape)

    def transform_grids(self, transform_obj):
        if not hasattr(transform_obj, "transform"):
            raise TypeError("transform_obj is not in correct data type")
        self.coords = transform_obj.transform(self.coords)
        self.transformed = True
        return


class LensTransform1(object):
    def __init__(self, gamma1, gamma2, kappa, non_affine=False,
                scale=None, hlr=None):
        """Initialize the transform object of 2D grids
        Args:
            gamma1 (float):     the first component of lensing shear field
            gamma2 (float):     the second component of lensing shear field
            kappa (float):      the lensing convergence field
        """
        
        self.non_affine = non_affine
        self.scale = scale
        self.hlr = hlr
        self.gamma1 = gamma1
        self.gamma2 = gamma2
        self.kappa = kappa
        
        return

    def transform(self, coords):
        """transform the center of pixels from lensed plane to pre-lensed plane
        Args:
            coords:   coordinates (x, y) of the pixel centers [arcsec]
        """
    
        if not self.non_affine:
            
            s2l_mat = np.array([[1 - self.kappa - self.gamma1, -self.gamma2], 
                                [-self.gamma2, 1 - self.kappa + self.gamma1]])
            
            return s2l_mat @ coords
        
        ### TO-DO: currently stuff is hardcoded since it was a quick mock-up ####
            ### Everything below here will need tidied up ### 
        x, y = coords
        
        x_prime = (1 - self.gamma1_coord(x,y))*x - self.gamma2_coord(x,y)*y
        y_prime = (1 + self.gamma1_coord(x,y))*y - self.gamma2_coord(x,y)*x
        
        return np.array([x_prime, y_prime])
        
    def gamma1_coord(self,x,y):
        
        radial_dist = np.sqrt((x - 33)**2 + (y - 33)**2)
        rwf = (radial_dist * self.scale) / self.hlr #hlr
        
        return self.gamma1 * rwf**1
    
    def gamma2_coord(self,x,y):
        
        radial_dist = np.sqrt((x - 33)**2 + (y - 33)**2)
        rwf = (radial_dist * self.scale) / self.hlr #hlr
        
        return self.gamma2 * rwf**1
        
        
        
    

    
        
    


        