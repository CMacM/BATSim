import numpy as np 
import galsim

class IaTransform(object):
    def __init__(self, g1, g2, center, scale, hlr):
        """Initialize the transform object of 2D grids
        Args:
            gamma1 (float):     the first component of lensing shear field
            gamma2 (float):     the second component of lensing shear field
            kappa (float):      the lensing convergence field
        """
        self.g1 = g1
        self.g2 = g2
        self.scale = scale
        self.hlr = hlr
        self.xcen = center[0]
        self.ycen = center[1]
        
        return

    def transform(self, coords):
        """transform the center of pixels from lensed plane to pre-lensed plane
        Args:
            coords:   coordinates (x, y) of the pixel centers [arcsec]
        """
        
        ### TO-DO: currently stuff is hardcoded since it was a quick mock-up ####
            ### Everything below here will need tidied up ### 
        x, y = coords
        
        norm = np.sqrt(1 - self.get_g1(x,y)**2 - self.get_g2(x,y)**2)
        
        x_prime = norm * ((1 - self.get_g1(x,y))*x - self.get_g2(x,y)*y)
        y_prime = norm * ((1 + self.get_g1(x,y))*y - self.get_g2(x,y)*x)
        
        
        return np.array([x_prime, y_prime])
        
    def get_g1(self,x,y):
        
        radial_dist = np.sqrt(abs(x - self.xcen)**2 + abs(y - self.ycen)**2)
        rwf = (radial_dist * self.scale) / self.hlr #hlr
        
        return self.g1 * rwf**1
    
    def get_g2(self,x,y):
        
        radial_dist = np.sqrt((x - self.xcen)**2 + (y - self.ycen)**2)
        rwf = (radial_dist * self.scale) / self.hlr #hlr
        
        return self.g2 * rwf**1