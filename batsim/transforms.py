import numpy as np 
import galsim

class AffineLensTransform(object):
    def __init__(self, g1, g2, kappa):
        """Initialize the transform object of 2D grids
        Args:
            gamma1 (float):     the first component of lensing shear field
            gamma2 (float):     the second component of lensing shear field
            kappa (float):      the lensing convergence field
        """
        
        self.s2l_mat = np.array([[1 - kappa - g1, -g2], 
                                 [-g2, 1 - kappa + g1]])
        
        return
    
    def transform(self, coords):
        return self.s2l_mat @ coords

class IaTransform(object):
    def __init__(self, scale, hlr, A1=-0.00136207, A2=0, 
                 beta=0.82404653, center=[0,0]):
        """
            Class to apply IA shear transform to a galsim image
            as a function of distance from the center of a galaxy.
            
            Args:
            g1 : The g1 component of shear (float)
            g2 : The g2 component of shear (float)
            center : Coordinates which 
                     define the image center (Lists | Tuple | Array)
            scale : The scale of the pixels in arcsec (float)
            hlr : The half light radius of the galaxy to transform
        """

        self.A1 = A1
        self.A2 = A2
        self.beta = beta
        self.scale = scale
        self.hlr = hlr
        self.xcen = center[0] * scale
        self.ycen = center[1] * scale
         
        return

    def transform(self, coords):
        """
            Transform each coordinate with a different
            shear value depending on its distance from the 
            center of the image.
        """
        x, y = coords
        
        x_prime = ((1 - self.get_g1(x,y))*x - self.get_g2(x,y)*y)
        y_prime = ((1 + self.get_g1(x,y))*y - self.get_g2(x,y)*x)
        
        return np.array([x_prime, y_prime])
        
    def get_g1(self,x,y):
        """
            Get g1 term for a set of image coordinates.
        """
        radial_dist = np.sqrt(abs(x - self.xcen)**2 + abs(y - self.ycen)**2)
        rwf = (radial_dist) / self.hlr

        return self.A1 * rwf**self.beta
    
    def get_g2(self,x,y):
        """
            Get g2 term for set of image coordinates.
        """
        radial_dist = np.sqrt(abs(x - self.xcen)**2 + abs(y - self.ycen)**2)
        rwf = (radial_dist) / self.hlr

        return self.A2 * rwf**self.beta