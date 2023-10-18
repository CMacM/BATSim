import numpy as np 
import galsim

class FlexionTransform(object):
    
    """Feel Free to merge this method with the next one. I wanted
    to be a little more careful because I don't know how to manuever with centers."""
    def __init__(self, gamma1, gamma2, kappa, F1=0, F2=0, G1=0, G2=0):
        """Initialize the transform object of 2D grids
        Args:
            gamma1 (float):     the first component of lensing shear field
            gamma2 (float):     the second component of lensing shear field
            kappa (float):      the lensing convergence field
            F1,F2,G1,G2 (float):  Flexion components
        """
        self.s2l_mat = np.array(
            [[1 - kappa - gamma1, -gamma2], [-gamma2, 1 - kappa + gamma1]]
        )
        self.s2l_mat_inv = np.linalg.inv(self.s2l_mat)
        D1 = -1/2*np.array([[3*F1+G1,F2+G2],[F2+G2,F1-G1]])
        D2 = -1/2*np.array([[F2+G2,F1-G1],[F1-G1,3*F2-G2]])
        self.D = np.stack([D1,D2],axis=2)
        return

    def transform(self, coords):
        """transform the center of pixels from lensed plane to pre-lensed plane
        Args:
            coords:   coordinates (x, y) of the pixel centers [arcsec]
        """
        return self.s2l_mat @ coords + np.einsum('ijk,jl,kl->il',self.D,coords,coords)    
    
    def inverse_transform(self, coords):
        """Details about this inverse transformation can be found here:
        https://github.com/garyang3/Notes/blob/main/Flexion_inverse_transform.pdf"""
        theta_0 = np.einsum('ij,jk',self.s2l_mat_inv,coords)
        theta_1 = -1/2*np.einsum('in,ijk,jl,lo,km,mo->no',self.s2l_mat_inv,self.D,self.s2l_mat_inv,coords,self.s2l_mat_inv,coords)
        return theta_0 + theta_1
    
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