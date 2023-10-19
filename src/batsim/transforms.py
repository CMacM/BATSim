import numpy as np
import galsim

<<<<<<< HEAD:src/batsim/transforms.py
=======
<<<<<<< HEAD
class AffineLensTransform(object):
    def __init__(self, g1, g2, kappa):
=======
>>>>>>> main:batsim/transforms.py
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
    
<<<<<<< HEAD:src/batsim/transforms.py
=======
class IaTransform(object):
    def __init__(self, g1, g2, center, scale, hlr):
>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
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
>>>>>>> main:batsim/transforms.py

class IaTransform(object):
    def __init__(self, scale, hlr, A=0.00136207, phi=0, 
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
        
        # intialise important class variables
        self.A = A
        self.phi = phi
        self.beta = beta
        self.scale = scale
        self.hlr = hlr
<<<<<<< HEAD:src/batsim/transforms.py
        self.xcen = center[0]
        self.ycen = center[1]
        
=======
<<<<<<< HEAD
        self.xcen = center[0] * scale
        self.ycen = center[1] * scale
         
=======
        self.xcen = center[0]
        self.ycen = center[1]

>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
>>>>>>> main:batsim/transforms.py
        return

    def transform(self, coords):
        """
<<<<<<< HEAD:src/batsim/transforms.py
=======
<<<<<<< HEAD
>>>>>>> main:batsim/transforms.py
            Transform each coordinate with a different
            shear value depending on its distance from the 
            center of the image.
        """
<<<<<<< HEAD:src/batsim/transforms.py
        
        # unpack x and y coordinates
        x, y = coords
        
        # transform coordinates with raidal dependence
        x_prime = ((1 - self.get_e1(x,y))*x - self.get_e2(x,y)*y)
        y_prime = ((1 + self.get_e1(x,y))*y - self.get_e2(x,y)*x)
        
        return np.array([x_prime, y_prime])

    def get_e1(self,x,y):
        """
            Get e1 term for a set of image coordinates.
        """
        
        # find distance from image center as ratio to hlr
=======
        x, y = coords
        
        x_prime = ((1 - self.get_g1(x,y))*x - self.get_g2(x,y)*y)
        y_prime = ((1 + self.get_g1(x,y))*y - self.get_g2(x,y)*x)
        
=======

        ### TO-DO: currently stuff is hardcoded since it was a quick mock-up ####
            ### Everything below here will need tidied up ###
        x, y = coords

        norm = np.sqrt(1 - self.get_g1(x,y)**2 - self.get_g2(x,y)**2)

        x_prime = norm * ((1 - self.get_g1(x,y))*x - self.get_g2(x,y)*y)
        y_prime = norm * ((1 + self.get_g1(x,y))*y - self.get_g2(x,y)*x)


>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
        return np.array([x_prime, y_prime])

    def get_g1(self,x,y):
<<<<<<< HEAD
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
=======

>>>>>>> main:batsim/transforms.py
        radial_dist = np.sqrt(abs(x - self.xcen)**2 + abs(y - self.ycen)**2)
        rwf = (radial_dist) / self.hlr
        
        # compute alignment amplitude at radius
        A_rwf = self.A * rwf**self.beta
        
        # get shear component for corresponding alignment amplitude
        e1 = A_rwf * np.cos(2*self.phi)
        
        # convert to g1 componet of shear for matrix application
        g1 = np.asarray([galsim.Shear(e1=e).g1 for e in e1])
        
        return g1
    
    def get_e2(self,x,y):
        """
            Get e2 term for set of image coordinates.
        """
        # find distance from image center as ratio to hlr
        radial_dist = np.sqrt(abs(x - self.xcen)**2 + abs(y - self.ycen)**2)
        rwf = (radial_dist) / self.hlr
        
        # compute alignment amplitude at radius
        A_rwf = self.A * rwf**self.beta
        
        # get shear component for corresponding alignment amplitude
        e2 = A_rwf * np.sin(2*self.phi)
        
        # convert to g1 componet of shear for matrix application
        #g2 = np.empty
        g2 = np.asarray([galsim.Shear(e2=e).g2 for e in e2])

        return g2


class LensTransform(object):
    def __init__(self, gamma1, gamma2, kappa, center=[-0.,-0.]):
        """Initialize the transform object of 2D grids
        Args:
            gamma1 (float):     the first component of lensing shear field
            gamma2 (float):     the second component of lensing shear field
            kappa (float):      the lensing convergence field
            xref (float):       reference coordinate x [in units of pixels]
            xref (float):       reference coordinate y [in units of pixels]
        """
        self.ref_vec = np.array([[center[0]],[center[1]]])
        self.s2l_mat = np.array(
            [[1 - kappa - gamma1, -gamma2], [-gamma2, 1 - kappa + gamma1]]
        )
        return

    def transform(self, coords):
        """transform the center of pixels from lensed plane to pre-lensed plane
        Args:
            coords:   coordinates (x, y) of the pixel centers [arcsec]
        """
        coords_relative = coords - self.ref_vec
        return self.s2l_mat @ coords_relative + self.ref_vec
<<<<<<< HEAD:src/batsim/transforms.py
=======
>>>>>>> 7910f724af76c1a1a61e60076c219bc877335c4d
>>>>>>> main:batsim/transforms.py
