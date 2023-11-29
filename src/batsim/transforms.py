import numpy as np
import galsim
import warnings

class FlexionTransform(object):

    """
        Feel Free to merge this method with the next one. I wanted
        to be a little more careful because I don't know how to manuever 
        with centers.
    """

    def __init__(self, gamma1, gamma2, kappa, F1=0, F2=0, G1=0, G2=0):
        """Initialize the transform object of 2D grids.
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
        """
            Transform the center of pixels from lensed plane to 
            pre-lensed plane.
            Args:
                coords: coordinates (x, y) of the pixel centers [arcsec]
        """
        return self.s2l_mat @ coords + np.einsum('ijk,jl,kl->il',self.D,coords,coords)

    def inverse_transform(self, coords):
        """
            Details about this inverse transformation can be found 
            here:
            https://github.com/garyang3/Notes/blob/main/Flexion_inverse_transform.pdf
        """
        theta_0 = np.einsum('ij,jk',self.s2l_mat_inv,coords)
        theta_1 = -1/2*np.einsum('in,ijk,jl,lo,km,mo->no',
                                 self.s2l_mat_inv,self.D,
                                 self.s2l_mat_inv,coords,
                                 self.s2l_mat_inv,coords
                                )
        return theta_0 + theta_1


class IaTransform(object):
    """
        Class to apply IA shear transform to a galsim image
        as a function of distance from the center of a galaxy.
    """
    def __init__(self, scale, hlr, A=0.00136207, phi=0,
                 beta=0.82404653, center=None):
        """
            Args:
                scale : The scale of the pixels in arcsec (float)
                hlr : The half light radius of the galaxy to transform
                A : Intrinsic alignment amplitude, this is the shear
                    applied at the half light radius in the distortion
                    definition of shear. Defaults to best fit of Georgiou+19
                    (float)
                beta : Index of the power law used to scale the alignment
                       strength with radius. Defaults to best fit of 
                       Georgiou+19 (float)
                phi : Angle in rads at which alignment should occur. Defaults to
                      zero, which is equivalent to alignment along the horizontal
                      axis. (float)
                center : Coordinates which define the image center from which
                         radius is calculated. (Lists | Tuple | Array)
        """
        # If a center has not been provided, default to [0,0]
        if center == None:
            center = [0,0]
        
        self.ref_vec = np.array([[center[0]],[center[1]]])
        
        # intialise important class variables
        self.A = A
        self.phi = phi
        self.beta = beta
        self.scale = scale
        self.hlr = hlr
        self.xcen = center[0]
        self.ycen = center[1]

        return

    def transform(self, coords):
        """
            Transforms each coordinate with a different shear 
            value depending on its distance from the center 
            of the image.
        """
        npix = np.sqrt(len(coords[0]))
        size_ratio = (npix/2 * self.scale) / self.hlr
        
        if size_ratio < 2.5:
            warnings.simplefilter("always")
            warning_message = ("The stamp provided is only %1.2f"
                               " times larger than the galaxy. To ensure" 
                               " accurate results, the stamp needs to be at"
                               " least 2.5 times larger.")%size_ratio
            warnings.warn(warning_message)
            
        # unpack x and y coordinates
        coords_relative = coords - self.ref_vec
        x, y = coords_relative
        
        g1, g2 = self.get_g1g2(x,y) 

        # transform coordinates with raidal dependence
        x_prime = ((1 - g1)*x - g2*y)
        y_prime = ((1 + g1)*y - g2*x)

        return np.array([x_prime, y_prime]) + self.ref_vec

    def get_g1g2(self,x,y):
        """
            Scales the amplitude according to power law, then
            gets the g1 and g2 components to construct the shear
            matrix.
        """

        # find distance from image center as ratio to hlr
        radial_dist = np.sqrt(abs(x - self.xcen)**2 + abs(y - self.ycen)**2)
        rwf = (radial_dist) / self.hlr

        # compute alignment amplitude at radius
        A_rwf = self.A * rwf**self.beta

        # get shear components for corresponding alignment amplitude
        e1 = A_rwf * np.cos(2*self.phi)
        e2 = A_rwf * np.sin(2*self.phi)
        
        # convert to g in same way galsim does
        absesq = e1**2 + e2**2

        # test to make sure level of distortion is reasonable
        if type(absesq) == np.float64 and absesq > 1:
            raise ValueError("Requested distortion exceeds 1.",
                                   np.sqrt(absesq), 0., 1.)
        elif type(absesq) ==  np.ndarray and any(absesq) > 1:
            raise ValueError("Requested distortion exceeds 1.",
                                   np.sqrt(absesq), 0., 1.)
            
        # define g from e1 and e2
        g = (e1 + 1j * e2) * self.e2g(absesq)
        
        # return real (g1) and imaginary (g2) components
        return g.real, g.imag
    
    # conversion used in galsim source code
    # modified to use binary arrays to speed up condition checking
    def e2g(self, absesq):
        if type(absesq) == np.ndarray:
            # this section of code is designed to avoid the use of a for loop
            # to maximise speed. Stable is a vector which will == 1 if a value
            # in absesq is big enough to use the simple calculation, and a
            # 0 if the Taylor expansion is needed for stability.
            stable = 2*np.ones(len(absesq))
            drill = 1e-10
            # if there are values greater than 1, continue with a deeper drill
            while stable.max() > 1:
                stable = np.floor((absesq/1.e-4)**(drill)).astype(int)
                drill = drill/10 
            # unstable values set to zero, stable values included    
            e2g = stable * (1. / (1. + np.sqrt(1.-absesq)))
            # now we invert to have unstable values as 1's
            unstable = abs(stable - 1).astype(int)
            # we add the unstable values to the array now, with the stable set to zero
            e2g += unstable * (0.5 + absesq*(0.125 
                                             + absesq*(0.0625 
                                                       + absesq*0.0390625)))
            # finally, return the array of conversion values
            return e2g
        # for if we just want a single shear value
        elif type(absesq) == np.float64:
            if absesq > 1.e-4:
                #return (1. - np.sqrt(1.-absesq)) / absesq
                return 1. / (1. + np.sqrt(1.-absesq))
            else:
                # Avoid numerical issues near e=0 using Taylor expansion
                return 0.5 + absesq*(0.125 + absesq*(0.0625 + absesq*0.0390625))

class LensTransform(object):
    def __init__(self, gamma1, gamma2, kappa, center=None):
        """
            Initialize the transform object of 2D grids.
            Args:
                gamma1 (float):   the first component of lensing shear field
                gamma2 (float):   the second component of lensing shear field
                kappa (float):    the lensing convergence field
                xref (float):     reference coordinate x [in units of pixels]
                xref (float):     reference coordinate y [in units of pixels]
        """
        
        if center==None:
            center=[0,0]
        
        self.ref_vec = np.array([[center[0]],[center[1]]])
        self.s2l_mat = np.array(
            [[1 - kappa - gamma1, -gamma2], [-gamma2, 1 - kappa + gamma1]]
        )
        return

    def transform(self, coords):
        """
            Transform the center of pixels from lensed plane to 
            pre-lensed plane.
            Args:
                coords:   coordinates (x, y) of the pixel centers [arcsec]
        """
        coords_relative = coords - self.ref_vec
        return self.s2l_mat @ coords_relative + self.ref_vec
