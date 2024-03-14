import warnings

import galsim
import numpy as np
from scipy.fft import fft2, fftshift, ifft2, ifftshift

from . import _gsinterface


class Stamp(object):
    def __init__(self, coords=None, nn=32, scale=0.2, centering="fpfs"):
        """Initialize the 2D stamp object. This class enables distorting
        an image by changing the samplinng position with non-affine
        transformation

        Args:
            nn (int):      number of grids on x and y direction
            scale (float): pixel scale in units of arcsec
        """
        self.scale = scale
        if centering == "galsim":
            self.centering = 0.5 * scale
        elif centering == "fpfs":
            self.centering = 0
        else:
            self.centering = centering

        # Set up the grids
        if coords is None:
            self.coords = self.get_coords(nn, scale)
        else:
            self.coords = coords

        if self.coords.shape[0] != 2:
            raise ValueError("Something wrong with the coordinate shape")
        self.shape = (nn, nn)

        self.pixel_values = None
        self.convolved = False
        self.flux_array = None
        self.transform_obj = None
        return

    def get_coords(self, nn, scale):
        indx = np.arange(-int(nn / 2), int((nn + 1) / 2), 1) * scale
        indy = np.arange(-int(nn / 2), int((nn + 1) / 2), 1) * scale
        inds = np.meshgrid(indy, indx, indexing="ij")
        coords = np.vstack([np.ravel(_) for _ in inds[::-1]])
        return coords

    def sample_galaxy(self, gal_obj, normalise=False):
        """Sample the surface density field of a galaxy at the grids
        This function only conduct sampling; PSF and pixel response are
        not included.

        Args:
            gal_obj (galsim):   Galsim galaxy object to sample on the grids
        Returns:
            outcome (ndarray):  2D galaxy image on the grids
        """

        # pixel_values = np.array([gal_obj.xValue(cc + self.centering) for cc in self.coords.T])
        if self.transform_obj is not None:
            coords = self.transform_obj.transform(self.coords)
            coords = coords.T + self.centering
        else:
            coords = self.coords.T + self.centering
        x_coords = coords[:, 0].tolist()
        y_coords = coords[:, 1].tolist()
        pixel_values = _gsinterface.getFluxVec(gal_obj._sbp, x_coords, y_coords)
        self.flux_array = np.reshape(pixel_values, self.shape)

        if normalise:
            self.flux_array /= np.sum(self.flux_array) / gal_obj.flux

        return self.flux_array

    def set_transform(self, transform_obj):
        if not hasattr(transform_obj, "transform"):
            raise TypeError("transform_obj is not in correct data type")
        self.transform_obj = transform_obj
        return

    def sample_convolution(self, gal_obj, psf_obj, normalise=True):
        # Get required size and scale for convolution to be accurate
        gal_scale = gal_obj.nyquist_scale
        gal_nn = gal_obj.getGoodImageSize(gal_scale)
        psf_nn = psf_obj.getGoodImageSize(gal_scale)

        # Get the galaxy coordinates and apply the transform if one is present
        gal_coords = self.get_coords(gal_nn, gal_scale)
        if self.transform_obj is not None:
            # If transform depends on pixel scale, update it
            local_transform = self.transform_obj
            if hasattr(self.transform_obj, "scale"):
                local_transform.scale = gal_scale

            gal_coords = local_transform.transform(gal_coords)

        gal_coords = gal_coords.T + self.centering
        gal_x = gal_coords[:, 0].tolist()
        gal_y = gal_coords[:, 1].tolist()
        galprof = _gsinterface.getFluxVec(gal_obj._sbp, gal_x, gal_y)
        galprof = np.reshape(galprof, (gal_nn, gal_nn))

        # Get the PSF coordinates
        psf_coords = self.get_coords(psf_nn, gal_scale).T + self.centering
        psf_x = psf_coords[:, 0].tolist()
        psf_y = psf_coords[:, 1].tolist()
        psfprof = _gsinterface.getFluxVec(psf_obj._sbp, psf_x, psf_y)
        psfprof = np.reshape(psfprof, (psf_nn, psf_nn))

        # Normalise the resulting images
        if normalise:
            galprof /= np.sum(galprof) / gal_obj.flux
            psfprof /= np.sum(psfprof) / psf_obj.flux

        # Pad smaller image with zeros to match the larger image
        if gal_nn > psf_nn:
            pad_width = ((gal_nn - psf_nn) // 2,)
            psfprof = np.pad(psfprof, pad_width, mode="constant", constant_values=0)
        elif psf_nn > gal_nn:
            pad_width = ((psf_nn - gal_nn) // 2,)
            galprof = np.pad(galprof, pad_width, mode="constant", constant_values=0)
        else:
            pass

        # Fourier transform the galaxy and PSF profiles
        galKprof = fft2(galprof)
        psfKprof = fft2(psfprof)

        # Convolve the galaxy and PSF profiles in Fourier space
        convKprof = galKprof * psfKprof

        # Get factor by which to rescale convolved image
        rescale = self.scale / gal_scale
        nn_cut = int(gal_nn / rescale)
        maxN = gal_nn // 2 + nn_cut // 2
        minN = gal_nn // 2 - nn_cut // 2
        convKprof = convKprof[minN:maxN, minN:maxN]

        # Inverse Fourier transform to obtain convolution result in real space
        conv_im = ifftshift(ifft2(convKprof))
        conv_im = conv_im.real  # Throw away imaginary part
        conv_nn = np.shape(conv_im)[0]

        # crop or pad the realspace imaghe to the original size
        if conv_nn > self.shape[0]:
            max_nn = conv_nn // 2 + self.shape[0] // 2
            min_nn = conv_nn // 2 - self.shape[0] // 2
            conv_im = conv_im[min_nn:max_nn, min_nn:max_nn]
        elif conv_nn < self.shape[0]:
            pad_width = ((self.shape[0] - conv_nn) // 2,)
            conv_im = np.pad(conv_im, pad_width, "constant", constant_values=0)

        # Update stored flux array and properties
        self.flux_array = conv_im
        self.convolved = True

        return self.flux_array

    def transform_grids(self, transform_obj):
        if not hasattr(transform_obj, "transform"):
            raise TypeError("transform_obj is not in correct data type")
        self.coords = transform_obj.transform(self.coords)
        self.transformed = True
        return
