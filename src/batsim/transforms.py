import warnings

import numpy as np


_ARRAY_BACKEND = None
_CPU_FALLBACK_WARNED = False


def _coords_backend(coords):
    """Return the array backend for a coordinate array."""
    module = type(coords).__module__.split(".", 1)[0]
    if module == "cupy":
        import cupy as cp

        return cp

    return np


def _coords_dtype(coords):
    """Return a floating dtype matching coords."""
    return getattr(coords, "dtype", None)


def _get_array_backend():
    """
    Return CuPy if available, otherwise NumPy.

    This mirrors the backend-selection behaviour used in the simulation code.
    """
    global _ARRAY_BACKEND
    global _CPU_FALLBACK_WARNED

    if _ARRAY_BACKEND is not None:
        return _ARRAY_BACKEND

    try:
        import cupy as cp

        cp.cuda.runtime.getDeviceCount()
        _ARRAY_BACKEND = cp
    except Exception:
        if not _CPU_FALLBACK_WARNED:
            warnings.warn(
                "CuPy unavailable; falling back to NumPy stamp coordinates.",
                RuntimeWarning,
                stacklevel=2,
            )
            _CPU_FALLBACK_WARNED = True
        _ARRAY_BACKEND = np

    return _ARRAY_BACKEND



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
        D1 = -1 / 2 * np.array([[3 * F1 + G1, F2 + G2], [F2 + G2, F1 - G1]])
        D2 = -1 / 2 * np.array([[F2 + G2, F1 - G1], [F1 - G1, 3 * F2 - G2]])
        self.D = np.stack([D1, D2], axis=2)
        return

    def transform(self, coords):
        """
        Transform the center of pixels from lensed plane to
        pre-lensed plane.

        Args:
        coords: coordinates (x, y) of the pixel centers [arcsec]
        """
        xp = _coords_backend(coords)
        dtype = _coords_dtype(coords)
        coords = xp.asarray(coords)
        s2l_mat = xp.asarray(self.s2l_mat, dtype=dtype)
        d_tensor = xp.asarray(self.D, dtype=dtype)

        return s2l_mat @ coords + xp.einsum(
            "ijk,jl,kl->il", d_tensor, coords, coords
        )

    def inverse_transform(self, coords):
        """
        Details about this inverse transformation can be found
        here:
        https://github.com/garyang3/Notes/blob/main/Flexion_inverse_transform.pdf
        """
        xp = _coords_backend(coords)
        dtype = _coords_dtype(coords)
        coords = xp.asarray(coords)
        s2l_mat_inv = xp.asarray(self.s2l_mat_inv, dtype=dtype)
        d_tensor = xp.asarray(self.D, dtype=dtype)

        theta_0 = xp.einsum("ij,jk", s2l_mat_inv, coords)
        theta_1 = (
            -1
            / 2
            * xp.einsum(
                "in,ijk,jl,lo,km,mo->no",
                s2l_mat_inv,
                d_tensor,
                s2l_mat_inv,
                coords,
                s2l_mat_inv,
                coords,
            )
        )
        return theta_0 + theta_1


class IaTransform(object):
    """
    Class to apply IA shear transform to a galsim image
    as a function of distance from the center of a galaxy.
    """

    def __init__(
        self,
        scale,
        hlr,
        A=0.00136207,
        phi=0,
        beta=0.82404653,
        center=None,
        clip_radius=5,
    ):
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
            center = [0, 0]

        self.ref_vec = np.array([[center[0]], [center[1]]])

        # intialise important class variables
        self.A = A
        self.phi = phi
        self.c2phi = np.cos(2 * self.phi)
        self.s2phi = np.sin(2 * self.phi)
        self.beta = beta
        self.scale = scale
        self.hlr = hlr
        self.xcen = center[0]
        self.ycen = center[1]
        self.clip_radius = clip_radius

        return

    def transform(self, coords):
        """
        Transforms each coordinate with a different shear
        value depending on its distance from the center
        of the image.
        """
        xp = _coords_backend(coords)
        dtype = _coords_dtype(coords)
        coords = xp.asarray(coords)
        ref_vec = xp.asarray(self.ref_vec, dtype=dtype)

        npix = xp.sqrt(len(coords[0]))

        # if size_ratio < 2.5:
        #     warnings.simplefilter("always")
        #     warning_message = ("The stamp provided is only %1.2f"
        #                        " times larger than the galaxy. To ensure"
        #                        " accurate results, the stamp needs to be at"
        #                        " least 2.5 times larger.")%size_ratio
        #     warnings.warn(warning_message)

        # unpack x and y coordinates
        coords_relative = coords - ref_vec

        x, y = coords_relative

        g1, g2 = self.get_g1g2(x, y)

        # transform coordinates with raidal dependence
        x_prime = (1 - g1) * x - g2 * y
        y_prime = (1 + g1) * y - g2 * x

        coords_realtive_transformed = xp.stack([x_prime, y_prime], axis=0)

        return coords_realtive_transformed + ref_vec

    def get_g1g2(self, x, y):
        """
        Scales the amplitude according to power law, then
        gets the g1 and g2 components to construct the shear
        matrix.
        """
        xp = _coords_backend(x)
        dtype = _coords_dtype(x)
        hlr = xp.asarray(self.hlr, dtype=dtype)
        amp = xp.asarray(self.A, dtype=dtype)
        c2phi = xp.asarray(self.c2phi, dtype=dtype)
        s2phi = xp.asarray(self.s2phi, dtype=dtype)

        # find distance from image center as ratio to hlr
        radial_dist = xp.sqrt(abs(x) ** 2 + abs(y) ** 2)
        rwf = radial_dist / hlr

        # fix shear beyond rfw >= clip_radius
        rwf = xp.clip(rwf, 0, self.clip_radius)

        # compute alignment amplitude at radius
        A_rwf = amp * rwf**self.beta
        absesq = A_rwf * A_rwf

        if bool(xp.any(absesq > 1)):
            raise ValueError(
                "Requested distortion exceeds 1.", xp.sqrt(absesq), 0.0, 1.0
                        )

        # factor to convert e1, e2 to g1, g2
        fac = self.e2g(absesq)

        g1 = A_rwf * c2phi * fac
        g2 = A_rwf * s2phi * fac

        # return real (g1) and imaginary (g2) components
        return g1, g2

    # conversion used in galsim source code
    # modified to use binary arrays to speed up condition checking
    def e2g(self, absesq):
        xp = _coords_backend(absesq)

        if hasattr(absesq, "shape"):
            # if absesq is big enough to use the simple calculation, and a
            # 0 if the Taylor expansion is needed for stability.
            # if there are values greater than 1, continue with a deeper drill
            stable = absesq > 1e-4
            # unstable values set to False, stable values included
            e2g = stable * (1.0 / (1.0 + xp.sqrt(1.0 - absesq)))
            # now we invert to have unstable values as True
            unstable = absesq <= 1e-4
            # we add the unstable values to the array now, with the stable set to zero
            e2g += unstable * (
                0.5 + absesq * (0.125 + absesq * (0.0625 + absesq * 0.0390625))
            )
            # finally, return the array of conversion values
            return e2g
        # for if we just want a single shear value
        elif type(absesq) == np.float64:
            if absesq > 1.0e-4:
                # return (1. - np.sqrt(1.-absesq)) / absesq
                return 1.0 / (1.0 + np.sqrt(1.0 - absesq))
            else:
                # Avoid numerical issues near e=0 using Taylor expansion
                return 0.5 + absesq * (0.125 + absesq * (0.0625 + absesq * 0.0390625))


class LensTransform:
    """
    Affine lensing transform supporting either NumPy or CuPy arrays.

    If CuPy is available and the input coordinates are CuPy arrays, the
    transform runs natively on GPU.
    """

    def __init__(
        self,
        gamma1,
        gamma2,
        kappa,
        center=None,
        backend=None,
        dtype=None,
    ):
        """
        Parameters
        ----------
        gamma1, gamma2 : float
            Components of lensing shear.
        kappa : float
            Lensing convergence.
        center : sequence of float, optional
            Reference coordinate [x, y].
        backend : module, optional
            NumPy or CuPy. If None, NumPy is used. Simulation code moves
            transforms to the requested backend before applying them.
        dtype : dtype, optional
            Floating dtype for transform arrays.
        """
        self.xp = np if backend is None else backend

        if dtype is None:
            dtype = self.xp.float64

        if center is None:
            center = [0.0, 0.0]

        self.dtype = dtype

        self.ref_vec = self.xp.asarray(
            [[center[0]], [center[1]]],
            dtype=self.dtype,
        )

        self.s2l_mat = self.xp.asarray(
            [
                [1.0 - kappa - gamma1, -gamma2],
                [-gamma2, 1.0 - kappa + gamma1],
            ],
            dtype=self.dtype,
        )

    def transform(self, coords):
        """
        Transform pixel-centre coordinates from lensed plane to pre-lensed plane.

        Parameters
        ----------
        coords : array-like
            Coordinate array with shape (2, npoints). Can be NumPy or CuPy.

        Returns
        -------
        array-like
            Transformed coordinates using the configured backend.
        """
        coords = self.xp.asarray(coords)
        dtype = _coords_dtype(coords)

        ref_vec = self.xp.asarray(self.ref_vec, dtype=dtype)
        s2l_mat = self.xp.asarray(self.s2l_mat, dtype=dtype)

        coords_relative = coords - ref_vec
        return s2l_mat @ coords_relative + ref_vec

    def to_backend(self, backend=None, dtype=None):
        """
        Return a copy of this transform on another backend/dtype.
        """
        xp = np if backend is None else backend

        if dtype is None:
            dtype = self.dtype

        ref_vec = self._to_numpy(self.ref_vec)
        s2l_mat = self._to_numpy(self.s2l_mat)

        new = object.__new__(LensTransform)
        new.xp = xp
        new.dtype = dtype
        new.ref_vec = xp.asarray(ref_vec, dtype=dtype)
        new.s2l_mat = xp.asarray(s2l_mat, dtype=dtype)

        return new

    @staticmethod
    def _to_numpy(array):
        """Convert NumPy/CuPy array to NumPy."""
        if isinstance(array, np.ndarray):
            return array

        return array.get()
