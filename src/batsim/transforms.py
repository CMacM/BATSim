import numpy as np


def _coords_backend(coords):
    """Return the array backend for a coordinate array."""
    module = type(coords).__module__.split(".", 1)[0]
    if module == "cupy":
        import cupy as cp

        return cp

    return np


def _coords_dtype(coords):
    """Return the dtype attached to a coordinate array, if present."""
    return getattr(coords, "dtype", None)


class FlexionTransform(object):
    """
    Coordinate transform for affine lensing with first-order flexion terms.

    The transform maps lensed-plane pixel-centre coordinates into the
    corresponding pre-lensed/source-plane coordinates. Inputs may be NumPy or
    CuPy arrays; the backend is inferred from the coordinate array passed to
    ``transform`` or ``inverse_transform``.

    Parameters
    ----------
    gamma1, gamma2 : float
        Components of the reduced shear field.
    kappa : float
        Lensing convergence.
    F1, F2 : float, optional
        Components of first flexion. Defaults to zero.
    G1, G2 : float, optional
        Components of third flexion. Defaults to zero.
    """

    def __init__(self, gamma1, gamma2, kappa, F1=0, F2=0, G1=0, G2=0):
        """
        Initialize the affine and flexion tensors.

        Parameters
        ----------
        gamma1, gamma2 : float
            Components of the shear field.
        kappa : float
            Lensing convergence.
        F1, F2, G1, G2 : float, optional
            Flexion components. Defaults to zero.
        """
        self.s2l_mat = np.array([[1 - kappa - gamma1, -gamma2], [-gamma2, 1 - kappa + gamma1]])
        self.s2l_mat_inv = np.linalg.inv(self.s2l_mat)
        D1 = -1 / 2 * np.array([[3 * F1 + G1, F2 + G2], [F2 + G2, F1 - G1]])
        D2 = -1 / 2 * np.array([[F2 + G2, F1 - G1], [F1 - G1, 3 * F2 - G2]])
        self.D = np.stack([D1, D2], axis=2)
        return

    def transform(self, coords):
        """
        Transform pixel-centre coordinates from lensed to pre-lensed plane.

        Parameters
        ----------
        coords : array-like
            Coordinate array with shape ``(2, npoints)``. The first row holds
            x coordinates and the second row holds y coordinates, in a
            consistent angular unit such as arcseconds.

        Returns
        -------
        array-like
            Transformed coordinate array using the same backend as ``coords``.
        """
        xp = _coords_backend(coords)
        dtype = _coords_dtype(coords)
        coords = xp.asarray(coords)
        s2l_mat = xp.asarray(self.s2l_mat, dtype=dtype)
        d_tensor = xp.asarray(self.D, dtype=dtype)

        return s2l_mat @ coords + 0.5 * xp.einsum("ijk,jl,kl->il", d_tensor, coords, coords)

    def inverse_transform(self, coords):
        """
        Approximate the inverse flexion transform for coordinate arrays.

        Parameters
        ----------
        coords : array-like
            Coordinate array with shape ``(2, npoints)`` in the transformed
            plane.

        Returns
        -------
        array-like
            Approximate inverse-transformed coordinates using the same backend
            as ``coords``.

        Notes
        -----
        The approximation follows the derivation in
        https://github.com/garyang3/Notes/blob/main/Flexion_inverse_transform.pdf.
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
    Radius-dependent intrinsic-alignment shear transform.

    The transform applies a reduced-shear coordinate mapping whose amplitude is
    set by a power law in radius from ``center``. The radial coordinate is
    measured in units of ``hlr`` and clipped at ``clip_radius`` before the
    amplitude law is evaluated.

    Parameters
    ----------
    scale : float
        Pixel scale in arcseconds. Stored for compatibility with existing
        callers.
    hlr : float
        Half-light radius used to normalize the radial profile.
    A : float, optional
        Distortion amplitude at one half-light radius. Defaults to the
        Georgiou+19 best-fit value used by the existing pipeline.
    phi : float, optional
        Alignment angle in radians. A value of zero aligns the distortion with
        the horizontal axis.
    beta : float, optional
        Power-law index for the radial amplitude profile.
    center : sequence of float, optional
        Coordinate origin ``[x, y]`` for the radial profile. Defaults to
        ``[0, 0]``.
    clip_radius : float, optional
        Maximum radius, in units of ``hlr``, used when evaluating the amplitude
        profile.
    """

    def __init__(
        self,
        scale,
        hlr,
        A=3.54e-4,
        phi=0,
        beta=1.11,
        center=None,
        clip_radius=5,
    ):
        """
        Initialize the radial IA transform parameters.

        Parameters
        ----------
        scale : float
            Pixel scale in arcseconds. Retained as transform metadata.
        hlr : float
            Half-light radius used to normalize radial distances.
        A : float, optional
            Distortion amplitude at ``r = hlr``.
        phi : float, optional
            Alignment angle in radians.
        beta : float, optional
            Power-law index controlling the radial amplitude profile.
        center : sequence of float, optional
            Coordinate origin ``[x, y]`` for the radial profile. Defaults to
            ``[0, 0]``.
        clip_radius : float, optional
            Maximum normalized radius used in the amplitude law.
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
        Transform coordinates using the local IA shear at each radius.

        Parameters
        ----------
        coords : array-like
            Coordinate array with shape ``(2, npoints)``. NumPy and CuPy inputs
            are both supported.

        Returns
        -------
        array-like
            Coordinate array transformed by the radius-dependent reduced shear.
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
        Compute reduced-shear components for coordinates relative to center.

        Parameters
        ----------
        x, y : array-like
            Coordinates relative to the transform center. Integer arrays are
            promoted to floating point before the radial amplitude law is
            evaluated.

        Returns
        -------
        tuple of array-like
            ``(g1, g2)`` reduced-shear components using the same backend as the
            input arrays.

        Raises
        ------
        ValueError
            Raised if the requested distortion amplitude exceeds one.
        """
        xp = _coords_backend(x)
        dtype = np.result_type(_coords_dtype(x), _coords_dtype(y))
        if not np.issubdtype(dtype, np.floating):
            dtype = np.float64

        x = xp.asarray(x, dtype=dtype)
        y = xp.asarray(y, dtype=dtype)
        hlr = xp.asarray(self.hlr, dtype=dtype)
        amp = xp.asarray(self.A, dtype=dtype)
        c2phi = xp.asarray(self.c2phi, dtype=dtype)
        s2phi = xp.asarray(self.s2phi, dtype=dtype)

        # Galactocentric radius in units of the half-light radius.
        radial_dist = xp.sqrt(x * x + y * y)
        radial_ratio = radial_dist / hlr

        # Prevent extrapolation beyond the adopted radial range.
        radial_ratio = xp.clip(
            radial_ratio,
            0.0,
            self.clip_radius,
        )

        # HLR-normalised reduced-shear amplitude.
        g_abs = amp * radial_ratio**self.beta

        if bool(xp.any(xp.abs(g_abs) >= 1.0)):
            raise ValueError(
                "Requested reduced-shear amplitude must satisfy |g| < 1.",
                xp.abs(g_abs),
                0.0,
                1.0,
            )

        g1 = g_abs * c2phi
        g2 = g_abs * s2phi

        return g1, g2


IATransform = IaTransform


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

        Parameters
        ----------
        backend : module, optional
            Array backend to use for the copied transform. Use NumPy by default
            or pass CuPy for GPU-backed arrays.
        dtype : dtype, optional
            Floating dtype for the copied transform arrays. Defaults to this
            transform's current dtype.

        Returns
        -------
        LensTransform
            New transform with equivalent parameters stored on the requested
            backend and dtype.
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


AffineLensingTransform = LensTransform


def e_to_g(absesq):
    """
    Convert distortion amplitude squared to reduced shear scale factor.

    This is a helper function not currently used in the transforms provided.

    Parameters
    ----------
    absesq : float or array-like
        Squared distortion amplitude. Values near zero use a Taylor
        expansion for numerical stability.

    Returns
    -------
    float or array-like
        Multiplicative factor that maps distortion components to reduced
        shear components.
    """
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
        e2g += unstable * (0.5 + absesq * (0.125 + absesq * (0.0625 + absesq * 0.0390625)))
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
