"""Coordinate transforms used during image rendering."""

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


class Transform:
    """
    Base class for BATSim coordinate transforms.

    Transform convention
    --------------------
    BATSim transforms implement backward coordinate mappings. Given coordinates
    ``x`` on the regular output-image grid, ``transform(x)`` returns coordinates
    ``u`` at which the source surface-brightness profile should be evaluated:

        I_output(x) = I_source(u),    u = transform(x).

    For a locally affine transformation with forward geometric matrix ``M``, a source
    feature at ``u`` appears at

        x = M @ u.

    The corresponding BATSim coordinate transform must therefore return

        u = inv(M) @ x.

    Although the inverse matrix is used to sample the source profile, visible
    features undergo the forward transformation. For example, a source feature
    at ``u0`` appears where

        inv(M) @ x = u0,

    which implies

        x = M @ u0.

    This backward, or pull-based, convention evaluates the source profile once
    for every output location and avoids the gaps and overlaps that can occur
    when source samples are pushed forwards onto an output grid.

    All transform implementations in this module must follow this convention:
    their inputs are output-plane coordinates and their return values are
    source-plane sampling coordinates.

    Subclasses implement only the transform-specific coordinate logic in
    ``_transform_relative``. The base class handles backend selection, dtype
    promotion, coordinate coercion, reference-centre subtraction/addition, and
    moving transform arrays between NumPy and CuPy via ``to_backend``.

    Parameters
    ----------
    center : sequence of float, optional
        Reference coordinate ``[x, y]``. Defaults to ``[0, 0]``.
    backend : module, optional
        NumPy or CuPy. If None, NumPy is used for stored transform arrays.
        Calls to ``transform`` still preserve the backend of the input
        coordinate array.
    dtype : dtype, optional
        Floating dtype for stored transform arrays.
    """

    _backend_array_attributes = ("ref_vec",)

    def __init__(self, center=None, backend=None, dtype=None):
        self.xp = np if backend is None else backend

        if dtype is None:
            dtype = self.xp.float64

        if center is None:
            center = [0.0, 0.0]

        self.dtype = dtype
        self.center = (float(center[0]), float(center[1]))
        self.ref_vec = self.xp.asarray(
            [[self.center[0]], [self.center[1]]],
            dtype=self.dtype,
        )

    def transform(self, coords):
        """
        Transform output-plane coordinates to source-plane sampling coordinates.

        Parameters
        ----------
        coords : array-like
            Coordinate array with shape ``(2, npoints)``.

        Returns
        -------
        array-like
            Transformed coordinate array using the same backend as ``coords``.
        """
        xp, dtype, coords, ref_vec = self._prepare_coords(coords)
        transformed = self._transform_relative(
            coords - ref_vec,
            xp=xp,
            dtype=dtype,
        )
        return transformed + ref_vec

    def _prepare_coords(self, coords):
        """
        Coerce coordinate arrays and return backend-aware transform metadata.

        Returns
        -------
        tuple
            ``(xp, dtype, coords, ref_vec)`` where ``xp`` is the coordinate
            backend and ``ref_vec`` is stored on that backend.
        """
        xp = _coords_backend(coords)
        dtype = _coords_dtype(coords)

        if dtype is None or not np.issubdtype(dtype, np.floating):
            dtype = self.dtype

        coords = xp.asarray(coords, dtype=dtype)
        ref_vec = self._as_backend_array(self.ref_vec, xp=xp, dtype=dtype)
        return xp, dtype, coords, ref_vec

    def _transform_relative(self, coords_relative, xp, dtype):
        """
        Transform coordinates relative to ``self.center``.

        Subclasses should override this method and return an array with shape
        ``(2, npoints)`` using the supplied backend ``xp``.
        """
        raise NotImplementedError("Transform subclasses must implement _transform_relative.")

    def to_backend(self, backend=None, dtype=None):
        """
        Return a copy of this transform on another backend/dtype.

        Parameters
        ----------
        backend : module, optional
            Array backend to use for stored transform arrays. Use NumPy by
            default or pass CuPy for GPU-backed arrays.
        dtype : dtype, optional
            Floating dtype for stored transform arrays. Defaults to this
            transform's current dtype.

        Returns
        -------
        Transform
            New transform with equivalent parameters stored on the requested
            backend and dtype.
        """
        xp = np if backend is None else backend

        if dtype is None:
            dtype = self.dtype

        new = object.__new__(type(self))
        new.__dict__ = self.__dict__.copy()
        new.xp = xp
        new.dtype = dtype

        for attr in self._backend_array_attributes:
            if hasattr(self, attr):
                setattr(
                    new,
                    attr,
                    xp.asarray(self._to_numpy(getattr(self, attr)), dtype=dtype),
                )

        return new

    @staticmethod
    def _to_numpy(array):
        """Convert a NumPy or CuPy array to NumPy."""
        if isinstance(array, np.ndarray):
            return array

        return array.get()

    def _as_backend_array(self, array, xp, dtype=None):
        """Return ``array`` on ``xp`` without implicit CuPy-to-NumPy conversion."""
        if xp is np:
            return np.asarray(self._to_numpy(array), dtype=dtype)

        if _coords_backend(array) is xp:
            return xp.asarray(array, dtype=dtype)

        return xp.asarray(self._to_numpy(array), dtype=dtype)


class FlexionTransform(Transform):
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

    _backend_array_attributes = ("ref_vec", "s2l_mat", "s2l_mat_inv", "D")

    def __init__(
        self,
        gamma1,
        gamma2,
        kappa,
        F1=0,
        F2=0,
        G1=0,
        G2=0,
        center=None,
        backend=None,
        dtype=None,
    ):
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
        super().__init__(center=center, backend=backend, dtype=dtype)

        self.s2l_mat = self.xp.asarray(
            [[1 - kappa - gamma1, -gamma2], [-gamma2, 1 - kappa + gamma1]],
            dtype=self.dtype,
        )
        self.s2l_mat_inv = self.xp.linalg.inv(self.s2l_mat)
        D1 = (
            -1
            / 2
            * self.xp.asarray(
                [[3 * F1 + G1, F2 + G2], [F2 + G2, F1 - G1]],
                dtype=self.dtype,
            )
        )
        D2 = (
            -1
            / 2
            * self.xp.asarray(
                [[F2 + G2, F1 - G1], [F1 - G1, 3 * F2 - G2]],
                dtype=self.dtype,
            )
        )
        self.D = self.xp.stack([D1, D2], axis=2)
        return

    def _transform_relative(self, coords_relative, xp, dtype):
        s2l_mat = self._as_backend_array(self.s2l_mat, xp=xp, dtype=dtype)
        d_tensor = self._as_backend_array(self.D, xp=xp, dtype=dtype)

        return s2l_mat @ coords_relative + 0.5 * xp.einsum(
            "ijk,jl,kl->il",
            d_tensor,
            coords_relative,
            coords_relative,
        )

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
        return super().transform(coords)

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
        xp, dtype, coords, ref_vec = self._prepare_coords(coords)
        coords_relative = coords - ref_vec
        s2l_mat_inv = self._as_backend_array(self.s2l_mat_inv, xp=xp, dtype=dtype)
        d_tensor = self._as_backend_array(self.D, xp=xp, dtype=dtype)

        theta_0 = xp.einsum("ij,jk", s2l_mat_inv, coords_relative)
        theta_1 = (
            -1
            / 2
            * xp.einsum(
                "in,ijk,jl,lo,km,mo->no",
                s2l_mat_inv,
                d_tensor,
                s2l_mat_inv,
                coords_relative,
                s2l_mat_inv,
                coords_relative,
            )
        )
        return theta_0 + theta_1 + ref_vec


class IaTransform(Transform):
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
        backend=None,
        dtype=None,
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
        backend : module, optional
            NumPy or CuPy. If None, NumPy is used for stored transform arrays.
        dtype : dtype, optional
            Floating dtype for stored transform arrays.
        """
        super().__init__(center=center, backend=backend, dtype=dtype)

        # initialise important class variables
        self.A = A
        self.phi = phi
        self.c2phi = np.cos(2 * self.phi)
        self.s2phi = np.sin(2 * self.phi)
        self.beta = beta
        self.scale = scale
        self.hlr = hlr
        self.xcen = self.center[0]
        self.ycen = self.center[1]
        self.clip_radius = clip_radius

        return

    def _transform_relative(self, coords_relative, xp, dtype):
        x, y = coords_relative

        # Get g1 and g2 shear components as a function of image position
        g1, g2 = self.get_g1g2(x, y)

        # Normalisation ensures magnification is zero
        inv_norm = 1.0 / xp.sqrt(1.0 - g1 * g1 - g2 * g2)

        # Coordinate transform uses the INVERSE of the shear matrix
        x_prime = inv_norm * ((1.0 - g1) * x - g2 * y)
        y_prime = inv_norm * (-g2 * x + (1.0 + g1) * y)
        return xp.stack([x_prime, y_prime], axis=0)

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
        return super().transform(coords)

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
        dtype = np.result_type(_coords_dtype(x) or self.dtype, _coords_dtype(y) or self.dtype)
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


class LensTransform(Transform):
    """
    Affine lensing transform supporting either NumPy or CuPy arrays.

    If CuPy is available and the input coordinates are CuPy arrays, the
    transform runs natively on GPU.
    """

    _backend_array_attributes = ("ref_vec", "lens_mat")

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
        super().__init__(center=center, backend=backend, dtype=dtype)

        self.lens_mat = self.xp.asarray(
            [
                [1.0 - kappa - gamma1, -gamma2],
                [-gamma2, 1.0 - kappa + gamma1],
            ],
            dtype=self.dtype,
        )

    def _transform_relative(self, coords_relative, xp, dtype):
        lens_mat = self._as_backend_array(self.lens_mat, xp=xp, dtype=dtype)
        return lens_mat @ coords_relative

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
        return super().transform(coords)


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
