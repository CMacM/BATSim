"""Array backend selection and memory-management helpers."""

import gc
import types
import numpy as np


def _get_array_backend():
    """
    Return BATSim's default array backend.

    NumPy is the default backend. CuPy is used only when explicitly requested
    through ``_resolve_array_backend("cp")`` or supplied as a module object.
    """
    return np


def _resolve_array_backend(backend):
    """
    Resolve a user-supplied backend selector to an array module.

    ``None`` means NumPy, string aliases request a specific backend, and
    module objects allow tests or callers to pass an already-imported array
    module.
    """
    if backend is None:
        return np

    if backend is np or backend in ("np", "numpy"):
        return np

    if backend in ("cp", "cupy"):
        try:
            import cupy as cp

            cp.cuda.runtime.getDeviceCount()
            return cp
        except Exception as exc:
            raise RuntimeError("backend='cp' was requested, but CuPy is unavailable.") from exc

    if isinstance(backend, types.ModuleType):
        return backend

    raise ValueError(f"Invalid backend '{backend}'; must be 'np', 'cp', a module, or None.")


def _to_numpy(xp, array):
    """Convert a NumPy/CuPy array to a NumPy array."""
    if xp is np:
        return np.asarray(array)
    return xp.asnumpy(array)


def _precision_dtypes(xp, precision):
    """Return real and complex dtypes for the requested FFT precision."""
    if precision == "single":
        return xp.float32, xp.complex64

    if precision == "double":
        return xp.float64, xp.complex128

    raise ValueError(f"Invalid precision '{precision}'; must be 'single' or 'double'.")


def _release_backend_memory(xp):
    """
    Release cached CuPy allocations after a simulation.

    NumPy has no comparable memory pools, so this is intentionally a no-op for
    CPU-backed renders.
    """
    if xp is np:
        return

    gc.collect()

    try:
        xp.cuda.Stream.null.synchronize()
    except Exception:
        pass

    try:
        xp.fft.config.get_plan_cache().clear()
    except Exception:
        pass

    gc.collect()

    try:
        xp.get_default_memory_pool().free_all_blocks()
    except Exception:
        pass

    try:
        xp.get_default_pinned_memory_pool().free_all_blocks()
    except Exception:
        pass

    gc.collect()


def clear_backend_memory(backend=None):
    """
    Clear cached memory held by the selected array backend.

    This is primarily useful after CuPy out-of-memory errors in notebooks or
    long-running sessions.  Calling it for the NumPy backend is a no-op.

    Parameters
    ----------
    backend : {"np", "numpy", "cp", "cupy"} or module or None, optional
        Backend to clear.  ``None`` selects NumPy, strings select NumPy or
        CuPy explicitly, and module objects are used directly.

    Returns
    -------
    None
    """
    xp = _resolve_array_backend(backend)
    _release_backend_memory(xp)


def sync_if_gpu(xp):
    """
    Synchronise the default CuPy stream when using a GPU backend.

    This is only used for timing and profiling and is a no-op for the NumPy backend.
    """
    if xp is not np:
        xp.cuda.Stream.null.synchronize()
