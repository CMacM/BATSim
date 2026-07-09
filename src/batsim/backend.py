"""Array backend selection and memory-management helpers."""

import gc
import types
import warnings

import numpy as np

_ARRAY_BACKEND = None
_CPU_FALLBACK_WARNED = False


def _get_array_backend(
    warning_message="CuPy unavailable; falling back to NumPy FFTs on CPU.",
):
    """
    Return CuPy if available, otherwise NumPy.

    The selected module is cached after the first call.
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
                warning_message,
                RuntimeWarning,
                stacklevel=2,
            )
            _CPU_FALLBACK_WARNED = True
        _ARRAY_BACKEND = np

    return _ARRAY_BACKEND


def _resolve_array_backend(backend):
    """
    Resolve a user-supplied backend selector to an array module.
    """
    if backend is None:
        return _get_array_backend()

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
    Clear cached backend allocations, primarily useful after CuPy OOMs.
    """
    xp = _resolve_array_backend(backend)
    _release_backend_memory(xp)


def sync_if_gpu(xp):
    """Synchronise the default CuPy stream when using a GPU backend."""
    if xp is not np:
        xp.cuda.Stream.null.synchronize()
