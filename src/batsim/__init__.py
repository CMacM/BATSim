from ._version import __version__

version = __version__

from . import WCS, pltutil, stamp

try:
    from . import _gsinterface
except ImportError:
    _gsinterface = None

from .sim import clear_backend_memory, simulate_galaxy
from .stamp import Stamp
from .transforms import (
    AffineLensingTransform,
    FlexionTransform,
    IATransform,
    IaTransform,
    LensTransform,
)

__all__ = [
    "__version__",
    "version",
    "simulate_galaxy",
    "clear_backend_memory",
    "Stamp",
    "LensTransform",
    "AffineLensingTransform",
    "IaTransform",
    "IATransform",
    "FlexionTransform",
    "WCS",
    "pltutil",
    "stamp",
]
