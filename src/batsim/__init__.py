# Basically copied from the galsim __init__.py
import os
from pathlib import Path

from ._version import __version__

version = __version__

batsim_dir = Path(__file__).parent
galsim_dir = Path(batsim_dir).parent.parent / 'extern' / 'GalSim'
include_dir = Path(galsim_dir) / 'include'
from . import _gsinterface
lib_file = os.path.abspath(_gsinterface.__file__)

from .transforms import *
from .stamp import *
from .pltutil import *
#from . import WCS
