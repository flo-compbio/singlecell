import os
import pkg_resources

__version__ = pkg_resources.require('singlecell')[0].version

_root = os.path.abspath(os.path.dirname(__file__))

from . import util
from . import qc