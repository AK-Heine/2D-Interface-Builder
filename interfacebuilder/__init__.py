"""2D interface builder."""

import sys
from distutils.version import LooseVersion

if sys.version_info[0] == 2:
    raise ImportError("Requires Python3. This is Python2.")

__version__ = "1.0.0"
