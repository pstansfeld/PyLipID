##############################################################################
# PyLipID: A python module for analysing protein-lipid interactions
#
# Author: Wanling Song
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
##############################################################################

# --- BEGIN: macOS-safe bootstrap -------------------------------------------
# Ensures PyLipID runs safely on macOS by:
#   1) forcing multiprocessing to use 'spawn'
#   2) preventing GUI backends from loading in worker processes
import multiprocessing as _mp
import sys

try:
    # only force spawn on macOS
    if sys.platform == "darwin":
        _mp.set_start_method("spawn", force=True)
except RuntimeError:
    # Already set or protected — ignore
    pass

# Use a non-interactive backend so importing PyLipID never loads macOS AppKit
try:
    import matplotlib
    matplotlib.use("Agg", force=True)
except Exception:
    pass
# --- END: macOS-safe bootstrap ---------------------------------------------

from ._version import __version__
from .api import LipidInteraction
