"""Analysis utilities for nuclear measurement data."""

from analysis_utilities.init_utils import set_root_preferences
from analysis_utilities.simultaneous_fit import (FitResult,
                                                   SimultaneousFit)
from analysis_utilities.unbinned_fit import (UnbinnedFit,
                                               UnbinnedFitResult)

__version__ = "@VERSION@"

__all__ = [
    "load_cpp_library",
    "set_root_preferences",
    "SimultaneousFit",
    "FitResult",
    "UnbinnedFit",
    "UnbinnedFitResult",
]

_cpp_loaded = False


def load_cpp_library():
    """Load the C++ Analysis-Utilities library into ROOT.

    After calling this, ROOT.PlottingUtils, ROOT.InitUtils,
    ROOT.PlotSaveFormat, and ROOT.PlotSaveOptions are available.

    Returns:
        ROOT module (for convenience)
    """
    global _cpp_loaded
    import ROOT

    if not _cpp_loaded:
        if ROOT.gSystem.Load("lib-analysis-utils") < 0:
            raise RuntimeError(
                "Could not load lib-analysis-utils.so. "
                "Make sure LD_LIBRARY_PATH includes the library directory.")
        ROOT.gInterpreter.Declare('#include "PlottingUtils.hpp"')
        ROOT.gInterpreter.Declare('#include "FittingUtils.hpp"')
        ROOT.gInterpreter.Declare('#include "InitUtils.hpp"')
        _cpp_loaded = True

    return ROOT
