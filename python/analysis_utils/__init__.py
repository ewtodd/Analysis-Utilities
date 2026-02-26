"""Analysis utilities for nuclear measurement data."""

__version__ = "@VERSION@"

#from analysis_utils.fitting import FitResult, PeakFitResult, UnbinnedFit

_cpp_loaded = False


def load_cpp_library():
    """Load the C++ Analysis-Utilities library into ROOT.

    After calling this, ROOT.PlottingUtils, ROOT.PlotSaveFormat,
    and ROOT.PlotSaveOptions are available.

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
        _cpp_loaded = True

    return ROOT
