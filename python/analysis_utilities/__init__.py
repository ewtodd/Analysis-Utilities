"""Analysis utilities for nuclear measurement data."""

from analysis_utilities.init_utils import set_root_preferences
from analysis_utilities.io import (get_root_files_base_dir, open_for_reading,
                                   open_for_writing, set_root_files_base_dir)

__version__ = "@VERSION@"

__all__ = [
    "load_cpp_library",
    "set_root_preferences",
    "set_root_files_base_dir",
    "get_root_files_base_dir",
    "open_for_reading",
    "open_for_writing",
]

# Every public header, in dependency order. Loading the shared library is not
# enough to make a class visible to PyROOT: Cling resolves a name only if it
# has parsed a declaration for it, or if the class carries a dictionary entry
# in libanalysis-utils.rootmap. Only the RooFit photopeak PDFs have the latter,
# so everything else has to be declared here or it is simply absent.
_CPP_HEADERS = (
    "IOUtils.hpp",
    "PlottingUtils.hpp",
    "InitUtils.hpp",
    "BinaryUtils.hpp",
    "WaveformProcessingUtils.hpp",
    "FittingUtils.hpp",
    "RooFitPhotopeakKernels.hpp",
    "RooFitPhotopeakPdfs.hpp",
    "RooFitUtils.hpp",
    "InteractiveEditorX11Guard.hpp",
    "InteractiveFitEditor.hpp",
    "InteractiveRooFitEditor.hpp",
    "InteractiveSimultaneousFitEditor.hpp",
)

_cpp_loaded = False


def load_cpp_library():
    """Load the C++ Analysis-Utilities library into ROOT.

    Declares every public header, so the whole C++ API is reachable from
    Python afterwards:

    - ``ROOT.PlottingUtils``, ``ROOT.InitUtils``, ``ROOT.PlotSaveFormat``,
      ``ROOT.PlotSaveOptions``
    - ``ROOT.WaveformProcessingUtils``, ``ROOT.FileProcessingConfig``,
      ``ROOT.WaveformFeatures``, ``ROOT.ProcessingStats``, ``ROOT.InputFormat``
    - ``ROOT.FittingUtils``, ``ROOT.RooFitUtils``, and the photopeak PDFs
      ``ROOT.RooLowExpTail``, ``ROOT.RooHighExpTail``, ``ROOT.RooLowLinTail``,
      ``ROOT.RooStepShelf``
    - ``ROOT.InteractiveFitEditor``, ``ROOT.InteractiveRooFitEditor``,
      ``ROOT.InteractiveSimultaneousFitEditor`` (these need a display to be
      useful, but declaring them costs nothing in batch mode)

    Two headers do not export a symbol matching their filename, which is worth
    knowing before guessing at a name:

    - ``IOUtils.hpp`` declares ``namespace IO``, so it is
      ``ROOT.IO.OpenForReading``, not ``ROOT.IOUtils``.
    - ``BinaryUtils.hpp`` declares the reader types directly:
      ``ROOT.CoMPASSReader``, ``ROOT.WaveDump742Reader``, ``ROOT.SOLReader``,
      and their data/hit companions. There is no ``ROOT.BinaryUtils``.

    Returns:
        ROOT module (for convenience)
    """
    global _cpp_loaded
    import ROOT

    if not _cpp_loaded:
        if ROOT.gSystem.Load("libanalysis-utils") < 0:
            raise RuntimeError(
                "Could not load libanalysis-utils.so. "
                "Make sure LD_LIBRARY_PATH includes the library directory.")
        failed = []
        for header in _CPP_HEADERS:
            if not ROOT.gInterpreter.Declare('#include "%s"' % header):
                failed.append(header)
        if failed:
            raise RuntimeError(
                "Could not declare these Analysis-Utilities headers: "
                + ", ".join(failed)
                + ". Make sure ROOT_INCLUDE_PATH includes the directory "
                "holding the installed headers.")
        _cpp_loaded = True

    return ROOT
