from analysis_utils import load_cpp_library

def set_root_preferences(save_format=None):
    """Initialize ROOT environment: batch mode, style, output directories.

    Parameters
    ----------
    save_format : ROOT.PlotSaveFormat, optional
        kPNG (default) or kPDF. Controls line widths in PlottingUtils.
    """
    ROOT = load_cpp_library()
    if save_format is None:
        save_format = ROOT.PlotSaveFormat.kPNG
    ROOT.InitUtils.SetROOTPreferences(save_format)
