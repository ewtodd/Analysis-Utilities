"""Python port of :cpp:class:`InitUtils` (the parts applicable to Python).

Binary file converters stay on the C++ side; only the ROOT environment
setup (styling, batch mode, standard output directories) is mirrored here.
"""

from __future__ import annotations

import os
from typing import Any, Optional


def set_root_preferences(save_format: Optional[Any] = None) -> Any:
    """Configure ROOT environment and plotting defaults.

    Mirrors :cpp:func:`InitUtils::SetROOTPreferences`: loads the C++
    library, calls :cpp:func:`PlottingUtils::SetStylePreferences`, forces
    the ROOT style, enables batch mode (no on-screen popups; also ensures
    consistent canvas sizing), and creates ``plots/`` and ``root_files/``
    in the current working directory if they do not exist.

    Parameters
    ----------
    save_format : ROOT.PlotSaveFormat, optional
        Output format for plots. Defaults to
        :cpp:enumerator:`PlotSaveFormat::kPNG` if not provided.

    Returns
    -------
    ROOT
        The ROOT module, for convenience.
    """
    from analysis_utilities import load_cpp_library
    ROOT = load_cpp_library()
    if save_format is None:
        save_format = ROOT.PlotSaveFormat.kPNG
    ROOT.PlottingUtils.SetStylePreferences(save_format)
    ROOT.gROOT.ForceStyle(True)
    ROOT.gROOT.SetBatch(True)
    for sub in ("plots", "root_files"):
        if not os.path.isdir(sub):
            os.makedirs(sub, exist_ok=True)
    return ROOT
