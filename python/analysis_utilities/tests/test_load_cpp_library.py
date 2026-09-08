"""Tests for analysis_utilities.load_cpp_library."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

ROOT = pytest.importorskip("ROOT")

try:
    from analysis_utilities import _CPP_HEADERS, load_cpp_library

    load_cpp_library()
except (RuntimeError, OSError) as exc:
    pytest.skip(f"C++ library not loadable: {exc}", allow_module_level=True)


# One representative symbol per header, so a failure names the header that
# stopped being declared rather than just a missing class. Two of these do not
# match their filename and are the ones most likely to be guessed wrong:
# IOUtils.hpp declares namespace IO, and BinaryUtils.hpp declares the reader
# classes directly with no BinaryUtils of its own.
#
# RooFitPhotopeakKernels.hpp is absent on purpose: it declares only CUDA launch
# wrappers, whose symbols exist solely in a CUDA build, so there is nothing here
# that a CPU build could resolve.
HEADER_SYMBOLS = (
    ("IOUtils.hpp", "IO"),
    ("PlottingUtils.hpp", "PlottingUtils"),
    ("InitUtils.hpp", "InitUtils"),
    ("BinaryUtils.hpp", "CoMPASSReader"),
    ("WaveformProcessingUtils.hpp", "WaveformProcessingUtils"),
    ("FittingUtils.hpp", "FittingUtils"),
    ("RooFitPhotopeakPdfs.hpp", "RooLowExpTail"),
    ("RooFitUtils.hpp", "RooFitUtils"),
    ("InteractiveEditorX11Guard.hpp", "AUXErrorHandlerSave"),
    ("InteractiveFitEditor.hpp", "InteractiveFitEditor"),
    ("InteractiveRooFitEditor.hpp", "InteractiveRooFitEditor"),
    ("InteractiveSimultaneousFitEditor.hpp", "InteractiveSimultaneousFitEditor"),
)


def _header_dir() -> Path | None:
    """Locate the directory holding the .hpp files, or None if not found.

    Works from a source checkout and from an installed package, which do not
    keep the headers in the same place.
    """
    source_tree = Path(__file__).resolve().parents[3] / "include"
    if (source_tree / "PlottingUtils.hpp").is_file():
        return source_tree
    for entry in os.environ.get("ROOT_INCLUDE_PATH", "").split(os.pathsep):
        if entry and (Path(entry) / "PlottingUtils.hpp").is_file():
            return Path(entry)
    return None


def test_every_installed_header_is_declared() -> None:
    """_CPP_HEADERS must match the headers on disk, exactly.

    CMake installs include/*.hpp by glob, so a new header ships without anyone
    touching the Python side. When that happens its classes are missing from
    PyROOT with no error at all, which is the failure this list prevents.
    """
    header_dir = _header_dir()
    if header_dir is None:
        pytest.skip("Could not locate the Analysis-Utilities headers")

    on_disk = {path.name for path in header_dir.glob("*.hpp")}
    assert on_disk, f"No headers found in {header_dir}"

    undeclared = sorted(on_disk - set(_CPP_HEADERS))
    stale = sorted(set(_CPP_HEADERS) - on_disk)
    assert not undeclared, (
        "These headers are installed but missing from _CPP_HEADERS, so their "
        f"classes are invisible to PyROOT: {undeclared}")
    assert not stale, (
        "These headers are in _CPP_HEADERS but not installed, so "
        f"load_cpp_library() will raise: {stale}")


@pytest.mark.parametrize("header,symbol", HEADER_SYMBOLS)
def test_header_symbol_resolves(header: str, symbol: str) -> None:
    assert hasattr(ROOT, symbol), (
        f"{symbol} is missing, so {header} was not declared")


def test_io_utils_is_reachable_as_a_namespace() -> None:
    """IOUtils.hpp declares namespace IO, not a class named IOUtils."""
    assert hasattr(ROOT.IO, "OpenForReading")
    assert hasattr(ROOT.IO, "OpenForWriting")
    assert hasattr(ROOT.IO, "ScopedRootLock")


def test_waveform_processing_types_resolve_together() -> None:
    """The processing entry points are useless without their companion types."""
    for symbol in ("FileProcessingConfig", "WaveformFeatures",
                   "ProcessingStats", "InputFormat"):
        assert hasattr(ROOT, symbol), f"{symbol} is missing"
