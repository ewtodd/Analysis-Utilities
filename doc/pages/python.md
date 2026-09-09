Python package {#python}
==============

The `analysis-utilities` Python package provides two things: a bridge that makes
the whole C++ API usable from Python scripts, and a loader for efficiently
reading ROOT TTrees into numpy arrays and pandas DataFrames for use with machine
learning libraries.

[TOC]

## Setup

To use the Python package in a downstream project, add the `pythonPackage`
output to your flake and include it in a Python environment:

```nix
# In your project's flake.nix
let
  pkgs = nixpkgs.legacyPackages.${system};
  analysis-utils = utils.packages.${system}.default;
  analysis-utilities-py = utils.packages.${system}.pythonPackage;
in
{
  devShells.default = pkgs.mkShell {
    buildInputs = [
      analysis-utils
      (pkgs.python3.withPackages (ps: [
        analysis-utilities-py
        ps.numpy
        ps.pandas
        # add ML libraries here, e.g. ps.scikit-learn
      ]))
      pkgs.root
    ];
    shellHook = ''
      export LD_LIBRARY_PATH="${analysis-utils}/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
      export ROOT_INCLUDE_PATH="${analysis-utils}/include''${ROOT_INCLUDE_PATH:+:$ROOT_INCLUDE_PATH}"
    '';
  };
}
```

## C++ bridge

`load_cpp_library()` loads the C++ shared library into ROOT and declares every
public header, making the whole C++ API available through PyROOT:

```python
from analysis_utilities import load_cpp_library

ROOT = load_cpp_library()

ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)
c = ROOT.PlottingUtils.GetConfiguredCanvas(False)

cfg = ROOT.FileProcessingConfig()
cfg.polarity = -1
processor = ROOT.WaveformProcessingUtils(cfg)
processor.ProcessFile("run.root", "run_features")
```

Two headers do not export a symbol matching their filename:

- `IOUtils.hpp` declares `namespace IO`, so it is `ROOT.IO.OpenForReading`, not
  `ROOT.IOUtils`.
- `BinaryUtils.hpp` declares the reader classes directly — `ROOT.CoMPASSReader`,
  `ROOT.WaveDump742Reader`, `ROOT.SOLReader` and their data/hit companions —
  with no `BinaryUtils` of its own.

@warning A class is visible to PyROOT only when Cling has parsed a declaration
for it, so a header added to `include/` must also be added to `_CPP_HEADERS` in
`python/analysis_utilities/__init__.py` or its classes are silently absent. The
test suite enforces that the two lists match.

## Project-rooted output paths

`set_root_preferences` configures the ROOT environment and pins both the plot
output base and the ROOT-files I/O base to absolute paths, so downstream scripts
produce identical output regardless of CWD. The same call wires
`PlottingUtils::SaveFigure` and `IO::OpenForReading` / `IO::OpenForWriting`
(Python: `analysis_utilities.io.open_for_reading` / `open_for_writing`) to the
configured locations on both the Python and C++ sides.

```python
from pathlib import Path
from analysis_utilities import (set_root_preferences, open_for_reading,
                                open_for_writing)

PROJECT_ROOT = Path(__file__).resolve().parent.parent
ROOT = set_root_preferences(plots_dir=PROJECT_ROOT / "plots" / "raw",
                            root_files_dir=PROJECT_ROOT / "root_files")

# Read: opens <root_files_base>/<subpath>
fin = open_for_reading("filtered/run42.root")

# Write: creates parent dirs under <root_files_base> before opening
fout = open_for_writing("raw/calibration_2026-05-01.root")   # default RECREATE
fupd = open_for_writing("filtered/run42.root", mode="UPDATE")

# Absolute subpaths bypass the base entirely.
```

Omitting `plots_dir` or `root_files_dir` prints a warning and falls back to the
CWD-relative defaults (`"plots"` and `"root_files"`). The IO helpers use only
`ROOT.gSystem` for filesystem operations (path joining via `ConcatFileName`,
parent creation via `mkdir`); no `pathlib` is imported on the open code path.

## TTree loader

`load_tree_data()` reads ROOT TTrees into pandas DataFrames and numpy arrays. It
handles type detection, TChain construction for multiple files, and optional
waveform array branches. Results are automatically cached to disk (as
pickle/npy files) so that subsequent loads skip the ROOT I/O entirely. The cache
is invalidated when any source ROOT file is newer than the cached file.

```python
from analysis_utilities.io import load_tree_data

# Load scalar branches into a DataFrame (cached to df_cache/ by default)
df = load_tree_data("output.root", tree_name="features")

# Load from multiple files with event limit
df = load_tree_data(
    ["run1.root", "run2.root"],
    tree_name="features",
    max_events=50000,
)

# Load waveforms alongside scalar data
df, waveforms = load_tree_data(
    "output.root",
    tree_name="features",
    array_branch="waveform",
)
# waveforms is a 2-D numpy array with shape (n_events, n_samples)

# Disable caching
df = load_tree_data("output.root", cache_dir=None)
```

**Parameters**:

- `root_files` — path or list of paths to ROOT files (combined via TChain)
- `tree_name` — TTree name (default: `"features"`)
- `scalar_branches` — branch names to load, or `None` to auto-detect all scalar
  branches. Caching requires `None` (all branches); pass `cache_dir=None` if you
  need a subset.
- `array_branch` — name of a `TArrayF` / `TArrayS` branch to load as a 2-D numpy
  array
- `max_events` — cap on number of events to read
- `cache_dir` — directory for cached files (default: `"df_cache"`). Set to
  `None` to disable caching.

**Returns** a `pandas.DataFrame` of scalar data. If `array_branch` is specified,
returns a tuple of `(DataFrame, numpy.ndarray)`.

Supported branch types: `Float_t`, `Double_t`, `Int_t`, `UInt_t`, `Short_t`,
`UShort_t`, `Long64_t`, `ULong64_t`, `UChar_t`.

## Fixed-size array loader

`load_leaf_array_data()` reads C-style fixed-size array leaves (branches created
as e.g. `tree->Branch("adc", buf, "adc[18]/s")`) directly into 2-D numpy arrays,
in the leaf's native dtype. It is the companion to `load_tree_data`, which
handles scalar branches plus `TArrayF` / `TArrayS` object branches; object
branches are rejected here. Results are cached as one `.npz` per (file set,
tree) inside `cache_dir`, invalidated when any source ROOT file is newer.

```python
from analysis_utilities.io import load_leaf_array_data

arrays = load_leaf_array_data(
    "output.root",
    tree_name="Data_R",
    array_branches=["Trace0", "Dig0"],
    max_events=50000,
)
# arrays is a dict: {"Trace0": (n_events, trace_len), "Dig0": (n_events, trace_len)}
```

**Parameters**:

- `root_files` — path or list of paths to ROOT files (combined via TChain)
- `tree_name` — TTree name
- `array_branches` — names of fixed-size array leaf branches to load
- `max_events` — cap on number of events to read
- `cache_dir` — directory for the cached `.npz` (default: `"df_cache"`). Set to
  `None` to disable caching.
