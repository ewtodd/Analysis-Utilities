# Analysis-Utilities
<!---->
C++ utilities for analysis of nuclear measurement data, built on [ROOT](https://root.cern/), with a focus on ergonomics and performance.
Supports CAEN digitizers via CoMPASS and WaveDump acquisition software.
The flake also exposes a Python package, with a wrapper around the PlottingUtils as well as a method for efficiently loading TTrees into Python native data types (numpy array, pandas df) so that machine learning libraries in Python can be used.
Warning: Currently in active development; breaking changes are all but guaranteed...
<!---->
## Installation
<!---->
### Prerequisites
<!---->
This project uses [Nix](https://nixos.org/) to manage dependencies.
Install it by following the instructions [here](https://nixos.org/download/).
Ensure [flakes are enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently) in your Nix configuration.
<!---->
### Setup
<!---->
In a new project directory:
<!---->
```bash
# C++ template (default) — ROOT + compiled libraries
nix flake init -t github:ewtodd/Analysis-Utilities --refresh
nix develop
```
<!---->
This creates a development environment with access to the libraries.
Compile local code with the included Makefile and run macros with `root -l macro.cpp+`.
<!---->
```bash
# Python template — includes analysis-utilities Python package and ML libraries
nix flake init -t github:ewtodd/Analysis-Utilities#python --refresh
nix develop
```
<!---->
## Modules
<!---->
### BinaryUtils
<!---->
Read binary data files from CAEN acquisition software.
<!---->
**CoMPASSReader** - Reads CoMPASS binary files (`.bin`):
- Parses event headers, timestamps, energy values, and waveforms
- Supports all CoMPASS waveform codes (INPUT, RC_CR, TRAPEZOID, CFD, etc.)
- Decodes status flags (pileup, saturation, deadtime, trigger lost, etc.)
<!---->
**WaveDump742Reader** - Reads WaveDump binary files for DT5742 family digitizers:
- Parses event structure including DC offset and start index cell
- Optional timing corrections support
<!---->
### WaveformProcessingUtils
<!---->
Process raw waveforms and extract physical parameters.
<!---->
- **Baseline subtraction** - Configurable number of pre-trigger samples
- **Trigger finding** - Fraction-of-peak threshold with configurable polarity
- **Waveform cropping** - Extract region of interest around trigger
- **Feature extraction** - Pulse height, peak position, short/long integrals, PSD ratio
- **Quality cuts** - Reject clipped signals, baseline issues, negative integrals
- **Parallel file processing** - Process multiple files concurrently using `std::async` with configurable worker count (defaults to 4).
Files are dispatched in batches, with each worker getting its own `WaveformProcessingUtils` instance.
Requires ROOT thread safety (`ROOT::EnableThreadSafety()`), which is handled automatically.
<!---->
Outputs processed data to ROOT TTrees with optional waveform storage.
<!---->
### FittingUtils
<!---->
Fit gamma-ray spectral photopeaks with a composable model.
<!---->
**Base model**:
- Gaussian peak + linear or flat background
<!---->
**Optional components** (individually togglable, for semiconductor detectors with segmented electrodes):
- **Step function** - Models events where part of the photon energy escapes the active volume, producing a step-like distribution on the left side of the peak, smeared by the detector resolution.
- **Low-energy exponential tail** - Exponential tail below the photopeak from incomplete charge collection, charge trapping, etc.
- **Low-energy linear tail** - Linear tail component that can capture asymmetric tailing not well described by a single exponential.
- **High-energy exponential tail** - Exponential tail above the photopeak from pileup effects.
<!---->
Components are tested using a hybrid group-and-prune approach: low-side components (step + both low-energy tails) are enabled as a group, then individually pruned if they do not improve the fit.
The high-energy tail is tested independently.
A component is kept only if it improves the reduced chi-squared.
<!---->
**Multi-peak fitting**:
- Double and triple peak variants with all components enabled by default
- Constrained fitting using results from previous fits
<!---->
All fits produce structured results (`FitResult` containing `PeakFitResult` entries) with parameter values, errors, and reduced chi-squared.
Failed fits return -1 for all parameters.
Reduced chi-squared is displayed on fit plots by default.
Individual peak components are plotted summed with the background for readability, and multi-peak fits use distinct line styles per peak.
<!---->
![Fit example](assets/FitExample.png)
<!---->
**Interactive fitting**:
- `SetInteractive()` enables a GUI editor that opens after the automated fit completes.
The editor shows the histogram, total fit, and individual components with a residual panel, all updating in real-time as parameters are adjusted via sliders and number entries.
A fit range slider allows adjusting the fit range visually.
Parameters can be fixed/freed via checkboxes, and a Refit button runs Minuit with the current values as initial guesses.
Live chi-squared is displayed on the plot.
<!---->
![Interactive fit editor GUI](assets/GUI.png)
<!---->
- On accept, parameters and fit range are saved to a `.fits` file in the `plots/fits/` directory.
On subsequent runs with `SetInteractive()`, saved parameters are loaded automatically and the editor is skipped, so interactive results only need to be produced once.
Delete the `.fits` file to redo the interactive fit.
- It is recommended to re-run the fitting macro after an interactive session, as the saved parameters serve as much better starting points for Minuit and often produce a lower chi-squared on the second pass.
- Batch mode is temporarily disabled for the editor window and restored afterward, so plot popups are not affected.
<!---->
**References**:
- Boggs SE, Pike SN.
Analytical fitting of gamma-ray photopeaks in germanium cross strip detectors.
*Experimental Astronomy*.
2023;56(2-3):403-420.
doi: [10.1007/s10686-023-09914-8](https://doi.org/10.1007/s10686-023-09914-8).
- Longoria LC, Naboulsi AH, Gray PW, MacMahon TD.
Analytical peak fitting for gamma-ray spectrum analysis with Ge detectors.
*Nuclear Instruments and Methods in Physics Research A*.
1990;299(1-3):308-312.
doi: [10.1016/0168-9002(90)90797-A](https://doi.org/10.1016/0168-9002(90)90797-A).
<!---->
### PlottingUtils
<!---->
![Light output spectrum example with subplot text](assets/LightOutputSpectrumExample.png)
All-static utility class for publication quality ROOT graphics with consistent styling.
No object instantiation required.
<!---->
**Initialization**:
- `SetStylePreferences(PlotSaveFormat)` - Must be called before using other methods (warns if not).
Sets global ROOT style and chooses output format (`PlotSaveFormat::kPNG` or `PlotSaveFormat::kPDF`, defaults to PNG).
Calling `InitUtils::SetROOTPreferences()` takes care of this and is recommended.
<!---->
**Canvas**:
- `GetConfiguredCanvas(Bool_t logy)` - Returns a pre-configured 1200x800 `TCanvas` with grid and tick marks
<!---->
**Object configuration**:
- `ConfigureGraph` / `ConfigureAndDrawGraph` - Set line color/width, axis label/title sizes, and offsets on a `TGraph`
- `ConfigureHistogram` / `ConfigureAndDrawHistogram` - Same for `TH1`, plus fill style and axis division settings
- `Configure2DHistogram` / `ConfigureAndDraw2DHistogram` - Same for `TH2`, enables log-z and adjusts right margin for color axis
<!---->
**Annotations**:
- `AddLegend(x1, x2, y1, y2)` - Returns a drawn `TLegend` with consistent font/border styling.
Note: argument order is the extremely sane `(x1, x2, y1, y2)`, not the ROOT default `(x1, y1, x2, y2)`.
- `AddText(label, x, y, angle)` - Returns a drawn `TLatex` in NDC coordinates for arbitrary annotations (e.g. subplot labels like "(a)", "(b)", or any other text).
Optional `angle` (default 0) sets text rotation in degrees.
<!---->
**Output**:
- `SaveFigure(canvas, name, subdirectory, PlotSaveOptions)` - Saves to `plots/` directory (or `plots/<subdirectory>/` if specified) using the format set in `SetStylePreferences`.
`PlotSaveOptions` controls linear (`kLINEAR`), log (`kLOG`), or both (`kBOTH`, default).
Log variants are prefixed with `log_`.
<!---->
**Utilities**:
- `GetDefaultColors()` - Returns a 24-color palette of distinct ROOT colors
- `GetRandomName()` - Generates a random canvas name to avoid ROOT name collisions
<!---->
### InitUtils
<!---->
Initialization and file conversion utilities.
<!---->
- `SetROOTPreferences()` - Configure ROOT environment and plotting defaults
- `ConvertCoMPASSBinToROOT()` - Convert CoMPASS binary files to ROOT format
<!---->
## Python Package
<!---->
The `analysis-utilities` Python package provides two things: a bridge to use PlottingUtils from Python scripts, and a loader for efficiently reading ROOT TTrees into numpy arrays and pandas DataFrames for use with machine learning libraries.
<!---->
### Setup
<!---->
To use the Python package in a downstream project, add the `pythonPackage` output to your flake and include it in a Python environment:
<!---->
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
<!---->
### PlottingUtils bridge
<!---->
`load_cpp_library()` loads the C++ shared library into ROOT and declares the PlottingUtils header, making the full PlottingUtils API available through PyROOT:
<!---->
```python
from analysis_utilities import load_cpp_library
#
ROOT = load_cpp_library()
#
ROOT.PlottingUtils.SetStylePreferences(ROOT.PlotSaveFormat.kPNG)
c = ROOT.PlottingUtils.GetConfiguredCanvas(False)
# Use ROOT.PlottingUtils.ConfigureGraph, ConfigureHistogram, etc.
```
<!---->
### TTree loader
<!---->
`load_tree_data()` reads ROOT TTrees into pandas DataFrames and numpy arrays.
It handles type detection, TChain construction for multiple files, and optional waveform array branches.
Results are automatically cached to disk (as pickle/npy files) so that subsequent loads skip the ROOT I/O entirely.
The cache is invalidated when any source ROOT file is newer than the cached file.
<!---->
```python
from analysis_utilities.io import load_tree_data
#
# Load scalar branches into a DataFrame (cached to df_cache/ by default)
df = load_tree_data("output.root", tree_name="features")
#
# Load from multiple files with event limit
df = load_tree_data(
    ["run1.root", "run2.root"],
    tree_name="features",
    max_events=50000,
)
#
# Load waveforms alongside scalar data
df, waveforms = load_tree_data(
    "output.root",
    tree_name="features",
    array_branch="waveform",
)
# waveforms is a 2-D numpy array with shape (n_events, n_samples)
#
# Disable caching
df = load_tree_data("output.root", cache_dir=None)
```
<!---->
**Parameters**:
- `root_files` - Path or list of paths to ROOT files (combined via TChain)
- `tree_name` - TTree name (default: `"features"`)
- `scalar_branches` - Branch names to load, or `None` to auto-detect all scalar branches. Caching requires `None` (all branches); pass `cache_dir=None` if you need a subset.
- `array_branch` - Name of a `TArrayF`/`TArrayS` branch to load as a 2-D numpy array
- `max_events` - Cap on number of events to read
- `cache_dir` - Directory for cached files (default: `"df_cache"`). Set to `None` to disable caching.
<!---->
**Returns** a `pandas.DataFrame` of scalar data. If `array_branch` is specified, returns a tuple of `(DataFrame, numpy.ndarray)`.
<!---->
Supported branch types: `Float_t`, `Double_t`, `Int_t`, `UInt_t`, `Short_t`, `UShort_t`, `Long64_t`, `ULong64_t`, `UChar_t`.
<!---->
## Roadmap
<!---->
- [x] Implement support for converting CoMPASS binary files to ROOT
- [x] Implement support for converting WaveDump binary files to ROOT (742 family digitizers only) - implemented.
- [ ] Implement support for converting CoMPASS CSV files to ROOT
