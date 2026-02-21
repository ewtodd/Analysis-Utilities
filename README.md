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
nix flake init -t github:ewtodd/Analysis-Utilities --refresh
nix develop
```
<!---->
This creates a development environment with access to the libraries.
Compile local code with the included Makefile and run macros with `root -l macro.cpp+`.
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
- `AddSubplotLabel(label, x, y)` - Returns a drawn `TText` in NDC coordinates for subplot labels (e.g. "(a)", "(b)")
<!---->
**Output**:
- `SaveFigure(canvas, name, PlotSaveOptions)` - Saves to `plots/` directory using the format set in `SetStylePreferences`.
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
## Roadmap
<!---->
- [x] Implement support for converting CoMPASS binary files to ROOT
- [ ] Implement support for converting WaveDump binary files to ROOT (742 family digitizers only) - implemented, but as of now untested.
- [ ] Implement support for converting CoMPASS CSV files to ROOT
