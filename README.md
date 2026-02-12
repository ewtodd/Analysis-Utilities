# Analysis-Utilities
<!---->
C++ utilities for analysis of nuclear measurement data, built on [ROOT](https://root.cern/).
Supports CAEN digitizers via CoMPASS and WaveDump acquisition software.
Warning: Currently in active development; breaking changes are all but guaranteed...
<!---->
## Installation
<!---->
### Prerequisites
<!---->
This project uses [Nix](https://nixos.org/) to manage dependencies. Install it with:
<!---->
```bash
curl -L https://nixos.org/nix/install | sh
```
<!---->
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
- **Parallel file processing** - Process multiple files concurrently using `std::async` with configurable worker count (defaults to 4). Files are dispatched in batches, with each worker getting its own `WaveformProcessingUtils` instance. Requires ROOT thread safety (`ROOT::EnableThreadSafety()`), which is handled automatically.
<!---->
Outputs processed data to ROOT TTrees with optional waveform storage.
<!---->
### FittingUtils
<!---->
Fit spectral peaks with configurable models.
<!---->
**Standard model** (5 parameters):
- Gaussian peak + linear background
<!---->
**Detailed model** (10 parameters):
Includes standard:
- Gaussian peak + linear background
<!---->
Details based on
L.C.
Longoria, A.H.
Naboulsi, P.W.
Gray, T.D.
MacMahon,
Analytical peak fitting for gamma-ray spectrum analysis with Ge detectors,
Nuclear Instruments and Methods in Physics Research Section A: Accelerators, Spectrometers, Detectors and Associated Equipment,
Volume 299, Issues 1–3,
1990,
Pages 308-312,
ISSN 0168-9002,
https://doi.org/10.1016/0168-9002(90)90797-A.
Quotes refer to the aforementioned article.
Reading the article reveals that the models discussed are not specific to Ge detectors!
<!---->
- Step function ("When part of the photon energy escapes from the
active volume of the detector, an event is recorded in the left-hand side of the peak, creating a tail with a steplike shape. It is assumed that if the detector had no tailing and an extremely narrow resolution, then a step-like distribution would result, with a cutoff at the centroid." In practice of course this means smearing a step by the resolution of the detector as determined by other parameters from the fit.)
- Low-energy tail (incomplete charge collection, instability, etc.)
- High-energy tail (pileup)
<!---->
**Multi-peak fitting**:
- Double and triple peak variants of both models
- Constrained fitting using results from previous fits
<!---->
All fits produce structured results with parameter values, errors, and reduced chi-squared.
<!---->
### PlottingUtils
<!---->
Configure ROOT graphics with consistent styling.
<!---->
- Style presets for publication-quality figures
- Graph and histogram configuration helpers
- 2D histogram support with color palettes
- Legend and subplot label utilities
- Save figures in multiple formats (with optional log scale variants)
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
