# Digital Signal Processing Toolkit for CAEN digitizers + CoMPASS + wavedump
## Features: Baseline subtraction, waveform inversion, charge integration, calibration, histogramming, plotting utilities, etc...
Machine learning utilities via sklearn in Python are in this repository but are not yet integrated into the flake.
<!---->
Example usage in nix flake, i.e.
to be used via `nix develop`:
```
{
  inputs = {
    waveform-analysis.url = "github:ewtodd/Waveform-DSP-Toolkit";
  };
  outputs = { self, nixpkgs, waveform-analysis }:
    # ...
    devShells.default = pkgs.mkShell {
      buildInputs = [
        waveform-analysis.packages.${system}.default
      ];
    };
}
```
This will make all the header and library files available in the result development shell.
Examples of actually using the headers/libraries are in the example-macros folder.
