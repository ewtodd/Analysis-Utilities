# Waveform Digital Signal Processing Toolkit for CAEN digitizers + CoMPASS
## Features: Baseline subtraction, waveform inversion, charge integration, calibration, histogramming, plotting utilities, etc... Machine learning utilities via sklearn in Python are here but are not yet integrated into the flake. 

Example usage in nix flake, i.e. to be used via `nix develop`: 
```
{
  inputs = {
    waveform-analysis.url = "github:yourusername/Waveform-DSP-Toolkit";
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
