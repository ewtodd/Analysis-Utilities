# Analysis Toolkit for Nuclear Measurements 
## Currently supporting CAEN digitizers (CoMPASS/wavedump)
## Features: Baseline subtraction, waveform inversion, pulse height analysis, charge integration, calibration, histogramming, plotting utilities, etc...
Usage in a new project directory:

```
nix flake init -t github:ewtodd/Nuclear-Measurement-Toolkit
```

This will create a flake.nix file containing a development environment that has access to the libraries. It also has example macros showing use with a fast inorganic scintillator and common radioactive check sources.

# Roadmap
## Short term:
- [ ] Implement "true" digital constant fraction discrimination for triggering. 
- [ ] Implement support for converting binary files (including those from wavedump) to ROOT so that the libraries can be used. 
- [ ] Implement support for converting CoMPASS csv files to ROOT so that the libraries can be used.
- [ ] Switch to using control files instead of coding the acquisition parameters. 
## Long term:
- [ ] Test on a platform that is not x86_64-linux.
- [ ] Support other industry standard digitizers.
