# C++ Utilities for Analysis of Nuclear Measurements
## Currently supporting CAEN digitizers (CoMPASS/wavedump)
## Features: Raw waveform processing and plotting utilities
Usage in a new project directory:

```
nix flake init -t github:ewtodd/Analysis-Utilities --refresh
```

This will create a flake.nix file containing a development environment that has access to the libraries. 

# Roadmap
- [ ] Implement "true" digital constant fraction discrimination for triggering. 
- [ ] Implement support for converting binary files (including those from wavedump) to ROOT so that the libraries can be used. 
- [ ] Implement support for converting CoMPASS csv files to ROOT so that the libraries can be used.
