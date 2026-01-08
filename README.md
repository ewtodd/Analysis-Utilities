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
- [x] Implement support for converting CoMPASS binary files to ROOT.
- [ ] Implement support for converting wavedump binary files to ROOT.
- [ ] Implement support for converting CoMPASS csv files to ROOT.
