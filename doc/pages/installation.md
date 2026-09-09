Installation {#installation}
============

[TOC]

## Prerequisites

This project uses [Nix](https://nixos.org/) to manage dependencies. Install it
by following the instructions [here](https://nixos.org/download/), and ensure
[flakes are enabled](https://nixos.wiki/wiki/Flakes#Enable_flakes_permanently)
in your Nix configuration.

## Setup

In a new project directory:

```bash
# C++ template (default) — ROOT + compiled libraries
nix flake init -t github:ewtodd/Analysis-Utilities --refresh
nix develop
```

This creates a development environment with access to the libraries. Compile
local code with the included Makefile and run macros with `root -l macro.cpp+`.

```bash
# Python template — includes analysis-utilities Python package and ML libraries
nix flake init -t github:ewtodd/Analysis-Utilities#python --refresh
nix develop
```

## Binary cache

The CUDA-overlaid ROOT and the `packages.cuda` library are large to build from
source. **e-desktop** (the primary development host for this project) re-serves
its nix store as a binary cache so other machines can fetch the pre-built
artifacts directly instead of rebuilding.

- **URL:** `https://cache.ethanwtodd.com`
- **Use on other hosts:** add the URL as a substituter (and trust its signing
  key) in your `nix.conf` / `nixos-rebuild` config. Once configured, `nix build`
  and `nix develop` will transparently pull `packages.cuda` and
  `packages.rootCuda` from the cache instead of compiling.

@warning Do not override the `nixpkgs` input (e.g. via
`inputs.utils.inputs.nixpkgs.follows`). Overriding it changes the derivation
hash and forces a local rebuild that the cache cannot satisfy.

## Building the documentation

`doxygen` and `graphviz` are already in both dev shells:

```bash
make docs                              # or:
cmake --build build --target docs
```

Output lands in `docs/html/index.html`, warnings in
`docs/doxygen-warnings.log`. Both paths are thin wrappers around `Doxyfile.in`,
whose `@VARIABLE@` placeholders are substituted at build time — do not run
`doxygen` against it directly.
