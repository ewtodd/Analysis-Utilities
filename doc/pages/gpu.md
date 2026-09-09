GPU acceleration {#cuda}
================

Opt-in CUDA evaluation for the RooFit photopeak backend.

[TOC]

## The two library variants

The flake exposes two builds:

- **`packages.default`** — CPU build. `doEval()` runs an OpenMP-parallelized
  loop on the host; the backend dispatches through
  `RooFit::EvalBackend::Cpu()`.
- **`packages.cuda`** — CUDA build. Each custom PDF ships a `__global__` kernel
  under `gpu/`; `doEval()` queries the output buffer with
  `cudaPointerGetAttributes` and, when RooFit hands it device memory, launches
  the kernel instead of running the host loop. The backend dispatches through
  `RooFit::EvalBackend::Cuda()`. Built against a `-Dcuda=ON` override of ROOT
  (re-exposed as `packages.rootCuda` so downstream consumers don't have to
  duplicate the override).

## Build details

The project is CMake-based. The CUDA build sets `-DAU_USE_CUDA=ON`, which
enables CUDA as a CMake language, picks up `gpu/*.cu` alongside `src/*.cpp`,
defines `AU_ROOFIT_BACKEND_CUDA=1` for the host code, and links `CUDA::cudart`.
The CPU build never compiles the kernel files and never touches the CUDA
toolchain.

## Opting in from a downstream flake

Swap `default` → `cuda` and `pkgs.root` → `utils.packages.${system}.rootCuda`,
and export the matching compile-define so downstream code that calls
`BestAvailableBackend()` agrees with the prebuilt library:

```nix
analysis-utils = utils.packages.${system}.cuda;
root           = utils.packages.${system}.rootCuda;
# ...
shellHook = ''
  export NIX_CFLAGS_COMPILE="-DAU_ROOFIT_BACKEND_CUDA=1''${NIX_CFLAGS_COMPILE:+ $NIX_CFLAGS_COMPILE}"
  export LD_LIBRARY_PATH="/run/opengl-driver/lib''${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
'';
```

The `LD_LIBRARY_PATH` line ensures nix-built binaries find the host's
`libcuda.so` (NixOS exposes it under `/run/opengl-driver/lib`) rather than a
stub from the nix store.

@note The `cudaCapabilities` list in the flake's
`pkgs = import nixpkgs { config = { ... }; }` block defaults to `[ "8.9" ]`
(Ada / RTX 40-series, e.g. RTX 4090). Adjust to your GPU's compute capability
and update `-DCMAKE_CUDA_ARCHITECTURES` in the `rootWithCuda` overlay to match.

@warning Do not override the `nixpkgs` input (e.g. via
`inputs.utils.inputs.nixpkgs.follows`). The CUDA-built ROOT is large and would
otherwise be rebuilt from source against your nixpkgs revision. Leaving the
input alone lets your machine pull the pre-built artifact from the
@ref installation "binary cache" instead.

## Performance, measured

On a 4 M-event unbinned extended maximum-likelihood fit (73mGe calibration
workload, RTX 50-series GPU):

| Path | Single fit | Simultaneous fit |
|------|-----------:|-----------------:|
| Scalar `evaluate()` (legacy) | ~minutes | ~minutes |
| OpenMP batched `doEval` | ~2 min | ~3 min |
| CUDA kernels (`packages.cuda`) | ~6 s | ~36 s |

The CUDA win for simultaneous fits is smaller than for single fits because once
`doEval` is effectively free, the remaining time is dominated by Minuit2 itself,
the analytic normalization integrals (computed on the host per-step), and
RooFit's evaluator orchestration — none of which benefit from GPU work.

@note The table above was measured on a 50-series (Blackwell) GPU, while the
flake's default `cudaCapabilities` targets 8.9.
