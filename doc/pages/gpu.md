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

The kernels launch on the CUDA stream RooFit's evaluator assigns to each node
of the computation graph, which the custom PDFs read through two of ROOT's
private batch-compute headers, `RooBatchCompute.h` and `CudaInterface.h`. No
ROOT install ships them, so the flake extracts both from the same source
tarball `rootCuda` is built from, and passes their directory to CMake as
`AU_ROOT_BATCHCOMPUTE_INCLUDE_DIR`; the CUDA dev shell exports it too. A CUDA
build outside the flake must set that variable to a directory holding
`roofit_batchcompute/RooBatchCompute.h` (from `roofit/batchcompute/res/`) and
`roofit_batchcompute/CudaInterface.h` (from `roofit/batchcompute/src/`) taken
from the source of the ROOT it links against.

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

Unbinned extended maximum-likelihood fits (Minuit2 Migrad, strategy 1) of the
toy photopeak model in `tests/CudaRegression.cpp`: a Gaussian, all four custom
tails and a linear background per channel. Measured on an RTX 4090 with NVIDIA
driver 595.104.02 and ROOT 6.40.00:

| Path | Single fit, 4 M events | Simultaneous fit, 4 × 1 M events |
|------|-----------:|-----------------:|
| OpenMP batched `doEval` (`packages.default`) | 67 s | 236 s |
| CUDA, without the allocator workaround | 83 s | 612 s |
| CUDA (`packages.cuda`) | 7.5 s | 32 s |

All three paths reach the same minimum. With the workaround the GPU is 88–93%
busy during the fit; without it, under 10%.

### The allocator workaround

ROOT 6.40's CUDA batch-compute backend allocates and frees a small descriptor
buffer in every kernel dispatch and NLL reduction. Once RooFit's per-node
buffers are larger than the driver's small-allocation chunk (from roughly 250k
events per node), no small allocation stays live between those calls. Each
`cudaMalloc` / `cudaFree` pair then maps and unmaps a fresh chunk, about 0.9 ms
on the hardware above, and that dominated the time per likelihood evaluation.
KeepCudaAllocatorChunkMapped() holds one 256-byte allocation for the life of
the process, which keeps the chunk mapped. BestAvailableBackend() calls it, so
fits set up through it, including every fit in this library, need nothing
more. Code that passes `RooFit::EvalBackend::Cuda()` itself should call
KeepCudaAllocatorChunkMapped() once before building the likelihood.

The simultaneous fit still has a smaller gap. ROOT 6.40's CUDA path recomputes
every channel on every evaluation, while the CPU path skips unchanged ones.
ROOT master fixes both problems: PR #23263 removes the per-call allocations and
PR #23753 skips clean nodes on the CUDA path.

The earlier figures of about 6 s and 36 s came from the 73mGe calibration
workload on an RTX 50-series GPU with ROOT 6.38.02. That is a different
workload, GPU and driver, so they are not directly comparable with the table.

## Regression check

`nix run .#cuda-regression` builds the toy models above with 1 M events per
channel and checks, for a single-channel and a four-channel simultaneous model,
that:

- the CPU and CUDA backends agree on the NLL to 1e-10 relative at 40 parameter
  points;
- the CPU NLL at the start point matches a recorded reference to 1e-9
  relative, which catches value changes made by both backends at once, as from a
  ROOT bump;
- the CUDA backend's median time per evaluation is below a fixed fraction of the
  CPU backend's (0.3 single, 0.6 simultaneous), so the check does not depend on
  the hardware.

It needs a GPU, which `nix flake check`'s sandbox does not have, so run it by
hand after changing the PDFs, the kernels, the flake's ROOT or nixpkgs, or the
NVIDIA driver. It exits non-zero if any check fails. After a deliberate change
to the model or the PDFs, update the reference values in
`tests/CudaRegression.cpp` with the ones it prints.
