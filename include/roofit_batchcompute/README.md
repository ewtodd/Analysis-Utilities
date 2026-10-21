# roofit_batchcompute — vendored ROOT batch-compute headers

Pinned verbatim copies of three headers from ROOT's `roofit/batchcompute`
tree, taken from the exact version this project builds against
(ROOT v6-40-00, the version pinned by flake.nix's nixpkgs input):

- `RooBatchCompute.h` — from `roofit/batchcompute/res/`. Defines
  `RooBatchCompute::Config`, whose `useCuda()` / `cudaStream()` the photopeak
  PDFs' `doEval` reads from the `EvalContext` so their CUDA kernels launch on
  ROOT's per-node stream.
- `CudaInterface.h` — from `roofit/batchcompute/src/`. Defines
  `RooBatchCompute::CudaInterface::CudaStream`. Its `operator cudaStream_t()`
  is declared only under `__CUDACC__`, so the conversion happens in the
  `gpu/*.cu` launchers, never in host code.
- `RooBatchComputeTypes.h` — from `roofit/batchcompute/res/`. The
  `__roodevice__` / `__roohost__` / `__rooglobal__` macro shim.

ROOT's nix package does not install its batch-compute headers, so without
these the CUDA build cannot name the `Config` / `CudaStream` types. They are
copied verbatim and should not be edited; bump them together with the pinned
ROOT version. Their only ROOT dependencies (`<ROOT/RSpan.hxx>`,
`<DllImport.h>`) are installed by the normal ROOT build.
