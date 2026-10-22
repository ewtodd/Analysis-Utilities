# CUDA kernels on ROOT's per-node stream — review checklist, 2026-10-22
<!---->
Branch `fix/cuda-per-node-stream`, based on `main` at `7515835`.
Not merged and not pushed.
It holds commit `074bd37` plus follow-up changes that were uncommitted when this was
written.
The photopeak PDFs' CUDA kernels used to launch on the legacy default stream; they
now launch on the stream RooFit's evaluator assigns to each graph node, as ROOT's
own batch-compute kernels do.
That needs two of ROOT's private headers. `074bd37` committed hand-copied
(and reformatted) versions of them; the follow-up deletes those copies, and
`flake.nix` instead extracts them from the source tarball the CUDA ROOT is built
from.
The kernel math and the CPU path are unchanged.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. ROOT line numbers
refer to tag `v6-40-00`. Checks indicate review by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
### Stream selection in `doEval`
`074bd37`, simplified in follow-up
<!---->
- [x] `src/RooFitPhotopeakPdfs.cpp:8-11` — inside `#ifdef AU_ROOFIT_BACKEND_CUDA`,
      includes `<roofit_batchcompute/RooBatchCompute.h>` for the full
      `RooBatchCompute::Config` type that `ctx.config(this)` returns by value.
      Check: the CPU build never sees this include, and the `Config` it defines is
      the one `libRooFitCore` was compiled with. The header comes from
      `rootCuda`'s own source tarball (`flake.nix:43-58`).
- [x] `src/RooFitPhotopeakPdfs.cpp:168-176` (`RooStepShelf`),
      `src/RooFitPhotopeakPdfs.cpp:249-257` (`RooLowExpTail`),
      `src/RooFitPhotopeakPdfs.cpp:329-337` (`RooLowLinTail`),
      `src/RooFitPhotopeakPdfs.cpp:410-418` (`RooHighExpTail`) — the device-pointer
      branch reads `ctx.config(this).cudaStream()` into `node_stream` and passes it
      as the launcher's new last argument.
      Check: at each site the other arguments, the `IsDevicePointer` gate, the early
      `return` and the host loop below are unchanged from `main`.
- [x] The same four sites, read against ROOT's
      `roofit/roofitcore/src/RooFit/Evaluator.cxx`.
      Lines 225-231 give every node its
      own stream and store it in that node's `Config`. Lines 558-561 make that stream
      wait on each input node's event, line 576 calls `doEval`, and line 577 records
      the node's event on the same stream. Line 495 marks a node finished only when
      `cudaStreamIsActive` on that stream is false.
      Check: a kernel launched on this stream therefore cannot start before its
      inputs are written, and its clients cannot start before it finishes. This is
      the ordering every fit result on the CUDA backend now depends on.
- [x] Same sites — `ctx.config(this)` is keyed on `this`, the node being evaluated.
      Check: the stream belongs to the PDF component whose output buffer is being
      written, not to its parent `RooAddPdf` or to the NLL.
- [x] Same sites — `074bd37` wrote `cfg.useCuda() ? cfg.cudaStream() : nullptr`;
      the follow-up passes `cudaStream()` directly.
      Check: these are the same value, because ROOT's `Config::useCuda()` is defined
      as `_cudaStream != nullptr` (`roofit/batchcompute/res/RooBatchCompute.h:57`).
<!---->
## Tooling
<!---->
### Kernel launchers
`074bd37`, comments trimmed in follow-up
<!---->
- [x] `gpu/RooHighExpTail.cu:39-48`, `gpu/RooLowExpTail.cu:39-48`,
      `gpu/RooLowLinTail.cu:26-35`, `gpu/RooStepShelf.cu:29-37` — each launcher
      takes a `CudaStream *`, converts it to `cudaStream_t` (null becomes stream 0)
      and launches with `<<<blocks, threads, 0, s>>>`.
      Check: the kernel bodies above each launcher are untouched, and the only
      change to the launch is the added shared-memory size `0` (the previous
      default) and the stream.
- [x] Line 5 of each `gpu/*.cu` — `#include <roofit_batchcompute/CudaInterface.h>`.
      Check: this is the only file that defines `CudaStream`'s conversion to
      `cudaStream_t`, and it is compiled only by nvcc.
<!---->
### Kernel header documentation
follow-up
<!---->
- [x] `include/RooFitPhotopeakKernels.hpp:6-32` — `@file` block restored from
      `main`, with its `@warning` (device-resident pointers). Its `@note` is
      rewritten to say that ordering now comes from the node's stream.
      Check: the note matches the evaluator behaviour described under core above.
- [x] `include/RooFitPhotopeakKernels.hpp:34-44` — forward declaration of
      `CudaStream`, wrapped in `@cond`/`@endcond`.
      Check: this installed header needs no ROOT-private header, and the published
      docs carry no `RooBatchCompute` namespace pages (none in the follow-up docs
      build).
- [x] `include/RooFitPhotopeakKernels.hpp:46-143` — the four launcher declarations
      in `/** */` blocks, back in `main`'s order, with `main`'s parameter wording
      plus the new `@param stream` and a one-line formula for each kernel.
      Check: each formula against its kernel: `gpu/RooLowExpTail.cu:17-35`,
      `gpu/RooStepShelf.cu:9-25`, `gpu/RooLowLinTail.cu:9-22`,
      `gpu/RooHighExpTail.cu:17-35`.
      `074bd37` had the high-energy tail's `erfc`
      argument without its minus sign; the kernel computes `z = mu - x`.
<!---->
### ROOT's private headers, from the flake
follow-up (replaces the copies `074bd37` committed under `include/roofit_batchcompute/`)
<!---->
- [x] `flake.nix:43-58` — `rootBatchComputeHeaders` unpacks
      `roofit/batchcompute/res/RooBatchCompute.h` and
      `roofit/batchcompute/src/CudaInterface.h` from `rootWithCuda.src` into
      `$out/include/roofit_batchcompute/`.
      The member paths are spelled out, so a
      ROOT version that moves either file fails this derivation.
      Check: `rootWithCuda.src` is the tarball the CUDA ROOT is compiled from.
      `pkgs.root` has no `patches`, and its `postPatch` does not touch `roofit/`, so
      the unpacked files are exactly what that build used.
      The output was compared
      with `cmp` against the files at tag `v6-40-00` on GitHub: byte-identical.
- [x] `flake.nix:91-95` — the CUDA package passes
      `-DAU_ROOT_BATCHCOMPUTE_INCLUDE_DIR` pointing at that derivation.
      Check: the CPU package gets no such flag, and `rootWithCuda` itself is
      unchanged, so the CUDA ROOT is not rebuilt.
- [x] `flake.nix:260-261` — the CUDA dev shell exports
      `AU_ROOT_BATCHCOMPUTE_INCLUDE_DIR` and adds it to `CPLUS_INCLUDE_PATH`.
      Check: this is what lets editors (clangd) and direct compiles find the headers,
      just as `$PWD/include` did for the copies.
      `gpu/.clangd` is deliberately left
      without a `/nix/store` path, so entering the shell does not modify that
      tracked file.
- [x] `CMakeLists.txt:45-61` — with `AU_USE_CUDA`, the cache variable defaults to
      the environment variable, and configuring stops with a clear message if
      either header is missing beneath it.
      Check: the message names both source locations in ROOT's tree. It was tested by
      building the CUDA package with the flag filtered out.
- [x] `CMakeLists.txt:100` and `CMakeLists.txt:120` — the directory is added
      `PRIVATE` to the shared and static targets.
      Check: it does not propagate to consumers, and nothing under
      `roofit_batchcompute/` is installed.
- [x] `doc/pages/gpu.md:31-40` — "Build details" now explains where the two headers
      come from, and what a CUDA build outside the flake must provide.
      Check: no root-absolute links were added.
<!---->
## Test evidence
<!---->
- [x] `nix flake check --no-build` passes. `nix build .#default`, `.#cuda` and
      `.#docs` all succeed on the follow-up working tree, and the docs build's
      warning log is empty.
- [x] GPU test on the RTX 4090: an extended 1,000,000-event fit of a Gaussian plus
      all four custom PDFs on `x` in [0, 200].
      The CUDA NLL was compared with the CPU
      NLL at 400 parameter points; the maximum difference was 1.9e-9 on |NLL| of
      about 9.5e6.
      Check: a stream-ordering race would show up here as a mismatch.
- [x] Same test, Minuit2 Migrad fits: three CUDA fits are identical to each other
      and to `main`'s CUDA build in every printed digit (12 significant figures).
      That held for `074bd37` itself, for the copied-header follow-up, and for the
      final flake-header build. They agree with the CPU fit to within 1e-8 relative
      (largest: `nShelf` at 6e-9 and `slope` at 4e-9).
      The fit status is 3 on every backend and on `main`, so that comes from the
      toy model and not from this change.
      Check: the result is unchanged from `main`, which is the claim for a change
      that only moves work between streams.
- [x] Wall time is unchanged: about 16.0–16.7 s per CUDA fit on both this branch
      and `main`.
      Check: the speed-up the commit message claims for simultaneous fits was not
      measured. No simultaneous fit was run.
- [x] In the CUDA dev shell, the CUDA path of `src/RooFitPhotopeakPdfs.cpp` compiles
      (`g++ -fsyntax-only`) with the headers found only through the shell's
      `CPLUS_INCLUDE_PATH`.
<!---->
## Checked and excluded
<!---->
- The old default-stream launches were not a race.
ROOT creates node streams with
  `cudaStreamCreate` (blocking), and blocking streams synchronise implicitly with the
  legacy default stream.
  This branch changes how work is scheduled; it does not fix
  wrong results on `main`.
- Scalar outputs (`nOut == 1`, Evaluator.cxx:566-568) are host buffers, so
  `IsDevicePointer` is false and the stream is never read.
- `Doxyfile.in` is back to `main`'s version.
The `EXCLUDE` that `074bd37` added
  pointed at the deleted copies.
- The `Makefile` build is CPU-only and never reads these headers.
- Pre-existing and not touched here: neither dev shell provides `cmake`.
With
  `cmake` added, a CUDA configure in the dev shell already fails earlier at
  `find_package(CUDAToolkit)` (`CMakeLists.txt:42`, "missing: CUDA_CUDART"), before
  reaching the new check.
  CUDA builds go through `nix build .#cuda`.
<!---->
## Nits
<!---->
- None open.
The three found during the first review pass are fixed in the
  follow-up: the comment repeated in each `gpu/*.cu` launcher (which also mentioned
  an impossible "CPU-built" case) is gone; the redundant `useCuda()` ternary is
  gone; and the unused `RooBatchComputeTypes.h` is no longer copied.
