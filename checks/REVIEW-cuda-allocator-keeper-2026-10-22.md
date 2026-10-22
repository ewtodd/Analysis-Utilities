# CUDA allocator workaround, simultaneous coefficient ranges and CUDA regression check — review checklist, 2026-10-22
<!---->
Branch `perf/cuda-allocator-keeper`, based on `main` at `73a8384`.
Not merged and not pushed.
It contains two core changes and one tooling addition.
<!---->
1.
CUDA allocator workaround.
ROOT 6.40's CUDA batch-compute backend runs `cudaMalloc` / `cudaFree` on a small
descriptor in every kernel dispatch and NLL reduction.
From about 250k events per graph node, no other small device allocation is live
between those calls, so each pair maps and unmaps a fresh driver chunk (about
0.9 ms on driver 595.104.02 with an RTX 4090).
That made CUDA fits slower than CPU fits: 83 s against 67 s for a single 4M-event
fit, and 612 s against 236 s for a 4 × 1M simultaneous fit.
Holding one small allocation for the life of the process keeps the chunk mapped.
With it, the same fits take 7.5 s and 32 s, and the GPU is 88–93% busy.
ROOT master fixes the cause upstream (PR #23263); the workaround is harmless once
the pinned ROOT includes it.
<!---->
2. Simultaneous fits fix each channel's coefficient range.
In an extended fit with `Range()`, RooFit adds a projection ratio per component
that should be 1.
It is computed as the normalised component evaluated at `x`'s
current value, so a `RooGaussian` that underflows there makes it 0/0 and the NLL
NaN.
With channel windows far enough apart that is unavoidable, and `main`'s
`FitSimultaneous` then fails with "all function calls during minimization gave
invalid NLL values".
Telling each channel's `RooAddPdf` which range its yields refer to (the range
its fit normalises over) removes that ratio.
The NLL and the fit result are unchanged wherever `main` works. Per-channel
curves, CSV dumps and reduced chi² now use the window normalisation the fit
used, so those reported numbers change.
<!---->
3. A `nix run .#cuda-regression` check that would have caught the slowdown.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
### Backend selection holds the allocation
this commit
<!---->
- [x] `include/RooFitUtils.hpp:43-67` — new public function
      `KeepCudaAllocatorChunkMapped()` and its Doxygen block.
      Check: the block states the side effects a caller needs to know: 256 bytes
      of device memory held for the process lifetime, CUDA context creation on the
      first call, and a no-op in the CPU build.
- [x] `include/RooFitUtils.hpp:85-92` — the CUDA branch of
      `BestAvailableBackend()` calls it before returning
      `RooFit::EvalBackend::Cuda()`. The CPU branch is unchanged.
      Check: every library fit selects its backend here (`RooFitUtils.cpp:572`,
      `:3005`, `:3102`, `InteractiveRooFitEditor.cpp:801`,
      `InteractiveSimultaneousFitEditor.cpp:875`). So the allocation exists before
      any likelihood is built, including for models with no custom tails, which is
      the default (`use_step_` and the tail flags start `kFALSE`).
- [x] `src/CudaAllocatorKeeper.cpp:1-29` — `AllocateChunkKeeper()` runs
      `cudaMalloc(256)`. On failure it clears the runtime's last-error state and
      returns null.
      `KeepCudaAllocatorChunkMapped()` stores the result in a
      function-local `static`, so it is allocated once, initialised thread-safely,
      and never freed.
      Without `AU_ROOFIT_BACKEND_CUDA` the function body is empty.
      Check: a failed allocation cannot surface later as an unrelated CUDA error,
      and nothing frees the pointer.
- [x] The behaviour change underneath: before this branch, every `cudaFree` inside
      ROOT's `compute()` / `reduceNLL()` / `reduceSum()` also synchronised the
      whole device. With the chunk kept mapped those frees no longer reach the
      driver, so that implicit synchronisation disappears.
      Check: results do not depend on it. Evidence under "Test evidence": CPU and
      CUDA agree to about 2e-16, and full fits match exactly, with and without the
      workaround.
<!---->
### Simultaneous fits fix each channel's coefficient range
this commit
<!---->
- [x] `include/RooFitUtils.hpp:379` — private declaration of
      `FixChannelCoefRanges(Bool_t split_range)`.
- [x] `src/RooFitUtils.cpp:2874-2893` — for every channel it `dynamic_cast`s
      `sim_channel_pdfs_[name]` to `RooAddPdf` and calls `fixCoefRange` with either
      the channel's `fitrange_<name>` window (`split_range`) or the shared
      `kFitRangeName`.
      Check: each channel sum is a `RooAddPdf` (built at `:2399-2400`), so the
      cast never skips a channel.
- [x] `src/RooFitUtils.cpp:2994` — refit after loading saved parameters: windows,
      immediately before the `SplitRange` fit at `:2997-3005`.
- [x] `src/RooFitUtils.cpp:3041-3042` — interactive editor branch: the shared
      `kFitRangeName`, before the editor opens (`:3046`). The editor's `DoRefit`
      fits over that shared range without `SplitRange`
      (`InteractiveSimultaneousFitEditor.cpp:870-875`), and draws with `x_`'s
      default range set to the same bounds (`:819-820`, `:928-929`).
      Check: the editor's fits and curves are unchanged, because its yields
      already referred to the shared range.
- [x] `src/RooFitUtils.cpp:3071` — non-interactive fit: windows, before both the
      `fit_debug_` seed NLL and the `SplitRange` fit at `:3095-3102`.
- [x] The loaded-parameters path without `refit_after_load_` sets nothing, so it
      behaves as before.
      Check: if an earlier `FitSimultaneous` call in the same process set the
      windows, they stay set for this path.
      That is only relevant if a process
      calls `FitSimultaneous` more than once. 
      ewt: I don't think this is something that happens? 
      son-of-anton-bot: confirmed it cannot. `FitSimultaneous` rebuilds every
      channel sum through `BuildChannelModel` (`:2911`) on each call, so a new
      `RooAddPdf` never inherits a coefficient range from an earlier call.
- [x] Reported numbers that change: the per-channel chi² (`:2535-2565`),
      `PlotChannel` (`:2738+`) and `DumpChannelCSV` (`:2787+`) all draw
      `expectedEvents(&nset) * getVal(&nset)` with `nset` over `x_`'s default
      range, which is the union of all channel windows (`:2927-2928`). With the
      coefficient range fixed, each channel's curve is normalised over its own
      window, which is what its fit used.
      Before, components reaching outside the
      window (background, tails, step) were drawn too low by
      `∫window / ∫union`.
      Check: this correction is wanted.
      On a toy channel with 200,000 events in its
      window, the curve integrates to 193,299 before and 200,000.0003 after.
      In
      the end-to-end run below the reduced chi² moves from 6.20 / 17.65 / 22.17 to
      3.72 / 17.01 / 9.40.
      Simultaneous fits whose windows nearly equal the union
      barely move.
<!---->
## Tooling
<!---->
### Regression check
this commit
<!---->
- [x] `tests/CudaRegression.cpp:45-64` — limits and reference values: agreement
      to 1e-10 relative at 40 points; start NLL within 1e-9 of the reference; a
      CUDA-to-CPU time ratio of at most 0.3 (single) and 0.6 (simultaneous).
      Check: the 0.6 simultaneous limit reflects ROOT 6.40 recomputing every
      channel on CUDA, and is to be tightened after the ROOT master pin.
      The
      reference values were recorded with ROOT 6.40.00.
- [x] `tests/CudaRegression.cpp:81-151` — `BuildChannel()`: one toy channel
      (Gaussian, all four custom tails, `RooPolynomial` background, six yields)
      with 1M events from a fixed-seed generator.
      Check: the seeds are fixed, so the reference NLL values are reproducible.
- [x] `tests/CudaRegression.cpp:153-193` — `BuildModel()`: a shared `x` with a
      `fitrange_<channel>` window per channel and a union `fitrange`.
      Lines 182-183
      fix each simultaneous channel's coefficient range to its window, as
      `FitSimultaneous` now does.
      Check: the reference values did not change when this replaced the earlier
      "park `x` mid-range" workaround, and ROOT no longer prints
      `checkForDup` errors.
- [x] `tests/CudaRegression.cpp:195-202` — `CreateNll()` uses
      `Extended + Range("fitrange")`, plus `SplitRange` for the simultaneous model,
      mirroring `RooFitUtils.cpp:3095-3102`.
- [x] `tests/CudaRegression.cpp:262-301` — `CheckModel()` builds the CPU NLL with
      `EvalBackend::Cpu()` and the CUDA NLL with `BestAvailableBackend()` (line
      269), so the check exercises the workaround the way library fits do.
      Check: the three reported checks are the ones listed in `doc/pages/gpu.md`.
- [x] `flake.nix:111-147` — `cudaRegression` builds `tests/` against `utilsCuda`
      and `rootWithCuda`. It wraps the binary with `/run/opengl-driver/lib` (host
      `libcuda`) and with `ROOT_INCLUDE_PATH` (cling finds the PDF dictionary
      header).
      Check: the link line uses `root-config --glibs`, not `--libs`, because
      `libanalysis-utils.so` leaves its GUI symbols to the executable (see Nits).
- [x] `flake.nix:238-243` — `packages.cudaRegression` and
      `apps.cuda-regression`.
      Check: it is an app and not a `checks` entry because the build sandbox has
      no GPU.
<!---->
### Documentation
this commit
<!---->
- [x] `doc/pages/gpu.md:73-131` — the performance table is replaced with the
      measurements from this branch, with a new section on the allocator
      workaround and a new "Regression check" section. The old 6 s / 36 s figures
      are kept in one paragraph, marked as a different workload, GPU and driver.
      Check: the numbers match "Test evidence" below, and no root-absolute links
      were added.
<!---->
## Test evidence
<!---->
- [x] `nix flake check --no-build` passes. `nix build .#default`, `.#docs` and
      `.#cudaRegression` (and so `.#cuda`) succeed, and the docs warning log is
      empty.
- [x] `nix run .#cuda-regression` on this branch exits 0, with no `checkForDup`
      output. Agreement is 1.9e-16 (single) and 2.2e-16 (simultaneous). Both
      reference values match. The CUDA-to-CPU ratios were 0.101–0.106 and
      0.400–0.412 across four runs.
- [x] The same `tests/CudaRegression.cpp`, compiled against `main`'s CUDA library
      (`73a8384`, no workaround), exits 1.
      Agreement and reference values pass;
      both speed checks fail at ratios 1.907 and 7.863.
      Check: the check catches exactly this regression and nothing else.
- [x] End-to-end `RooFitUtils::FitSimultaneous` (scratch driver: 3 channels,
      200k events each, ±20 windows, flat background, step, low and high
      exponential tails; CUDA build, non-interactive):
      - windows 40 apart (union three times a window): `main` fails with "all
        function calls during minimization gave invalid NLL values" (status 302,
        `min_nll` 0, every μ left at its start value).
        This branch fits:
        `min_nll` −4435984.272, μ 40.069 / 80.077 / 120.077.
      - windows 5 apart (union close to each window): `min_nll` −4435983.549 and
        μ 40.06528854 / 45.06578445 / 50.07823181 on both, identical in every
        printed digit.
        Reduced chi² 6.20 / 17.65 / 22.17 on `main`, 3.72 / 17.01
        / 9.40 on this branch.
      Both runs end with Migrad status 3 (EDM 64 on the first), because the
      driver's event generator does not match the model's shapes exactly.
      Check: the fit result is unchanged wherever `main` works.
- [x] Coefficient-range toy experiments (scratch, CPU): the NLL is unchanged by
      `fixCoefRange` with `SplitRange` (−4358073.87483) and with one shared
      range (−4294351.33316).
      A mismatched coefficient range (the windows, in a
      shared-range fit) brings the duplicate projection nodes back without
      changing the NLL; that is why each fit path sets the range matching its own
      normalisation.
- [x] Full Migrad fits of the regression toy (scratch benchmark, workaround
      applied through a direct `cudaMalloc` in the benchmark rather than this
      function):
      - single 4M: CPU 67.5 s, CUDA 83.1 s, CUDA + workaround 7.5 s.
      - simultaneous 4 × 1M: CPU 236 s, CUDA 612 s, CUDA + workaround 31.6 s.
      The simultaneous fits reach the same `minNll` (−34117699.1986) in the same
      7354 evaluations on all three paths, with and without `fixCoefRange`.
      The
      two single CUDA fits match each other exactly; the single CPU fit stops at a
      slightly different point (`minNll` differs by 0.017 after 1617 instead of
      1947 evaluations).
      The toy's single fit ends with Migrad status −1 on every
      path.
      Check: the CUDA results are unchanged by the workaround.
<!---->
## Checked and excluded
<!---->
- The single-peak model (`total_pdf_`) is not changed. Every path sets `x_`'s
  default range to the same bounds as `kFitRangeName` (for example `:863-864`).
  In a toy with that setup, a σ = 0.3 peak and `x` left 65σ away, the NLL is
  finite and identical with and without `fixCoefRange`.
- The `checkForDup` errors never reached production logs: the library removes
  the `InputArguments` topic from every RooFit stream (`RooFitUtils.cpp:69-71`).
  They appeared only in the standalone regression check, and are gone with the
  fixed coefficient ranges.
- The 1 ms polling sleep in ROOT's evaluator was measured and essentially never
  fires on these graphs; it is not a cause of the slowdown.
- Our own kernels take 0.1–0.14 ms each at 1M events, and memcpy time is
  negligible. Neither is a cause.
- `Range` + `SplitRange` (from `391e96b`) costs nothing measurable: 76 against
  78 ms per evaluation before the workaround. `fixCoefRange` does not change
  speed either (4.21 against 4.23 ms per CUDA evaluation).
<!---->
## Nits
<!---->
- Pre-existing: `libanalysis-utils.so` does not record ROOT's GUI library as a
  dependency, so an executable linking it must add it itself (with
  `root-config --libs` instead of `--glibs`, the link fails on undefined `TG*`
  symbols).
  That is why the regression app uses `--glibs`.
  Left unchanged; this
  branch does not touch the library's link line.
- `2ce3aaa`'s message describes the 0/0 as "the component's integral over the fit
  range".
  It is actually the normalised component evaluated at `x`'s current
  value.
  Its 1e-300 floors still help, but they cover only the custom tails, not
  `RooGaussian`.
