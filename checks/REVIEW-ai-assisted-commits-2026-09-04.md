# AI-assisted commits — review checklist, 2026-09-04
<!---->
Per the AI-assisted-development disclosure section of the README, parts of this
codebase (interactive GUI editor, routine boilerplate) were written with the
help of Claude Code and/or son-of-anton, but **all changes to the core analysis
logic (fitting models, signal processing, result extraction) are human-reviewed
and approved before being committed.** This doc is that review.
Line numbers are as of HEAD (`63f7fbc`).
No code is reproduced here; open the
file at the line and review.
Checks indicate review by ewt.
<!---->
AI-assisted commits (identified by the `Co-Authored-By: Claude` trailer),
newest first:
<!---->
| Commit | Title | Model | Date |
| --- | --- | --- | --- |
| `78e0a0b` | fix(RooFitUtils): saved state must not override model configuration on single fits | Claude Opus 5 | 2026-08-27 |
| `8bbdb1b` | feat(RooFitUtils): constrain a peak separation to an external value | Claude Opus 5 | 2026-08-27 |
| `88cec44` | feat(RooFitUtils): cap tail decay lengths; saved files carry values, not constness | Claude Opus 5 | 2026-08-27 |
| `0d8711d` | fix(RooFitUtils): saved state must not override model configuration | Claude Opus 5 | 2026-08-27 |
| `3523967` | fix(RooFitUtils): apply the EDM validity criterion to the non-interactive path | Claude Opus 5 | 2026-08-27 |
| `26ff1b4` | fix(RooFitUtils): fail cleanly on unresolvable parameter links; relax edm gate | Claude Opus 5 | 2026-08-26 |
| `48bb04f` | feat(RooFitUtils): let saved simultaneous-fit params seed a real fit | Claude Opus 5 | 2026-08-26 |
| `7d4360a` | fix: NaN in simultaneous photopeak fit from zero-integral components | Claude Opus 4.8 (1M context) | 2026-06-01 |
| `222f155` | fix: invalid NLL in simultaneous photopeak fits at high statistics | Claude Opus 4.8 (1M context) | 2026-05-31 |
| `08c983e` | feat: per-peak step/shape locking for sim fits + X11 teardown guard | Claude Opus 4.8 (1M context) | 2026-05-31 |
| `abba503` | fix: segfault on interactive editor teardown | Claude Opus 4.7 (1M context) | 2026-05-25 |
| `52a4f4a` | perf: GPU-accelerated RooFit photopeak fitting | Claude Opus 4.7 | 2026-05-19 |
| `92df438` | perf: bulk-read fast path for CoMPASS binary conversion | Claude Opus 4.7 | 2026-05-19 |
| `0da3794` | tail decay: parametrize as ratio over sigma, enforce tau >= sigma | Claude Opus 4.7 (1M context) | 2026-05-11 |
| `71bd648` | roofit save: 17-digit precision, persist Hesse errors; silence eval-error spam | Claude Opus 4.7 (1M context) | 2026-05-11 |
| `ebaf64c` | switch RooFitUtils to unbinned likelihood; live-rebin display | Claude Opus 4.7 (1M context) | 2026-05-11 |
| `bf45a6e` | add interactive simultaneous fit editor | Claude Opus 4.7 (1M context) | 2026-05-10 |
| `703179e` | sim fit: shared union range, tighter sigma bounds | Claude Opus 4.7 (1M context) | 2026-05-10 |
| `df06f5a` | analytical integrals for step+tail PDFs | Claude Opus 4.7 (1M context) | 2026-05-10 |
| `231f577` | add simultaneous fit support to RooFitUtils | Claude Opus 4.7 (1M context) | 2026-05-10 |
| `db3c39c` | add RooFit backend (RooFitUtils + InteractiveRooFitEditor), drop python sim fit | Claude Opus 4.7 (1M context) | 2026-05-10 |
| `62217b4` | document hyper-EMG peak shape in README | Claude Opus 4.6 (1M context) | 2026-04-11 |
| `ea8a1da` | add hyper-EMG peak shape for CZT gamma-ray spectroscopy | Claude Opus 4.6 (1M context) | 2026-04-11 |
| `9649e02` | refactor: amplitude params as ratios, fix constrained fits and canvas cleanup | Claude Opus 4.6 (1M context) | 2026-03-26 |
| `20e8672` | fix: interactive editor segfault on close and ignoring use_* flags | Claude Opus 4.6 (1M context) | 2026-03-26 |
<!---->
The seven newest commits (`78e0a0b` … `48bb04f`) all touch the RooFit fit
backend's saved-state handling, EDM validity gate, and constraint machinery —
i.e. core result-extraction logic — and deserve line-by-line review. The rest
are the older foundations: the unbinned-likelihood rewrite, the custom PDFs and
their analytical integrals, the simultaneous-fit machinery, the GPU kernels,
the binary-reader fast path, and the GUI teardown fixes.
<!---->
## Core analysis logic — review these lines
<!---->
### Saved state must not override model configuration (single + simultaneous)
`78e0a0b`, `0d8711d`
<!---->
- [ ] `src/RooFitUtils.cpp:846-924` (`LoadInteractiveParams`) — the single-fit
      `.roofits` loader now matches saved parameters **by NAME** and refuses to
      let saved constness override the live model. Check: confirm a saved
      `.roofits` written when a component was OFF cannot re-enable it, and that
      a component ON in the model but missing from the file keeps its live
      (unfixed) value rather than being zeroed.
- [ ] `src/RooFitUtils.cpp:981-1060` (`LoadSimInteractiveParams`) — the
      simultaneous loader applies the same by-name / no-constness-override rule
      (`0d8711d`), then `48bb04f`'s locks run after the load.
      Check: the
      comment at `:1060` says locks run after the load — confirm the load
      cannot leave a parameter fixed that a later `lock_shape_after_seed`
      should free, and vice-versa.
- [x] `src/RooFitUtils.cpp:867-869` — comment notes the single path "never got"
      the fix the sim path had.
      Check: the two loaders now share the invariant
      (name-match, values-only).
      Diff them side by side and confirm no
      single-path-only branch still applies a saved `isConstant()`.
NOTE: While looking at this, the below were also noted and fixed:
      1.
      `LoadSimInteractiveParams` returned `kTRUE` after matching **zero**
         parameters.
         A `.simroofits` whose names had all drifted — a renamed
         channel, a changed peak count — printed "Loaded sim interactive params
         from …" and silently seeded nothing, which is the same class of quiet
         failure `78e0a0b` fixed on the single path, moved from the constness
         path to the arrival path.
         The behavior is fixed to be identical to the single loader: an
         `n_set == 0` guard warns and returns `kFALSE`, which drops the caller
         at `:2940` through to the interactive editor exactly as a missing file
         does.
      2.
      `kFitRangeName` was **private** in `include/RooFitUtils.hpp` and used
         in only 7 places, while the raw `"fitrange"` literal appeared 29 times
         across `RooFitUtils.cpp`, `InteractiveRooFitEditor.cpp` and
         `InteractiveSimultaneousFitEditor.cpp` — including inside this loader
         and in the `RooFit::Range(...)` arguments that drive the actual fits.
         The constant is now public (`:286`) and all 29 literals go through it.
<!---->
### EDM validity gate (interactive + non-interactive)
`3523967`, `26ff1b4`
<!---->
- [x] `src/RooFitUtils.cpp:2997-3022` (refit-after-load path) — validity is
      gated on EDM (`edm >= 0 && edm < 10 && covQual >= 2`) **or** `status()==0`,
      because Minuit2 returns a nonzero status on covariance grounds for fits
      that genuinely converged.
      Check: the threshold `10.0` is an absolute
      distance-to-minimum in NLL units — confirm it is appropriate for the
      extended-NLL magnitude here (~1e7) and that a fit with `status()!=0`
      and `edm >= 10` is still rejected.
- [x] `src/RooFitUtils.cpp:3103-3115` (non-interactive path, `3523967`) — the
      same EDM criterion is now applied to the cold-start path so a good fit is
      not discarded just because it landed on a covariance boundary.
      Check:
      this path and the refit path use the **same** `edm`/`covQual` test; a
      difference here is the exact bug `3523967` was fixing, so confirm parity.
- [x] `src/RooFitUtils.cpp:2008` — comment warns Minuit can return `covQual 0`
      with `edm` exactly `0`, which would discard every fit.
      Check: the gate
      treats `edm < 10` (not `edm == 0`) as converged, so an `edm==0` success is
      not rejected.
NOTE: This one could use some more thought.
<!---->
### Saved simultaneous params seed a real fit (not adopted as the result)
`48bb04f`
<!---->
- [x] `src/RooFitUtils.cpp:2966-3022` — with a saved `.simroofits` present the
      params now **seed** a minimisation instead of being adopted as the result
      with no minimisation.
      Check: after the load the code actually runs a
      minimiser and reports the *minimised* chi2/ndf, not the loaded values'
      chi2; confirm `sim_valid` is set from the fit, not from the load.
- [x] `include/RooFitUtils.hpp:331-343` — the long comment documents the old
      "load-and-adopt" trap and the new seed behaviour (cites 73mGe precal:
      ~1.4 from a saved seed vs ~8 cold, edm 24000x tolerance).
      Check: the
      behaviour matches the comment — a re-run is reproducible and responds to
      upstream changes, and the saved file is deleted-only-to-redo as the README
      states.
<!---->
### Peak-separation constraint (external value)
`8bbdb1b`
<!---->
- [x] `src/RooFitUtils.cpp:1950-2010` (`ConstrainPeakSeparation`) — penalises
      `(mu_hi - mu_lo)` toward `delta` via a RooFit external constraint rather
      than fixing either centroid.
      Check: the penalty is a Gaussian on the
      **difference** with width `sigma`, so the pair can slide together (the
      absolute scale stays free) but is forbidden from stretching; confirm it
      guards `sigma > 0` (`:1954`) and the identically-zero-spacing no-op
      (`:1996`).
- [x] `include/RooFitUtils.hpp:131-141` — the comment distinguishes locking the
      two mus (pins absolute scale + stretch) from constraining separation
      (measurement, not equality).
      Check: `ConstrainPeakSeparation` adds an
      external RooFit constraint and does **not** `setConstant` either mu.
<!---->
### Tail decay-length cap
`88cec44`
<!---->
- [x] `src/RooFitUtils.cpp:333,344,2238,2246` — `tail_ratio_max_` is passed into
      the PDF constructors to bound the exponential-tail decay length.
      Check:
      the cap bounds the **ratio** `tau/sigma` (default 100.0 preserves historic
      behaviour per `include/RooFitUtils.hpp:169-170`), not an absolute keV
      length, so a wide peak is not over-constrained.
- [x] `src/RooFitPhotopeakPdfs.cpp:223` — `Double_t tau = ratio * sigma;`
      reconstructs the decay length from the ratio.
      Check: the ratio is
      floored at `>= 1` (enforced in `0da3794`) so `tau >= sigma` always holds
      after the cap is applied.
- [x] `include/RooFitUtils.hpp:167-170` — comment: "Upper bound on the
      exponential tail DECAY LENGTHS, in keV." Check: the comment says *keV*
      but the mechanism is a dimensionless ratio cap — confirm which is right
      and fix the comment if it is misleading (nit, but it describes a physics
      knob).
<!---->
### Invalid-NLL and NaN fixes (simultaneous fits)
`7d4360a`, `222f155`
<!---->
- [x] `src/RooFitPhotopeakPdfs.cpp:69` — a "tiny positive floor keeps every
      component's normalized integral well-defined." Check: with a zero-integral
      component (e.g. a tail that was effectively zeroed) the floor prevents the
      X/X `0/0` that made `RooRatio::rangeProj` produce NaN; confirm the floor
      is small enough not to bias the normalisation at normal statistics.
- [x] `src/RooFitUtils.cpp:3074-3092` — `SetFitDebug` un-suppresses RooFit eval
      errors and prints the seed NLL so an invalid-NLL failure names the
      offending PDF.
      Check: debug is OFF by default (`fit_debug_`), and the
      seed-NLL print only fires when debug is on.
- [x] `src/RooFitUtils.cpp:3092` — comment: high statistics populate
      out-of-range bins and invalidate the NLL once.
      Check: the fix (clipping /
      `enableSilentClipping`) means values outside the current range are
      clipped, not thrown, so a high-stat simultaneous fit does not fail at
      high statistics.
NOTE: This is more of a question of the RooFit function enableSilentClipping than it is this code.
<!---->
### Tail decay parametrised as ratio over sigma
`0da3794`
<!---->
- [ ] `src/RooFitPhotopeakPdfs.cpp:223` and the `ratio` initialisation in
      `BuildChannelModel` — the tail decay is `tau = ratio * sigma` with
      `ratio >= 1`.
      Check: the lower bound on `ratio` is `1` (i.e. `tau >=
      sigma`) and the fit parameter is the dimensionless ratio, so the result
      is scale-invariant and the cap in `88cec44` acts on the same quantity.
- [ ] `include/FittingUtils.hpp:36-44` — the `PeakFitResult` tail fields are
      `ratio`/`amplitude` (1M-context commit refactored both backends).
      Check:
      the TF1 (`FittingUtils`) and RooFit (`RooFitUtils`) backends report the
      **same** ratio quantity so a constrained fit's `component_amplitude /
      gaus_amplitude` is identical across backends (README states this).
<!---->
### Unbinned-likelihood rewrite + analytical integrals
`ebaf64c`, `df06f5a`
<!---->
- [ ] `src/RooFitUtils.cpp:78-160` (`LoadEventsFromTree`, `RooDataSet` build) —
      the fit is unbinned extended maximum likelihood over an event-level
      `std::vector<Double_t>`; binning is display-only. Check: the fit never
      bins, and `LoadEventsFromTree` extracts the branch into the event vector
      without loss.
- [ ] `src/RooFitPhotopeakPdfs.cpp:182-360` — `getAnalyticalIntegral` /
      `analyticalIntegral` for `RooStepShelf`, `RooLowExpTail`, `RooLowLinTail`.
      Check: each custom PDF provides an analytic integral so the per-Minuit
      numerical re-normalisation is cheap, and the integrals are correct
      (spot-check one against a numeric integral at a test parameter set).
- [ ] `src/RooFitPhotopeakPdfs.cpp:49-83` — the `ExpErfc` derivation (comment at
      `:49-50`: for the tail, `a` and `b` share the sign of `y` and
      `b = a*tau/(sqrt2*sigma)`, so `exp(a)` overflow is cancelled by `erfc(b)`
      underflow). Check: the CPU `doEval` and the CUDA kernel use the **same**
      `ExpErfc` form so CPU and GPU agree to machine precision.
<!---->
### Amplitude params as ratios (TF1 backend)
`9649e02`
<!---->
- [ ] `src/FittingUtils.cpp:77-91` — `gaus_amplitude` is a free parameter and
      `exp_tail_ratio` scales the tail relative to it. Check: the tail
      amplitude in the returned `PeakFitResult` is consistent with the
      ratio-times-gaussian convention, so a constrained fit reusing a prior
      ratio is meaningful.
- [ ] `src/FittingUtils.cpp:1222-1240` (`FitDoublePeak`) — multi-peak variant.
      Check: the double/triple peak functions sum the per-peak components plus
      one shared background, and each peak's amplitudes use the same ratio
      convention as the single-peak path.
<!---->
### Hyper-EMG peak shape (added, then apparently superseded)
`ea8a1da`, `62217b4`
<!---->
- [x] `src/FittingUtils.cpp` — `ea8a1da` added `FittingFunctions::HyperEMGPeak`,
      `DoublePeakHyperEMG`, `TriplePeakHyperEMG`, `EMGLowComponent`,
      `EMGHighComponent`.
      **These symbols are NOT present at HEAD** — a later
      commit (likely `9649e02`'s large `FittingUtils.cpp` rewrite) removed or
      renamed them. Check: confirm whether the hyper-EMG shape is still
      available (search `grep -rn "EMG" src include` returns nothing), and if it
      was intentionally dropped, confirm the README (which `62217b4` documented)
      no longer promises it.
<!---->
## Signal processing — review these lines
<!---->
### Bulk-read fast path for CoMPASS conversion
`92df438`
<!---->
- [x] `src/InitUtils.cpp:362-557` (`ConvertCoMPASSBinToROOT`) — the fixed-width
      (no-waveform) records are parsed from a bulk buffer with `std::memcpy`
      field extraction (`:509-555`) and `tree->Fill()` (`:555`), instead of the
      old per-record `CoMPASSReader` loop. Check: the `off_*` field offsets and
      `memcpy` widths match the on-disk record layout, and the variable-length
      (waveform) `else` branch still uses the `CoMPASSReader`-based loop.
- [x] `src/InitUtils.cpp:509-556` — each memcpy-extracted record fills the tree
      once. Check: the fast path fills every field the old path filled (energy
      variants, flags, board, channel, timestamp) so the output tree is
      identical for the same input; the parallel in-memory path
      `ConvertCoMPASSBinToHits` (`:641`) does the same.
<!---->
## Tooling / IO / GUI
<!---->
### GPU-accelerated photopeak fitting
`52a4f4a`
<!---->
NOTE: these are mechanical, and I have confirmed that it works, so they're not going to get as much attention.
- [x] `gpu/Roo*.cu` (all four) + `src/RooFitPhotopeakPdfs.cpp` `doEval` — each
      custom PDF ships a `__global__` kernel; `doEval` uses
      `cudaPointerGetAttributes` to detect a device buffer and launches the
      kernel, else runs the OpenMP host loop. Check: the kernel and host loop
      compute the same value (see the ExpErfc parity item above), and the CUDA
      build is fully opt-in (`-DAU_USE_CUDA=ON`; the CPU build never compiles
      `gpu/*.cu`).
- [x] `include/RooFitPhotopeakPdfs.hpp:24-123` — `canComputeBatchWithCuda()`
      returns `true` only under `AU_ROOFIT_BACKEND_CUDA`. Check: on a CPU build
      this returns `false` for every PDF, so RooFit never routes device buffers
      to a kernel that was not compiled.
- [x] `flake.nix:29,72,139` + `CMakeLists.txt` — `rootWithCuda` override,
      `packages.cuda`, `-DCMAKE_CUDA_ARCHITECTURES` from `cudaCapabilities`
      (default `8.9`). Check: the CUDA capability matches the target GPU and the
      README's warning about not overriding the `nixpkgs` input (cache
      invalidation) is honoured.
<!---->
### X11 teardown guard + per-peak locking
`08c983e`, `abba503`
<!---->
- [x] `include/InteractiveEditorX11Guard.hpp:40-60` — installs a no-op X11 error
      handler for the editor's lifetime, resolved via `dlsym` so there is no
      link-time libX11 dependency. Check: the previous handler is saved and
      restored (`AURestoreXErrorHandler`), and a failed `dlsym` makes install a
      harmless no-op.
- [x] `src/InteractiveRooFitEditor.cpp:1027`, `src/InteractiveFitEditor.cpp:1048`
      — the guard is engaged for the editor loop. Check: the handler is restored
      on every exit path (normal close, error, exception) so a recoverable
      `BadWindow`/`BadDrawable` during redraw does not SIGSEGV the whole session.
- [x] `src/InteractiveFitEditor.cpp` / `RooFitEditor` / `SimultaneousFitEditor`
      destructors (`abba503`) — the segfault-on-teardown fix. Check: each editor
      destroys its own drawing primitives and restores global state (batch
      mode, directory) before ROOT objects are torn down, avoiding the double
      free / use-after-free.
<!---->
### Interactive simultaneous-fit editor
`bf45a6e`
<!---->
- [x] `src/InteractiveSimultaneousFitEditor.cpp:701-808` (`OnSliderMoved`,
      `OnEntryChanged`, `OnBoundsChanged`, `OnFixToggled`, `OnRangeChanged`) —
      the multi-channel GUI callbacks. Check: a parameter change in one channel
      that is linked to another updates the linked value, and a range change
      rebuilds the display histogram from the events (the fit stays unbinned).
<!---->
### Simultaneous-fit machinery
`231f577`, `703179e`
<!---->
NOTE: these two seem redundant
- [x] `src/RooFitUtils.cpp:2075-2163` (`AddChannel`, `SeedChannel`) — channel
      registration and seeding. Check: per-channel locks (`mu_fixed`,
      `bkg_yield_fixed`, `lock_shape_after_seed`, `shape_lock_per_peak`,
      `use_step_per_peak`) are stored in the channel config and applied at fit
      time, and a seed result is only used as starting values.
- [x] `src/RooFitUtils.cpp:2449-2480` (`SeedChannel` application) — the seed
      result's parameters are copied into the channel model as initial values.
      Check: seeding does not fix the parameters (they still float unless a
      separate lock says so) and a missing seed is a clean no-op.
NOTE: this is not redundant
- [x] `src/RooFitUtils.cpp` (shared union range, `703179e`) — the simultaneous
      fit uses a shared union range across channels and tighter sigma bounds.
      Check: `LinkParameter` / `LinkPeakShape` tie the shared degree of freedom
      correctly, and an unresolvable link fails cleanly (see the `26ff1b4`
      unresolvable-link item) rather than segfaulting.
NOTE: this behavior is purely to establish a consistent convention
<!---->
## Checked and excluded (not core analysis logic)
<!---->
- `db3c39c` — adds the whole `RooFitUtils` backend + `InteractiveRooFitEditor`
  and deletes the old Python `fit_models.py`/`simultaneous_fit.py`.
  The new
  C++ backend's *math* is covered by the items above; the Python file removal
  itself is bookkeeping.
- `71bd648` — 17-digit save precision + persisting Hesse errors + silencing
  eval-error spam. Covered by the save/load items; the precision constant is a
  format choice, not physics.
- `92df438` flake/`InitUtils.hpp` signature change — IO plumbing, not the
  analysis result.
- `62217b4` — README-only documentation.
- `20e8672` — the earlier interactive-editor segfault/use-flag fix; the
  remaining behaviour is exercised by the GUI teardown items above.
<!---->
Nits:
- `include/RooFitUtils.hpp:169-170` — the `SetTailRatioMax` doc says "DECAY
  LENGTHS, in keV" but the mechanism caps a dimensionless `tau/sigma` ratio.
  Reconcile the wording with the implementation before the reviewer is misled.
- `src/RooFitUtils.cpp:2650-2651` and the two `extract` helpers were refactored
  out of this session's style-audit fix (now `RooFitExtractParamDiagnostic`);
  confirm the refactor preserved the exact near-limit margin logic
  (`abs_tol 1e-4`, `frac_tol 1e-3`).
