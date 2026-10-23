# Linked doublet tails and linked-parameter reporting — review checklist, 2026-10-22
<!---->
Working tree on `main` at `7dc3ae8`.
Not committed.
`FitDoublePeak` gains `link_tails` (default `kFALSE`): the second peak uses the first
peak's `LowExpTailAmplitude`, `LowExpTailRatio`, `HighExpTailAmplitude` and
`HighExpTailRatio`, so the two lines carry one physical tail on each side. Linking now
works the way `LinkParameter` does in a simultaneous fit: the second peak's
`RooFitPeakModel` points at the first peak's `RooRealVar`, for `link_sigma` as well.
That makes a `link_sigma` or `link_tails` result report the shared value and error,
and the saved state write each shared parameter once. Option 4 of the request (fixing
the second Gaussian yield as a ratio of the first) is not part of this change; see
"Checked and excluded".
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Why
<!---->
In the 73mGe Ba-133 01/15 anchor (79.61 keV next to 80.998 keV, 1.4 keV apart at
`sigma` 0.82) two free tails are degenerate: refit from the production saved state,
the weak line's `R` rails at 2.997 (limit 3) and the 80.998 keV tail ends at
`k = 0.026`, `tau/sigma` 13.7 ± 15.3. A `link_sigma` result also reported
`sigma = 0.51` for the second line, `Sigma2`'s unused starting value, and saved it as
`Sigma2 0.51 0 1`.
<!---->
## Core analysis logic — review these lines
<!---->
### Linking
<!---->
- [x] `src/RooFitUtils.cpp:518-559` — `LinkSecondPeakToFirst`: records the second
      peak's own names (`Sigma2`; the four `*ExpTail*2`), points its struct at the
      first peak's variables, and rebuilds all five PDFs and both exponential-tail
      yields on the variables the struct now holds. The PDF and yield names
      (`*_pdf2_linked`, suffix `2_linked`) and the objects built for `link_sigma`
      alone are the ones the old block built, so that graph is unchanged.
      Check: the window factors (`BuildExpTailYields`, `:409-427`) take the second
      peak's own `mu` with the shared `sigma` and tail ratio, which is what the
      request asks for.
- [x] `src/RooFitUtils.cpp:1614-1615` — `FitDoublePeak` calls it after both peaks are
      built and before `BuildTotalModel`, so the `RooAddPdf` only sees the rebuilt
      PDFs and yields. The old `link_sigma` block also set `Sigma2` to `Sigma1`'s start
      and made it constant; `Sigma2` is now in no PDF of the model and in no parameter
      list, so neither is needed.
      Check: no path still reaches `Sigma2` or the second peak's own tail variables
      (`CollectAllParams` `:655-672`, the editor's rows, `ExtractPeakResult` `:864-890`
      all go through the struct).
- [x] `src/RooFitUtils.cpp:1693-1697`, `:1698-1705`, `:1724-1726`, `:1752-1757` — the
      automatic search with `link_tails`: the inter-peak group no longer releases,
      prunes or fixes peak 1's high tail or peak 2's low exponential tail, because
      they are the shared tails already decided by `TestLowSideGroup(0)` and
      `TestHighTailIndependent(1)`; releasing them would restart them at 0.15 / 1.5.
      With `link_tails` off both flags equal the old `use_*` flags.
      Check: you agree each shared tail should be tested once, and that an inter-peak
      group left with nothing to release (73mGe's configuration: low exponential tail
      only) still runs its one refit, as it already does when no component is
      offered.
<!---->
### Reported results
<!---->
- [x] `src/RooFitUtils.cpp:864-890` (unchanged) — `ExtractPeakResult` reads the struct,
      so peak 2's `sigma`, `sigma_error`, tail ratios and their errors are now the
      shared ones, and its tail amplitudes are the shared `R` times its own
      `gaus_amplitude`.
      Check: this is the reporting you asked for in point 2. A downstream that took a
      `link_sigma` result's `peaks[1].sigma` (0.51) now gets the fitted width; 73mGe's
      only `link_sigma` call (`CalibrateAndFit.cpp:1004`) uses `mu` and `mu_error`.
- [x] `src/RooFitUtils.cpp:2759-2783` (unchanged) — `ExtractPeakResultFor`, the
      simultaneous path: `ResolveOrCreate` (`:2190-2210`) already returns the source
      variable for a link target, so no placeholder exists there and a linked peak
      already reports the source's value and error (check F below).
      Check: you are content that point 2 needs no change in the simultaneous path.
      What does change there is a `SeedChannel` from a `link_sigma` `FitDoublePeak`
      result (`ApplySeedToChannel` `:2600-2636`), which now seeds `Sigma2` with the
      fitted width instead of 0.51.
<!---->
### Saved states
<!---->
- [x] `src/RooFitUtils.cpp:674-683`, `:686` — `DistinctParams` keeps the first
      occurrence of each variable; `CollectFloatingParams` uses it, so the
      post-editor reduced chi-squared (`ComputeReducedChi2(nullptr, …)`) counts the
      shared width once, as it did when `Sigma2` was a constant placeholder.
- [x] `src/RooFitUtils.cpp:974-975` — `SaveInteractiveParams` writes each variable
      once, under the first peak's name: no `Sigma2` line for `link_sigma`, no
      `*ExpTail*2` lines for `link_tails`.
      Check: you accept that a `link_tails` `.roofits` file holds the shared tail only
      as `LowExpTailAmplitude1` / `LowExpTailRatio1`.
- [x] `src/RooFitUtils.cpp:1015`, `:1023-1029`, `:1069-1072` — `LoadInteractiveParams`
      skips the names `LinkSecondPeakToFirst` recorded, before the name lookup, and
      counts them. An old `link_sigma` file's `Sigma2 0.51 0 1` is skipped; so are
      the second peak's tail lines when a separate-tail state is loaded into a
      `link_tails` model. A saved fixed flag is never applied by this loader in any
      case (`:1042-1043`, unchanged).
      Check: skipping is the behaviour you want rather than "treat as linked" (copying
      the line's value to the first peak's variable). The old `Sigma2` line holds a
      starting value, not a fit result, so copying it would be wrong.
- [x] `src/RooFitUtils.cpp:442-469` (unchanged) — legacy conversion with `link_tails`:
      the shared ratio is converted once, with peak 1's window factor, because peak 1
      is the first peak that holds it.
      Check: see "Migration" below. Converting the weak line's old tail gives a shared
      start the refit cannot leave.
<!---->
## Tooling
<!---->
- [x] `include/RooFitUtils.hpp:809-856` — `FitDoublePeak` documentation: `link_sigma`
      reworded, `link_tails`, what a linked parameter reports, the component search,
      saved states.
- [x] `include/RooFitUtils.hpp:192-198` — `RooFitPeakModel`: two peaks can hold one
      parameter pointer; walkers must de-duplicate.
- [x] `include/RooFitUtils.hpp:41`, `:409-411`, `:474`, `:482-486` — `<string>`, the
      name set, the helper declarations.
- [x] `src/InteractiveRooFitEditor.cpp:592-604` — the editor's chi2/ndf counts a
      parameter shared by two rows once, so the live label for a `link_sigma` doublet
      is unchanged now that row 11 is `Sigma1` instead of a fixed `Sigma2`.
- [x] `src/InteractiveRooFitEditor.cpp:667-678`, `:693`, `:727`, `:760`, `:779`,
      `include/InteractiveRooFitEditor.hpp:132-134` — `SyncAliasedWidgets`: after a
      slider, value, bound or fix change, the other row bound to the same variable is
      re-synced (bounds copied, then `SyncWidget`), so the two tabs never disagree.
      Check: try it on a `link_sigma` and a `link_tails` doublet; it is GUI code and
      was compiled but not driven.
- [x] `doc/pages/fitting.md:172-213` — new section "Linked doublets", including how to
      seed a `link_tails` fit from a separate-tail state.
<!---->
## Test evidence
<!---->
The check programs are scratch code outside the repository; they reach the private
parameter list to print every `RooRealVar` at 17 digits.
<!---->
- [x] `nix build` of `.#default`, `.#cuda` and `.#docs` succeeds; the docs warning log
      is empty; the compiler warnings are the pre-existing ones (two `-Wformat` in
      `RooFitUtils.cpp`, two in `BinaryUtils.cpp`).
- [x] Default calls unchanged: one program, compiled against `7dc3ae8` and against
      this tree, CPU and CUDA builds, runs `FitSinglePeak`; `FitDoublePeak` plain,
      with `link_sigma`, and constrained; constrained `FitTriplePeak`;
      `FitSimultaneous` with `LinkParameter` and `LinkPeakShape` links; a `link_sigma`
      doublet loaded from a pre-convention state with the `Sigma2 0.51 0 1` line; and
      73mGe's production Ba-133 01/15 fit from its saved state. Every `FitResult` field
      and every parameter's value, error and constness are identical to 17 digits in
      both backends, except peak 2's `sigma` and `sigma_error` in the three
      `link_sigma` fits, which now equal peak 1's (0.8244585 ± 0.0037586 instead of
      0.51 ± 0 for the production Ba fit on CPU).
      Check: identical floating-point results are the "bit-for-bit" you asked for.
      The loader's summary line changes from "fixed by model configuration" to
      "linked to the first peak's" for the `Sigma2` line.
- [x] Synthetic doublet, 200 toys (CUDA build): Ba-like lines at 79.6142 and
      80.9979 keV, `sigma` 0.82, Gaussian counts 13363 and 166000, one tail with
      `R = 0.589`, `tau/sigma = 10.5` (`k = 0.0700`) on both, flat background, window
      75-90 keV. Each toy is fitted from the same neutral state (`R` 0.3, `tau/sigma`
      5) through the saved-state path, so no component pruning, with the same
      Strategy 1 pass and library `RunFit` after it for both modes (next item).
      Shared tail: 193 of 200 status 0; `R` 0.5923 ± 0.0015, `tau/sigma`
      10.62 ± 0.06, `k` 0.0699; `mu1` bias +6.8 ± 7.7 eV, `mu2` +2.0 ± 0.9 eV, pull
      widths 1.09 and 1.07.
      Separate tails: 50 of 200 status 0 (the other 150 status 1); the weak line's
      `R` at its limit in 36; `k` of the 80.998 keV tail 0.046; `mu1` bias
      +52.8 ± 8.4 eV. `mu2` is not biased on average here (−0.3 ± 1.1 eV), but its
      spread is 16.1 eV against 12.6 eV, with pull width 1.64.
      Check: the residual +2.0 ± 0.9 eV on `mu2` with the shared tail (0.18 of the
      per-fit error) is acceptable or worth chasing.
- [x] Found while testing, not changed: the library's single-channel refit
      (`RunFit`, Strategy 2) can stall at a poor start. MIGRAD's seed Hessian is not
      positive definite there, the diagonal shift MnPosDef adds produces a step no
      line search recovers from ("No improvement in line search" twice, status 3),
      and every parameter is returned at its starting value. The saved-state path
      (`FitDoublePeak` `:1635-1642`; `:1494`, `:1866`, `:2015` in the other `Fit*`
      methods) does not check the status, so the start is returned as a valid
      result. It happened in all 200 shared-tail toys from the neutral start, in 0 of
      200 separate-tail ones, in the baseline build too (the pre-convention
      `link_sigma` state above, identical in both), and on real data from the legacy
      state below. A Strategy 1 pass before `RunFit` brought 193 of the 200 toys and
      the real-data case to status 0; that pass is in the toy harness, not in the
      library.
      Check: whether `RunFit` should retry (Strategy 1 then 2) or at least report a
      non-zero status. That would change today's results for fits that stall, so it
      is a separate decision from this change.
- [x] 73mGe Ba-133 01/15 (`filtered/PostReactor_Ba133_20260115`, 711821 events,
      75-90 keV, flat background, low exponential tail, `tau/sigma` ≤ 20,
      `link_sigma`), CUDA build:
      production state, separate tails: `R` 2.997 (railed) and 0.286, `k` of the
      80.998 keV tail 0.026, `mu2` 81.264336, chi2/ndf 1.558.
      `link_tails` from a neutral state: `R` 0.5983, `tau/sigma` 11.02 ± 0.93,
      `k` 0.0677, `mu2` 81.274704 ± 0.0083, chi2/ndf 1.505, status 0, `cov_qual` 3.
      `link_tails` from the production state with the shared tail started at the
      80.998 keV line's old tail: `mu2` 81.274780, `k` 0.0678.
      CPU build, from the unedited production state after a Strategy 1 pass:
      `mu2` 81.274713.
      So `mu2` moves +10.4 eV from production, and lies 2.9 eV below the shared-tail
      binned refit's 81.2776.
      Check: 2.9 eV is "a few eV", and `tau/sigma` 11.0 ± 0.9 against your 10.5 is
      fine.
- [x] Migration: the unedited production state loaded into a `link_tails` model
      stalls as described above. Its shared tail starts at the weak line's converted
      `R` = 0.935 (`k` 0.145), the starting chi2/ndf is 355, and the state comes back
      unchanged. Before 73mGe switches this fit to `link_tails`, the state needs a new
      start: delete it and refit in the editor, or copy the 80.998 keV line's
      `LowExpTail*2` values into the `*1` lines.
      Check: 73mGe is updated that way before it takes this build.
- [x] Automatic path (no saved state), for reference: it fails on this data with or
      without `link_tails`. Separate tails: the weak line's `R` rails at 3, the
      80.998 keV tail is pruned, `mu2` 81.329. Shared: the initial tail-less fit puts
      `mu1` at 76.94 keV, the shared tail is rejected, `mu2` 81.036. 73mGe uses the
      saved-state path.
- [x] Saved states: the `link_tails` endpoint saved through `SaveInteractiveParams`
      has no `Sigma2` and no `*ExpTail*2` lines, `LowExpTailAmplitude1` and
      `LowExpTailRatio1` carry the shared tail, and reloading it ignores nothing and
      returns the same point. Loading the old production state reports 1 line ignored
      as linked for `link_sigma` (`Sigma2`) and 5 for `link_sigma` + `link_tails`.
- [x] Reporting: `link_sigma` doublets report peak 2's `sigma` and `sigma_error`
      equal to peak 1's (synthetic: 0.8234853 ± 0.0027182 on both; production Ba as
      above). In the simultaneous fit with `c1:Sigma2` linked to `c1:Sigma1`, both
      peaks already reported 0.8058816 ± 0.0041394 in the baseline, and the shared
      tail ratio, 9.133 ± 0.475, on both.
<!---->
## Checked and excluded
<!---->
- Option 4, fixing the second Gaussian yield as a ratio of the first, is not
  implemented. It needs the tied yield to be a derived quantity, while
  `RooFitPeakModel::gaus_yield` is a free `RooRealVar` that the coefficient list, the
  tail yields, `ExtractPeakResult`, `AppendPeakGraphs`, the editor's component curves
  and the editor's "Yield" row all read; as an inert placeholder it would recreate the
  problem points 2 and 3 remove. A separate change.
- `CollectAllParams` keeps ten entries per peak, with the shared variable repeated,
  so `SetManualParameters` (`:250-272`), `SnapshotParams` / `RestoreParams`
  (`:746-769`) and the editor's fixed row layout keep their indices. Snapshot and
  restore write the same values to a repeated variable twice.
- A model without `link_sigma` loading an old `link_sigma` file still takes
  `Sigma2 = 0.51` as its starting value (floating, not fixed), as before.
- `ExtractParameterDiagnosticsSingle` (`:2882`) would list the shared variable under
  both peaks' names; it has no caller.
- The simultaneous editor (`InteractiveSimultaneousFitEditor`) shows linked variables
  in several rows and counts each row in its chi2/ndf, as before this change; the
  single-channel fixes above are not ported to it.
- `SortPeaksByMu` swaps whole structs, which keeps the links; `FitDoublePeak`'s
  centroid bounds stop it swapping in practice.
- FittingUtils (TF1 backend) has no `link_sigma` or `link_tails`.
<!---->
## Nits
<!---->
- The library version stays 26.10.22. 73mGe cannot pin a build with `link_tails` by
  version.
