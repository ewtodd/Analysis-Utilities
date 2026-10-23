# Full-range exponential-tail amplitudes — review checklist, 2026-10-22
<!---->
Working tree on `main` at `da6489d`.
Not committed.
`LowExpTailAmplitude` and `HighExpTailAmplitude` change from in-window yield ratios `r`
to full-range yield ratios `R` (tail integral over the whole line divided by the
Gaussian's), so that the same value is the same physical tail on every peak. The PDFs
stay normalised over the fit range; each exponential tail's yield gains a window
coefficient `f_T / f_G` (`RooExpTailWindowFactor`):
`tail yield in window = GausAmplitude * R * f_T / f_G`. This is the default, with no
switch. `StepAmplitude` and `LowLinTailAmplitude` are unchanged: their shapes have no
finite full-range integral.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Why
<!---->
With `r` shared through `LinkParameter`, a peak whose tail is cut off by the window
(the Ge line at 68.75 keV in 65–80 keV) is forced to carry a taller tail than one that
is not (Pb Ka1 at 75 keV). For one physical tail (`k = 0.06`, `rho = 11`) the in-window
ratios are 0.18 (Ge), 0.30 (Ka2) and 0.35 (Ka1); sharing `r` gives Ge about twice its
true tail and pulls its centroid up (73mGe synthetic testbed: +32 eV). `R` is 0.529 for
all three.
<!---->
## Core analysis logic — review these lines
<!---->
### The window factor
<!---->
- [x] `src/RooFitPhotopeakPdfs.cpp:89-114` — `ExpTailWindowFraction`: `f_T` as the
      difference of `ExpTailAntideriv` (the antiderivative `RooLowExpTail` and
      `RooHighExpTail` already normalise with) over `2 tau exp(sigma^2 / (2 tau^2))`,
      that shape's whole-line integral. `GaussWindowFraction`: `f_G`, from `erfc`
      differences when both edges are on one side of the centroid.
      Check: the whole-line integral is `A(+inf) - A(-inf)` of that antiderivative.
- [x] `src/RooFitPhotopeakPdfs.cpp:483-529` — the class: proxies for `mu`, `sigma`,
      `tau_ratio` and the two window edges, and `evaluate()`. The high-side tail uses
      the window mirrored as `RooHighExpTail::analyticalIntegral` does; `f_G` is
      floored at `1e-300`; a non-positive result returns 0.
      Check: the mirroring (`y_lo = mu - x_hi`, `y_hi = mu - x_lo`) matches lines
      462-481 of the same file.
- [x] `include/RooFitPhotopeakPdfs.hpp:387-472` — the declaration and the formula as
      documented. No `doEval`: RooFit's default evaluates the scalar on the CPU in
      both backends, as it does the `RooFormulaVar` yields.
      Check: the documented `f_G` and `f_T` are the ones in the request.
<!---->
### Model construction
<!---->
- [x] `src/RooFitUtils.cpp:340-353` (`BuildPeak`) and `src/RooFitUtils.cpp:2346-2357`
      (`BuildChannelModel`) — the exponential-tail amplitude limit goes from 0.5 to
      `kExpTailAmplitudeMax` = 3 (`include/RooFitUtils.hpp:516-534`). Real data has
      `R` of 0.45–0.57 (and 0.9 implied by the old linked Ge tail), so 0.5 would rail;
      a long tail reaches 1.6. 3 covers every old `r` up to 0.5 while at least a sixth
      of the tail is inside the window (the Ge peak at `tau/sigma = 20`: `R` = 2.5).
      It is your floor rather than something larger because the wider the range, the
      deeper a typical 0.5 sits in the curved part of Minuit's bounded-parameter
      transform: in the linked closure (below) HESSE returned `cov_qual` 2 in 7 of 40
      fits with a limit of 5 and 3 of 40 with 3; centroids were identical.
      Check: 3 is acceptable.
- [x] `src/RooFitUtils.cpp:386-404` — `BuildPeak`: step and linear-tail yields stay
      `GausAmplitude * ratio`; the single-channel window variables are created once,
      at the fit range, and both exponential-tail yields come from
      `BuildExpTailYields`.
      Check: the fit range is what the single-channel yields count events in (`RunFit`
      fits `Range("fitrange")` and `x`'s default range is set to it after the data
      are built).
- [x] `src/RooFitUtils.cpp:409-427` — `BuildExpTailYields`: one low-side and one
      high-side factor per peak, and the yields `GausAmplitude * R * factor`.
      Check: the `sigma` passed in is the one the tail PDFs use.
- [x] `src/RooFitUtils.cpp:1575-1576` — `FitDoublePeak(link_sigma)`: peak 2's factors
      and yields are rebuilt on peak 1's sigma along with its PDFs; otherwise they
      would follow the frozen `Sigma2`.
- [x] `src/RooFitUtils.cpp:2302-2309` and `src/RooFitUtils.cpp:2421-2431` —
      `BuildChannelModel`: one window pair per channel, at the channel's own range,
      shared by its peaks; exponential-tail yields as above.
- [x] `src/RooFitUtils.cpp:3018-3031` — `FixChannelCoefRanges` moves each channel's
      window to the range its yields are fixed to: the channel's own range for
      `SplitRange`, the shared `fitrange` for the editor.
      Check: every path that changes what a channel's yields count also moves the
      window.
      The paths are this function, the two loaders below and the editors'
      range handlers. The replay path (no refit) never calls this function; there the
      window is the file's `RANGE`, which is also `x`'s default range.
<!---->
### Editors
<!---->
- [x] `src/InteractiveRooFitEditor.cpp:774`, `:791-795`, `:835`, `:857` and
      `src/InteractiveSimultaneousFitEditor.cpp:821`, `:842-850`, `:916`, `:941` —
      the range slider, Cancel and Reset move every peak's window with `fitrange`.
      Check: no other place in either editor sets `fitrange`.
<!---->
### Saved states
<!---->
- [x] `src/RooFitUtils.cpp:8-16` — the marker `ExpTailAmplitudeConvention 1 0 1`,
      written as a four-field line so that older builds and 73mGe's
      `InjectionTest::LoadSavedState` read it as one unknown parameter instead of
      losing their place.
- [x] `src/RooFitUtils.cpp:435-469` — `ConvertLegacyExpTailAmplitudes`:
      `R = r / factor` (that is, `r * f_G / f_T`) at the loaded `mu`, `sigma`, `rho`
      and window; the error is scaled by the same factor; a value above the limit is
      clipped with a warning. A ratio shared through a link is converted once, with
      the first peak that uses it, which is its source.
      Check: you accept that a linked old state cannot start exactly where it was
      saved. The source peak does; the peaks linked to it start with the same
      physical tail instead of the old shared in-window ratio. 73mGe refits all of
      these.
- [x] `src/RooFitUtils.cpp:918-920` and `src/RooFitUtils.cpp:1030-1032` — both savers
      write the marker after `RANGE`.
- [x] `src/RooFitUtils.cpp:950`, `:962-968`, `:985`, `:995-1003` —
      `LoadInteractiveParams`: the window follows the file's `RANGE`; the marker is
      recognised; files without it are converted after every value is loaded.
- [x] `src/RooFitUtils.cpp:1112-1115`, `:1121-1127`, `:1143`, `:1153-1166` —
      `LoadSimInteractiveParams`: the same for simultaneous files; every channel's
      window is set to `RANGE`, which `FitSimultaneous` then gives every channel.
      Check: converting with the window from the file's `RANGE` is right for your
      `.simroofits` files (all channels share one window there).
<!---->
### Reported results
<!---->
- [x] `src/RooFitUtils.cpp:823`, `:831`, `:2728`, `:2736` (unchanged lines) —
      `low_exp_tail_amplitude` and `high_exp_tail_amplitude` are still
      `ratio * GausAmplitude`, so their meaning changes to `GausAmplitude * R`.
      Seeding
      (`ApplySeedToChannel`, `src/RooFitUtils.cpp:2578`) and the constrained fits
      (`:1799`, `:1950`) divide by `GausAmplitude` again and now carry `R`, which
      transfers exactly between fits with different windows.
      A constrained
      `FitDoublePeak` / `FitTriplePeak` whose source fit used a different window
      therefore changes result.
- [x] `include/FittingUtils.hpp:144-196` — `PeakFitResult`: a table of what each
      amplitude holds per backend, `R`, the conversion to the FittingUtils height ratio
      `k = R sqrt(2 pi) / (2 rho exp(1 / (2 rho^2)))` and to the in-window tail yield.
      Check: this is the documentation you asked for, and the old claim that
      `amplitude / gaus_amplitude` means the same in both backends (false before this
      change too) is gone everywhere: also `include/RooFitUtils.hpp:344-351` and
      `doc/pages/fitting.md:96-101`.
- [x] `include/RooFitUtils.hpp:169-246` — `RooFitPeakModel`: which ratios are full
      range and which are in window, the window variables and factors, and
      `SetWindow()`.
- [x] Downstream, 73mGe (you are updating these): `GainTransferFit.hpp`
      `TailHeightRatio()` becomes the closed form above and no longer needs the
      window; `InjectionTest.cpp` `ExtractGeShape()` reads `R` from new state files
      (marker present) and `r` from old ones.
      Check: both are updated before 73mGe takes this build.
<!---->
## Tooling
<!---->
- [x] `include/RooFitUtils.hpp:388-392`, `:426-433` — private members and helpers.
- [x] `include/LinkDef.h:11` — dictionary entry for the new class.
- [x] `src/InteractiveRooFitEditor.cpp:997-1001`,
      `src/InteractiveSimultaneousFitEditor.cpp:1190-1194`,
      `include/InteractiveRooFitEditor.hpp:138-140`,
      `include/InteractiveSimultaneousFitEditor.hpp:177-180` — the editors' fallback
      bounds for parameters without limits (all have limits, so not reached) follow
      the new limit; helper declarations.
- [x] `include/RooFitUtils.hpp:889-896`, `:922-923` — `LinkParameter` and
     `LinkPeakShape` notes: linked exponential tails are the same physical tail;
     linked step and linear-tail ratios are not.
- [x] `tests/CudaRegression.cpp` — the toy models' exponential-tail yields are built
      as the library builds them (`ExpTailYield`, `:113-132`), windows follow the
      yield range (`:351-363`), and the second parameter point also moves `R`
      (`:434-435`).
      `R` starts at the value that reproduces the old start yields in
      each channel's own window, so the single and `SplitRange` reference NLLs are
      unchanged; the shared-range one is re-recorded (`:83-92`) because over the union
      range a fixed `R` puts more of each tail inside the window.
      Check: re-recording the shared-range reference is justified; the direct
      likelihood check passed for it first.
- [x] `doc/pages/fitting.md:114-170` — new section "Tail amplitudes": the
      convention, conversions, the limit, saved-state compatibility.
<!---->
## Test evidence
<!---->
- [x] `nix build` of `.#default`, `.#cuda` and `.#docs` succeeds; the docs warning
      log is empty; the compiler warnings are the pre-existing ones (two `-Wformat` in
      `RooFitUtils.cpp`, two in `BinaryUtils.cpp`).
- [x] Check 1, single-peak fits unchanged: the three 73mGe Am-241 59.5 keV runs
      (`calibration_low_tau20` states, 73mGe's `FitCalPeak` configuration), baseline
      `da6489d` against this tree.
      Automatic path: `mu`, `sigma` and yields identical
      to 0.01 events, NLL equal to 6e-9 absolute, `mu_error` equal to 0.1%.
      Production path (saved legacy state loaded, converted, refit): `mu` within
      0.02 eV, yields within 0.004 sigma, NLL within 1.5e-4 (the new fit ends lower
      each time; both stop within Minuit's EDM tolerance), and 01/16 equal to 9
      digits.
- [x] Check 1, centroid errors: on the production path `mu_error` agrees to 0.6%
      for 01/13 and 01/15, but 01/16 gives 11.15 eV against the baseline's 11.70 eV
      with `mu` equal to 8 digits. MINOS on 01/16 is identical in both
      parametrisations (−11.684 / +11.891 eV to 1e-4 eV), so the likelihood is
      unchanged; HESSE on that fit moves by ±4% with incidental settings alone
      (11.65–12.54 eV for `R` limits of 1, 1.5, 2, 3, 5 and 50), and the baseline's
      own automatic and loaded fits differ by 0.8–1.9%.
      Check: if the 73mGe calibration error budget is sensitive at the few-percent
      level, the Am/Ba centroid errors should come from MINOS rather than HESSE,
      independently of this change.
- [x] Check 2, `nix run .#cuda-regression`: all checks pass. CPU/CUDA NLL agreement
      ≤ 2.2e-16 relative over 40 points in all three configurations; the direct
      likelihood check passes in all three (the evaluators compute the factor as the
      scalar path does); single and `SplitRange` start NLLs match the old references
      to 1e-14; offset fits converge with CPU/CUDA parameters equal to 8e-4 sigma;
      CUDA at 0.04 of the CPU time per evaluation, as before.
- [x] Saved states, simultaneous (linked testbed after `RunBkgGeSim`, legacy
      `.simroofits` holding the old shared `r`): the new build converts `r = 0.29802`
      to `R = 0.52857` (truth 0.52878); on replay the source peak's in-window tail is
      7450.4986 events against the old 7450.5, and the Ge tail drops from 5960 to
      3496. Refit and save writes the marker; reloading that file does not convert
      and refits to the same point. The baseline build reading the new file warns
      about `ExpTailAmplitudeConvention` and parses everything else.
- [x] Saved states, single channel: the 01/16 Am state converted by hand to `R` with
      the marker added loads without conversion and gives a result identical to the
      legacy file's.
- [x] Linked closure, 40 testbeds at Cu 01/14 scale (bkg: Ka2 125k, Ka1 216k on
      1.34M flat; sig: 258k, 450k and Ge 53k on 3.5M), one physical tail
      (`k = 0.06`, `rho = 11`), `RunBkgGeSim`'s links and `FitSimultaneous` defaults,
      CUDA builds, the same testbeds for both (data checksums compared). Baseline:
      Ge `mu` bias +34.6 ± 3.0 eV, fitted `r` 0.323, `tau/sigma` 10.0. This tree: Ge
      `mu` bias +1.2 ± 3.1 eV, Ka1 −0.7 ± 0.6 eV, fitted `R` 0.531 (truth 0.529),
      `tau/sigma` 11.2 (truth 11); spread of the Ge bias 19.4 eV, median HESSE error
      20.9 eV (baseline 20.2 eV).
      Check: you still run SynthSim "physical" from 73mGe (check 3).
- [x] Linked closure, centroid errors: with this tree, HESSE inflates the Ge error in
      a few fits at the default strategy 1: 1 of 40 (75.6 eV, testbed 4, `cov_qual`
      3) with the limit of 3, 3 of 40 (42–58 eV, `cov_qual` 2) with a limit of 5;
      the baseline had none (max 22.5 eV). `cov_qual` is not 3 in 3 of 40 fits (0 for
      the baseline). The centroids are unaffected. Every such fit rerun with
      73mGe's `FIT_HARD` settings (strategy 2, final HESSE, 200k calls), seven in all,
      gave `cov_qual` 3 and 20–24 eV. The baseline is not immune either: with
      `FIT_HARD`, its testbed 8 failed (status 3, EDM 650). The coupling the new model
      adds (the Ge tail's in-window yield depends on the Ge `mu` through the window
      factor) is physical; the factor is analytic over the fit region (none of its
      branch switches can trigger with `mu` inside the window).
      Check: the Ge centroid errors you quote come from `FIT_HARD` fits, or are
      checked against the testbed spread.
<!---->
## Checked and excluded
<!---->
- Starting values: `ConfigureComponentFlagsForPeak` and `ReleaseComponent` still set
  0.1 / 0.15 and `tau/sigma = 1.5`.
  At `rho = 1.5` the factor is 1 to about 1e-14 for
  a peak inside its window, so the starting model is unchanged.
- `SetManualParameters` / `SetManualParameter` take values in the model's parameter
  order, so manual values for the exponential-tail amplitudes are now `R`.
- Plots, chi-squared, CSV dumps and the editors' component curves use the yield
  `RooFormulaVar`s, so they show in-window yields as before.
- The window variables are constant `RooRealVar`s outside every parameter list
  (`CollectAllParams`, the editors, the diagnostics, the save files), and
  `ComputeChannelChi2` counts only floating parameters.
- FittingUtils (TF1 backend) is untouched.
<!---->
## Nits
<!---->
- `LinkPeakShape`'s documentation says position stays free, but it links `Mu` too
  (`src/RooFitUtils.cpp`, `shape_params`). Not changed here.
- The library's version is still 26.10.22, the same as the build before this change,
  so the documentation describes the old convention as "until this convention was
  introduced" rather than by version.
  A version bump would let 73mGe pin it.
