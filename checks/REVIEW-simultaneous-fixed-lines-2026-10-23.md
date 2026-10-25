# Fixed lines in simultaneous fits — review checklist, 2026-10-23
<!---->
Branch `review/simultaneous-fixed-lines`, on top of `5d7032b`. **Merged into `main` on
2026-10-24 before ewt's review**, which is in progress; the boxes below are still open.
73mGe needs the K X-ray escapes of the 95.902 keV 113Cd(n,g) line in its in-situ
fit: eight lines at 64-73 keV whose positions and yields follow from that line, which
lies outside the 65-80 keV window. They sit under Pb Ka2 and beside and under Ge, so
they cannot be fitted, and a simultaneous channel had no way to carry a component of
fixed yield. `AddFixedLines(channel, shape_peak, mus, gaus_yields)` adds lines at fixed
centroids with fixed whole-line Gaussian yields, drawn with one of the channel's
peaks' shape. Each line is a `RooFitPeakModel` that owns only a constant centroid and
yield and holds the shape peak's own `Sigma`, step and tail variables, so it follows
that shape through the fit and adds no free parameter. The lines are kept apart from
the peaks: they enter the channel's `RooAddPdf` but not its results, saved states or
editor rows.
<!---->
Line numbers are from `main` at `69098de`, where every open checklist points since the
review branches went onto `main` on 2026-10-24. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
- [x] `src/RooFitUtils.cpp:2475-2512` — `AddFixedLines`: refuses an unknown channel,
      a shape peak outside the channel's peaks, mismatched lengths and negative or NaN
      yields; stores the lines in the channel's config, replacing any earlier call's.
      Check: replacing rather than appending is what you want for a second call.
- [x] `src/RooFitUtils.cpp:2781-2782` — `BuildChannelModel` adds the lines after every
      peak's components and before the background, so the existing coefficients keep
      their order.
- [x] `src/RooFitUtils.cpp:2810-2846` — `BuildChannelFixedLines`: each line starts as a
      copy of the shape peak's `RooFitPeakModel` (its resolved, possibly linked,
      `Sigma` and step and tail variables and the channel's window edges), then gets a
      constant `FixedLineMu` and `FixedLineYield`. `gaus_yield_scaled` is the yield
      times the Gaussian's share of the window, `erf`-difference over 2, as a
      `RooFormulaVar` of the yield, centroid, the shape peak's sigma and the window
      edges, so a line near an edge is counted as a fitted peak's `GausAmplitude` is.
      Check: the conversion from the whole-line yield the caller gives to the in-window
      coefficient.
- [x] `src/RooFitUtils.cpp:2847-2878` — the five component PDFs on the line's own
      centroid; the step and linear-tail yields are the in-window coefficient times the
      shape peak's ratios (the fitted peaks' convention), and `BuildExpTailYields`
      builds both exponential tails on it with the line's own window factors, so a
      tail's in-window yield is the whole-line yield times `R` times `f_T`.
- [x] `src/RooFitUtils.cpp:2880-2895` — the Gaussian always enters the sum; each
      other component only when the shape peak's ratio is not held at zero when the
      channel is built, which saves two to four PDFs per line in the usual
      configurations. Check: a ratio held at zero at build time cannot be released
      later in simultaneous mode (configuration-fixed parameters are skipped by the
      state loader, and the shape locks only fix).
      ewt: roofit code is so weird 
<!---->
## Tooling
<!---->
- [x] `include/RooFitUtils.hpp:323-328` — the channel config's three fields.
- [x] `include/RooFitUtils.hpp:463-467`, `:487-489`, `:585-587` — the per-channel line
      store, the plotted channel's lines, and the private declarations.
- [x] `include/RooFitUtils.hpp:1083-1122` — `AddFixedLines` documentation: what a line
      borrows from the shape peak, the whole-line yield against `GausAmplitude`'s
      in-window one, switched-off components, where the lines appear and where not.
- [x] `src/RooFitUtils.cpp:1424-1453` — `AppendFixedLinesGraph`: all of the plotted
      channel's lines as one dashed dark cyan curve on the background, every component
      on its in-window coefficient; nothing when there are none, so single-channel and
      line-free plots are unchanged.
- [x] `src/RooFitUtils.cpp:1490`, `:1540`, `:1592` — the three plot functions call it
      after the peaks' components.
- [x] `src/RooFitUtils.cpp:3225`, `:3242` — `PlotChannel` hands the channel's lines to
      the plot functions and clears them afterwards.
- [x] `doc/pages/fitting.md:297-299`, `:333-366` — the method in the simultaneous-fit
      list and a "Fixed lines" section with an example.
<!---->
## Test evidence
<!---->
- [x] `nix build` of `.#default`, `.#cuda` and `.#docs` succeeds; the docs warning log
      is empty; `clang-format` (LLVM, `BreakStringLiterals: false`) leaves both touched
      sources unchanged.
- [x] Testbed (CUDA, 40 spectra, `/tmp/73mGe-escape/FixedLinesTestbed.C`): a Pb
      Ka-like doublet at 72.80 / 74.97 keV (20000 / 34000 Gaussian events, sigma 0.80,
      one low tail with R 0.322 and tau/sigma 8, generated with the TF1 lineshape
      through `TailHeightRatio`'s relation) on a flat continuum, plus three lines of
      the same shape at 72.73, 69.50 and 64.90 keV (1600, 600, 400; the last half
      outside the 65-80 keV window). With the three as fixed lines at their true
      yields: mu1 bias −0.1 ± 2.1 eV, mu2 −0.3 ± 1.4 eV, mu1 pull width 0.95; sigma
      0.7992, R 0.3203 (0.3217), yields 19992 / 33996. Without them: +12.4 ± 2.0 and
      +7.6 ± 1.4 eV, R 0.368.
      Check: that recovering the yields and R with an edge line among the three is
      enough to settle the yield convention.
- [x] Default calls unchanged: 73mGe's production on this branch (CUDA) with its
      escapes switched off reproduces every row of its committed `tau20_sigpb` result,
      made on `5d7032b`, character for character (five Ge values, errors and chi2 to
      the file's six to seven digits).
- [x] 73mGe production with the escapes (CUDA, `tau20_sigpb_esc`): eight lines in
      each channel of the six in-situ fits, Gaussian yields 150-6000 events. Channel
      chi2/ndf as before or lower (Cu 01/13 signal 1.37 → 1.36, Cd 25% 1.51 → 1.47);
      Ge moves −11.7 eV (Cu, captures inside the crystal) and −0.7 / −0.8 eV (Cd
      shield pairs, external photons), against −12.7 / +0.1 / −1.2 eV from a TF1
      replica of the in-situ fit with the same lines; fit errors within 0.4 eV.
      Check: the plots (`log_PostCal_*`) show the lines' curve where expected.
<!---->
## Checked and excluded
<!---->
- `SaveSimInteractiveParams` and `LoadSimInteractiveParams` walk
  `sim_channel_peaks_` only; the lines have no free parameter and are not saved, so
  states saved without them load unchanged.
- `ComputeChannelChi2` evaluates the channel's `RooAddPdf`, so it includes the lines,
  and counts free parameters from it, which the lines' constants do not add to.
- `FitSimultaneous` extracts `cfg.num_peaks` peaks per channel; the lines are not in
  `FitResult::peaks`.
- The simultaneous editor draws its total curve from the channel PDF (lines
  included) and its component curves and rows from the peaks (lines not shown). GUI
  code, compiled but not driven.
- `GausYield()` of a line would return its whole-line yield; nothing calls it on a
  line (the plot reads `GausYieldCoefficient()`).
- CUDA: the in-window coefficient is a `RooFormulaVar` scalar, evaluated on the CPU
  like the exponential-tail yields.
<!---->
## Nits
<!---->
- `include/FittingUtils.hpp:249`, `src/RooFitUtils.cpp:2939` — `FitResult::lin_bkg_slope`
  defaults to −1 and `ApplySeedToChannel` applies any value above −1e6, so a
  hand-built seed that leaves it unset sets a held (flat-background) slope to −1,
  clamped at its lower bound. Found while writing the testbed; seeds that come from
  fits are unaffected. Not changed here.
