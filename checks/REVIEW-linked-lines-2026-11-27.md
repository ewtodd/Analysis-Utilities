# Linked lines in simultaneous fits — review checklist, 2026-11-27
<!---->
On `main` before ewt's review; review in progress, the boxes below are open.
73Ge measures the K X-ray escapes of the 95.902 keV 113Cd(n,g) line on the 01/15
active background run: eight lines at 64-73 keV whose positions and relative
intensities are known (Geant4's series ratios and shares) and whose common
intensity is the measurement. `AddFixedLines` can only carry lines of fixed yield,
and a peak's yield is its own, so that fit was a hand-written likelihood outside
AU, without the editor or saved states. `AddLinkedLines(channel, shape_peak, mus,
yield_ratios)` adds lines at fixed centroids whose whole-line Gaussian yields are
fixed multiples of one of the channel's peaks', drawn with its shape: fit one
escape line as the peak and link the other seven to it. The lines go through the
same build, model, plot and editor paths as `AddFixedLines` lines and add no free
parameter.
<!---->
Line numbers are from `main` at this commit. No code is reproduced here; open the
file at the line and review. Checks indicate review by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
- [ ] `src/RooFitUtils.cpp:2524-2563` — `AddLinkedLines` (declared at
      `include/RooFitUtils.hpp:1284`): refuses an unknown channel, a shape peak
      outside the channel's peaks, mismatched lengths and negative, NaN or
      infinite ratios; **appends** the lines to the channel's config, each with
      its own shape peak, so two peaks fitted apart can each carry a set.
      Check: appending (unlike `AddFixedLines`, which replaces) is what you want.
- [ ] `src/RooFitUtils.cpp:2875-2884` — the comment on what a fixed and a linked
      line's `gaus_yield` and coefficient are.
- [ ] `src/RooFitUtils.cpp:2885-2913` — `BuildChannelFixedLines` now loops over the
      fixed lines then the linked lines; each line takes its own shape peak (the
      channel's one `fixed_line_shape_peak` for fixed lines, its own entry for
      linked ones), checks it exists, copies its `RooFitPeakModel` and gets a
      constant centroid. Fixed lines keep their old names (`_fixed<i>`); linked
      ones are `_linked<i>`.
- [ ] `src/RooFitUtils.cpp:2915-2928` — the linked line's coefficient: a constant
      `LinkedLineRatio` and the `RooFormulaVar`
      `ratio * GausYieldCoefficient(peak) * f_line / f_peak`, each `f` the
      erf-difference window share at its own centroid with the shared sigma and
      the channel's window (the 1/2 cancels). Check: the peak's in-window yield is
      taken to the whole line through the peak's share and back into the window
      through the line's, so a line near an edge is counted as a fixed line is.
- [ ] `src/RooFitUtils.cpp:2929-2939` — the fixed-line branch, unchanged in
      substance (whole-line yield times its window share).
- [ ] `src/RooFitUtils.cpp:2940-2989` — the PDFs, step/tail yields and the
      switched-off-component rule, shared by both kinds as before; the inner
      component index is renamed `c` so it no longer shadows the line index.
- [ ] `src/RooFitUtils.cpp:2991-2998` — one line printed per kind present.
<!---->
## Tooling
<!---->
- [ ] `include/RooFitUtils.hpp:359-364` — the channel config's three per-line
      fields.
- [ ] `include/RooFitUtils.hpp:502-504` — the per-channel line store's comment.
- [ ] `include/RooFitUtils.hpp:1246-1286` — `AddLinkedLines` documentation: the
      whole-line ratio and its in-window form, appending, where the lines appear
      and where not, the editor's `GausAmplitude` moving them.
- [ ] `doc/pages/fitting.md:346-348`, `:421-451` — the method in the
      simultaneous-fit list and a "Linked lines" section with an example.
<!---->
## Test evidence
<!---->
- [ ] `nix build` of `.#default`, `.#cuda` and `.#docs` succeeds; the docs log's
      only warning is the dirty git tree; `clang-format` (LLVM,
      `BreakStringLiterals: false`) leaves both touched sources unchanged.
- [ ] `nix run .#cuda-regression`: all checks pass, the start NLLs equal to their
      references.
- [ ] Regression of the fixed-line path: UM-ANSG 73Ge production
      (`CalibrateAndFit.cpp`, 14 simultaneous fits with fixed escape lines) on
      this AU through `--override-input utils path:...` against the pinned
      26.11.26 run of the same day: the 599-line logs (every fitted parameter,
      chi2, capture line and Ge centroid) and the result file are identical.
- [ ] Testbed of the linked lines: see below.
<!---->
### Testbed
<!---->
- [ ] `LinkedLinesTestbed.C` (CPU build, in the session scratchpad, unbinned): eight
      escape lines of a 95.79 keV parent at `ReferenceLines::ESCAPE_XRAY_KEV` offsets
      with Geant4's shares, Cd K-beta 14000 events, Cd Ka/Kb 3.67, Te Ka/Cd Ka 0.053,
      Te Kb/Te Ka 0.30, sigma 0.837 keV and a low exponential tail (R 0.5,
      tau/sigma 9.6) generated from AU's own tail density, on a ~39000/bin linear
      continuum (~1.36 M events per spectrum). Fit: one peak (Cd K-beta1', centroid
      fixed), `HoldParameter` on its Sigma and tail, the other seven linked to it.
      Pulls of the whole-line peak yield against the true 11732:
      66-74.5 keV, 40 spectra: mean -0.04 +- 0.15, RMS 0.92, bias -0.11% +- 0.41%,
      covQual 3 in all; window 69.5-78.0 keV with the peak 0.25 sigma above the edge
      and Cd Kb2 and Te Ka straddling it, 80 spectra: mean -0.10 +- 0.12, RMS 1.05,
      bias -0.20% +- 0.23%. The linked lines are drawn in the fixed-line curve.
- [ ] Not this change, found by the testbed: a simultaneous channel bounds Sigma at
      0.1 x the window width (0.85 keV for 66-74.5 keV, against the measured 0.837),
      and `ComputeChannelChi2` counts a display bin the window edge cuts through.
