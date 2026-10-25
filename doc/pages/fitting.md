Photopeak fitting {#fitting}
=================
<!---->
Fit gamma-ray spectral photopeaks with a composable model.
Two backends are
provided with identical APIs and result types so call sites can swap one for the
other.
<!---->
[TOC]
<!---->
## Backends
<!---->
- **`FittingUtils`** — `TF1` + ROOT::Math::Minimizer (Minuit2), binned
  chi-squared.
  Fast, pointwise un-normalized model evaluation.
- **`RooFitUtils`** — `RooAddPdf` of custom `RooAbsPdf` components, **unbinned
  extended maximum likelihood**.
  The custom PDFs (`RooStepShelf`,
  `RooLowExpTail`, `RooLowLinTail`, `RooHighExpTail`) implement `doEval()` so
  every NLL pass is vectorized over events rather than dispatching scalar
  `evaluate()` calls per event.
  The CPU `doEval` paths are OpenMP-parallelized;
  when built with CUDA support each PDF also includes a `__global__` kernel
  under `gpu/` and `doEval` dispatches to it when RooFit hands it
  device-resident buffers — see @ref cuda.
  Each component PDF still re-normalizes
  numerically over the fit range per Minuit step, but the per-step cost is
  amortized.
  Provides RooFit's full machinery for downstream extensions (e.g.
  simultaneous fits).
<!---->
Both classes take the same constructor signature, expose the same `Set*` flag
setters, and return the same `FitResult` struct. The interactive fit editor is
implemented per-backend (`InteractiveFitEditor` for TF1,
`InteractiveRooFitEditor` for RooFit, and `InteractiveSimultaneousFitEditor` for
multi-channel simultaneous RooFit fits). Saved-parameter files use the `.fits`
extension for the TF1 backend and `.roofits` for the RooFit backend, so the two
can coexist for the same input.
<!---->
## The model
<!---->
**Base model**: Gaussian peak + linear or flat background.
<!---->
**Optional components** (individually togglable, for semiconductor detectors
with segmented electrodes):
<!---->
- **Step function** — models events where part of the photon energy escapes the
  active volume, producing a step-like distribution on the left side of the
  peak, smeared by the detector resolution.
- **Low-energy exponential tail** — exponential tail below the photopeak from
  incomplete charge collection, charge trapping, etc.
- **Low-energy linear tail** — linear tail component that can capture asymmetric
  tailing not well described by a single exponential.
- **High-energy exponential tail** — exponential tail above the photopeak from
  pileup effects.
<!---->
Components are tested using a hybrid group-and-prune approach: low-side
components (step + both low-energy tails) are enabled as a group, then
individually pruned if they do not improve the fit. The high-energy tail is
tested independently. A component is kept only if it improves the reduced
chi-squared.
<!---->
**Multi-peak fitting**: double and triple peak variants with all components
enabled by default, and constrained fitting using results from previous fits.
<!---->
![Fit example](FitExample.png)
<!---->
## Results and diagnostics
<!---->
All fits produce structured results (`FitResult` containing `PeakFitResult`
entries) with parameter values, errors, and reduced chi-squared. Failed fits
return -1 for all parameters.
<!---->
RooFitUtils additionally populates Minuit2 / `RooFitResult` diagnostics on each
returned `FitResult` so downstream code can distinguish a true minimum from a
local or failed solution:
<!---->
- `fit_status` (Minuit2 status, 0 = success; from
  `RooFitUtils::kFirstInvalidMinimumStatus`, 3, up there is no valid minimum and
  the parameters may still be the starting values) and `cov_qual` (0 = Exact,
  1 = NotPositiveDefinite, 2 = Approximate, 3 = External)
- `edm` (estimated distance to minimum) and `min_nll` (minimized negative
  log-likelihood), with `has_fit_diagnostics` flagging whether they are
  meaningful. Every fit minimizes with `RooFit::Offset("initial")`: the NLL
  at the starting point is subtracted, so Minuit2 works with values near zero
  instead of ~1e7 and can resolve EDM below 1e-3. `min_nll` is therefore the
  improvement over the start, not the absolute NLL, and values from separately
  run fits are not comparable.
- `parameter_diagnostics` — one `FitParameterDiagnostic` per fitted
  `RooRealVar` carrying the value, error, limits, and near-limit flags (set when
  a parameter sits within a small margin of its bound)
<!---->
`FitSimultaneous` fills all of these, plus `parameter_correlations`. A
single-channel RooFitUtils fit fills them from the minimisation its result comes
from: the refit of a saved state, or the final fit of the automatic search (the
initial fit, when that fails and ends the search). A result accepted in the
interactive editor has none, and FittingUtils leaves every field at its sentinel
default.
<!---->
A single-channel minimisation that ends without a valid minimum is retried once:
a Strategy 1 pass, then the Strategy 2 fit again. Strategy 2 seeds MIGRAD with
the full Hessian at the starting point, and far from the minimum that matrix can
be indefinite; the fit then stops where it started, which a saved-state refit
used to return as its result without comment. A fit still without a valid
minimum after the retry prints a warning, and its `fit_status` says so. A
retried fit sets `fit_retried`: the Strategy 1 pass usually recovers the stall,
but from a poor start it can also reach a different local minimum with status
0, which only the reduced chi-squared shows.
<!---->
Reduced chi-squared is displayed on fit plots by default. Individual
peak components are plotted summed with the background for readability, and
multi-peak fits use distinct line styles per peak.
<!---->
## Backend-specific notes for RooFitUtils
<!---->
- Input is an event-level `std::vector<Double_t>` plus a display bin width (used
  only for plotting and the post-hoc reduced chi-squared used for component
  pruning). A static `LoadEventsFromTree` helper extracts events from a ROOT
  `TTree` branch. The fit itself never bins.
- Amplitude/yield semantics: gaussian, step and tail "amplitude" fields in
  `PeakFitResult` carry RooFit yields (event counts) rather than peak heights,
  and the ratio `component_amplitude / gaus_amplitude` is a yield ratio, not
  the height ratio FittingUtils reports. See
  @ref tail-amplitudes "Tail amplitudes" below for what each ratio means and how
  to convert.
- The low-energy linear tail factor `(1 + s·(x-μ))` is floored at zero to keep
  the PDF non-negative as required by RooFit normalization. The TF1 backend has
  no such constraint.
- `RooRealVar::enableSilentClipping()` is enabled so observable values outside
  the current range are silently clipped rather than throwing.
- **CPU vs CUDA backend**: by default `fitTo()` calls dispatch through
  `RooFit::EvalBackend::Cpu()`, and the OpenMP-parallelized `doEval` path
  handles batched evaluation. When compiled with `-DAU_ROOFIT_BACKEND_CUDA=1`,
  `BestAvailableBackend()` (in `RooFitUtils.hpp`) returns
  `RooFit::EvalBackend::Cuda()` instead, and all four custom PDFs report
  `canComputeBatchWithCuda() = true`. See @ref gpu for the nix-flake opt-in.
<!---->
### Tail amplitudes {#tail-amplitudes}
<!---->
Every component is normalised over the fit window, and its yield is the
Gaussian yield `GausAmplitude` (Gaussian events inside the window) times a
ratio parameter.
<!---->
The two exponential tails use the **full-range yield ratio** `R`
(`LowExpTailAmplitude`, `HighExpTailAmplitude`): the tail's integral over the
whole line divided by the Gaussian's. Their yield inside the window is
<!---->
    GausAmplitude * R * f_T / f_G
<!---->
where `f_G` and `f_T` are the fractions of the Gaussian and of the tail that
fall inside the window for that peak's `mu`, `sigma` and `tau = rho * sigma`.
`RooExpTailWindowFactor` computes `f_T / f_G`.
For a tail
`exp(y/tau) erfc(y/(sqrt2 sigma))` with height ratio `k` to its Gaussian,
<!---->
    R = k * 2 rho * exp(1 / (2 rho^2)) / sqrt(2 pi)
<!---->
which depends on the tail's shape only, not on the window or on where the peak
sits in it. So sharing `R` and `rho` through `LinkParameter()` gives every
linked peak the same physical tail, and a constrained fit or a `SeedChannel()`
seed carries the tail over exactly even when the earlier fit used a different
window. The bound on `R` is `RooFitUtils::kExpTailAmplitudeMax` (3); measured
values are typically around 0.5, and a long tail can exceed 1. A wider bound
makes HESSE less reliable in linked simultaneous fits; the constant's
documentation has the numbers.
<!---->
The step and the linear tail (`StepAmplitude`, `LowLinTailAmplitude`) are still
**in-window yield ratios**: their yield is `GausAmplitude * ratio`.
Both shapes
extend to `-infinity` without falling to zero, so they have no full-range
integral to form an `R` from.
Linking them between peaks at different positions
in the window, or between windows of different width, ties their in-window
yields rather than their shapes.
<!---->
In a `PeakFitResult` from RooFitUtils, `low_exp_tail_amplitude / gaus_amplitude`
is `R`.
To compare with FittingUtils, which reports heights, the height ratio is
`k = R * sqrt(2 pi) / (2 rho exp(1 / (2 rho^2)))`.
The tail's events inside the
window are `gaus_amplitude * R * f_T / f_G`, not `low_exp_tail_amplitude`.
<!---->
Until this convention was introduced, the exponential-tail amplitudes were
in-window ratios `r = R * f_T / f_G` as well. A single, unlinked peak fits to
the same model either way; only the reported ratio changes. Where `r` was shared
between peaks at different positions, it forced the peak with less of its tail
inside the window to carry a taller tail.
<!---->
**Saved states.** `.roofits` and `.simroofits` files now carry a line
`ExpTailAmplitudeConvention 1 0 1` after `RANGE`. A file without it is from
before the change; on load, each exponential-tail amplitude is converted with
`R = r * f_G / f_T` at the file's window and the loaded `mu`, `sigma` and
`rho`, so an unlinked peak starts exactly where it was saved.
A ratio linked
between peaks is converted with the peak it belongs to (the link's source), and
the peaks linked to it then start with the same physical tail, not the old
in-window ratio, so their start differs from the saved fit; a refit
(`SetRefitAfterLoad()` for simultaneous fits) recovers them. An older build of
this library reading a new file warns that `ExpTailAmplitudeConvention` is not
a model parameter and reads `R` as `r`.
<!---->
### Linked doublets {#linked-doublets}
<!---->
`FitDoublePeak(input_name, peak_name, mu1_init, mu2_init, link_sigma,
link_tails, gaus_yield_ratio)` can make the two peaks share parameters, the
single-channel counterpart of `LinkParameter()`:
<!---->
- **`link_sigma`** — the second peak uses `Sigma1` as its Gaussian width.
- **`link_tails`** — the second peak uses `LowExpTailAmplitude1`,
  `LowExpTailRatio1`, `HighExpTailAmplitude1` and `HighExpTailRatio1`. Because
  the amplitude is the full-range ratio `R`, that is one physical tail: each
  peak's tail is `R` times its own Gaussian yield, scaled into the window with
  its own `mu` and its width.
  Use it when the lines overlap enough that two free
  tails are degenerate, for example Ba-133 79.61 keV next to 80.998 keV
  (1.4 keV apart at `sigma` around 0.8 keV), where the weak line's tail
  otherwise runs to its limit and takes the region below both lines, or for two
  lines seen in one detector that should carry the same tail. The step and the
  linear tail are not shared: their amplitudes are in-window ratios, so sharing
  them would tie yields rather than shapes.
- **`gaus_yield_ratio`** — when positive, the first peak's Gaussian yield is
  fixed at that multiple of the second's, `GausAmplitude1 = gaus_yield_ratio *
  GausAmplitude2`, and only `GausAmplitude2` is fitted. The peaks are in
  centroid order, so it is the lower line's yield over the higher line's: for
  Ba-133, 79.61 / 80.998 keV = 2.65 / 32.9 = 0.0805. It is a ratio of Gaussian
  yields, equal to the intensity ratio times the efficiency ratio when both
  lines carry the same tail (`link_tails`) and the same window fraction. With a
  known branching ratio it removes the remaining trade-off between the weak
  line's amplitude and the shared tail.
  Start it from a state with a realistic
  tail, such as a saved `link_tails` fit: from a tail much too short, the weak
  line can settle in the tail region instead (a local minimum with reduced
  chi-squared around 10), because its yield can no longer adjust.
<!---->
All three default to off, and with all three off the fit is unchanged.
A linked
parameter is reported with the value and error of the parameter it follows:
with `link_sigma` both peaks' `sigma` and `sigma_error` are `Sigma1`'s, and with
`link_tails` both report the shared tail ratio, with
`low_exp_tail_amplitude = R * gaus_amplitude` for each peak's own
`gaus_amplitude`. With `gaus_yield_ratio` the first peak's `gaus_amplitude` and
its error are the ratio times the second's. In the editor, the first peak's
yield row then edits `GausAmplitude2` and is labelled with the ratio.
<!---->
In the automatic component search the shared low-side tail is tested in the
first peak's low-side group and the shared high-side tail as the second peak's
high tail; the inter-peak group then tests only the second peak's step and
linear tail. In the interactive editor both peaks' tabs show the shared
parameter, and editing either row changes it.
<!---->
A saved `.roofits` state writes a linked parameter once, under the name of the
peak it belongs to: the first peak's for a width or tail, `GausAmplitude2` for
a tied yield. On load, lines naming the replaced parameters are ignored.
States written before this behaviour carry a `Sigma2` line for a `link_sigma`
fit holding its unused starting value, marked as fixed; it is ignored too, and
the loader never takes a parameter's fixed flag from a file in any case.
A state
saved from a fit with separate tails starts the shared tail from the first
peak's old tail; when the first peak is the weaker line, that tail is usually
the poorly constrained one, so delete the state and refit, or copy the stronger
line's `LowExpTailAmplitude2` and `LowExpTailRatio2` values into the `1` lines.
<!---->
## Interactive fitting
<!---->
`SetInteractive()` enables a GUI editor that opens after the automated fit
completes.
The editor shows the histogram, total fit, and individual components
with a residual (pull) panel marked by dashed ±3σ guide lines, all updating in
real time as parameters are adjusted via sliders and number entries.
A fit range
slider allows adjusting the fit range visually; for the RooFit backend the
display histogram is rebuilt from the underlying events whenever the range
changes, so binning follows the zoom in real time (the fit itself remains
unbinned).
Parameters can be fixed/freed via checkboxes, and a Refit button runs
the underlying minimizer with the current values as initial guesses.
Live
chi-squared is displayed on the plot.
<!---->
![Interactive fit editor GUI](GUI.png)
<!---->
On accept, parameters and fit range are saved to a `.fits` (TF1) or `.roofits`
(RooFit) file in the `plots/fits/` directory.
On subsequent runs with
`SetInteractive()`, saved parameters are loaded automatically and the editor is
skipped, so interactive results only need to be produced once.
Delete the saved
file to redo the interactive fit.
<!---->
@note It is recommended to re-run the fitting macro after an interactive
session, as the saved parameters serve as much better starting points for the
minimizer and often produce a lower chi-squared on the second pass.
<!---->
Batch mode is temporarily disabled for the editor window and restored
afterward, so plot popups are not affected.
<!---->
## Simultaneous multi-channel fits
<!---->
Construct a `RooFitUtils` with the default constructor to put it in simultaneous
mode, then register channels and run one joint extended-likelihood fit across
all of them.
A channel is one spectrum (its own events, fit range, peaks, and
background) sharing a single observable with the others.
<!---->
- **`AddChannel(name, events, fit_range_low, fit_range_high, display_bin_width_kev, num_peaks, mu_inits, ...)`**
  — adds a channel.
  The component flags (`use_flat_background`, `use_step`,
  `use_low_exp_tail`, `use_low_lin_tail`, `use_high_exp_tail`) match the
  single-channel model, plus the per-channel locking options below.
- **`LinkParameter(target, source)`** /
  **`LinkPeakShape(target_channel, target_peak, source_channel, source_peak)`**
  — tie a parameter (or a peak's whole shape) in one channel to another so they
  are fit as a single shared degree of freedom. A peak's shape is its width and
  every step and tail parameter; its centroid and Gaussian yield stay free, so
  the same line can sit at a different position in each run (call
  `LinkParameter()` on `Mu` to tie the centroids too).
  Linked exponential tails
  are the same physical tail wherever the peaks sit; see
  @ref tail-amplitudes "Tail amplitudes".
  `LinkPeakShape()` used to link `Mu`
  as well; a `.simroofits` state saved then has no line for the target
  centroid, which on load starts from its source's value.
- **`AddFixedLines(channel, shape_peak, mus, gaus_yields)`** — add lines of
  known position and yield to a channel, drawn with one of its peaks' shapes;
  see @ref fixed-lines "Fixed lines".
- **`SeedChannel(channel_name, result)`** — supply a prior single-channel
  `FitResult` as starting values for that channel.
- **`FitSimultaneous(input_name, base_label)`** — builds the `RooSimultaneous`,
  applies the locks, and returns one `FitResult` per channel.
  With
  `SetInteractive()` the simultaneous editor opens after the automated fit, and
  accepted parameters are saved/reloaded as for the single-channel editors.
- **`DumpChannelCSV(channel, base_path, npts = 1000)`** — read-only export of a
  converged channel for external plotting: writes `<base>_spectrum.csv`
  (per-bin energy, data counts, fit total/background, residual pull) and
  `<base>_fitcurve.csv` (smooth `npts`-sample fit overlay).
  Evaluates the
  current parameter state and changes no fit state.
- **`SetFitDebug(kTRUE)`** — for simultaneous fits, un-suppresses RooFit
  evaluation errors and prints the seed NLL so an invalid-NLL failure names the
  offending PDF.
  Off by default.
<!---->
### Per-channel locking options
<!---->
All optional, on `AddChannel`:
<!---->
- **`mu_fixed`** — per-peak `std::vector<Bool_t>`; where `kTRUE`, that peak's
  centroid is held at its `mu_inits` value instead of floating.
- **`bkg_yield_fixed`** / **`bkg_slope_fixed`** — hold the channel's background
  yield and/or slope constant.
- **`lock_shape_after_seed`** — after seeding, fix sigma, the gaussian yield,
  and every tail/step shape parameter so only peak positions and normalizations
  float.
- **`shape_lock_per_peak`** — per-peak `std::vector<Bool_t>` refining
  `lock_shape_after_seed`: when non-empty, entry `i` controls whether peak `i`'s
  shape is locked (when empty, all peaks lock, as before).
- **`use_step_per_peak`** — per-peak `std::vector<Bool_t>` overriding the
  channel-level `use_step` for individual peaks (e.g. a step on some peaks but
  not on another peak sharing the channel).
<!---->
### Fixed lines {#fixed-lines}
<!---->
Some contributions to a window are known but cannot be fitted: the X-ray escape
peaks of a line outside the window, for example, sit under or beside the fitted
peaks with yields that follow from that line's yield.
`AddFixedLines(channel,
shape_peak, mus, gaus_yields)` adds them to a channel as lines at fixed
centroids with fixed Gaussian yields, drawn with the shape of one of the
channel's peaks:
<!---->
```cpp
RooFitUtils sim;
sim.AddChannel("run", events, 65, 80, 0.255, 2, pb_mus, kTRUE, kFALSE, kTRUE);
// Two escapes of a 95.9 keV line seen with 12000 events in its Gaussian,
// at 1.5% and 0.8% of it, with the shape of the Pb Ka1 peak (index 1).
std::vector<Double_t> mus = {95.9 - 23.17, 95.9 - 26.09};
std::vector<Double_t> yields = {0.015 * 12000, 0.008 * 12000};
sim.AddFixedLines("run", 1, mus, yields);
```
<!---->
- The lines take the shape peak's own `Sigma`, step and tail parameters,
  through any links on them, so they follow its shape as the fit moves it and
  add no free parameter.
- Their yields are Gaussian yields over the whole line, not inside the window
  as a fitted peak's `GausAmplitude` is; each line's in-window share follows
  from its centroid, the shape peak's sigma and the window. The tails are
  scaled as a fitted peak's are (see @ref tail-amplitudes "Tail amplitudes"),
  and a component the shape peak has switched off is left out.
- They are part of the channel's model, its chi-squared, its plots (one dashed
  dark cyan curve on the background) and the simultaneous editor's total curve.
  They are not peaks: they do not appear in the channel's `FitResult`, in saved
  `.simroofits` states or in the editor's parameter rows, so states saved
  without them load unchanged.
- Call it once per channel, after `AddChannel()` and before
  `FitSimultaneous()`; a second call for the same channel replaces the first.
<!---->
## References
<!---->
- Boggs SE, Pike SN.
Analytical fitting of gamma-ray photopeaks in germanium
  cross strip detectors.
  *Experimental Astronomy*.
  2023;56(2-3):403-420.
  doi: [10.1007/s10686-023-09914-8](https://doi.org/10.1007/s10686-023-09914-8).
- Longoria LC, Naboulsi AH, Gray PW, MacMahon TD.
Analytical peak fitting for
  gamma-ray spectrum analysis with Ge detectors.
  *Nuclear Instruments and
  Methods in Physics Research A*.
  1990;299(1-3):308-312.
  doi: [10.1016/0168-9002(90)90797-A](https://doi.org/10.1016/0168-9002\(90\)90797-A).
