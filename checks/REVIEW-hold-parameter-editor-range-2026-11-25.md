# Held parameters and the simultaneous editor's range — review checklist, 2026-11-25
<!---->
On `main` before ewt's review, which is in progress; the boxes below are still open.
Two changes to simultaneous fits, both from 73mGe's direct fit of its 01/15 no-shield
run (one channel: the Ge line with fixed escape lines on a linear continuum).
<!---->
**`HoldParameter(target, value)`.** That fit holds the Ge line's low tail at the
day's source-run tail while its width and yield float. Simultaneous mode could hold a
centroid (`mu_fixed`), the background (`bkg_yield_fixed`, `bkg_slope_fixed`) or a
peak's whole shape (`lock_shape_after_seed`, which also fixes `Sigma` and the Gaussian
yield), but not one parameter at a value set in code. The editor's Fix does not
survive a reload either: `LoadSimInteractiveParams` drops the saved flag on purpose
("Constness is model configuration, not saved state"). `SetManualParameter` does fix,
but on the single-channel model. A hold is stored by full name and applied after the
locks, at the cold start and again after a saved state loads.
<!---->
**The editor's curves after a range change.** Moving the fit range in the
simultaneous editor left the drawn total, its pulls and its chi2 on the old window:
the channel's `RooAddPdf` keeps integrals cached for the range it was first evaluated
in, so its `expectedEvents()` and `getVal()` stop matching its own components, which
follow the new range. Neither `fixCoefRange` under a new name, nor
`setShapeDirty()` on the observable, nor a copy of the sum (which shares the
components) cleared it. The fit itself was right: `fitTo` compiles a fresh graph. The
editor, and `ComputeChannelChi2` (the chi2 reported after Accept), now sum the leaf
components instead, which is what the component curves already drew.
<!---->
Line numbers are from `main` at the commit that adds this checklist. No code is
reproduced here; open the file at the line and review. Checks indicate review by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
- [ ] `src/RooFitUtils.cpp:2370-2377` — `HoldParameter`: refuses a target without
      the `channel:Param` form; a second call for one target replaces the first.
- [ ] `src/RooFitUtils.cpp:2379-2401` — `ApplyParameterHolds`: looks each target up
      in the build's name registry (a linked name maps to the variable it follows, so
      holding either name holds both), refuses a missing target or a value outside the
      variable's range rather than letting `setVal` clamp it, then sets the value and
      makes it constant. Check: failing rather than clamping.
- [ ] `src/RooFitUtils.cpp:3359-3363` — applied after the mu, background and shape
      locks at the cold start; a failure aborts the fit with an empty result, as a
      failed channel build does.
- [ ] `src/RooFitUtils.cpp:3416-3417` — applied again after a saved state loads. The
      loader skips constant parameters, so a held one is not loaded in the first place;
      the second call restates the hold. Its return is ignored because the targets
      were validated before the load.
- [ ] `src/RooFitUtils.cpp:434-453` — `PeakModelDensity` and `ChannelModelDensity`:
      each peak's and fixed line's Gaussian (on `GausYieldCoefficient()`, the
      coefficient the sum uses), step and three tails, plus the background, each
      normalised over the observable's current range. Check: that this is the same
      sum as the channel's `RooAddPdf` built in `BuildChannelModel`, component for
      component.
- [ ] `src/RooFitUtils.cpp:2967-3004` — `ComputeChannelChi2` takes the expected bin
      contents from `ChannelModelDensity` (its `peaks` and `bkg` arguments, unused
      until now, plus the channel's fixed lines) instead of the sum's
      `expectedEvents() * getVal()`. The parameter count still comes from the sum.
- [ ] `src/InteractiveSimultaneousFitEditor.cpp:529-539` — the editor's
      `ModelDensity`, and its three uses: the total curve (`:552`), the pulls (`:607`)
      and the chi2 (`:643`).
- [ ] `src/RooFitUtils.cpp:3468` — the editor's channel view now carries the
      channel's fixed lines.
<!---->
## Tooling
<!---->
- [ ] `include/RooFitUtils.hpp:1086-1102` — `HoldParameter` documentation.
- [ ] `include/RooFitUtils.hpp:689-711` — `ChannelModelDensity` documentation: why
      the components and not the sum.
- [ ] `include/RooFitUtils.hpp:457`, `:505` — the hold store and `ApplyParameterHolds`'
      declaration.
- [ ] `include/InteractiveSimultaneousFitEditor.hpp:51`, `:167` — the view's
      `fixed_lines` (null when there are none) and `ModelDensity`.
- [ ] `doc/pages/fitting.md:331` — `HoldParameter` in the simultaneous-fit list.
<!---->
## Verification
<!---->
- **The editor's stale total, reproduced and fixed.** A scratch replay of the
  editor's steps on AU's own objects (one channel, 400,000 flat events over
  55.25-92 keV plus 20,000 in a Gaussian at 68.75 keV, fitted over 65-80 keV): after
  the range moves to 61.75-80 keV and the editor's `fitTo` runs, the yields match the
  events in the new window (background 198,383 against 198,639 expected), the
  components give 13,852 per keV at 70 keV against the generated 13,827, and the sum
  gives 11,917, 14% low. `ChannelModelDensity` is that component sum.
- **`HoldParameter` in use.** 73mGe's 01/15 direct fit holds `LowExpTailAmplitude1`
  and `LowExpTailRatio1` at 0.5424 and 9.575: the editor shows both fixed, the saved
  state records them fixed, and a batch refit from the saved state returns them
  unchanged (Ge 68.7705 ± 0.0166 keV against the editor's 68.7704).
- **Builds.** `nix build` and the 73mGe shell built against this tree
  (`--override-input utils path:...`).
