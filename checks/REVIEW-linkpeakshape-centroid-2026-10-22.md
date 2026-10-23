# LinkPeakShape leaves the centroid free — review checklist, 2026-10-22
<!---->
Branch `review/linkpeakshape-fitstatus-yieldratio`, based on `main` at `13c67d8`.
**Pushed to `main` on 2026-10-22 without ewt's review**, so that 73mGe can check it
in production. ewt reviews it on 2026-10-23; the boxes below are still open.
`LinkPeakShape` linked `Mu` along with the shape since `bb1721d` (2026-05-10),
contrary to its documentation ("leaving position and normalisation free"). It now
links the eight shape parameters only. `.simroofits` states saved while `Mu` was tied
hold only the source peak's centroid; on load, a target centroid with no line of its
own starts from its source's loaded value.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Why
<!---->
73mGe's `RunBkgGeSim` calls `LinkPeakShape("sig", i, "bkg", i)` to share the Pb Ka
shape between the background and the signal run, meaning to leave the signal run's
Pb centroids free. In every production pair the two channels' raw Pb centroids are
equal to 0.1 eV (01/13 Cu10: 73.1497 / 75.3414 in both), while independent fits of
the same runs differ by −50 to +30 eV because the reading moves with count rate. The
in-situ fit has been imposing one compromise centroid on both runs. 73mGe works
around it behind an opt-in switch with eight `LinkParameter` calls.
<!---->
## Core analysis logic — review these lines
<!---->
- [ ] `src/RooFitUtils.cpp:2320-2347` — `LinkPeakShape`: the list is `Sigma`,
      `StepAmplitude`, `LowExpTailAmplitude`, `LowExpTailRatio`,
      `LowLinTailAmplitude`, `LowLinTailSlope`, `HighExpTailAmplitude`,
      `HighExpTailRatio`; `Mu` is gone. The target/source centroid pair is recorded
      (`sim_shape_mu_pairs_`, `include/RooFitUtils.hpp:415-417`) for the loader and
      not linked.
      Check: no caller in this repository relies on the tie (none calls
      `LinkPeakShape`; see "Checked and excluded"), and 73mGe's `RunBkgGeSim` is the
      fit you want to change result.
- [ ] `src/RooFitUtils.cpp:1187`, `:1193` — `LoadSimInteractiveParams` records every
      parameter name the file has a line for, whether or not the line is used.
- [ ] `src/RooFitUtils.cpp:1219-1246` — after the values are loaded and before the
      legacy tail conversion: for each recorded pair, a target `Mu` that is its own
      variable, not constant (`mu_fixed` locks are applied before the load), and has
      no line, while its source's `Mu` was loaded, takes the source's value and
      error. A target `Mu` that the caller also tied with `LinkParameter` is the
      source variable and is skipped.
      Check: copying the error along with the value is fine; it only seeds the
      refit.
- [ ] `src/RooFitUtils.cpp:1256-1259` — one line, after "Loaded sim interactive
      params", naming each centroid started this way.
<!---->
## Tooling
<!---->
- [ ] `include/RooFitUtils.hpp:966-990` — `LinkPeakShape` documentation: which
      parameters are linked, that `Mu` and the Gaussian yield stay free, how to tie
      the centroids, and the old-state behaviour.
- [ ] `doc/pages/fitting.md:253-263` — the simultaneous-fit list says what a peak's
      shape is and how old states load.
<!---->
## Test evidence
<!---->
The check programs are scratch code outside the repository.
<!---->
- [ ] `nix build .#default` and `.#docs` succeed; the docs warning log is empty.
- [ ] Displacement test (requested item 4), CPU build: two channels with a Pb
      Ka2/Ka1-like doublet (72.80 and 75.00 keV, `sigma` 0.80, 120k and 200k Gaussian
      counts, one tail `R` 0.53, `tau/sigma` 11, 150k flat counts in 65-80 keV); the
      target channel's lines sit exactly 30 eV higher. Links as in `RunBkgGeSim`
      (one low tail tied to `bkg:` peak 1, the Ka widths linked, then
      `LinkPeakShape("sig", i, "bkg", i)`). Before the fix both channels' centroids
      are one variable and the fitted shift is 0 (5.5 and 9.1 standard errors from
      30 eV). After it, over 20 toys, the fitted shifts are 27.3 ± 1.9 eV (Ka2) and
      29.3 ± 1.0 eV (Ka1), with spreads 8.4 and 4.3 eV against median HESSE errors
      7.2 and 4.4 eV; 18 of 20 fits end with status 0. In the reference toy, every
      linked shape parameter of each target peak is the same `RooRealVar` pointer as
      its source's (`sigma`, `ratio_step`, both exponential tails' amplitude and ratio,
      the linear tail's amplitude and slope), while `mu` and `gaus_yield` are not;
      before the fix `mu` was shared as well.
- [ ] Old state: a `.simroofits` written as the tied version wrote it (one centroid
      line per Pb line, under `bkg:`), replayed without a refit, so the result is the
      start. Both `sig:` centroids equal the `bkg:` values from the file (72.8123 and
      75.0123), and the loader prints "centroids saved while LinkPeakShape tied them
      start from their source peak's: sig:Mu1 <- bkg:Mu1, sig:Mu2 <- bkg:Mu2".
- [ ] Every other fit path unchanged: the regression program from the linked-tails
      change (single, double, `link_sigma`, constrained double and triple,
      saved-state loads, 73mGe's production Ba fit), built against `13c67d8` and this
      tree, gives identical results to 17 digits, except its simultaneous case, which
      uses `LinkPeakShape`: there the two channels' centroids were equal (79.8004 and
      81.0208 keV in both) and are now 79.8695 / 79.7819 and 81.0220 / 81.0264 keV,
      with the minimum NLL 4.53 lower for the two freed centroids.
      Check: 73mGe's in-situ results change when it takes this build, by the
      compromise-centroid effect it has measured, and you want that.
<!---->
## Checked and excluded
<!---->
- Every use of `LinkPeakShape` in this repository: none in `src/`, `python/`,
  `templates/` or `tests/`; the documentation in `include/RooFitUtils.hpp` and
  `doc/pages/fitting.md` is updated. Older checklists mention it as history and are
  left as they are (`REVIEW-full-range-tail-amplitudes-2026-10-22.md` lists this bug
  as a nit).
- `SaveSimInteractiveParams` writes each variable once by pointer, so a state saved
  after this change has a line for each target centroid; nothing changes there.
- `SeedChannel` / `ApplySeedToChannel` set each peak's own `Mu`; a seeded target
  channel's centroids now start from its own seed instead of following the source.
- `ConstrainPeakSeparation` and the per-channel `mu_fixed` locks act on each
  channel's own `Mu` and are unaffected.
