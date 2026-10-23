# Fixed Gaussian yield ratio in FitDoublePeak — review checklist, 2026-10-22
<!---->
Branch `review/linkpeakshape-fitstatus-yieldratio`, on top of `b5b0cf3`.
**Pushed to `main` on 2026-10-22 without ewt's review**, so that 73mGe can check it
in production. ewt reviews it on 2026-10-23; the boxes below are still open.
Item 4 of the linked-tails request: `FitDoublePeak` gains `gaus_yield_ratio`. When
positive, `GausAmplitude1 = gaus_yield_ratio * GausAmplitude2` and only
`GausAmplitude2` is fitted. It works like the links: the first peak's
`RooFitPeakModel` holds the second peak's variable plus a fixed multiplier, and its
Gaussian coefficient and component yields are rebuilt on the product, so reporting,
saving and the editor follow without a dead placeholder.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Direction of the ratio
<!---->
The request says "fix peak 2's Gaussian yield as a ratio of peak 1's" and gives
"Ba 79.61/80.998 = 2.65/32.9 = 0.0805", which is peak 1 over peak 2 in
`FitDoublePeak`'s centroid order. The parameter follows the example: it is the lower
line's yield over the higher line's, so 73mGe passes `2.65 / 32.9`.
<!---->
## Core analysis logic — review these lines
<!---->
- [ ] `src/RooFitUtils.cpp:564-592` — `TieFirstGausYieldToSecond`: records
      `GausAmplitude1` for the loader; a constant `GausAmplitudeRatio`; the first
      peak's `gaus_yield` becomes `GausAmplitude2` with `gaus_yield_ratio` set;
      `gaus_yield_scaled` = `GausAmplitude2 * ratio`; the step, linear-tail and both
      exponential-tail yields are rebuilt on it (suffix `1_tied`), the exponential
      ones with the first peak's own `mu`, `sigma` and tail ratio.
      Check: the direction (above) and that the second peak's yields are untouched.
- [ ] `src/RooFitUtils.cpp:1697`, `:1729-1730` — `FitDoublePeak` calls it after the
      `link_sigma` / `link_tails` links and before `BuildTotalModel`, for any
      positive ratio.
- [ ] `src/RooFitUtils.cpp:418-425`, `:504` — `BuildExpTailYields` and
      `BuildTotalModel` build on `GausYieldCoefficient()`, which is `gaus_yield`
      itself unless the multiplier is set, so every other model is the same graph.
- [ ] `src/RooFitUtils.cpp:945-946`, `:2897-2898`, `:1356` — `ExtractPeakResult`,
      `ExtractPeakResultFor` and `AppendPeakGraphs` read `GausYield()` /
      `GausYieldError()`: `getVal()` and `getError()` times the multiplier, which is
      exactly 1.0 for every other peak.
      Check: the first peak's reported error is the second's error times the ratio,
      the right propagation for a fixed ratio.
- [ ] `src/RooFitUtils.cpp:1145` — the loader's "ignored" line now says "linked to
      the other peak's", since a tied yield is the second peak's. An old state's
      `GausAmplitude1` line is skipped; its `GausAmplitude2` line seeds the fit.
<!---->
## Tooling
<!---->
- [ ] `include/RooFitUtils.hpp:200-204`, `:242-261` — `RooFitPeakModel`: the
      multiplier, the scaled coefficient, `GausYieldCoefficient()`, `GausYield()`,
      `GausYieldError()`, and the rule to read yields through them.
- [ ] `include/RooFitUtils.hpp:444-447`, `:514` — name-set comment; private
      declaration.
- [ ] `include/RooFitUtils.hpp:894-935` — `FitDoublePeak` documentation: the
      parameter, its direction with the Ba example, when it equals an intensity
      ratio, the start it needs, what is reported, the saved state.
- [ ] `src/InteractiveRooFitEditor.cpp:257-263`, `:514`;
      `src/InteractiveSimultaneousFitEditor.cpp:555` — both editors draw the
      Gaussian curve from `GausYield()`; the single-channel editor labels the tied
      peak's yield row "Yield x <ratio>", since that row edits `GausAmplitude2`.
      Check: the label is clear enough; GUI code, compiled but not driven.
- [ ] `doc/pages/fitting.md:190-246` — "Linked doublets": the third option, its
      direction and start, reporting and saved states.
<!---->
## Test evidence
<!---->
- [ ] `nix build` of `.#default`, `.#cuda` and `.#docs` succeeds; the docs warning
      log is empty; `clang-format` leaves every touched file unchanged.
- [ ] Default calls unchanged: the regression program against the status change
      (built at `3f39025`; `b5b0cf3` only adds the `fit_retried` field, which the
      program does not print) gives identical results to 17 digits in all 298 result
      and parameter lines.
- [ ] Synthetic Ba-like toys (CUDA): one shared tail, Gaussian counts 13363 and
      166000 (true ratio 0.0805), fitted with `link_sigma` + `link_tails` +
      `gaus_yield_ratio` = 0.0805, against the free ratio on the same toys.
      - From a start with a realistic tail (`R` 0.55, `tau/sigma` 10), 100 toys: all
        100 in the right minimum, 95 status 0. `mu2` bias +0.1 ± 0.3 eV with spread
        3.4 eV and median HESSE error 4.1 eV, against +2.5 ± 1.2 eV, 12.5 eV and
        11.8 eV with the ratio free; `mu1` spread 57 eV against 107; `k` 0.0701
        (truth 0.0700).
      - From the neutral start used before (`R` 0.3, `tau/sigma` 5), 200 toys: 134
        in the right minimum, with the same gains (`mu2` bias −0.1 ± 0.3 eV, spread
        3.8 eV, against 12.7 eV free). 66 converge directly, without the retry, to a
        local minimum with the first peak at 76.25 keV and reduced chi-squared
        11-12: the weak line, whose yield can no longer adjust, slides into the tail
        region the too-short start leaves unexplained. All 68 retried fits found the
        right minimum.
      Check: documenting the start requirement (realistic tail, for example a saved
      `link_tails` fit) is enough, with chi2/ndf flagging the wrong minimum, rather
      than constraining `mu1` or fitting the free ratio first inside `FitDoublePeak`.
- [ ] 73mGe Ba-133 01/15 with `gaus_yield_ratio` = 2.65 / 32.9 (CUDA build):
      - From the production state: `mu2` 81.269188 ± 0.0043 keV, `k` 0.0679,
        `tau/sigma` 10.99 ± 0.91, chi2/ndf 1.488 (1.503 with the ratio free), status 0,
        `cov_qual` 3. Against the free ratio, `mu2` moves −5.5 eV and its error halves
        (8.3 → 4.3 eV); against the shared-tail binned refit's 81.2776, it is 8.4 eV
        lower (the free-ratio fit is 2.9 eV lower).
      - From the neutral start: the wrong minimum (first peak at 76.35 keV,
        chi2/ndf 14.9), reached through the retry, with `fit_retried` set.
      Check: which `mu2` the binned comparison should be against: 81.2776 came from
      your shared-tail refit, and the request does not say whether the ratio was
      fixed there.
- [ ] Reporting and saved states: the first peak's `gaus_amplitude` / `gaus_amplitude`
      ratio is 0.080547 (= 2.65 / 32.9) and its error is exactly 0.0805 times the
      second's (242.42 both ways). The saved state has a `GausAmplitude2` line and
      no `GausAmplitude1`; loading the production state skips 6 lines as linked (5
      from `link_sigma` + `link_tails`, plus `GausAmplitude1`).
<!---->
## Checked and excluded
<!---->
- The constrained `FitDoublePeak` overload, `FitTriplePeak` and simultaneous fits do
  not take a ratio; their peaks keep a multiplier of 1.
- `CollectAllParams` keeps ten entries per peak; the first peak's yield slot is
  `GausAmplitude2`, written once by the saver and counted once in chi2/ndf (the
  linked-tails change's `DistinctParams`).
- `ConfigureComponentFlagsForPeak`, `FixComponent` and `ReleaseComponent` do not touch
  yields; the automatic search is unchanged.
- Peaks never swap in `FitDoublePeak` (each centroid is bounded at the midpoint of
  the starting values), so the multiplier stays with the lower line.
