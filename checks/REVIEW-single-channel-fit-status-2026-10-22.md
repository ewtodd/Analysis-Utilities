# Single-channel fit status and retry — review checklist, 2026-10-22
<!---->
Branch `review/linkpeakshape-fitstatus-yieldratio`, on top of `b6140c2`.
**Pushed to `main` on 2026-10-22 without ewt's review**, so that 73mGe can check it
in production. ewt is reviewing it (from 2026-10-24); the boxes below are still open.
A single-channel RooFitUtils refit can stop at its starting point (status 3), and the
saved-state path returned that start as a valid result, with no way for a downstream
macro to notice. Both remedies you offered are in: `RunFit` retries a minimisation
that ends without a valid minimum, and every single-channel `FitResult` now carries
the Minuit2 status, covariance quality, EDM, minimum NLL and per-parameter
diagnostics of the fit it comes from, and whether that fit needed the retry.
<!---->
Line numbers are from `main` at `5b82845` (the code as at `69098de`, where every open
checklist pointed since the review branches went onto `main` on 2026-10-24). No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Why both
<!---->
The retry removes the stall where it occurs: in the linked-tails tests the library's
own refit stalled in 200 of 200 shared-tail toys from a neutral start and on 73mGe's
Ba state, and a Strategy 1 pass recovered nearly all of them. The stored status is
still needed, because a retry can also fail, the editor's Refit does not go through
`RunFit`, and a downstream macro should not have to infer convergence from chi2/ndf.
A retry can also succeed in the wrong place: from a poor start the Strategy 1 pass
can reach another local minimum, returned with status 0.
`fit_retried` marks every
retried result so that its chi2/ndf gets looked at.
<!---->
## Core analysis logic — review these lines
<!---->
- [x] `src/RooFitUtils.cpp:729-741` — `Minimize`: the `fitTo` call `RunFit` made
      before, with the strategy as an argument.
      Strategy 2 reproduces the old call
      argument for argument.
- [x] `src/RooFitUtils.cpp:743-764` — `RunFit`: Strategy 2 as before; if that ends
      with status `>= kFirstInvalidMinimumStatus` (3: EDM above maximum, which
      includes stopping at the start; 4: call limit; 5: other) or returns nothing, it
      prints a warning, runs a Strategy 1 pass from where it stopped, and runs the
      Strategy 2 fit again, whose result it returns, and records that it retried
      (`last_fit_retried_`, reset at every call and in `InitState`, `:68`).
      Statuses 1 and 2 (a valid minimum whose covariance was forced positive
      definite, or whose HESSE failed) are not retried.
      Check: the threshold.
      Retrying 1 and 2 as well would change many converged
      fits (150 of 200 separate-tail toys end with status 1) for little gain.
- [x] `src/RooFitUtils.cpp:766-784` — `StoreFitDiagnostics`: status, `covQual`, EDM,
      `minNll`, whether the fit was retried, and `ExtractParameterDiagnosticsSingle`
      into the `FitResult`, and a warning when the status is still 3 or more.
      It is
      called right after the `RunFit` whose result it stores, so the retry flag is
      that fit's.
- [x] `src/RooFitUtils.cpp:1642`, `:1789`, `:2023`, `:2175` (saved-state refits),
      `:1682`, `:1825`, `:2053`, `:2204` (initial fit failed, ending the search),
      `:1703`, `:1921`, `:2074`, `:2225` (final fit) — the four `Fit*` methods store
      the diagnostics of the minimisation their result comes from. `valid` is not
      changed: a saved-state refit is still `valid` whatever its status.
      Check: you want `valid` left as it is, with `fit_status` the indicator. Making
      a saved-state refit invalid from status 3 would also blank its result (every
      field -1) for any downstream that does not look.
      ewt: that sure is a lot of repetitive code.
- [x] `src/RooFitUtils.cpp:3186-3195` — `ExtractParameterDiagnosticsSingle`, which
      had no caller, now lists each variable once under its own name, so a linked
      parameter is not listed twice and the names survive `SortPeaksByMu`.
<!---->
## Tooling
<!---->
- [x] `include/RooFitUtils.hpp:615-625` — public `kFirstInvalidMinimumStatus` = 3,
      for downstream macros to compare `fit_status` against.
- [x] `include/RooFitUtils.hpp:396-403` — class note: what a single-channel result
      carries, the retry, and why a retried result's chi2/ndf needs a look.
- [x] `include/RooFitUtils.hpp:540-546` — private `Minimize`, `last_fit_retried_`,
      `StoreFitDiagnostics`.
- [x] `include/FittingUtils.hpp:210-211`, `:242-244`, `:256-271` — `FitResult` and
      `FitParameterDiagnostic` documentation: which fits fill which fields;
      `fit_status`'s meaning.
- [x] `include/FittingUtils.hpp:287-295` — new `FitResult::fit_retried`, default
      `kFALSE`, never set by simultaneous or FittingUtils fits.
      Check: adding a field to `FitResult` is fine for downstream code (it is a
      plain struct; 73mGe rebuilds against each library version).
- [x] `doc/pages/fitting.md:74-108` — "Results and diagnostics": what each fit fills,
      the retry, and `fit_retried`.
<!---->
## Test evidence
<!---->
- [x] `nix build` of `.#default`, `.#cuda` and `.#docs` succeeds; the docs warning
      log is empty.
- [x] Regression program (CPU) against `b6140c2`: the single, `link_sigma` double,
      constrained double and production Ba fits (A, C, D, H) gain only their status
      line, with every value identical to 17 digits; one of A's component tests
      retried and was rejected as before.
      Changed:
      - G, a `link_sigma` doublet refit from a pre-convention state, which stalled:
        chi2/ndf 22.06 → 1.147, `mu2` 81.0000 (the start) → 80.9914 keV (truth
        80.9979), status 1, `cov_qual` 3.
      - B, the automatic double fit with free widths on a doublet with heavy tails:
        its initial fit failed (status 3) and the fit returned `valid = 0`.
        Now the
        retry converges, but to a wrong minimum: peak 1 a 2.28 keV-wide Gaussian at
        77.74 keV standing in for the tails, final chi2/ndf 15.4, status 0.
        E (the
        constrained triple seeded from B) and F (the simultaneous fit seeded from
        B) change with it.
      Check: you accept that the retry can turn an outright failure into a valid but
      wrong local minimum, which chi2/ndf and `fit_retried` flag and `valid` and
      `fit_status` do not.
      The yield-ratio change on top of this one shows the same
      on real data from a neutral start.
- [x] `fit_retried` (CPU): set on G, whose refit stalled and was retried (status 1,
      chi2/ndf 1.147); clear on C, whose final fit converged at once.
- [x] 73mGe Ba-133 01/15 with `link_sigma` + `link_tails` from the unedited
      production state (CUDA build): the library's refit now converges (it stalled
      before): `mu2` 81.274712, status 0, `cov_qual` 3, matching the neutral and
      hand-edited starts to 0.08 eV. The production fit (separate tails) is
      unchanged at `mu2` 81.264336 and reports status 0, `cov_qual` 3.
- [x] Synthetic Ba-like toys (50, CUDA build): the library's own shared-tail refit
      stalled in 0 of 50 (200 of 200 before); the retry fired in all 50 and reached a
      valid minimum in all 50.
- [x] Reporting: every single-channel result above has `has_fit_diagnostics` set,
      and `parameter_diagnostics` has 21 entries for the production Ba fit (22
      variables, `Sigma1` shared) and 17 with `link_tails`.
<!---->
## Found while testing, not changed
<!---->
- The automatic component search accepts a group only when its fit has status
  exactly 0 (`TestLowSideGroup`, `TestHighTailIndependent`, the inter-peak group).
  A group fit that ends with status 1 is rejected however much it improves chi2/ndf.
  In the regression program's single peak with a large low tail, the low-side group
  is rejected and the result keeps chi2/ndf 64; releasing that group after the
  search and refitting gives status 1, `cov_qual` 2, chi2/ndf 1.26.
- `near_limit` fires for disabled components held constant at their bounds (step,
  linear tail and high tail at 0, tail ratio at 1), as it does in the simultaneous
  path; filter on `is_constant`. Its margin is 0.1% of the range, so `R` = 2.997
  against a limit of 3 is not flagged.
- The editor's Refit runs its own `fitTo` without the retry.
- `FitSimultaneous` already stores its status; it has no retry.
