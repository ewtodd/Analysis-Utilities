# Fit convergence controls and correlation diagnostics — review checklist, 2026-10-19
<!---->
Branch `fit-convergence-controls`, based on `main` at `781932d`. Not merged and not
pushed. This change exists so the 73mGe convergence study can (a) restart a
simultaneous fit from its own endpoint, (b) switch Minuit2 strategy, call limit and
final HESSE, and (c) read the fitted parameters' limits and strong correlations.
Every new control defaults to the exact previous behaviour, so existing callers are
unaffected.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
### Simultaneous-fit minimizer controls
this commit
<!---->
- [ ] `include/RooFitUtils.hpp:297-305` — new private members `save_after_refit_ =
      kFALSE`, `fit_strategy_ = 1`, `fit_max_calls_ = 0`, `fit_final_hesse_ =
      kFALSE`.
      Check: these four defaults together reproduce the previous `fitTo`
      call exactly — Strategy 1, no MaxCalls argument, no Hesse argument, and no
      write-back after a refit.
- [ ] `include/RooFitUtils.hpp:560-608` — `SetSaveAfterRefit`, `SetFitStrategy`,
      `SetFitMaxCalls`, `SetFitFinalHesse` and their Doxygen blocks.
      Check: the setters only assign the members above; none of them can change a fit
      unless explicitly called.
- [ ] `src/RooFitUtils.cpp:2943-2950` — `max_calls_arg` and `hesse_arg` are built
      with `RooCmdArg::none()` unless the corresponding control is set.
      Check: `MaxCalls(0)` is never passed, and `Hesse(false)` is never passed; an
      unset control passes exactly nothing to `fitTo`.
- [ ] `src/RooFitUtils.cpp:2974-2981` (refit-after-load path) and
      `src/RooFitUtils.cpp:3067-3076` (cold path) — both `fitTo` calls now pass
      `Strategy(fit_strategy_)`, `max_calls_arg` and `hesse_arg`.
      Check: the two paths still use identical fit configuration, which was the
      invariant established by `3523967`. With the defaults set, the argument list
      must be the same as before this change.
- [ ] `src/RooFitUtils.cpp:2997-3000` — with `SetSaveAfterRefit(kTRUE)`, a refit from
      a saved `.simroofits` writes its endpoint back to that same file.
      Check: with the default `kFALSE` there is no write; with `kTRUE` the write
      happens after `fit_result` exists, including a fit the EDM gate rejects. That
      is deliberate: the restart loop records a non-converged endpoint so the next
      iteration can start from it.
<!---->
### Parameter diagnostics and correlations
this commit
<!---->
- [ ] `include/FittingUtils.hpp:181-192` — `FitParameterDiagnostic::is_constant`.
      Check: the field is additive and initialised `kFALSE`; existing consumers see
      the same values in all pre-existing fields.
- [ ] `include/FittingUtils.hpp:194-207` (new `FitParameterCorrelation`) and
      `include/FittingUtils.hpp:239-243` (`FitResult::parameter_correlations`).
      Check: an empty vector is the sentinel for "no pair above threshold", which is
      also what non-simultaneous fits return.
- [ ] `src/RooFitUtils.cpp:2628` — `RooFitExtractParamDiagnostic` records
      `rv->isConstant()`.
      Check: constness is reported only; the near-limit flags and the parameter list
      are unchanged.
- [ ] `src/RooFitUtils.cpp:3098-3123` — pairs are read from
      `RooFitResult::correlationMatrix()` at indices matching `floatParsFinal()`,
      keeping `|rho| >= 0.95`; the matrix dimension is compared with the parameter
      list size first.
      Check: the index order of `correlationMatrix()` is indeed `floatParsFinal()`
      (the RooFit convention), and a failed or empty covariance cannot index out of
      bounds.
- [ ] `src/RooFitUtils.cpp:3163` — the same correlation list is attached to every
      channel's `FitResult`.
      Check: this is a global (single combined minimisation) diagnostic, so
      duplicating it per channel is intended rather than a stale per-channel copy.
<!---->
## Test evidence
<!---->
- [ ] `nix build .#cuda` succeeds at the commit under review.
- [ ] The 73mGe convergence macro's first iteration, run against this build with all
      controls at their defaults, reproduces the production pre-cal fit's Minuit
      diagnostics (status, covQual, EDM, NLL) and raw Ge centroid.
      Check: if this does not reproduce exactly, the default-preservation claim is
      false and the convergence study's baseline is not the production fit.
<!---->
## Checked and excluded
<!---->
- The 73mGe-side macro (`ConvergenceTest.cpp`) and the combiner changes are not part
  of this repository and are not part of this review.
- `FitResult` gains fields, so this is an ABI change. The library and all consumers
  are rebuilt together by the flake; no serialized `FitResult` is persisted.
<!---->
## Nits
<!---->
- `include/FittingUtils.hpp:214-217` (pre-existing) documents `cov_qual` as
  "`0` Exact, `1` NotPositiveDefinite, `2` Approximate, `3` External", while the
  EDM gate treats `covQual >= 2` as acceptable and the 73mGe convergence criterion
  asks for 3. Left unchanged because this branch does not touch that contract.
