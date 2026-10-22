# ROOT master pin and likelihood offsetting — review checklist, 2026-10-22
<!---->
Branch `build/root-master-cuda`, based on `main` at `66c895d`.
Not merged and not pushed.
Both library builds move from nixpkgs' ROOT 6.40.00 to ROOT master pinned at
`e0e5723692`, and every fit now minimises with `RooFit::Offset("initial")`.
As agreed with ewt, the mechanics of building ROOT master are summarised in one
item; the items with a `Check:` are the ones that change results or interfaces.
<!---->
Line numbers are from the working tree at the commit this checklist accompanies. No
code is reproduced here; open the file at the line and review. Checks indicate review
by ewt.
<!---->
## The RooFit bug this pin fixes
<!---->
In an extended `RooSimultaneous` fit with `Range(...)`, each channel's `RooAddPdf`
should interpret its yields relative to a defined range (the observable's default
range, or the range set with `fixCoefRange`). With RooFit's vectorized backends
(`EvalBackend::Cpu()` / `Cuda()`, which `BestAvailableBackend()` returns), the channels
lose that setting and fall back to RooFit's "automatic" interpretation, so the
minimised likelihood differs from the intended one by a parameter-dependent amount and
the minimum moves. `fixCoefRange` on the channel sums is silently ignored.
<!---->
Cause: `RooSimultaneous::compileForNormSet()` clones and prefixes each channel's
observables (`x` → `_ch0_x`) and only then fixes the coefficient normalisation, with
the unprefixed set, so the components cannot find their observables by name
([v6-40-00 L1253-L1255](https://github.com/root-project/root/blob/v6-40-00/roofit/roofitcore/src/RooSimultaneous.cxx#L1253-L1255),
[v6-40-04 L1317-L1319](https://github.com/root-project/root/blob/v6-40-04/roofit/roofitcore/src/RooSimultaneous.cxx#L1317-L1319),
[v6-38-02 L1253-L1255](https://github.com/root-project/root/blob/v6-38-02/roofit/roofitcore/src/RooSimultaneous.cxx#L1253-L1255)).
Introduced by [87ce8bfe9c](https://github.com/root-project/root/commit/87ce8bfe9c55dd4cac731d8f35edb85191a45dfc)
(PR #15795) and first released in 6.34.00; fixed in master by
[96d8760a77](https://github.com/root-project/root/commit/96d8760a776eaade609d4ad71e352b08fdc4b1c9)
(PR #23452, issue #23444). Backported to `v6-40-00-patches`, but in no release yet;
6.40.04 still has it.
<!---->
Affected here: `FitSimultaneous` fits (the automatic and refit-after-load paths, which
use `SplitRange`) in which at least one channel's window is narrower than the union of
all windows, which is `x_`'s default range. That is every such result since
simultaneous fitting was added on 2026-05-10 (`bb1721d`, ROOT 6.38.02). Not affected:
the interactive simultaneous editor, which sets the default range to its shared fit
range; simultaneous fits whose channels all use the same window; single-peak fits. On a
4 × 1M toy the bug moved the yields by +10–16% (60–90σ for the peak yields) and the
centroids by up to 0.5σ.
<!---->
## Core analysis logic — review these lines
<!---->
### Results change with the ROOT version
this commit
<!---->
- [x] The bug above: `FitSimultaneous` results whose channel windows differ from the
      union change with this pin.
      Check: you know which of your results this covers, and plan to re-run them.
- [x] `src/RooFitUtils.cpp:2879-2898` — `FixChannelCoefRanges` (from `c1f73f9`) takes
      effect only now; 6.40 ignored it inside the `RooSimultaneous`. Simultaneous yields
      therefore mean "events in the channel's window", and `master` without it would
      mean "events over the full default range".
      Check: downstream code that uses simultaneous yields as counts expects window
      counts.
<!---->
### Likelihood offsetting
this commit
<!---->
- [x] `src/RooFitUtils.cpp:567-577` (`RunFit`, with the explanation),
      `src/RooFitUtils.cpp:3002-3010` (refit after load), `src/RooFitUtils.cpp:3100-3108`
      (`FitSimultaneous`), `src/InteractiveRooFitEditor.cpp:797-802`,
      `src/InteractiveSimultaneousFitEditor.cpp:870-875` — `RooFit::Offset("initial")`
      added to every `fitTo`. The `AU_ROOFIT_FIT_DEBUG` seed-NLL print stays absolute,
      and the `FitSimultaneous` EDM gate is unchanged.
      Check: all five minimisations are covered and nothing else changed in them.
- [x] `include/FittingUtils.hpp:233-240` and `doc/pages/fitting.md:74-80` — `min_nll`
      is now the improvement over the fit's own start, not the absolute NLL.
      Check: nothing downstream compares `min_nll` between separate fits (the 73mGe
      convergence restart loop is the candidate).
<!---->
## Tooling
<!---->
- [x] `flake.nix:4-22, 46-140` — building ROOT master, summarised: the `root-src` and
      `clad-src` inputs; nixpkgs' derivation with the bundled LLVM/Clang (master needs
      LLVM 22) and nixpkgs' LLVM-20 clang dropped; all `builtins/` kept; clad patched
      in a writable copy, because master's patch step needs a git clone; `http` and
      `webgui` off, because master downloads civetweb and openui5 at build time, so
      `TBrowser`/`rootbrowse` use the classic X11 interface; `C_INCLUDE_DIRS` pointing
      the bundled Clang at glibc; no separate debug output; `rootWithCuda` adds only
      CUDA to `rootMaster`.
      Check: this is verified by the builds and the regression check below rather
      than line by line.
- [x] `flake.nix:320-321, 329-353`, `templates/standard/flake.nix:6,32,43`,
      `templates/python/flake.nix:6,32,44`, `doc/pages/python.md`,
      `doc/pages/installation.md` — new `packages.root`; the CPU dev shell and both
      templates use the flake's ROOT instead of `pkgs.root`; the templates no longer
      set `utils.inputs.nixpkgs.follows`.
      Check: downstream flakes must take `utils.packages.${system}.root` (or `rootCuda`)
      and drop that `follows`; a different ROOT is a different ABI.
- [x] `tests/CudaRegression.cpp` — the regression check now covers all three fit
      configurations, compares RooFit's NLL with a directly computed likelihood at
      three parameter points, re-records the references (only after that check passed),
      adds converged CPU/CUDA fits with offsetting on a model generated from itself,
      tightens the speed limit to 0.15 and mutes `InputArguments` as the library does.
      Check: compiled against 6.40 it fails the direct-likelihood check for
      `SplitRange`, so it catches the bug above.
- [x] `doc/pages/gpu.md:21-40` ("Which ROOT"), `:92-136` (performance, allocator
      note), `:137+` (regression check).
<!---->
## Test evidence
<!---->
- [x] `nix flake check --no-build` passes; `nix build` of `.#root`, `.#rootCuda`,
      `.#default`, `.#cuda`, `.#docs` and `.#cudaRegression` succeeds; the docs warning
      log is empty.
- [x] `nix run .#cuda-regression` exits 0 in 50 s: CPU/CUDA NLL agreement ≤ 2.2e-16;
      `(RooFit − direct)` constant to ≤ 6.4e-6 in all three configurations; references
      match; CUDA at 0.041 of the CPU's time per evaluation; offset fits converge
      (status 0, EDM < 1e-3) with CPU/CUDA parameters equal to 2e-5σ (single) and
      5e-5σ (`SplitRange`). The same source built against 6.40 fails the direct check
      for `SplitRange` (`(RooFit − direct)` varies by 1.96e4).
- [x] Offset check 1, library `FitSimultaneous` (3 × 100k events generated from the
      library's model), CPU vs CUDA: identical in every printed digit, status 0, EDM
      1.1e-5.
- [x] Offset check 2, library `FitSinglePeak` with and without offsetting (same ROOT):
      where the fit without offsetting reaches the full solution (2,000 and 20,000
      events) the results agree (20k: μ 59.98693 both, high tail 977.8 vs 978.4).
      Between 4,000 and 6,000 events both versions sometimes land in a second minimum
      with the tails at zero. From 7,000 events up, the fit without offsetting leaves
      the high tail at zero while the offset fit recovers it (the data has a 4% tail).
- [x] Offset check 3, library `FitSimultaneous`, 3 × 1M events: without offsetting
      status 3, EDM 0.091, "machine limit 1.65" (the 73mGe symptom); with it status 0,
      EDM 0.0019, machine limit 0.034. EDM < 1e-3 is not reached here, because the
      offset is taken at the start and this fit improves on it by 5.6e5; see Nits.
      Check: re-run the 73mGe Cu 01/14 post-cal simultaneous fit with
      `TAIL_RATIO_MAX = 20`, which is not available here.
<!---->
## Checked and excluded
<!---->
- Single-peak fits: correct under 6.40 and `master` alike (the direct check passes on
  both).
- The `KeepCudaAllocatorChunkMapped()` workaround from `c1f73f9` has no effect with
  master's batch compute, which no longer allocates per call. Kept for builds against
  released ROOT.
- ROOT master prints `RooArgSet::checkForDup` for the range-projection nodes of a
  `fixCoefRange`'d simultaneous model. The direct-likelihood check shows they do not
  change the NLL, and the library mutes that topic.
<!---->
## Nits
<!---->
- A second Migrad pass started from the first minimum would take the offset there,
  so that Minuit sees O(1) values and EDM < 1e-3 becomes reachable even when a fit
  starts far from its minimum. That is a further core change and is not made here.
