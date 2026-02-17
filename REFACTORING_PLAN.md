# FittingUtils Refactoring Plan
<!---->
## Overview
<!---->
Unify the Standard/Detailed split into a single composable fitting model where you toggle individual components (step, low tail, high tail).
Add an exponential background model for fitting photopeaks on Compton edges.
Multi-peak methods always use full components (detailed) with linear background by default.
<!---->
---
<!---->
## 1.
Header File (`include/FittingUtils.hpp`)
<!---->
### New Types
<!---->
- [x] Add `BackgroundModel` enum:
  ```cpp
   ```
  - `kFlat`: constant B0 only
  - `kLinear`: B0 + B1*x
  - `kExponential`: B0 + A*exp(-lambda*x)
<!---->
- [x] Replace all 6 result structs with unified structs:
  ```cpp
  struct PeakFitResult {
      Float_t mu, mu_error;
      Float_t sigma, sigma_error;
      Float_t gaus_amplitude, gaus_amplitude_error;
      Float_t step_amplitude, step_amplitude_error;
      Float_t low_exp_tail_amplitude, low_exp_tail_amplitude_error;
      Float_t low_exp_tail_slope, low_exp_tail_slope_error;
      Float_t low_lin_tail_amplitude, low_lin_tail_amplitude_error;
      Float_t high_exp_tail_amplitude, high_exp_tail_amplitude_error;
      Float_t high_exp_tail_slope, high_exp_tail_slope_error;
  };
  ```
<!---->
- [x] Remove `FitResultStandard`, `FitResultDetailed`, `FitResultDoublePeakStandard`, `FitResultDoublePeakDetailed`, `FitResultTriplePeakStandard`, `FitResultTriplePeakDetailed`
<!---->
### FittingFunctions Namespace
<!---->
- [x] Add `ExponentialBackground` declaration
- [x] Replace `Standard` / `Detailed` with `PeakFunction`
- [x] Replace `DoublePeakStandard` / `DoublePeakDetailed` with `DoublePeak`
- [x] Replace `TriplePeakStandard` / `TriplePeakDetailed` with `TriplePeak`
- [x] Keep individual component functions: `Gaussian`, `LinearBackground`, `Step`, `LowTail`, `HighTail`
<!---->
### Class Members
<!---->
- [x] Remove `Bool_t isDetailed_`
- [x] Remove `Bool_t use_flat_background_`
- [x] Add `BackgroundModel bkg_model_`
- [x] Add `use_low_exp_tail_`, `use_low_lin_tail_`, `use_high_exp_tail_` (replaces `use_low_tail_`, `use_high_tail_`)
<!---->
### Constructor
<!---->
- [x] Change signature:
  ```cpp
<!---->
### Setters
<!---->
- [x] Replace `UseFlatBackground()` with `SetBackgroundModel(BackgroundModel)`
- [x] Replace `UseStep/UseLowTail/UseHighTail` with `SetStep()`, `SetLowExpTail()`, `SetLowLinTail()`, `SetHighExpTail()`
<!---->
### Public Fit Methods
<!---->
- [x] Replace `FitPeakStandard` + `FitPeakDetailed` with:
  ```cpp
  FitResult FitPeak(const TString input_name, const TString peak_name);
  ```
<!---->
- [x] Replace all `FitDoublePeakStandard` / `FitDoublePeakDetailed` overloads with:
  ```cpp
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakFitResult &constrained_peak, Double_t mu2_init);
  ```
<!---->
- [x] Replace `FitTriplePeakStandard` / `FitTriplePeakDetailed` with:
  ```cpp
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks, Double_t mu3_init);
  ```
<!---->
### Private Plot Methods
<!---->
- [x] Replace 6 plot methods with 3:
  ```cpp
  void PlotFit(const TString input_name, const TString peak_name);
  void PlotFitDoublePeak(const TString input_name, const TString peak_name);
  void PlotFitTriplePeak(const TString input_name, const TString peak_name);
  ```
<!---->
### Private Swap Methods
<!---->
- [x] Replace `SwapDoublePeakStandardParameters` + `SwapDoublePeakDetailedParameters` with:
  ```cpp
  void SwapDoublePeakParameters();  // always swaps 9-param blocks
  ```
<!---->
---
<!---->
## 2.
Source File (`src/FittingUtils.cpp`)
<!---->
### ExponentialBackground Function
<!---->
- [ ] Implement `ExponentialBackground`:
  ```cpp
  Double_t ExponentialBackground(Double_t *x, Double_t *par) {
      Double_t bkg_const = par[0];
      Double_t exp_amplitude = par[1];
      Double_t exp_decay = par[2];
      return bkg_const + exp_amplitude * TMath::Exp(-exp_decay * x[0]);
  }
  ```
<!---->
### Unified PeakFunction (replaces Standard + Detailed)
<!---->
- [ ] Implement `PeakFunction` with 14 parameters:
<!---->
  | Index | Parameter | Notes |
  |-------|-----------|-------|
  | 0 | Mu | |
  | 1 | Sigma | |
  | 2 | GausAmplitude | |
  | 3 | StepAmplitude | Fixed to 0 if step disabled |
  | 4 | LowExpTailAmplitude | Fixed to 0 if low exp tail disabled |
  | 5 | LowExpTailArg | Fixed if low exp tail disabled |
  | 6 | LowLinTailAmplitude | Fixed to 0 if low lin tail disabled |
  | 7 | LowLinTailSlope | Fixed if low lin tail disabled |
  | 8 | HighExpTailAmplitude | Fixed to 0 if high exp tail disabled |
  | 9 | HighExpTailArg | Fixed if high exp tail disabled |
  | 10 | BkgConst | Always present |
  | 11 | BkgSlope | Fixed to 0 if not kLINEAR |
  | 12 | ExpBkgAmplitude | Fixed to 0 if not kEXPONENTIAL |
  | 13 | ExpBkgDecay | Fixed if not kEXPONENTIAL |
<!---->
  Per-peak params (indices 0-8, 9 total): Mu, Sigma, GausAmp, StepAmp, LowExpTailAmp, LowExpTailArg, LowLinTailAmp, HighExpTailAmp, HighExpTailArg
<!---->
  Implementation sums: `Gaussian + Step + LowExpTail + LowLinTail + HighExpTail + LinearBackground + ExponentialBackground`
  (disabled components contribute 0 via fixed amplitudes)
<!---->
### Unified DoublePeak (replaces DoublePeakStandard + DoublePeakDetailed)
<!---->
- [ ] Implement `DoublePeak` with 20 parameters:
  - Params 0-8: Peak 1 (mu, sigma, gaus_amp, step_amp, low_exp_tail_amp, low_exp_tail_slope, low_lin_tail_amp, high_exp_tail_amp, high_exp_tail_slope)
  - Params 9-17: Peak 2 (same layout)
  - Params 18-19: BkgConst, BkgSlope (shared, always linear)
<!---->
### Unified TriplePeak (replaces TriplePeakStandard + TriplePeakDetailed)
<!---->
- [ ] Implement `TriplePeak` with 29 parameters:
  - Params 0-8: Peak 1
  - Params 9-17: Peak 2
  - Params 18-26: Peak 3
  - Params 27-28: BkgConst, BkgSlope (shared)
<!---->
### Constructor
<!---->
- [ ] Always build a 13-param TF1 using `PeakFunction`
- [ ] Set parameter names for all 13 params
- [ ] Initialize and set limits for Gaussian + background params
- [ ] Fix disabled component amplitudes to 0:
  - If `!use_step_`: fix param 3 to 0
  - If `!use_low_exp_tail_`: fix params 4,5
  - If `!use_low_lin_tail_`: fix param 6
  - If `!use_high_exp_tail_`: fix params 7,8
  - If `kFLAT`: fix param 10 (slope) to 0, fix params 11,12 to 0
  - If `kLINEAR`: fix params 11,12 to 0
  - If `kEXPONENTIAL`: fix param 10 (slope) to 0, set limits for params 11,12
<!---->
### FitPeak (replaces FitPeakStandard + FitPeakDetailed)
<!---->
- [ ] Implement sequential component-enabling logic:
  1.
  Start with only Gaussian + background (all optional component amplitudes fixed to 0)
  2.
  Fit to get baseline chi2
  3.
  If `use_step_`: release step param (3), fit, keep if chi2 improves
  4.
  If `use_low_exp_tail_`: release low exp tail params (4,5), fit, keep if chi2 improves
  5.
  If `use_low_lin_tail_`: release low lin tail param (6), fit, keep if chi2 improves
  6.
  If `use_high_exp_tail_`: release high exp tail params (7,8), fit, keep if chi2 improves
  6.
  When no components are enabled, steps 3-5 are skipped (equivalent to old Standard behavior)
- [ ] Return `FitResult` with one `PeakFitResult` entry
<!---->
### FitDoublePeak
<!---->
- [ ] Merge the 4 existing double peak methods into 2 overloads
- [ ] Always use full per-peak components (9 params per peak)
- [ ] Apply sequential component-enabling for each peak
- [ ] Return `FitResult` with 2 `PeakFitResult` entries
<!---->
### FitTriplePeak
<!---->
- [ ] Merge the 2 existing triple peak methods into 1
- [ ] Always use full per-peak components
- [ ] Return `FitResult` with 3 `PeakFitResult` entries
<!---->
### PlotFit (replaces PlotFitStandard + PlotFitDetailed)
<!---->
- [ ] Single method that conditionally draws each component
- [ ] Check if component amplitude != 0 before drawing it
- [ ] Draw exponential background curve when `bkg_model_ == kExponential`
- [ ] Same pad1/pad2 layout (histogram+fit on top, residuals on bottom)
<!---->
### PlotFitDoublePeak (replaces Standard + Detailed variants)
<!---->
- [ ] Merge into one method
- [ ] Always draw all per-peak components (Gaussian, step, low exp tail, low lin tail, high exp tail for each peak)
<!---->
### PlotFitTriplePeak (replaces Standard + Detailed variants)
<!---->
- [ ] Merge into one method
<!---->
### SwapDoublePeakParameters (replaces Standard + Detailed variants)
<!---->
- [ ] Always swap 9-parameter blocks (indices 0-8 with 9-17)
<!---->
### RegisterCustomFunctions
<!---->
- [ ] Update to register `PeakFunction` instead of `Standard`/`Detailed`
<!---->
### Keep As-Is
<!---->
- [ ] `EstimateBackground()`
- [ ] `ClampToBounds()`
- [ ] `SetManualParameters()`, `SetManualParameter()`, `ClearManualParameters()`
- [ ] Individual component functions: `Gaussian`, `LinearBackground`, `Step`, `LowTail`, `HighTail`
<!---->
---
<!---->
## 3.
Caller Migration Reference
<!---->
For code like the Calibration example:
```cpp
// Old:
fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
FitResultStandard result = fitter->FitPeakStandard(input_name, peak_name);
result.mu;
//
// New (equivalent — just Gaussian + linear background):
fitter = new FittingUtils(hist, fit_low, fit_high);
FitResult result = fitter->FitPeak(input_name, peak_name);
result.peaks[0].mu;
//
// New with exponential background:
fitter = new FittingUtils(hist, fit_low, fit_high, BackgroundModel::kEXPONENTIAL);
FitResult result = fitter->FitPeak(input_name, peak_name);
//
// New with step + all tails enabled:
fitter = new FittingUtils(hist, fit_low, fit_high, BackgroundModel::kLINEAR,
                          true, true, true, true);
FitResult result = fitter->FitPeak(input_name, peak_name);
```
<!---->
---
<!---->
## 4.
Verification
<!---->
- [ ] Project compiles cleanly
- [ ] Basic Gaussian+linear fit produces same results as old `FitPeakStandard`
- [ ] Full component fit produces same results as old `FitPeakDetailed`
- [ ] Exponential background fits converge on test data
- [ ] Multi-peak fits work correctly
- [ ] Plot output looks correct for all configurations
