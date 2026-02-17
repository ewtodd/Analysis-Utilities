# FittingUtils Refactoring Plan

## Overview

Unify the Standard/Detailed split into a single composable fitting model where you toggle individual components (step, low tail, high tail). Add an exponential background model for fitting photopeaks on Compton edges. Multi-peak methods always use full components (detailed) with linear background by default.

---

## 1. Header File (`include/FittingUtils.hpp`)

### New Types

- [ ] Add `BackgroundModel` enum:
  ```cpp
  enum class BackgroundModel { kFlat, kLinear, kExponential };
  ```
  - `kFlat`: constant B0 only
  - `kLinear`: B0 + B1*x
  - `kExponential`: B0 + A*exp(-lambda*x)

- [ ] Replace all 6 result structs with unified structs:
  ```cpp
  struct PeakResult {
      Float_t mu, mu_error;
      Float_t sigma, sigma_error;
      Float_t gaus_amplitude, gaus_amplitude_error;
      Float_t step_amplitude, step_amplitude_error;
      Float_t low_tail_amplitude, low_tail_amplitude_error;
      Float_t low_tail_slope, low_tail_slope_error;
      Float_t high_tail_amplitude, high_tail_amplitude_error;
      Float_t high_tail_slope, high_tail_slope_error;
  };

  struct FitResult {
      std::vector<PeakResult> peaks;  // 1, 2, or 3 entries
      Float_t bkg_const, bkg_const_error;
      Float_t bkg_slope, bkg_slope_error;
      Float_t exp_bkg_amplitude, exp_bkg_amplitude_error;
      Float_t exp_bkg_decay, exp_bkg_decay_error;
      Float_t reduced_chi2;
      Bool_t valid;
  };
  ```

- [ ] Remove `FitResultStandard`, `FitResultDetailed`, `FitResultDoublePeakStandard`, `FitResultDoublePeakDetailed`, `FitResultTriplePeakStandard`, `FitResultTriplePeakDetailed`

### FittingFunctions Namespace

- [ ] Add `ExponentialBackground` declaration
- [ ] Replace `Standard` / `Detailed` with `PeakFunction`
- [ ] Replace `DoublePeakStandard` / `DoublePeakDetailed` with `DoublePeak`
- [ ] Replace `TriplePeakStandard` / `TriplePeakDetailed` with `TriplePeak`
- [ ] Keep individual component functions: `Gaussian`, `LinearBackground`, `Step`, `LowTail`, `HighTail`

### Class Members

- [ ] Remove `Bool_t isDetailed_`
- [ ] Remove `Bool_t use_flat_background_`
- [ ] Add `BackgroundModel bkg_model_`

### Constructor

- [ ] Change signature:
  ```cpp
  // Before:
  FittingUtils(TH1 *hist, Float_t low, Float_t high,
               Bool_t use_flat_background, Bool_t isDetailed,
               Bool_t use_step = kTRUE, Bool_t use_low_tail = kTRUE,
               Bool_t use_high_tail = kTRUE);

  // After:
  FittingUtils(TH1 *hist, Float_t low, Float_t high,
               BackgroundModel bkg_model = BackgroundModel::kLinear,
               Bool_t use_step = kFALSE, Bool_t use_low_tail = kFALSE,
               Bool_t use_high_tail = kFALSE);
  ```

### Setters

- [ ] Replace `UseFlatBackground()` with `SetBackgroundModel(BackgroundModel)`
- [ ] Keep `UseStep()`, `UseLowTail()`, `UseHighTail()`

### Public Fit Methods

- [ ] Replace `FitPeakStandard` + `FitPeakDetailed` with:
  ```cpp
  FitResult FitPeak(const TString input_name, const TString peak_name);
  ```

- [ ] Replace all `FitDoublePeakStandard` / `FitDoublePeakDetailed` overloads with:
  ```cpp
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakResult &constrained_peak, Double_t mu2_init);
  ```

- [ ] Replace `FitTriplePeakStandard` / `FitTriplePeakDetailed` with:
  ```cpp
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks, Double_t mu3_init);
  ```

### Private Plot Methods

- [ ] Replace 6 plot methods with 3:
  ```cpp
  void PlotFit(const TString input_name, const TString peak_name);
  void PlotFitDoublePeak(const TString input_name, const TString peak_name);
  void PlotFitTriplePeak(const TString input_name, const TString peak_name);
  ```

### Private Swap Methods

- [ ] Replace `SwapDoublePeakStandardParameters` + `SwapDoublePeakDetailedParameters` with:
  ```cpp
  void SwapDoublePeakParameters();  // always swaps 8-param blocks
  ```

---

## 2. Source File (`src/FittingUtils.cpp`)

### ExponentialBackground Function

- [ ] Implement `ExponentialBackground`:
  ```cpp
  Double_t ExponentialBackground(Double_t *x, Double_t *par) {
      Double_t bkg_const = par[0];
      Double_t exp_amplitude = par[1];
      Double_t exp_decay = par[2];
      return bkg_const + exp_amplitude * TMath::Exp(-exp_decay * x[0]);
  }
  ```

### Unified PeakFunction (replaces Standard + Detailed)

- [ ] Implement `PeakFunction` with 12 parameters:

  | Index | Parameter | Notes |
  |-------|-----------|-------|
  | 0 | Mu | |
  | 1 | Sigma | |
  | 2 | GausAmplitude | |
  | 3 | StepAmplitude | Fixed to 0 if step disabled |
  | 4 | LowTailAmplitude | Fixed to 0 if low tail disabled |
  | 5 | LowTailSlope | Fixed if low tail disabled |
  | 6 | HighTailAmplitude | Fixed to 0 if high tail disabled |
  | 7 | HighTailSlope | Fixed if high tail disabled |
  | 8 | BkgConst | Always present |
  | 9 | BkgSlope | Fixed to 0 if not kLinear |
  | 10 | ExpBkgAmplitude | Fixed to 0 if not kExponential |
  | 11 | ExpBkgDecay | Fixed if not kExponential |

  Implementation sums: `Gaussian + Step + LowTail + HighTail + LinearBackground + ExponentialBackground`
  (disabled components contribute 0 via fixed amplitudes)

### Unified DoublePeak (replaces DoublePeakStandard + DoublePeakDetailed)

- [ ] Implement `DoublePeak` with 18 parameters:
  - Params 0-7: Peak 1 (mu, sigma, gaus_amp, step_amp, low_tail_amp, low_tail_slope, high_tail_amp, high_tail_slope)
  - Params 8-15: Peak 2 (same layout)
  - Params 16-17: BkgConst, BkgSlope (shared, always linear)

### Unified TriplePeak (replaces TriplePeakStandard + TriplePeakDetailed)

- [ ] Implement `TriplePeak` with 26 parameters:
  - Params 0-7: Peak 1
  - Params 8-15: Peak 2
  - Params 16-23: Peak 3
  - Params 24-25: BkgConst, BkgSlope (shared)

### Constructor

- [ ] Always build a 12-param TF1 using `PeakFunction`
- [ ] Set parameter names for all 12 params
- [ ] Initialize and set limits for Gaussian + background params
- [ ] Fix disabled component amplitudes to 0:
  - If `!use_step_`: fix param 3 to 0
  - If `!use_low_tail_`: fix params 4,5
  - If `!use_high_tail_`: fix params 6,7
  - If `kFlat`: fix param 9 (slope) to 0, fix params 10,11 to 0
  - If `kLinear`: fix params 10,11 to 0
  - If `kExponential`: fix param 9 (slope) to 0, set limits for params 10,11

### FitPeak (replaces FitPeakStandard + FitPeakDetailed)

- [ ] Implement sequential component-enabling logic:
  1. Start with only Gaussian + background (all optional component amplitudes fixed to 0)
  2. Fit to get baseline chi2
  3. If `use_step_`: release step param (3), fit, keep if chi2 improves
  4. If `use_low_tail_`: release low tail params (4,5), fit, keep if chi2 improves
  5. If `use_high_tail_`: release high tail params (6,7), fit, keep if chi2 improves
  6. When no components are enabled, steps 3-5 are skipped (equivalent to old Standard behavior)
- [ ] Return `FitResult` with one `PeakResult` entry

### FitDoublePeak

- [ ] Merge the 4 existing double peak methods into 2 overloads
- [ ] Always use full per-peak components (8 params per peak)
- [ ] Apply sequential component-enabling for each peak
- [ ] Return `FitResult` with 2 `PeakResult` entries

### FitTriplePeak

- [ ] Merge the 2 existing triple peak methods into 1
- [ ] Always use full per-peak components
- [ ] Return `FitResult` with 3 `PeakResult` entries

### PlotFit (replaces PlotFitStandard + PlotFitDetailed)

- [ ] Single method that conditionally draws each component
- [ ] Check if component amplitude != 0 before drawing it
- [ ] Draw exponential background curve when `bkg_model_ == kExponential`
- [ ] Same pad1/pad2 layout (histogram+fit on top, residuals on bottom)

### PlotFitDoublePeak (replaces Standard + Detailed variants)

- [ ] Merge into one method
- [ ] Always draw all per-peak components (Gaussian, step, tails for each peak)

### PlotFitTriplePeak (replaces Standard + Detailed variants)

- [ ] Merge into one method

### SwapDoublePeakParameters (replaces Standard + Detailed variants)

- [ ] Always swap 8-parameter blocks (indices 0-7 with 8-15)

### RegisterCustomFunctions

- [ ] Update to register `PeakFunction` instead of `Standard`/`Detailed`

### Keep As-Is

- [ ] `EstimateBackground()`
- [ ] `ClampToBounds()`
- [ ] `SetManualParameters()`, `SetManualParameter()`, `ClearManualParameters()`
- [ ] Individual component functions: `Gaussian`, `LinearBackground`, `Step`, `LowTail`, `HighTail`

---

## 3. Caller Migration Reference

For code like the Calibration example:
```cpp
// Old:
fitter = new FittingUtils(hist, fit_low, fit_high, kFALSE, kFALSE);
FitResultStandard result = fitter->FitPeakStandard(input_name, peak_name);
result.mu;

// New (equivalent — just Gaussian + linear background):
fitter = new FittingUtils(hist, fit_low, fit_high);
FitResult result = fitter->FitPeak(input_name, peak_name);
result.peaks[0].mu;

// New with exponential background:
fitter = new FittingUtils(hist, fit_low, fit_high, BackgroundModel::kExponential);
FitResult result = fitter->FitPeak(input_name, peak_name);

// New with step + tails enabled:
fitter = new FittingUtils(hist, fit_low, fit_high, BackgroundModel::kLinear,
                          true, true, true);
FitResult result = fitter->FitPeak(input_name, peak_name);
```

---

## 4. Verification

- [ ] Project compiles cleanly
- [ ] Basic Gaussian+linear fit produces same results as old `FitPeakStandard`
- [ ] Full component fit produces same results as old `FitPeakDetailed`
- [ ] Exponential background fits converge on test data
- [ ] Multi-peak fits work correctly
- [ ] Plot output looks correct for all configurations
