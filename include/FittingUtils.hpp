#ifndef FITTINGUTILS_H
#define FITTINGUTILS_H

#include "PlottingUtils.hpp"
#include <TROOT.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTree.h>

// Forward declaration (defined in InteractiveFitEditor.hpp)
Bool_t LaunchInteractiveFitEditor(TH1 *hist, TF1 *fit_func,
                                  Double_t range_low, Double_t range_high,
                                  Int_t num_peaks = 1,
                                  const TString &info_label = "");

struct HyperEMGPeakFitResult {
  Float_t mu = -1, mu_error = -1;
  Float_t sigma = -1, sigma_error = -1;
  Float_t amplitude = -1, amplitude_error = -1;
  Float_t w_low1 = -1, w_low1_error = -1;
  Float_t tau_low1 = -1, tau_low1_error = -1;
  Float_t w_low2 = -1, w_low2_error = -1;
  Float_t tau_low2 = -1, tau_low2_error = -1;
  Float_t w_high1 = -1, w_high1_error = -1;
  Float_t tau_high1 = -1, tau_high1_error = -1;
  Float_t w_high2 = -1, w_high2_error = -1;
  Float_t tau_high2 = -1, tau_high2_error = -1;
};

struct HyperEMGFitResult {
  std::vector<HyperEMGPeakFitResult> peaks;
  Float_t bkg_constant = -1, bkg_constant_error = -1;
  Float_t lin_bkg_slope = -1, lin_bkg_slope_error = -1;
  Float_t reduced_chi2 = -1;
  Bool_t valid = kFALSE;
};

namespace FittingFunctions {
Double_t Gaussian(Double_t *x, Double_t *par);
Double_t LinearBackground(Double_t *x, Double_t *par);
Double_t Step(Double_t *x, Double_t *par);
Double_t LowTail(Double_t *x, Double_t *par);
Double_t HighTail(Double_t *x, Double_t *par);
Double_t PeakFunction(Double_t *x, Double_t *par);
Double_t DoublePeakFunction(Double_t *x, Double_t *par);
Double_t TriplePeakFunction(Double_t *x, Double_t *par);

// Hyper-EMG peak shape: sum of exponentially-modified Gaussians
// 2 low-side + 2 high-side components (each toggleable via weight=0)
// Per-peak params (11):
//   [0] mu  [1] sigma  [2] amplitude
//   [3] w_low1   [4] tau_low1   (fast low-side tail)
//   [5] w_low2   [6] tau_low2   (slow low-side tail, fix w=0 to disable)
//   [7] w_high1  [8] tau_high1  (primary high-side tail)
//   [9] w_high2  [10] tau_high2 (secondary high-side, fix w=0 to disable)
// Weights are normalized inside the function: w_i / sum(w_i)
static const Int_t kHyperEMGParamsPerPeak = 11;
static const Int_t kHyperEMGNLow = 2;
static const Int_t kHyperEMGNHigh = 2;

// Single EMG component evaluators (for component drawing in the editor)
// Returns w/(2*tau) * exp(...) * erfc(...) with overflow protection
Double_t EMGLowComponent(Double_t x, Double_t mu, Double_t sigma,
                         Double_t w, Double_t tau);
Double_t EMGHighComponent(Double_t x, Double_t mu, Double_t sigma,
                          Double_t w, Double_t tau);

Double_t HyperEMGPeak(Double_t *x, Double_t *par);
Double_t DoublePeakHyperEMG(Double_t *x, Double_t *par);
Double_t TriplePeakHyperEMG(Double_t *x, Double_t *par);

void SetHyperEMGParLimits(TF1 *func, Int_t peak_offset);
void SetHyperEMGInitialParams(TF1 *func, Int_t peak_offset,
                              Double_t mu_init, Double_t sigma_init,
                              Double_t amp_init);
void SetHyperEMGParNames(TF1 *func, Int_t peak_offset,
                         const TString &peak_label = "");

// Hyper-EMG fit orchestration (save/load, interactive editor, plotting)
void SaveHyperEMGParams(TF1 *func, Double_t rlo, Double_t rhi,
                        const TString &input, const TString &peak);
Bool_t LoadHyperEMGParams(TF1 *func, Double_t &rlo, Double_t &rhi,
                          const TString &input, const TString &peak);
HyperEMGPeakFitResult ExtractHyperEMGPeak(TF1 *func, Int_t offset);
void SetupHyperEMGBackground(TF1 *func, Int_t bkg_offset,
                             Double_t peak_height, Bool_t use_flat_bkg);

void PlotHyperEMGFit(TH1 *hist, TF1 *func, Double_t range_low,
                     Double_t range_high, Int_t num_peaks,
                     const TString &input_name, const TString &peak_name,
                     const TString &label = "");

HyperEMGFitResult FitSingleHyperEMG(TH1 *hist, Double_t fit_low,
                                     Double_t fit_high, Bool_t use_flat_bkg,
                                     const TString &label,
                                     const TString &peak_name,
                                     Bool_t interactive);
HyperEMGFitResult FitDoubleHyperEMG(TH1 *hist, Double_t fit_low,
                                     Double_t fit_high, Double_t mu1_init,
                                     Double_t mu2_init, Bool_t use_flat_bkg,
                                     const TString &label,
                                     const TString &peak_name,
                                     Bool_t interactive);
HyperEMGFitResult FitTripleHyperEMG(TH1 *hist, Double_t fit_low,
                                     Double_t fit_high, Double_t mu1_init,
                                     Double_t mu2_init, Double_t mu3_init,
                                     Bool_t use_flat_bkg,
                                     const TString &label,
                                     const TString &peak_name,
                                     Bool_t interactive);
} // namespace FittingFunctions

struct PeakFitResult {
  Float_t mu = -1, mu_error = -1;
  Float_t sigma = -1, sigma_error = -1;
  Float_t gaus_amplitude = -1, gaus_amplitude_error = -1;
  Float_t step_amplitude = -1, step_amplitude_error = -1;
  Float_t low_exp_tail_amplitude = -1, low_exp_tail_amplitude_error = -1;
  Float_t low_exp_tail_decay = -1, low_exp_tail_decay_error = -1;
  Float_t low_lin_tail_amplitude = -1, low_lin_tail_amplitude_error = -1;
  Float_t low_lin_tail_slope = -1, low_lin_tail_slope_error = -1;
  Float_t high_exp_tail_amplitude = -1, high_exp_tail_amplitude_error = -1;
  Float_t high_exp_tail_decay = -1, high_exp_tail_decay_error = -1;
};

struct FitResult {
  std::vector<PeakFitResult> peaks; // 1-3 entries supported
  Float_t bkg_constant = -1, bkg_constant_error = -1;
  Float_t lin_bkg_slope = -1, lin_bkg_slope_error = -1;
  Float_t reduced_chi2 = -1;
  Bool_t valid = kFALSE;
};

class FittingUtils {
private:
  TF1 *fit_function_;
  TH1 *working_hist_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;

  Bool_t use_flat_background_;
  Bool_t use_step_;
  Bool_t use_low_exp_tail_;
  Bool_t use_low_lin_tail_;
  Bool_t use_high_exp_tail_;

  Bool_t use_manual_init_;
  Bool_t interactive_;
  std::vector<Double_t> manual_params_;

  Double_t EstimateBackground();
  Double_t ClampToBounds(Int_t param_index, Double_t value);

  void SaveInteractiveParams(const TString &input_name,
                             const TString &peak_name);
  Bool_t LoadInteractiveParams(const TString &input_name,
                               const TString &peak_name);

  void SortPeaksByMu(Int_t num_peaks);
  void PlotResidualHistogram(TGraph *residuals, const TString &input_name,
                             const TString &peak_name);

public:
  FittingUtils(TH1 *working_hist, Float_t fit_range_low, Float_t fit_range_high,
               Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
               Bool_t use_low_exp_tail = kFALSE,
               Bool_t use_low_lin_tail = kFALSE,
               Bool_t use_high_exp_tail = kFALSE);
  ~FittingUtils();

  void SetBackgroundModel(Bool_t use_flat_background) {
    use_flat_background_ = use_flat_background;
  }
  void SetStep(Bool_t use_step = kTRUE) { use_step_ = use_step; }
  void SetLowExpTail(Bool_t use_low_exp_tail = kTRUE) {
    use_low_exp_tail_ = use_low_exp_tail;
  }
  void SetLowLinTail(Bool_t use_low_lin_tail = kTRUE) {
    use_low_lin_tail_ = use_low_lin_tail;
  }
  void SetHighExpTail(Bool_t use_high_exp_tail = kTRUE) {
    use_high_exp_tail_ = use_high_exp_tail;
  }
  void SetInteractive(Bool_t interactive = kTRUE) {
    interactive_ = interactive;
  }

  void SetManualParameters(const std::vector<Double_t> &params);
  void SetManualParameter(Int_t index, Double_t value);
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  TF1 *GetFitFunction() { return fit_function_; }
  void SetFitFunction(TF1 *func) { fit_function_ = func; }

  void PlotFitSinglePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitDoublePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitTriplePeak(const TString input_name, const TString peak_name,
                         const TString label = "");

  FitResult FitSinglePeak(const TString input_name, const TString peak_name);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakFitResult &constrained_peak,
                          Double_t mu2_init);
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks,
                          Double_t mu3_init);
};

#endif
