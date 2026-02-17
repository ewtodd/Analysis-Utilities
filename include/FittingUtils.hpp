#ifndef FITTINGUTILS_H
#define FITTINGUTILS_H

#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TMath.h>
#include <TPad.h>
#include <TSystem.h>
#include <TTree.h>

enum class BackgroundModel { kFLAT, kLINEAR };

namespace FittingFunctions {
Double_t Gaussian(Double_t *x, Double_t *par);
Double_t LinearBackground(Double_t *x, Double_t *par);
Double_t Step(Double_t *x, Double_t *par);
Double_t LowTail(Double_t *x, Double_t *par);
Double_t HighTail(Double_t *x, Double_t *par);
Double_t PeakFunction(Double_t *x, Double_t *par);
Double_t DoublePeakFunction(Double_t *x, Double_t *par);
Double_t TriplePeakFunction(Double_t *x, Double_t *par);
} // namespace FittingFunctions

struct PeakFitResult {
  Float_t mu, mu_error;
  Float_t sigma, sigma_error;
  Float_t gaus_amplitude, gaus_amplitude_error;
  Float_t step_amplitude, step_amplitude_error;
  Float_t low_exp_tail_amplitude, low_exp_tail_amplitude_error;
  Float_t low_exp_tail_decay, low_exp_tail_decay_error;
  Float_t low_lin_tail_amplitude, low_lin_tail_amplitude_error;
  Float_t low_lin_tail_slope, low_lin_tail_slope_error;
  Float_t high_exp_tail_amplitude, high_exp_tail_amplitude_error;
  Float_t high_exp_tail_decay, high_exp_tail_decay_error;
};

struct FitResult {
  std::vector<PeakFitResult> peaks; // 1, 2, or 3 entries
  Float_t bkg_constant, bkg_constant_error;
  Float_t lin_bkg_slope, lin_bkg_slope_error;
  Float_t reduced_chi2;
  Bool_t valid;
};

class FittingUtils {
private:
  TF1 *fit_function_;
  TH1 *working_hist_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;

  BackgroundModel bkg_model_;
  Bool_t use_step_;
  Bool_t use_low_exp_tail_;
  Bool_t use_low_lin_tail_;
  Bool_t use_high_exp_tail_;

  Bool_t use_manual_init_;
  std::vector<Double_t> manual_params_;

  void PlotFit(const TString input_name, const TString peak_name);
  void PlotFitDoublePeak(const TString input_name, const TString peak_name);
  void PlotFitTriplePeak(const TString input_name, const TString peak_name);

  Double_t EstimateBackground();
  Double_t ClampToBounds(Int_t param_index, Double_t value);

  void SwapDoublePeakParameters();

public:
  FittingUtils(TH1 *working_hist, Float_t fit_range_low, Float_t fit_range_high,
               BackgroundModel bkg_model = BackgroundModel::kLINEAR,
               Bool_t use_step = kFALSE, Bool_t use_low_exp_tail = kFALSE,
               Bool_t use_low_lin_tail = kFALSE,
               Bool_t use_high_exp_tail = kFALSE);
  ~FittingUtils();

  void SetBackgroundModel(BackgroundModel bkg_model) { bkg_model_ = bkg_model; }
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

  void SetManualParameters(const std::vector<Double_t> &params);
  void SetManualParameter(Int_t index, Double_t value);
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  TF1 *GetFitFunction() { return fit_function_; }

  FitResult FitPeak(const TString input_name, const TString peak_name);
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
