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

namespace FittingFunctions {
Double_t Gaussian(Double_t *x, Double_t *par);
Double_t LinearBackground(Double_t *x, Double_t *par);
Double_t Step(Double_t *x, Double_t *par);
Double_t LowTail(Double_t *x, Double_t *par);
Double_t HighTail(Double_t *x, Double_t *par);
Double_t Standard(Double_t *x, Double_t *par);
Double_t Detailed(Double_t *x, Double_t *par);
} // namespace FittingFunctions

struct FitResultStandard {
  Float_t gaus_amplitude;
  Float_t gaus_amplitude_error;
  Float_t mu;
  Float_t mu_error;
  Float_t sigma;
  Float_t sigma_error;
  Float_t bkg_const;
  Float_t bkg_const_error;
  Float_t bkg_slope;
  Float_t bkg_slope_error;
};

struct FitResultDetailed {
  Float_t gaus_amplitude;
  Float_t gaus_amplitude_error;
  Float_t mu;
  Float_t mu_error;
  Float_t sigma;
  Float_t sigma_error;
  Float_t bkg_const;
  Float_t bkg_const_error;
  Float_t bkg_slope;
  Float_t bkg_slope_error;
  Float_t step_amplitude;
  Float_t step_amplitude_error;
  Float_t low_tail_amplitude;
  Float_t low_tail_amplitude_error;
  Float_t low_tail_range;
  Float_t low_tail_range_error;
  Float_t high_tail_amplitude;
  Float_t high_tail_amplitude_error;
  Float_t high_tail_range;
  Float_t high_tail_range_error;
};

class FittingUtils {
private:
  TF1 *fit_function_;
  TH1 *working_hist_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;
  Bool_t use_flat_background_;
  Bool_t isDetailed_;
  Bool_t use_step_;
  Bool_t use_low_tail_;
  Bool_t use_high_tail_;
  Bool_t use_manual_init_;
  std::vector<Double_t> manual_params_;

  void PlotFitStandard(const TString input_name, const TString peak_name);
  void PlotFitDetailed(const TString input_name, const TString peak_name);
  Double_t EstimateBackground();
  Double_t ClampToBounds(Int_t param_index, Double_t value);

public:
  FittingUtils(TH1 *working_hist, Float_t fit_range_low, Float_t fit_range_high,
               Bool_t use_flat_background, Bool_t isDetailed,
               Bool_t use_step = kTRUE, Bool_t use_low_tail = kTRUE,
               Bool_t use_high_tail = kTRUE);
  ~FittingUtils();

  void UseFlatBackground(Bool_t use_flat_background = kTRUE) {
    use_flat_background_ = use_flat_background;
  }
  void UseStep(Bool_t use_step = kTRUE) { use_step_ = use_step; }
  void UseLowTail(Bool_t use_low_tail = kTRUE) { use_low_tail_ = use_low_tail; }
  void UseHighTail(Bool_t use_high_tail = kTRUE) {
    use_high_tail_ = use_high_tail;
  }

  void SetManualParameters(const std::vector<Double_t> &params);
  void SetManualParameter(Int_t index, Double_t value);
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  TF1 *GetFitFunction() { return fit_function_; }

  FitResultStandard FitPeakStandard(const TString input_name,
                                    const TString peak_name);
  FitResultDetailed FitPeakDetailed(const TString input_name,
                                    const TString peak_name);

  static void RegisterCustomFunctions();
};

#endif
