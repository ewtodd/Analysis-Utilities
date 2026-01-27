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

struct FitResultGaussianLinear {
  Float_t mu;
  Float_t mu_error;
  Float_t sigma;
  Float_t sigma_error;
  Float_t bkg_const;
  Float_t bkg_const_error;
  Float_t bkg_slope;
  Float_t bkg_slope_error;
};

struct FitResultGaussianTailStep {
  Float_t mu;
  Float_t mu_error;
  Float_t sigma;
  Float_t sigma_error;
  Float_t tau;
  Float_t tau_error;
  Float_t tail_amplitude;
  Float_t tail_amplitude_error;
  Float_t step_amplitude;
  Float_t step_amplitude_error;
};

class FittingUtils {
private:
  Bool_t isGaussianLinear_;
  Bool_t isGaussianLowTailLowStep_;
  Bool_t isGaussianHighTailHighStep_;

  TF1 *fit_function_;
  TTree *working_tree_;
  TH1 *working_hist_;
  Float_t working_value_;
  Int_t num_hist_bins_;
  Int_t min_hist_value_;
  Int_t max_hist_value_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;

  void PlotFitGaussianLinear(TCanvas *canvas, Int_t color,
                             const TString peak_name);
  void PlotFitGaussianLowTailLowStep(TCanvas *canvas, Int_t color,
                                     const TString peak_name);
  void PlotFitGaussianHighTailHighStep(TCanvas *canvas, Int_t color,
                                       const TString peak_name);

public:
  FittingUtils(Bool_t isGaussianLinear, Bool_t isGaussianLowTailLowStep,
               Bool_t isGaussianHighTailHighStep);
  ~FittingUtils();
  void SetNumHistBins(Int_t num_hist_bins) { num_hist_bins_ = num_hist_bins; }

  void SetMinHistValue(Float_t min_hist_value) {
    min_hist_value_ = min_hist_value;
  }

  void SetMaxHistValue(Float_t max_hist_value) {
    max_hist_value_ = max_hist_value;
  }

  void SetExpectedMu(Double_t expected_mu) {
    fit_function_->SetParameter(1, expected_mu);
  }

  void SetExpectedSigma(Double_t expected_sigma) {
    fit_function_->SetParameter(2, expected_sigma);
  }

  void SetExpectedAmplitude(Double_t expected_amplitude) {
    fit_function_->SetParameter(0, expected_amplitude);
  }

  void SetExpectedBackground(Double_t expected_background) {
    fit_function_->SetParameter(3, expected_background);
  }

  void SetExpectedTail(Double_t expected_tail) {
    fit_function_->SetParameter(3, expected_tail);
  }

  void SetExpectedTailAmplitude(Double_t expected_tail_amplitude) {
    fit_function_->SetParameter(4, expected_tail_amplitude);
  }

  void SetExpectedStepAmplitude(Double_t expected_step_amplitude) {
    fit_function_->SetParameter(5, expected_step_amplitude);
  }

  void SetFitRange(Double_t fit_range_low, Double_t fit_range_high) {
    fit_function_->SetRange(fit_range_low, fit_range_high);
    fit_range_low_ = fit_range_low;
    fit_range_high_ = fit_range_high;
  }

  TF1 *GetFitFunction() { return fit_function_; }

  FitResultGaussianLinear
  FitPeakGaussianLinear(TCanvas *canvas, Int_t color, const TString peak_name,
                        const TString formatted_branch_name_with_units);

  FitResultGaussianTailStep
  FitPeakGaussianLowTailLowStep(TCanvas *canvas, Int_t color,
                                const TString peak_name,
                                const TString formatted_branch_name_with_units);

  FitResultGaussianTailStep FitPeakGaussianHighTailHighStep(
      TCanvas *canvas, Int_t color, const TString peak_name,
      const TString formatted_branch_name_with_units);

  Bool_t LoadProcessed(const TString input_name, const TString branch_name,
                       const TString tree_name = "features");

  static void RegisterCustomFunctions();
};

#endif
