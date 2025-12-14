#ifndef FITUTILS_H
#define FITUTILS_H

#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1F.h>
#include <TSystem.h>
#include <TTree.h>

class FitUtils {
private:
  TF1 *fit_function_;
  TTree *working_tree_;
  Float_t working_value_;
  Bool_t fit_successful_;
  Int_t num_hist_bins_;
  Int_t max_hist_value_;

public:
  FitUtils();
  ~FitUtils();
  void SetNumHistBins(Int_t num_hist_bins) { num_hist_bins_ = num_hist_bins; }

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
    fit_function_->SetParameter(4, expected_background);
  }

  void SetFitRange(Double_t fit_range_low, Double_t fit_range_high) {
    fit_function_->SetRange(fit_range_low, fit_range_high);
  }

  Bool_t FitPeak(TCanvas *canvas, const TString input_name,
                 const TString formatted_branch_name);

  Bool_t LoadProcessed(const TString input_name, const TString branch_name);
};

#endif
