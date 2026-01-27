#ifndef FITTINGUTILS_H
#define FITTINGUTILS_H

#include "PlottingUtils.hpp"
#include <RooMinimizer.h>
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
  Bool_t isDetailed_;
  TF1 *fit_function_;
  TH1 *working_hist_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;

  void PlotFitStandard(TCanvas *canvas, Int_t color, const TString peak_name);
  void PlotFitDetailed(TCanvas *canvas, Int_t color, const TString peak_name);

public:
  FittingUtils(TH1 *working_hist, Bool_t isDetailed);
  ~FittingUtils();

  void SetMu(Double_t expected_mu) {
    fit_function_->SetParameter(0, expected_mu);
  }

  void SetSigma(Double_t expected_sigma) {
    fit_function_->SetParameter(1, expected_sigma);
  }

  void SetGausAmplitude(Double_t expected_amplitude) {
    fit_function_->SetParameter(2, expected_amplitude);
  }

  void SetBkgConst(Double_t expected_bkg_const) {
    fit_function_->SetParameter(3, expected_bkg_const);
  }

  void SetBkgSlope(Double_t expected_bkg_slope) {
    fit_function_->SetParameter(4, expected_bkg_slope);
  }

  void SetStepAmplitude(Double_t expected_step_amplitude) {
    fit_function_->SetParameter(5, expected_step_amplitude);
  }

  void SetLowTailAmplitude(Double_t expected_low_tail_amplitude) {
    fit_function_->SetParameter(6, expected_low_tail_amplitude);
  }

  void SetLowTailRange(Double_t expected_low_tail_range) {
    fit_function_->SetParameter(7, expected_low_tail_range);
  }

  void SetHighTailAmplitude(Double_t expected_high_tail_amplitude) {
    fit_function_->SetParameter(8, expected_high_tail_amplitude);
  }

  void SetHighTailRange(Double_t expected_high_tail_range) {
    fit_function_->SetParameter(9, expected_high_tail_range);
  }

  void SetFitRange(Double_t fit_range_low, Double_t fit_range_high) {
    fit_function_->SetRange(fit_range_low, fit_range_high);
    fit_range_low_ = fit_range_low;
    fit_range_high_ = fit_range_high;
  }

  void FixMu(Double_t value) { fit_function_->FixParameter(0, value); }

  void ReleaseMu() { fit_function_->ReleaseParameter(0); }

  void FixSigma(Double_t value) { fit_function_->FixParameter(1, value); }

  void ReleaseSigma() { fit_function_->ReleaseParameter(1); }

  void FixGausAmplitude(Double_t value) {
    fit_function_->FixParameter(2, value);
  }

  void ReleaseGausAmplitude() { fit_function_->ReleaseParameter(2); }

  void FixBkgConst(Double_t value) { fit_function_->FixParameter(3, value); }

  void ReleaseBkgConst() { fit_function_->ReleaseParameter(3); }

  void FixBkgSlope(Double_t value) { fit_function_->FixParameter(4, value); }

  void ReleaseBkgSlope() { fit_function_->ReleaseParameter(4); }

  void FixStepAmplitude(Double_t value) {
    if (isDetailed_)
      fit_function_->FixParameter(5, value);
  }

  void ReleaseStepAmplitude() {
    if (isDetailed_)
      fit_function_->ReleaseParameter(5);
  }

  void FixLowTailAmplitude(Double_t value) {
    if (isDetailed_)
      fit_function_->FixParameter(6, value);
  }

  void ReleaseLowTailAmplitude() {
    if (isDetailed_)
      fit_function_->ReleaseParameter(6);
  }

  void FixLowTailRange(Double_t value) {
    if (isDetailed_)
      fit_function_->FixParameter(7, value);
  }

  void ReleaseLowTailRange() {
    if (isDetailed_)
      fit_function_->ReleaseParameter(7);
  }

  void FixHighTailAmplitude(Double_t value) {
    if (isDetailed_)
      fit_function_->FixParameter(8, value);
  }

  void ReleaseHighTailAmplitude() {
    if (isDetailed_)
      fit_function_->ReleaseParameter(8);
  }

  void FixHighTailRange(Double_t value) {
    if (isDetailed_)
      fit_function_->FixParameter(9, value);
  }

  void ReleaseHighTailRange() {
    if (isDetailed_)
      fit_function_->ReleaseParameter(9);
  }

  void ReleaseAllParameters() {
    Int_t npar = isDetailed_ ? 10 : 5;
    for (Int_t i = 0; i < npar; i++) {
      fit_function_->ReleaseParameter(i);
    }
  }

  TF1 *GetFitFunction() { return fit_function_; }

  FitResultStandard FitPeakStandard(TCanvas *canvas, Int_t color,
                                    const TString peak_name);

  FitResultDetailed FitPeakDetailed(TCanvas *canvas, Int_t color,
                                    const TString peak_name);

  static void RegisterCustomFunctions();
};

#endif
