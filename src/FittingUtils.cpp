#include "FittingUtils.hpp"

Double_t FittingUtils::Gaussian(Double_t x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x - mu) / sigma;
  Double_t gaus_amplitude = par[2];
  return gaus_amplitude * TMath::Exp(-0.5 * z * z);
}

Double_t FittingUtils::LinearBackground(Double_t x, Double_t *par) {
  Double_t bkg_const = par[0];
  Double_t bkg_slope = par[1];
  return bkg_slope * x + bkg_const;
}

Double_t FittingUtils::LowTail(Double_t x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x - mu) / sigma;
  Double_t low_tail_amplitude = par[2];
  Double_t low_tail_range = par[3];
  return low_tail_amplitude * TMath::Exp(low_tail_range * z) /
         TMath::Power(1 + TMath::Exp(z), 4);
}

Double_t FittingUtils::HighTail(Double_t x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x - mu) / sigma;
  Double_t high_tail_amplitude = par[2];
  Double_t high_tail_range = par[3];
  return high_tail_amplitude * TMath::Exp(high_tail_range * z) /
         TMath::Power(1 + TMath::Exp(-z), 4);
}

Double_t FittingUtils::Step(Double_t x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x - mu) / sigma;
  Double_t step_amplitude = par[2];
  return step_amplitude / TMath::Power(1 + TMath::Exp(z), 2);
}

Double_t FittingUtils::Standard(Double_t x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t gaus_amplitude = par[2];
  Double_t bkg_const = par[3];
  Double_t bkg_slope = par[4];

  Double_t gaus_par[3] = {mu, sigma, gaus_amplitude};
  Double_t bkg_par[2] = {bkg_const, bkg_slope};
  return Gaussian(x, gaus_par) + LinearBackground(x, bkg_par);
}

Double_t FittingUtils::Detailed(Double_t x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t gaus_amplitude = par[2];
  Double_t bkg_const = par[3];
  Double_t bkg_slope = par[4];
  Double_t step_amplitude = par[5];
  Double_t low_tail_amplitude = par[6];
  Double_t low_tail_range = par[7];
  Double_t high_tail_amplitude = par[8];
  Double_t high_tail_range = par[9];

  Double_t gaus_par[3] = {mu, sigma, gaus_amplitude};
  Double_t bkg_par[2] = {bkg_const, bkg_slope};
  Double_t step_par[3] = {mu, sigma, step_amplitude};
  Double_t low_tail_par[4] = {mu, sigma, low_tail_amplitude, low_tail_range};
  Double_t high_tail_par[4] = {mu, sigma, high_tail_amplitude, high_tail_range};

  return Gaussian(x, gaus_par) + LinearBackground(x, bkg_par) +
         Step(x, step_par) + LowTail(x, low_tail_par) +
         HighTail(x, high_tail_par);
}

FittingUtils::FittingUtils(TH1 *working_hist, Bool_t isDetailed)
    : fit_range_low_(0), fit_range_high_(1e6) {
  isDetailed_ = isDetailed;
  working_hist_ = working_hist;
  if (!isDetailed_)
    fit_function_ =
        new TF1("Standard", Standard, fit_range_low_, fit_range_high_, 5);
  else
    fit_function_ =
        new TF1("Detailed", Detailed, fit_range_low_, fit_range_high_, 10);
}

FittingUtils::~FittingUtils() {
  fit_function_ = nullptr;
  working_hist_ = nullptr;
}

void FittingUtils::PlotFitStandard(TCanvas *canvas, Int_t color,
                                   const TString peak_name) {
  canvas->Clear();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad1->SetBottomMargin(0.02);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  pad1->SetTopMargin(0.12);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Float_t min_hist_value = 0.9 * fit_range_low_;
  Float_t max_hist_value = 1.1 * fit_range_high_;

  working_hist_->GetXaxis()->SetRangeUser(min_hist_value, max_hist_value);
  working_hist_->Draw();

  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);
  TF1 *peak = new TF1("gaussian", Gaussian, fit_range_low_, fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");

  TF1 *background = new TF1("background", LinearBackground, fit_range_low_,
                            fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(3));
  background->SetParameter(1, fit_function_->GetParameter(4));
  background->SetLineColor(kGreen);
  background->Draw("same");
  pad2->cd();

  Int_t nbins = working_hist_->GetNbinsX();
  TGraph *residuals = new TGraph();
  Int_t point_counter = 0;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t x = working_hist_->GetBinCenter(i);
    if (x < fit_range_low_ || x > fit_range_high_)
      continue;
    Double_t data = working_hist_->GetBinContent(i);
    Double_t fit_val = fit_function_->Eval(x);
    Double_t error = working_hist_->GetBinError(i);

    if (error > 0 && data > 0) {
      Double_t pull = (data - fit_val) / error;
      residuals->SetPoint(point_counter, x, pull);
      point_counter++;
    }
  }

  residuals->SetMarkerStyle(20);
  residuals->SetMarkerSize(0.8);
  residuals->SetMarkerColor(color);
  residuals->SetLineColor(color);
  residuals->SetTitle("");
  residuals->GetXaxis()->SetLimits(min_hist_value, max_hist_value);
  residuals->GetYaxis()->SetTitle("#delta/#sigma");
  residuals->GetXaxis()->SetTitle(working_hist_->GetXaxis()->GetTitle());
  residuals->GetXaxis()->SetTitleSize(0.13);
  residuals->GetYaxis()->SetTitleSize(0.13);
  residuals->GetXaxis()->SetLabelSize(0.12);
  residuals->GetYaxis()->SetLabelSize(0.12);
  residuals->GetXaxis()->SetTitleOffset(1.0);
  residuals->GetYaxis()->SetTitleOffset(0.3);
  residuals->GetYaxis()->SetNdivisions(505);
  residuals->GetXaxis()->SetNdivisions(510);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", min_hist_value, max_hist_value);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  canvas->cd();
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
}

void FittingUtils::PlotFitDetailed(TCanvas *canvas, Int_t color,
                                   const TString peak_name) {
  canvas->Clear();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad1->SetBottomMargin(0.02);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  pad1->SetTopMargin(0.12);
  pad2->SetTopMargin(0.02);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Float_t min_hist_value = 0.9 * fit_range_low_;
  Float_t max_hist_value = 1.1 * fit_range_high_;

  working_hist_->GetXaxis()->SetRangeUser(min_hist_value, max_hist_value);
  working_hist_->Draw();

  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);

  TF1 *peak = new TF1("gaussian", Gaussian, fit_range_low_, fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");

  TF1 *background = new TF1("background", LinearBackground, fit_range_low_,
                            fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(3));
  background->SetParameter(1, fit_function_->GetParameter(4));
  background->SetLineColor(kGreen);
  background->Draw("same");

  TF1 *step = new TF1("step", Step, fit_range_low_, fit_range_high_, 3);
  step->SetParameter(0, fit_function_->GetParameter(0));
  step->SetParameter(1, fit_function_->GetParameter(1));
  step->SetParameter(2, fit_function_->GetParameter(5));
  step->SetLineColor(kMagenta);
  step->Draw("same");

  TF1 *low_tail =
      new TF1("lowtail", LowTail, fit_range_low_, fit_range_high_, 4);
  low_tail->SetParameter(0, fit_function_->GetParameter(0));
  low_tail->SetParameter(1, fit_function_->GetParameter(1));
  low_tail->SetParameter(2, fit_function_->GetParameter(6));
  low_tail->SetParameter(3, fit_function_->GetParameter(7));
  low_tail->SetLineColor(kRed);
  low_tail->Draw("same");

  TF1 *high_tail =
      new TF1("hightail", HighTail, fit_range_low_, fit_range_high_, 4);
  high_tail->SetParameter(0, fit_function_->GetParameter(0));
  high_tail->SetParameter(1, fit_function_->GetParameter(1));
  high_tail->SetParameter(2, fit_function_->GetParameter(8));
  high_tail->SetParameter(3, fit_function_->GetParameter(9));
  high_tail->SetLineColor(kOrange);
  high_tail->Draw("same");

  pad2->cd();

  Int_t nbins = working_hist_->GetNbinsX();
  TGraph *residuals = new TGraph();
  Int_t point_counter = 0;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t x = working_hist_->GetBinCenter(i);
    if (x < fit_range_low_ || x > fit_range_high_)
      continue;
    Double_t data = working_hist_->GetBinContent(i);
    Double_t fit_val = fit_function_->Eval(x);
    Double_t error = working_hist_->GetBinError(i);

    if (error > 0 && data > 0) {
      Double_t pull = (data - fit_val) / error;
      residuals->SetPoint(point_counter, x, pull);
      point_counter++;
    }
  }

  residuals->SetMarkerStyle(20);
  residuals->SetMarkerSize(0.8);
  residuals->SetMarkerColor(color);
  residuals->SetLineColor(color);
  residuals->SetTitle("");
  residuals->GetXaxis()->SetLimits(min_hist_value, max_hist_value);
  residuals->GetYaxis()->SetTitle("#delta/#sigma");
  residuals->GetXaxis()->SetTitle(working_hist_->GetXaxis()->GetTitle());
  residuals->GetXaxis()->SetTitleSize(0.13);
  residuals->GetYaxis()->SetTitleSize(0.13);
  residuals->GetXaxis()->SetLabelSize(0.12);
  residuals->GetYaxis()->SetLabelSize(0.12);
  residuals->GetXaxis()->SetTitleOffset(1.0);
  residuals->GetYaxis()->SetTitleOffset(0.3);
  residuals->GetYaxis()->SetNdivisions(505);
  residuals->GetXaxis()->SetNdivisions(510);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", min_hist_value, max_hist_value);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  canvas->cd();
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
}

FitResultStandard FittingUtils::FitPeakStandard(TCanvas *canvas, Int_t color,
                                                const TString peak_name) {
  FitResultStandard results;

  fit_function_->SetParName(0, "Mu");
  fit_function_->SetParName(1, "Sigma");
  fit_function_->SetParName(2, "GausAmplitude");
  fit_function_->SetParName(3, "BkgConst");
  fit_function_->SetParName(4, "BkgSlope");

  Float_t min_hist_value = 0.9 * fit_range_low_;
  Float_t max_hist_value = 1.1 * fit_range_high_;

  fit_function_->SetParLimits(0, min_hist_value, max_hist_value);
  fit_function_->SetParLimits(1, 0, 0.1 * max_hist_value);
  fit_function_->SetParLimits(2, 0, 1e6);
  fit_function_->SetParLimits(3, 0, 1e6);

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSRN+");
  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFitStandard(canvas, color, peak_name);
    results.mu = fit_function_->GetParameter(0);
    results.mu_error = fit_function_->GetParError(0);
    results.sigma = fit_function_->GetParameter(1);
    results.sigma_error = fit_function_->GetParError(1);
    results.gaus_amplitude = fit_function_->GetParameter(2);
    results.gaus_amplitude_error = fit_function_->GetParError(2);
    results.bkg_const = fit_function_->GetParameter(3);
    results.bkg_const_error = fit_function_->GetParError(3);
    results.bkg_slope = fit_function_->GetParameter(4);
    results.bkg_slope_error = fit_function_->GetParError(4);
  }
  return results;
}

FitResultDetailed FittingUtils::FitPeakDetailed(TCanvas *canvas, Int_t color,
                                                const TString peak_name) {
  FitResultDetailed results;

  fit_function_->SetParName(0, "Mu");
  fit_function_->SetParName(1, "Sigma");
  fit_function_->SetParName(2, "GausAmplitude");
  fit_function_->SetParName(3, "BkgConst");
  fit_function_->SetParName(4, "BkgSlope");
  fit_function_->SetParName(5, "StepAmplitude");
  fit_function_->SetParName(6, "LowTailAmplitude");
  fit_function_->SetParName(7, "LowTailRange");
  fit_function_->SetParName(8, "HighTailAmplitude");
  fit_function_->SetParName(9, "HighTailRange");

  Float_t min_hist_value = 0.9 * fit_range_low_;
  Float_t max_hist_value = 1.1 * fit_range_high_;

  fit_function_->SetParLimits(0, min_hist_value, max_hist_value);
  fit_function_->SetParLimits(1, 0, 0.1 * max_hist_value);
  fit_function_->SetParLimits(2, 0, 1e6);
  fit_function_->SetParLimits(3, 0, 1e6);
  fit_function_->SetParLimits(5, 0, 1e6);
  fit_function_->SetParLimits(6, 0, 1e6);
  fit_function_->SetParLimits(7, -10, 10);
  fit_function_->SetParLimits(8, 0, 1e6);
  fit_function_->SetParLimits(9, -10, 10);

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSRN+");
  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFitDetailed(canvas, color, peak_name);
    results.mu = fit_function_->GetParameter(0);
    results.mu_error = fit_function_->GetParError(0);
    results.sigma = fit_function_->GetParameter(1);
    results.sigma_error = fit_function_->GetParError(1);
    results.gaus_amplitude = fit_function_->GetParameter(2);
    results.gaus_amplitude_error = fit_function_->GetParError(2);
    results.bkg_const = fit_function_->GetParameter(3);
    results.bkg_const_error = fit_function_->GetParError(3);
    results.bkg_slope = fit_function_->GetParameter(4);
    results.bkg_slope_error = fit_function_->GetParError(4);
    results.step_amplitude = fit_function_->GetParameter(5);
    results.step_amplitude_error = fit_function_->GetParError(5);
    results.low_tail_amplitude = fit_function_->GetParameter(6);
    results.low_tail_amplitude_error = fit_function_->GetParError(6);
    results.low_tail_range = fit_function_->GetParameter(7);
    results.low_tail_range_error = fit_function_->GetParError(7);
    results.high_tail_amplitude = fit_function_->GetParameter(8);
    results.high_tail_amplitude_error = fit_function_->GetParError(8);
    results.high_tail_range = fit_function_->GetParameter(9);
    results.high_tail_range_error = fit_function_->GetParError(9);
  }
  return results;
}

void FittingUtils::RegisterCustomFunctions() {
  TF1 *f_standard = new TF1("Standard", Standard, 0, 1000, 5);
  f_standard->SetParNames("Mu", "Sigma", "GausAmp", "BkgConst", "BkgSlope");
  gROOT->GetListOfFunctions()->Add(f_standard);

  TF1 *f_detailed = new TF1("Detailed", Detailed, 0, 1000, 10);
  f_detailed->SetParNames("Mu", "Sigma", "GausAmp", "BkgConst", "BkgSlope",
                          "StepAmp", "LowTailAmp", "LowTailRange",
                          "HighTailAmp", "HighTailRange");
  gROOT->GetListOfFunctions()->Add(f_detailed);

  std::cout << "Custom fitting functions registered and available in FitPanel!"
            << std::endl;
}
