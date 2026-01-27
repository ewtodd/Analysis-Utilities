#include "FittingUtils.hpp"

FittingUtils::FittingUtils(TH1 *working_hist, Bool_t isDetailed)
    : fit_range_low_(0), fit_range_high_(1e6) {
  isDetailed_ = isDetailed;
  working_hist_ = static_cast<TH1F *>(working_hist->Clone());

  ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(Int_t(1e10));
  ROOT::Math::MinimizerOptions::SetDefaultMaxIterations(Int_t(1e6));
  ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-6);
  ROOT::Math::MinimizerOptions::SetDefaultStrategy(2);

  if (!isDetailed_) {
    fit_function_ = new TF1("Standard", &FittingFunctions::Standard,
                            fit_range_low_, fit_range_high_, 5);

    fit_function_->SetParName(0, "Mu");
    fit_function_->SetParName(1, "Sigma");
    fit_function_->SetParName(2, "GausAmplitude");
    fit_function_->SetParName(3, "BkgConst");
    fit_function_->SetParName(4, "BkgSlope");

    fit_function_->SetParLimits(0, 0, 1e7);
    fit_function_->SetParLimits(1, 0, 1e7);
    fit_function_->SetParLimits(2, 0, 1e7);
    fit_function_->SetParLimits(3, 0, 1e6);
    fit_function_->SetParLimits(4, -1e3, 1e3);

  } else {
    fit_function_ = new TF1("Detailed", &FittingFunctions::Detailed,
                            fit_range_low_, fit_range_high_, 10);
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

    fit_function_->SetParLimits(0, 0, 1e7);
    fit_function_->SetParLimits(1, 0, 1e7);
    fit_function_->SetParLimits(2, -1e-6, 1e7);
    fit_function_->SetParLimits(5, -1e-6, 1e7);
    fit_function_->SetParLimits(6, -1e-6, 1e7);
    fit_function_->SetParLimits(7, 0, 1e7);
    fit_function_->SetParLimits(8, -1e-6, 1e7);
    fit_function_->SetParLimits(9, 0, 1e7);
  }
}

FittingUtils::~FittingUtils() {
  fit_function_ = nullptr;
  working_hist_ = nullptr;
}

Double_t FittingFunctions::Gaussian(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x[0] - mu) / sigma;
  Double_t gaus_amplitude = par[2];
  return gaus_amplitude * TMath::Exp(-0.5 * z * z);
}

Double_t FittingFunctions::LinearBackground(Double_t *x, Double_t *par) {
  Double_t bkg_const = par[0];
  Double_t bkg_slope = par[1];
  return bkg_slope * x[0] + bkg_const;
}

Double_t FittingFunctions::LowTail(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x[0] - mu) / sigma;
  Double_t low_tail_amplitude = par[2];
  Double_t low_tail_range = par[3];
  return low_tail_amplitude * TMath::Exp(-low_tail_range * z) /
         TMath::Power(1 + TMath::Exp(z), 4);
}

Double_t FittingFunctions::HighTail(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x[0] - mu) / sigma;
  Double_t high_tail_amplitude = par[2];
  Double_t high_tail_range = par[3];
  return high_tail_amplitude * TMath::Exp(high_tail_range * z) /
         TMath::Power(1 + TMath::Exp(-z), 4);
}

Double_t FittingFunctions::Step(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t z = (x[0] - mu) / sigma;
  Double_t step_amplitude = par[2];
  return step_amplitude / TMath::Power(1 + TMath::Exp(z), 2);
}

Double_t FittingFunctions::Standard(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  Double_t gaus_amplitude = par[2];
  Double_t bkg_const = par[3];
  Double_t bkg_slope = par[4];

  Double_t gaus_par[3] = {mu, sigma, gaus_amplitude};
  Double_t bkg_par[2] = {bkg_const, bkg_slope};
  return Gaussian(x, gaus_par) + LinearBackground(x, bkg_par);
}

Double_t FittingFunctions::Detailed(Double_t *x, Double_t *par) {
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

void FittingUtils::PlotFitStandard(TCanvas *canvas, Int_t color,
                                   const TString peak_name) {
  canvas->Clear();
  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad1->SetBottomMargin(0.04);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  pad1->SetTopMargin(0.12);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Float_t min_hist_value = 0.9 * fit_range_low_;
  Float_t max_hist_value = 1.1 * fit_range_high_;

  working_hist_->GetXaxis()->SetRangeUser(min_hist_value, max_hist_value);
  working_hist_->GetXaxis()->SetLabelSize(0);
  working_hist_->GetXaxis()->SetTitleSize(0);
  working_hist_->SetLineColor(kViolet);
  working_hist_->SetLineWidth(2);
  working_hist_->Draw();

  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);

  TF1 *peak = new TF1("gaussian", FittingFunctions::Gaussian, fit_range_low_,
                      fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
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
  residuals->GetYaxis()->CenterTitle(kTRUE);
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
  pad1->SetBottomMargin(0.04);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  pad1->SetTopMargin(0.12);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Float_t min_hist_value = 0.9 * fit_range_low_;
  Float_t max_hist_value = 1.1 * fit_range_high_;

  working_hist_->GetXaxis()->SetRangeUser(min_hist_value, max_hist_value);
  working_hist_->GetXaxis()->SetLabelSize(0);
  working_hist_->GetXaxis()->SetTitleSize(0);
  working_hist_->SetLineColor(kViolet);
  working_hist_->SetLineWidth(2);
  working_hist_->Draw();

  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);

  TF1 *peak = new TF1("gaussian", FittingFunctions::Gaussian, fit_range_low_,
                      fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(3));
  background->SetParameter(1, fit_function_->GetParameter(4));
  background->SetLineColor(kGreen);
  background->Draw("same");

  TF1 *step = new TF1("step", FittingFunctions::Step, fit_range_low_,
                      fit_range_high_, 3);
  step->SetParameter(0, fit_function_->GetParameter(0));
  step->SetParameter(1, fit_function_->GetParameter(1));
  step->SetParameter(2, fit_function_->GetParameter(5));
  step->SetLineColor(kMagenta);
  step->Draw("same");

  TF1 *low_tail = new TF1("lowtail", FittingFunctions::LowTail, fit_range_low_,
                          fit_range_high_, 4);
  low_tail->SetParameter(0, fit_function_->GetParameter(0));
  low_tail->SetParameter(1, fit_function_->GetParameter(1));
  low_tail->SetParameter(2, fit_function_->GetParameter(6));
  low_tail->SetParameter(3, fit_function_->GetParameter(7));
  low_tail->SetLineColor(kRed);
  low_tail->Draw("same");

  TF1 *high_tail = new TF1("hightail", FittingFunctions::HighTail,
                           fit_range_low_, fit_range_high_, 4);
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
  residuals->GetYaxis()->CenterTitle(kTRUE);
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

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMEN+");
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

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMEN+");
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
  TF1 *f_standard =
      new TF1("Standard", &FittingFunctions::Standard, 0, 1000, 5);
  f_standard->SetParNames("Mu", "Sigma", "GausAmp", "BkgConst", "BkgSlope");
  gROOT->GetListOfFunctions()->Add(f_standard);

  TF1 *f_detailed =
      new TF1("Detailed", &FittingFunctions::Detailed, 0, 1000, 10);
  f_detailed->SetParNames("Mu", "Sigma", "GausAmp", "BkgConst", "BkgSlope",
                          "StepAmp", "LowTailAmp", "LowTailRange",
                          "HighTailAmp", "HighTailRange");
  gROOT->GetListOfFunctions()->Add(f_detailed);

  std::cout << "Custom fitting functions registered and available in FitPanel!"
            << std::endl;
}
