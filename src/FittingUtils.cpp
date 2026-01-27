#include "FittingUtils.hpp"

Double_t Gaus(Double_t *x, Double_t *par) {
  Double_t A = par[0];
  Double_t m = par[1];
  Double_t s = par[2];
  Double_t y = x[0] - m;
  return A * (1. / s) * (1. / TMath::Sqrt(2 * TMath::Pi())) *
         TMath::Exp(-0.5 * y * y / (s * s));
}

Double_t LowTail(Double_t *x, Double_t *par) {
  Double_t A = par[0];
  Double_t m = par[1];
  Double_t s = par[2];
  Double_t t = par[3];
  Double_t y = x[0] - m;
  return A * TMath::Exp(y / t) * (TMath::Erfc(y / (TMath::Sqrt(2) * s)) / 2);
}

Double_t HighTail(Double_t *x, Double_t *par) {
  Double_t A = par[0];
  Double_t m = par[1];
  Double_t s = par[2];
  Double_t t = par[3];
  Double_t y = -(x[0] - m);
  return A * TMath::Exp(y / t) * (TMath::Erfc(y / (TMath::Sqrt(2) * s)) / 2);
}

Double_t LowStep(Double_t *x, Double_t *par) {
  Double_t A = par[0];
  Double_t m = par[1];
  Double_t s = par[2];
  Double_t y = x[0] - m;
  return A * (TMath::Erfc(y / (TMath::Sqrt(2) * s)) / 2);
}

Double_t HighStep(Double_t *x, Double_t *par) {
  Double_t A = par[0];
  Double_t m = par[1];
  Double_t s = par[2];
  Double_t y = -(x[0] - m);
  return A * (TMath::Erfc(y / (TMath::Sqrt(2) * s)) / 2);
}

Double_t GaussianLowTailLowStep(Double_t *x, Double_t *par) {
  Double_t gaus_par[3] = {par[0], par[1], par[2]};
  Double_t tail_par[4] = {par[4], par[1], par[2], par[3]};
  Double_t step_par[3] = {par[5], par[1], par[2]};
  return Gaus(x, gaus_par) + LowTail(x, tail_par) + LowStep(x, step_par);
}

Double_t GaussianHighTailHighStep(Double_t *x, Double_t *par) {
  Double_t gaus_par[3] = {par[0], par[1], par[2]};
  Double_t tail_par[4] = {par[4], par[1], par[2], par[3]};
  Double_t step_par[3] = {par[5], par[1], par[2]};
  return Gaus(x, gaus_par) + HighTail(x, tail_par) + HighStep(x, step_par);
}

FittingUtils::FittingUtils(Bool_t isGaussianLinear,
                           Bool_t isGaussianLowTailLowStep,
                           Bool_t isGaussianHighTailHighStep)
    : working_tree_(nullptr), working_value_(0), num_hist_bins_(0),
      min_hist_value_(0), max_hist_value_(0), fit_range_low_(0),
      fit_range_high_(1e6) {
  isGaussianLinear_ = isGaussianLinear;
  isGaussianLowTailLowStep_ = isGaussianLowTailLowStep;
  isGaussianHighTailHighStep_ = isGaussianHighTailHighStep;
  if (isGaussianLinear_)
    fit_function_ = new TF1("gaussian_plus_linear", "gaus(0)+pol1(3)",
                            fit_range_low_, fit_range_high_);
  if (isGaussianLowTailLowStep_)
    fit_function_ = new TF1("gaussian_lowtail_lowstep", GaussianLowTailLowStep,
                            fit_range_low_, fit_range_high_, 6);
  if (isGaussianHighTailHighStep_)
    fit_function_ =
        new TF1("gaussian_hightail_highstep", GaussianHighTailHighStep,
                fit_range_low_, fit_range_high_, 6);
}

FittingUtils::~FittingUtils() {
  fit_function_ = nullptr;
  working_hist_ = nullptr;
  working_tree_ = nullptr;
}

Bool_t FittingUtils::LoadProcessed(const TString input_name,
                                   const TString branch_name,
                                   const TString tree_name) {
  if (gSystem->AccessPathName("root_files")) {
    gSystem->mkdir("root_files", kTRUE);
    std::cout << "Subdirectory root_files not found, creating..." << std::endl;
    std::cout << "Place input files in root_files subdirectory. If you process "
                 "with WaveformProcessingUtils, this is done automatically."
              << std::endl;
    return kFALSE;
  }
  TString input_filename = "root_files/" + input_name + ".root";
  TFile *input_file = new TFile(input_filename, "READ");
  if (!input_file || input_file->IsZombie()) {
    std::cout << "Error: Could not open input file " << input_filename
              << std::endl;
    return kFALSE;
  }
  working_tree_ = static_cast<TTree *>(input_file->Get(tree_name));
  working_tree_->SetBranchAddress(branch_name, &working_value_);
  return kTRUE;
}

void FittingUtils::PlotFitGaussianLinear(TCanvas *canvas, Int_t color,
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
  working_hist_->GetXaxis()->SetRangeUser(min_hist_value_, max_hist_value_);
  PlottingUtils::ConfigureAndDrawHistogram(working_hist_, color);
  working_hist_->GetXaxis()->SetLabelSize(0);
  working_hist_->GetXaxis()->SetTitleSize(0);
  working_hist_->GetYaxis()->SetTitleSize(0.06);
  working_hist_->GetYaxis()->SetLabelSize(0.06);
  working_hist_->GetYaxis()->SetTitleOffset(1.2);
  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);
  TF1 *peak = new TF1("gaussian", "gaus", fit_range_low_, fit_range_high_);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");
  TF1 *background =
      new TF1("background", "pol1", fit_range_low_, fit_range_high_);
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
  residuals->GetXaxis()->SetLimits(min_hist_value_, max_hist_value_);
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
  TF1 *zero_line = new TF1("zero_line", "0", min_hist_value_, max_hist_value_);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");
  canvas->cd();
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
}

void FittingUtils::PlotFitGaussianLowTailLowStep(TCanvas *canvas, Int_t color,
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
  working_hist_->GetXaxis()->SetRangeUser(min_hist_value_, max_hist_value_);
  PlottingUtils::ConfigureAndDrawHistogram(working_hist_, color);
  working_hist_->GetXaxis()->SetLabelSize(0);
  working_hist_->GetXaxis()->SetTitleSize(0);
  working_hist_->GetYaxis()->SetTitleSize(0.06);
  working_hist_->GetYaxis()->SetLabelSize(0.06);
  working_hist_->GetYaxis()->SetTitleOffset(1.2);
  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);
  TF1 *peak = new TF1("gaussian", Gaus, fit_range_low_, fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");
  TF1 *tail = new TF1("lowtail", LowTail, fit_range_low_, fit_range_high_, 4);
  tail->SetParameter(0, fit_function_->GetParameter(4));
  tail->SetParameter(1, fit_function_->GetParameter(1));
  tail->SetParameter(2, fit_function_->GetParameter(2));
  tail->SetParameter(3, fit_function_->GetParameter(3));
  tail->SetLineColor(kRed);
  tail->Draw("same");
  TF1 *step = new TF1("lowstep", LowStep, fit_range_low_, fit_range_high_, 3);
  step->SetParameter(0, fit_function_->GetParameter(5));
  step->SetParameter(1, fit_function_->GetParameter(1));
  step->SetParameter(2, fit_function_->GetParameter(2));
  step->SetLineColor(kGreen);
  step->Draw("same");
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
  residuals->GetXaxis()->SetLimits(min_hist_value_, max_hist_value_);
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
  TF1 *zero_line = new TF1("zero_line", "0", min_hist_value_, max_hist_value_);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");
  canvas->cd();
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
}

void FittingUtils::PlotFitGaussianHighTailHighStep(TCanvas *canvas, Int_t color,
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
  working_hist_->GetXaxis()->SetRangeUser(min_hist_value_, max_hist_value_);
  PlottingUtils::ConfigureAndDrawHistogram(working_hist_, color);
  working_hist_->GetXaxis()->SetLabelSize(0);
  working_hist_->GetXaxis()->SetTitleSize(0);
  working_hist_->GetYaxis()->SetTitleSize(0.06);
  working_hist_->GetYaxis()->SetLabelSize(0.06);
  working_hist_->GetYaxis()->SetTitleOffset(1.2);
  pad1->SetTickx(0);
  fit_function_->Draw("same");
  fit_function_->SetLineColor(kAzure);
  TF1 *peak = new TF1("gaussian", Gaus, fit_range_low_, fit_range_high_, 3);
  peak->SetParameter(0, fit_function_->GetParameter(0));
  peak->SetParameter(1, fit_function_->GetParameter(1));
  peak->SetParameter(2, fit_function_->GetParameter(2));
  peak->SetLineColor(kBlack);
  peak->Draw("same");
  TF1 *tail = new TF1("hightail", HighTail, fit_range_low_, fit_range_high_, 4);
  tail->SetParameter(0, fit_function_->GetParameter(4));
  tail->SetParameter(1, fit_function_->GetParameter(1));
  tail->SetParameter(2, fit_function_->GetParameter(2));
  tail->SetParameter(3, fit_function_->GetParameter(3));
  tail->SetLineColor(kRed);
  tail->Draw("same");
  TF1 *step = new TF1("highstep", HighStep, fit_range_low_, fit_range_high_, 3);
  step->SetParameter(0, fit_function_->GetParameter(5));
  step->SetParameter(1, fit_function_->GetParameter(1));
  step->SetParameter(2, fit_function_->GetParameter(2));
  step->SetLineColor(kGreen);
  step->Draw("same");
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
  residuals->GetXaxis()->SetLimits(min_hist_value_, max_hist_value_);
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
  TF1 *zero_line = new TF1("zero_line", "0", min_hist_value_, max_hist_value_);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");
  canvas->cd();
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
}

FitResultGaussianLinear FittingUtils::FitPeakGaussianLinear(
    TCanvas *canvas, Int_t color, const TString peak_name,
    const TString formatted_branch_name_with_units) {
  FitResultGaussianLinear results;
  working_hist_ =
      new TH1F("", Form(";%s; Counts", formatted_branch_name_with_units.Data()),
               num_hist_bins_, min_hist_value_, max_hist_value_);
  Int_t num_entries = working_tree_->GetEntries();
  for (Int_t i = 0; i < num_entries; i++) {
    working_tree_->GetEntry(i);
    working_hist_->Fill(working_value_);
  }

  fit_function_->SetParName(0, "Amplitude");
  fit_function_->SetParName(1, "Mu");
  fit_function_->SetParName(2, "Sigma");
  fit_function_->SetParName(3, "BkgConst");
  fit_function_->SetParName(4, "BkgSlope");
  fit_function_->SetParLimits(0, 0, 1e6);
  fit_function_->SetParLimits(1, min_hist_value_, max_hist_value_);
  fit_function_->SetParLimits(2, 0, 0.1 * max_hist_value_);
  fit_function_->SetParLimits(3, 0, 1e6);
  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSRN+");
  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFitGaussianLinear(canvas, color, peak_name);
    results.mu = fit_function_->GetParameter("Mu");
    results.mu_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Mu"));
    results.sigma = fit_function_->GetParameter("Sigma");
    results.sigma_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Sigma"));
    results.bkg_const = fit_function_->GetParameter("BkgConst");
    results.bkg_const_error =
        fit_function_->GetParError(fit_function_->GetParNumber("BkgConst"));
    results.bkg_slope = fit_function_->GetParameter("BkgSlope");
    results.bkg_slope_error =
        fit_function_->GetParError(fit_function_->GetParNumber("BkgSlope"));
    return results;
  }
  return results;
}

FitResultGaussianTailStep FittingUtils::FitPeakGaussianLowTailLowStep(
    TCanvas *canvas, Int_t color, const TString peak_name,
    const TString formatted_branch_name_with_units) {
  FitResultGaussianTailStep results;
  working_hist_ =
      new TH1F("", Form(";%s; Counts", formatted_branch_name_with_units.Data()),
               num_hist_bins_, min_hist_value_, max_hist_value_);
  Int_t num_entries = working_tree_->GetEntries();
  for (Int_t i = 0; i < num_entries; i++) {
    working_tree_->GetEntry(i);
    working_hist_->Fill(working_value_);
  }

  fit_function_->SetParName(0, "Amplitude");
  fit_function_->SetParName(1, "Mu");
  fit_function_->SetParName(2, "Sigma");
  fit_function_->SetParName(3, "Tau");
  fit_function_->SetParName(4, "TailAmplitude");
  fit_function_->SetParName(5, "StepAmplitude");

  fit_function_->SetParLimits(1, min_hist_value_, max_hist_value_);
  fit_function_->SetParLimits(2, 0, 0.1 * max_hist_value_);
  fit_function_->SetParLimits(3, 1, max_hist_value_);
  fit_function_->SetParLimits(4, 1, 1e6);
  fit_function_->SetParLimits(5, 10, 1e6);

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSRN+");
  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFitGaussianLowTailLowStep(canvas, color, peak_name);
    results.mu = fit_function_->GetParameter("Mu");
    results.mu_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Mu"));
    results.sigma = fit_function_->GetParameter("Sigma");
    results.sigma_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Sigma"));
    results.tau = fit_function_->GetParameter("Tau");
    results.tau_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Tau"));
    results.tail_amplitude = fit_function_->GetParameter("TailAmplitude");
    results.tail_amplitude_error = fit_function_->GetParError(
        fit_function_->GetParNumber("TailAmplitude"));
    results.step_amplitude = fit_function_->GetParameter("StepAmplitude");
    results.step_amplitude_error = fit_function_->GetParError(
        fit_function_->GetParNumber("StepAmplitude"));
    return results;
  }
  return results;
}

FitResultGaussianTailStep FittingUtils::FitPeakGaussianHighTailHighStep(
    TCanvas *canvas, Int_t color, const TString peak_name,
    const TString formatted_branch_name_with_units) {
  FitResultGaussianTailStep results;
  working_hist_ =
      new TH1F("", Form(";%s; Counts", formatted_branch_name_with_units.Data()),
               num_hist_bins_, min_hist_value_, max_hist_value_);
  Int_t num_entries = working_tree_->GetEntries();
  for (Int_t i = 0; i < num_entries; i++) {
    working_tree_->GetEntry(i);
    if (working_value_ > min_hist_value_ && working_value_ < max_hist_value_)
      working_hist_->Fill(working_value_);
  }

  fit_function_->SetParName(0, "Amplitude");
  fit_function_->SetParName(1, "Mu");
  fit_function_->SetParName(2, "Sigma");
  fit_function_->SetParName(3, "Tau");
  fit_function_->SetParName(4, "TailAmplitude");
  fit_function_->SetParName(5, "StepAmplitude");
  fit_function_->SetParLimits(1, min_hist_value_, max_hist_value_);
  fit_function_->SetParLimits(2, 0, 0.1 * max_hist_value_);
  fit_function_->SetParLimits(3, 0, max_hist_value_);
  fit_function_->SetParLimits(4, 0, 1e6);
  fit_function_->SetParLimits(5, 0, 1e6);

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSRN+");
  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFitGaussianHighTailHighStep(canvas, color, peak_name);
    results.mu = fit_function_->GetParameter("Mu");
    results.mu_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Mu"));
    results.sigma = fit_function_->GetParameter("Sigma");
    results.sigma_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Sigma"));
    results.tau = fit_function_->GetParameter("Tau");
    results.tau_error =
        fit_function_->GetParError(fit_function_->GetParNumber("Tau"));
    results.tail_amplitude = fit_function_->GetParameter("TailAmplitude");
    results.tail_amplitude_error = fit_function_->GetParError(
        fit_function_->GetParNumber("TailAmplitude"));
    results.step_amplitude = fit_function_->GetParameter("StepAmplitude");
    results.step_amplitude_error = fit_function_->GetParError(
        fit_function_->GetParNumber("StepAmplitude"));
    return results;
  }
  return results;
}

void FittingUtils::RegisterCustomFunctions() {
  TF1 *f_gaus = new TF1("CustomGaus", Gaus, 0, 1000, 3);
  f_gaus->SetParNames("Amplitude", "Mean", "Sigma");
  gROOT->GetListOfFunctions()->Add(f_gaus);

  TF1 *f_gltls =
      new TF1("GausLowTailLowStep", GaussianLowTailLowStep, 0, 1000, 6);
  f_gltls->SetParNames("Amplitude", "Mean", "Sigma", "Tau", "TailAmp",
                       "StepAmp");
  gROOT->GetListOfFunctions()->Add(f_gltls);

  TF1 *f_ghths =
      new TF1("GausHighTailHighStep", GaussianHighTailHighStep, 0, 1000, 6);
  f_ghths->SetParNames("Amplitude", "Mean", "Sigma", "Tau", "TailAmp",
                       "StepAmp");
  gROOT->GetListOfFunctions()->Add(f_ghths);

  TF1 *f_lowtail = new TF1("LowTail", LowTail, 0, 1000, 4);
  f_lowtail->SetParNames("Amplitude", "Mean", "Sigma", "Tau");
  gROOT->GetListOfFunctions()->Add(f_lowtail);

  TF1 *f_hightail = new TF1("HighTail", HighTail, 0, 1000, 4);
  f_hightail->SetParNames("Amplitude", "Mean", "Sigma", "Tau");
  gROOT->GetListOfFunctions()->Add(f_hightail);

  TF1 *f_lowstep = new TF1("LowStep", LowStep, 0, 1000, 3);
  f_lowstep->SetParNames("Amplitude", "Mean", "Sigma");
  gROOT->GetListOfFunctions()->Add(f_lowstep);

  TF1 *f_highstep = new TF1("HighStep", HighStep, 0, 1000, 3);
  f_highstep->SetParNames("Amplitude", "Mean", "Sigma");
  gROOT->GetListOfFunctions()->Add(f_highstep);

  std::cout << "Custom fitting functions registered and available in FitPanel!"
            << std::endl;
}
