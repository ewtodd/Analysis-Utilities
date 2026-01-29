#include "FittingUtils.hpp"

FittingUtils::FittingUtils(TH1 *working_hist, Float_t fit_range_low,
                           Float_t fit_range_high, Bool_t isDetailed,
                           Bool_t use_step, Bool_t use_low_tail,
                           Bool_t use_high_tail) {
  working_hist_ = static_cast<TH1F *>(working_hist->Clone());
  fit_range_low_ = fit_range_low;
  fit_range_high_ = fit_range_high;
  isDetailed_ = isDetailed;
  use_step_ = use_step;
  use_low_tail_ = use_low_tail;
  use_high_tail_ = use_high_tail;

  if (!isDetailed_) {
    fit_function_ = new TF1("Standard", &FittingFunctions::Standard,
                            fit_range_low_, fit_range_high_, 5);

    fit_function_->SetParName(0, "Mu");
    fit_function_->SetParName(1, "Sigma");
    fit_function_->SetParName(2, "GausAmplitude");
    fit_function_->SetParName(3, "BkgConst");
    fit_function_->SetParName(4, "BkgSlope");

    Double_t mu_init = (fit_range_low_ + fit_range_high_) / 2;
    Double_t range_width = fit_range_high_ - fit_range_low_;
    Double_t sigma_init = range_width * 0.01;

    Double_t max_bin = working_hist_->GetMaximumBin();
    Double_t peak_height = working_hist_->GetBinContent(max_bin);
    Double_t bkg_estimate = EstimateBackground();

    fit_function_->SetParLimits(0, fit_range_low_, fit_range_high_);
    fit_function_->SetParLimits(1, range_width * 0.001, range_width * 0.5);
    fit_function_->SetParLimits(2, 0, peak_height * 1.5);
    fit_function_->SetParLimits(3, 0, peak_height * 0.5);
    fit_function_->SetParLimits(4, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);

    fit_function_->SetParameter(0, mu_init);
    fit_function_->SetParameter(1, sigma_init);
    fit_function_->SetParameter(2, peak_height * 0.8);
    fit_function_->SetParameter(3, bkg_estimate);
    fit_function_->SetParameter(4, 0);
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
    fit_function_->SetParName(7, "LowTailSlope");
    fit_function_->SetParName(8, "HighTailAmplitude");
    fit_function_->SetParName(9, "HighTailSlope");

    Double_t mu_init = (fit_range_low_ + fit_range_high_) / 2;
    Double_t range_width = fit_range_high_ - fit_range_low_;
    Double_t sigma_init = range_width * 0.01;

    Double_t max_bin = working_hist_->GetMaximumBin();
    Double_t peak_height = working_hist_->GetBinContent(max_bin);
    Double_t bkg_estimate = EstimateBackground();

    fit_function_->SetParLimits(0, fit_range_low_, fit_range_high_);
    fit_function_->SetParLimits(1, range_width * 0.001, range_width * 0.5);
    fit_function_->SetParLimits(2, 0, peak_height * 1.5);
    fit_function_->SetParLimits(3, 0, peak_height * 0.5);
    fit_function_->SetParLimits(4, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);

    fit_function_->SetParameter(0, mu_init);
    fit_function_->SetParameter(1, sigma_init);
    fit_function_->SetParameter(2, peak_height * 0.8);
    fit_function_->SetParameter(3, bkg_estimate);
    fit_function_->SetParameter(4, 0);

    fit_function_->SetParLimits(5, 0, peak_height * 0.3);
    fit_function_->SetParameter(5, 0);

    fit_function_->SetParLimits(6, 0, peak_height);
    fit_function_->SetParLimits(7, 0, 100);
    fit_function_->SetParameter(6, peak_height * 0.10);
    fit_function_->SetParameter(7, 0.2);

    fit_function_->SetParLimits(8, 0, peak_height);
    fit_function_->SetParLimits(9, 0, 100);
    fit_function_->SetParameter(8, peak_height * 0.10);
    fit_function_->SetParameter(9, 0.2);

    std::cout << "Detailed fit configuration:" << std::endl;
    std::cout << "  Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
              << std::endl;
    std::cout << "  Low tail: " << (use_low_tail_ ? "ENABLED" : "DISABLED")
              << std::endl;
    std::cout << "  High tail: " << (use_high_tail_ ? "ENABLED" : "DISABLED")
              << std::endl;
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
  if (sigma <= 0)
    return 0;

  Double_t tail_amplitude = par[2];
  Double_t tail_slope = par[3];

  Double_t dx = (x[0] - mu) / (TMath::Sqrt(2) * sigma);

  Double_t exp_arg = (x[0] - mu) / tail_slope;

  Double_t transition = 0.5 * (1.0 - TMath::Erf(dx));

  return tail_amplitude * TMath::Exp(exp_arg) * transition;
}

Double_t FittingFunctions::HighTail(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  if (sigma <= 0)
    return 0;

  Double_t tail_amplitude = par[2];
  Double_t tail_slope = par[3];

  Double_t dx = (mu - x[0]) / (TMath::Sqrt(2) * sigma);

  Double_t exp_arg = (mu - x[0]) / tail_slope;

  Double_t transition = 0.5 * (1.0 - TMath::Erf(dx));

  return tail_amplitude * TMath::Exp(exp_arg) * transition;
}

Double_t FittingFunctions::Step(Double_t *x, Double_t *par) {
  Double_t mu = par[0];
  Double_t sigma = par[1];
  if (sigma <= 0)
    return 0;

  Double_t z = (x[0] - mu) / sigma;
  Double_t step_amplitude = par[2];

  Double_t denominator = TMath::Power(1 + TMath::Exp(z), 2);
  if (denominator < 1e-100)
    return 0;

  return step_amplitude / denominator;
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

void FittingUtils::PlotFitStandard(const TString peak_name) {
  TCanvas *canvas = new TCanvas(PlottingUtils::GetRandomName(), "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);

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
  residuals->SetMarkerColor(kAzure);
  residuals->SetLineColor(kAzure);
  residuals->SetTitle("");
  Double_t actual_min = working_hist_->GetXaxis()->GetBinLowEdge(
      working_hist_->GetXaxis()->GetFirst());
  Double_t actual_max = working_hist_->GetXaxis()->GetBinUpEdge(
      working_hist_->GetXaxis()->GetLast());
  residuals->GetXaxis()->SetLimits(actual_min, actual_max);
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

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  canvas->cd();
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(canvas, "log_" + peak_name + ".png", kFALSE);
}

void FittingUtils::PlotFitDetailed(const TString peak_name) {
  TCanvas *canvas = new TCanvas(PlottingUtils::GetRandomName(), "", 1200, 800);
  PlottingUtils::ConfigureCanvas(canvas, kFALSE);

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
  residuals->SetMarkerColor(kAzure);
  residuals->SetLineColor(kAzure);
  residuals->SetTitle("");
  Double_t actual_min = working_hist_->GetXaxis()->GetBinLowEdge(
      working_hist_->GetXaxis()->GetFirst());
  Double_t actual_max = working_hist_->GetXaxis()->GetBinUpEdge(
      working_hist_->GetXaxis()->GetLast());
  residuals->GetXaxis()->SetLimits(actual_min, actual_max);
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

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kFALSE);
  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(canvas, "log_" + peak_name + ".png", kFALSE);
}

Double_t FittingUtils::EstimateBackground() {
  Int_t left_bin = working_hist_->FindBin(fit_range_low_);
  Int_t right_bin = working_hist_->FindBin(fit_range_high_);

  Int_t n_sideband = (right_bin - left_bin) / 10;
  Double_t left_avg = 0, right_avg = 0;

  for (Int_t i = 0; i < n_sideband; i++) {
    left_avg += working_hist_->GetBinContent(left_bin + i);
    right_avg += working_hist_->GetBinContent(right_bin - i);
  }

  return (left_avg + right_avg) / (2.0 * n_sideband);
}

Double_t FittingUtils::ClampToBounds(Int_t param_index, Double_t value) {
  Double_t low, high;
  fit_function_->GetParLimits(param_index, low, high);

  if (low < high) {
    return TMath::Max(low, TMath::Min(value, high));
  }
  return value;
}

FitResultStandard FittingUtils::FitPeakStandard(const TString peak_name) {
  FitResultStandard results;

  TF1 *bkg_only = new TF1("bkg_temp", FittingFunctions::LinearBackground,
                          fit_range_low_, fit_range_high_, 2);

  Double_t exclude_low =
      fit_function_->GetParameter(0) - 3 * fit_function_->GetParameter(1);
  Double_t exclude_high =
      fit_function_->GetParameter(0) + 3 * fit_function_->GetParameter(1);
  working_hist_->Fit(bkg_only, "QN0+", "", fit_range_low_, exclude_low);

  fit_function_->SetParameter(3, ClampToBounds(3, bkg_only->GetParameter(0)));
  fit_function_->SetParameter(4, ClampToBounds(4, bkg_only->GetParameter(1)));

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMEN+");

  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFitStandard(peak_name);
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

FitResultDetailed FittingUtils::FitPeakDetailed(const TString peak_name) {
  FitResultDetailed results;

  fit_function_->FixParameter(5, 0);
  fit_function_->FixParameter(6, 0);
  fit_function_->FixParameter(7, 0.5);
  fit_function_->FixParameter(8, 0);
  fit_function_->FixParameter(9, 0.5);

  TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMNQ0+");

  if (!initial_fit.Get() || !initial_fit->IsValid()) {
    std::cout << "ERROR: Initial fit failed" << std::endl;
    results.mu_error = -1;
    return results;
  }

  Double_t chi2_standard = initial_fit->Chi2() / initial_fit->Ndf();
  std::cout << "Initial chi2/ndf = " << chi2_standard << std::endl;

  Double_t gaus_amp = TMath::Abs(fit_function_->GetParameter(2));
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());

  std::vector<Double_t> best_params(fit_function_->GetNpar());
  std::vector<Double_t> best_errors(fit_function_->GetNpar());
  for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
    best_params[i] = fit_function_->GetParameter(i);
    best_errors[i] = fit_function_->GetParError(i);
  }
  Double_t best_chi2 = chi2_standard;

  std::cout << "Peak height: " << peak_height << std::endl;
  std::cout << "Gaussian amplitude: " << gaus_amp << std::endl;

  if (use_low_tail_) {
    std::cout << "Testing low tail..." << std::endl;

    fit_function_->ReleaseParameter(6);
    fit_function_->ReleaseParameter(7);

    Double_t tail_amp_init = TMath::Min(gaus_amp * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(6, tail_amp_init);
    fit_function_->SetParameter(7, 0.3);
    fit_function_->SetParLimits(6, 0, peak_height * 0.30);
    fit_function_->SetParLimits(7, 0.01, 2.0);

    TFitResultPtr tail_fit = working_hist_->Fit(fit_function_, "LSMNQ0+");

    if (tail_fit.Get() && tail_fit->IsValid()) {
      Double_t chi2_with_tail = tail_fit->Chi2() / tail_fit->Ndf();
      Double_t improvement = (best_chi2 - chi2_with_tail) / best_chi2;

      Double_t fitted_amp = fit_function_->GetParameter(6);
      Double_t fitted_slope = fit_function_->GetParameter(7);

      std::cout << "Fitted tail amplitude: " << fitted_amp << std::endl;
      std::cout << "Fitted tail slope: " << fitted_slope << std::endl;
      std::cout << "Chi2/ndf: " << chi2_with_tail << " vs " << best_chi2
                << std::endl;
      std::cout << "Improvement: " << improvement * 100 << " percent"
                << std::endl;

      Bool_t amp_ok = (fitted_amp > 1e-6 && fitted_amp <= peak_height * 0.30);
      Bool_t slope_ok = (fitted_slope > 0.01 && fitted_slope < 2.0);
      Bool_t chi2_ok = (chi2_with_tail < best_chi2);
      Bool_t improvement_ok = (improvement > 0.002);

      if (amp_ok && slope_ok && chi2_ok && improvement_ok) {
        std::cout << "Low tail ACCEPTED" << std::endl;
        best_chi2 = chi2_with_tail;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Low tail REJECTED" << std::endl;
        fit_function_->FixParameter(6, 0);
        fit_function_->FixParameter(7, 0.5);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Low tail fit FAILED" << std::endl;
      fit_function_->FixParameter(6, 0);
      fit_function_->FixParameter(7, 0.5);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_high_tail_) {
    std::cout << "Testing high tail..." << std::endl;

    fit_function_->ReleaseParameter(8);
    fit_function_->ReleaseParameter(9);

    Double_t tail_amp_init = TMath::Min(gaus_amp * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(8, tail_amp_init);
    fit_function_->SetParameter(9, 0.3);
    fit_function_->SetParLimits(8, 0, peak_height * 0.30);
    fit_function_->SetParLimits(9, 0.01, 2.0);

    TFitResultPtr htail_fit = working_hist_->Fit(fit_function_, "LSMNQ0+");

    if (htail_fit.Get() && htail_fit->IsValid()) {
      Double_t chi2_with_htail = htail_fit->Chi2() / htail_fit->Ndf();
      Double_t improvement = (best_chi2 - chi2_with_htail) / best_chi2;

      Double_t fitted_amp = fit_function_->GetParameter(8);
      Double_t fitted_slope = fit_function_->GetParameter(9);

      std::cout << "Fitted tail amplitude: " << fitted_amp << std::endl;
      std::cout << "Fitted tail slope: " << fitted_slope << std::endl;
      std::cout << "Chi2/ndf: " << chi2_with_htail << " vs " << best_chi2
                << std::endl;
      std::cout << "Improvement: " << improvement * 100 << " percent"
                << std::endl;

      Bool_t amp_ok = (fitted_amp > 1e-6 && fitted_amp <= peak_height * 0.30);
      Bool_t slope_ok = (fitted_slope > 0.01 && fitted_slope < 2.0);
      Bool_t chi2_ok = (chi2_with_htail < best_chi2);
      Bool_t improvement_ok = (improvement > 0.002);

      if (amp_ok && slope_ok && chi2_ok && improvement_ok) {
        std::cout << "High tail ACCEPTED" << std::endl;
        best_chi2 = chi2_with_htail;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "High tail REJECTED" << std::endl;
        fit_function_->FixParameter(8, 0);
        fit_function_->FixParameter(9, 0.5);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "High tail fit FAILED" << std::endl;
      fit_function_->FixParameter(8, 0);
      fit_function_->FixParameter(9, 0.5);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_step_) {
    std::cout << "Testing step function..." << std::endl;

    fit_function_->ReleaseParameter(5);
    fit_function_->SetParameter(5, gaus_amp * 0.05);

    TFitResultPtr step_fit = working_hist_->Fit(fit_function_, "LSMNQ0+");

    if (step_fit.Get() && step_fit->IsValid()) {
      Double_t chi2_with_step = step_fit->Chi2() / step_fit->Ndf();
      Double_t improvement = (best_chi2 - chi2_with_step) / best_chi2;

      std::cout << "Chi2/ndf: " << chi2_with_step << " vs " << best_chi2
                << std::endl;
      std::cout << "Improvement: " << improvement * 100 << " percent"
                << std::endl;

      if (improvement > 0.002 && chi2_with_step < best_chi2) {
        std::cout << "Step ACCEPTED" << std::endl;
        best_chi2 = chi2_with_step;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Step REJECTED" << std::endl;
        fit_function_->FixParameter(5, 0);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Step fit FAILED" << std::endl;
      fit_function_->FixParameter(5, 0);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  std::cout << "Final fit with selected components..." << std::endl;
  for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
    fit_function_->SetParameter(i, best_params[i]);
    fit_function_->SetParError(i, best_errors[i]);
  }

  TFitResultPtr final_fit = working_hist_->Fit(fit_function_, "LSMEN+");

  if (final_fit.Get() && final_fit->IsValid()) {
    Double_t final_chi2 = final_fit->Chi2() / final_fit->Ndf();
    std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;

    std::cout << "Gaussian: YES" << std::endl;
    std::cout << "Background: YES" << std::endl;
    std::cout << "Low tail: "
              << (fit_function_->GetParameter(6) > 1e-6 ? "YES" : "NO")
              << std::endl;
    std::cout << "Step: "
              << (fit_function_->GetParameter(5) > 1e-6 ? "YES" : "NO")
              << std::endl;
    std::cout << "High tail: "
              << (fit_function_->GetParameter(8) > 1e-6 ? "YES" : "NO")
              << std::endl;

    PlotFitDetailed(peak_name);

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
  } else {
    std::cout << "ERROR: Final fit did not converge" << std::endl;
    results.mu_error = -1;
  }

  return results;
}
void FittingUtils::RegisterCustomFunctions() {
  TF1 *f_standard =
      new TF1("Standard", &FittingFunctions::Standard, 0, 1000, 5);
  f_standard->SetParName(0, "Mu");
  f_standard->SetParName(1, "Sigma");
  f_standard->SetParName(2, "GausAmplitude");
  f_standard->SetParName(3, "BkgConst");
  f_standard->SetParName(4, "BkgSlope");

  f_standard->SetParLimits(0, 0, 1e7);
  f_standard->SetParLimits(1, 0, 1e7);
  f_standard->SetParLimits(2, 0, 1e7);
  f_standard->SetParLimits(3, 0, 1e6);
  f_standard->SetParLimits(4, -1e3, 1e3);

  TF1 *f_detailed =
      new TF1("Detailed", &FittingFunctions::Detailed, 0, 1000, 10);
  f_detailed->SetParName(0, "Mu");
  f_detailed->SetParName(1, "Sigma");
  f_detailed->SetParName(2, "GausAmplitude");
  f_detailed->SetParName(3, "BkgConst");
  f_detailed->SetParName(4, "BkgSlope");
  f_detailed->SetParName(5, "StepAmplitude");
  f_detailed->SetParName(6, "LowTailAmplitude");
  f_detailed->SetParName(7, "LowTailRange");
  f_detailed->SetParName(8, "HighTailAmplitude");
  f_detailed->SetParName(9, "HighTailRange");

  f_detailed->SetParLimits(0, 0, 1e7);
  f_detailed->SetParLimits(1, 0, 1e7);
  f_detailed->SetParLimits(2, -1e-6, 1e7);
  f_detailed->SetParLimits(5, -1e-6, 1e7);
  f_detailed->SetParLimits(6, -1e-6, 1e7);
  f_detailed->SetParLimits(7, 0, 1e7);
  f_detailed->SetParLimits(8, -1e-6, 1e7);
  f_detailed->SetParLimits(9, 0, 1e7);

  std::cout << "Custom fitting functions registered and available in FitPanel!"
            << std::endl;
}
