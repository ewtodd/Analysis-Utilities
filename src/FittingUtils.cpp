#include "FittingUtils.hpp"

FittingUtils::FittingUtils(TH1 *working_hist, Float_t fit_range_low,
                           Float_t fit_range_high, Bool_t use_flat_background,
                           Bool_t isDetailed, Bool_t use_step,
                           Bool_t use_low_tail, Bool_t use_high_tail) {

  working_hist_ = static_cast<TH1F *>(working_hist->Clone());
  fit_range_low_ = fit_range_low;
  fit_range_high_ = fit_range_high;
  use_flat_background_ = use_flat_background;
  isDetailed_ = isDetailed;
  use_step_ = use_step;
  use_low_tail_ = use_low_tail;
  use_high_tail_ = use_high_tail;
  use_manual_init_ = kFALSE;

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
    if (!use_flat_background_) {
      fit_function_->SetParLimits(4, -0.1 * bkg_estimate / range_width,
                                  0.1 * bkg_estimate / range_width);
    } else
      fit_function_->SetParLimits(4, 0, 0);

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
    if (!use_flat_background_) {
      fit_function_->SetParLimits(4, -0.1 * bkg_estimate / range_width,
                                  0.1 * bkg_estimate / range_width);
    } else
      fit_function_->SetParLimits(4, 0, 0);

    fit_function_->SetParameter(0, mu_init);
    fit_function_->SetParameter(1, sigma_init);
    fit_function_->SetParameter(2, peak_height * 0.8);
    fit_function_->SetParameter(3, bkg_estimate);
    fit_function_->SetParameter(4, 0);

    fit_function_->SetParLimits(5, 0, peak_height * 0.5);
    fit_function_->SetParameter(5, 0);

    fit_function_->SetParLimits(6, 0, peak_height);
    fit_function_->SetParLimits(7, 1, 10);
    fit_function_->SetParameter(6, peak_height * 0.10);
    fit_function_->SetParameter(7, 0.2);

    fit_function_->SetParLimits(8, 0, peak_height);
    fit_function_->SetParLimits(9, 1, 10);
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

Double_t FittingFunctions::DoublePeakStandard(Double_t *x, Double_t *par) {
  Double_t mu1 = par[0];
  Double_t sigma1 = par[1];
  Double_t gaus_amplitude1 = par[2];
  Double_t mu2 = par[3];
  Double_t sigma2 = par[4];
  Double_t gaus_amplitude2 = par[5];
  Double_t bkg_const = par[6];
  Double_t bkg_slope = par[7];

  Double_t gaus1_par[3] = {mu1, sigma1, gaus_amplitude1};
  Double_t gaus2_par[3] = {mu2, sigma2, gaus_amplitude2};
  Double_t bkg_par[2] = {bkg_const, bkg_slope};

  return Gaussian(x, gaus1_par) + Gaussian(x, gaus2_par) +
         LinearBackground(x, bkg_par);
}

Double_t FittingFunctions::DoublePeakDetailed(Double_t *x, Double_t *par) {
  Double_t mu1 = par[0];
  Double_t sigma1 = par[1];
  Double_t gaus_amplitude1 = par[2];
  Double_t step_amplitude1 = par[3];
  Double_t low_tail_amplitude1 = par[4];
  Double_t low_tail_range1 = par[5];
  Double_t high_tail_amplitude1 = par[6];
  Double_t high_tail_range1 = par[7];

  Double_t mu2 = par[8];
  Double_t sigma2 = par[9];
  Double_t gaus_amplitude2 = par[10];
  Double_t step_amplitude2 = par[11];
  Double_t low_tail_amplitude2 = par[12];
  Double_t low_tail_range2 = par[13];
  Double_t high_tail_amplitude2 = par[14];
  Double_t high_tail_range2 = par[15];

  Double_t bkg_const = par[16];
  Double_t bkg_slope = par[17];

  Double_t gaus1_par[3] = {mu1, sigma1, gaus_amplitude1};
  Double_t step1_par[3] = {mu1, sigma1, step_amplitude1};
  Double_t low_tail1_par[4] = {mu1, sigma1, low_tail_amplitude1,
                               low_tail_range1};
  Double_t high_tail1_par[4] = {mu1, sigma1, high_tail_amplitude1,
                                high_tail_range1};

  Double_t gaus2_par[3] = {mu2, sigma2, gaus_amplitude2};
  Double_t step2_par[3] = {mu2, sigma2, step_amplitude2};
  Double_t low_tail2_par[4] = {mu2, sigma2, low_tail_amplitude2,
                               low_tail_range2};
  Double_t high_tail2_par[4] = {mu2, sigma2, high_tail_amplitude2,
                                high_tail_range2};

  Double_t bkg_par[2] = {bkg_const, bkg_slope};

  return Gaussian(x, gaus1_par) + Step(x, step1_par) +
         LowTail(x, low_tail1_par) + HighTail(x, high_tail1_par) +
         Gaussian(x, gaus2_par) + Step(x, step2_par) +
         LowTail(x, low_tail2_par) + HighTail(x, high_tail2_par) +
         LinearBackground(x, bkg_par);
}

Double_t FittingFunctions::TriplePeakStandard(Double_t *x, Double_t *par) {
  Double_t mu1 = par[0];
  Double_t sigma1 = par[1];
  Double_t gaus_amplitude1 = par[2];
  Double_t mu2 = par[3];
  Double_t sigma2 = par[4];
  Double_t gaus_amplitude2 = par[5];
  Double_t mu3 = par[6];
  Double_t sigma3 = par[7];
  Double_t gaus_amplitude3 = par[8];
  Double_t bkg_const = par[9];
  Double_t bkg_slope = par[10];

  Double_t gaus1_par[3] = {mu1, sigma1, gaus_amplitude1};
  Double_t gaus2_par[3] = {mu2, sigma2, gaus_amplitude2};
  Double_t gaus3_par[3] = {mu3, sigma3, gaus_amplitude3};
  Double_t bkg_par[2] = {bkg_const, bkg_slope};

  return Gaussian(x, gaus1_par) + Gaussian(x, gaus2_par) +
         Gaussian(x, gaus3_par) + LinearBackground(x, bkg_par);
}

Double_t FittingFunctions::TriplePeakDetailed(Double_t *x, Double_t *par) {
  Double_t mu1 = par[0];
  Double_t sigma1 = par[1];
  Double_t gaus_amplitude1 = par[2];
  Double_t step_amplitude1 = par[3];
  Double_t low_tail_amplitude1 = par[4];
  Double_t low_tail_range1 = par[5];
  Double_t high_tail_amplitude1 = par[6];
  Double_t high_tail_range1 = par[7];

  Double_t mu2 = par[8];
  Double_t sigma2 = par[9];
  Double_t gaus_amplitude2 = par[10];
  Double_t step_amplitude2 = par[11];
  Double_t low_tail_amplitude2 = par[12];
  Double_t low_tail_range2 = par[13];
  Double_t high_tail_amplitude2 = par[14];
  Double_t high_tail_range2 = par[15];

  Double_t mu3 = par[16];
  Double_t sigma3 = par[17];
  Double_t gaus_amplitude3 = par[18];
  Double_t step_amplitude3 = par[19];
  Double_t low_tail_amplitude3 = par[20];
  Double_t low_tail_range3 = par[21];
  Double_t high_tail_amplitude3 = par[22];
  Double_t high_tail_range3 = par[23];

  Double_t bkg_const = par[24];
  Double_t bkg_slope = par[25];

  Double_t gaus1_par[3] = {mu1, sigma1, gaus_amplitude1};
  Double_t step1_par[3] = {mu1, sigma1, step_amplitude1};
  Double_t low_tail1_par[4] = {mu1, sigma1, low_tail_amplitude1,
                               low_tail_range1};
  Double_t high_tail1_par[4] = {mu1, sigma1, high_tail_amplitude1,
                                high_tail_range1};

  Double_t gaus2_par[3] = {mu2, sigma2, gaus_amplitude2};
  Double_t step2_par[3] = {mu2, sigma2, step_amplitude2};
  Double_t low_tail2_par[4] = {mu2, sigma2, low_tail_amplitude2,
                               low_tail_range2};
  Double_t high_tail2_par[4] = {mu2, sigma2, high_tail_amplitude2,
                                high_tail_range2};

  Double_t gaus3_par[3] = {mu3, sigma3, gaus_amplitude3};
  Double_t step3_par[3] = {mu3, sigma3, step_amplitude3};
  Double_t low_tail3_par[4] = {mu3, sigma3, low_tail_amplitude3,
                               low_tail_range3};
  Double_t high_tail3_par[4] = {mu3, sigma3, high_tail_amplitude3,
                                high_tail_range3};

  Double_t bkg_par[2] = {bkg_const, bkg_slope};

  return Gaussian(x, gaus1_par) + Step(x, step1_par) +
         LowTail(x, low_tail1_par) + HighTail(x, high_tail1_par) +
         Gaussian(x, gaus2_par) + Step(x, step2_par) +
         LowTail(x, low_tail2_par) + HighTail(x, high_tail2_par) +
         Gaussian(x, gaus3_par) + Step(x, step3_par) +
         LowTail(x, low_tail3_par) + HighTail(x, high_tail3_par) +
         LinearBackground(x, bkg_par);
}

void FittingUtils::SwapDoublePeakStandardParameters() {
  std::cout << "Swapping double peak standard parameters to enforce mu1 < mu2"
            << std::endl;

  Double_t temp_mu = fit_function_->GetParameter(0);
  Double_t temp_mu_err = fit_function_->GetParError(0);
  fit_function_->SetParameter(0, fit_function_->GetParameter(3));
  fit_function_->SetParError(0, fit_function_->GetParError(3));
  fit_function_->SetParameter(3, temp_mu);
  fit_function_->SetParError(3, temp_mu_err);

  Double_t temp_sigma = fit_function_->GetParameter(1);
  Double_t temp_sigma_err = fit_function_->GetParError(1);
  fit_function_->SetParameter(1, fit_function_->GetParameter(4));
  fit_function_->SetParError(1, fit_function_->GetParError(4));
  fit_function_->SetParameter(4, temp_sigma);
  fit_function_->SetParError(4, temp_sigma_err);

  Double_t temp_amp = fit_function_->GetParameter(2);
  Double_t temp_amp_err = fit_function_->GetParError(2);
  fit_function_->SetParameter(2, fit_function_->GetParameter(5));
  fit_function_->SetParError(2, fit_function_->GetParError(5));
  fit_function_->SetParameter(5, temp_amp);
  fit_function_->SetParError(5, temp_amp_err);
}

void FittingUtils::SwapDoublePeakDetailedParameters() {
  std::cout << "Swapping double peak detailed parameters to enforce mu1 < mu2"
            << std::endl;

  for (Int_t i = 0; i < 8; i++) {
    Double_t temp_val = fit_function_->GetParameter(i);
    Double_t temp_err = fit_function_->GetParError(i);

    fit_function_->SetParameter(i, fit_function_->GetParameter(8 + i));
    fit_function_->SetParError(i, fit_function_->GetParError(8 + i));

    fit_function_->SetParameter(8 + i, temp_val);
    fit_function_->SetParError(8 + i, temp_err);
  }

  std::cout << "  Swapped: Mu, Sigma, GausAmp, StepAmp, LowTailAmp, "
               "LowTailSlope, HighTailAmp, HighTailSlope"
            << std::endl;
}

void FittingUtils::SetManualParameters(const std::vector<Double_t> &params) {
  if (params.size() != (size_t)fit_function_->GetNpar()) {
    std::cerr << "ERROR: Manual parameters size (" << params.size()
              << ") doesn't match number of fit parameters ("
              << fit_function_->GetNpar() << ")" << std::endl;
    return;
  }

  manual_params_ = params;
  use_manual_init_ = kTRUE;

  // Apply the parameters immediately
  for (size_t i = 0; i < params.size(); i++) {
    fit_function_->SetParameter(i, params[i]);
  }

  std::cout << "Manual parameters set:" << std::endl;
  for (size_t i = 0; i < params.size(); i++) {
    std::cout << "  Par[" << i << "] " << fit_function_->GetParName(i) << " = "
              << params[i] << std::endl;
  }
}

void FittingUtils::SetManualParameter(Int_t index, Double_t value) {
  if (index < 0 || index >= fit_function_->GetNpar()) {
    std::cerr << "ERROR: Parameter index " << index << " out of range [0, "
              << fit_function_->GetNpar() - 1 << "]" << std::endl;
    return;
  }

  if (!use_manual_init_) {
    manual_params_.resize(fit_function_->GetNpar(), 0.0);
    use_manual_init_ = kTRUE;
  }

  manual_params_[index] = value;
  fit_function_->SetParameter(index, value);

  std::cout << "Set Par[" << index << "] " << fit_function_->GetParName(index)
            << " = " << value << std::endl;
}

void FittingUtils::PlotFitStandard(const TString input_name,
                                   const TString peak_name) {
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
  peak->SetNpx(1000);
  peak->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(3));
  background->SetParameter(1, fit_function_->GetParameter(4));
  background->SetLineColor(kGreen);
  background->SetNpx(1000);
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
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(
      canvas, "log_" + peak_name + "_" + input_name + ".png", kFALSE);
}

void FittingUtils::PlotFitDetailed(const TString input_name,
                                   const TString peak_name) {
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
  peak->SetNpx(1000);
  peak->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(3));
  background->SetParameter(1, fit_function_->GetParameter(4));
  background->SetLineColor(kGreen);
  background->SetNpx(1000);
  background->Draw("same");

  TF1 *step = new TF1("step", FittingFunctions::Step, fit_range_low_,
                      fit_range_high_, 3);
  step->SetParameter(0, fit_function_->GetParameter(0));
  step->SetParameter(1, fit_function_->GetParameter(1));
  step->SetParameter(2, fit_function_->GetParameter(5));
  step->SetLineColor(kGray);
  step->SetNpx(1000);
  step->Draw("same");

  TF1 *low_tail = new TF1("lowtail", FittingFunctions::LowTail, fit_range_low_,
                          fit_range_high_, 4);
  low_tail->SetParameter(0, fit_function_->GetParameter(0));
  low_tail->SetParameter(1, fit_function_->GetParameter(1));
  low_tail->SetParameter(2, fit_function_->GetParameter(6));
  low_tail->SetParameter(3, fit_function_->GetParameter(7));
  low_tail->SetLineColor(kRed);
  low_tail->SetNpx(1000);
  low_tail->Draw("same");

  TF1 *high_tail = new TF1("hightail", FittingFunctions::HighTail,
                           fit_range_low_, fit_range_high_, 4);
  high_tail->SetParameter(0, fit_function_->GetParameter(0));
  high_tail->SetParameter(1, fit_function_->GetParameter(1));
  high_tail->SetParameter(2, fit_function_->GetParameter(8));
  high_tail->SetParameter(3, fit_function_->GetParameter(9));
  high_tail->SetLineColor(kOrange);
  high_tail->SetNpx(1000);
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
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(
      canvas, "log_" + peak_name + "_" + input_name + ".png", kFALSE);
}

void FittingUtils::PlotFitDoublePeakStandard(const TString input_name,
                                             const TString peak_name) {
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

  TF1 *peak1 = new TF1("gaussian1", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak1->SetParameter(0, fit_function_->GetParameter(0));
  peak1->SetParameter(1, fit_function_->GetParameter(1));
  peak1->SetParameter(2, fit_function_->GetParameter(2));
  peak1->SetLineColor(kBlack);
  peak1->SetNpx(1000);
  peak1->Draw("same");

  TF1 *peak2 = new TF1("gaussian2", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak2->SetParameter(0, fit_function_->GetParameter(3));
  peak2->SetParameter(1, fit_function_->GetParameter(4));
  peak2->SetParameter(2, fit_function_->GetParameter(5));
  peak2->SetLineColor(kBlack);
  peak2->SetNpx(1000);
  peak2->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(6));
  background->SetParameter(1, fit_function_->GetParameter(7));
  background->SetLineColor(kGreen);
  background->SetNpx(1000);
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
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(
      canvas, "log_" + peak_name + "_" + input_name + ".png", kFALSE);
}

void FittingUtils::PlotFitDoublePeakDetailed(const TString input_name,
                                             const TString peak_name) {
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

  TF1 *peak1 = new TF1("gaussian1", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak1->SetParameter(0, fit_function_->GetParameter(0));
  peak1->SetParameter(1, fit_function_->GetParameter(1));
  peak1->SetParameter(2, fit_function_->GetParameter(2));
  peak1->SetLineColor(kBlack);
  peak1->SetNpx(1000);
  peak1->Draw("same");

  TF1 *peak2 = new TF1("gaussian2", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak2->SetParameter(0, fit_function_->GetParameter(8));
  peak2->SetParameter(1, fit_function_->GetParameter(9));
  peak2->SetParameter(2, fit_function_->GetParameter(10));
  peak2->SetLineColor(kBlack);
  peak2->SetNpx(1000);
  peak2->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(16));
  background->SetParameter(1, fit_function_->GetParameter(17));
  background->SetLineColor(kGreen);
  background->SetNpx(1000);
  background->Draw("same");

  TF1 *step1 = new TF1("step1", FittingFunctions::Step, fit_range_low_,
                       fit_range_high_, 3);
  step1->SetParameter(0, fit_function_->GetParameter(0));
  step1->SetParameter(1, fit_function_->GetParameter(1));
  step1->SetParameter(2, fit_function_->GetParameter(3));
  step1->SetLineColor(kGray);
  step1->SetNpx(1000);
  step1->Draw("same");

  TF1 *low_tail1 = new TF1("lowtail1", FittingFunctions::LowTail,
                           fit_range_low_, fit_range_high_, 4);
  low_tail1->SetParameter(0, fit_function_->GetParameter(0));
  low_tail1->SetParameter(1, fit_function_->GetParameter(1));
  low_tail1->SetParameter(2, fit_function_->GetParameter(4));
  low_tail1->SetParameter(3, fit_function_->GetParameter(5));
  low_tail1->SetLineColor(kRed);
  low_tail1->SetNpx(1000);
  low_tail1->Draw("same");

  TF1 *high_tail1 = new TF1("hightail1", FittingFunctions::HighTail,
                            fit_range_low_, fit_range_high_, 4);
  high_tail1->SetParameter(0, fit_function_->GetParameter(0));
  high_tail1->SetParameter(1, fit_function_->GetParameter(1));
  high_tail1->SetParameter(2, fit_function_->GetParameter(6));
  high_tail1->SetParameter(3, fit_function_->GetParameter(7));
  high_tail1->SetLineColor(kOrange);
  high_tail1->SetNpx(1000);
  high_tail1->Draw("same");

  TF1 *step2 = new TF1("step2", FittingFunctions::Step, fit_range_low_,
                       fit_range_high_, 3);
  step2->SetParameter(0, fit_function_->GetParameter(8));
  step2->SetParameter(1, fit_function_->GetParameter(9));
  step2->SetParameter(2, fit_function_->GetParameter(11));
  step2->SetLineColor(kGray);
  step2->SetNpx(1000);
  step2->Draw("same");

  TF1 *low_tail2 = new TF1("lowtail2", FittingFunctions::LowTail,
                           fit_range_low_, fit_range_high_, 4);
  low_tail2->SetParameter(0, fit_function_->GetParameter(8));
  low_tail2->SetParameter(1, fit_function_->GetParameter(9));
  low_tail2->SetParameter(2, fit_function_->GetParameter(12));
  low_tail2->SetParameter(3, fit_function_->GetParameter(13));
  low_tail2->SetLineColor(kRed);
  low_tail2->SetNpx(1000);
  low_tail2->Draw("same");

  TF1 *high_tail2 = new TF1("hightail2", FittingFunctions::HighTail,
                            fit_range_low_, fit_range_high_, 4);
  high_tail2->SetParameter(0, fit_function_->GetParameter(8));
  high_tail2->SetParameter(1, fit_function_->GetParameter(9));
  high_tail2->SetParameter(2, fit_function_->GetParameter(14));
  high_tail2->SetParameter(3, fit_function_->GetParameter(15));
  high_tail2->SetLineColor(kOrange);
  high_tail2->SetNpx(1000);
  high_tail2->Draw("same");

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
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(
      canvas, "log_" + peak_name + "_" + input_name + ".png", kFALSE);
}

void FittingUtils::PlotFitTriplePeakStandard(const TString input_name,
                                             const TString peak_name) {
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

  TF1 *peak1 = new TF1("gaussian1", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak1->SetParameter(0, fit_function_->GetParameter(0));
  peak1->SetParameter(1, fit_function_->GetParameter(1));
  peak1->SetParameter(2, fit_function_->GetParameter(2));
  peak1->SetLineColor(kBlack);
  peak1->SetNpx(1000);
  peak1->Draw("same");

  TF1 *peak2 = new TF1("gaussian2", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak2->SetParameter(0, fit_function_->GetParameter(3));
  peak2->SetParameter(1, fit_function_->GetParameter(4));
  peak2->SetParameter(2, fit_function_->GetParameter(5));
  peak2->SetLineColor(kBlack);
  peak2->SetNpx(1000);
  peak2->Draw("same");

  TF1 *peak3 = new TF1("gaussian3", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak3->SetParameter(0, fit_function_->GetParameter(6));
  peak3->SetParameter(1, fit_function_->GetParameter(7));
  peak3->SetParameter(2, fit_function_->GetParameter(8));
  peak3->SetLineColor(kMagenta);
  peak3->SetNpx(1000);
  peak3->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(9));
  background->SetParameter(1, fit_function_->GetParameter(10));
  background->SetLineColor(kGreen);
  background->SetNpx(1000);
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
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(
      canvas, "log_" + peak_name + "_" + input_name + ".png", kFALSE);
}

void FittingUtils::PlotFitTriplePeakDetailed(const TString input_name,
                                             const TString peak_name) {
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

  TF1 *peak1 = new TF1("gaussian1", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak1->SetParameter(0, fit_function_->GetParameter(0));
  peak1->SetParameter(1, fit_function_->GetParameter(1));
  peak1->SetParameter(2, fit_function_->GetParameter(2));
  peak1->SetLineColor(kBlack);
  peak1->SetNpx(1000);
  peak1->Draw("same");

  TF1 *peak2 = new TF1("gaussian2", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak2->SetParameter(0, fit_function_->GetParameter(8));
  peak2->SetParameter(1, fit_function_->GetParameter(9));
  peak2->SetParameter(2, fit_function_->GetParameter(10));
  peak2->SetLineColor(kBlack);
  peak2->SetNpx(1000);
  peak2->Draw("same");

  TF1 *peak3 = new TF1("gaussian3", FittingFunctions::Gaussian, fit_range_low_,
                       fit_range_high_, 3);
  peak3->SetParameter(0, fit_function_->GetParameter(16));
  peak3->SetParameter(1, fit_function_->GetParameter(17));
  peak3->SetParameter(2, fit_function_->GetParameter(18));
  peak3->SetLineColor(kMagenta);
  peak3->SetNpx(1000);
  peak3->Draw("same");

  TF1 *background = new TF1("background", FittingFunctions::LinearBackground,
                            fit_range_low_, fit_range_high_, 2);
  background->SetParameter(0, fit_function_->GetParameter(24));
  background->SetParameter(1, fit_function_->GetParameter(25));
  background->SetLineColor(kGreen);
  background->SetNpx(1000);
  background->Draw("same");

  TF1 *step1 = new TF1("step1", FittingFunctions::Step, fit_range_low_,
                       fit_range_high_, 3);
  step1->SetParameter(0, fit_function_->GetParameter(0));
  step1->SetParameter(1, fit_function_->GetParameter(1));
  step1->SetParameter(2, fit_function_->GetParameter(3));
  step1->SetLineColor(kGray);
  step1->SetNpx(1000);
  step1->Draw("same");

  TF1 *low_tail1 = new TF1("lowtail1", FittingFunctions::LowTail,
                           fit_range_low_, fit_range_high_, 4);
  low_tail1->SetParameter(0, fit_function_->GetParameter(0));
  low_tail1->SetParameter(1, fit_function_->GetParameter(1));
  low_tail1->SetParameter(2, fit_function_->GetParameter(4));
  low_tail1->SetParameter(3, fit_function_->GetParameter(5));
  low_tail1->SetLineColor(kRed);
  low_tail1->SetNpx(1000);
  low_tail1->Draw("same");

  TF1 *high_tail1 = new TF1("hightail1", FittingFunctions::HighTail,
                            fit_range_low_, fit_range_high_, 4);
  high_tail1->SetParameter(0, fit_function_->GetParameter(0));
  high_tail1->SetParameter(1, fit_function_->GetParameter(1));
  high_tail1->SetParameter(2, fit_function_->GetParameter(6));
  high_tail1->SetParameter(3, fit_function_->GetParameter(7));
  high_tail1->SetLineColor(kOrange);
  high_tail1->SetNpx(1000);
  high_tail1->Draw("same");

  TF1 *step2 = new TF1("step2", FittingFunctions::Step, fit_range_low_,
                       fit_range_high_, 3);
  step2->SetParameter(0, fit_function_->GetParameter(8));
  step2->SetParameter(1, fit_function_->GetParameter(9));
  step2->SetParameter(2, fit_function_->GetParameter(11));
  step2->SetLineColor(kGray);
  step2->SetNpx(1000);
  step2->Draw("same");

  TF1 *low_tail2 = new TF1("lowtail2", FittingFunctions::LowTail,
                           fit_range_low_, fit_range_high_, 4);
  low_tail2->SetParameter(0, fit_function_->GetParameter(8));
  low_tail2->SetParameter(1, fit_function_->GetParameter(9));
  low_tail2->SetParameter(2, fit_function_->GetParameter(12));
  low_tail2->SetParameter(3, fit_function_->GetParameter(13));
  low_tail2->SetLineColor(kRed);
  low_tail2->SetNpx(1000);
  low_tail2->Draw("same");

  TF1 *high_tail2 = new TF1("hightail2", FittingFunctions::HighTail,
                            fit_range_low_, fit_range_high_, 4);
  high_tail2->SetParameter(0, fit_function_->GetParameter(8));
  high_tail2->SetParameter(1, fit_function_->GetParameter(9));
  high_tail2->SetParameter(2, fit_function_->GetParameter(14));
  high_tail2->SetParameter(3, fit_function_->GetParameter(15));
  high_tail2->SetLineColor(kOrange);
  high_tail2->SetNpx(1000);
  high_tail2->Draw("same");

  TF1 *step3 = new TF1("step3", FittingFunctions::Step, fit_range_low_,
                       fit_range_high_, 3);
  step3->SetParameter(0, fit_function_->GetParameter(16));
  step3->SetParameter(1, fit_function_->GetParameter(17));
  step3->SetParameter(2, fit_function_->GetParameter(19));
  step3->SetLineColor(kGray);
  step3->SetNpx(1000);
  step3->Draw("same");

  TF1 *low_tail3 = new TF1("lowtail3", FittingFunctions::LowTail,
                           fit_range_low_, fit_range_high_, 4);
  low_tail3->SetParameter(0, fit_function_->GetParameter(16));
  low_tail3->SetParameter(1, fit_function_->GetParameter(17));
  low_tail3->SetParameter(2, fit_function_->GetParameter(20));
  low_tail3->SetParameter(3, fit_function_->GetParameter(21));
  low_tail3->SetLineColor(kRed);
  low_tail3->SetNpx(1000);
  low_tail3->Draw("same");

  TF1 *high_tail3 = new TF1("hightail3", FittingFunctions::HighTail,
                            fit_range_low_, fit_range_high_, 4);
  high_tail3->SetParameter(0, fit_function_->GetParameter(16));
  high_tail3->SetParameter(1, fit_function_->GetParameter(17));
  high_tail3->SetParameter(2, fit_function_->GetParameter(22));
  high_tail3->SetParameter(3, fit_function_->GetParameter(23));
  high_tail3->SetLineColor(kOrange);
  high_tail3->SetNpx(1000);
  high_tail3->Draw("same");

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
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(kTRUE);
  PlottingUtils::SaveFigure(
      canvas, "log_" + peak_name + "_" + input_name + ".png", kFALSE);
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

FitResultStandard FittingUtils::FitPeakStandard(const TString input_name,
                                                const TString peak_name) {
  FitResultStandard results;

  if (!use_manual_init_) {
    if (!use_flat_background_) {
      TF1 *bkg_only = new TF1("bkg_temp", FittingFunctions::LinearBackground,
                              fit_range_low_, fit_range_high_, 2);

      Double_t exclude_low =
          fit_function_->GetParameter(0) - 3 * fit_function_->GetParameter(1);
      Double_t exclude_high =
          fit_function_->GetParameter(0) + 3 * fit_function_->GetParameter(1);
      working_hist_->Fit(bkg_only, "QN0R", "", fit_range_low_, exclude_low);

      fit_function_->SetParameter(3,
                                  ClampToBounds(3, bkg_only->GetParameter(0)));
      fit_function_->SetParameter(4,
                                  ClampToBounds(4, bkg_only->GetParameter(1)));

      delete bkg_only;
    } else {
      fit_function_->FixParameter(4, 0);
    }
  } else {
    std::cout << "Using manually initialized parameters for Standard fit"
              << std::endl;
    for (size_t i = 0; i < manual_params_.size(); i++) {
      fit_function_->SetParameter(i, manual_params_[i]);
      std::cout << "  Par[" << i << "] " << fit_function_->GetParName(i)
                << " = " << manual_params_[i] << std::endl;
    }

    if (use_flat_background_) {
      fit_function_->FixParameter(4, 0);
      std::cout << "  Fixing BkgSlope to 0 (flat background mode)" << std::endl;
    }
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    std::cout << "Standard fit converged successfully" << std::endl;
    std::cout << "Chi2/ndf = " << fit_result->Chi2() / fit_result->Ndf()
              << std::endl;

    PlotFitStandard(input_name, peak_name);

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
    results.reduced_chi2 = fit_result->Chi2() / fit_result->Ndf();
  } else {
    std::cout << "ERROR: Standard fit failed to converge" << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.mu_error = -1;
  }

  return results;
}

FitResultDetailed FittingUtils::FitPeakDetailed(const TString input_name,
                                                const TString peak_name) {
  FitResultDetailed results;

  fit_function_->FixParameter(5, 0);
  fit_function_->FixParameter(6, 0);
  fit_function_->FixParameter(7, 1);
  fit_function_->FixParameter(8, 0);
  fit_function_->FixParameter(9, 1);

  if (use_flat_background_) {
    fit_function_->FixParameter(4, 0);
  }

  if (use_manual_init_) {
    std::cout << "Using manually initialized parameters" << std::endl;
    for (size_t i = 0; i < manual_params_.size(); i++) {
      fit_function_->SetParameter(i, manual_params_[i]);
    }

    std::cout << "Skipping auto-initialization, using provided values"
              << std::endl;
  } else {
    TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (!initial_fit.Get() || !initial_fit->IsValid()) {
      std::cout << "ERROR: Initial fit failed" << std::endl;
      results.mu_error = -1;
      return results;
    }

    Double_t chi2_standard = initial_fit->Chi2() / initial_fit->Ndf();
    std::cout << "Initial chi2/ndf = " << chi2_standard << std::endl;
  }

  Double_t gaus_amp = TMath::Abs(fit_function_->GetParameter(2));
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());

  std::vector<Double_t> best_params(fit_function_->GetNpar());
  std::vector<Double_t> best_errors(fit_function_->GetNpar());
  for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
    best_params[i] = fit_function_->GetParameter(i);
    best_errors[i] = fit_function_->GetParError(i);
  }

  TFitResultPtr current_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

  Double_t best_chi2 = (current_fit.Get() && current_fit->IsValid())
                           ? current_fit->Chi2() / current_fit->Ndf()
                           : 1e9;

  std::cout << "Starting chi2/ndf = " << best_chi2 << std::endl;
  std::cout << "Peak height: " << peak_height << std::endl;
  std::cout << "Gaussian amplitude: " << gaus_amp << std::endl;
  std::cout << "Background mode: " << (use_flat_background_ ? "FLAT" : "LINEAR")
            << std::endl;

  if (use_step_) {
    std::cout << "Testing step function..." << std::endl;

    fit_function_->ReleaseParameter(5);
    fit_function_->SetParLimits(5, 0, peak_height);

    if (!use_manual_init_) {
      fit_function_->SetParameter(5, gaus_amp);
    }

    TFitResultPtr step_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (step_fit.Get() && step_fit->IsValid()) {
      Double_t chi2_with_step = step_fit->Chi2() / step_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_step << " vs " << best_chi2
                << std::endl;

      if (chi2_with_step < best_chi2) {
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

  if (use_low_tail_) {
    std::cout << "Testing low tail..." << std::endl;

    fit_function_->ReleaseParameter(6);
    fit_function_->ReleaseParameter(7);
    fit_function_->SetParLimits(6, 0, peak_height * 1.2);
    fit_function_->SetParLimits(7, 1, 10);

    if (!use_manual_init_) {
      Double_t tail_amp_init = TMath::Min(gaus_amp * 0.15, peak_height * 0.25);
      fit_function_->SetParameter(6, tail_amp_init);
      fit_function_->SetParameter(7, 1);
    }

    TFitResultPtr tail_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (tail_fit.Get() && tail_fit->IsValid()) {
      Double_t chi2_with_tail = tail_fit->Chi2() / tail_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_tail << " vs " << best_chi2
                << std::endl;

      if (chi2_with_tail < best_chi2) {
        std::cout << "Low tail ACCEPTED" << std::endl;
        best_chi2 = chi2_with_tail;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Low tail REJECTED" << std::endl;
        fit_function_->FixParameter(6, 0);
        fit_function_->FixParameter(7, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Low tail fit FAILED" << std::endl;
      fit_function_->FixParameter(6, 0);
      fit_function_->FixParameter(7, 1);
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
    fit_function_->SetParLimits(8, 0, peak_height * 1.2);
    fit_function_->SetParLimits(9, 1, 10);

    if (!use_manual_init_) {
      Double_t tail_amp_init = TMath::Min(gaus_amp * 0.15, peak_height * 0.25);
      fit_function_->SetParameter(8, tail_amp_init);
      fit_function_->SetParameter(9, 1);
    }

    TFitResultPtr htail_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (htail_fit.Get() && htail_fit->IsValid()) {
      Double_t chi2_with_htail = htail_fit->Chi2() / htail_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_htail << " vs " << best_chi2
                << std::endl;

      if (chi2_with_htail < best_chi2) {
        std::cout << "High tail ACCEPTED" << std::endl;
        best_chi2 = chi2_with_htail;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "High tail REJECTED" << std::endl;
        fit_function_->FixParameter(8, 0);
        fit_function_->FixParameter(9, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "High tail fit FAILED" << std::endl;
      fit_function_->FixParameter(8, 0);
      fit_function_->FixParameter(9, 1);
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

  if (use_flat_background_) {
    fit_function_->FixParameter(4, 0);
  }

  TFitResultPtr final_fit = working_hist_->Fit(fit_function_, "LSMRBENR+");

  if (final_fit.Get() && final_fit->IsValid()) {
    Double_t final_chi2 = final_fit->Chi2() / final_fit->Ndf();
    std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;

    std::cout << "Gaussian: YES" << std::endl;
    std::cout << "Background: " << (use_flat_background_ ? "FLAT" : "LINEAR")
              << std::endl;
    std::cout << "Step: "
              << (fit_function_->GetParameter(5) > 1e-6 ? "YES" : "NO")
              << std::endl;
    std::cout << "Low tail: "
              << (fit_function_->GetParameter(6) > 1e-6 ? "YES" : "NO")
              << std::endl;
    std::cout << "High tail: "
              << (fit_function_->GetParameter(8) > 1e-6 ? "YES" : "NO")
              << std::endl;

    PlotFitDetailed(input_name, peak_name);

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
    results.reduced_chi2 = final_chi2;
  } else {
    std::cout << "ERROR: Final fit did not converge" << std::endl;
    results.mu_error = -1;
  }

  return results;
}

FitResultDoublePeakStandard
FittingUtils::FitDoublePeakStandard(const TString input_name,
                                    const TString peak_name, Double_t mu1_init,
                                    Double_t mu2_init) {
  FitResultDoublePeakStandard results;

  if (mu1_init > mu2_init) {
    std::cout << "Warning: mu1_init > mu2_init, swapping initial values"
              << std::endl;
    Double_t temp = mu1_init;
    mu1_init = mu2_init;
    mu2_init = temp;
  }

  fit_function_ =
      new TF1("DoublePeakStandard", &FittingFunctions::DoublePeakStandard,
              fit_range_low_, fit_range_high_, 8);

  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "Mu2");
  fit_function_->SetParName(4, "Sigma2");
  fit_function_->SetParName(5, "GausAmplitude2");
  fit_function_->SetParName(6, "BkgConst");
  fit_function_->SetParName(7, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  fit_function_->SetParLimits(0, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(1, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(2, 0, peak_height * 1.5);
  fit_function_->SetParLimits(3, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(4, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(5, 0, peak_height * 1.5);
  fit_function_->SetParLimits(6, 0, peak_height * 0.5);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(7, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  } else {
    fit_function_->SetParLimits(7, 0, 0);
  }

  fit_function_->SetParameter(0, mu1_init);
  fit_function_->SetParameter(1, sigma_init);
  fit_function_->SetParameter(2, peak_height * 0.5);
  fit_function_->SetParameter(3, mu2_init);
  fit_function_->SetParameter(4, sigma_init);
  fit_function_->SetParameter(5, peak_height * 0.5);
  fit_function_->SetParameter(6, bkg_estimate);
  fit_function_->SetParameter(7, 0);

  if (use_flat_background_) {
    fit_function_->FixParameter(7, 0);
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    Double_t fitted_mu1 = fit_function_->GetParameter(0);
    Double_t fitted_mu2 = fit_function_->GetParameter(3);

    if (fitted_mu1 > fitted_mu2) {
      std::cout << "Warning: Fitted mu1 (" << fitted_mu1 << ") > mu2 ("
                << fitted_mu2 << "), swapping peaks" << std::endl;
      SwapDoublePeakStandardParameters();
    }

    std::cout << "Double peak standard fit converged successfully" << std::endl;
    std::cout << "Chi2/ndf = " << fit_result->Chi2() / fit_result->Ndf()
              << std::endl;

    PlotFitDoublePeakStandard(input_name, peak_name);

    results.peak1.mu = fit_function_->GetParameter(0);
    results.peak1.mu_error = fit_function_->GetParError(0);
    results.peak1.sigma = fit_function_->GetParameter(1);
    results.peak1.sigma_error = fit_function_->GetParError(1);
    results.peak1.gaus_amplitude = fit_function_->GetParameter(2);
    results.peak1.gaus_amplitude_error = fit_function_->GetParError(2);
    results.peak1.bkg_const = fit_function_->GetParameter(6);
    results.peak1.bkg_const_error = fit_function_->GetParError(6);
    results.peak1.bkg_slope = fit_function_->GetParameter(7);
    results.peak1.bkg_slope_error = fit_function_->GetParError(7);
    results.peak1.reduced_chi2 = -1;

    results.peak2.mu = fit_function_->GetParameter(3);
    results.peak2.mu_error = fit_function_->GetParError(3);
    results.peak2.sigma = fit_function_->GetParameter(4);
    results.peak2.sigma_error = fit_function_->GetParError(4);
    results.peak2.gaus_amplitude = fit_function_->GetParameter(5);
    results.peak2.gaus_amplitude_error = fit_function_->GetParError(5);
    results.peak2.bkg_const = fit_function_->GetParameter(6);
    results.peak2.bkg_const_error = fit_function_->GetParError(6);
    results.peak2.bkg_slope = fit_function_->GetParameter(7);
    results.peak2.bkg_slope_error = fit_function_->GetParError(7);
    results.peak2.reduced_chi2 = -1;

    results.reduced_chi2 = fit_result->Chi2() / fit_result->Ndf();
  } else {
    std::cout << "ERROR: Double peak standard fit failed to converge"
              << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
  }

  return results;
}

FitResultDoublePeakDetailed
FittingUtils::FitDoublePeakDetailed(const TString input_name,
                                    const TString peak_name, Double_t mu1_init,
                                    Double_t mu2_init) {
  FitResultDoublePeakDetailed results;

  if (mu1_init > mu2_init) {
    std::cout << "Warning: mu1_init > mu2_init, swapping initial values"
              << std::endl;
    Double_t temp = mu1_init;
    mu1_init = mu2_init;
    mu2_init = temp;
  }

  fit_function_ =
      new TF1("DoublePeakDetailed", &FittingFunctions::DoublePeakDetailed,
              fit_range_low_, fit_range_high_, 18);

  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "StepAmplitude1");
  fit_function_->SetParName(4, "LowTailAmplitude1");
  fit_function_->SetParName(5, "LowTailSlope1");
  fit_function_->SetParName(6, "HighTailAmplitude1");
  fit_function_->SetParName(7, "HighTailSlope1");
  fit_function_->SetParName(8, "Mu2");
  fit_function_->SetParName(9, "Sigma2");
  fit_function_->SetParName(10, "GausAmplitude2");
  fit_function_->SetParName(11, "StepAmplitude2");
  fit_function_->SetParName(12, "LowTailAmplitude2");
  fit_function_->SetParName(13, "LowTailSlope2");
  fit_function_->SetParName(14, "HighTailAmplitude2");
  fit_function_->SetParName(15, "HighTailSlope2");
  fit_function_->SetParName(16, "BkgConst");
  fit_function_->SetParName(17, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  fit_function_->SetParLimits(0, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(1, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(2, 0, peak_height * 1.5);
  fit_function_->SetParLimits(3, 0, peak_height * 0.5);
  fit_function_->SetParLimits(4, 0, peak_height);
  fit_function_->SetParLimits(5, 1, 10);
  fit_function_->SetParLimits(6, 0, peak_height);
  fit_function_->SetParLimits(7, 1, 10);

  fit_function_->SetParLimits(8, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(9, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(10, 0, peak_height * 1.5);
  fit_function_->SetParLimits(11, 0, peak_height * 0.5);
  fit_function_->SetParLimits(12, 0, peak_height);
  fit_function_->SetParLimits(13, 1, 10);
  fit_function_->SetParLimits(14, 0, peak_height);
  fit_function_->SetParLimits(15, 1, 10);

  fit_function_->SetParLimits(16, 0, peak_height * 0.5);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(17, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  } else {
    fit_function_->SetParLimits(17, 0, 0);
  }

  fit_function_->SetParameter(0, mu1_init);
  fit_function_->SetParameter(1, sigma_init);
  fit_function_->SetParameter(2, peak_height * 0.5);
  fit_function_->FixParameter(3, 0);
  fit_function_->FixParameter(4, 0);
  fit_function_->FixParameter(5, 1);
  fit_function_->FixParameter(6, 0);
  fit_function_->FixParameter(7, 1);

  fit_function_->SetParameter(8, mu2_init);
  fit_function_->SetParameter(9, sigma_init);
  fit_function_->SetParameter(10, peak_height * 0.5);
  fit_function_->FixParameter(11, 0);
  fit_function_->FixParameter(12, 0);
  fit_function_->FixParameter(13, 1);
  fit_function_->FixParameter(14, 0);
  fit_function_->FixParameter(15, 1);

  fit_function_->SetParameter(16, bkg_estimate);
  fit_function_->SetParameter(17, 0);

  if (use_flat_background_) {
    fit_function_->FixParameter(17, 0);
  }

  std::cout << "Detailed fit configuration:" << std::endl;
  std::cout << "  Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "  Low tail: " << (use_low_tail_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "  High tail: " << (use_high_tail_ ? "ENABLED" : "DISABLED")
            << std::endl;

  TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

  if (!initial_fit.Get() || !initial_fit->IsValid()) {
    std::cout << "ERROR: Initial double peak fit failed" << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
    return results;
  }

  Double_t gaus_amp1 = TMath::Abs(fit_function_->GetParameter(2));
  Double_t gaus_amp2 = TMath::Abs(fit_function_->GetParameter(10));

  std::vector<Double_t> best_params(fit_function_->GetNpar());
  std::vector<Double_t> best_errors(fit_function_->GetNpar());
  for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
    best_params[i] = fit_function_->GetParameter(i);
    best_errors[i] = fit_function_->GetParError(i);
  }

  Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
  std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

  if (use_step_) {
    std::cout << "Testing step function for peak1..." << std::endl;
    fit_function_->ReleaseParameter(3);
    fit_function_->SetParLimits(3, 0, peak_height);
    fit_function_->SetParameter(3, gaus_amp1);

    TFitResultPtr step1_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (step1_fit.Get() && step1_fit->IsValid()) {
      Double_t chi2_with_step1 = step1_fit->Chi2() / step1_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_step1 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_step1 < best_chi2) {
        std::cout << "Step1 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_step1;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Step1 REJECTED" << std::endl;
        fit_function_->FixParameter(3, 0);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Step1 fit FAILED" << std::endl;
      fit_function_->FixParameter(3, 0);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }

    std::cout << "Testing step function for peak2..." << std::endl;
    fit_function_->ReleaseParameter(11);
    fit_function_->SetParLimits(11, 0, peak_height);
    fit_function_->SetParameter(11, gaus_amp2);

    TFitResultPtr step2_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (step2_fit.Get() && step2_fit->IsValid()) {
      Double_t chi2_with_step2 = step2_fit->Chi2() / step2_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_step2 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_step2 < best_chi2) {
        std::cout << "Step2 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_step2;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Step2 REJECTED" << std::endl;
        fit_function_->FixParameter(11, 0);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Step2 fit FAILED" << std::endl;
      fit_function_->FixParameter(11, 0);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_low_tail_) {
    std::cout << "Testing low tail for peak1..." << std::endl;
    fit_function_->ReleaseParameter(4);
    fit_function_->ReleaseParameter(5);
    fit_function_->SetParLimits(4, 0, peak_height * 1.2);
    fit_function_->SetParLimits(5, 1, 10);
    Double_t tail_amp_init = TMath::Min(gaus_amp1 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(4, tail_amp_init);
    fit_function_->SetParameter(5, 1);

    TFitResultPtr lowtail1_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (lowtail1_fit.Get() && lowtail1_fit->IsValid()) {
      Double_t chi2_with_lowtail1 = lowtail1_fit->Chi2() / lowtail1_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_lowtail1 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_lowtail1 < best_chi2) {
        std::cout << "Low tail1 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_lowtail1;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Low tail1 REJECTED" << std::endl;
        fit_function_->FixParameter(4, 0);
        fit_function_->FixParameter(5, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Low tail1 fit FAILED" << std::endl;
      fit_function_->FixParameter(4, 0);
      fit_function_->FixParameter(5, 1);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }

    std::cout << "Testing low tail for peak2..." << std::endl;
    fit_function_->ReleaseParameter(12);
    fit_function_->ReleaseParameter(13);
    fit_function_->SetParLimits(12, 0, peak_height * 1.2);
    fit_function_->SetParLimits(13, 1, 10);
    tail_amp_init = TMath::Min(gaus_amp2 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(12, tail_amp_init);
    fit_function_->SetParameter(13, 1);

    TFitResultPtr lowtail2_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (lowtail2_fit.Get() && lowtail2_fit->IsValid()) {
      Double_t chi2_with_lowtail2 = lowtail2_fit->Chi2() / lowtail2_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_lowtail2 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_lowtail2 < best_chi2) {
        std::cout << "Low tail2 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_lowtail2;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Low tail2 REJECTED" << std::endl;
        fit_function_->FixParameter(12, 0);
        fit_function_->FixParameter(13, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Low tail2 fit FAILED" << std::endl;
      fit_function_->FixParameter(12, 0);
      fit_function_->FixParameter(13, 1);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_high_tail_) {
    std::cout << "Testing high tail for peak1..." << std::endl;
    fit_function_->ReleaseParameter(6);
    fit_function_->ReleaseParameter(7);
    fit_function_->SetParLimits(6, 0, peak_height * 1.2);
    fit_function_->SetParLimits(7, 1, 10);
    Double_t tail_amp_init = TMath::Min(gaus_amp1 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(6, tail_amp_init);
    fit_function_->SetParameter(7, 1);

    TFitResultPtr hightail1_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (hightail1_fit.Get() && hightail1_fit->IsValid()) {
      Double_t chi2_with_hightail1 =
          hightail1_fit->Chi2() / hightail1_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_hightail1 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_hightail1 < best_chi2) {
        std::cout << "High tail1 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_hightail1;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "High tail1 REJECTED" << std::endl;
        fit_function_->FixParameter(6, 0);
        fit_function_->FixParameter(7, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "High tail1 fit FAILED" << std::endl;
      fit_function_->FixParameter(6, 0);
      fit_function_->FixParameter(7, 1);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }

    std::cout << "Testing high tail for peak2..." << std::endl;
    fit_function_->ReleaseParameter(14);
    fit_function_->ReleaseParameter(15);
    fit_function_->SetParLimits(14, 0, peak_height * 1.2);
    fit_function_->SetParLimits(15, 1, 10);
    tail_amp_init = TMath::Min(gaus_amp2 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(14, tail_amp_init);
    fit_function_->SetParameter(15, 1);

    TFitResultPtr hightail2_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (hightail2_fit.Get() && hightail2_fit->IsValid()) {
      Double_t chi2_with_hightail2 =
          hightail2_fit->Chi2() / hightail2_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_hightail2 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_hightail2 < best_chi2) {
        std::cout << "High tail2 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_hightail2;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "High tail2 REJECTED" << std::endl;
        fit_function_->FixParameter(14, 0);
        fit_function_->FixParameter(15, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "High tail2 fit FAILED" << std::endl;
      fit_function_->FixParameter(14, 0);
      fit_function_->FixParameter(15, 1);
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

  if (use_flat_background_) {
    fit_function_->FixParameter(17, 0);
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMRBENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    Double_t fitted_mu1 = fit_function_->GetParameter(0);
    Double_t fitted_mu2 = fit_function_->GetParameter(8);

    if (fitted_mu1 > fitted_mu2) {
      std::cout << "Warning: Fitted mu1 (" << fitted_mu1 << ") > mu2 ("
                << fitted_mu2
                << "), swapping ALL peak parameters including tails and steps"
                << std::endl;
      SwapDoublePeakDetailedParameters();
    }

    Double_t final_chi2 = fit_result->Chi2() / fit_result->Ndf();
    std::cout << "Double peak detailed fit converged successfully" << std::endl;
    std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;

    PlotFitDoublePeakDetailed(input_name, peak_name);

    results.peak1.mu = fit_function_->GetParameter(0);
    results.peak1.mu_error = fit_function_->GetParError(0);
    results.peak1.sigma = fit_function_->GetParameter(1);
    results.peak1.sigma_error = fit_function_->GetParError(1);
    results.peak1.gaus_amplitude = fit_function_->GetParameter(2);
    results.peak1.gaus_amplitude_error = fit_function_->GetParError(2);
    results.peak1.step_amplitude = fit_function_->GetParameter(3);
    results.peak1.step_amplitude_error = fit_function_->GetParError(3);
    results.peak1.low_tail_amplitude = fit_function_->GetParameter(4);
    results.peak1.low_tail_amplitude_error = fit_function_->GetParError(4);
    results.peak1.low_tail_range = fit_function_->GetParameter(5);
    results.peak1.low_tail_range_error = fit_function_->GetParError(5);
    results.peak1.high_tail_amplitude = fit_function_->GetParameter(6);
    results.peak1.high_tail_amplitude_error = fit_function_->GetParError(6);
    results.peak1.high_tail_range = fit_function_->GetParameter(7);
    results.peak1.high_tail_range_error = fit_function_->GetParError(7);
    results.peak1.bkg_const = fit_function_->GetParameter(16);
    results.peak1.bkg_const_error = fit_function_->GetParError(16);
    results.peak1.bkg_slope = fit_function_->GetParameter(17);
    results.peak1.bkg_slope_error = fit_function_->GetParError(17);
    results.peak1.reduced_chi2 = -1;

    results.peak2.mu = fit_function_->GetParameter(8);
    results.peak2.mu_error = fit_function_->GetParError(8);
    results.peak2.sigma = fit_function_->GetParameter(9);
    results.peak2.sigma_error = fit_function_->GetParError(9);
    results.peak2.gaus_amplitude = fit_function_->GetParameter(10);
    results.peak2.gaus_amplitude_error = fit_function_->GetParError(10);
    results.peak2.step_amplitude = fit_function_->GetParameter(11);
    results.peak2.step_amplitude_error = fit_function_->GetParError(11);
    results.peak2.low_tail_amplitude = fit_function_->GetParameter(12);
    results.peak2.low_tail_amplitude_error = fit_function_->GetParError(12);
    results.peak2.low_tail_range = fit_function_->GetParameter(13);
    results.peak2.low_tail_range_error = fit_function_->GetParError(13);
    results.peak2.high_tail_amplitude = fit_function_->GetParameter(14);
    results.peak2.high_tail_amplitude_error = fit_function_->GetParError(14);
    results.peak2.high_tail_range = fit_function_->GetParameter(15);
    results.peak2.high_tail_range_error = fit_function_->GetParError(15);
    results.peak2.bkg_const = fit_function_->GetParameter(16);
    results.peak2.bkg_const_error = fit_function_->GetParError(16);
    results.peak2.bkg_slope = fit_function_->GetParameter(17);
    results.peak2.bkg_slope_error = fit_function_->GetParError(17);
    results.peak2.reduced_chi2 = -1;

    results.reduced_chi2 = final_chi2;
  } else {
    std::cout << "ERROR: Double peak detailed fit failed to converge"
              << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
  }

  return results;
}

FitResultDoublePeakStandard FittingUtils::FitDoublePeakStandard(
    const TString input_name, const TString peak_name,
    const FitResultStandard &constrained_peak, Double_t mu2_init) {
  FitResultDoublePeakStandard results;

  fit_function_ =
      new TF1("DoublePeakStandard", &FittingFunctions::DoublePeakStandard,
              fit_range_low_, fit_range_high_, 8);

  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "Mu2");
  fit_function_->SetParName(4, "Sigma2");
  fit_function_->SetParName(5, "GausAmplitude2");
  fit_function_->SetParName(6, "BkgConst");
  fit_function_->SetParName(7, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  fit_function_->FixParameter(0, constrained_peak.mu);
  fit_function_->FixParameter(1, constrained_peak.sigma);
  fit_function_->SetParLimits(2, 0, peak_height * 1.5);
  fit_function_->SetParameter(2, constrained_peak.gaus_amplitude);

  fit_function_->SetParLimits(3, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(4, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(5, 0, peak_height * 1.5);
  fit_function_->SetParameter(3, mu2_init);
  fit_function_->SetParameter(4, sigma_init);
  fit_function_->SetParameter(5, peak_height * 0.5);

  fit_function_->SetParLimits(6, 0, peak_height * 0.5);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(7, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  } else {
    fit_function_->SetParLimits(7, 0, 0);
  }
  fit_function_->SetParameter(6, bkg_estimate);
  fit_function_->SetParameter(7, 0);

  if (use_flat_background_) {
    fit_function_->FixParameter(7, 0);
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    std::cout << "Double peak standard fit (with constrained peak1) converged "
                 "successfully"
              << std::endl;
    std::cout << "Chi2/ndf = " << fit_result->Chi2() / fit_result->Ndf()
              << std::endl;

    PlotFitDoublePeakStandard(input_name, peak_name);

    results.peak1.mu = fit_function_->GetParameter(0);
    results.peak1.mu_error = fit_function_->GetParError(0);
    results.peak1.sigma = fit_function_->GetParameter(1);
    results.peak1.sigma_error = fit_function_->GetParError(1);
    results.peak1.gaus_amplitude = fit_function_->GetParameter(2);
    results.peak1.gaus_amplitude_error = fit_function_->GetParError(2);
    results.peak1.bkg_const = fit_function_->GetParameter(6);
    results.peak1.bkg_const_error = fit_function_->GetParError(6);
    results.peak1.bkg_slope = fit_function_->GetParameter(7);
    results.peak1.bkg_slope_error = fit_function_->GetParError(7);
    results.peak1.reduced_chi2 = -1;

    results.peak2.mu = fit_function_->GetParameter(3);
    results.peak2.mu_error = fit_function_->GetParError(3);
    results.peak2.sigma = fit_function_->GetParameter(4);
    results.peak2.sigma_error = fit_function_->GetParError(4);
    results.peak2.gaus_amplitude = fit_function_->GetParameter(5);
    results.peak2.gaus_amplitude_error = fit_function_->GetParError(5);
    results.peak2.bkg_const = fit_function_->GetParameter(6);
    results.peak2.bkg_const_error = fit_function_->GetParError(6);
    results.peak2.bkg_slope = fit_function_->GetParameter(7);
    results.peak2.bkg_slope_error = fit_function_->GetParError(7);
    results.peak2.reduced_chi2 = -1;

    results.reduced_chi2 = fit_result->Chi2() / fit_result->Ndf();
  } else {
    std::cout << "ERROR: Double peak standard fit (with constrained peak1) "
                 "failed to converge"
              << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
  }

  return results;
}

FitResultDoublePeakDetailed FittingUtils::FitDoublePeakDetailed(
    const TString input_name, const TString peak_name,
    const FitResultDetailed &constrained_peak, Double_t mu2_init) {
  FitResultDoublePeakDetailed results;

  fit_function_ =
      new TF1("DoublePeakDetailed", &FittingFunctions::DoublePeakDetailed,
              fit_range_low_, fit_range_high_, 18);

  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "StepAmplitude1");
  fit_function_->SetParName(4, "LowTailAmplitude1");
  fit_function_->SetParName(5, "LowTailSlope1");
  fit_function_->SetParName(6, "HighTailAmplitude1");
  fit_function_->SetParName(7, "HighTailSlope1");
  fit_function_->SetParName(8, "Mu2");
  fit_function_->SetParName(9, "Sigma2");
  fit_function_->SetParName(10, "GausAmplitude2");
  fit_function_->SetParName(11, "StepAmplitude2");
  fit_function_->SetParName(12, "LowTailAmplitude2");
  fit_function_->SetParName(13, "LowTailSlope2");
  fit_function_->SetParName(14, "HighTailAmplitude2");
  fit_function_->SetParName(15, "HighTailSlope2");
  fit_function_->SetParName(16, "BkgConst");
  fit_function_->SetParName(17, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  fit_function_->FixParameter(0, constrained_peak.mu);
  fit_function_->FixParameter(1, constrained_peak.sigma);
  fit_function_->SetParLimits(2, 0, peak_height * 1.5);
  fit_function_->SetParameter(2, constrained_peak.gaus_amplitude);
  fit_function_->SetParameter(3, constrained_peak.step_amplitude);
  fit_function_->SetParLimits(3, 0, constrained_peak.step_amplitude * 1.1);
  fit_function_->SetParameter(4, constrained_peak.low_tail_amplitude);
  fit_function_->FixParameter(5, constrained_peak.low_tail_range);
  fit_function_->SetParameter(6, constrained_peak.high_tail_amplitude);
  fit_function_->FixParameter(7, constrained_peak.high_tail_range);

  fit_function_->SetParLimits(8, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(9, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(10, 0, peak_height * 1.5);
  fit_function_->SetParLimits(11, 0, peak_height * 0.5);
  fit_function_->SetParLimits(12, 0, peak_height);
  fit_function_->SetParLimits(13, 1, 10);
  fit_function_->SetParLimits(14, 0, peak_height);
  fit_function_->SetParLimits(15, 1, 10);

  fit_function_->SetParameter(8, mu2_init);
  fit_function_->SetParameter(9, sigma_init);
  fit_function_->SetParameter(10, peak_height * 0.5);
  fit_function_->FixParameter(11, 0);
  fit_function_->FixParameter(12, 0);
  fit_function_->FixParameter(13, 1);
  fit_function_->FixParameter(14, 0);
  fit_function_->FixParameter(15, 1);

  fit_function_->SetParLimits(16, 0, peak_height * 0.5);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(17, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  } else {
    fit_function_->SetParLimits(17, 0, 0);
  }

  fit_function_->SetParameter(16, bkg_estimate);
  fit_function_->SetParameter(17, 0);

  if (use_flat_background_) {
    fit_function_->FixParameter(17, 0);
  }

  std::cout << "Detailed fit configuration (with constrained peak1):"
            << std::endl;
  std::cout << "  Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "  Low tail: " << (use_low_tail_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "  High tail: " << (use_high_tail_ ? "ENABLED" : "DISABLED")
            << std::endl;

  TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

  if (!initial_fit.Get() || !initial_fit->IsValid()) {
    std::cout << "ERROR: Initial double peak fit (with constrained peak1) "
                 "failed"
              << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
    return results;
  }

  Double_t gaus_amp2 = TMath::Abs(fit_function_->GetParameter(10));

  std::vector<Double_t> best_params(fit_function_->GetNpar());
  std::vector<Double_t> best_errors(fit_function_->GetNpar());
  for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
    best_params[i] = fit_function_->GetParameter(i);
    best_errors[i] = fit_function_->GetParError(i);
  }

  Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
  std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

  if (use_step_) {
    std::cout << "Testing step function for peak2..." << std::endl;
    fit_function_->ReleaseParameter(11);
    fit_function_->SetParLimits(11, 0, peak_height);
    fit_function_->SetParameter(11, gaus_amp2);

    TFitResultPtr step2_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (step2_fit.Get() && step2_fit->IsValid()) {
      Double_t chi2_with_step2 = step2_fit->Chi2() / step2_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_step2 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_step2 < best_chi2) {
        std::cout << "Step2 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_step2;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Step2 REJECTED" << std::endl;
        fit_function_->FixParameter(11, 0);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Step2 fit FAILED" << std::endl;
      fit_function_->FixParameter(11, 0);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_low_tail_) {
    std::cout << "Testing low tail for peak2..." << std::endl;
    fit_function_->ReleaseParameter(12);
    fit_function_->ReleaseParameter(13);
    fit_function_->SetParLimits(12, 0, peak_height * 1.2);
    fit_function_->SetParLimits(13, 1, 10);
    Double_t tail_amp_init = TMath::Min(gaus_amp2 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(12, tail_amp_init);
    fit_function_->SetParameter(13, 1);

    TFitResultPtr lowtail2_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (lowtail2_fit.Get() && lowtail2_fit->IsValid()) {
      Double_t chi2_with_lowtail2 = lowtail2_fit->Chi2() / lowtail2_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_lowtail2 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_lowtail2 < best_chi2) {
        std::cout << "Low tail2 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_lowtail2;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Low tail2 REJECTED" << std::endl;
        fit_function_->FixParameter(12, 0);
        fit_function_->FixParameter(13, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Low tail2 fit FAILED" << std::endl;
      fit_function_->FixParameter(12, 0);
      fit_function_->FixParameter(13, 1);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_high_tail_) {
    std::cout << "Testing high tail for peak2..." << std::endl;
    fit_function_->ReleaseParameter(14);
    fit_function_->ReleaseParameter(15);
    fit_function_->SetParLimits(14, 0, peak_height * 1.2);
    fit_function_->SetParLimits(15, 1, 10);
    Double_t tail_amp_init = TMath::Min(gaus_amp2 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(14, tail_amp_init);
    fit_function_->SetParameter(15, 1);

    TFitResultPtr hightail2_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (hightail2_fit.Get() && hightail2_fit->IsValid()) {
      Double_t chi2_with_hightail2 =
          hightail2_fit->Chi2() / hightail2_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_hightail2 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_hightail2 < best_chi2) {
        std::cout << "High tail2 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_hightail2;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "High tail2 REJECTED" << std::endl;
        fit_function_->FixParameter(14, 0);
        fit_function_->FixParameter(15, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "High tail2 fit FAILED" << std::endl;
      fit_function_->FixParameter(14, 0);
      fit_function_->FixParameter(15, 1);
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

  if (use_flat_background_) {
    fit_function_->FixParameter(17, 0);
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMRBENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    Double_t final_chi2 = fit_result->Chi2() / fit_result->Ndf();
    std::cout << "Double peak detailed fit (with constrained peak1) converged "
                 "successfully"
              << std::endl;
    std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;

    PlotFitDoublePeakDetailed(input_name, peak_name);

    results.peak1.mu = fit_function_->GetParameter(0);
    results.peak1.mu_error = fit_function_->GetParError(0);
    results.peak1.sigma = fit_function_->GetParameter(1);
    results.peak1.sigma_error = fit_function_->GetParError(1);
    results.peak1.gaus_amplitude = fit_function_->GetParameter(2);
    results.peak1.gaus_amplitude_error = fit_function_->GetParError(2);
    results.peak1.step_amplitude = fit_function_->GetParameter(3);
    results.peak1.step_amplitude_error = fit_function_->GetParError(3);
    results.peak1.low_tail_amplitude = fit_function_->GetParameter(4);
    results.peak1.low_tail_amplitude_error = fit_function_->GetParError(4);
    results.peak1.low_tail_range = fit_function_->GetParameter(5);
    results.peak1.low_tail_range_error = fit_function_->GetParError(5);
    results.peak1.high_tail_amplitude = fit_function_->GetParameter(6);
    results.peak1.high_tail_amplitude_error = fit_function_->GetParError(6);
    results.peak1.high_tail_range = fit_function_->GetParameter(7);
    results.peak1.high_tail_range_error = fit_function_->GetParError(7);
    results.peak1.bkg_const = fit_function_->GetParameter(16);
    results.peak1.bkg_const_error = fit_function_->GetParError(16);
    results.peak1.bkg_slope = fit_function_->GetParameter(17);
    results.peak1.bkg_slope_error = fit_function_->GetParError(17);
    results.peak1.reduced_chi2 = -1;

    results.peak2.mu = fit_function_->GetParameter(8);
    results.peak2.mu_error = fit_function_->GetParError(8);
    results.peak2.sigma = fit_function_->GetParameter(9);
    results.peak2.sigma_error = fit_function_->GetParError(9);
    results.peak2.gaus_amplitude = fit_function_->GetParameter(10);
    results.peak2.gaus_amplitude_error = fit_function_->GetParError(10);
    results.peak2.step_amplitude = fit_function_->GetParameter(11);
    results.peak2.step_amplitude_error = fit_function_->GetParError(11);
    results.peak2.low_tail_amplitude = fit_function_->GetParameter(12);
    results.peak2.low_tail_amplitude_error = fit_function_->GetParError(12);
    results.peak2.low_tail_range = fit_function_->GetParameter(13);
    results.peak2.low_tail_range_error = fit_function_->GetParError(13);
    results.peak2.high_tail_amplitude = fit_function_->GetParameter(14);
    results.peak2.high_tail_amplitude_error = fit_function_->GetParError(14);
    results.peak2.high_tail_range = fit_function_->GetParameter(15);
    results.peak2.high_tail_range_error = fit_function_->GetParError(15);
    results.peak2.bkg_const = fit_function_->GetParameter(16);
    results.peak2.bkg_const_error = fit_function_->GetParError(16);
    results.peak2.bkg_slope = fit_function_->GetParameter(17);
    results.peak2.bkg_slope_error = fit_function_->GetParError(17);
    results.peak2.reduced_chi2 = -1;

    results.reduced_chi2 = final_chi2;
  } else {
    std::cout << "ERROR: Double peak detailed fit (with constrained peak1) "
                 "failed to converge"
              << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
  }

  return results;
}

FitResultTriplePeakStandard FittingUtils::FitTriplePeakStandard(
    const TString input_name, const TString peak_name,
    const FitResultDoublePeakStandard &constrained_peaks, Double_t mu3_init) {
  FitResultTriplePeakStandard results;

  fit_function_ =
      new TF1("TriplePeakStandard", &FittingFunctions::TriplePeakStandard,
              fit_range_low_, fit_range_high_, 11);

  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "Mu2");
  fit_function_->SetParName(4, "Sigma2");
  fit_function_->SetParName(5, "GausAmplitude2");
  fit_function_->SetParName(6, "Mu3");
  fit_function_->SetParName(7, "Sigma3");
  fit_function_->SetParName(8, "GausAmplitude3");
  fit_function_->SetParName(9, "BkgConst");
  fit_function_->SetParName(10, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  fit_function_->FixParameter(0, constrained_peaks.peak1.mu);
  fit_function_->FixParameter(1, constrained_peaks.peak1.sigma);
  fit_function_->SetParLimits(2, 0, peak_height * 1.5);
  fit_function_->SetParameter(2, constrained_peaks.peak1.gaus_amplitude);

  fit_function_->FixParameter(3, constrained_peaks.peak2.mu);
  fit_function_->FixParameter(4, constrained_peaks.peak2.sigma);
  fit_function_->SetParLimits(5, 0, peak_height * 1.5);
  fit_function_->SetParameter(5, constrained_peaks.peak2.gaus_amplitude);

  fit_function_->SetParLimits(6, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(7, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(8, 0, peak_height * 1.5);
  fit_function_->SetParameter(6, mu3_init);
  fit_function_->SetParameter(7, sigma_init);
  fit_function_->SetParameter(8, peak_height * 0.5);

  fit_function_->SetParLimits(9, 0, peak_height * 0.5);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(10, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  } else {
    fit_function_->SetParLimits(10, 0, 0);
  }
  fit_function_->SetParameter(9, bkg_estimate);
  fit_function_->SetParameter(10, 0);

  if (use_flat_background_) {
    fit_function_->FixParameter(10, 0);
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    std::cout << "Triple peak standard fit converged successfully" << std::endl;
    std::cout << "Chi2/ndf = " << fit_result->Chi2() / fit_result->Ndf()
              << std::endl;

    PlotFitTriplePeakStandard(input_name, peak_name);

    results.peak1.mu = fit_function_->GetParameter(0);
    results.peak1.mu_error = fit_function_->GetParError(0);
    results.peak1.sigma = fit_function_->GetParameter(1);
    results.peak1.sigma_error = fit_function_->GetParError(1);
    results.peak1.gaus_amplitude = fit_function_->GetParameter(2);
    results.peak1.gaus_amplitude_error = fit_function_->GetParError(2);
    results.peak1.bkg_const = fit_function_->GetParameter(9);
    results.peak1.bkg_const_error = fit_function_->GetParError(9);
    results.peak1.bkg_slope = fit_function_->GetParameter(10);
    results.peak1.bkg_slope_error = fit_function_->GetParError(10);
    results.peak1.reduced_chi2 = -1;

    results.peak2.mu = fit_function_->GetParameter(3);
    results.peak2.mu_error = fit_function_->GetParError(3);
    results.peak2.sigma = fit_function_->GetParameter(4);
    results.peak2.sigma_error = fit_function_->GetParError(4);
    results.peak2.gaus_amplitude = fit_function_->GetParameter(5);
    results.peak2.gaus_amplitude_error = fit_function_->GetParError(5);
    results.peak2.bkg_const = fit_function_->GetParameter(9);
    results.peak2.bkg_const_error = fit_function_->GetParError(9);
    results.peak2.bkg_slope = fit_function_->GetParameter(10);
    results.peak2.bkg_slope_error = fit_function_->GetParError(10);
    results.peak2.reduced_chi2 = -1;

    results.peak3.mu = fit_function_->GetParameter(6);
    results.peak3.mu_error = fit_function_->GetParError(6);
    results.peak3.sigma = fit_function_->GetParameter(7);
    results.peak3.sigma_error = fit_function_->GetParError(7);
    results.peak3.gaus_amplitude = fit_function_->GetParameter(8);
    results.peak3.gaus_amplitude_error = fit_function_->GetParError(8);
    results.peak3.bkg_const = fit_function_->GetParameter(9);
    results.peak3.bkg_const_error = fit_function_->GetParError(9);
    results.peak3.bkg_slope = fit_function_->GetParameter(10);
    results.peak3.bkg_slope_error = fit_function_->GetParError(10);
    results.peak3.reduced_chi2 = -1;

    results.reduced_chi2 = fit_result->Chi2() / fit_result->Ndf();

  } else {
    std::cout << "ERROR: Triple peak standard fit failed to converge"
              << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
    results.peak3.mu_error = -1;
  }

  return results;
}

FitResultTriplePeakDetailed FittingUtils::FitTriplePeakDetailed(
    const TString input_name, const TString peak_name,
    const FitResultDoublePeakDetailed &constrained_peaks, Double_t mu3_init) {
  FitResultTriplePeakDetailed results;

  fit_function_ =
      new TF1("TriplePeakDetailed", &FittingFunctions::TriplePeakDetailed,
              fit_range_low_, fit_range_high_, 26);

  fit_function_->SetParName(0, "Mu1");
  fit_function_->SetParName(1, "Sigma1");
  fit_function_->SetParName(2, "GausAmplitude1");
  fit_function_->SetParName(3, "StepAmplitude1");
  fit_function_->SetParName(4, "LowTailAmplitude1");
  fit_function_->SetParName(5, "LowTailSlope1");
  fit_function_->SetParName(6, "HighTailAmplitude1");
  fit_function_->SetParName(7, "HighTailSlope1");
  fit_function_->SetParName(8, "Mu2");
  fit_function_->SetParName(9, "Sigma2");
  fit_function_->SetParName(10, "GausAmplitude2");
  fit_function_->SetParName(11, "StepAmplitude2");
  fit_function_->SetParName(12, "LowTailAmplitude2");
  fit_function_->SetParName(13, "LowTailSlope2");
  fit_function_->SetParName(14, "HighTailAmplitude2");
  fit_function_->SetParName(15, "HighTailSlope2");
  fit_function_->SetParName(16, "Mu3");
  fit_function_->SetParName(17, "Sigma3");
  fit_function_->SetParName(18, "GausAmplitude3");
  fit_function_->SetParName(19, "StepAmplitude3");
  fit_function_->SetParName(20, "LowTailAmplitude3");
  fit_function_->SetParName(21, "LowTailSlope3");
  fit_function_->SetParName(22, "HighTailAmplitude3");
  fit_function_->SetParName(23, "HighTailSlope3");
  fit_function_->SetParName(24, "BkgConst");
  fit_function_->SetParName(25, "BkgSlope");

  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  fit_function_->FixParameter(0, constrained_peaks.peak1.mu);
  fit_function_->FixParameter(1, constrained_peaks.peak1.sigma);
  fit_function_->SetParLimits(2, 0, peak_height * 1.5);
  fit_function_->SetParameter(2, constrained_peaks.peak1.gaus_amplitude);
  fit_function_->SetParameter(3, constrained_peaks.peak1.step_amplitude);
  fit_function_->SetParLimits(3, 0,
                              constrained_peaks.peak1.step_amplitude * 1.1);
  fit_function_->SetParameter(4, constrained_peaks.peak1.low_tail_amplitude);
  fit_function_->FixParameter(5, constrained_peaks.peak1.low_tail_range);
  fit_function_->SetParameter(6, constrained_peaks.peak1.high_tail_amplitude);
  fit_function_->FixParameter(7, constrained_peaks.peak1.high_tail_range);

  fit_function_->FixParameter(8, constrained_peaks.peak2.mu);
  fit_function_->FixParameter(9, constrained_peaks.peak2.sigma);
  fit_function_->SetParLimits(10, 0, peak_height * 1.5);
  fit_function_->SetParameter(10, constrained_peaks.peak2.gaus_amplitude);
  fit_function_->SetParameter(11, constrained_peaks.peak2.step_amplitude);
  fit_function_->SetParLimits(11, 0,
                              constrained_peaks.peak2.step_amplitude * 1.1);
  fit_function_->SetParameter(12, constrained_peaks.peak2.low_tail_amplitude);
  fit_function_->FixParameter(13, constrained_peaks.peak2.low_tail_range);
  fit_function_->SetParameter(14, constrained_peaks.peak2.high_tail_amplitude);
  fit_function_->FixParameter(15, constrained_peaks.peak2.high_tail_range);

  fit_function_->SetParLimits(16, fit_range_low_, fit_range_high_);
  fit_function_->SetParLimits(17, range_width * 0.001, range_width * 0.5);
  fit_function_->SetParLimits(18, 0, peak_height * 1.5);
  fit_function_->SetParLimits(19, 0, peak_height * 0.5);
  fit_function_->SetParLimits(20, 0, peak_height);
  fit_function_->SetParLimits(21, 1, 10);
  fit_function_->SetParLimits(22, 0, peak_height);
  fit_function_->SetParLimits(23, 1, 10);

  fit_function_->SetParameter(16, mu3_init);
  fit_function_->SetParameter(17, sigma_init);
  fit_function_->SetParameter(18, peak_height * 0.5);
  fit_function_->FixParameter(19, 0);
  fit_function_->FixParameter(20, 0);
  fit_function_->FixParameter(21, 1);
  fit_function_->FixParameter(22, 0);
  fit_function_->FixParameter(23, 1);

  fit_function_->SetParLimits(24, 0, peak_height * 0.5);
  if (!use_flat_background_) {
    fit_function_->SetParLimits(25, -0.1 * bkg_estimate / range_width,
                                0.1 * bkg_estimate / range_width);
  } else {
    fit_function_->SetParLimits(25, 0, 0);
  }
  fit_function_->SetParameter(24, bkg_estimate);
  fit_function_->SetParameter(25, 0);

  if (use_flat_background_) {
    fit_function_->FixParameter(25, 0);
  }

  std::cout << "Detailed fit configuration:" << std::endl;
  std::cout << "  Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "  Low tail: " << (use_low_tail_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "  High tail: " << (use_high_tail_ ? "ENABLED" : "DISABLED")
            << std::endl;

  TFitResultPtr initial_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

  if (!initial_fit.Get() || !initial_fit->IsValid()) {
    std::cout << "ERROR: Initial triple peak fit failed" << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
    results.peak3.mu_error = -1;
    return results;
  }

  Double_t gaus_amp3 = TMath::Abs(fit_function_->GetParameter(18));

  std::vector<Double_t> best_params(fit_function_->GetNpar());
  std::vector<Double_t> best_errors(fit_function_->GetNpar());
  for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
    best_params[i] = fit_function_->GetParameter(i);
    best_errors[i] = fit_function_->GetParError(i);
  }

  Double_t best_chi2 = initial_fit->Chi2() / initial_fit->Ndf();
  std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;

  if (use_step_) {
    std::cout << "Testing step function for peak3..." << std::endl;
    fit_function_->ReleaseParameter(19);
    fit_function_->SetParLimits(19, 0, peak_height);
    fit_function_->SetParameter(19, gaus_amp3);

    TFitResultPtr step3_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (step3_fit.Get() && step3_fit->IsValid()) {
      Double_t chi2_with_step3 = step3_fit->Chi2() / step3_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_step3 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_step3 < best_chi2) {
        std::cout << "Step3 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_step3;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Step3 REJECTED" << std::endl;
        fit_function_->FixParameter(19, 0);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Step3 fit FAILED" << std::endl;
      fit_function_->FixParameter(19, 0);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_low_tail_) {
    std::cout << "Testing low tail for peak3..." << std::endl;
    fit_function_->ReleaseParameter(20);
    fit_function_->ReleaseParameter(21);
    fit_function_->SetParLimits(20, 0, peak_height * 1.2);
    fit_function_->SetParLimits(21, 1, 10);
    Double_t tail_amp_init = TMath::Min(gaus_amp3 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(20, tail_amp_init);
    fit_function_->SetParameter(21, 1);

    TFitResultPtr lowtail3_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (lowtail3_fit.Get() && lowtail3_fit->IsValid()) {
      Double_t chi2_with_lowtail3 = lowtail3_fit->Chi2() / lowtail3_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_lowtail3 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_lowtail3 < best_chi2) {
        std::cout << "Low tail3 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_lowtail3;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "Low tail3 REJECTED" << std::endl;
        fit_function_->FixParameter(20, 0);
        fit_function_->FixParameter(21, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "Low tail3 fit FAILED" << std::endl;
      fit_function_->FixParameter(20, 0);
      fit_function_->FixParameter(21, 1);
      for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
        fit_function_->SetParameter(i, best_params[i]);
        fit_function_->SetParError(i, best_errors[i]);
      }
    }
  }

  if (use_high_tail_) {
    std::cout << "Testing high tail for peak3..." << std::endl;
    fit_function_->ReleaseParameter(22);
    fit_function_->ReleaseParameter(23);
    fit_function_->SetParLimits(22, 0, peak_height * 1.2);
    fit_function_->SetParLimits(23, 1, 10);
    Double_t tail_amp_init = TMath::Min(gaus_amp3 * 0.15, peak_height * 0.25);
    fit_function_->SetParameter(22, tail_amp_init);
    fit_function_->SetParameter(23, 1);

    TFitResultPtr hightail3_fit = working_hist_->Fit(fit_function_, "LSMBNQ0R");

    if (hightail3_fit.Get() && hightail3_fit->IsValid()) {
      Double_t chi2_with_hightail3 =
          hightail3_fit->Chi2() / hightail3_fit->Ndf();
      std::cout << "Chi2/ndf: " << chi2_with_hightail3 << " vs " << best_chi2
                << std::endl;

      if (chi2_with_hightail3 < best_chi2) {
        std::cout << "High tail3 ACCEPTED" << std::endl;
        best_chi2 = chi2_with_hightail3;
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          best_params[i] = fit_function_->GetParameter(i);
          best_errors[i] = fit_function_->GetParError(i);
        }
      } else {
        std::cout << "High tail3 REJECTED" << std::endl;
        fit_function_->FixParameter(22, 0);
        fit_function_->FixParameter(23, 1);
        for (Int_t i = 0; i < fit_function_->GetNpar(); i++) {
          fit_function_->SetParameter(i, best_params[i]);
          fit_function_->SetParError(i, best_errors[i]);
        }
      }
    } else {
      std::cout << "High tail3 fit FAILED" << std::endl;
      fit_function_->FixParameter(22, 0);
      fit_function_->FixParameter(23, 1);
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

  if (use_flat_background_) {
    fit_function_->FixParameter(25, 0);
  }

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSMRBENR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    Double_t final_chi2 = fit_result->Chi2() / fit_result->Ndf();
    std::cout << "Triple peak detailed fit converged successfully" << std::endl;
    std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;

    PlotFitTriplePeakDetailed(input_name, peak_name);

    results.peak1.mu = fit_function_->GetParameter(0);
    results.peak1.mu_error = fit_function_->GetParError(0);
    results.peak1.sigma = fit_function_->GetParameter(1);
    results.peak1.sigma_error = fit_function_->GetParError(1);
    results.peak1.gaus_amplitude = fit_function_->GetParameter(2);
    results.peak1.gaus_amplitude_error = fit_function_->GetParError(2);
    results.peak1.step_amplitude = fit_function_->GetParameter(3);
    results.peak1.step_amplitude_error = fit_function_->GetParError(3);
    results.peak1.low_tail_amplitude = fit_function_->GetParameter(4);
    results.peak1.low_tail_amplitude_error = fit_function_->GetParError(4);
    results.peak1.low_tail_range = fit_function_->GetParameter(5);
    results.peak1.low_tail_range_error = fit_function_->GetParError(5);
    results.peak1.high_tail_amplitude = fit_function_->GetParameter(6);
    results.peak1.high_tail_amplitude_error = fit_function_->GetParError(6);
    results.peak1.high_tail_range = fit_function_->GetParameter(7);
    results.peak1.high_tail_range_error = fit_function_->GetParError(7);
    results.peak1.bkg_const = fit_function_->GetParameter(24);
    results.peak1.bkg_const_error = fit_function_->GetParError(24);
    results.peak1.bkg_slope = fit_function_->GetParameter(25);
    results.peak1.bkg_slope_error = fit_function_->GetParError(25);
    results.peak1.reduced_chi2 = -1;

    results.peak2.mu = fit_function_->GetParameter(8);
    results.peak2.mu_error = fit_function_->GetParError(8);
    results.peak2.sigma = fit_function_->GetParameter(9);
    results.peak2.sigma_error = fit_function_->GetParError(9);
    results.peak2.gaus_amplitude = fit_function_->GetParameter(10);
    results.peak2.gaus_amplitude_error = fit_function_->GetParError(10);
    results.peak2.step_amplitude = fit_function_->GetParameter(11);
    results.peak2.step_amplitude_error = fit_function_->GetParError(11);
    results.peak2.low_tail_amplitude = fit_function_->GetParameter(12);
    results.peak2.low_tail_amplitude_error = fit_function_->GetParError(12);
    results.peak2.low_tail_range = fit_function_->GetParameter(13);
    results.peak2.low_tail_range_error = fit_function_->GetParError(13);
    results.peak2.high_tail_amplitude = fit_function_->GetParameter(14);
    results.peak2.high_tail_amplitude_error = fit_function_->GetParError(14);
    results.peak2.high_tail_range = fit_function_->GetParameter(15);
    results.peak2.high_tail_range_error = fit_function_->GetParError(15);
    results.peak2.bkg_const = fit_function_->GetParameter(24);
    results.peak2.bkg_const_error = fit_function_->GetParError(24);
    results.peak2.bkg_slope = fit_function_->GetParameter(25);
    results.peak2.bkg_slope_error = fit_function_->GetParError(25);
    results.peak2.reduced_chi2 = -1;

    results.peak3.mu = fit_function_->GetParameter(16);
    results.peak3.mu_error = fit_function_->GetParError(16);
    results.peak3.sigma = fit_function_->GetParameter(17);
    results.peak3.sigma_error = fit_function_->GetParError(17);
    results.peak3.gaus_amplitude = fit_function_->GetParameter(18);
    results.peak3.gaus_amplitude_error = fit_function_->GetParError(18);
    results.peak3.step_amplitude = fit_function_->GetParameter(19);
    results.peak3.step_amplitude_error = fit_function_->GetParError(19);
    results.peak3.low_tail_amplitude = fit_function_->GetParameter(20);
    results.peak3.low_tail_amplitude_error = fit_function_->GetParError(20);
    results.peak3.low_tail_range = fit_function_->GetParameter(21);
    results.peak3.low_tail_range_error = fit_function_->GetParError(21);
    results.peak3.high_tail_amplitude = fit_function_->GetParameter(22);
    results.peak3.high_tail_amplitude_error = fit_function_->GetParError(22);
    results.peak3.high_tail_range = fit_function_->GetParameter(23);
    results.peak3.high_tail_range_error = fit_function_->GetParError(23);
    results.peak3.bkg_const = fit_function_->GetParameter(24);
    results.peak3.bkg_const_error = fit_function_->GetParError(24);
    results.peak3.bkg_slope = fit_function_->GetParameter(25);
    results.peak3.bkg_slope_error = fit_function_->GetParError(25);
    results.peak3.reduced_chi2 = -1;

    results.reduced_chi2 = final_chi2;
  } else {
    std::cout << "ERROR: Triple peak detailed fit failed to converge"
              << std::endl;
    std::cout << "Fit status: " << fit_result->Status() << std::endl;
    results.peak1.mu_error = -1;
    results.peak2.mu_error = -1;
    results.peak3.mu_error = -1;
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
  f_detailed->SetParLimits(2, 0, 1e7);
  f_detailed->SetParLimits(5, 0, 1e7);
  f_detailed->SetParLimits(6, 0, 1e7);
  f_detailed->SetParLimits(7, 0, 1e7);
  f_detailed->SetParLimits(8, 0, 1e7);
  f_detailed->SetParLimits(9, 0, 1e7);

  std::cout << "Custom fitting functions registered and available in FitPanel!"
            << std::endl;
}
