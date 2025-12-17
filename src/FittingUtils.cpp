#include "FittingUtils.hpp"

FittingUtils::FittingUtils()
    : working_tree_(nullptr), working_value_(0), num_hist_bins_(0),
      max_hist_value_(0) {
  fit_function_ = new TF1("gaussian_plus_linear", "gaus(0)+pol1(3)");
  fit_function_->SetParameter(3, 0);
}

FittingUtils::~FittingUtils() {
  fit_function_ = nullptr;
  working_hist_ = nullptr;
  working_tree_ = nullptr;
}

Bool_t FittingUtils::LoadProcessed(const TString input_name,
                                   const TString branch_name) {

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
  working_tree_ = static_cast<TTree *>(input_file->Get("features"));
  working_tree_->SetBranchAddress(branch_name, &working_value_);

  return kTRUE;
}

void FittingUtils::PlotFit(TCanvas *canvas, Int_t color,
                           const TString peak_name) {
  Float_t range_low = fit_range_low_ - 0.1 * fit_range_low_;
  Float_t range_up = fit_range_high_ + 0.1 * fit_range_high_;
  working_hist_->GetXaxis()->SetRangeUser(range_low, range_up);
  PlottingUtils::ConfigureAndDrawHistogram(working_hist_, color);
  fit_function_->Draw("same");
  TF1 *peak = new TF1("gaussian", "gaus", fit_range_low_, fit_range_high_);
  peak->SetParameter(0, fit_function_->GetParameter("Amplitude"));
  peak->SetParameter(1, fit_function_->GetParameter("Mu"));
  peak->SetParameter(2, fit_function_->GetParameter("Sigma"));
  peak->SetLineColor(kBlack);
  peak->Draw("same");
  TF1 *background =
      new TF1("background", "pol1", fit_range_low_, fit_range_high_);
  background->SetParameter(0, fit_function_->GetParameter("BkgConst"));
  background->SetParameter(1, fit_function_->GetParameter("BkgSlope"));
  background->SetLineColor(kGreen);
  background->Draw("same");
  PlottingUtils::SaveFigure(canvas, peak_name + ".png", kTRUE);
}

FitResult
FittingUtils::FitPeak(TCanvas *canvas, Int_t color, const TString peak_name,
                      const TString formatted_branch_name_with_units) {
  FitResult results;
  working_hist_ =
      new TH1F("", Form(";%s; Counts", formatted_branch_name_with_units.Data()),
               num_hist_bins_, 0, max_hist_value_);

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
  fit_function_->SetParLimits(1, 0, max_hist_value_);
  fit_function_->SetParLimits(2, 0, 0.1 * max_hist_value_);
  fit_function_->SetParLimits(3, 0, 1e6);

  TFitResultPtr fit_result = working_hist_->Fit(fit_function_, "LSRN+");

  if (fit_result.Get() && fit_result->IsValid()) {
    PlotFit(canvas, color, peak_name);
    results.mu = fit_function_->GetParameter("Mu");
    results.mu_error = fit_function_->GetParError("Mu");
    results.sigma = fit_function_->GetParameter("Sigma");
    results.sigma_error = fit_function_->GetParError("Sigma");
    return results;
  }

  return results;
}
