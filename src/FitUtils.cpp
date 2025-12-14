#include "FitUtils.hpp"

FitUtils::FitUtils()
    : working_tree_(nullptr), working_value_(0), fit_successful_(kFALSE),
      num_hist_bins_(0), max_hist_value_(0) {
  fit_function_ = new TF1("gaussian_plus_linear", "gaus(0)+pol1(3)");
  fit_function_->SetParameter(3, 0);
}

FitUtils::~FitUtils() {
  fit_function_ = nullptr;
  working_tree_ = nullptr;
}

Bool_t FitUtils::LoadProcessed(const TString input_name,
                               const TString branch_name) {
  TString input_filename = input_name + ".root";
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

Bool_t FitUtils::FitPeak(TCanvas *canvas, const TString input_name,
                         const TString formatted_branch_name) {
  canvas->SetLogy(kFALSE);
  TH1F *hist =
      new TH1F("", Form(";%s [a.u.]; Counts", formatted_branch_name.Data()),
               num_hist_bins_, 0, max_hist_value_);

  fit_function_->SetParName(0, "Amplitude");
  fit_function_->SetParName(1, "Mean");
  fit_function_->SetParName(2, "Sigma");
  fit_function_->SetParName(3, "Bkg_Const");
  fit_function_->SetParName(4, "Bkg_Slope");

  fit_function_->SetParLimits(0, 0, 1e5);
  fit_function_->SetParLimits(3, 0, 0);
  fit_function_->SetParLimits(4, 0, 0);

  TFitResultPtr fit_result = hist->Fit(fit_function_, "LSRN+");

  if (fit_result.Get() && fit_result->IsValid()) {
    fit_result->Draw();
    PlottingUtils::SaveFigure(canvas, input_name, kFALSE);
    fit_successful_ = kTRUE;
    return kTRUE;
  }

  fit_successful_ = kFALSE;
  return kFALSE;
}
