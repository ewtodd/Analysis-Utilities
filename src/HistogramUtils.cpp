#include "HistogramUtils.hpp"
#include <TMath.h>
#include <glob.h>
#include <iostream>
#include <string>

HistogramUtils::HistogramUtils() {}

HistogramUtils::~HistogramUtils() {
  for (auto &pair : calibrated_spectra_) {
    delete pair.second;
  }
  for (auto &pair : integral_spectra_) {
    delete pair.second;
  }
  for (auto &pair : pulse_height_spectra_) {
    delete pair.second;
  }
  for (auto &pair : charge_comparison_spectra_) {
    delete pair.second;
  }
}

void HistogramUtils::SetConfig(const HistogramConfig &config) {
  config_ = config;
}

void HistogramUtils::AddSource(Int_t source_id, const std::string &name) {
  source_names_[source_id] = name;
}

void HistogramUtils::CreateAllHistograms() {
  for (const auto &pair : source_names_) {
    CreateHistogramForSource(pair.first);
  }
}

void HistogramUtils::CreateHistogramForSource(Int_t source_id) {
  std::string source_name = source_names_[source_id];

  std::string clean_name = source_name;
  std::replace(clean_name.begin(), clean_name.end(), ' ', '_');
  std::replace(clean_name.begin(), clean_name.end(), '-', '_');
  std::replace(clean_name.begin(), clean_name.end(), '&', '_');

  std::string cal_hist_name = "h_calibrated_" + clean_name;
  std::string cal_hist_title =
      source_name +
      " Calibrated Energy Deposition;Deposited Energy [keVee];Counts / " +
      std::to_string(Int_t(config_.calibrated_bin_width)) + " keV";

  calibrated_spectra_[source_id] = new TH1F(
      cal_hist_name.c_str(), cal_hist_title.c_str(), config_.calibrated_bins,
      config_.calibrated_min, config_.calibrated_max);

  std::string integral_hist_name = "h_integral_" + clean_name;
  std::string integral_hist_title =
      source_name + " Pulse Integral;Pulse Integral [a.u.];Counts";

  integral_spectra_[source_id] = new TH1F(
      integral_hist_name.c_str(), integral_hist_title.c_str(),
      config_.integral_bins, config_.integral_min, config_.integral_max);

  std::string ph_hist_name = "h_pulse_height_" + clean_name;
  std::string ph_hist_title =
      source_name + " Pulse Height;Pulse Height [a.u.];Counts";

  pulse_height_spectra_[source_id] =
      new TH1F(ph_hist_name.c_str(), ph_hist_title.c_str(), config_.ph_bins,
               config_.ph_min, config_.ph_max);

  // Charge comparison spectrum
  std::string cc_hist_name = "h_charge_comparison_" + clean_name;
  std::string cc_hist_title =
      source_name + " Charge Comparison;CC Parameter;Counts";

  charge_comparison_spectra_[source_id] =
      new TH1F(cc_hist_name.c_str(), cc_hist_title.c_str(), config_.cc_bins,
               config_.cc_min, config_.cc_max);

  // Charge comparison spectrum
  std::string si_hist_name = "h_shape_indicator_" + clean_name;
  std::string si_hist_title =
      source_name + " Shape Indicator;SI Parameter;Counts";

  // Set ownership to avoid ROOT auto-deletion
  calibrated_spectra_[source_id]->SetDirectory(0);
  integral_spectra_[source_id]->SetDirectory(0);
  pulse_height_spectra_[source_id]->SetDirectory(0);
  charge_comparison_spectra_[source_id]->SetDirectory(0);
}

void HistogramUtils::FillFromTree(const std::string &filename,
                                  const std::string &treename) {
  TFile *file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open file " << filename << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get(treename.c_str()));
  if (!tree) {
    std::cout << "Error: Could not find tree " << treename << " in file "
              << filename << std::endl;
    file->Close();
    return;
  }

  if (!tree) {
    std::cout << "Error: Null tree pointer!" << std::endl;
    return;
  }

  // Set up branch addresses
  Float_t pulse_height, long_integral;
  Int_t source_id;
  Bool_t passes_cuts;

  tree->SetBranchAddress("pulse_height", &pulse_height);
  tree->SetBranchAddress("long_integral", &long_integral);
  tree->SetBranchAddress("source_id", &source_id);
  tree->SetBranchAddress("passes_cuts", &passes_cuts);

  Long64_t n_entries = tree->GetEntries();
  std::cout << "Filling histograms from " << n_entries << " entries..."
            << std::endl;

  for (Long64_t i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);

    if (!passes_cuts)
      continue;

    if (integral_spectra_.find(source_id) == integral_spectra_.end()) {
      std::cout << "Warning: No histograms found for source_id " << source_id
                << std::endl;
      continue;
    }

    integral_spectra_[source_id]->Fill(long_integral);
    pulse_height_spectra_[source_id]->Fill(pulse_height);
  }

  std::cout << "Histogram filling complete." << std::endl;

  file->Close();
}

TH1F *HistogramUtils::GetCalibratedSpectrum(Int_t source_id) const {
  auto it = calibrated_spectra_.find(source_id);
  return (it != calibrated_spectra_.end()) ? it->second : nullptr;
}

TH1F *HistogramUtils::GetIntegralSpectrum(Int_t source_id) const {
  auto it = integral_spectra_.find(source_id);
  return (it != integral_spectra_.end()) ? it->second : nullptr;
}

TH1F *HistogramUtils::GetPulseHeightSpectrum(Int_t source_id) const {
  auto it = pulse_height_spectra_.find(source_id);
  return (it != pulse_height_spectra_.end()) ? it->second : nullptr;
}

TH1F *HistogramUtils::GetChargeComparisonSpectrum(Int_t source_id) const {
  auto it = charge_comparison_spectra_.find(source_id);
  return (it != charge_comparison_spectra_.end()) ? it->second : nullptr;
}

Bool_t HistogramUtils::LoadFromFile(const std::string &filename) {
  TFile *file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open file " << filename << std::endl;
    return kFALSE;
  }

  std::cout << "Loading histograms from " << filename << "..." << std::endl;

  for (const auto &source_pair : source_names_) {
    Int_t source_id = source_pair.first;
    std::string source_name = source_pair.second;

    // Clean name for histogram lookup
    std::string clean_name = source_name;
    std::replace(clean_name.begin(), clean_name.end(), ' ', '_');
    std::replace(clean_name.begin(), clean_name.end(), '-', '_');
    std::replace(clean_name.begin(), clean_name.end(), '&', '_');

    // Try to load each histogram type
    std::string cal_hist_name = "calibrated_spectra/h_calibrated_" + clean_name;
    std::string integral_hist_name =
        "integral_spectra/h_integral_" + clean_name;
    std::string ph_hist_name =
        "pulse_height_spectra/h_pulse_height_" + clean_name;
    std::string cc_hist_name =
        "charge_comparison_spectra/h_charge_comparison_" + clean_name;
    std::string si_hist_name =
        "shape_indicator_spectra/h_shape_indicator_" + clean_name;

    TH1F *cal_hist = static_cast<TH1F *>(file->Get(cal_hist_name.c_str()));
    TH1F *integral_hist =
        static_cast<TH1F *>(file->Get(integral_hist_name.c_str()));
    TH1F *ph_hist = static_cast<TH1F *>(file->Get(ph_hist_name.c_str()));
    TH1F *cc_hist = static_cast<TH1F *>(file->Get(cc_hist_name.c_str()));
    TH1F *si_hist = static_cast<TH1F *>(file->Get(si_hist_name.c_str()));

    if (cal_hist) {
      cal_hist->SetDirectory(0);
      calibrated_spectra_[source_id] = cal_hist;
    }
    if (integral_hist) {
      integral_hist->SetDirectory(0);
      integral_spectra_[source_id] = integral_hist;
    }
    if (ph_hist) {
      ph_hist->SetDirectory(0);
      pulse_height_spectra_[source_id] = ph_hist;
    }
    if (cc_hist) {
      cc_hist->SetDirectory(0);
      charge_comparison_spectra_[source_id] = cc_hist;
    }
  }

  file->Close();
  std::cout << "Histograms loaded successfully." << std::endl;
  return kTRUE;
}

void HistogramUtils::SaveToFile(const std::string &filename) {
  TFile *file = new TFile(filename.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not create output file " << filename
              << std::endl;
    return;
  }

  std::cout << "Saving histograms to " << filename << "..." << std::endl;

  // Create directories for organization
  TDirectory *calibrated_dir = file->mkdir("calibrated_spectra");
  TDirectory *integral_dir = file->mkdir("integral_spectra");
  TDirectory *pulse_height_dir = file->mkdir("pulse_height_spectra");
  TDirectory *cc_dir = file->mkdir("charge_comparison_spectra");
  TDirectory *si_dir = file->mkdir("shape_indicator_spectra");

  // Save light output spectra
  calibrated_dir->cd();
  for (const auto &pair : calibrated_spectra_) {
    pair.second->Write("", TObject::kOverwrite);
  }

  // Save integral spectra
  integral_dir->cd();
  for (const auto &pair : integral_spectra_) {
    pair.second->Write("", TObject::kOverwrite);
  }

  // Save pulse height spectra
  pulse_height_dir->cd();
  for (const auto &pair : pulse_height_spectra_) {
    pair.second->Write("", TObject::kOverwrite);
  }

  // Save charge comparison spectra
  cc_dir->cd();
  for (const auto &pair : charge_comparison_spectra_) {
    pair.second->Write("", TObject::kOverwrite);
  }

  file->Close();
  delete file;

  std::cout << "Histograms saved successfully." << std::endl;
}
