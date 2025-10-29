#include "HistogramUtils.h"
#include <TMath.h>
#include <fstream> // For std::ifstream
#include <glob.h>  // For glob, glob_t, GLOB_TILDE, globfree
#include <iostream>
#include <sstream> // For std::istringstream
#include <string>  // For std::string operations
#include <vector>  // For std::vector

HistogramUtils::HistogramUtils() {}

HistogramUtils::~HistogramUtils() {
  // Clean up histograms
  for (auto &pair : light_output_spectra_) {
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
  for (auto &pair : shape_indicator_spectra_) {
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

  // Clean histogram name (remove spaces, special characters)
  std::string clean_name = source_name;
  std::replace(clean_name.begin(), clean_name.end(), ' ', '_');
  std::replace(clean_name.begin(), clean_name.end(), '-', '_');
  std::replace(clean_name.begin(), clean_name.end(), '&', '_');

  Int_t light_output_bins = config_.GetLightOutputBins();
  std::string light_hist_name = "h_light_output_" + clean_name;
  std::string light_hist_title =
      source_name + " Light Output;Light Output [keVee];Counts / " +
      std::to_string(Int_t(config_.light_output_bin_width)) + " keV";

  light_output_spectra_[source_id] = new TH1F(
      light_hist_name.c_str(), light_hist_title.c_str(), light_output_bins,
      config_.light_output_min, config_.light_output_max);

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

  shape_indicator_spectra_[source_id] =
      new TH1F(si_hist_name.c_str(), si_hist_title.c_str(), config_.si_bins,
               config_.si_min, config_.si_max);

  // Set ownership to avoid ROOT auto-deletion
  light_output_spectra_[source_id]->SetDirectory(0);
  integral_spectra_[source_id]->SetDirectory(0);
  pulse_height_spectra_[source_id]->SetDirectory(0);
  charge_comparison_spectra_[source_id]->SetDirectory(0);
  shape_indicator_spectra_[source_id]->SetDirectory(0);
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

TH1F *HistogramUtils::GetLightOutputSpectrum(Int_t source_id) const {
  auto it = light_output_spectra_.find(source_id);
  return (it != light_output_spectra_.end()) ? it->second : nullptr;
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

TH1F *HistogramUtils::GetShapeIndicatorSpectrum(Int_t source_id) const {
  auto it = shape_indicator_spectra_.find(source_id);
  return (it != shape_indicator_spectra_.end()) ? it->second : nullptr;
}

Bool_t HistogramUtils::LoadFromFile(const std::string &filename) {
  TFile *file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open file " << filename << std::endl;
    return kFALSE;
  }

  std::cout << "Loading histograms from " << filename << "..." << std::endl;

  // Load histograms from each directory using the OLD reliable approach
  for (const auto &source_pair : source_names_) {
    Int_t source_id = source_pair.first;
    std::string source_name = source_pair.second;

    // Clean name for histogram lookup
    std::string clean_name = source_name;
    std::replace(clean_name.begin(), clean_name.end(), ' ', '_');
    std::replace(clean_name.begin(), clean_name.end(), '-', '_');
    std::replace(clean_name.begin(), clean_name.end(), '&', '_');

    // Try to load each histogram type
    std::string light_hist_name =
        "light_output_spectra/h_light_output_" + clean_name;
    std::string integral_hist_name =
        "integral_spectra/h_integral_" + clean_name;
    std::string ph_hist_name =
        "pulse_height_spectra/h_pulse_height_" + clean_name;
    std::string cc_hist_name =
        "charge_comparison_spectra/h_charge_comparison_" + clean_name;
    std::string si_hist_name =
        "shape_indicator_spectra/h_shape_indicator_" + clean_name;

    TH1F *light_hist = static_cast<TH1F *>(file->Get(light_hist_name.c_str()));
    TH1F *integral_hist =
        static_cast<TH1F *>(file->Get(integral_hist_name.c_str()));
    TH1F *ph_hist = static_cast<TH1F *>(file->Get(ph_hist_name.c_str()));
    TH1F *cc_hist = static_cast<TH1F *>(file->Get(cc_hist_name.c_str()));
    TH1F *si_hist = static_cast<TH1F *>(file->Get(si_hist_name.c_str()));

    if (light_hist) {
      light_hist->SetDirectory(0);
      light_output_spectra_[source_id] = light_hist;
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
    if (si_hist) {
      si_hist->SetDirectory(0);
      shape_indicator_spectra_[source_id] = si_hist;
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
  TDirectory *light_output_dir = file->mkdir("light_output_spectra");
  TDirectory *integral_dir = file->mkdir("integral_spectra");
  TDirectory *pulse_height_dir = file->mkdir("pulse_height_spectra");
  TDirectory *cc_dir = file->mkdir("charge_comparison_spectra");
  TDirectory *si_dir = file->mkdir("shape_indicator_spectra");

  // Save light output spectra
  light_output_dir->cd();
  for (const auto &pair : light_output_spectra_) {
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

  si_dir->cd();
  for (const auto &pair : shape_indicator_spectra_) {
    pair.second->Write("", TObject::kOverwrite);
  }

  file->Close();
  delete file;

  std::cout << "Histograms saved successfully." << std::endl;
}

void HistogramUtils::LoadMeasurementStatistics(
    const std::vector<std::string> &directories,
    const std::vector<Int_t> &source_ids) {
  std::cout << "Loading measurement statistics..." << std::endl;

  for (size_t i = 0; i < directories.size() && i < source_ids.size(); ++i) {
    std::string stats_file = FindStatsFile(directories[i]);

    if (!stats_file.empty()) {
      if (ParseStatsFile(stats_file, source_ids[i])) {
        std::cout << "Loaded stats for source " << source_ids[i] << " from "
                  << stats_file << std::endl;
      }
    } else {
      std::cout << "Warning: No stats file found in " << directories[i]
                << std::endl;
    }
  }

  background_source_id_ = -1;
}

std::string HistogramUtils::FindStatsFile(const std::string &directory) {
  // Look for .txt files in the directory
  glob_t glob_result;
  std::string pattern = directory + "/*.txt";

  if (glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result) == 0) {
    if (glob_result.gl_pathc > 0) {
      std::string result = glob_result.gl_pathv[0];
      globfree(&glob_result);
      return result;
    }
    globfree(&glob_result);
  }

  return "";
}

Double_t HistogramUtils::ParseTimeString(const std::string &time_str) {
  // Parse format "0:09:53.037" to seconds
  std::istringstream iss(time_str);
  std::string token;

  // Split by ':'
  std::vector<std::string> time_parts;
  while (std::getline(iss, token, ':')) {
    time_parts.push_back(token);
  }

  if (time_parts.size() >= 3) {
    Double_t hours = std::stod(time_parts[0]);
    Double_t minutes = std::stod(time_parts[1]);
    Double_t seconds = std::stod(time_parts[2]);

    return hours * 3600 + minutes * 60 + seconds;
  }

  return 0.0;
}

Bool_t HistogramUtils::ParseStatsFile(const std::string &filepath,
                                      Int_t source_id) {
  std::ifstream file(filepath);
  if (!file.is_open()) {
    std::cout << "Error: Could not open stats file " << filepath << std::endl;
    return kFALSE;
  }

  MeasurementStats stats;
  std::string line;

  while (std::getline(file, line)) {
    if (line.find("Run ID = ") != std::string::npos) {
      stats.run_id = line.substr(line.find("= ") + 2);
    } else if (line.find("Live time = ") != std::string::npos) {
      std::string time_str = line.substr(line.find("= ") + 2);
      stats.live_time_seconds = ParseTimeString(time_str);
    } else if (line.find("Output counts = ") != std::string::npos) {
      std::string counts_line = line.substr(line.find("= ") + 2);
      std::istringstream iss(counts_line);
      iss >> stats.output_counts;

      // Extract rate
      size_t rate_pos = counts_line.find("rate (cps) = ");
      if (rate_pos != std::string::npos) {
        std::string rate_str = counts_line.substr(rate_pos + 13);
        std::istringstream rate_iss(rate_str);
        rate_iss >> stats.output_rate_cps;
      }
    }
  }

  source_stats_[source_id] = stats;
  return kTRUE;
}

void HistogramUtils::SubtractSpectrumType(Int_t source_id,
                                          const std::string &spectrum_type,
                                          Double_t source_norm_factor,
                                          Double_t bg_norm_factor) {
  TH1F *source_hist = nullptr;
  TH1F *bg_hist = nullptr;

  // Get the appropriate histograms
  if (spectrum_type == "light_output") {
    auto source_it = light_output_spectra_.find(source_id);
    auto bg_it = light_output_spectra_.find(background_source_id_);
    if (source_it != light_output_spectra_.end() &&
        bg_it != light_output_spectra_.end()) {
      source_hist = source_it->second;
      bg_hist = bg_it->second;
    }
  } else if (spectrum_type == "integral") {
    auto source_it = integral_spectra_.find(source_id);
    auto bg_it = integral_spectra_.find(background_source_id_);
    if (source_it != integral_spectra_.end() &&
        bg_it != integral_spectra_.end()) {
      source_hist = source_it->second;
      bg_hist = bg_it->second;
    }
  } else if (spectrum_type == "pulse_height") {
    auto source_it = pulse_height_spectra_.find(source_id);
    auto bg_it = pulse_height_spectra_.find(background_source_id_);
    if (source_it != pulse_height_spectra_.end() &&
        bg_it != pulse_height_spectra_.end()) {
      source_hist = source_it->second;
      bg_hist = bg_it->second;
    }
  }
  if (!source_hist || !bg_hist) {
    std::cout << "DEBUG: Missing histogram for " << spectrum_type
              << " source_id " << source_id << std::endl;
    return;
  }

  // Debug: Print scaling factors
  std::cout << "DEBUG " << spectrum_type << " source " << source_id << ":"
            << std::endl;
  std::cout << "  Source time: " << source_norm_factor << " s" << std::endl;
  std::cout << "  Background time: " << bg_norm_factor << " s" << std::endl;
  std::cout << "  Scaling factor: " << (source_norm_factor / bg_norm_factor)
            << std::endl;

  Double_t original_integral = source_hist->Integral();
  Double_t bg_integral = bg_hist->Integral();

  for (Int_t bin = 1; bin <= source_hist->GetNbinsX(); ++bin) {
    Double_t source_counts = source_hist->GetBinContent(bin);
    Double_t bg_counts = bg_hist->GetBinContent(bin);

    Double_t bg_rate = bg_counts / bg_norm_factor;

    Double_t bg_scaled_to_source = bg_rate * source_norm_factor;

    Double_t subtracted_counts = source_counts - bg_scaled_to_source;

    if (subtracted_counts < 0)
      subtracted_counts = 0;

    source_hist->SetBinContent(bin, subtracted_counts);

    // Propagate uncertainties
    Double_t source_error = source_hist->GetBinError(bin);
    Double_t bg_error = bg_hist->GetBinError(bin);

    Double_t bg_rate_error = bg_error / bg_norm_factor;

    Double_t bg_scaled_error = bg_rate_error * source_norm_factor;

    // Combined error (assuming independent uncertainties)
    Double_t combined_error = TMath::Sqrt(source_error * source_error +
                                          bg_scaled_error * bg_scaled_error);

    source_hist->SetBinError(bin, combined_error);
  }

  Double_t final_integral = source_hist->Integral();
  std::cout << "  Original integral: " << original_integral << std::endl;
  std::cout << "  Background integral: " << bg_integral << std::endl;
  std::cout << "  Final integral: " << final_integral << std::endl;
  std::cout << "  Change: " << (final_integral - original_integral)
            << std::endl;
  std::string current_title = source_hist->GetTitle();
  std::string new_title = current_title + " (Background Subtracted)";
  source_hist->SetTitle(new_title.c_str());

  source_hist->SetEntries(source_hist->Integral());
}

TH1F *HistogramUtils::GetBackgroundSubtractedSpectrum(
    Int_t source_id, const std::string &spectrum_type) const {
  // This returns the spectrum after background subtraction has been applied
  if (spectrum_type == "light_output") {
    return GetLightOutputSpectrum(source_id);
  } else if (spectrum_type == "integral") {
    return GetIntegralSpectrum(source_id);
  } else if (spectrum_type == "pulse_height") {
    return GetPulseHeightSpectrum(source_id);
  } else if (spectrum_type == "charge_comparison") {
    return GetChargeComparisonSpectrum(source_id);
  }

  return nullptr;
}

void HistogramUtils::ApplyBackgroundSubtraction() {
  if (background_source_id_ < 0) {
    std::cout << "Warning: No background source set for subtraction"
              << std::endl;
    return;
  }

  std::cout << "Applying background subtraction using source ID "
            << background_source_id_ << "..." << std::endl;

  // Get background normalization factor (live time)
  auto bg_it = source_stats_.find(background_source_id_);
  if (bg_it == source_stats_.end()) {
    std::cout << "Error: No statistics found for background source"
              << std::endl;
    return;
  }

  Double_t bg_norm_factor = bg_it->second.live_time_seconds;
  std::cout << "Background live time: " << bg_norm_factor << " seconds"
            << std::endl;

  for (const auto &source_pair : source_names_) {
    Int_t source_id = source_pair.first;

    if (source_id == background_source_id_)
      continue; // Skip background itself

    std::string source_name = source_pair.second;

    auto source_stats_it = source_stats_.find(source_id);
    if (source_stats_it == source_stats_.end()) {
      std::cout << "Warning: No statistics found for source " << source_name
                << std::endl;
      continue;
    }

    Double_t source_norm_factor = source_stats_it->second.live_time_seconds;
    std::cout << "Processing " << source_name
              << " (live time: " << source_norm_factor << " s)..." << std::endl;

    // Subtract background from each spectrum type
    SubtractSpectrumType(source_id, "light_output", source_norm_factor,
                         bg_norm_factor);
    SubtractSpectrumType(source_id, "integral", source_norm_factor,
                         bg_norm_factor);
    SubtractSpectrumType(source_id, "pulse_height", source_norm_factor,
                         bg_norm_factor);
    SubtractSpectrumType(source_id, "charge_comparison", source_norm_factor,
                         bg_norm_factor);
  }

  std::cout << "Background subtraction complete." << std::endl;
}
