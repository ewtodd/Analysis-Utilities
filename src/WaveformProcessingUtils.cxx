#include "WaveformProcessingUtils.h"
#include <TMath.h>
#include <algorithm>
#include <glob.h>
#include <iostream>
#include <numeric>

WaveformProcessingUtils::WaveformProcessingUtils()
    : polarity_("positive"), trigger_threshold_(0.1), pre_samples_(20),
      post_samples_(200), max_events_(-1), verbose_(kFALSE),
      output_file_(nullptr), output_tree_(nullptr), store_waveforms_(kFALSE),
      current_waveform_(nullptr) {}

WaveformProcessingUtils::~WaveformProcessingUtils() {
  if (current_waveform_) {
    delete current_waveform_;
  }
  if (output_file_) {
    output_file_->Close();
    delete output_file_;
  }
}
void WaveformProcessingUtils::AddSource(const std::string &name, Int_t id) {
  source_map_[name] = id;
  stats_[id] = ProcessingStats();
}

Bool_t WaveformProcessingUtils::ProcessDataSets(
    const std::vector<std::string> &filepaths,
    const std::vector<std::string> &source_names,
    const std::string &output_filename) {
  if (filepaths.size() != source_names.size()) {
    std::cout << "Error: Number of filepaths and source names must match!"
              << std::endl;
    return kFALSE;
  }

  // Create output file and tree
  output_file_ = new TFile(output_filename.c_str(), "RECREATE");
  if (!output_file_ || output_file_->IsZombie()) {
    std::cout << "Error: Could not create output file " << output_filename
              << std::endl;
    return kFALSE;
  }

  output_tree_ = new TTree("features", "Waveform Features");

  // Set up feature branches
  output_tree_->Branch("pulse_height", &current_features_.pulse_height,
                       "pulse_height/F");
  output_tree_->Branch("peak_position", &current_features_.peak_position,
                       "peak_position/I");
  output_tree_->Branch("trigger_position", &current_features_.trigger_position,
                       "trigger_position/I");
  output_tree_->Branch("long_integral", &current_features_.long_integral,
                       "long_integral/F");
  output_tree_->Branch("source_id", &current_features_.source_id,
                       "source_id/I");
  output_tree_->Branch("passes_cuts", &current_features_.passes_cuts,
                       "passes_cuts/O");
  output_tree_->Branch("negative_fraction",
                       &current_features_.negative_fraction,
                       "negative_fraction/F");

  // Add waveform branch if requested
  if (store_waveforms_) {
    current_waveform_ = nullptr;
    output_tree_->Branch("Samples", &current_waveform_);
    std::cout << "Waveform storage enabled - storing TArrayS data for events "
                 "that pass cuts"
              << std::endl;
  }

  // Process each dataset
  for (size_t i = 0; i < filepaths.size(); ++i) {
    if (verbose_) {
      std::cout << "Processing: " << filepaths[i] << " (" << source_names[i]
                << ")" << std::endl;
    }

    Int_t source_id = source_map_[source_names[i]];
    if (!ProcessSingleFile(filepaths[i], source_id)) {
      std::cout << "Warning: Failed to process " << filepaths[i] << std::endl;
    }
  }

  output_file_->cd();
  // Write and close
  output_tree_->Write("", TObject::kOverwrite);
  output_file_->Close();

  if (verbose_) {
    PrintAllStatistics();
  }

  return kTRUE;
}

Bool_t WaveformProcessingUtils::ProcessSingleFile(const std::string &filepath,
                                                  Int_t source_id) {
  std::vector<std::string> root_files = FindROOTFiles(filepath);

  if (root_files.empty()) {
    if (verbose_) {
      std::cout << "No ROOT files found in " << filepath << std::endl;
    }
    return kFALSE;
  }

  if (verbose_) {
    std::cout << "  Found " << root_files.size() << " ROOT files" << std::endl;
  }

  Int_t total_processed_this_source = 0;

  for (const std::string &root_file : root_files) {
    TFile *file = TFile::Open(root_file.c_str(), "READ");
    if (!file || file->IsZombie()) {
      if (verbose_) {
        std::cout << "  Error opening file: " << root_file << std::endl;
      }
      continue;
    }

    TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
    if (!tree) {
      if (verbose_) {
        std::cout << "  Error: TTree 'Data_R' not found in " << root_file
                  << std::endl;
      }
      file->Close();
      continue;
    }

    TArrayS *samples = new TArrayS();
    tree->SetBranchAddress("Samples", &samples);

    Long64_t n_entries = tree->GetEntries();
    Int_t processed_this_file = 0;
    tree->GetEntry(0);
    for (Long64_t entry = 0; entry < n_entries; ++entry) {
      if (max_events_ > 0 && stats_[source_id].accepted >= max_events_) {
        break;
      }

      if (tree->GetEntry(entry) <= 0)
        continue;

      stats_[source_id].total_processed++;

      // Convert TArrayS to std::vector
      std::vector<Short_t> waveform_data;
      waveform_data.reserve(samples->GetSize());
      for (Int_t i = 0; i < samples->GetSize(); ++i) {
        waveform_data.push_back(samples->At(i));
      }
      if (ProcessWaveform(waveform_data, source_id)) {
        processed_this_file++;
      }
    }

    total_processed_this_source += processed_this_file;

    if (verbose_) {
      std::cout << "    " << root_file << ": " << processed_this_file
                << " events" << std::endl;
    }

    delete samples;
    file->Close();

    if (max_events_ > 0 && stats_[source_id].accepted >= max_events_) {
      break;
    }
  }

  if (verbose_) {
    std::cout << "  Total processed for this source: "
              << total_processed_this_source << std::endl;
  }

  return kTRUE;
}

Bool_t
WaveformProcessingUtils::ProcessWaveform(const std::vector<Short_t> &samples,
                                         Int_t source_id) {
  // Baseline subtraction
  std::vector<Float_t> processed_wf = SubtractBaseline(samples);

  // Find trigger
  Int_t trigger_pos = FindTrigger(processed_wf);
  if (trigger_pos < 0) {
    stats_[source_id].rejected_no_trigger++;
    return kFALSE;
  }

  // Check sufficient samples
  if (trigger_pos < pre_samples_ ||
      (Int_t(processed_wf.size()) - trigger_pos) <= post_samples_) {
    stats_[source_id].rejected_insufficient_samples++;
    return kFALSE;
  }

  // Crop waveform
  std::vector<Float_t> cropped_wf = CropWaveform(processed_wf, trigger_pos);

  // Extract features
  WaveformFeatures features = ExtractFeatures(cropped_wf, source_id);

  // Apply quality cuts
  Bool_t passes_cuts = ApplyQualityCuts(features);
  features.passes_cuts = passes_cuts;

  // ONLY process further if cuts are passed
  if (!passes_cuts) {
    return kFALSE;
  }

  // Store the CROPPED waveform if requested (ONLY for events that pass cuts)
  if (store_waveforms_) {
    // Create TArrayS from cropped, baseline-subtracted waveform
    if (current_waveform_)
      delete current_waveform_;
    current_waveform_ = new TArrayS(cropped_wf.size());
    for (size_t i = 0; i < cropped_wf.size(); ++i) {
      current_waveform_->SetAt(Short_t(cropped_wf[i]), i);
    }
  }

  current_features_ = features;
  output_tree_->Fill();

  stats_[source_id].accepted++;
  return kTRUE;
}

std::vector<Float_t>
WaveformProcessingUtils::SubtractBaseline(const std::vector<Short_t> &samples) {
  // Calculate baseline from first 10 samples
  Float_t baseline = 0;
  Int_t baseline_samples = TMath::Min(10, Int_t(samples.size()));

  for (Int_t i = 0; i < baseline_samples; ++i) {
    baseline += samples[i];
  }
  baseline /= baseline_samples;

  std::vector<Float_t> processed;
  processed.reserve(samples.size());

  for (size_t i = 0; i < samples.size(); ++i) {
    if (polarity_ == "negative") {
      processed.push_back(baseline - samples[i]);
    } else {
      processed.push_back(samples[i] - baseline);
    }
  }

  return processed;
}

Int_t WaveformProcessingUtils::FindTrigger(
    const std::vector<Float_t> &waveform) {
  // Find peak value
  Float_t peak_value = *std::max_element(waveform.begin(), waveform.end());
  Float_t trigger_level = peak_value * trigger_threshold_;

  // Find first point above trigger level
  for (size_t i = 0; i < waveform.size(); ++i) {
    if (waveform[i] >= trigger_level) {
      return Int_t(i);
    }
  }

  return -1; // No trigger found
}

std::vector<Float_t>
WaveformProcessingUtils::CropWaveform(const std::vector<Float_t> &waveform,
                                      Int_t trigger_pos) {
  Int_t start = trigger_pos - pre_samples_;
  Int_t end = trigger_pos + post_samples_;

  std::vector<Float_t> cropped;
  cropped.reserve(pre_samples_ + post_samples_);

  for (Int_t i = start; i < end && i < Int_t(waveform.size()); ++i) {
    cropped.push_back(waveform[i]);
  }

  return cropped;
}

WaveformFeatures
WaveformProcessingUtils::ExtractFeatures(const std::vector<Float_t> &cropped_wf,
                                         Int_t source_id) {
  WaveformFeatures features;
  features.source_id = source_id;
  features.trigger_position = pre_samples_;

  // Find peak in cropped waveform
  auto max_it = std::max_element(cropped_wf.begin(), cropped_wf.end());
  features.pulse_height = *max_it;
  features.peak_position = std::distance(cropped_wf.begin(), max_it);

  features.long_integral = 0;

  Int_t negative_samples = 0;
  Int_t long_end = Int_t(cropped_wf.size());

  for (Int_t i = pre_samples_; i < long_end; ++i) {
    Float_t sample_value = cropped_wf[i];
    features.long_integral += sample_value;
    if (sample_value < 0)
      negative_samples++;
  }

  features.passes_cuts = kTRUE; // Will be updated in ApplyQualityCuts
  features.negative_fraction =
      Float_t(negative_samples) / Float_t(long_end - pre_samples_);

  return features;
}

Bool_t
WaveformProcessingUtils::ApplyQualityCuts(const WaveformFeatures &features) {
  Int_t source_id = features.source_id;

  // 1. Check for negative raw long integral (before clipping)
  if (features.long_integral <= 0) {
    stats_[source_id].rejected_negative_integral++;
    return kFALSE;
  }

  if (features.negative_fraction > 0.4) {
    stats_[source_id].rejected_baseline++;
    return kFALSE;
  }

  if (features.long_integral <= 0) {
    stats_[source_id].rejected_negative_integral++;
    return kFALSE;
  }
  return kTRUE;
}

std::vector<std::string>
WaveformProcessingUtils::FindROOTFiles(const std::string &basepath,
                                       const std::string &pattern) {
  std::vector<std::string> files;
  std::string full_pattern = basepath + "/RAW/" + pattern;

  glob_t glob_result;
  int glob_status = glob(full_pattern.c_str(), GLOB_TILDE, NULL, &glob_result);

  if (glob_status == 0) {
    for (unsigned int i = 0; i < glob_result.gl_pathc; ++i) {
      files.push_back(std::string(glob_result.gl_pathv[i]));
    }
  }

  globfree(&glob_result);

  if (files.empty() && verbose_) {
    std::cout << "  No files found matching pattern: " << full_pattern
              << std::endl;
  }

  return files;
}

void WaveformProcessingUtils::PrintAllStatistics() const {
  std::cout << "Waveform processing statistics..." << std::endl;

  for (const auto &pair : stats_) {
    Int_t source_id = pair.first;
    const ProcessingStats &stats = pair.second;

    // Find source name
    std::string source_name = "Unknown";
    for (const auto &source_pair : source_map_) {
      if (source_pair.second == source_id) {
        source_name = source_pair.first;
        break;
      }
    }

    stats.Print(source_name);
  }
}

ProcessingStats WaveformProcessingUtils::GetStats(Int_t source_id) const {
  auto it = stats_.find(source_id);
  if (it != stats_.end()) {
    return it->second;
  }
  return ProcessingStats();
}

// ProcessingStats methods
void ProcessingStats::Print(const std::string &source_name) const {
  std::cout << "Source: " << source_name << std::endl;
  std::cout << "  Total processed: " << total_processed << std::endl;
  std::cout << "  Accepted: " << accepted << std::endl;
  std::cout << "  Rejected breakdown:" << std::endl;
  std::cout << "    No trigger: " << rejected_no_trigger << std::endl;
  std::cout << "    Insufficient samples: " << rejected_insufficient_samples
            << std::endl;
  std::cout << "    Negative integral: " << rejected_negative_integral
            << std::endl;
  std::cout << "    Bad baseline: " << rejected_baseline << std::endl;

  if (total_processed > 0) {
    Double_t acceptance_rate = GetAcceptanceRate();
    std::cout << "  Acceptance rate: " << acceptance_rate << "%" << std::endl;
  }
  std::cout << std::endl;
}

Double_t ProcessingStats::GetAcceptanceRate() const {
  if (total_processed > 0) {
    return (100.0 * accepted) / total_processed;
  }
  return 0.0;
}
