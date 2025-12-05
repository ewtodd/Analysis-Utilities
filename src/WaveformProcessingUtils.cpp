#include "WaveformProcessingUtils.hpp"
#include <TMath.h>
#include <algorithm>
#include <glob.h>
#include <iostream>
#include <numeric>

WaveformProcessingUtils::WaveformProcessingUtils()
    : polarity_(1), trigger_threshold_(0.15), pre_samples_(10),
      post_samples_(100), max_events_(-1), verbose_(kFALSE),
      output_file_(nullptr), output_tree_(nullptr), store_waveforms_(kTRUE),
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

Bool_t WaveformProcessingUtils::ProcessFile(const TString filepath,
                                            const TString output_filename) {

  // Create output file and tree
  output_file_ = new TFile(output_filename, "RECREATE");
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
  output_tree_->Branch("passes_cuts", &current_features_.passes_cuts,
                       "passes_cuts/O");
  output_tree_->Branch("negative_fraction",
                       &current_features_.negative_fraction,
                       "negative_fraction/F");
  output_tree_->Branch("timestamp", &current_features_.timestamp,
                       "timestamp/l");

  if (store_waveforms_) {
    current_waveform_ = nullptr;
    output_tree_->Branch("Samples", &current_waveform_);
    std::cout << "Storing waveforms that pass cuts." << std::endl;
  }

  TFile *file = TFile::Open(filepath, "READ");
  if (!file || file->IsZombie()) {
    if (verbose_) {
      std::cout << "Error opening file: " << filepath << std::endl;
    }
  }

  TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
  if (!tree) {
    if (verbose_) {
      std::cout << "Error: TTree 'Data_R' not found in " << filepath
                << std::endl;
    }
    file->Close();
  }

  TArrayS *samples = new TArrayS();
  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("Timestamp", &current_timestamp_);

  Long64_t n_entries = tree->GetEntries();
  tree->GetEntry(0);

  for (Long64_t entry = 0; entry < n_entries; ++entry) {
    if (max_events_ > 0 && stats_.accepted >= max_events_) {
      break;
    }

    if (tree->GetEntry(entry) <= 0)
      continue;

    // Convert TArrayS to std::vector
    std::vector<Short_t> waveform_data;
    waveform_data.reserve(samples->GetSize());
    for (Int_t i = 0; i < samples->GetSize(); ++i) {
      waveform_data.push_back(samples->At(i));
    }
    stats_.total_processed++;
    if (ProcessWaveform(waveform_data)) {
    }
  }

  delete samples;
  file->Close();

  output_file_->cd();
  // Write and close
  output_tree_->Write("", TObject::kOverwrite);
  output_file_->Close();

  if (verbose_) {
    PrintAllStatistics();
  }

  return kTRUE;
}

Bool_t
WaveformProcessingUtils::ProcessWaveform(const std::vector<Short_t> &samples) {
  // Baseline subtraction
  std::vector<Float_t> processed_wf = SubtractBaseline(samples);

  // Find trigger
  Int_t trigger_pos = FindTrigger(processed_wf);
  if (trigger_pos < 0) {
    stats_.rejected_no_trigger++;
    return kFALSE;
  }

  // Check sufficient samples
  if (trigger_pos < pre_samples_ ||
      (Int_t(processed_wf.size()) - trigger_pos) <= post_samples_) {
    stats_.rejected_insufficient_samples++;
    return kFALSE;
  }

  // Crop waveform
  std::vector<Float_t> cropped_wf = CropWaveform(processed_wf, trigger_pos);

  // Extract features
  WaveformFeatures features = ExtractFeatures(cropped_wf);

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

  stats_.accepted++;
  return kTRUE;
}

std::vector<Float_t>
WaveformProcessingUtils::SubtractBaseline(const std::vector<Short_t> &samples) {
  Float_t baseline = 0;
  Int_t baseline_samples = TMath::Min(10, Int_t(samples.size()));

  for (Int_t i = 0; i < baseline_samples; ++i) {
    baseline += samples[i];
  }
  baseline /= baseline_samples;

  std::vector<Float_t> processed;
  processed.reserve(samples.size());

  for (size_t i = 0; i < samples.size(); ++i) {
    if (polarity_ == -1) {
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

WaveformFeatures WaveformProcessingUtils::ExtractFeatures(
    const std::vector<Float_t> &cropped_wf) {
  WaveformFeatures features;
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
  features.timestamp = current_timestamp_;

  features.passes_cuts = kTRUE; // Will be updated in ApplyQualityCuts
  features.negative_fraction =
      Float_t(negative_samples) / Float_t(long_end - pre_samples_);

  return features;
}

Bool_t
WaveformProcessingUtils::ApplyQualityCuts(const WaveformFeatures &features) {

  if (features.long_integral <= 0) {
    stats_.rejected_negative_integral++;
    return kFALSE;
  }

  if (features.negative_fraction > 0.4) {
    stats_.rejected_baseline++;
    return kFALSE;
  }

  if (features.long_integral <= 0) {
    stats_.rejected_negative_integral++;
    return kFALSE;
  }
  return kTRUE;
}

void WaveformProcessingUtils::PrintAllStatistics() const {
  std::cout << "Waveform processing statistics..." << std::endl;
  std::cout << "Total processed: " << stats_.total_processed << std::endl;
  std::cout << std::endl;
  std::cout << "Accepted: " << stats_.accepted << std::endl;
  std::cout << std::endl;
  std::cout << "Rejected breakdown:" << std::endl;
  std::cout << "No trigger: " << stats_.rejected_no_trigger << std::endl;
  std::cout << "Insufficient samples: " << stats_.rejected_insufficient_samples
            << std::endl;
  std::cout << "Negative integral: " << stats_.rejected_negative_integral
            << std::endl;
  std::cout << "Bad baseline: " << stats_.rejected_baseline << std::endl;
  std::cout << std::endl;

  if (stats_.total_processed > 0) {
    std::cout << "  Acceptance rate: "
              << 100 * Float_t(stats_.accepted) /
                     Float_t(stats_.total_processed)
              << "%" << std::endl;
  }
  std::cout << std::endl;
}
