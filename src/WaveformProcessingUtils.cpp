#include "WaveformProcessingUtils.hpp"

WaveformProcessingUtils::WaveformProcessingUtils()
    : polarity_(1), trigger_threshold_(0.15), num_samples_baseline_(10),
      pre_samples_(10), post_samples_(100), max_events_(-1), verbose_(kFALSE),
      sample_waveforms_to_save_(0), sample_waveforms_saved_(0),
      output_file_(nullptr), output_tree_(nullptr), store_waveforms_(kTRUE),
      current_waveform_(nullptr) {}

WaveformProcessingUtils::WaveformProcessingUtils(
    const FileProcessingConfig &config)
    : polarity_(config.polarity), trigger_threshold_(config.trigger_threshold),
      num_samples_baseline_(config.num_samples_baseline),
      pre_samples_(config.pre_samples), post_samples_(config.post_samples),
      pre_gate_(config.pre_gate), short_gate_(config.short_gate),
      long_gate_(config.long_gate), max_events_(config.max_events),
      verbose_(config.verbose),
      sample_waveforms_to_save_(config.sample_waveforms_to_save),
      sample_waveforms_saved_(0), output_file_(nullptr), output_tree_(nullptr),
      store_waveforms_(config.store_waveforms), current_waveform_(nullptr) {}

WaveformProcessingUtils::~WaveformProcessingUtils() {
  if (current_waveform_) {
    delete current_waveform_;
    current_waveform_ = nullptr;
  }
  if (output_file_) {
    if (output_file_->IsOpen()) {
      output_file_->Close();
    }
    delete output_file_;
    output_file_ = nullptr;
  }
}

Bool_t
WaveformProcessingUtils::ProcessWaveform(const std::vector<Short_t> &samples) {
  Short_t raw_max = (polarity_ == 1)
                        ? *std::max_element(samples.begin(), samples.end())
                        : *std::min_element(samples.begin(), samples.end());

  std::vector<Float_t> processed_wf = SubtractBaseline(samples);

  Float_t trigger_pos = FindTrigger(processed_wf);
  if (trigger_pos < 0) {
    stats_.rejected_no_trigger++;
    return kFALSE;
  }

  Int_t trigger_pos_int = Int_t(trigger_pos);
  if (trigger_pos_int < pre_samples_ ||
      (Int_t(processed_wf.size()) - trigger_pos_int) <= post_samples_) {
    stats_.rejected_insufficient_samples++;
    return kFALSE;
  }

  std::vector<Float_t> cropped_wf = CropWaveform(processed_wf, trigger_pos);

  WaveformFeatures features = ExtractFeatures(cropped_wf);
  features.raw_pulse_height = std::abs(raw_max);
  features.trigger_position = trigger_pos;
  Bool_t passes_cuts = ApplyQualityCuts(features);
  features.passes_cuts = passes_cuts;

  if (!passes_cuts) {
    return kFALSE;
  }

  if (store_waveforms_) {
    if (current_waveform_)
      delete current_waveform_;
    current_waveform_ = new TArrayS(cropped_wf.size());
    for (size_t i = 0; i < cropped_wf.size(); ++i) {
      current_waveform_->SetAt(Short_t(cropped_wf[i]), i);
    }
  }

  if (sample_waveforms_saved_ < sample_waveforms_to_save_) {
    SaveSampleWaveform(cropped_wf);
  }

  current_features_ = features;
  output_tree_->Fill();

  stats_.accepted++;
  return kTRUE;
}

std::mutex WaveformProcessingUtils::canvas_mutex_;

void WaveformProcessingUtils::SaveSampleWaveform(
    const std::vector<Float_t> &waveform) {

  std::lock_guard<std::mutex> lock(canvas_mutex_);

  if (gSystem->AccessPathName("plots/samplewaveforms")) {
    gSystem->mkdir("plots/samplewaveforms", kTRUE);
  }

  Int_t n = waveform.size();
  std::vector<Double_t> x(n), y(n);
  for (Int_t i = 0; i < n; ++i) {
    x[i] = i;
    y[i] = waveform[i];
  }

  TGraph *graph = new TGraph(n, x.data(), y.data());
  TCanvas *canvas =
      new TCanvas(PlottingUtils::GetRandomName(), "Waveform", 1200, 800);
  PlottingUtils::SetStylePreferences();
  PlottingUtils::ConfigureCanvas(canvas);
  PlottingUtils::ConfigureGraph(graph, kBlue + 1, ";Sample;Amplitude [ADC]");
  graph->Draw("AL");

  TString filename = Form("plots/samplewaveforms/%s_waveform_%04d.png",
                          current_output_name_.Data(), sample_waveforms_saved_);
  canvas->SaveAs(filename);

  delete graph;
  delete canvas;

  sample_waveforms_saved_++;
}

std::vector<Float_t>
WaveformProcessingUtils::SubtractBaseline(const std::vector<Short_t> &samples) {
  Float_t baseline = 0;
  Int_t baseline_samples =
      TMath::Min(num_samples_baseline_, Int_t(samples.size()));

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

Float_t
WaveformProcessingUtils::FindTrigger(const std::vector<Float_t> &waveform) {

  Float_t peak_value = *std::max_element(waveform.begin(), waveform.end());
  Float_t trigger_level = peak_value * trigger_threshold_;

  for (size_t i = 0; i < waveform.size(); ++i) {
    if (waveform[i] >= trigger_level) {
      return Float_t(i);
    }
  }

  return -1.0;
}

std::vector<Float_t>
WaveformProcessingUtils::CropWaveform(const std::vector<Float_t> &waveform,
                                      Float_t trigger_pos) {
  Int_t trigger_pos_int = Int_t(trigger_pos);
  Int_t start = trigger_pos_int - pre_samples_;
  Int_t end = trigger_pos_int + post_samples_;

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
  Int_t integration_start = pre_samples_ - pre_gate_;

  auto max_it = std::max_element(cropped_wf.begin(), cropped_wf.end());
  features.pulse_height = *max_it;
  features.peak_position = std::distance(cropped_wf.begin(), max_it);

  features.short_integral = 0;
  features.long_integral = 0;

  Int_t negative_samples = 0;
  Int_t short_end =
      TMath::Min(integration_start + short_gate_, Int_t(cropped_wf.size()));
  Int_t long_end =
      TMath::Min(integration_start + long_gate_, Int_t(cropped_wf.size()));

  for (Int_t i = integration_start; i < long_end; ++i) {
    Float_t sample_value = cropped_wf[i];
    features.long_integral += sample_value;
    if (i < short_end) {
      features.short_integral += sample_value;
    }
    if (sample_value < 0)
      negative_samples++;
  }
  features.timestamp = current_timestamp_;

  features.passes_cuts = kTRUE;
  features.negative_fraction =
      Float_t(negative_samples) / Float_t(long_end - integration_start);

  return features;
}

Bool_t
WaveformProcessingUtils::ApplyQualityCuts(const WaveformFeatures &features) {

  if (((features.raw_pulse_height == 16384) && (polarity_ == 1)) ||
      ((features.raw_pulse_height == 0) && (polarity_ == -1))) {
    stats_.rejected_clipped++;
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
  std::cout << "Rejected no trigger: " << stats_.rejected_no_trigger
            << std::endl;
  std::cout << "Rejected clipped ADC: " << stats_.rejected_clipped << std::endl;
  std::cout << "Rejected insufficient samples: "
            << stats_.rejected_insufficient_samples << std::endl;
  std::cout << "Rejected negative integral: "
            << stats_.rejected_negative_integral << std::endl;
  std::cout << "Rejected bad baseline: " << stats_.rejected_baseline
            << std::endl;
  std::cout << std::endl;

  if (stats_.total_processed > 0) {
    std::cout << "Acceptance rate: "
              << 100 * Float_t(stats_.accepted) /
                     Float_t(stats_.total_processed)
              << "%" << std::endl;
  }
  std::cout << std::endl;

  std::ofstream stats_file(
      "root_files/" + std::string(current_output_name_.Data()) + ".stats",
      std::ios::app);
  if (stats_file.is_open()) {
    stats_file << "Waveform processing statistics..." << std::endl;
    stats_file << "Total processed: " << stats_.total_processed << std::endl;
    stats_file << std::endl;
    stats_file << "Accepted: " << stats_.accepted << std::endl;
    stats_file << std::endl;
    stats_file << "Rejected no trigger: " << stats_.rejected_no_trigger
               << std::endl;
    stats_file << "Rejected clipped ADC: " << stats_.rejected_clipped
               << std::endl;
    stats_file << "Rejected insufficient samples: "
               << stats_.rejected_insufficient_samples << std::endl;
    stats_file << "Rejected negative integral: "
               << stats_.rejected_negative_integral << std::endl;
    stats_file << "Rejected bad baseline: " << stats_.rejected_baseline
               << std::endl;
    stats_file << std::endl;

    if (stats_.total_processed > 0) {
      stats_file << "Acceptance rate: "
                 << 100 * Float_t(stats_.accepted) /
                        Float_t(stats_.total_processed)
                 << "%" << std::endl;
    }
    stats_file << std::endl;
  }
}

Bool_t WaveformProcessingUtils::ProcessFile(const TString filepath,
                                            const TString output_name) {
  current_output_name_ = output_name;
  sample_waveforms_saved_ = 0;

  if (gSystem->AccessPathName("root_files")) {
    gSystem->mkdir("root_files", kTRUE);
  }

  // clear file
  std::ofstream("root_files/" + std::string(output_name.Data()) + ".stats",
                std::ios::trunc);

  TString output_filename = "root_files/" + output_name + ".root";
  output_file_ = new TFile(output_filename, "RECREATE");
  if (!output_file_ || output_file_->IsZombie()) {
    std::cout << "Error: Could not create output file " << output_filename
              << std::endl;
    return kFALSE;
  }

  output_tree_ = new TTree("features", "Waveform Features");

  output_tree_->Branch("pulse_height", &current_features_.pulse_height,
                       "pulse_height/F");
  output_tree_->Branch("trigger_position", &current_features_.trigger_position,
                       "trigger_position/F");
  output_tree_->Branch("short_integral", &current_features_.short_integral,
                       "short_integral/F");
  output_tree_->Branch("long_integral", &current_features_.long_integral,
                       "long_integral/F");
  output_tree_->Branch("passes_cuts", &current_features_.passes_cuts,
                       "passes_cuts/O");
  output_tree_->Branch("timestamp", &current_features_.timestamp,
                       "timestamp/l");

  if (store_waveforms_) {
    current_waveform_ = nullptr;
    output_tree_->Branch("Samples", &current_waveform_);
    std::cout << "Storing waveforms that pass cuts." << std::endl;
  }

  TFile *file = TFile::Open(filepath, "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error opening file: " << filepath << std::endl;
    return kFALSE;
  }

  TTree *tree = static_cast<TTree *>(file->Get("Data_R"));
  if (!tree) {
    std::cout << "Error: TTree 'Data_R' not found in " << filepath << std::endl;
    file->Close();
    return kFALSE;
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

    std::vector<Short_t> waveform_data;
    waveform_data.reserve(samples->GetSize());
    for (Int_t i = 0; i < samples->GetSize(); ++i) {
      waveform_data.push_back(samples->At(i));
    }
    stats_.total_processed++;
    ProcessWaveform(waveform_data);
  }

  delete samples;
  file->Close();

  output_file_->cd();
  output_tree_->Write("", TObject::kOverwrite);
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;
  output_tree_ = nullptr;
  current_waveform_ = nullptr;

  if (verbose_) {
    PrintAllStatistics();
  }

  return kTRUE;
}

void WaveformProcessingUtils::ProcessFilesParallel(
    const std::vector<TString> &filepaths,
    const std::vector<TString> &output_names,
    const FileProcessingConfig &config, Int_t max_workers) {

  ROOT::EnableThreadSafety();

  Int_t n_files = Int_t(filepaths.size());
  Int_t n_workers = max_workers > 0
                        ? max_workers
                        : Int_t(std::thread::hardware_concurrency());
  n_workers = TMath::Min(n_workers, n_files);

  std::cout << "Processing " << n_files << " files with " << n_workers
            << " workers." << std::endl;

  std::function<Bool_t(const TString &, const TString &)> process_one =
      [&config](const TString &filepath, const TString &output_name) -> Bool_t {
    WaveformProcessingUtils *processor = new WaveformProcessingUtils(config);
    Bool_t result = processor->ProcessFile(filepath, output_name);
    delete processor;
    return result;
  };

  for (Int_t i = 0; i < n_files; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_files);

    for (Int_t j = i; j < batch_end; ++j) {
      futures.push_back(std::async(std::launch::async, process_one,
                                   std::cref(filepaths[j]),
                                   std::cref(output_names[j])));
    }

    for (size_t j = 0; j < futures.size(); ++j) {
      Bool_t result = futures[j].get();
      std::cout << "Finished: " << output_names[i + j]
                << (result ? " [OK]" : " [FAILED]") << std::endl;
    }
  }
}
