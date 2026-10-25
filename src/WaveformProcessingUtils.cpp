#include "WaveformProcessingUtils.hpp"
#include "IOUtils.hpp"
#include <TArrayS.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TLine.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <future>
#include <thread>

// ------------------------------------------------------------ ROOT sources

RootTreeWaveformSource::RootTreeWaveformSource()
    : file_(nullptr), tree_(nullptr), entry_(0) {}

RootTreeWaveformSource::~RootTreeWaveformSource() {
  RootTreeWaveformSource::Close();
}

Bool_t RootTreeWaveformSource::Open(const TString &path) {
  Close();
  file_ = TFile::Open(path, "READ");
  if (!file_ || file_->IsZombie()) {
    std::cout << "ERROR opening file: " << path << std::endl;
    delete file_;
    file_ = nullptr;
    return kFALSE;
  }
  tree_ = static_cast<TTree *>(file_->Get("Data_R"));
  if (!tree_) {
    std::cout << "ERROR: TTree 'Data_R' not found in " << path << std::endl;
    Close();
    return kFALSE;
  }
  if (!BindBranches()) {
    std::cout << "ERROR: 'Data_R' in " << path
              << " has no waveform branch of this format" << std::endl;
    Close();
    return kFALSE;
  }
  entry_ = 0;
  return kTRUE;
}

void RootTreeWaveformSource::Close() {
  if (file_) {
    file_->Close();
    delete file_;
    file_ = nullptr;
  }
  tree_ = nullptr;
  entry_ = 0;
}

Double_t RootTreeWaveformSource::LeafValue(const char *name,
                                           Double_t fallback) const {
  TLeaf *leaf = tree_->GetLeaf(name);
  return leaf ? leaf->GetValue() : fallback;
}

CoMPASSRootSource::CoMPASSRootSource() : samples_(new TArrayS()) {}

CoMPASSRootSource::~CoMPASSRootSource() {
  Close();
  delete samples_;
}

Bool_t CoMPASSRootSource::BindBranches() {
  if (!tree_->GetBranch("Samples"))
    return kFALSE;
  tree_->SetBranchAddress("Samples", &samples_);
  return kTRUE;
}

Bool_t CoMPASSRootSource::Next(WaveformRecord &record) {
  if (!tree_)
    return kFALSE;
  const Long64_t n_entries = tree_->GetEntries();
  while (entry_ < n_entries) {
    tree_->GetEntry(entry_++);
    Int_t n = samples_->GetSize();
    if (n == 0)
      continue;
    record.samples.Set(n);
    for (Int_t i = 0; i < n; ++i)
      record.samples.SetAt(samples_->At(i), i);
    // A 64-bit picosecond timestamp does not survive a trip through double.
    TLeaf *timestamp = tree_->GetLeaf("Timestamp");
    record.timestamp =
        timestamp ? *static_cast<ULong64_t *>(timestamp->GetValuePointer()) : 0;
    record.board = Int_t(LeafValue("Board", -1));
    record.channel = Int_t(LeafValue("Channel", -1));
    record.energy = LeafValue("Energy", -1);
    record.flags = UInt_t(LeafValue("Flags", 0));
    return kTRUE;
  }
  return kFALSE;
}

WaveDumpRootSource::WaveDumpRootSource() : samples_(new TArrayS()) {}

WaveDumpRootSource::~WaveDumpRootSource() {
  Close();
  delete samples_;
}

Bool_t WaveDumpRootSource::BindBranches() {
  if (!tree_->GetBranch("Samples"))
    return kFALSE;
  tree_->SetBranchAddress("Samples", &samples_);
  return kTRUE;
}

Bool_t WaveDumpRootSource::Next(WaveformRecord &record) {
  if (!tree_)
    return kFALSE;
  const Long64_t n_entries = tree_->GetEntries();
  while (entry_ < n_entries) {
    tree_->GetEntry(entry_++);
    Int_t n = samples_->GetSize();
    if (n == 0)
      continue;
    record.samples.Set(n);
    for (Int_t i = 0; i < n; ++i)
      record.samples.SetAt(samples_->At(i), i);
    record.timestamp = ULong64_t(LeafValue("TriggerTimeTag", 0));
    record.board = -1;
    record.channel = Int_t(LeafValue("Channel", -1));
    record.energy = -1;
    record.flags = 0;
    return kTRUE;
  }
  return kFALSE;
}

SolarisRootSource::SolarisRootSource() : trace_(new TArrayI()) {}

SolarisRootSource::~SolarisRootSource() {
  Close();
  delete trace_;
}

Bool_t SolarisRootSource::BindBranches() {
  if (!tree_->GetBranch("Trace0"))
    return kFALSE;
  tree_->SetBranchAddress("Trace0", &trace_);
  return kTRUE;
}

Bool_t SolarisRootSource::Next(WaveformRecord &record) {
  if (!tree_)
    return kFALSE;
  const Long64_t n_entries = tree_->GetEntries();
  while (entry_ < n_entries) {
    tree_->GetEntry(entry_++);
    if (trace_->GetSize() == 0)
      continue;
    record.samples = *trace_;
    TLeaf *timestamp = tree_->GetLeaf("Timestamp");
    record.timestamp =
        timestamp ? *static_cast<ULong64_t *>(timestamp->GetValuePointer()) : 0;
    record.board = -1;
    record.channel = Int_t(LeafValue("Channel", -1));
    record.energy = LeafValue("Energy", -1);
    record.flags = (UInt_t(LeafValue("FlagsHigh", 0)) << 16) |
                   UInt_t(LeafValue("FlagsLow", 0));
    return kTRUE;
  }
  return kFALSE;
}

// ---------------------------------------------------------- binary sources

Bool_t CoMPASSBinarySource::Open(const TString &path) {
  return reader_.Open(path.Data());
}

Bool_t CoMPASSBinarySource::Next(WaveformRecord &record) {
  while (reader_.ReadEvent()) {
    const CoMPASSData &event = reader_.GetCurrentEvent();
    Int_t n = event.samples.GetSize();
    if (!event.hasWaveform() || n == 0)
      continue;
    record.samples.Set(n);
    for (Int_t i = 0; i < n; ++i)
      record.samples.SetAt(event.samples.At(i), i);
    record.timestamp = event.timestamp;
    record.board = event.board;
    record.channel = event.channel;
    if (event.hasEnergyCh())
      record.energy = event.energy_ch;
    else if (event.hasEnergyCal())
      record.energy = event.energy_cal;
    else
      record.energy = -1;
    record.flags = event.flags;
    return kTRUE;
  }
  return kFALSE;
}

void CoMPASSBinarySource::Close() { reader_.Close(); }

Bool_t WaveDump742BinarySource::Open(const TString &path) {
  return reader_.Open(path.Data());
}

Bool_t WaveDump742BinarySource::Next(WaveformRecord &record) {
  while (reader_.ReadEvent()) {
    const WaveDump742Data &event = reader_.GetCurrentEvent();
    Int_t n = event.samples.GetSize();
    if (n == 0)
      continue;
    record.samples.Set(n);
    for (Int_t i = 0; i < n; ++i)
      record.samples.SetAt(event.samples.At(i), i);
    record.timestamp = event.group_trigger_time_tag;
    record.board = Int_t(event.board_id);
    record.channel = Int_t(event.channel);
    record.energy = -1;
    record.flags = 0;
    return kTRUE;
  }
  return kFALSE;
}

void WaveDump742BinarySource::Close() { reader_.Close(); }

Bool_t SolarisBinarySource::Open(const TString &path) {
  return reader_.Open(path.Data());
}

Bool_t SolarisBinarySource::Next(WaveformRecord &record) {
  while (reader_.ReadEvent()) {
    const SOLData &event = reader_.GetCurrentEvent();
    const Int_t *trace = event.getAnalog0();
    Int_t n = Int_t(event.getSamples());
    if (!event.hasTraces() || !trace || n == 0)
      continue;
    record.samples.Set(n, trace);
    record.timestamp = event.timestamp;
    record.board = -1;
    record.channel = event.channel;
    record.energy = event.energy;
    record.flags = (UInt_t(event.flags_high) << 16) | UInt_t(event.flags_low);
    return kTRUE;
  }
  return kFALSE;
}

void SolarisBinarySource::Close() { reader_.Close(); }

WaveformSource *CreateWaveformSource(InputFormat format) {
  switch (format) {
  case InputFormat::kCOMPASS:
    return new CoMPASSRootSource();
  case InputFormat::kWAVEDUMP:
    return new WaveDumpRootSource();
  case InputFormat::kSOLARIS:
    return new SolarisRootSource();
  case InputFormat::kCOMPASS_BINARY:
    return new CoMPASSBinarySource();
  case InputFormat::kWAVEDUMP742_BINARY:
    return new WaveDump742BinarySource();
  case InputFormat::kSOLARIS_BINARY:
    return new SolarisBinarySource();
  }
  return nullptr;
}

// ------------------------------------------------------------- base pipeline

std::mutex WaveformProcessor::canvas_mutex_;

WaveformProcessor::WaveformProcessor(const FileProcessingConfig &config)
    : config_(config), sample_waveforms_saved_(0), output_file_(nullptr),
      output_tree_(nullptr), baseline_subtracted_(new TArrayF()),
      cropped_(new TArrayF()), current_baseline_rms_(0),
      current_baseline_rms_valid_(kFALSE), out_timestamp_(0), out_channel_(-1),
      out_board_(-1), out_digitizer_energy_(-1), out_flags_(0),
      out_trigger_position_(-1), out_pulse_height_(0), out_baseline_rms_(-1) {}

WaveformProcessor::~WaveformProcessor() {
  if (output_file_) {
    if (output_file_->IsOpen())
      output_file_->Close();
    delete output_file_;
    output_file_ = nullptr;
  }
  delete baseline_subtracted_;
  delete cropped_;
}

Bool_t WaveformProcessor::ProcessFile(const TString &filepath,
                                      const TString &output_name) {
  WaveformSource *source = CreateWaveformSource(config_.input_format);
  if (!source || !source->Open(filepath)) {
    delete source;
    return kFALSE;
  }
  Bool_t result = ProcessSource(*source, output_name);
  source->Close();
  delete source;
  return result;
}

Bool_t WaveformProcessor::ProcessSource(WaveformSource &source,
                                        const TString &output_name) {
  current_output_name_ = output_name;
  sample_waveforms_saved_ = 0;
  stats_ = ProcessingStats();
  if (!CheckSettings())
    return kFALSE;
  if (!OpenOutput(output_name))
    return kFALSE;

  WaveformRecord record;
  while (source.Next(record)) {
    if (config_.max_events > 0 && stats_.accepted >= config_.max_events)
      break;
    if (config_.channel >= 0 && record.channel != config_.channel)
      continue;
    stats_.total_processed++;
    ProcessRecord(record);
  }

  CloseOutput();
  if (config_.verbose)
    PrintAllStatistics();
  return kTRUE;
}

Bool_t WaveformProcessor::OpenOutput(const TString &output_name) {
  const TString base_dir = IO::GetRootFilesBaseDir();
  if (gSystem->AccessPathName(base_dir))
    gSystem->mkdir(base_dir, kTRUE);

  TString stats_path = base_dir + "/" + output_name + ".stats";
  std::ofstream(stats_path.Data(), std::ios::trunc);

  TString output_subpath = output_name + ".root";
  output_file_ = IO::OpenForWriting(output_subpath);
  if (!output_file_ || output_file_->IsZombie()) {
    std::cout << "ERROR: Could not create output file " << base_dir << "/"
              << output_subpath << std::endl;
    return kFALSE;
  }

  output_tree_ = new TTree("features", "Waveform Features");
  output_tree_->Branch("timestamp", &out_timestamp_, "timestamp/l");
  output_tree_->Branch("channel", &out_channel_, "channel/I");
  output_tree_->Branch("board", &out_board_, "board/I");
  output_tree_->Branch("digitizer_energy", &out_digitizer_energy_,
                       "digitizer_energy/D");
  output_tree_->Branch("flags", &out_flags_, "flags/i");
  output_tree_->Branch("trigger_position", &out_trigger_position_,
                       "trigger_position/I");
  output_tree_->Branch("pulse_height", &out_pulse_height_, "pulse_height/F");
  output_tree_->Branch("baseline_rms", &out_baseline_rms_, "baseline_rms/F");
  BookBranches(output_tree_);
  if (config_.store_waveforms)
    output_tree_->Branch("Samples", &cropped_);
  return kTRUE;
}

void WaveformProcessor::CloseOutput() {
  output_file_->cd();
  output_tree_->Write("", TObject::kOverwrite);
  output_file_->Close();
  delete output_file_;
  output_file_ = nullptr;
  output_tree_ = nullptr;
}

Bool_t WaveformProcessor::ProcessRecord(const WaveformRecord &record) {
  const TArrayI &samples = record.samples;
  Int_t n = samples.GetSize();
  if (n == 0)
    return kFALSE;

  Int_t raw_extremum = samples.At(0);
  for (Int_t i = 1; i < n; ++i) {
    if (config_.polarity == 1 && samples[i] > raw_extremum)
      raw_extremum = samples[i];
    if (config_.polarity != 1 && samples[i] < raw_extremum)
      raw_extremum = samples[i];
  }

  SubtractBaseline(samples);

  Float_t trigger_pos = FindTrigger(*baseline_subtracted_);
  if (trigger_pos < 0) {
    stats_.rejected_no_trigger++;
    return kFALSE;
  }
  if (trigger_pos < config_.pre_samples ||
      (baseline_subtracted_->GetSize() - trigger_pos) <= config_.post_samples) {
    stats_.rejected_insufficient_samples++;
    return kFALSE;
  }

  CropWaveform(Int_t(trigger_pos));
  const Float_t *cropped = cropped_->GetArray();
  out_pulse_height_ = *std::max_element(cropped, cropped + cropped_->GetSize());
  out_trigger_position_ = Int_t(FindTrigger(*cropped_));

  Int_t raw_pulse_height = std::abs(raw_extremum);
  if ((raw_pulse_height == config_.adc_saturation_code &&
       config_.polarity == 1) ||
      (raw_pulse_height == 0 && config_.polarity == -1)) {
    stats_.rejected_clipped++;
    return kFALSE;
  }

  if (!ExtractFeatures(*cropped_))
    return kFALSE;

  out_timestamp_ = record.timestamp;
  out_channel_ = record.channel;
  out_board_ = record.board;
  out_digitizer_energy_ = record.energy;
  out_flags_ = record.flags;
  out_baseline_rms_ = current_baseline_rms_valid_ ? current_baseline_rms_ : -1;

  if (sample_waveforms_saved_ < config_.sample_waveforms_to_save) {
    std::lock_guard<std::mutex> lock(canvas_mutex_);
    DrawSampleWaveform(*cropped_);
    sample_waveforms_saved_++;
  }

  output_tree_->Fill();
  stats_.accepted++;
  if (current_baseline_rms_valid_) {
    stats_.sum_baseline_rms_accepted += current_baseline_rms_;
    stats_.baseline_rms_count_accepted++;
  }
  return kTRUE;
}

void WaveformProcessor::SubtractBaseline(const TArrayI &samples) {
  Int_t n = samples.GetSize();

  Float_t baseline = 0;
  Int_t baseline_samples = TMath::Min(config_.num_samples_baseline, n);
  for (Int_t i = 0; i < baseline_samples; ++i)
    baseline += samples.GetAt(i);
  baseline /= baseline_samples;

  current_baseline_rms_valid_ = kFALSE;
  if (baseline_samples > 1) {
    Double_t sum_sq = 0.0;
    for (Int_t i = 0; i < baseline_samples; ++i) {
      Double_t d = samples.GetAt(i) - baseline;
      sum_sq += d * d;
    }
    Float_t rms = TMath::Sqrt(sum_sq / (baseline_samples - 1));
    stats_.sum_baseline_rms += rms;
    stats_.baseline_rms_count++;
    current_baseline_rms_ = rms;
    current_baseline_rms_valid_ = kTRUE;
  }

  baseline_subtracted_->Set(n);
  if (config_.polarity == -1) {
    for (Int_t i = 0; i < n; ++i)
      baseline_subtracted_->SetAt(baseline - samples.GetAt(i), i);
  } else {
    for (Int_t i = 0; i < n; ++i)
      baseline_subtracted_->SetAt(samples.GetAt(i) - baseline, i);
  }
}

Float_t WaveformProcessor::FindTrigger(const TArrayF &waveform) const {
  Int_t n = waveform.GetSize();
  const Float_t *arr = waveform.GetArray();
  Float_t peak_value = *std::max_element(arr, arr + n);
  // The level in double, then stored as float, as before the variants
  // existed, so the trigger sample is unchanged to the last bit.
  Double_t threshold = config_.trigger_threshold;
  Float_t trigger_level = peak_value * threshold;
  for (Int_t i = 0; i < n; ++i) {
    if (arr[i] >= trigger_level)
      return Float_t(i);
  }
  return -1.0;
}

void WaveformProcessor::CropWaveform(Int_t trigger_pos) {
  Int_t start = trigger_pos - config_.pre_samples;
  Int_t end = TMath::Min(trigger_pos + config_.post_samples,
                         baseline_subtracted_->GetSize());
  Int_t crop_size = end - start;
  cropped_->Set(crop_size);
  const Float_t *source = baseline_subtracted_->GetArray();
  for (Int_t i = 0; i < crop_size; ++i)
    cropped_->SetAt(source[start + i], i);
}

void WaveformProcessor::DrawSampleWaveform(const TArrayF &cropped) {
  Int_t n = cropped.GetSize();
  std::vector<Double_t> x(n), y(n);
  for (Int_t i = 0; i < n; ++i) {
    x[i] = i;
    y[i] = cropped[i];
  }
  PlottingUtils::SetStylePreferences(PlotSaveFormat::kPNG);
  TGraph *graph = new TGraph(n, x.data(), y.data());
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  PlottingUtils::ConfigureGraph(graph, kBlue + 1, ";Sample;Amplitude [ADC]");
  graph->Draw("AL");
  TString output_name = Form("%s_waveform_%04d", current_output_name_.Data(),
                             sample_waveforms_saved_);
  PlottingUtils::SaveFigure(canvas, output_name, "samplewaveforms",
                            PlotSaveOptions::kLINEAR);
  delete graph;
  delete canvas;
}

void WaveformProcessor::PrintAllStatistics() const {
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
  if (stats_.baseline_rms_count > 0) {
    std::cout << "Mean baseline RMS (all processed): "
              << stats_.sum_baseline_rms / stats_.baseline_rms_count
              << " ADC counts (over " << stats_.baseline_rms_count
              << " waveforms)" << std::endl;
  }
  if (stats_.baseline_rms_count_accepted > 0) {
    std::cout << "Mean baseline RMS (accepted only): "
              << stats_.sum_baseline_rms_accepted /
                     stats_.baseline_rms_count_accepted
              << " ADC counts (over " << stats_.baseline_rms_count_accepted
              << " waveforms)" << std::endl;
  }
  std::cout << std::endl;

  TString stats_path =
      IO::GetRootFilesBaseDir() + "/" + current_output_name_ + ".stats";
  std::ofstream stats_file(stats_path.Data(), std::ios::app);
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
    if (stats_.baseline_rms_count > 0) {
      stats_file << "Mean baseline RMS (all processed): "
                 << stats_.sum_baseline_rms / stats_.baseline_rms_count
                 << " ADC counts (over " << stats_.baseline_rms_count
                 << " waveforms)" << std::endl;
    }
    if (stats_.baseline_rms_count_accepted > 0) {
      stats_file << "Mean baseline RMS (accepted only): "
                 << stats_.sum_baseline_rms_accepted /
                        stats_.baseline_rms_count_accepted
                 << " ADC counts (over " << stats_.baseline_rms_count_accepted
                 << " waveforms)" << std::endl;
    }
    stats_file << std::endl;
  }
}

Bool_t WaveformProcessor::ProcessWithClone(const WaveformProcessor *prototype,
                                           const TString &filepath,
                                           const TString &output_name) {
  WaveformProcessor *processor = prototype->Clone();
  Bool_t result = processor->ProcessFile(filepath, output_name);
  delete processor;
  return result;
}

void WaveformProcessor::ProcessFilesParallel(
    const std::vector<TString> &filepaths,
    const std::vector<TString> &output_names, Int_t max_workers) const {
  ROOT::EnableThreadSafety();
  IO::SetThreadSafe(kTRUE);

  Int_t n_files = Int_t(filepaths.size());
  Int_t n_workers = max_workers > 0
                        ? max_workers
                        : Int_t(std::thread::hardware_concurrency());
  n_workers = TMath::Min(n_workers, n_files);

  std::cout << "Processing " << n_files << " files with " << n_workers
            << " workers." << std::endl;

  for (Int_t i = 0; i < n_files; i += n_workers) {
    std::vector<std::future<Bool_t>> futures;
    Int_t batch_end = TMath::Min(i + n_workers, n_files);
    for (Int_t j = i; j < batch_end; ++j) {
      futures.push_back(std::async(
          std::launch::async, &WaveformProcessor::ProcessWithClone, this,
          std::cref(filepaths[j]), std::cref(output_names[j])));
    }
    for (size_t j = 0; j < futures.size(); ++j) {
      Bool_t result = futures[j].get();
      std::cout << "Finished: " << output_names[i + j]
                << (result ? " [OK]" : " [FAILED]") << std::endl;
    }
  }
}

// ------------------------------------------------------------ pulse integral

PulseIntegralProcessor::PulseIntegralProcessor(
    const FileProcessingConfig &config, const PulseIntegralSettings &settings)
    : WaveformProcessor(config), settings_(settings), short_integral_(0),
      long_integral_(0) {}

WaveformProcessor *PulseIntegralProcessor::Clone() const {
  return new PulseIntegralProcessor(config_, settings_);
}

void PulseIntegralProcessor::BookBranches(TTree *tree) {
  tree->Branch("short_integral", &short_integral_, "short_integral/F");
  tree->Branch("long_integral", &long_integral_, "long_integral/F");
}

Bool_t PulseIntegralProcessor::ExtractFeatures(const TArrayF &cropped) {
  Int_t integration_start = config_.pre_samples - settings_.pre_gate;
  Int_t n = cropped.GetSize();
  const Float_t *arr = cropped.GetArray();

  short_integral_ = 0;
  long_integral_ = 0;
  Int_t negative_samples = 0;
  Int_t short_end = TMath::Min(integration_start + settings_.short_gate, n);
  Int_t long_end = TMath::Min(integration_start + settings_.long_gate, n);
  for (Int_t i = integration_start; i < long_end; ++i) {
    Float_t sample_value = arr[i];
    long_integral_ += sample_value;
    if (i < short_end)
      short_integral_ += sample_value;
    if (sample_value < 0)
      negative_samples++;
  }
  Float_t negative_fraction =
      Float_t(negative_samples) / Float_t(long_end - integration_start);

  if (negative_fraction > 0.50) {
    stats_.rejected_baseline++;
    return kFALSE;
  }
  if (long_integral_ <= 0) {
    stats_.rejected_negative_integral++;
    return kFALSE;
  }
  return kTRUE;
}

// --------------------------------------------------------------- pulse height

PulseHeightProcessor::PulseHeightProcessor(const FileProcessingConfig &config,
                                           const PulseHeightSettings &settings)
    : WaveformProcessor(config), settings_(settings), trapezoid_height_(0),
      rise_time_(-1) {}

WaveformProcessor *PulseHeightProcessor::Clone() const {
  return new PulseHeightProcessor(config_, settings_);
}

void PulseHeightProcessor::BookBranches(TTree *tree) {
  tree->Branch("trapezoid_height", &trapezoid_height_, "trapezoid_height/F");
  tree->Branch("rise_time", &rise_time_, "rise_time/F");
}

TArrayF PulseHeightProcessor::Deconvolve(const TArrayF &waveform,
                                         Double_t decay_time_samples) {
  Int_t n = waveform.GetSize();
  TArrayF step(n);
  if (decay_time_samples <= 0) {
    for (Int_t i = 0; i < n; ++i)
      step[i] = waveform[i];
    return step;
  }
  const Double_t decay_per_sample = TMath::Exp(-1.0 / decay_time_samples);
  Double_t value = 0;
  Double_t previous_input = 0;
  for (Int_t i = 0; i < n; ++i) {
    value += waveform[i] - decay_per_sample * previous_input;
    previous_input = waveform[i];
    step[i] = value;
  }
  return step;
}

TArrayF PulseHeightProcessor::Trapezoid(const TArrayF &step, Int_t rise_samples,
                                        Int_t flat_samples) {
  Int_t n = step.GetSize();
  TArrayF trapezoid(n);
  // cumulative[i] is the sum of the first i samples.
  std::vector<Double_t> cumulative(n + 1, 0.0);
  for (Int_t i = 0; i < n; ++i)
    cumulative[i + 1] = cumulative[i] + step[i];
  const Int_t k = rise_samples;
  const Int_t delay = rise_samples + flat_samples;
  // Samples before the start count as zero, the baseline, so a window that
  // reaches before it sums only the samples it has.
  for (Int_t i = 0; i < n; ++i) {
    Int_t late_first = TMath::Max(i - k + 1, 0);
    Int_t early_last = i - delay;
    Int_t early_first = TMath::Max(early_last - k + 1, 0);
    Double_t late = cumulative[i + 1] - cumulative[late_first];
    Double_t early = 0;
    if (early_last >= 0)
      early = cumulative[early_last + 1] - cumulative[early_first];
    trapezoid[i] = (late - early) / k;
  }
  return trapezoid;
}

Int_t PulseHeightProcessor::PeakingStart() const {
  return config_.pre_samples + settings_.rise_samples - 1 +
         Int_t(TMath::Nint(settings_.peaking_fraction *
                           settings_.flat_top_samples));
}

Bool_t PulseHeightProcessor::CheckSettings() const {
  const PulseHeightSettings &s = settings_;
  if (s.rise_samples < 1 || s.flat_top_samples < 0 ||
      s.peak_average_samples < 1 || s.peaking_fraction < 0 ||
      s.peaking_fraction > 1 || s.rise_plateau_samples < 1) {
    std::cout << "ERROR: PulseHeightSettings out of range: rise >= 1, flat "
                 "top >= 0, peak average >= 1, peaking fraction in [0, 1], "
                 "rise plateau >= 1"
              << std::endl;
    return kFALSE;
  }
  const Int_t crop_size = config_.pre_samples + config_.post_samples;
  const Int_t earliest =
      PeakingStart() - s.rise_samples - s.flat_top_samples - s.rise_samples + 1;
  const Int_t latest = PeakingStart() + s.peak_average_samples - 1;
  if (earliest < 0 || latest >= crop_size) {
    std::cout << "ERROR: the trapezoid at its peaking point needs samples "
              << earliest << " to " << latest << " of a crop of " << crop_size
              << " (pre_samples " << config_.pre_samples << ", post_samples "
              << config_.post_samples << ")" << std::endl;
    return kFALSE;
  }
  const Int_t plateau_end =
      config_.pre_samples + s.rise_plateau_start + s.rise_plateau_samples;
  if (plateau_end > crop_size || config_.pre_samples - s.rise_search_back < 1) {
    std::cout << "ERROR: the rise-time windows do not fit a crop of "
              << crop_size << " samples" << std::endl;
    return kFALSE;
  }
  return kTRUE;
}

Double_t PulseHeightProcessor::CrossingTime(const TArrayF &step, Int_t i,
                                            Double_t level) {
  if (step[i - 1] >= level || step[i] == step[i - 1])
    return i;
  return i - 1 + (level - step[i - 1]) / (step[i] - step[i - 1]);
}

Float_t PulseHeightProcessor::RiseTime(const TArrayF &step) const {
  const Int_t trigger = config_.pre_samples;
  const Int_t plateau_first = trigger + settings_.rise_plateau_start;
  Double_t plateau = 0;
  for (Int_t i = plateau_first;
       i < plateau_first + settings_.rise_plateau_samples; ++i)
    plateau += step[i];
  plateau /= settings_.rise_plateau_samples;
  if (!(plateau > 0))
    return -1;
  Double_t t10 = -1, t90 = -1;
  for (Int_t i = trigger - settings_.rise_search_back; i < plateau_first; ++i) {
    if (t10 < 0 && step[i] >= 0.1 * plateau)
      t10 = CrossingTime(step, i, 0.1 * plateau);
    if (t90 < 0 && step[i] >= 0.9 * plateau) {
      t90 = CrossingTime(step, i, 0.9 * plateau);
      break;
    }
  }
  if (t10 < 0 || t90 < 0)
    return -1;
  return Float_t(t90 - t10);
}

Bool_t PulseHeightProcessor::ExtractFeatures(const TArrayF &cropped) {
  const Int_t peaking_start = PeakingStart();
  if (peaking_start + settings_.peak_average_samples > cropped.GetSize()) {
    stats_.rejected_insufficient_samples++;
    return kFALSE;
  }
  TArrayF step = Deconvolve(cropped, settings_.decay_time_samples);
  TArrayF trapezoid =
      Trapezoid(step, settings_.rise_samples, settings_.flat_top_samples);
  Double_t sum = 0;
  for (Int_t i = 0; i < settings_.peak_average_samples; ++i)
    sum += trapezoid[peaking_start + i];
  trapezoid_height_ = sum / settings_.peak_average_samples;
  rise_time_ = RiseTime(step);
  return kTRUE;
}

void PulseHeightProcessor::DrawSampleWaveform(const TArrayF &cropped) {
  Int_t n = cropped.GetSize();
  TArrayF step = Deconvolve(cropped, settings_.decay_time_samples);
  TArrayF trapezoid =
      Trapezoid(step, settings_.rise_samples, settings_.flat_top_samples);
  std::vector<Double_t> x(n), y_input(n), y_step(n), y_trapezoid(n);
  Double_t y_lo = 0, y_hi = 0;
  for (Int_t i = 0; i < n; ++i) {
    x[i] = i;
    y_input[i] = cropped[i];
    y_step[i] = step[i];
    y_trapezoid[i] = trapezoid[i];
    y_lo = TMath::Min(y_lo, TMath::Min(y_input[i], y_step[i]));
    y_hi = TMath::Max(y_hi, TMath::Max(y_input[i], y_step[i]));
  }
  PlottingUtils::SetStylePreferences(PlotSaveFormat::kPNG);
  TCanvas *canvas = PlottingUtils::GetConfiguredCanvas(kFALSE);
  TGraph *input = new TGraph(n, x.data(), y_input.data());
  TGraph *deconvolved = new TGraph(n, x.data(), y_step.data());
  TGraph *filtered = new TGraph(n, x.data(), y_trapezoid.data());
  PlottingUtils::ConfigureGraph(input, kBlue + 1, ";Sample;Amplitude [ADC]");
  PlottingUtils::ConfigureGraph(deconvolved, kRed + 1, "");
  PlottingUtils::ConfigureGraph(filtered, kGreen + 2, "");
  input->SetMinimum(y_lo - 0.05 * (y_hi - y_lo));
  input->SetMaximum(y_hi + 0.35 * (y_hi - y_lo));
  input->Draw("AL");
  deconvolved->Draw("L SAME");
  filtered->Draw("L SAME");
  Double_t peaking_start = PeakingStart();
  Double_t peaking_end = peaking_start + settings_.peak_average_samples - 1;
  TLine *start_line =
      new TLine(peaking_start, y_lo, peaking_start, trapezoid_height_);
  TLine *end_line =
      new TLine(peaking_end, y_lo, peaking_end, trapezoid_height_);
  start_line->SetLineStyle(2);
  end_line->SetLineStyle(2);
  start_line->Draw();
  end_line->Draw();
  TLegend *legend = PlottingUtils::AddLegend(0.55, 0.9, 0.7, 0.9);
  legend->AddEntry(input, "Waveform", "l");
  legend->AddEntry(deconvolved, "Pole-zero deconvolved", "l");
  legend->AddEntry(filtered,
                   Form("Trapezoid, height %.1f ADC", trapezoid_height_), "l");
  legend->Draw();
  TString output_name = Form("%s_waveform_%04d", current_output_name_.Data(),
                             sample_waveforms_saved_);
  PlottingUtils::SaveFigure(canvas, output_name, "samplewaveforms",
                            PlotSaveOptions::kLINEAR);
  delete legend;
  delete start_line;
  delete end_line;
  delete input;
  delete deconvolved;
  delete filtered;
  delete canvas;
}
