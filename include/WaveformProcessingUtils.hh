#ifndef WAVEFORMPROCESSOR_H
#define WAVEFORMPROCESSOR_H

#include <TArrayS.h>
#include <TFile.h>
#include <TTree.h>
#include <map>
#include <string>
#include <vector>

struct WaveformFeatures {
  Float_t pulse_height;
  Int_t peak_position;
  Int_t trigger_position;
  Float_t long_integral;
  Float_t negative_fraction;
  Int_t source_id;
  Bool_t passes_cuts;
};

struct ProcessingStats {
  Int_t total_processed = 0;
  Int_t accepted = 0;
  Int_t rejected_no_trigger = 0;
  Int_t rejected_insufficient_samples = 0;
  Int_t rejected_negative_integral = 0;
  Int_t rejected_baseline = 0;

  void Print(const std::string &source_name) const;
  Double_t GetAcceptanceRate() const;
};

class WaveformProcessingUtils {
private:
  // Analysis parameters - set by user
  std::string polarity_;
  Double_t trigger_threshold_;
  Int_t pre_samples_;
  Int_t post_samples_;
  Int_t max_events_;
  Bool_t verbose_;

  std::map<std::string, Int_t> source_map_;
  std::map<Int_t, ProcessingStats> stats_;

  TFile *output_file_;
  TTree *output_tree_;
  WaveformFeatures current_features_;
  Bool_t store_waveforms_;
  TArrayS *current_waveform_; // For ROOT tree branch
public:
  WaveformProcessingUtils();
  ~WaveformProcessingUtils();

  TTree *GetOutputTree() { return output_tree_; }
  TFile *GetOutputFile() { return output_file_; }
  // Parameter setters (called from analysis macro)
  void SetPolarity(const std::string &polarity) { polarity_ = polarity; }
  void SetTriggerThreshold(Double_t threshold) {
    trigger_threshold_ = threshold;
  }

  void SetSampleWindows(Int_t pre_samples, Int_t post_samples) {
    pre_samples_ = pre_samples;
    post_samples_ = post_samples;
  }
  void SetMaxEvents(Int_t max_events) { max_events_ = max_events; }
  void SetVerbose(Bool_t verbose) { verbose_ = verbose; }
  void SetStoreWaveforms(Bool_t store = kTRUE) { store_waveforms_ = store; }
  void AddSource(const std::string &name, Int_t id);

  // Main processing
  Bool_t ProcessDataSets(const std::vector<std::string> &filepaths,
                         const std::vector<std::string> &source_names,
                         const std::string &output_filename);

  Bool_t ProcessSingleFile(const std::string &filepath, Int_t source_id);
  Bool_t ProcessWaveform(const std::vector<Short_t> &samples, Int_t source_id);
  // Core waveform analysis
  std::vector<Float_t> SubtractBaseline(const std::vector<Short_t> &samples);
  Int_t FindTrigger(const std::vector<Float_t> &waveform);
  std::vector<Float_t> CropWaveform(const std::vector<Float_t> &waveform,
                                    Int_t trigger_pos);
  WaveformFeatures ExtractFeatures(const std::vector<Float_t> &cropped_wf,
                                   Int_t source_id);
  Bool_t ApplyQualityCuts(const WaveformFeatures &features);

  // Utilities
  std::vector<std::string>
  FindROOTFiles(const std::string &basepath,
                const std::string &pattern = "*DT5730*.root");

  // Statistics
  void PrintAllStatistics() const;
  ProcessingStats GetStats(Int_t source_id) const;
};

#endif
