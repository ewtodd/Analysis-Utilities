#ifndef WAVEFORMPROCESSOR_H
#define WAVEFORMPROCESSOR_H

#include <TArrayS.h>
#include <TFile.h>
#include <TROOT.h>
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
  Bool_t passes_cuts;
  ULong64_t timestamp;
};

struct ProcessingStats {
  Int_t total_processed = 0;
  Int_t accepted = 0;
  Int_t rejected_no_trigger = 0;
  Int_t rejected_insufficient_samples = 0;
  Int_t rejected_negative_integral = 0;
  Int_t rejected_baseline = 0;
};

class WaveformProcessingUtils {
private:
  // Analysis parameters - set by user
  Int_t polarity_;
  Double_t trigger_threshold_;
  Int_t pre_samples_;
  Int_t post_samples_;
  Int_t short_gate_;
  Int_t long_gate_;
  Int_t max_events_;
  Bool_t verbose_;

  ProcessingStats stats_;

  TFile *output_file_;
  TTree *output_tree_;
  WaveformFeatures current_features_;
  Bool_t store_waveforms_;
  TArrayS *current_waveform_;
  ULong64_t current_timestamp_;

public:
  WaveformProcessingUtils();
  ~WaveformProcessingUtils();

  TTree *GetOutputTree() { return output_tree_; }
  TFile *GetOutputFile() { return output_file_; }

  void SetPolarity(const Int_t polarity) { polarity_ = polarity; }
  void SetTriggerThreshold(Double_t threshold) {
    trigger_threshold_ = threshold;
  }
  void SetSampleWindows(Int_t pre_samples, Int_t post_samples) {
    pre_samples_ = pre_samples;
    post_samples_ = post_samples;
  }
  void SetPSDGates(Int_t short_gate, Int_t long_gate) {
    short_gate_ = short_gate;
    long_gate_ = long_gate;
  }
  void SetMaxEvents(Int_t max_events) { max_events_ = max_events; }
  void SetVerbose(Bool_t verbose) { verbose_ = verbose; }
  void SetStoreWaveforms(Bool_t store = kTRUE) { store_waveforms_ = store; }

  // Main processing
  Bool_t ProcessFile(const TString filepath, const TString output_filename);

  Bool_t ProcessWaveform(const std::vector<Short_t> &samples);

  // Core waveform analysis
  std::vector<Float_t> SubtractBaseline(const std::vector<Short_t> &samples);
  Int_t FindTrigger(const std::vector<Float_t> &waveform);
  std::vector<Float_t> CropWaveform(const std::vector<Float_t> &waveform,
                                    Int_t trigger_pos);
  WaveformFeatures ExtractFeatures(const std::vector<Float_t> &cropped_wf);
  Bool_t ApplyQualityCuts(const WaveformFeatures &features);

  // Statistics
  void PrintAllStatistics() const;
  ProcessingStats GetStats() const { return stats_; };
};

#endif
