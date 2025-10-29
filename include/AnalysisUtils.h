#ifndef ANALYSIS_UTILS_H
#define ANALYSIS_UTILS_H

#include "PlottingUtils.h"
#include <TArrayS.h>
#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

class AnalysisUtils {
public:
  struct SpectralCuts {
    Float_t min_light_output = 0.0;
    Float_t max_light_output = 2000.0;
    Int_t max_events_per_source = 5000;
  };

  std::map<Int_t, std::vector<Float_t>>
  CalculateAverageWaveforms(const std::string &filename,
                            const SpectralCuts &cuts);

  void SaveAverageWaveforms(
      const std::map<Int_t, std::vector<Float_t>> &avgWaveforms,
      const std::string &filename = "average_waveforms.root");
  std::map<Int_t, std::vector<Float_t>>
  LoadAverageWaveforms(const std::string &filename = "average_waveforms.root");

  std::vector<Float_t> CalculateSIWeightingFactor(
      const std::map<Int_t, std::vector<Float_t>> &waveforms,
      Int_t alpha_source_id = 0, Int_t gamma_source_id = 2);

  void
  SaveSIWeightingFactor(const std::vector<Float_t> &SIWeightingFactor,
                        const std::string &filename = "average_waveforms.root");

  std::vector<Float_t>
  LoadSIWeightingFactor(const std::string &filename = "average_waveforms.root");

  Bool_t CalculateChargeComparisonPSD(
      Int_t short_gate, Int_t long_gate,
      const std::string &waveforms_file = "processed_waveforms.root",
      const std::string &histograms_file = "histograms.root",
      Bool_t force_recalculate = kFALSE);

  Bool_t CalculateSIPSD(
      const std::string &waveforms_file = "processed_waveforms.root",
      const std::string &histograms_file = "histograms.root",
      const std::string &average_waveforms_file = "average_waveforms.root",
      Bool_t force_recalculate = kFALSE);

  TH1F *GetFilteredPSDHistogram(
      const std::string &psd_type, // "charge_comparison_psd" or "si_psd"
      Float_t min_light_output = 500.0, Float_t max_light_output = 1250.0,
      Int_t source_id_calc = 0, // Am-241
      const std::string &waveforms_file = "processed_waveforms.root");

  std::pair<std::vector<Float_t>, std::vector<Float_t>>
  ExtractWaveformsBySourceId(
      const std::map<Int_t, std::vector<Float_t>> &waveforms,
      Int_t alpha_source_id = 0, Int_t gamma_source_id = 2);

  struct GateOptimizationResult {
    Int_t short_gate;
    Int_t long_gate;
    Double_t fom;
    Int_t alpha_events;
    Int_t gamma_events;
  };

  std::vector<GateOptimizationResult> OptimizeChargeComparisonGates(
      Float_t min_light_output = 400.0, Float_t max_light_output = 800.0,
      Int_t alpha_source_id = 0, Int_t gamma_source_id = 2,
      Int_t min_short_gate = 7, Int_t max_short_gate = 20,
      Int_t short_gate_step = 1, // Separate step for short gate
      Int_t min_long_gate = 110, Int_t max_long_gate = 200,
      Int_t long_gate_step = 10, // Separate step for long gate
      const std::string &waveforms_file = "processed_waveforms.root");

  TH1F *CalculateChargeComparisonPSDOnTheFly(
      Int_t short_gate, Int_t long_gate, Float_t min_light_output,
      Float_t max_light_output, Int_t source_id_calc,
      const std::string &waveforms_file = "processed_waveforms.root");
};

#endif
