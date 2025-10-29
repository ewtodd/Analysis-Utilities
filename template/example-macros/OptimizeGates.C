#include "AnalysisUtils.h"
#include "PlottingUtils.h"
#include <TSystem.h>
#include <iostream>

void OptimizeGates() {
  // Load the analysis library
  if (gSystem->Load("../lib/libWaveformAnalysis.so") < 0) {
    std::cerr << "Error: Could not load libWaveformAnalysis.so" << std::endl;
    return;
  }

  // Add include path
  gROOT->ProcessLine(".I ../include");

  // Set ROOT style
  PlottingUtils::SetROOTStyle();

  std::cout << "Starting gate optimization for Charge Comparison PSD..."
            << std::endl;

  // Create analysis utility instance
  AnalysisUtils analyzer;

  // Optimization parameters
  Float_t min_light = 1000.0; // keVee
  Float_t max_light = 1200.0; // keVee
  Int_t alpha_source_id = 0;  // Am-241
  Int_t gamma_source_id = 2;  // Na-22

  Int_t min_short_gate = 7;
  Int_t max_short_gate = 20;
  Int_t short_gate_step = 1;
  Int_t min_long_gate = 110;
  Int_t max_long_gate = 200;
  Int_t long_gate_step = 10;

  // Run optimization
  auto results = analyzer.OptimizeChargeComparisonGates(
      min_light, max_light, alpha_source_id, gamma_source_id, min_short_gate,
      max_short_gate, short_gate_step, min_long_gate, max_long_gate,
      long_gate_step);

  if (results.empty()) {
    std::cout << "No valid gate combinations found!" << std::endl;
    return;
  }

  // Print top 10 results
  std::cout << "\nTop 10 gate combinations:" << std::endl;
  std::cout << "Rank | Short | Long | FOM    | alpha Events | gamma Events"
            << std::endl;
  std::cout << "-----|-------|------|--------|----------|----------"
            << std::endl;

  Int_t max_results = std::min(10, (Int_t)results.size());
  for (Int_t i = 0; i < max_results; ++i) {
    const auto &result = results[i];
    printf("%4d | %5d | %4d | %6.3f | %8d | %8d\n", i + 1, result.short_gate,
           result.long_gate, result.fom, result.alpha_events,
           result.gamma_events);
  }

  // Create plot with optimal gates
  std::cout << "\nCreating PSD plot with optimal gates..." << std::endl;
  const auto &best = results[0];

  TH1F *best_alpha = analyzer.CalculateChargeComparisonPSDOnTheFly(
      best.short_gate, best.long_gate, min_light, max_light, alpha_source_id);
  TH1F *best_gamma = analyzer.CalculateChargeComparisonPSDOnTheFly(
      best.short_gate, best.long_gate, min_light, max_light, gamma_source_id);

  if (best_alpha && best_gamma) {
    std::string output_name = "plots/optimal_gates_" +
                              std::to_string(best.short_gate) + "_" +
                              std::to_string(best.long_gate);
    PlottingUtils::PlotPSDWithFOM(best_alpha, best_gamma, "PSP_{CC}",
                                  output_name);
  }

  std::cout << "\nOptimal gates: Short=" << best.short_gate
            << ", Long=" << best.long_gate << ", FOM=" << best.fom << std::endl;
  std::cout << "Gate optimization completed!" << std::endl;
}
