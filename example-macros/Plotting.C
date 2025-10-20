#include "AnalysisUtils.h"
#include "HistogramUtils.h"
#include "PlottingUtils.h"
#include <TArrayS.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <vector>

void Plotting() {
  // Load the analysis library
  if (gSystem->Load("../lib/libWaveformAnalysis.so") < 0) {
    std::cerr << "Error: Could not load libWaveformAnalysis.so" << std::endl;
    return;
  }

  // Add include path
  gROOT->ProcessLine(".I ../include");

  // Set ROOT style
  PlottingUtils::SetROOTStyle();

  // Step 1: Plot Light Output Spectra (no cuts needed)
  std::cout << "Creating light output spectra plots..." << std::endl;
  std::map<Int_t, std::string> source_names = {
      {0, "Am-241"}, {1, "Cs-137"}, {2, "Na-22"}, {3, "Am-241 & Cs-137"}};

  HistogramUtils histMgr;
  std::vector<std::string> labels = {"Am-241", "Cs-137", "Na-22",
                                     "Am-241 & Cs-137"};

  for (size_t i = 0; i < labels.size(); ++i) {
    histMgr.AddSource(i, labels[i]);
  }

  histMgr.LoadFromFile("histograms.root"); // Loads EVERYTHING

  auto light_output_spectra = histMgr.GetAllLightOutputSpectra();

  PlottingUtils::PlotIndividualSpectra(light_output_spectra, source_names,
                                       "plots/individual_light_output");

  // Step 2: Try to load existing average waveforms, or calculate if not found
  std::cout << "Loading or calculating average waveforms..." << std::endl;

  // Create analysis utility instance
  AnalysisUtils analyzer;
  std::map<Int_t, std::vector<Float_t>> average_waveforms;

  // Try to load existing average waveforms first
  average_waveforms = analyzer.LoadAverageWaveforms("average_waveforms.root");

  if (average_waveforms.empty()) {
    std::cout
        << "No existing average waveforms found. Calculating from scratch..."
        << std::endl;

    // Define spectral cuts for waveform averaging
    AnalysisUtils::SpectralCuts cuts;
    cuts.min_light_output = 1000.0; // keVee
    cuts.max_light_output = 1750.0; // keVee
    cuts.max_events_per_source = 20000;

    // Calculate average waveforms (this handles file I/O internally)
    average_waveforms =
        analyzer.CalculateAverageWaveforms("processed_waveforms.root", cuts);
  } else {
    std::cout << "Loaded existing average waveforms from file." << std::endl;
  }

  // Step 3: Plot average waveforms
  if (!average_waveforms.empty()) {
    auto [f_alpha, f_gamma] =
        analyzer.ExtractWaveformsBySourceId(average_waveforms, 0, 2);

    if (!f_alpha.empty() && !f_gamma.empty()) {
      PlottingUtils::PlotAverageWaveforms(f_alpha, f_gamma,
                                          "plots/average_waveforms");

      // Step 4: Try to load existing SI weighting factor, or calculate if not
      // found
      std::cout << "Loading or calculating SI weighting factor..." << std::endl;

      std::vector<Float_t> weighting_factor =
          analyzer.LoadSIWeightingFactor("average_waveforms.root");

      if (weighting_factor.empty()) {
        std::cout << "No existing SI weighting factor found. Calculating..."
                  << std::endl;
        weighting_factor = analyzer.CalculateSIWeightingFactor(
            average_waveforms, 0, 2); // alpha_source_id=0, gamma_source_id=2
      } else {
        std::cout << "Loaded existing SI weighting factor from file."
                  << std::endl;
      }

      if (!weighting_factor.empty()) {
        PlottingUtils::PlotSIWeightingFactor(f_alpha, f_gamma, weighting_factor,
                                             "plots/weighting_factor");
        std::cout << "SI weighting factor plotting completed." << std::endl;
      } else {
        std::cout << "Could not load or calculate SI weighting factor!"
                  << std::endl;
      }
    } else {
      std::cout << "Could not extract alpha and gamma waveforms for plotting!"
                << std::endl;
    }
  } else {
    std::cout << "No waveforms available for plotting!" << std::endl;
  }

  PlottingUtils::PlotAverageWaveformsByLightOutput(
      "processed_waveforms.root", // input file
      "plots/avg_wf_by_light",    // output prefix
      250.0,                      // bin width (keVee)
      0.0,                        // min light output
      1750.0,                     // max light output
      50                          // minimum events per bin
  );

  std::cout << "Plotting macro completed." << std::endl;
}
