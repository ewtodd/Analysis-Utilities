#include "AnalysisUtils.h"
#include "CalibrationUtils.h"
#include "HistogramUtils.h"
#include "PlottingUtils.h"
#include <TArrayS.h>
#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <vector>

void PSD() {
  // Load the analysis library
  if (gSystem->Load("../lib/libWaveformAnalysis.so") < 0) {
    std::cerr << "Error: Could not load libWaveformAnalysis.so" << std::endl;
    return;
  }

  // Add include path
  gROOT->ProcessLine(".I ../include");

  // Set ROOT style
  PlottingUtils::SetROOTStyle();

  std::cout << "Starting PSD analysis..." << std::endl;

  // Create analysis utility instance
  AnalysisUtils analyzer;

  // PSD parameters
  Int_t short_gate = 11;
  Int_t long_gate = 150;

  std::cout << "Using gates: short=" << short_gate << ", long=" << long_gate
            << std::endl;

  // Step 1: Calculate Charge Comparison PSD
  std::cout << "\nStep 1: Calculating Charge Comparison PSD..." << std::endl;
  Bool_t cc_success = analyzer.CalculateChargeComparisonPSD(
      short_gate, long_gate, "processed_waveforms.root", "histograms.root");
  if (!cc_success) {
    std::cout << "Error: Failed to calculate charge comparison PSD!"
              << std::endl;
    return;
  } else {
    std::cout << "Charge comparison PSD calculation completed successfully."
              << std::endl;
  }

  // Step 2: Calculate SI PSD
  std::cout << "\nStep 2: Calculating SI PSD..." << std::endl;
  Bool_t si_success = analyzer.CalculateSIPSD("processed_waveforms.root",
                                              "partial_histograms.root",
                                              "average_waveforms.root");

  if (!si_success) {
    std::cout << "Error: Failed to calculate SI PSD!" << std::endl;
    return;
  } else {
    std::cout << "SI PSD calculation completed successfully." << std::endl;
  }

  // Step 3: Load histograms and create plots
  std::cout << "\nStep 3: Creating PSD plots..." << std::endl;

  // Create HistogramUtils instance and load all histograms
  HistogramUtils histMgr;
  std::vector<std::string> labels = {"Am-241", "Cs-137", "Na-22",
                                     "Am-241 & Cs-137"};

  for (size_t i = 0; i < labels.size(); ++i) {
    histMgr.AddSource(i, labels[i]);
  }

  // Load all histograms (will load from the directories created by
  // AnalysisUtils)
  histMgr.LoadFromFile("histograms.root");

  // Get source names
  auto source_names = histMgr.GetSourceNames();

  auto cc_psd_spectra =
      histMgr.GetAllChargeComparisonSpectra(); // This maps to
                                               // charge_collection_spectra_
  if (!cc_psd_spectra.empty()) {
    PlottingUtils::PlotPSDSpectra(cc_psd_spectra, source_names,
                                  "plots/charge_comparison_psd", "PSP_{CC}");
  } else {
    std::cout << "Warning: No charge comparison PSD spectra found!"
              << std::endl;
  }

  // Plot SI PSD 1D spectra (using shape_indicator_spectra_)
  auto si_psd_spectra =
      histMgr.GetAllShapeIndicatorSpectra(); // This maps to
                                             // shape_indicator_spectra_
  if (!si_psd_spectra.empty()) {
    PlottingUtils::PlotPSDSpectra(si_psd_spectra, source_names, "plots/si_psd",
                                  "PSP_{SI}");
  } else {
    std::cout << "Warning: No SI PSD spectra found!" << std::endl;
  }

  // Create 2D plots
  PlottingUtils::PlotPSD2D("processed_waveforms.root", "charge_comparison_psd",
                           "plots/charge_comparison_psd_2d", "PSP_{CC}");

  PlottingUtils::PlotPSD2D("processed_waveforms.root", "si_psd",
                           "plots/si_psd_2d", "PSP_{SI}");

  std::cout << "\nPSD analysis completed!" << std::endl;
  std::cout << "Results saved to:" << std::endl;
  std::cout << "  - Branches: charge_comparison_psd and si_psd in "
               "processed_waveforms.root"
            << std::endl;
  std::cout << "  - Histograms: charge_comparison_psd_spectra and "
               "si_psd_spectra in histograms.root"
            << std::endl;
  std::cout << "  - Plots: Individual 1D and 2D PSD plots in plots/ directory "
            << std::endl;

  std::cout << "\nStep 4: Creating filtered PSD analysis with FOM..."
            << std::endl;

  // Define light output filter range
  Float_t min_light = 1000.0; // keVee
  Float_t max_light = 1200.0; // keVee

  TH1F *cc_alpha = analyzer.GetFilteredPSDHistogram("charge_comparison_psd",
                                                    min_light, max_light, 0);
  TH1F *cc_gamma = analyzer.GetFilteredPSDHistogram("charge_comparison_psd",
                                                    min_light, max_light, 2);

  if (cc_alpha && cc_gamma) {
    std::cout << "Alpha events: " << cc_alpha->GetEntries() << std::endl;
    std::cout << "Gamma events: " << cc_gamma->GetEntries() << std::endl;

    PlottingUtils::PlotPSDWithFOM(cc_alpha, cc_gamma, "PSP_{CC}",
                                  "plots/charge_comparison_psd_fom");
  } else {
    std::cout << "Failed to create CC PSD histograms" << std::endl;
  }

  // Analyze SI PSD
  std::cout << "\nAnalyzing SI PSD..." << std::endl;
  TH1F *si_alpha =
      analyzer.GetFilteredPSDHistogram("si_psd", min_light, max_light, 0);
  TH1F *si_gamma =
      analyzer.GetFilteredPSDHistogram("si_psd", min_light, max_light, 2);

  if (si_alpha && si_gamma) {
    std::cout << "Alpha events: " << si_alpha->GetEntries() << std::endl;
    std::cout << "Gamma events: " << si_gamma->GetEntries() << std::endl;

    PlottingUtils::PlotPSDWithFOM(si_alpha, si_gamma, "PSP_{SI}",
                                  "plots/si_psd_fom");
  } else {
    std::cout << "Failed to create SI PSD histograms" << std::endl;
  }

  std::cout << "FOM analysis completed!" << std::endl;

  // Optional: Print some statistics
  TFile *wf_file = TFile::Open("processed_waveforms.root", "READ");
  if (wf_file && !wf_file->IsZombie()) {
    TTree *tree = static_cast<TTree *>(wf_file->Get("features"));
    if (tree) {
      Long64_t n_entries = tree->GetEntries();
      std::cout << "Total events processed: " << n_entries << std::endl;

      // Quick validation - check if branches exist
      if (tree->GetBranch("charge_comparison_psd")) {
        std::cout << "✓ Charge comparison PSD branch exists" << std::endl;
      }
      if (tree->GetBranch("si_psd")) {
        std::cout << "✓ SI PSD branch exists" << std::endl;
      }
    }
    wf_file->Close();
  }

  std::cout << "All PSD plots created successfully!" << std::endl;
}
