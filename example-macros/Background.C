#include "CalibrationUtils.h"
#include "HistogramUtils.h"
#include "PlottingUtils.h"

void Background() {
  HistogramUtils histMgr;

  // Configure histogram ranges (same as main macro)
  HistogramConfig histConfig;
  histConfig.light_output_min = 0;
  histConfig.light_output_max = 2000;
  histConfig.light_output_bin_width = 10;
  histConfig.integral_min = 0;
  histConfig.integral_max = 120000;
  histConfig.integral_bins = 1000;
  histConfig.ph_min = 0;
  histConfig.ph_max = 5000;
  histConfig.ph_bins = 200;
  histMgr.SetConfig(histConfig);
  // Data paths and labels
  std::vector<std::string> filepaths = {
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/YAP_am241_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/YAP_am241_flipped_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/YAP_cs137_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/YAP_na22_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/YAP_am241_cs137_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/YAP_am241_na22_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/May12-14/"
      "YAP_bkg_50lsb_weekend_new_settings"};

  std::vector<std::string> labels = {
      "Am-241",          "Am-241 Gamma Only", "Cs-137",    "Na-22",
      "Am-241 & Cs-137", "Am-241 & Na-22",    "Background"};

  for (size_t i = 0; i < labels.size(); ++i) {
    histMgr.AddSource(i, labels[i]);
  }

  // Load histograms from existing file
  if (!histMgr.LoadFromFile("histograms.root")) {
    std::cout << "Error: Could not load histograms.root" << std::endl;
    std::cout << "Run the full analysis first!" << std::endl;
    return;
  }

  // Step 3: Load measurement statistics and apply background subtraction
  std::cout << "=== STEP 3: Background Subtraction ===" << std::endl;

  // Load statistics from the same directories
  std::vector<Int_t> source_ids = {0, 1, 2, 3, 4, 5, 6};
  histMgr.LoadMeasurementStatistics(filepaths, source_ids);

  // Set background source (assuming "Background" is the last source, index 6)
  histMgr.SetBackgroundSource(6);

  // Add after SetBackgroundSource(6):
  if (histMgr.GetIntegralSpectrum(6) == nullptr) {
    std::cout << "ERROR: Background histogram not found!" << std::endl;
  } else {
    std::cout << "Background histogram found with "
              << histMgr.GetIntegralSpectrum(6)->GetEntries() << " entries"
              << std::endl;
  }
  // Add this BEFORE applying background subtraction:
  std::cout << "=== BEFORE BACKGROUND SUBTRACTION ===" << std::endl;
  histMgr.PrintStatistics();

  // Apply background subtraction
  histMgr.ApplyBackgroundSubtraction();

  // Add this AFTER background subtraction:
  std::cout << "=== AFTER BACKGROUND SUBTRACTION ===" << std::endl;
  histMgr.PrintStatistics();

  // Save background-subtracted histograms
  histMgr.SaveToFile("histograms_background_subtracted.root");

  std::cout << "=== Analysis Complete ===" << std::endl;
}
