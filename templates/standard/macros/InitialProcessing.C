#include "CalibrationUtils.h"
#include "HistogramUtils.h"
#include "PlottingUtils.h"
#include "WaveformProcessingUtils.h"
#include <TSystem.h>

void InitialProcessing() {
  std::vector<std::string> filepaths = {
      "/home/e-work/LABDATA/ANSG/CeYAP/"
      "May12-14/YAP_am241_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/"
      "May12-14/YAP_cs137_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/"
      "May12-14/YAP_na22_new_settings",
      "/home/e-work/LABDATA/ANSG/CeYAP/"
      "May12-14/YAP_am241_cs137_new_settings",
  };

  std::vector<std::string> labels = {"Am-241", "Cs-137", "Na-22",
                                     "Am-241 & Cs-137"};

  std::cout << "Processing waveforms..." << std::endl;
  WaveformProcessingUtils *processor = new WaveformProcessingUtils();

  processor->SetPolarity("negative");
  processor->SetTriggerThreshold(0.15);
  processor->SetSampleWindows(17, 190); // pre_samples, post_samples
  processor->SetMaxEvents(-1);          // -1 for process all events
  processor->SetVerbose(kTRUE);
  processor->SetStoreWaveforms(kTRUE);
  for (size_t i = 0; i < labels.size(); ++i) {
    processor->AddSource(labels[i], i);
  }

  if (!processor->ProcessDataSets(filepaths, labels,
                                  "processed_waveforms.root")) {
    std::cout << "Error in waveform processing!" << std::endl;
    return;
  }

  // Step 2: Create histograms
  std::cout << "Creating histograms..." << std::endl;
  HistogramUtils *histMgr = new HistogramUtils();

  HistogramConfig histConfig;
  histConfig.light_output_min = 0;
  histConfig.light_output_max = 2000;

  histConfig.integral_min = 0;
  histConfig.integral_max = 120000;
  histConfig.integral_bins = 300;

  histConfig.ph_min = 0;
  histConfig.ph_max = 5000;
  histConfig.ph_bins = 200;

  histMgr->SetConfig(histConfig);

  for (size_t i = 0; i < labels.size(); ++i) {
    histMgr->AddSource(i, labels[i]);
  }
  histMgr->CreateAllHistograms();
  histMgr->FillFromTree("processed_waveforms.root");
  histMgr->SaveToFile("histograms.root");
  std::cout << "Initial processing complete." << std::endl;
}
