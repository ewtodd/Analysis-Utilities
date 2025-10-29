#include "CalibrationUtils.h"
#include "HistogramUtils.h"
#include "PlottingUtils.h"
#include <TSystem.h>

void Calibration() {
  // Load the analysis library
  if (gSystem->Load("../lib/libWaveformAnalysis.so") < 0) {
    std::cerr << "Error: Could not load libWaveformAnalysis.so" << std::endl;
    return;
  }

  gROOT->ProcessLine(".I ../include");

  // Step 1: Load existing histograms from file
  std::cout << "Loading histograms from file..." << std::endl;
  HistogramUtils *histMgr = new HistogramUtils();

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
  histMgr->SetConfig(histConfig);

  std::vector<std::string> labels = {"Am-241", "Cs-137", "Na-22",
                                     "Am-241 & Cs-137"};

  for (size_t i = 0; i < labels.size(); ++i) {
    histMgr->AddSource(i, labels[i]);
  }

  histMgr->LoadFromFile("histograms.root");

  std::cout << "ENERGY CALIBRATION" << std::endl;
  CalibrationUtils *calMgr = new CalibrationUtils();
  calMgr->SetIncludeZeroPoint(kTRUE);
  if (calMgr->LoadCalibration("calibration.root")) {
    std::cout << "Using existing calibration." << std::endl;
  } else {
    std::cout << "No existing calibration found. Creating new one..."
              << std::endl;

    calMgr->LoadIntegralSpectra(histMgr, "histograms.root");

    calMgr->AddCalibrationPeak("Am-241", 59.5409, 4100, 1500, 60000, 10, 3000,
                               7000, 0);
    calMgr->AddCalibrationPeak("Cs-137", 661.7, 39000, 1848, 18656, 25000,
                               34000, 45000, 1);
    calMgr->AddCalibrationPeak("Na-22_511", 511.0, 30251, 1635, 3176, 9059,
                               25000, 36000, 2);
    calMgr->AddCalibrationPeak("Na-22_1274", 1274.5, 76000, 2298, 135, 885,
                               71000, 82000, 2);

    if (!calMgr->FitAllPeaks()) {
      std::cout << "Error in peak fitting!" << std::endl;
      return;
    }

    if (!calMgr->CreateCalibrationCurve()) {
      std::cout << "Error creating quadratic calibration curve!" << std::endl;
      return;
    }

    if (!calMgr->CreateLinearCalibrationCurve()) {
      std::cout << "Error creating linear calibration curve!" << std::endl;
      return;
    }

    calMgr->SaveCalibration("calibration.root");
    calMgr->PrintResults();

    std::cout << "CREATING PLOTS" << std::endl;
    PlottingUtils::SetROOTStyle();

    PlottingUtils::PlotFittedPeaks(calMgr->GetPeaks(),
                                   histMgr->GetAllIntegralSpectra(),
                                   "plots/calibration_peaks");

    PlottingUtils::PlotCalibrationCurve(calMgr->GetCalibrationCurve(),
                                        calMgr->GetLinearCalibrationFunction(),
                                        "plots/calibration");

    std::cout << "Calibration and plotting complete." << std::endl;
  }

  std::cout << "Applying calibration to processed waveforms..." << std::endl;

  calMgr->ApplyCalibratedLightOutput("processed_waveforms.root",
                                     "histograms.root", histMgr);

  // Save calibrated spectra
  histMgr->SaveToFile("histograms.root");
  std::cout << "Calibrated spectra saved." << std::endl;

  std::cout << "Waveform data calibration complete!" << std::endl;
}
