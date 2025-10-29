#include "CalibrationUtils.h"
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TMath.h>
#include <TSpline.h>
#include <algorithm>
#include <iostream>
#include <ostream>

CalibrationUtils::CalibrationUtils()
    : calibration_curve_(nullptr), calibration_function_(nullptr),
      linear_calibration_function_(nullptr), linear_calibration_curve_(nullptr),
      include_zero_point_(kTRUE) {}

CalibrationUtils::~CalibrationUtils() {
  calibration_curve_ = nullptr;
  calibration_function_ = nullptr;
  linear_calibration_function_ = nullptr;
  linear_calibration_curve_ = nullptr;
  integral_spectra_.clear();
}

void CalibrationUtils::AddCalibrationPeak(
    const std::string &isotope, Double_t deposited_energy_kev,
    Double_t expected_integral, Double_t expected_sigma,
    Double_t expected_amplitude, Double_t expected_background,
    Double_t fit_range_low, Double_t fit_range_high, Int_t source_id) {

  CalibrationPeak peak;
  peak.isotope = isotope;
  peak.deposited_energy_kev = deposited_energy_kev;
  peak.expected_integral = expected_integral;
  peak.expected_sigma = expected_sigma;
  peak.expected_amplitude = expected_amplitude;
  peak.expected_background = expected_background;
  peak.fit_range_low = fit_range_low;
  peak.fit_range_high = fit_range_high;
  peak.fit_successful = kFALSE;
  peak.source_id = source_id;

  peaks_.push_back(peak);
  std::cout << "Added calibration peak: " << isotope << " ("
            << deposited_energy_kev << " keV) for source_id " << source_id
            << std::endl;
}

void CalibrationUtils::LoadIntegralSpectra(HistogramUtils *histMgr,
                                           const std::string &filename) {
  // Use HistogramUtils to load the spectra
  histMgr->LoadFromFile(filename);
  // Copy the loaded spectra to our local map for fitting
  integral_spectra_ = histMgr->GetAllIntegralSpectra();

  std::cout << "CalibrationUtils loaded " << integral_spectra_.size()
            << " integral spectra from HistogramUtils." << std::endl;
}

Bool_t CalibrationUtils::FitAllPeaks() {
  if (peaks_.empty()) {
    std::cout << "Error: No calibration peaks defined!" << std::endl;
    return kFALSE;
  }

  if (integral_spectra_.empty()) {
    std::cout << "Error: No integral spectra loaded!" << std::endl;
    return kFALSE;
  }

  std::cout << "Fitting calibration peaks..." << std::endl;
  Bool_t all_fits_successful = kTRUE;

  for (size_t i = 0; i < peaks_.size(); ++i) {
    CalibrationPeak &peak = peaks_[i];
    std::cout << "Fitting peak " << (i + 1) << "/" << peaks_.size() << ": "
              << peak.isotope << " (" << peak.deposited_energy_kev << " keV)"
              << " using source_id " << peak.source_id << std::endl;

    // Use the specific source for this peak
    Bool_t fit_result = FitSinglePeak(peak, peak.source_id);

    if (fit_result) {
      std::cout << " SUCCESS: Position = " << peak.fitted_position << " ± "
                << peak.fitted_position_error << ", FWHM = " << peak.fitted_fwhm
                << std::endl;
    } else {
      std::cout << " FAILED: Could not fit peak" << std::endl;
      all_fits_successful = kFALSE;
    }
  }

  return all_fits_successful;
}

Bool_t CalibrationUtils::FitSinglePeak(CalibrationPeak &peak, Int_t source_id) {
  auto it = integral_spectra_.find(source_id);
  if (it == integral_spectra_.end()) {
    std::cout << " Error: No spectrum found for source_id " << source_id
              << std::endl;
    return kFALSE;
  }

  TH1F *hist = it->second;

  // Create fit function with Gaussian + background + power law tail
  TF1 *fit_func = nullptr;

  fit_func =
      new TF1("gauss_plus_bkg_tail",
              "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*pow(x, -[6])",
              peak.fit_range_low, peak.fit_range_high);

  // fit_func->SetParLimits(2, peak.expected_sigma - 500,
  //                       peak.expected_sigma + 500);

  if (source_id == 0) {
    fit_func->SetParLimits(1, peak.expected_integral - 100,
                           peak.expected_integral + 100);
    fit_func->SetParLimits(4, -10, 10);
    fit_func->SetParLimits(3, 0, 1e6);
    fit_func->SetParLimits(5, 0, 1e6);
    fit_func->SetParLimits(6, -10, 10);
    // fit_func->SetParLimits(0, peak.expected_amplitude - 10000,
    //                       peak.expected_amplitude + 10000);

  } else {
    fit_func->SetParLimits(1, peak.expected_integral - 500,
                           peak.expected_integral + 500);
    fit_func->SetParLimits(3, 0, peak.expected_background + 500);
    fit_func->SetParLimits(4, -10, 0);
    fit_func->SetParLimits(0, peak.expected_amplitude - 10000,
                           peak.expected_amplitude + 10000);
    fit_func->SetParLimits(5, -1e-5, 1e-5);
    fit_func->SetParLimits(6, -1e-5, 1e-5);
  }

  // Set initial parameters (now 7 parameters instead of 5)
  fit_func->SetParameters(peak.expected_amplitude, peak.expected_integral,
                          peak.expected_sigma, peak.expected_background, -0.01,
                          peak.expected_amplitude * 0.1, 2.0);

  // Parameter names
  fit_func->SetParName(0, "Amplitude");
  fit_func->SetParName(1, "Mean");
  fit_func->SetParName(2, "Sigma");
  fit_func->SetParName(3, "Bkg_Const");
  fit_func->SetParName(4, "Bkg_Slope");
  fit_func->SetParName(5, "PowerLaw_Amp");
  fit_func->SetParName(6, "PowerLaw_Exp");

  // Perform fit
  TFitResultPtr fit_result = hist->Fit(fit_func, "LSRN+");

  if (fit_result.Get() && fit_result->IsValid()) {
    peak.fitted_position = fit_func->GetParameter(1);
    peak.fitted_position_error = fit_func->GetParError(1);
    peak.fitted_sigma = fit_func->GetParameter(2);
    peak.fitted_fwhm = 2.355 * peak.fitted_sigma;
    peak.fitted_amplitude = fit_func->GetParameter(0);
    peak.fitted_bkg_const = fit_func->GetParameter(3);
    peak.fitted_bkg_slope = fit_func->GetParameter(4);

    // Store power law parameters (you may need to add these to CalibrationPeak)
    peak.fitted_powerlaw_amp = fit_func->GetParameter(5);
    peak.fitted_powerlaw_exp = fit_func->GetParameter(6);

    peak.fit_successful = kTRUE;
    delete fit_func;
    return kTRUE;
  }

  peak.fit_successful = kFALSE;
  delete fit_func;
  return kFALSE;
}
Bool_t CalibrationUtils::CreateCalibrationCurve() {
  // First, ensure linear calibration exists
  if (!CreateLinearCalibrationCurve()) {
    std::cout << "Error: Need linear calibration first!" << std::endl;
    return kFALSE;
  }

  // Get linear fit parameters
  Double_t linear_slope = linear_calibration_function_->GetParameter(1);
  Double_t linear_intercept = linear_calibration_function_->GetParameter(0);

  std::cout << "Creating quadratic calibration..." << std::endl;
  std::cout << "Linear parameters from high-energy fit:" << std::endl;
  std::cout << "  Slope: " << linear_slope << std::endl;
  std::cout << "  Intercept: " << linear_intercept << std::endl;

  // Prepare ALL calibration points (including low energy)
  std::vector<CalibrationPeak> successful_peaks;
  for (const auto &peak : peaks_) {
    if (peak.fit_successful) {
      successful_peaks.push_back(peak);
    }
  }

  // Calculate number of points (add 1 if including zero)
  Int_t n_points = successful_peaks.size();
  if (include_zero_point_) {
    n_points += 1;
    std::cout << "Including zero point (0,0) in fit" << std::endl;
  }

  Double_t *x_values = new Double_t[n_points];
  Double_t *x_errors = new Double_t[n_points];
  Double_t *y_values = new Double_t[n_points];
  Double_t *zero_errors = new Double_t[n_points];

  // Fill arrays starting with zero point if included
  Int_t point_idx = 0;
  if (include_zero_point_) {
    x_values[0] = 0.0;
    x_errors[0] = 0.0;
    y_values[0] = 0.0;
    point_idx = 1;
  }

  // Add calibration peaks
  for (Int_t i = 0; i < successful_peaks.size(); ++i) {
    x_values[point_idx] = successful_peaks[i].fitted_position;
    x_errors[point_idx] = successful_peaks[i].fitted_position_error;
    y_values[point_idx] = successful_peaks[i].deposited_energy_kev;
    zero_errors[point_idx] = 0;
    point_idx++;
  }

  // Create TGraph for fitting
  TGraph *fit_graph = new TGraph(n_points, x_values, y_values);

  // Create constrained quadratic function: E = a*I² + b*I + c
  calibration_function_ =
      new TF1("calibration_function", "pol2",
              0, // Start from 0 to include zero point
              *std::max_element(x_values, x_values + n_points));

  if (include_zero_point_) {
    // If including zero point, we can fix the intercept to 0
    // E = a*I² + b*I (no constant term)
    calibration_function_->SetParameter(0, 0.0);          // c = 0 (intercept)
    calibration_function_->SetParameter(1, linear_slope); // b (slope)
    calibration_function_->SetParameter(2, 0.0);          // a (quadratic term)

    // Fix the intercept to zero
    // calibration_function_->FixParameter(0, 0.0);

  } else {
    // Original behavior without zero constraint
    Double_t tolerance = 1e-5;
    calibration_function_->SetParameter(0, linear_intercept);
    calibration_function_->SetParameter(1, linear_slope);
    calibration_function_->SetParameter(2, 0.0);

    calibration_function_->SetParLimits(0, linear_intercept - tolerance,
                                        linear_intercept + tolerance);
    calibration_function_->SetParLimits(1, linear_slope - tolerance,
                                        linear_slope + tolerance);
    calibration_function_->SetParLimits(2, -1e6, 1e6);
  }

  // Perform constrained fit
  TFitResultPtr fit_result = fit_graph->Fit(calibration_function_, "QSR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    Double_t a = calibration_function_->GetParameter(2);
    Double_t b = calibration_function_->GetParameter(1);
    Double_t c = calibration_function_->GetParameter(0);

    std::cout << "Constrained quadratic calibration:" << std::endl;
    std::cout << "Quadratic coefficient: " << a << " ± " 
              << calibration_function_->GetParError(2) << std::endl;
    std::cout << "Chi2/NDF: " << calibration_function_->GetChisquare() << "/"
              << calibration_function_->GetNDF() << " = "
              << calibration_function_->GetChisquare() /
                     calibration_function_->GetNDF()
              << std::endl;

    calibration_curve_ =
        new TGraphErrors(n_points, x_values, y_values, x_errors, zero_errors);
    calibration_curve_->SetName("calibration_curve");

    std::string title =
        include_zero_point_
            ? "Constrained Quadratic Calibration (with origin);Pulse Integral "
              "[a.u.];Deposited Energy [keV]"
            : "Constrained Quadratic Calibration;Pulse Integral "
              "[a.u.];Deposited Energy [keV]";
    calibration_curve_->SetTitle(title.c_str());
    calibration_curve_->GetListOfFunctions()->Add(calibration_function_);

    delete fit_graph;
    delete[] x_values;
    delete[] x_errors;
    delete[] y_values;
    delete[] zero_errors;

    std::cout << "Constrained quadratic calibration created successfully!"
              << std::endl;
    return kTRUE;
  }

  // Cleanup on failure
  delete fit_graph;
  delete[] x_values;
  delete[] x_errors;
  delete[] y_values;
  return kFALSE;
}
Bool_t CalibrationUtils::CreateLinearCalibrationCurve() {
  // Filter for high energy points only (exclude 59.5 keV Am-241)
  std::vector<CalibrationPeak> high_energy_peaks;
  for (const auto &peak : peaks_) {
    if (peak.fit_successful &&
        peak.deposited_energy_kev > 400) { // Exclude Am-241 at 59.5 keV
      high_energy_peaks.push_back(peak);
    }
  }

  if (high_energy_peaks.size() < 2) {
    std::cout
        << "Error: Need at least 2 high energy points for linear calibration!"
        << std::endl;
    return kFALSE;
  }

  std::cout << "Creating linear calibration curve, ignoring low light output "
               "points..."
            << std::endl;
  std::cout << "Using " << high_energy_peaks.size()
            << " high energy calibration points:" << std::endl;

  // Prepare data for TGraph
  Int_t n_points = high_energy_peaks.size();
  Double_t *x_values = new Double_t[n_points];
  Double_t *x_errors = new Double_t[n_points];
  Double_t *y_values = new Double_t[n_points];

  for (Int_t i = 0; i < n_points; ++i) {
    x_values[i] = high_energy_peaks[i].fitted_position;
    x_errors[i] = high_energy_peaks[i].fitted_position_error;
    y_values[i] = high_energy_peaks[i].deposited_energy_kev;
    std::cout << "  Point " << (i + 1) << ": " << high_energy_peaks[i].isotope
              << " - Integral: " << x_values[i] << ", Energy: " << y_values[i]
              << " keV" << std::endl;
  }

  // Create TGraph for fitting
  TGraph *fit_graph = new TGraph(n_points, x_values, y_values);

  // Create LINEAR calibration function with extended range
  Double_t max_x = *std::max_element(x_values, x_values + n_points);

  linear_calibration_function_ =
      new TF1("linear_calibration_function", "pol1", 0, max_x + 1000);

  // Perform linear fit
  TFitResultPtr fit_result =
      fit_graph->Fit(linear_calibration_function_, "QSR+");

  if (fit_result.Get() && fit_result->IsValid()) {
    Double_t slope = linear_calibration_function_->GetParameter(1);
    Double_t intercept = linear_calibration_function_->GetParameter(0);

    std::cout << "Linear calibration function: E = " << slope << " * I + "
              << intercept << std::endl;
    std::cout << "Chi2/NDF: " << linear_calibration_function_->GetChisquare()
              << "/" << linear_calibration_function_->GetNDF() << " = "
              << linear_calibration_function_->GetChisquare() /
                     linear_calibration_function_->GetNDF()
              << std::endl;

    // Create TGraphErrors for plotting (high energy points only)
    Double_t *zero_errors = new Double_t[n_points];
    for (Int_t i = 0; i < n_points; ++i) {
      zero_errors[i] = 0;
    }

    linear_calibration_curve_ =
        new TGraphErrors(n_points, x_values, y_values, x_errors, zero_errors);
    linear_calibration_curve_->SetName("linear_calibration_curve");
    linear_calibration_curve_->SetTitle(
        "Linear Energy Calibration;Pulse Integral [a.u.];Deposited Energy "
        "[keV]");

    delete fit_graph;
    delete[] x_values;
    delete[] x_errors;
    delete[] y_values;
    delete[] zero_errors;

    std::cout << "Linear calibration curve created successfully!" << std::endl;
    return kTRUE;
  } else {
    std::cout << "Error: Linear calibration fit failed!" << std::endl;
    delete fit_graph;
    delete[] x_values;
    delete[] x_errors;
    delete[] y_values;
    return kFALSE;
  }
}
Double_t CalibrationUtils::CalibrateToLightOutput(Double_t integral) const {
  if (!calibration_function_) {
    std::cout << "Warning: No calibration function available!" << std::endl;
    return -1;
  }

  return calibration_function_->Eval(integral);
}

void CalibrationUtils::SaveCalibration(const std::string &filename) {
  TFile *file = new TFile(filename.c_str(), "RECREATE");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not create calibration file " << filename
              << std::endl;
    return;
  }

  std::cout << "Saving calibration to " << filename << "..." << std::endl;

  // Save calibration curve and function
  if (calibration_curve_) {
    calibration_curve_->Write("calibration_curve");
  }
  if (calibration_function_) {
    calibration_function_->Write("calibration_function");
  }

  // Save peak information as histograms (for easy retrieval)
  TH1F *peak_positions = new TH1F("peak_positions", "Fitted Peak Positions",
                                  peaks_.size(), 0, peaks_.size());
  TH1F *peak_energies = new TH1F("peak_energies", "Peak Energies",
                                 peaks_.size(), 0, peaks_.size());

  for (size_t i = 0; i < peaks_.size(); ++i) {
    peak_positions->SetBinContent(i + 1, peaks_[i].fitted_position);
    peak_positions->SetBinError(i + 1, peaks_[i].fitted_position_error);
    peak_energies->SetBinContent(i + 1, peaks_[i].deposited_energy_kev);

    // Set bin labels with isotope names
    peak_positions->GetXaxis()->SetBinLabel(i + 1, peaks_[i].isotope.c_str());
    peak_energies->GetXaxis()->SetBinLabel(i + 1, peaks_[i].isotope.c_str());
  }

  peak_positions->Write();
  peak_energies->Write();

  delete peak_positions;
  delete peak_energies;

  file->Close();
  delete file;

  std::cout << "Calibration saved successfully." << std::endl;
}

Bool_t CalibrationUtils::LoadCalibration(const std::string &filename) {
  TFile *file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open calibration file " << filename
              << std::endl;
    return kFALSE;
  }

  std::cout << "Loading calibration from " << filename << "..." << std::endl;

  // Load calibration curve
  TGraphErrors *curve =
      static_cast<TGraphErrors *>(file->Get("calibration_curve"));
  if (curve) {
    if (calibration_curve_)
      delete calibration_curve_;
    calibration_curve_ = static_cast<TGraphErrors *>(curve->Clone());
  }

  // Load calibration function
  TF1 *function = static_cast<TF1 *>(file->Get("calibration_function"));
  if (function) {
    if (calibration_function_)
      delete calibration_function_;
    calibration_function_ = static_cast<TF1 *>(function->Clone());
  }

  file->Close();

  if (calibration_curve_ && calibration_function_) {
    std::cout << "Calibration loaded successfully." << std::endl;
    return kTRUE;
  } else {
    std::cout << "Error: Could not load calibration data." << std::endl;
    return kFALSE;
  }
}

void CalibrationUtils::PrintResults() const {
  std::cout << "Calibration results:" << std::endl;

  std::cout << "Peak fitting results:" << std::endl;
  for (size_t i = 0; i < peaks_.size(); ++i) {
    const CalibrationPeak &peak = peaks_[i];
    std::cout << "  " << peak.isotope << " (" << peak.deposited_energy_kev
              << " keV): ";

    if (peak.fit_successful) {
      std::cout << "Position = " << peak.fitted_position << " ± " 
                << peak.fitted_position_error << " a.u." << std::endl;
    } else {
      std::cout << "FIT FAILED" << std::endl;
    }
  }

  if (calibration_function_) {
    Double_t slope = calibration_function_->GetParameter(1);
    Double_t intercept = calibration_function_->GetParameter(0);

    std::cout << "\nCalibration function: E = " << slope << " * I + "
              << intercept << std::endl;
    std::cout << "Conversion factor: " << slope << " keV per unit integral"
              << std::endl;
  }
}

Bool_t
CalibrationUtils::ApplyCalibratedLightOutput(const std::string &waveforms_file,
                                             const std::string &histograms_file,
                                             HistogramUtils *histMgr) {

  if (!calibration_function_) {
    std::cout << "Error: No calibration function loaded!" << std::endl;
    return kFALSE;
  }

  // Open input files
  TFile *wffile = TFile::Open(waveforms_file.c_str(), "UPDATE");
  TFile *histfile = TFile::Open(histograms_file.c_str(), "UPDATE");

  if (!wffile || wffile->IsZombie()) {
    std::cout << "Error: Could not open " << waveforms_file << std::endl;
    return kFALSE;
  }
  if (!histfile || histfile->IsZombie()) {
    std::cout << "Error: Could not open " << histograms_file << std::endl;
    return kFALSE;
  }

  // Get the tree and check for light output branch
  TTree *wftree = static_cast<TTree *>(wffile->Get("features"));
  if (!wftree) {
    std::cout << "Error: No features tree found" << std::endl;
    wffile->Close();
    histfile->Close();
    return kFALSE;
  }

  // Check if light_output_keVee branch already exists
  TBranch *light_branch = wftree->GetBranch("light_output_keVee");
  if (!light_branch) {
    // Create the branch first
    std::cout << "Creating light_output_keVee branch..." << std::endl;
    Float_t long_integral;
    Float_t light_output_keVee;

    wftree->SetBranchAddress("long_integral", &long_integral);
    wftree->Branch("light_output_keVee", &light_output_keVee,
                   "light_output_keVee/F");

    Long64_t n_entries = wftree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
      wftree->GetEntry(i);
      light_output_keVee = CalibrateToLightOutput(long_integral);
      wftree->GetBranch("light_output_keVee")->Fill();
    }
    wftree->AutoSave("SaveSelf");
    std::cout << "Branch created and saved." << std::endl;
  }

  bool histograms_already_filled = false;
  for (const auto &source_pair : histMgr->GetSourceNames()) {
    Int_t sid = source_pair.first;
    TH1F *hist = histMgr->GetLightOutputSpectrum(sid);
    if (hist && hist->GetEntries() > 0) {
      histograms_already_filled = true;
      std::cout << "Light output histogram for source " << sid
                << " already has " << hist->GetEntries() << " entries."
                << std::endl;
      break;
    }
  }

  if (histograms_already_filled) {
    std::cout << "Light output histograms already filled. Skipping histogram "
                 "filling..."
              << std::endl;
  } else {
    // Now read the calibrated data and fill histograms
    std::cout << "Filling light output histograms..." << std::endl;

    Float_t light_output_keVee;
    Int_t source_id;

    wftree->SetBranchAddress("light_output_keVee", &light_output_keVee);
    wftree->SetBranchAddress("source_id", &source_id);

    Long64_t n_entries = wftree->GetEntries();

    for (Long64_t i = 0; i < n_entries; ++i) {
      wftree->GetEntry(i);

      // Get the histogram for this source_id and fill it
      TH1F *hist = histMgr->GetLightOutputSpectrum(source_id);
      if (hist) {
        hist->Fill(light_output_keVee);
      }
    }

    // Navigate to light output directory and save updated histograms
    TDirectory *light_output_dir =
        histfile->GetDirectory("light_output_spectra");
    if (!light_output_dir) {
      light_output_dir = histfile->mkdir("light_output_spectra");
    }

    light_output_dir->cd();

    // Save all light output histograms
    for (const auto &source_pair : histMgr->GetSourceNames()) {
      Int_t sid = source_pair.first;
      TH1F *hist = histMgr->GetLightOutputSpectrum(sid);
      if (hist) {
        hist->Write("", TObject::kOverwrite);
      }
    }

    std::cout << "Light output histograms filled and saved." << std::endl;
  }

  wffile->Close();
  histfile->Close();

  return kTRUE;
}
