#ifndef CALIBRATIONUTILS_H
#define CALIBRATIONUTILS_H

#include "HistogramUtils.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <map>
#include <string>
#include <vector>

struct CalibrationPeak {
  std::string isotope;
  Double_t deposited_energy_kev;
  Double_t expected_integral;
  Double_t expected_sigma;
  Double_t expected_amplitude;
  Double_t expected_background;
  Double_t fit_range_low;
  Double_t fit_range_high;
  Bool_t fit_successful;
  Double_t fitted_position;
  Double_t fitted_position_error;
  Double_t fitted_sigma;
  Double_t fitted_fwhm;
  Int_t source_id;
  Double_t fitted_amplitude; // Gaussian amplitude [0]
  Double_t fitted_bkg_const; // Background constant [3]
  Double_t fitted_bkg_slope; // Background slope [4]
  Double_t fitted_powerlaw_amp;
  Double_t fitted_powerlaw_exp;
};

class CalibrationUtils {
private:
  std::vector<CalibrationPeak> peaks_;
  TGraphErrors *calibration_curve_;
  TF1 *calibration_function_;
  TF1 *linear_calibration_function_;
  TGraphErrors *linear_calibration_curve_;
  std::map<Int_t, TH1F *> integral_spectra_;
  Bool_t include_zero_point_;

public:
  CalibrationUtils();
  ~CalibrationUtils();
  void SetIncludeZeroPoint(Bool_t include_zero = kTRUE) {
    include_zero_point_ = include_zero;
  }
  void AddCalibrationPeak(const std::string &isotope,
                          Double_t deposited_energy_kev,
                          Double_t expected_integral, Double_t expected_sigma,
                          Double_t expected_amplitude,
                          Double_t expected_background, Double_t fit_range_low,
                          Double_t fit_range_high, Int_t source_id);

  void LoadIntegralSpectra(HistogramUtils *histMgr,
                           const std::string &filename);

  Bool_t FitAllPeaks();
  Bool_t FitSinglePeak(CalibrationPeak &peak, Int_t source_id);

  Bool_t CreateCalibrationCurve();
  Bool_t CreateLinearCalibrationCurve();
  TF1 *GetLinearCalibrationFunction() const {
    return linear_calibration_function_;
  }
  TGraphErrors *GetLinearCalibrationCurve() const {
    return linear_calibration_curve_;
  }
  Double_t CalibrateToLightOutput(Double_t integral) const;
  TGraphErrors *GetCalibrationCurve() const { return calibration_curve_; }
  TF1 *GetCalibrationFunction() const { return calibration_function_; }
  std::vector<CalibrationPeak> GetPeaks() const { return peaks_; }

  Bool_t ApplyCalibratedLightOutput(const std::string &waveforms_file,
                                    const std::string &histograms_file,
                                    HistogramUtils *histMgr);
  void SaveCalibration(const std::string &filename);
  Bool_t LoadCalibration(const std::string &filename);

  void PrintResults() const;
};

#endif
