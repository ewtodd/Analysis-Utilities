#ifndef CALIBRATIONUTILS_H
#define CALIBRATIONUTILS_H

#include "HistogramUtils.hh"
#include <TCanvas.h>
#include <TF1.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <map>
#include <string>
#include <vector>

struct CalibrationPeak {
  // Assumes peak shape is Gaussian on top of a linear background
  std::string calibration_source_;
  Int_t source_id_;
  Double_t deposited_energy_keV_;

  Double_t guess_peak_mu_; // This is most likely pulse integral
                           // for scintillators, and pulse
                           // height for semiconductors
  Double_t guess_sigma_;
  Double_t guess_amplitude_;
  Double_t guess_bkg_const_;
  Double_t guess_bkg_slope_;
  Double_t fit_range_low_;
  Double_t fit_range_high_;

  Bool_t fit_successful_;
  Double_t fit_position_;
  Double_t fit_position_error_;
  Double_t fit_sigma_;
  Double_t fit_amplitude_;
  Double_t fit_bkg_const_;
  Double_t fit_bkg_slope_;
};

class CalibrationUtils {
private:
  std::vector<CalibrationPeak> peaks_;
  TGraphErrors *calibration_curve_;
  TF1 *calibration_function_;
  std::map<Int_t, TH1F *> measured_spectra_;
  Bool_t include_zero_point_;

public:
  CalibrationUtils();
  ~CalibrationUtils();
  void SetIncludeZeroPoint(Bool_t include_zero = kTRUE) {
    include_zero_point_ = include_zero;
  }
  void AddCalibrationPeak(const std::string &calibration_source,
                          Int_t source_id, Double_t deposited_energy_keV,
                          Double_t guess_peak_mu, Double_t guess_sigma,
                          Double_t guess_amplitude, Double_t guess_bkg_const,
                          Double_t guess_bkg_slope, Double_t fit_range_low,
                          Double_t fit_range_high);

  void LoadSpectra(HistogramUtils *histMgr, const std::string &filename);

  Bool_t FitAllPeaks();
  Bool_t FitSinglePeak(CalibrationPeak &peak, Int_t source_id);

  Bool_t CreateCalibrationCurve();
  Bool_t CreateLinearCalibrationCurve();

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
