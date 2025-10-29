#ifndef PLOTTINGUTILS_H
#define PLOTTINGUTILS_H

#include "CalibrationUtils.h"
#include <TCanvas.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <map>
#include <set>
#include <string>
#include <vector>

class PlottingUtils {
public:
  // Style configuration
  static void SetROOTStyle();
  static void ConfigureHistogram(TH1 *hist, Int_t color,
                                 const std::string &title = "");
  static void ConfigureCanvas(TCanvas *canvas, Bool_t logy = kFALSE);

  static void
  PlotIndividualSpectra(const std::map<Int_t, TH1F *> &spectra,
                        const std::map<Int_t, std::string> &source_names,
                        const std::string &output_prefix = "spectrum");

  static void
  PlotCalibrationCurve(TGraphErrors *calibration_curve, TF1 *linear_function,
                       const std::string &output_name = "calibration_curve");

  static void
  PlotFittedPeaks(const std::vector<CalibrationPeak> &peaks,
                  const std::map<Int_t, TH1F *> &spectra,
                  const std::string &output_prefix = "fitted_peaks");

  static void PlotSinglePeakFit(const CalibrationPeak &peak, TH1F *histogram,
                                const std::string &output_filename);
  static void
  PlotAverageWaveforms(const std::map<Int_t, std::vector<Float_t>> &waveforms,
                       const std::map<Int_t, std::string> &source_names,
                       const std::string &output_name = "average_waveforms");

  static void
  PlotAverageWaveforms(const std::vector<Float_t> &f_alpha,
                       const std::vector<Float_t> &f_gamma,
                       const std::string &output_name = ";average_waveforms");
  static void
  PlotSIWeightingFactor(const std::vector<Float_t> &f_alpha,
                        const std::vector<Float_t> &f_gamma,
                        const std::vector<Float_t> &weighting_factor,
                        const std::string &output_name = "weighting_factor");

  static void PlotPSDSpectra(const std::map<Int_t, TH1F *> &psd_spectra,
                             const std::map<Int_t, std::string> &source_names,
                             const std::string &output_prefix = "psd_spectrum",
                             const std::string &psd_type = "PSD");

  static void
  PlotPSD2D(const std::string &waveforms_file = "processed_waveforms.root",
            const std::string &psd_branch = "charge_comparison_psd",
            const std::string &output_name = "psd_2d",
            const std::string &psd_label = "Charge Comparison PSD");
  static void PlotAverageWaveformsByLightOutput(
      const std::string &waveforms_file = "processed_waveforms.root",
      const std::string &output_prefix = "avg_waveforms_by_light",
      Float_t bin_width = 250.0, Float_t min_light = 0.0,
      Float_t max_light = 1750.0, Int_t min_events_per_bin = 50);

  // Utility
  static TLegend *CreateLegend(Double_t x1 = 0.7, Double_t y1 = 0.7,
                               Double_t x2 = 0.9, Double_t y2 = 0.9);
  static void AddSubplotLabel(const std::string &label, Double_t x = 0.9,
                              Double_t y = 0.85);
  static std::string CleanSourceName(const std::string &source_name);
  static Double_t CalculateFigureOfMerit(TH1F *hist_alpha, TH1F *hist_gamma,
                                         Double_t &separation,
                                         Double_t &fwhm_alpha,
                                         Double_t &fwhm_gamma);

  static void PlotPSDWithFOM(TH1F *hist_alpha, TH1F *hist_gamma,
                             const std::string &psd_type = "PSP_{CC}",
                             const std::string &output_name = "psd_with_fom");

  static void
  PlotFilteredPSDHistograms(TH1F *hist_alpha, TH1F *hist_gamma,
                            const std::string &psd_type = "PSD",
                            const std::string &output_name = "filtered_psd");

  // Colors
  static std::vector<Int_t> GetDefaultColors();
  static Int_t GetSourceColor(Int_t source_id);
};

#endif
