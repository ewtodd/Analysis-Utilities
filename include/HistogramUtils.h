#ifndef HISTOGRAMUTILS_H
#define HISTOGRAMUTILS_H

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <map>
#include <string>
#include <vector>

struct HistogramConfig {
  // Light output spectra (keVee)
  Double_t light_output_min = 0;
  Double_t light_output_max = 2000;
  Double_t light_output_bin_width = 10;

  Int_t integral_bins = 300;
  Double_t integral_min = 0;
  Double_t integral_max = 100000;

  Int_t ph_bins = 200;
  Double_t ph_min = 0;
  Double_t ph_max = 4000;

  Int_t cc_bins = 100;
  Double_t cc_min = 0;
  Double_t cc_max = 1;

  Int_t si_bins = 100;
  Double_t si_min = 0;
  Double_t si_max = 1;

  Int_t GetLightOutputBins() const {
    return static_cast<Int_t>((light_output_max - light_output_min) /
                              light_output_bin_width);
  }
};

struct MeasurementStats {
  std::string run_id;
  Double_t live_time_seconds;
  Int_t output_counts;
  Double_t output_rate_cps;
};

class HistogramUtils {
private:
  HistogramConfig config_;
  std::map<Int_t, std::string> source_names_;

  std::map<Int_t, TH1F *> light_output_spectra_;
  std::map<Int_t, TH1F *> integral_spectra_;
  std::map<Int_t, TH1F *> pulse_height_spectra_;
  std::map<Int_t, TH1F *> charge_comparison_spectra_;
  std::map<Int_t, TH1F *> shape_indicator_spectra_;

  std::map<Int_t, MeasurementStats> source_stats_;
  Int_t background_source_id_;

  Bool_t ParseStatsFile(const std::string &filepath, Int_t source_id);
  std::string FindStatsFile(const std::string &directory);
  Double_t ParseTimeString(const std::string &time_str);
  void SubtractSpectrumType(Int_t source_id, const std::string &spectrum_type,
                            Double_t source_norm_factor,
                            Double_t bg_norm_factor);

public:
  HistogramUtils();
  ~HistogramUtils();

  void SetBackgroundSource(Int_t source_id) {
    background_source_id_ = source_id;
  }
  void LoadMeasurementStatistics(const std::vector<std::string> &directories,
                                 const std::vector<Int_t> &source_ids);
  TH1F *GetBackgroundSubtractedSpectrum(Int_t source_id,
                                        const std::string &spectrum_type) const;
  void ApplyBackgroundSubtraction();

  void SetConfig(const HistogramConfig &config);
  void AddSource(Int_t source_id, const std::string &name);

  void CreateAllHistograms();
  void CreateHistogramForSource(Int_t source_id);

  void FillFromTree(const std::string &filename,
                    const std::string &treename = "features");

  TH1F *GetLightOutputSpectrum(Int_t source_id) const;
  TH1F *GetIntegralSpectrum(Int_t source_id) const;
  TH1F *GetPulseHeightSpectrum(Int_t source_id) const;
  TH1F *GetChargeComparisonSpectrum(Int_t source_id) const;
  TH1F *GetShapeIndicatorSpectrum(Int_t source_id) const;
  const std::map<Int_t, std::string> GetSourceNames() const {
    return source_names_;
  };

  std::map<Int_t, TH1F *> GetAllLightOutputSpectra() const {
    return light_output_spectra_;
  }
  std::map<Int_t, TH1F *> GetAllIntegralSpectra() const {
    return integral_spectra_;
  }
  std::map<Int_t, TH1F *> GetAllChargeComparisonSpectra() const {
    return charge_comparison_spectra_;
  }
  std::map<Int_t, TH1F *> GetAllShapeIndicatorSpectra() const {
    return shape_indicator_spectra_;
  }
  void SaveToFile(const std::string &filename);
  Bool_t LoadFromFile(const std::string &filename);
};

#endif
