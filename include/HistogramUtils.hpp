#ifndef HISTOGRAMUTILS_H
#define HISTOGRAMUTILS_H

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TTree.h>
#include <map>
#include <string>

struct HistogramConfig {
  Int_t calibrated_bin_width = 10;
  Double_t calibrated_min = 0;
  Double_t calibrated_max = 100000;
  Int_t calibrated_bins =
      GetBins(calibrated_min, calibrated_max, calibrated_bin_width);

  Int_t integral_bin_width = 10;
  Double_t integral_min = 0;
  Double_t integral_max = 100000;
  Int_t integral_bins = GetBins(integral_min, integral_max, integral_bin_width);

  Int_t ph_bin_width = 10;
  Double_t ph_min = 0;
  Double_t ph_max = 4000;
  Int_t ph_bins = GetBins(ph_min, ph_max, ph_bin_width);

  Int_t cc_bins = 100;
  Double_t cc_min = 0;
  Double_t cc_max = 1;

  Int_t GetBins(Int_t min, Int_t max, Int_t bin_width) const {
    return static_cast<Int_t>((max - min) / bin_width);
  }
};

class HistogramUtils {
private:
  HistogramConfig config_;
  std::map<Int_t, std::string> source_names_;

  std::map<Int_t, TH1F *> calibrated_spectra_;
  std::map<Int_t, TH1F *> integral_spectra_;
  std::map<Int_t, TH1F *> pulse_height_spectra_;
  std::map<Int_t, TH1F *> charge_comparison_spectra_;

public:
  HistogramUtils();
  ~HistogramUtils();

  void SetConfig(const HistogramConfig &config);
  void AddSource(Int_t source_id, const std::string &name);

  void CreateAllHistograms();
  void CreateHistogramForSource(Int_t source_id);

  void FillFromTree(const std::string &filename,
                    const std::string &treename = "features");

  TH1F *GetCalibratedSpectrum(Int_t source_id) const;
  TH1F *GetIntegralSpectrum(Int_t source_id) const;
  TH1F *GetPulseHeightSpectrum(Int_t source_id) const;
  TH1F *GetChargeComparisonSpectrum(Int_t source_id) const;
  const std::map<Int_t, std::string> GetSourceNames() const {
    return source_names_;
  };

  std::map<Int_t, TH1F *> GetAllIntegralSpectra() const {
    return integral_spectra_;
  }
  std::map<Int_t, TH1F *> GetAllPulseHeightSpectra() const {
    return pulse_height_spectra_;
  }
  std::map<Int_t, TH1F *> GetAllCalibratedSpectra() const {
    return calibrated_spectra_;
  }
  std::map<Int_t, TH1F *> GetAllChargeComparisonSpectra() const {
    return charge_comparison_spectra_;
  }

  void SaveToFile(const std::string &filename);
  Bool_t LoadFromFile(const std::string &filename);
};

#endif
