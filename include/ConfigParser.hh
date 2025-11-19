#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <TROOT.h>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

class WaveformProcessingUtils;
class HistogramUtils;
class AnalysisUtils;
class CalibrationUtils;

struct SourceConfig {
  Int_t id;
  std::string name;
  std::string path;
};

struct PSDSourceConfig {
  Int_t first_source_id;
  Int_t second_source_id;
  std::string first_source_name;
  std::string second_source_name;
  std::string first_particle_type;
  std::string second_particle_type;
  bool valid;
};

enum class DetectorType { SCINTILLATOR, SEMICONDUCTOR, UNKNOWN };

class ConfigParser {
private:
  std::map<std::string, std::string> values_;
  DetectorType detector_type_;
  bool psd_capable_;

  std::string trim(const std::string &str) {
    size_t first = str.find_first_not_of(" \t\r\n");
    if (first == std::string::npos)
      return "";
    size_t last = str.find_last_not_of(" \t\r\n");
    return str.substr(first, last - first + 1);
  }

  void DetermineDetectorCapabilities() {
    std::string type_str = GetString("detector.type");
    std::transform(type_str.begin(), type_str.end(), type_str.begin(),
                   ::tolower);

    if (type_str == "scintillator") {
      detector_type_ = DetectorType::SCINTILLATOR;
      psd_capable_ = GetBool("detector.psd_capable", false);
    } else if (type_str == "semiconductor") {
      detector_type_ = DetectorType::SEMICONDUCTOR;
      psd_capable_ = false;
    } else {
      detector_type_ = DetectorType::UNKNOWN;
      psd_capable_ = false;
      if (!type_str.empty()) {
        std::cerr << "WARNING: Unknown detector type '" << type_str
                  << "'. Please specify 'scintillator' or 'semiconductor'."
                  << std::endl;
      }
    }
  }

public:
  ConfigParser() : detector_type_(DetectorType::UNKNOWN), psd_capable_(false) {}

  bool Load(const std::string &filename);
  bool Validate();
  void Print() const;

  std::string GetString(const std::string &key,
                        const std::string &default_val = "") const;
  Float_t GetFloat(const std::string &key, Float_t default_val = 0.0) const;
  Int_t GetInt(const std::string &key, Int_t default_val = 0) const;
  Bool_t GetBool(const std::string &key, Bool_t default_val = kFALSE) const;

  DetectorType GetDetectorType() const { return detector_type_; }
  PSDSourceConfig GetPSDSourceConfig() const;

  bool IsChargeComparisonEnabled() const {
    return GetBool("charge_comparison.enabled",
                   true); // Default: always enabled
  }
  bool IsPSDCapable() const { return psd_capable_; }
  bool IsScintillator() const {
    return detector_type_ == DetectorType::SCINTILLATOR;
  }
  bool IsSemiconductor() const {
    return detector_type_ == DetectorType::SEMICONDUCTOR;
  }

  std::vector<SourceConfig> GetSources() const;
  Int_t GetSourceID(const std::string &source_name) const;

  std::string GetProcessedWaveformsFile() const;
  std::string GetHistogramsFile() const;
  std::string GetAverageWaveformsFile() const;

  WaveformProcessingUtils *CreateWaveformProcessor() const;
  HistogramUtils *CreateHistogramManager() const;
  AnalysisUtils *CreateAnalysisUtils() const;
  CalibrationUtils *CreateCalibrationManager() const;
};

#endif
