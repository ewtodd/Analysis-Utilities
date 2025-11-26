#ifndef CONFIG_PARSER_H
#define CONFIG_PARSER_H

#include <TROOT.h>
#include <optional>
#include <string>
#include <toml++/toml.hpp>
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

enum class DetectorType { SCINTILLATOR, SEMICONDUCTOR };

class ConfigParser {
private:
  toml::table config_;
  DetectorType detector_type_;
  bool psd_capable_;

  void DetermineDetectorCapabilities();

public:
  ConfigParser() = default;

  bool Load(const std::string &filename);
  bool Validate();
  void Print() const;

  template <typename T> std::optional<T> Get(const std::string &key) const {
    auto keys = SplitKey(key);
    const toml::table *current = &config_;

    for (size_t i = 0; i < keys.size() - 1; i++) {
      auto node = current->get(keys[i]);
      if (!node || !node->is_table())
        return std::nullopt;
      current = node->as_table();
    }

    auto node = current->get(keys.back());
    if (!node)
      return std::nullopt;
    return node->value<T>();
  }

  template <typename T> T GetOr(const std::string &key, T default_value) const {
    return Get<T>(key).value_or(default_value);
  }

  std::vector<SourceConfig> GetSources() const;
  PSDSourceConfig GetPSDSources() const;
  DetectorType GetDetectorType() const { return detector_type_; }
  bool IsPSDCapable() const { return psd_capable_; }

  WaveformProcessingUtils *CreateWaveformProcessor() const;
  HistogramUtils *CreateHistogramManager() const;
  AnalysisUtils *CreateAnalysisUtils() const;
  CalibrationUtils *CreateCalibrationManager() const;

private:
  std::vector<std::string> SplitKey(const std::string &key) const;
};

#endif
