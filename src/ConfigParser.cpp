#include "ConfigParser.hpp"
#include "AnalysisUtils.hpp"
#include "CalibrationUtils.hpp"
#include "HistogramUtils.hpp"
#include "WaveformProcessingUtils.hpp"
#include <iostream>

bool ConfigParser::Load(const std::string &filename) {
  try {
    config_ = toml::parse_file(filename);
    DetermineDetectorCapabilities();
    return true;
  } catch (const toml::parse_error &err) {
    std::cerr << "TOML parsing error:\n" << err << std::endl;
    return false;
  }
}

std::vector<std::string> ConfigParser::SplitKey(const std::string &key) const {
  std::vector<std::string> result;
  std::stringstream ss(key);
  std::string item;
  while (std::getline(ss, item, '.')) {
    result.push_back(item);
  }
  return result;
}

std::vector<SourceConfig> ConfigParser::GetSources() const {
  std::vector<SourceConfig> sources;

  // TOML array of tables syntax: [[sources]]
  if (auto sources_array = config_["sources"].as_array()) {
    for (auto &&elem : *sources_array) {
      if (auto table = elem.as_table()) {
        SourceConfig src;
        src.id = table->get("id")->value_or(0);
        src.name = table->get("name")->value_or("");
        src.path = table->get("path")->value_or("");
        sources.push_back(src);
      }
    }
  }

  return sources;
}

PSDSourceConfig ConfigParser::GetPSDSources() const {
  PSDSourceConfig psd;
  psd.valid = false;

  if (!IsPSDCapable())
    return psd;

  auto psd_table = config_["psd"].as_table();
  if (!psd_table)
    return psd;

  psd.first_source_id = GetOr<int>("psd.first_source", -1);
  psd.second_source_id = GetOr<int>("psd.second_source", -1);
  psd.first_particle_type = GetOr<std::string>("psd.first_particle_type", "");
  psd.second_particle_type = GetOr<std::string>("psd.second_particle_type", "");

  // Get source names from source list
  auto sources = GetSources();
  for (const auto &src : sources) {
    if (src.id == psd.first_source_id)
      psd.first_source_name = src.name;
    if (src.id == psd.second_source_id)
      psd.second_source_name = src.name;
  }

  psd.valid =
      (psd.first_source_id >= 0 && psd.second_source_id >= 0 &&
       !psd.first_source_name.empty() && !psd.second_source_name.empty());

  return psd;
}

void ConfigParser::DetermineDetectorCapabilities() {
  std::string type_str = GetOr<std::string>("detector.type", "unknown");

  if (type_str == "scintillator") {
    detector_type_ = DetectorType::SCINTILLATOR;
    psd_capable_ = true;
  } else if (type_str == "semiconductor") {
    detector_type_ = DetectorType::SEMICONDUCTOR;
    psd_capable_ = false;
  } else {
    std::cerr << "The only two supported classes of detectors are "
                 "scintillators and semiconductors."
              << std::endl;
    psd_capable_ = false;
  }
}

bool ConfigParser::Validate() {
  bool valid = true;

  // Required fields
  std::vector<std::string> required = {"detector.type", "processing.polarity",
                                       "processing.trigger_threshold",
                                       "processing.output_file"};

  for (const auto &key : required) {
    if (!Get<std::string>(key) && !Get<double>(key) && !Get<bool>(key)) {
      std::cerr << "Missing required key: " << key << std::endl;
      valid = false;
    }
  }

  // PSD-specific validation
  if (IsPSDCapable()) {
    if (!Get<int>("psd.short_gate") || !Get<int>("psd.long_gate")) {
      std::cerr << "PSD-capable detector missing gate configuration"
                << std::endl;
      valid = false;
    }
  }

  return valid;
}

WaveformProcessingUtils *ConfigParser::CreateWaveformProcessor() const {
  auto processor = new WaveformProcessingUtils();

  processor->SetPolarity(GetOr<std::string>("processing.polarity", "positive"));
  processor->SetTriggerThreshold(
      GetOr<double>("processing.trigger_threshold", 0.1));
  processor->SetSampleWindows(GetOr<int>("processing.pre_samples", 20),
                              GetOr<int>("processing.post_samples", 200));
  processor->SetMaxEvents(GetOr<int>("processing.max_events", -1));
  processor->SetVerbose(GetOr<bool>("processing.verbose", false));

  if (IsPSDCapable()) {
    processor->SetPSDGates(GetOr<int>("psd.short_gate", 30),
                           GetOr<int>("psd.long_gate", 100));
  }

  // Add sources
  for (const auto &src : GetSources()) {
    processor->AddSource(src.name, src.id);
  }

  return processor;
}

void ConfigParser::Print() const {
  std::cout << "\n=== Configuration ===\n";
  std::cout << config_ << std::endl;
}
