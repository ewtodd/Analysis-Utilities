#include "ConfigParser.hh"
#include "AnalysisUtils.hh"
#include "CalibrationUtils.hh"
#include "HistogramUtils.hh"
#include "WaveformProcessingUtils.hh"

bool ConfigParser::Load(const std::string &filename) {
  std::ifstream file(filename);
  if (!file.is_open()) {
    std::cerr << "ERROR: Cannot open config file: " << filename << std::endl;
    return false;
  }

  std::string line;
  int line_num = 0;

  while (std::getline(file, line)) {
    line_num++;
    line = trim(line);

    if (line.empty() || line[0] == '#')
      continue;

    size_t pos = line.find('=');
    if (pos != std::string::npos) {
      std::string key = trim(line.substr(0, pos));
      std::string value = trim(line.substr(pos + 1));

      if (key.empty()) {
        std::cerr << "WARNING: Empty key on line " << line_num << std::endl;
        continue;
      }

      values_[key] = value;
    } else {
      std::cerr << "WARNING: Invalid line " << line_num << ": " << line
                << std::endl;
    }
  }

  file.close();

  DetermineDetectorCapabilities();

  return true;
}

bool ConfigParser::Validate() {
  bool valid = true;

  std::vector<std::string> required_common = {"detector.type",
                                              "processing.polarity",
                                              "processing.trigger_threshold",
                                              "processing.output_file",
                                              "psd.short_gate",
                                              "psd.long_gate"};

  for (const auto &key : required_common) {
    if (values_.find(key) == values_.end()) {
      std::cerr << "ERROR: Missing required key: " << key << std::endl;
      valid = false;
    }
  }

  if (detector_type_ == DetectorType::UNKNOWN) {
    std::cerr << "ERROR: Invalid or missing detector.type (must be "
                 "'scintillator' or 'semiconductor')"
              << std::endl;
    valid = false;
  }

  if (IsScintillator()) {
    if (values_.find("histogram.integral_min") == values_.end()) {
      std::cerr << "ERROR: Scintillator requires histogram.integral_min"
                << std::endl;
      valid = false;
    }
    if (values_.find("histogram.integral_max") == values_.end()) {
      std::cerr << "ERROR: Scintillator requires histogram.integral_max"
                << std::endl;
      valid = false;
    }

    if (IsPSDCapable()) {
      if (values_.find("psd.optimize_gates") == values_.end()) {
        std::cerr
            << "WARNING: PSD-capable detector but no gate optimization defined"
            << std::endl;
      }
    }
  }

  if (IsSemiconductor()) {
    if (values_.find("histogram.pulse_height_min") == values_.end()) {
      std::cerr << "ERROR: Semiconductor requires histogram.pulse_height_min"
                << std::endl;
      valid = false;
    }
    if (values_.find("histogram.pulse_height_max") == values_.end()) {
      std::cerr << "ERROR: Semiconductor requires histogram.pulse_height_max"
                << std::endl;
      valid = false;
    }
  }

  auto sources = GetSources();
  if (sources.empty()) {
    std::cerr << "ERROR: No sources defined (need at least source.0.name and "
                 "source.0.path)"
              << std::endl;
    valid = false;
  }

  return valid;
}

void ConfigParser::Print() const {

  std::cout << "Detector type: ";
  switch (detector_type_) {
  case DetectorType::SCINTILLATOR:
    std::cout << "scintillator";
    if (psd_capable_)
      std::cout << " (PSD capable)";
    break;
  case DetectorType::SEMICONDUCTOR:
    std::cout << "semiconductor";
    break;
  case DetectorType::UNKNOWN:
    std::cout << "UNKNOWN - Please specify!";
    break;
  }
  std::cout << std::endl;

  std::string det_name = GetString("detector.name");
  if (!det_name.empty()) {
    std::cout << "Name: " << det_name << std::endl;
  }

  std::string det_model = GetString("detector.model");
  if (!det_model.empty()) {
    std::cout << "Model: " << det_model << std::endl;
  }

  auto sources = GetSources();
  if (!sources.empty()) {
    std::cout << "Sources:" << std::endl;

    for (const auto &src : sources) {
      std::cout << "[" << src.id << "] " << src.name << std::endl;
      std::cout << src.path << std::endl;
    }
  }

  std::cout << "Enabled features: " << std::endl;

  if (IsScintillator()) {
    std::cout << "Integral-based light output" << std::endl;
    if (IsPSDCapable()) {
      std::cout << "Pulse Shape Discrimination" << std::endl;
    }
  }

  if (IsSemiconductor()) {
    std::cout << "Pulse height energy measurement" << std::endl;
  }

  if (GetBool("calibration.enabled", false)) {
    std::cout << "Calibrated energy" << std::endl;
  }

  std::cout << "Processing parameters: " << std::endl;

  std::cout << "Polarity: " << GetString("processing.polarity") << std::endl;
  std::cout << "Trigger threshold: " << GetFloat("processing.trigger_threshold")
            << std::endl;
  std::cout << "Sample window: " << GetInt("processing.pre_samples")
            << " (pre) + " << GetInt("processing.post_samples") << " (post)"
            << std::endl;
  std::cout << "Max events: " << GetInt("processing.max_events", -1)
            << " (-1 = all)" << std::endl;

  std::cout << "Output file names: " << std::endl;

  std::cout << "Waveforms: " << GetString("processing.output_file")
            << std::endl;
  std::cout << "Histograms: " << GetString("histogram.output_file")
            << std::endl;
}

std::string ConfigParser::GetString(const std::string &key,
                                    const std::string &default_val) const {
  auto it = values_.find(key);
  return (it != values_.end()) ? it->second : default_val;
}

Float_t ConfigParser::GetFloat(const std::string &key,
                               Float_t default_val) const {
  std::string val = GetString(key);
  if (val.empty())
    return default_val;
  try {
    return std::stof(val);
  } catch (...) {
    std::cerr << "ERROR: Cannot convert '" << val << "' to float for key '"
              << key << "'" << std::endl;
    return default_val;
  }
}

Int_t ConfigParser::GetInt(const std::string &key, Int_t default_val) const {
  std::string val = GetString(key);
  if (val.empty())
    return default_val;
  try {
    return std::stoi(val);
  } catch (...) {
    std::cerr << "ERROR: Cannot convert '" << val << "' to int for key '" << key
              << "'" << std::endl;
    return default_val;
  }
}

Bool_t ConfigParser::GetBool(const std::string &key, Bool_t default_val) const {
  std::string val = GetString(key);
  if (val.empty())
    return default_val;

  std::string val_lower = val;
  std::transform(val_lower.begin(), val_lower.end(), val_lower.begin(),
                 ::tolower);
  return (val_lower == "true" || val_lower == "1" || val_lower == "yes");
}

std::vector<SourceConfig> ConfigParser::GetSources() const {
  std::vector<SourceConfig> sources;

  for (int i = 0; i < 100; ++i) {
    std::string name_key = "source." + std::to_string(i) + ".name";
    std::string path_key = "source." + std::to_string(i) + ".path";

    std::string name = GetString(name_key);
    std::string path = GetString(path_key);

    if (!name.empty() && !path.empty()) {
      SourceConfig src;
      src.id = i;
      src.name = name;
      src.path = path;
      sources.push_back(src);
    } else if (!name.empty() || !path.empty()) {
      std::cerr << "WARNING: Incomplete source definition for source." << i
                << std::endl;
    }
  }

  return sources;
}

Int_t ConfigParser::GetSourceID(const std::string &source_name) const {
  auto sources = GetSources();
  for (const auto &src : sources) {
    if (src.name == source_name) {
      return src.id;
    }
  }
  std::cerr << "WARNING: Source '" << source_name << "' not found" << std::endl;
  return -1;
}

PSDSourceConfig ConfigParser::GetPSDSourceConfig() const {
  PSDSourceConfig config;

  std::string first_name = GetString("charge_comparison.first_source", "");
  std::string second_name = GetString("charge_comparison.second_source", "");

  config.first_source_name = first_name;
  config.second_source_name = second_name;

  // Get particle types
  config.first_particle_type =
      GetString("charge_comparison.first_particle_type", "alpha");
  config.second_particle_type =
      GetString("charge_comparison.second_particle_type", "gamma");

  // Get source IDs
  config.first_source_id = GetSourceID(first_name);
  config.second_source_id = GetSourceID(second_name);

  // Validate
  config.valid = (config.first_source_id >= 0 && config.second_source_id >= 0);

  if (!config.valid) {
    std::cerr << "ERROR: Could not find charge comparison source IDs"
              << std::endl;
    std::cerr << "  Looking for first source (" << config.first_particle_type
              << "): " << first_name << std::endl;
    std::cerr << "  Looking for second source (" << config.second_particle_type
              << "): " << second_name << std::endl;
  }

  return config;
}

std::string ConfigParser::GetProcessedWaveformsFile() const {
  return GetString("processing.output_file", "processed_waveforms.root");
}

std::string ConfigParser::GetHistogramsFile() const {
  return GetString("histogram.output_file", "histograms.root");
}

std::string ConfigParser::GetAverageWaveformsFile() const {
  return GetString("analysis.average_waveforms_file", "average_waveforms.root");
}

WaveformProcessingUtils *ConfigParser::CreateWaveformProcessor() const {
  WaveformProcessingUtils *processor = new WaveformProcessingUtils();

  processor->SetPolarity(GetString("processing.polarity"));
  processor->SetTriggerThreshold(GetFloat("processing.trigger_threshold"));
  processor->SetSampleWindows(GetInt("processing.pre_samples"),
                              GetInt("processing.post_samples"));
  processor->SetMaxEvents(GetInt("processing.max_events", -1));
  processor->SetPSDGates(GetInt("psd.short_gate"), GetInt("psd.long_gate"));
  processor->SetVerbose(GetBool("processing.verbose", kTRUE));
  processor->SetStoreWaveforms(GetBool("processing.store_waveforms", kTRUE));

  // Add all sources
  auto sources = GetSources();
  for (const auto &src : sources) {
    processor->AddSource(src.name, src.id);
  }

  std::cout << "✓ WaveformProcessor configured with " << sources.size() 
            << " sources" << std::endl;

  return processor;
}

HistogramUtils *ConfigParser::CreateHistogramManager() const {
  HistogramUtils *histMgr = new HistogramUtils();

  HistogramConfig histConfig;

  if (IsScintillator()) {
    histConfig.integral_min = GetFloat("histogram.integral_min", 0);
    histConfig.integral_max = GetFloat("histogram.integral_max", 120000);
    histConfig.integral_bins = GetInt("histogram.integral_bins", 300);
    histConfig.create_ph_histograms = false;
    histConfig.create_integral_histograms = true;

    std::cout << "Histogram mode: Scintillator (pulse integral)" << std::endl;

  } else if (IsSemiconductor()) {
    histConfig.ph_min = GetFloat("histogram.pulse_height_min", 0);
    histConfig.ph_max = GetFloat("histogram.pulse_height_max", 5000);
    histConfig.ph_bins = GetInt("histogram.pulse_height_bins", 200);
    histConfig.create_integral_histograms = false;
    histConfig.create_ph_histograms = true;

    std::cout << "Histogram mode: Semiconductor (pulse height)" << std::endl;
  }

  histMgr->SetConfig(histConfig);

  // Add all sources
  auto sources = GetSources();
  for (const auto &src : sources) {
    histMgr->AddSource(src.id, src.name);
  }

  std::cout << "HistogramManager configured for " << sources.size()
            << " sources" << std::endl;

  return histMgr;
}

AnalysisUtils *ConfigParser::CreateAnalysisUtils() const {
  AnalysisUtils *analyzer = new AnalysisUtils();
  std::cout << "AnalysisUtils created" << std::endl;
  return analyzer;
}

CalibrationUtils *ConfigParser::CreateCalibrationManager() const {
  if (!GetBool("calibration.enabled", false)) {
    std::cout << "⊘ Calibration disabled in config" << std::endl;
    return nullptr;
  }

  CalibrationUtils *calibMgr = new CalibrationUtils();

  calibMgr->SetIncludeZeroPoint(GetBool("calibration.include_zero", kTRUE));

  // Load calibration peaks from config
  // Format: calibration.source_name.peak_energy = deposited_energy_keV
  // Additional parameters: calibration.source_name.peak_energy.param = value
  auto sources = GetSources();

  int calibration_points = 0;

  for (const auto &src : sources) {
    // Look for calibration points for this source
    // e.g., calibration.Cs-137.662 = 661.7
    std::string prefix = "calibration." + src.name + ".";

    // Find all energy peaks for this source
    std::map<std::string, Double_t> peaks_for_source;

    for (const auto &pair : values_) {
      if (pair.first.find(prefix) == 0) {
        std::string remainder = pair.first.substr(prefix.length());

        size_t dot_pos = remainder.find('.');

        if (dot_pos == std::string::npos) {
          std::string peak_id = remainder;
          Double_t energy_kev = std::stof(pair.second);
          peaks_for_source[peak_id] = energy_kev;
        }
      }
    }

    for (const auto &peak_pair : peaks_for_source) {
      std::string peak_id = peak_pair.first;
      Double_t energy_kev = peak_pair.second;

      std::string peak_prefix = prefix + peak_id + ".";

      Double_t guess_mu =
          GetFloat(peak_prefix + "guess_mu", energy_kev * 100.0);
      Double_t guess_sigma =
          GetFloat(peak_prefix + "guess_sigma", guess_mu * 0.05);
      Double_t guess_amplitude =
          GetFloat(peak_prefix + "guess_amplitude", 1000.0);
      Double_t guess_bkg_const =
          GetFloat(peak_prefix + "guess_bkg_const", 10.0);
      Double_t guess_bkg_slope = GetFloat(peak_prefix + "guess_bkg_slope", 0.0);
      Double_t fit_range_low =
          GetFloat(peak_prefix + "fit_range_low", guess_mu - 3.0 * guess_sigma);
      Double_t fit_range_high = GetFloat(peak_prefix + "fit_range_high",
                                         guess_mu + 3.0 * guess_sigma);

      calibMgr->AddCalibrationPeak(src.name,        // calibration_source
                                   src.id,          // source_id
                                   energy_kev,      // deposited_energy_keV
                                   guess_mu,        // guess_peak_mu
                                   guess_sigma,     // guess_sigma
                                   guess_amplitude, // guess_amplitude
                                   guess_bkg_const, // guess_bkg_const
                                   guess_bkg_slope, // guess_bkg_slope
                                   fit_range_low,   // fit_range_low
                                   fit_range_high   // fit_range_high
      );

      std::cout << "  Added calibration peak: " << src.name << " @ "
                << energy_kev << " keV." << std::endl;

      calibration_points++;
    }
  }

  if (calibration_points > 0) {
    std::cout << "CalibrationManager configured with " << calibration_points
              << " calibration point(s)" << std::endl;
  } else {
    std::cout << "CalibrationManager created but no calibration points defined"
              << std::endl;
    std::cout << "  Add peaks like: calibration.Cs-137.662 = 661.7"
              << std::endl;
  }

  return calibMgr;
}
