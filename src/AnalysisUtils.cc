#include "AnalysisUtils.hh"
std::map<Int_t, std::vector<Float_t>>
AnalysisUtils::CalculateAverageWaveforms(const std::string &filename,
                                         const SpectralCuts &cuts) {

  TFile *wf_file = TFile::Open(filename.c_str(), "READ");
  if (!wf_file || wf_file->IsZombie()) {
    std::cout << "Error: Could not load processed_waveforms.root!" << std::endl;
    return {};
  }

  TTree *tree = static_cast<TTree *>(wf_file->Get("features"));
  if (!tree) {
    std::cout << "Error: Could not find features tree!" << std::endl;
    wf_file->Close();
    return {};
  }

  // Set up branch addresses
  TArrayS *samples = nullptr;
  Int_t source_id;
  Float_t light_output_keVee;

  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("source_id", &source_id);
  tree->SetBranchAddress("light_output_keVee", &light_output_keVee);

  std::map<Int_t, std::vector<std::vector<Float_t>>> waveforms_by_source;
  std::map<Int_t, Int_t> event_counts;

  Float_t min_light_output = cuts.min_light_output;
  Float_t max_light_output = cuts.max_light_output;
  Int_t max_events_per_source = cuts.max_events_per_source;

  Long64_t n_entries = tree->GetEntries();
  Int_t length = 0;
  tree->GetEntry(0);

  for (Long64_t i = 0; i < n_entries; ++i) {
    if (tree->GetEntry(i) <= 0)
      continue;

    if (event_counts[source_id] >= max_events_per_source)
      continue;

    if (light_output_keVee < min_light_output ||
        light_output_keVee > max_light_output)
      continue;

    if (i == 0) {
      length = samples->GetSize();
    }

    // Convert TArrayS to vector
    std::vector<Float_t> waveform;
    for (Int_t j = 0; j < samples->GetSize(); ++j) {
      waveform.push_back(Float_t(samples->At(j)));
    }

    waveforms_by_source[source_id].push_back(waveform);
    event_counts[source_id]++;
  }

  // Calculate average waveforms
  std::map<Int_t, std::vector<Float_t>> average_waveforms;
  for (const auto &source_pair : waveforms_by_source) {
    Int_t src_id = source_pair.first;
    const std::vector<std::vector<Float_t>> &waveforms = source_pair.second;

    if (waveforms.empty())
      continue;

    std::vector<Float_t> avg_waveform(length, 0.0);
    for (const auto &wf : waveforms) {
      for (Int_t i = 0; i < length; ++i) {
        avg_waveform[i] += wf[i];
      }
    }

    // Normalize by number of waveforms
    for (Int_t i = 0; i < length; ++i) {
      avg_waveform[i] /= Float_t(waveforms.size());
    }

    // Normalize to peak = 1 for better comparison
    Float_t peak = *std::max_element(avg_waveform.begin(), avg_waveform.end());
    if (peak > 0) {
      for (Float_t &val : avg_waveform) {
        val /= peak;
      }
    }
    average_waveforms[src_id] = avg_waveform;
  }

  SaveAverageWaveforms(average_waveforms);

  return average_waveforms;
}

void AnalysisUtils::SaveAverageWaveforms(
    const std::map<Int_t, std::vector<Float_t>> &avgWaveforms,
    const std::string &filename) {
  TFile *file = new TFile(filename.c_str(), "RECREATE");
  TTree *tree = new TTree("averages", "Average waveforms by source");

  Int_t source_id;
  std::vector<Float_t> avg_waveform;

  TBranch *source_id_branch = tree->GetBranch("source_id");
  TBranch *avg_waveform_branch = tree->GetBranch("avg_waveform");

  if (!source_id_branch || !avg_waveform_branch) {
    tree->Branch("source_id", &source_id);
    tree->Branch("avg_waveform", &avg_waveform);

    for (const auto &pair : avgWaveforms) {
      source_id = pair.first;
      avg_waveform = pair.second;
      tree->Fill();
    }
    tree->AutoSave("SaveSelf");
  } else {
    std::cout
        << "Source ID and average waveform branches already exist, skipping..."
        << std::endl;
  }
  file->Close();
}

std::map<Int_t, std::vector<Float_t>>
AnalysisUtils::LoadAverageWaveforms(const std::string &filename) {
  TFile *file = TFile::Open(filename.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not load " << filename << "!" << std::endl;
    return {};
  }

  TTree *tree = static_cast<TTree *>(file->Get("averages"));
  if (!tree) {
    std::cout << "Error: Could not find averages tree!" << std::endl;
    file->Close();
    return {};
  }

  Int_t source_id;
  std::vector<Float_t> *avg_waveform = nullptr;

  tree->SetBranchAddress("source_id", &source_id);
  tree->SetBranchAddress("avg_waveform", &avg_waveform);

  std::map<Int_t, std::vector<Float_t>> avgWaveforms;

  Long64_t n_entries = tree->GetEntries();
  for (Long64_t i = 0; i < n_entries; ++i) {
    if (tree->GetEntry(i) <= 0)
      continue;
    avgWaveforms[source_id] = *avg_waveform;
  }

  file->Close();
  return avgWaveforms;
}

void AnalysisUtils::SaveSIWeightingFactor(
    const std::vector<Float_t> &SIWeightingFactor,
    const std::string &filename) {

  TFile *file = TFile::Open(filename.c_str(), "UPDATE");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not load average_waveforms.root!" << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("averages"));

  if (!tree) {
    std::cout << "Error: Could not find averages tree!" << std::endl;
    return;
  }

  std::vector<Float_t> si_weighting_factor;

  TBranch *si_weighting_factor_branch = tree->GetBranch("si_weighting_factor");

  if (!si_weighting_factor_branch) {
    tree->Branch("si_weighting_factor", &si_weighting_factor);
    si_weighting_factor = SIWeightingFactor;
    tree->GetBranch("si_weighting_factor")->Fill();

    tree->AutoSave("SaveSelf");
  } else {
    std::cout << "SI weighting factor branch already exists, skipping..."
              << std::endl;
  }
  file->Close();
}

Bool_t AnalysisUtils::CalculateChargeComparisonPSD(
    Int_t short_gate, Int_t long_gate, const std::string &waveforms_file,
    const std::string &histograms_file, Bool_t force_recalculate) {

  // Open input files
  TFile *wffile = TFile::Open(waveforms_file.c_str(), "UPDATE");
  TFile *histfile = TFile::Open(histograms_file.c_str(), "UPDATE");

  if (!wffile || wffile->IsZombie()) {
    std::cout << "Error: Could not open " << waveforms_file << std::endl;
    return kFALSE;
  }
  if (!histfile || histfile->IsZombie()) {
    std::cout << "Error: Could not open " << histograms_file << std::endl;
    wffile->Close();
    return kFALSE;
  }

  // Get the tree
  TTree *wftree = static_cast<TTree *>(wffile->Get("features"));
  if (!wftree) {
    std::cout << "Error: No features tree found" << std::endl;
    wffile->Close();
    histfile->Close();
    return kFALSE;
  }

  // Check if charge_comparison_psd branch already exists
  TBranch *psd_branch = wftree->GetBranch("charge_comparison_psd");
  if (psd_branch && !force_recalculate) {
    std::cout << "Charge comparison PSD branch already exists. Use "
                 "force_recalculate=true to recreate."
              << std::endl;
  } else {
    if (psd_branch && force_recalculate) {
      std::cout << "Recreating charge_comparison_psd branch..." << std::endl;
    } else {
      std::cout << "Creating charge_comparison_psd branch..." << std::endl;
    }

    // Create the branch
    TArrayS *samples = nullptr;
    Int_t source_id;
    Int_t trigger_position;
    Float_t charge_comparison_psd;

    wftree->SetBranchAddress("Samples", &samples);
    wftree->SetBranchAddress("source_id", &source_id);
    wftree->SetBranchAddress("trigger_position",
                             &trigger_position); // Set branch address

    // Remove existing branch if force recalculating
    if (psd_branch && force_recalculate) {
      // ROOT doesn't easily delete branches, so we'll overwrite
      wftree->SetBranchAddress("charge_comparison_psd", &charge_comparison_psd);
    } else {
      wftree->Branch("charge_comparison_psd", &charge_comparison_psd,
                     "charge_comparison_psd/F");
    }
    wftree->GetEntry(0);
    Int_t n_samples = samples->GetSize();

    Long64_t n_entries = wftree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
      wftree->GetEntry(i);
      Float_t CC = 0.0;
      Float_t short_integral = 0.0;
      Float_t long_integral = 0.0;

      Int_t short_end = std::min(trigger_position + short_gate, n_samples);
      Int_t long_end = std::min(trigger_position + long_gate, n_samples);

      for (Int_t j = trigger_position; j < long_end; ++j) {
        Float_t sample_value = Float_t(samples->At(j));
        long_integral += sample_value;
        if (j < short_end) {
          short_integral += sample_value;
        }
      }
      CC = (long_integral - short_integral) / long_integral;
      if (0.0 < CC && CC < 1.0) {
        charge_comparison_psd = CC;
      } else {
        CC = -1;
        charge_comparison_psd = CC;
      }
      wftree->GetBranch("charge_comparison_psd")->Fill();
    }
    wftree->AutoSave("SaveSelf");
    std::cout << "Charge comparison PSD branch created and saved." << std::endl;
  }

  // Check if histograms already exist and are filled
  bool histograms_already_filled = false;
  TDirectory *psd_dir = histfile->GetDirectory("charge_comparison_spectra");
  if (psd_dir && !force_recalculate) {
    psd_dir->cd();

    // Simply check if directory has any histogram keys
    TList *keys = gDirectory->GetListOfKeys();
    if (keys && keys->GetEntries() > 0) {
      histograms_already_filled = true;
      std::cout << "Charge comparison PSD histograms already filled. Use "
                   "force_recalculate=true to recreate."
                << std::endl;
    }
  }
  if (!histograms_already_filled || force_recalculate) {
    if (force_recalculate) {
      std::cout << "Recreating charge comparison PSD histograms..."
                << std::endl;
    } else {
      std::cout << "Filling charge comparison PSD histograms..." << std::endl;
    }

    Float_t charge_comparison_psd;
    Int_t source_id;

    wftree->SetBranchAddress("charge_comparison_psd", &charge_comparison_psd);
    wftree->SetBranchAddress("source_id", &source_id);

    // Create source name mapping (like HistogramUtils does)
    std::map<Int_t, std::string> source_names = {
        {0, "Am-241"}, {1, "Cs-137"}, {2, "Na-22"}, {3, "Am-241 & Cs-137"}};

    // Create histograms for each source with proper names
    std::map<Int_t, TH1F *> psd_histograms;
    std::set<Int_t> source_ids;

    // First pass: collect all source IDs
    Long64_t n_entries = wftree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
      wftree->GetEntry(i);
      source_ids.insert(source_id);
    }

    // Create histograms with proper names (matching HistogramUtils pattern)
    for (Int_t sid : source_ids) {
      // Get source name and clean it
      auto name_it = source_names.find(sid);
      std::string source_name = (name_it != source_names.end())
                                    ? name_it->second
                                    : "Source_" + std::to_string(sid);

      std::string clean_name = source_name;
      std::replace(clean_name.begin(), clean_name.end(), ' ', '_');
      std::replace(clean_name.begin(), clean_name.end(), '-', '_');
      std::replace(clean_name.begin(), clean_name.end(), '&', '_');

      // Create histogram with proper name (matching HistogramUtils pattern)
      std::string hist_name = "h_charge_comparison_" + clean_name;
      std::string hist_title =
          source_name + " Charge Comparison PSD;PSD Parameter;Counts";

      TH1F *hist =
          new TH1F(hist_name.c_str(), hist_title.c_str(), 100, 0.0, 1.0);
      psd_histograms[sid] = hist;
    }

    // Second pass: fill histograms
    for (Long64_t i = 0; i < n_entries; ++i) {
      wftree->GetEntry(i);
      if (psd_histograms.count(source_id)) {
        if (0 < charge_comparison_psd && charge_comparison_psd < 1) {
          psd_histograms[source_id]->Fill(charge_comparison_psd);
        }
      }
    }

    // Save histograms
    if (!psd_dir) {
      psd_dir = histfile->mkdir("charge_comparison_spectra");
    }
    psd_dir->cd();

    for (auto &pair : psd_histograms) {
      pair.second->Write("", TObject::kOverwrite);
    }

    std::cout << "Charge comparison PSD histograms filled and saved."
              << std::endl;
  }

  wffile->Close();
  histfile->Close();
  return kTRUE;
}

TH1F *AnalysisUtils::GetFilteredPSDHistogram(
    Float_t min_light_output, Float_t max_light_output, Int_t source_id_calc,
    const std::string &waveforms_file) {

  TFile *file = TFile::Open(waveforms_file.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << waveforms_file << std::endl;
    return nullptr;
  }

  TTree *tree = static_cast<TTree *>(file->Get("features"));
  if (!tree) {
    std::cout << "Error: Could not find features tree" << std::endl;
    file->Close();
    return nullptr;
  }

  // Set up branch addresses
  Float_t psd_value;
  Float_t light_output_keVee;
  Int_t source_id;

  tree->SetBranchAddress("charge_comparison_psd", &psd_value);
  tree->SetBranchAddress("light_output_keVee", &light_output_keVee);
  tree->SetBranchAddress("source_id", &source_id);
  std::string hist_name =
      "hist_charge_comparison_psd_" + std::to_string(source_id_calc);
  TH1F *hist = new TH1F(hist_name.c_str(), ";PSD;Counts", 200, 0, 1);
  hist->SetDirectory(0);
  Long64_t n_entries = tree->GetEntries();

  for (Long64_t i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);

    if (light_output_keVee < min_light_output ||
        light_output_keVee > max_light_output) {
      continue;
    }

    if (source_id == source_id_calc) {
      hist->Fill(psd_value);
    }
  }

  file->Close();
  return hist;
}

TH1F *AnalysisUtils::CalculateChargeComparisonPSDOnTheFly(
    Int_t short_gate, Int_t long_gate, Float_t min_light_output,
    Float_t max_light_output, Int_t source_id_calc,
    const std::string &waveforms_file) {

  TFile *file = TFile::Open(waveforms_file.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << waveforms_file << std::endl;
    return nullptr;
  }

  TTree *tree = static_cast<TTree *>(file->Get("features"));
  if (!tree) {
    std::cout << "Error: Could not find features tree" << std::endl;
    file->Close();
    return nullptr;
  }

  // Set up branch addresses
  TArrayS *samples = nullptr;
  Float_t light_output_keVee;
  Int_t source_id;
  Int_t trigger_position;

  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("light_output_keVee", &light_output_keVee);
  tree->SetBranchAddress("source_id", &source_id);
  tree->SetBranchAddress("trigger_position", &trigger_position);

  // Create histogram
  std::string hist_name = "hist_cc_" + std::to_string(short_gate) + "_" +
                          std::to_string(long_gate) + "_" +
                          std::to_string(source_id_calc);
  TH1F *hist = new TH1F(hist_name.c_str(), ";CC PSD;Counts", 100, 0, 1);
  hist->SetDirectory(0);

  // Determine entry range based on source ID
  Long64_t start_entry = 0;
  Long64_t end_entry = 0;

  if (source_id_calc == 0) {
    start_entry = 0;
    end_entry = 500000;
  } else if (source_id_calc == 2) {
    start_entry = 1000000;
    end_entry = 1500000;
  } else {
    std::cout << "Error: Unsupported source_id: " << source_id_calc
              << std::endl;
    file->Close();
    delete hist;
    return nullptr;
  }

  // Process entries
  for (Long64_t i = start_entry; i < end_entry; ++i) {
    tree->GetEntry(i);

    // Apply light output filter
    if (light_output_keVee < min_light_output ||
        light_output_keVee > max_light_output) {
      continue;
    }

    // Calculate PSD on the fly
    Int_t n_samples = samples->GetSize();
    Float_t short_integral = 0.0;
    Float_t long_integral = 0.0;

    Int_t short_end = std::min(trigger_position + short_gate, n_samples);
    Int_t long_end = std::min(trigger_position + long_gate, n_samples);

    for (Int_t j = trigger_position; j < long_end; ++j) {
      Float_t sample_value = Float_t(samples->At(j));
      long_integral += sample_value;
      if (j < short_end) {
        short_integral += sample_value;
      }
    }

    if (long_integral > 0) {
      Float_t cc_psd = (long_integral - short_integral) / long_integral;
      if (cc_psd > 0.0 && cc_psd < 1.0) {
        hist->Fill(cc_psd);
      }
    }
  }

  file->Close();
  return hist;
}

std::vector<AnalysisUtils::GateOptimizationResult>
AnalysisUtils::OptimizeChargeComparisonGates(
    Float_t min_light_output, Float_t max_light_output, Int_t alpha_source_id,
    Int_t gamma_source_id, Int_t min_short_gate, Int_t max_short_gate,
    Int_t short_gate_step, Int_t min_long_gate, Int_t max_long_gate,
    Int_t long_gate_step, const std::string &waveforms_file) {

  std::cout << "Starting gate optimization..." << std::endl;
  std::cout << "  Light output range: " << min_light_output << "-"
            << max_light_output << " keVee" << std::endl;
  std::cout << "  Short gate range: " << min_short_gate << "-" << max_short_gate
            << " (step " << short_gate_step << ")" << std::endl;
  std::cout << "  Long gate range: " << min_long_gate << "-" << max_long_gate
            << " (step " << long_gate_step << ")" << std::endl;

  std::vector<GateOptimizationResult> results;

  Int_t total_short_gates =
      (max_short_gate - min_short_gate) / short_gate_step + 1;
  Int_t total_long_gates = (max_long_gate - min_long_gate) / long_gate_step + 1;
  Int_t total_combinations = total_short_gates * total_long_gates;
  Int_t current_combination = 0;

  std::cout << "  Total combinations to test: " << total_combinations
            << std::endl;

  for (Int_t short_gate = min_short_gate; short_gate <= max_short_gate;
       short_gate += short_gate_step) {
    for (Int_t long_gate = min_long_gate; long_gate <= max_long_gate;
         long_gate += long_gate_step) {
      current_combination++;

      // Skip invalid combinations (short gate must be < long gate)
      if (short_gate >= long_gate) {
        std::cout << "Skipping invalid combination: short=" << short_gate
                  << ", long=" << long_gate << " (short >= long)" << std::endl;
        continue;
      }

      std::cout << "Testing gates " << short_gate << "/" << long_gate << " ("
                << current_combination << "/" << total_combinations << ")"
                << std::endl;

      // Calculate PSD histograms for both sources
      TH1F *hist_alpha = CalculateChargeComparisonPSDOnTheFly(
          short_gate, long_gate, min_light_output, max_light_output,
          alpha_source_id, waveforms_file);

      TH1F *hist_gamma = CalculateChargeComparisonPSDOnTheFly(
          short_gate, long_gate, min_light_output, max_light_output,
          gamma_source_id, waveforms_file);

      if (!hist_alpha || !hist_gamma) {
        std::cout << "  Failed to create histograms" << std::endl;
        if (hist_alpha)
          delete hist_alpha;
        if (hist_gamma)
          delete hist_gamma;
        continue;
      }

      Int_t alpha_events = hist_alpha->GetEntries();
      Int_t gamma_events = hist_gamma->GetEntries();

      if (alpha_events < 100 || gamma_events < 100) {
        std::cout << "  Insufficient statistics (α:" << alpha_events 
                  << ", γ:" << gamma_events << ")" << std::endl;
        delete hist_alpha;
        delete hist_gamma;
        continue;
      }

      // Calculate FOM
      Double_t separation, fwhm_alpha, fwhm_gamma;
      Double_t fom = PlottingUtils::CalculateFigureOfMerit(
          hist_alpha, hist_gamma, separation, fwhm_alpha, fwhm_gamma);

      if (fom > 0) {
        GateOptimizationResult result;
        result.short_gate = short_gate;
        result.long_gate = long_gate;
        result.fom = fom;
        result.alpha_events = alpha_events;
        result.gamma_events = gamma_events;

        results.push_back(result);

        std::cout << "  ✓ FOM: " << fom << " (α:" << alpha_events 
                  << ", γ:" << gamma_events << ")" << std::endl;
      } else {
        std::cout << "  ✗ FOM calculation failed" << std::endl;
      }

      delete hist_alpha;
      delete hist_gamma;
    }
  }

  // Sort results by FOM (highest first)
  std::sort(results.begin(), results.end(),
            [](const GateOptimizationResult &a,
               const GateOptimizationResult &b) { return a.fom > b.fom; });

  std::cout << "Optimization completed. Found " << results.size()
            << " valid combinations." << std::endl;

  return results;
}
