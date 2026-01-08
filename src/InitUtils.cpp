#include "InitUtils.hpp"

void InitUtils::SetROOTPreferences(Bool_t setupPlotting) {
  if (setupPlotting) {
    PlottingUtils::SetStylePreferences();
  }
  gROOT->ForceStyle(kTRUE);
  gROOT->SetBatch(kTRUE);
  if (gSystem->AccessPathName("plots")) {
    gSystem->mkdir("plots", kTRUE);
  }
  if (gSystem->AccessPathName("root_files")) {
    gSystem->mkdir("root_files", kTRUE);
  }
}

UShort_t InitUtils::ConvertCoMPASSBinToROOT(const TString input_filename,
                                            const TString output_name,
                                            UShort_t global_header_override,
                                            Bool_t skip_bad_events) {
  if (gSystem->AccessPathName("root_files")) {
    gSystem->mkdir("root_files", kTRUE);
  }

  if (gSystem->AccessPathName(input_filename)) {
    std::cout << "Error: Input file does not exist: " << input_filename
              << std::endl;
    return 0;
  }

  TString output_filename = "root_files/" + output_name + ".root";

  CoMPASSReader reader;
  Bool_t open_success =
      (global_header_override != 0)
          ? reader.Open(input_filename.Data(), global_header_override)
          : reader.Open(input_filename.Data());

  if (!open_success) {
    std::cout << "Error: Failed to open CoMPASS binary file" << std::endl;
    return 0;
  }

  UShort_t global_header = reader.GetGlobalHeader();

  Bool_t has_energy_ch = (global_header & 0x0001);
  Bool_t has_energy_cal = (global_header & 0x0002);
  Bool_t has_energy_short = (global_header & 0x0004);
  Bool_t has_waveform = (global_header & 0x0008);

  TFile *outfile = new TFile(output_filename, "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Error: Could not create output file " << output_filename
              << std::endl;
    reader.Close();
    return 0;
  }

  TTree *tree = new TTree("Data_R", "CoMPASS Binary Data");

  UShort_t board, channel, energy, energy_short;
  ULong64_t timestamp;
  Double_t energy_cal;
  UInt_t flags, num_samples;
  UChar_t waveform_code;
  std::vector<UShort_t> *samples = nullptr;

  tree->Branch("Board", &board, "Board/s");
  tree->Branch("Channel", &channel, "Channel/s");
  tree->Branch("Timestamp", &timestamp, "Timestamp/l");

  if (has_energy_ch) {
    tree->Branch("Energy", &energy, "Energy/s");
    std::cout << "Energy type: Channel (ADC counts)" << std::endl;
  } else if (has_energy_cal) {
    tree->Branch("Energy", &energy_cal, "Energy/D");
    std::cout << "Energy type: Calibrated (keV/MeV)" << std::endl;
  }

  if (has_energy_short) {
    tree->Branch("EnergyShort", &energy_short, "EnergyShort/s");
  }

  tree->Branch("Flags", &flags, "Flags/i");

  if (has_waveform) {
    samples = new std::vector<UShort_t>();
    tree->Branch("WaveformCode", &waveform_code, "WaveformCode/b");
    tree->Branch("NumSamples", &num_samples, "NumSamples/i");
    tree->Branch("Samples", &samples);
  }

  Long64_t event_count = 0;
  Long64_t warning_fake = 0;
  Long64_t warning_saturated = 0;
  Long64_t warning_pileup = 0;
  Long64_t warning_memory_full = 0;
  Long64_t warning_trigger_lost = 0;
  Long64_t warning_pll_loss = 0;
  Long64_t warning_over_temp = 0;
  Long64_t warning_adc_shutdown = 0;

  std::cout << "Reading events..." << std::endl;
  if (skip_bad_events) {
    std::cout
        << "Filtering enabled: skipping fake, saturated, and pileup events"
        << std::endl;
  }

  while (reader.ReadEvent()) {
    const CoMPASSData &event = reader.GetCurrentEvent();
    if (event_count == 0) {
      event.PrintHeader();
    }

    if (event.isFakeEvent()) {
      warning_fake++;
      if (skip_bad_events)
        continue;
    }
    if (event.isInputSaturating() || event.hasSaturation()) {
      warning_saturated++;
      if (skip_bad_events)
        continue;
    }
    if (event.isPileup()) {
      warning_pileup++;
      if (skip_bad_events)
        continue;
    }

    if (event.hasMemoryFull())
      warning_memory_full++;
    if (event.hasTriggerLost())
      warning_trigger_lost++;
    if (event.hasPLLLockLoss())
      warning_pll_loss++;
    if (event.isOverTemperature())
      warning_over_temp++;
    if (event.isADCShutdown())
      warning_adc_shutdown++;

    board = event.board;
    channel = event.channel;
    timestamp = event.timestamp;
    flags = event.flags;

    if (has_energy_ch) {
      energy = event.energy_ch;
    }
    if (has_energy_cal) {
      energy_cal = event.energy_cal;
    }
    if (has_energy_short) {
      energy_short = event.energy_short_ch;
    }
    if (has_waveform) {
      waveform_code = event.waveform_code;
      num_samples = event.num_samples;
      *samples = event.samples;
    }

    tree->Fill();
    event_count++;
  }

  std::cout << "Conversion complete." << std::endl;
  std::cout << "Total events processed: " << event_count << std::endl;

  if (warning_fake > 0 || warning_saturated > 0 || warning_pileup > 0) {
    std::cout << "Events with rejection-quality flags:" << std::endl;
    if (warning_fake > 0) {
      std::cout << "  Fake events: " << warning_fake;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    if (warning_saturated > 0) {
      std::cout << "  Saturated: " << warning_saturated;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    if (warning_pileup > 0) {
      std::cout << "  Pileup: " << warning_pileup;
      if (skip_bad_events)
        std::cout << " (rejected)";
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  if (warning_memory_full > 0) {
    std::cout << "Warning: " << warning_memory_full
              << " events with memory full flag" << std::endl;
  }
  if (warning_trigger_lost > 0) {
    std::cout << "Warning: " << warning_trigger_lost
              << " events with trigger lost flag" << std::endl;
  }
  if (warning_pll_loss > 0) {
    std::cout << "Warning: " << warning_pll_loss << " events with PLL lock loss"
              << std::endl;
  }
  if (warning_over_temp > 0) {
    std::cout << "Warning: " << warning_over_temp
              << " events with over temperature" << std::endl;
  }
  if (warning_adc_shutdown > 0) {
    std::cout << "Warning: " << warning_adc_shutdown
              << " events with ADC shutdown" << std::endl;
  }

  std::cout << "Total bytes read: " << reader.GetBytesRead() << std::endl;

  outfile->cd();
  tree->Write("", TObject::kOverwrite);
  outfile->Close();
  reader.Close();

  delete outfile;
  if (samples) {
    delete samples;
  }

  std::cout << "Output saved to: " << output_filename << std::endl;

  return global_header;
}
