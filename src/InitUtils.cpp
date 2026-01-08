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

Bool_t InitUtils::ConvertCoMPASSBinToROOT(const TString input_filename,
                                          const TString output_name) {
  if (gSystem->AccessPathName("root_files")) {
    gSystem->mkdir("root_files", kTRUE);
  }

  if (gSystem->AccessPathName(input_filename)) {
    std::cout << "Error: Input file does not exist: " << input_filename
              << std::endl;
    return kFALSE;
  }

  TString output_filename = "root_files/" + output_name + ".root";

  CoMPASSReader reader;
  if (!reader.Open(input_filename.Data())) {
    std::cout << "Error: Failed to open CoMPASS binary file" << std::endl;
    return kFALSE;
  }

  TFile *outfile = new TFile(output_filename, "RECREATE");
  if (!outfile || outfile->IsZombie()) {
    std::cout << "Error: Could not create output file " << output_filename
              << std::endl;
    reader.Close();
    return kFALSE;
  }

  TTree *tree = new TTree("compass_data", "CoMPASS Binary Data");

  UShort_t header, board, channel, energy_ch, energy_short_ch;
  UShort_t global_header;
  ULong64_t timestamp;
  Double_t energy_cal;
  UInt_t flags, num_samples;
  UChar_t waveform_code;
  std::vector<UShort_t> *samples = new std::vector<UShort_t>();

  global_header = reader.GetGlobalHeader();

  tree->Branch("global_header", &global_header, "global_header/s");
  tree->Branch("header", &header, "header/s");
  tree->Branch("board", &board, "board/s");
  tree->Branch("channel", &channel, "channel/s");
  tree->Branch("timestamp", &timestamp, "timestamp/l");
  tree->Branch("energy_ch", &energy_ch, "energy_ch/s");
  tree->Branch("energy_cal", &energy_cal, "energy_cal/D");
  tree->Branch("energy_short_ch", &energy_short_ch, "energy_short_ch/s");
  tree->Branch("flags", &flags, "flags/i");
  tree->Branch("waveform_code", &waveform_code, "waveform_code/b");
  tree->Branch("num_samples", &num_samples, "num_samples/i");
  tree->Branch("samples", &samples);

  Long64_t event_count = 0;
  Long64_t print_interval = 10000;

  std::cout << "Reading events..." << std::endl;

  while (reader.ReadEvent()) {
    const CoMPASSData &event = reader.GetCurrentEvent();

    header = event.header;
    board = event.board;
    channel = event.channel;
    timestamp = event.timestamp;
    energy_ch = event.energy_ch;
    energy_cal = event.energy_cal;
    energy_short_ch = event.energy_short_ch;
    flags = event.flags;
    waveform_code = event.waveform_code;
    num_samples = event.num_samples;
    *samples = event.samples;

    tree->Fill();
    event_count++;

    if (event_count % print_interval == 0) {
      std::cout << "Processed " << event_count << " events..." << std::endl;
    }
  }

  std::cout << "Conversion complete." << std::endl;
  std::cout << "Total events processed: " << event_count << std::endl;
  std::cout << "Total bytes read: " << reader.GetBytesRead() << std::endl;

  outfile->cd();
  tree->Write();
  outfile->Close();
  reader.Close();

  delete outfile;
  delete samples;

  std::cout << "Output saved to: " << output_filename << std::endl;

  return kTRUE;
}
