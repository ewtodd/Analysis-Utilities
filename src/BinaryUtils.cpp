#include "BinaryUtils.hpp"

CoMPASSData::CoMPASSData()
    : header(0), board(0), channel(0), timestamp(0), energy_ch(0),
      energy_cal(0.0), energy_short_ch(0), flags(0), waveform_code(0),
      num_samples(0) {
  std::cout << "This version of CoMPASS binary conversion is based on manual "
               "revision 25.1..."
            << std::endl;
}

TString CoMPASSData::getWaveformCodeName() const {
  switch (waveform_code) {
  case INPUT:
    return "Input";
  case RC_CR:
    return "RC-CR (DPP-PHA)";
  case RC_CR2:
    return "RC-CR2 (DPP-PHA)";
  case TRAPEZOID:
    return "Trapezoid (DPP-PHA)";
  case BASELINE:
    return "Baseline";
  case THRESHOLD:
    return "Threshold";
  case CFD:
    return "CFD (DPP-PSD)";
  case TRAPEZOID_BASELINE:
    return "Trapezoid-Baseline (DPP-PHA)";
  case FAST_TRIANGLE:
    return "Fast Triangle (x27xx DPP-PHA)";
  case SMOOTHED_INPUT:
    return "Smoothed Input (DPP-PSD)";
  default:
    return "Unknown";
  }
}

std::vector<TString> CoMPASSData::getActiveFlags() const {
  std::vector<TString> active_flags;

  if (hasDeadtime())
    active_flags.push_back("Deadtime");
  if (hasTimestampRollover())
    active_flags.push_back("Timestamp Rollover");
  if (hasTimestampResetExt())
    active_flags.push_back("Timestamp Reset (Ext)");
  if (isFakeEvent())
    active_flags.push_back("Fake Event");
  if (hasMemoryFull())
    active_flags.push_back("Memory Full");
  if (hasTriggerLost())
    active_flags.push_back("Trigger Lost");
  if (hasNTriggersLost())
    active_flags.push_back("N Triggers Lost");
  if (hasSaturation())
    active_flags.push_back("Saturation");
  if (has1024Triggers())
    active_flags.push_back("1024 Triggers");
  if (isFirstAfterBusy())
    active_flags.push_back("First After Busy");
  if (isInputSaturating())
    active_flags.push_back("Input Saturating");
  if (hasNTriggersCounted())
    active_flags.push_back("N Triggers Counted");
  if (isNotMatchedTimeFilter())
    active_flags.push_back("Not Matched Time Filter");
  if (hasFineTimestamp())
    active_flags.push_back("Fine Timestamp");
  if (isPileup())
    active_flags.push_back("Pile-up");
  if (hasPLLLockLoss())
    active_flags.push_back("PLL Lock Loss");
  if (isOverTemperature())
    active_flags.push_back("Over Temperature");
  if (isADCShutdown())
    active_flags.push_back("ADC Shutdown");

  return active_flags;
}

void CoMPASSData::PrintHeader() const {
  std::cout << "CoMPASS event header..." << std::endl;
  std::cout << "Header:    0x" << std::hex << header << std::dec
            << " (Binary: " << std::bitset<16>(header) << ")" << std::endl;

  std::cout << std::endl;
  std::cout << "Control bits..." << std::endl;
  std::cout << "Energy (ch):    " << (hasEnergyCh() ? "YES" : "NO")
            << std::endl;
  std::cout << "Energy (cal):   " << (hasEnergyCal() ? "YES" : "NO")
            << std::endl;
  std::cout << "Energy (short): " << (hasEnergyShort() ? "YES" : "NO")
            << std::endl;
  std::cout << "Waveform:       " << (hasWaveform() ? "YES" : "NO")
            << std::endl;
}

void CoMPASSData::PrintFlags() const {
  std::cout << "Flags (0x" << std::hex << flags << std::dec << ")" << std::endl;
  std::cout << "Binary: " << std::bitset<32>(flags) << std::endl;

  std::vector<TString> active = getActiveFlags();
  if (active.empty()) {
    std::cout << "No flags set" << std::endl;
  } else {
    std::cout << "Active flags:" << std::endl;
    Int_t n_flags = active.size();
    for (Int_t flag = 0; flag < n_flags; flag++) {
      std::cout << "  - " << active.at(flag).Data() << std::endl;
    }
  }
}

void CoMPASSData::PrintWaveform() const {
  if (hasWaveform()) {
    std::cout << "Waveform..." << std::endl;
    std::cout << "Code:    " << static_cast<Int_t>(waveform_code) << " ("
              << getWaveformCodeName().Data() << ")" << std::endl;
    std::cout << "Samples: " << num_samples << std::endl;
    if (samples.GetSize() > 0) {
      std::cout << "First 5 samples: ";
      for (Int_t i = 0; i < TMath::Min(5, samples.GetSize()); i++) {
        std::cout << samples[i] << " ";
      }
      std::cout << std::endl;
    }
  }
}

void CoMPASSData::Print() const {
  PrintHeader();

  if (hasEnergyCh()) {
    std::cout << "Energy (ch):    " << energy_ch << std::endl;
  }
  if (hasEnergyCal()) {
    std::cout << "Energy (cal):   " << energy_cal << " keV/MeV" << std::endl;
  }
  if (hasEnergyShort()) {
    std::cout << "Energy (short): " << energy_short_ch << std::endl;
  }

  PrintFlags();
  PrintWaveform();
}

Bool_t CoMPASSReader::Open(const char *fname) {
  if (!BinaryReader::Open(fname)) {
    return kFALSE;
  }

  file.read(reinterpret_cast<char *>(&global_header), sizeof(UShort_t));
  bytes_read += sizeof(UShort_t);

  if ((global_header & 0xCAE0) != 0xCAE0) {
    std::cerr << "WARNING: Header does not match 0xCAEx pattern: 0x" << std::hex
              << global_header << std::dec << std::endl;
  }

  std::cout << "CoMPASS file opened: " << fname << std::endl;

  return kTRUE;
}

Bool_t CoMPASSReader::Open(const char *fname, UShort_t header_override) {
  if (!BinaryReader::Open(fname)) {
    return kFALSE;
  }

  if (header_override != 0) {
    global_header = header_override;
    std::cout << "CoMPASS continuation file opened: " << fname << std::endl;
    std::cout << "Using provided global header: 0x" << std::hex << global_header
              << std::dec << std::endl;
    std::cout << "Control bits: " << std::bitset<4>(global_header & 0x000F)
              << std::endl;
  } else {
    file.read(reinterpret_cast<char *>(&global_header), sizeof(UShort_t));
    bytes_read += sizeof(UShort_t);

    if ((global_header & 0xCAE0) != 0xCAE0) {
      std::cerr << "WARNING: Header does not match 0xCAEx pattern: 0x"
                << std::hex << global_header << std::dec << std::endl;
    }

    std::cout << "CoMPASS file opened: " << fname << std::endl;
    std::cout << "Global header: 0x" << std::hex << global_header << std::dec
              << std::endl;
    std::cout << "Control bits: " << std::bitset<4>(global_header & 0x000F)
              << std::endl;
  }

  return kTRUE;
}

Bool_t CoMPASSReader::ReadEvent() {
  if (IsEOF()) {
    return kFALSE;
  }

  current_event.header = global_header;

  file.read(reinterpret_cast<char *>(&current_event.board), sizeof(UShort_t));
  file.read(reinterpret_cast<char *>(&current_event.channel), sizeof(UShort_t));
  file.read(reinterpret_cast<char *>(&current_event.timestamp),
            sizeof(ULong64_t));
  bytes_read += sizeof(UShort_t) * 2 + sizeof(ULong64_t);

  if (current_event.hasEnergyCh()) {
    file.read(reinterpret_cast<char *>(&current_event.energy_ch),
              sizeof(UShort_t));
    bytes_read += sizeof(UShort_t);
  }

  if (current_event.hasEnergyCal()) {
    file.read(reinterpret_cast<char *>(&current_event.energy_cal),
              sizeof(Double_t));
    bytes_read += sizeof(Double_t);
  }

  if (current_event.hasEnergyShort()) {
    file.read(reinterpret_cast<char *>(&current_event.energy_short_ch),
              sizeof(UShort_t));
    bytes_read += sizeof(UShort_t);
  }

  file.read(reinterpret_cast<char *>(&current_event.flags), sizeof(UInt_t));
  bytes_read += sizeof(UInt_t);

  if (current_event.hasWaveform()) {
    file.read(reinterpret_cast<char *>(&current_event.waveform_code),
              sizeof(UChar_t));
    file.read(reinterpret_cast<char *>(&current_event.num_samples),
              sizeof(UInt_t));
    bytes_read += sizeof(UChar_t) + sizeof(UInt_t);

    // Sanity-cap num_samples; a corrupt header can read garbage and trigger
    // a multi-GB allocation in TArrayS::Set.
    const UInt_t kMaxSamples = 1u << 20;
    if (current_event.num_samples > kMaxSamples) {
      std::cerr << "ERROR: implausible num_samples "
                << current_event.num_samples << " at byte " << bytes_read
                << " (treating as truncated)" << std::endl;
      return kFALSE;
    }

    current_event.samples.Set(current_event.num_samples);
    std::vector<UShort_t> sample_buf(current_event.num_samples);
    file.read(reinterpret_cast<char *>(sample_buf.data()),
              current_event.num_samples * sizeof(UShort_t));
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete waveform at byte " << bytes_read
                << " (truncated file, event discarded)" << std::endl;
      return kFALSE;
    }
    for (UInt_t i = 0; i < current_event.num_samples; i++) {
      current_event.samples.SetAt(static_cast<Short_t>(sample_buf[i]), i);
    }
    bytes_read += current_event.num_samples * sizeof(UShort_t);
  }

  return kTRUE;
}

WaveDump742Data::WaveDump742Data()
    : event_size(0), board_id(0), pattern(0), channel(0), event_counter(0),
      group_trigger_time_tag(0), dc_offset(0), start_index_cell(0) {
  std::cout << "This version of wavedump binary conversion for 742 family "
               "digitizers is based on manual "
               "revision 21"
            << std::endl;
}

void WaveDump742Data::Print() const {
  std::cout << "WaveDump 742 Event..." << std::endl;
  std::cout << "Event size:        " << event_size << std::endl;
  std::cout << "Board ID:          " << board_id << std::endl;
  std::cout << "Channel:           " << channel << std::endl;
  std::cout << "Event counter:     " << event_counter << std::endl;
  std::cout << "Group trigger tag: " << group_trigger_time_tag << std::endl;
  std::cout << "DC offset:         " << dc_offset << std::endl;
  std::cout << "Start cell:        " << start_index_cell << std::endl;
  std::cout << "Samples:           " << samples.GetSize() << std::endl;
}

Bool_t WaveDump742Reader::ReadEvent() {
  if (IsEOF()) {
    return kFALSE;
  }

  UInt_t headers[8];
  file.read(reinterpret_cast<char *>(headers), 8 * sizeof(UInt_t));
  if (file.fail()) {
    if (file.gcount() > 0) {
      std::cerr << "WARNING: Incomplete event header at byte " << bytes_read
                << " (" << file.gcount() << " of " << 8 * sizeof(UInt_t)
                << " header bytes read, truncated file)" << std::endl;
    }
    return kFALSE;
  }
  bytes_read += 8 * sizeof(UInt_t);

  current_event.event_size = headers[0];
  current_event.board_id = headers[1];
  current_event.pattern = headers[2];
  current_event.channel = headers[3];
  current_event.event_counter = headers[4];
  current_event.group_trigger_time_tag = headers[5];
  current_event.dc_offset = headers[6];
  current_event.start_index_cell = headers[7];

  // event_size is in bytes and includes the 8-word (32-byte) header
  UInt_t header_bytes = 8 * sizeof(UInt_t);
  if (current_event.event_size <= header_bytes) {
    std::cerr << "ERROR: Invalid event_size " << current_event.event_size
              << " at event " << current_event.event_counter << std::endl;
    return kFALSE;
  }

  UInt_t sample_bytes = current_event.event_size - header_bytes;
  UInt_t sample_size = sample_bytes / sizeof(UInt_t);

  current_event.samples.Set(sample_size);

  if (corrections_enabled) {
    std::vector<Float_t> float_samples(sample_size);
    file.read(reinterpret_cast<char *>(float_samples.data()), sample_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete event at byte " << bytes_read
                << " (truncated file, event " << current_event.event_counter
                << " discarded)" << std::endl;
      return kFALSE;
    }
    for (UInt_t i = 0; i < sample_size; i++) {
      current_event.samples.SetAt(static_cast<Short_t>(float_samples[i]), i);
    }
  } else {
    std::vector<UInt_t> int_samples(sample_size);
    file.read(reinterpret_cast<char *>(int_samples.data()), sample_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete event at byte " << bytes_read
                << " (truncated file, event " << current_event.event_counter
                << " discarded)" << std::endl;
      return kFALSE;
    }
    for (UInt_t i = 0; i < sample_size; i++) {
      current_event.samples.SetAt(static_cast<Short_t>(int_samples[i] & 0xFFF),
                                  i);
    }
  }

  bytes_read += sample_bytes;

  return kTRUE;
}

// --- SOLData ---

SOLData::SOLData()
    : block_header(0), channel(0), energy(0), energy_short(0), timestamp(0),
      fine_timestamp(0), flags_high(0), flags_low(0), data_type(0),
      is_psd(kFALSE), down_sampling(0), board_fail(0), flush(0), trigger_thr(0),
      event_size(0), agg_counter(0), trace_len(0), block_id(0) {
  ana_probe_type[0] = 0xFF;
  ana_probe_type[1] = 0xFF;
  dig_probe_type[0] = 0xFF;
  dig_probe_type[1] = 0xFF;
  dig_probe_type[2] = 0xFF;
  dig_probe_type[3] = 0xFF;
}

TString SOLData::getDataTypeName() const {
  switch (data_type) {
  case ALL:
    return "ALL (full traces)";
  case OneTrace:
    return "OneTrace";
  case NoTrace:
    return "NoTrace";
  case Minimum:
    return "Minimum";
  case MiniWithFineTime:
    return "MiniWithFineTime";
  case Raw:
    return "Raw (FPGA)";
  default:
    return TString::Format("Unknown (0x%x)", data_type);
  }
}

void SOLData::Print() const {
  std::cout << "SOL block..." << std::endl;
  std::cout << "Block header:  0x" << std::hex << block_header << std::dec
            << std::endl;
  std::cout << "Data type:     " << static_cast<Int_t>(data_type) << " ("
            << getDataTypeName().Data() << ")" << std::endl;
  std::cout << "PSD:           " << (is_psd ? "YES" : "NO") << std::endl;
  std::cout << "Channel:       " << static_cast<Int_t>(channel) << std::endl;
  std::cout << "Energy:        " << energy << std::endl;
  if (is_psd) {
    std::cout << "Energy (short): " << energy_short << std::endl;
  }
  std::cout << "Timestamp:     " << timestamp << " ns" << std::endl;
  std::cout << "Fine ts:       " << fine_timestamp << " ps" << std::endl;
  std::cout << "Flags (high):  0x" << std::hex << flags_high << std::dec
            << std::endl;
  std::cout << "Flags (low):   0x" << std::hex << flags_low << std::dec
            << std::endl;

  if (data_type == ALL) {
    std::cout << "Down-sampling: " << static_cast<Int_t>(down_sampling)
              << std::endl;
    std::cout << "Board fail:    " << static_cast<Int_t>(board_fail)
              << std::endl;
    std::cout << "Flush:         " << static_cast<Int_t>(flush) << std::endl;
    std::cout << "Trigger thr:   " << trigger_thr << std::endl;
    std::cout << "Event size:    " << event_size << " bytes" << std::endl;
    std::cout << "Agg counter:   " << agg_counter << std::endl;
    std::cout << "Trace length:  " << trace_len << " samples" << std::endl;
    std::cout << "Analog probe 0: " << static_cast<Int_t>(ana_probe_type[0])
              << std::endl;
    std::cout << "Analog probe 1: " << static_cast<Int_t>(ana_probe_type[1])
              << std::endl;
    std::cout << "Digital probe 0: " << static_cast<Int_t>(dig_probe_type[0])
              << std::endl;
    std::cout << "Digital probe 1: " << static_cast<Int_t>(dig_probe_type[1])
              << std::endl;
    std::cout << "Digital probe 2: " << static_cast<Int_t>(dig_probe_type[2])
              << std::endl;
    std::cout << "Digital probe 3: " << static_cast<Int_t>(dig_probe_type[3])
              << std::endl;
  }

  if (hasTraces()) {
    std::cout << "Trace0 samples: " << trace0.size() << std::endl;
    std::cout << "Trace1 samples: " << trace1.size() << std::endl;
    if (!trace0.empty()) {
      std::cout << "Trace0 first 5: ";
      Int_t n = TMath::Min(5, static_cast<Int_t>(trace0.size()));
      for (Int_t i = 0; i < n; i++) {
        std::cout << trace0[i] << " ";
      }
      std::cout << std::endl;
    }
  }
}

// --- SOLReader ---

std::vector<TString> SOLReader::SplitSolFileByTime(const char *inputFile,
                                                   const char *outputDir,
                                                   Double_t chunkSeconds,
                                                   Int_t &totalBlocks,
                                                   Int_t &totalChunks) {
  std::vector<TString> outputFiles;
  totalBlocks = 0;
  totalChunks = 0;

  // Create output directory if it doesn't exist
  gSystem->mkdir(outputDir, kTRUE);

  // Derive base name from input file
  TString baseName = TString(inputFile);
  Int_t lastSlash = baseName.Last('/');
  if (lastSlash >= 0) {
    baseName = baseName(lastSlash + 1, baseName.Length() - lastSlash - 1);
  }
  Int_t dotPos = baseName.Last('.sol');
  if (dotPos >= 3) {
    baseName = baseName(0, dotPos - 3);
  }

  // Open input file
  std::ifstream inFile(inputFile, std::ios::binary);
  if (!inFile.is_open()) {
    std::cerr << "ERROR: Cannot open input file " << inputFile << std::endl;
    return outputFiles;
  }

  // Open first output file
  TString outPath = TString(outputDir) + "/" + baseName + "_chunk%03d.sol";
  TString currentOutPath = Form(outPath.Data(), 0);
  std::ofstream outFile(currentOutPath.Data(), std::ios::binary);

  if (!outFile.is_open()) {
    std::cerr << "ERROR: Cannot open output file " << currentOutPath
              << std::endl;
    inFile.close();
    return outputFiles;
  }

  // Read and split by time
  Long64_t chunkEndTs = 0; // ns
  Int_t chunkIndex = 0;
  Bool_t firstBlock = kTRUE;

  // Buffer for reading fixed headers (max ~100 bytes)
  std::vector<char> headerBuf(256);

  while (!inFile.eof()) {
    // Read 2-byte block header into buffer
    inFile.read(headerBuf.data(), 2);
    if (inFile.fail()) {
      break;
    }

    // Verify block start identifier
    UShort_t blockHeader;
    std::memcpy(&blockHeader, headerBuf.data(), 2);
    if ((blockHeader & 0xAA00) != 0xAA00) {
      std::cerr << "WARNING: Invalid block header 0x" << std::hex << blockHeader
                << std::dec << " at position " << inFile.tellg() << std::endl;
      break;
    }

    // Parse data type and PSD flag from header
    UChar_t dataType = blockHeader & 0xF;
    Bool_t isPsd = ((blockHeader >> 4) & 0xF) != 0;

    Long64_t timestamp = 0;
    ULong64_t traceLen = 0;
    Long64_t headerSize = 0; // size of fixed portion after the 2-byte header

    if (dataType == SOLData::ALL) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2 + 1 + 2 + 1 + 1 + 1 + 2 + 8 +
                   4 + 8 + 2 + 4;
      // Read remaining fixed header into buffer
      inFile.read(headerBuf.data() + 2, headerSize);
      if (inFile.fail()) {
        break;
      }

      // Extract timestamp from buffer (offset from start of header)
      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = 0;
      for (Int_t i = 0; i < 6; i++) {
        timestamp |= static_cast<Long64_t>(
                         static_cast<unsigned char>(headerBuf[tsOffset + i]))
                     << (static_cast<Long64_t>(i) * 8);
      }

      // Extract traceLen from buffer
      Int_t tlOffset = tsOffset + 6 + 2 + 1 + 2 + 1 + 1 + 1 + 2 + 8 + 4;
      traceLen = 0;
      for (Int_t i = 0; i < 8; i++) {
        traceLen |= static_cast<ULong64_t>(
                        static_cast<unsigned char>(headerBuf[tlOffset + i]))
                    << (static_cast<ULong64_t>(i) * 8);
      }

      // Read and write trace data if present
      if (traceLen > 0) {
        Long64_t traceBytes = traceLen * 12;
        std::vector<char> traceBuf(traceBytes);
        inFile.read(traceBuf.data(), traceBytes);
        if (inFile.fail()) {
          break;
        }

        // Write header + trace to current output
        outFile.write(headerBuf.data(), 2 + headerSize);
        outFile.write(traceBuf.data(), traceBytes);

        // Check if we need to start a new chunk
        if (firstBlock) {
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          firstBlock = kFALSE;
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        } else if (timestamp >= chunkEndTs) {
          outFile.close();
          chunkIndex++;
          currentOutPath = Form(outPath.Data(), chunkIndex);
          outFile.open(currentOutPath.Data(), std::ios::binary);
          if (!outFile.is_open()) {
            std::cerr << "ERROR: Cannot open output file " << currentOutPath
                      << std::endl;
            break;
          }
          outFile.write(headerBuf.data(), 2 + headerSize);
          outFile.write(traceBuf.data(), traceBytes);
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        }
      } else {
        // No trace data, just write header
        outFile.write(headerBuf.data(), 2 + headerSize);

        if (firstBlock) {
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          firstBlock = kFALSE;
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        } else if (timestamp >= chunkEndTs) {
          outFile.close();
          chunkIndex++;
          currentOutPath = Form(outPath.Data(), chunkIndex);
          outFile.open(currentOutPath.Data(), std::ios::binary);
          if (!outFile.is_open()) {
            std::cerr << "ERROR: Cannot open output file " << currentOutPath
                      << std::endl;
            break;
          }
          outFile.write(headerBuf.data(), 2 + headerSize);
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        }
      }

    } else if (dataType == SOLData::OneTrace) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2 + 1 + 2 + 8 + 1;
      inFile.read(headerBuf.data() + 2, headerSize);
      if (inFile.fail()) {
        break;
      }

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = 0;
      for (Int_t i = 0; i < 6; i++) {
        timestamp |= static_cast<Long64_t>(
                         static_cast<unsigned char>(headerBuf[tsOffset + i]))
                     << (static_cast<Long64_t>(i) * 8);
      }

      Int_t tlOffset = tsOffset + 6 + 2 + 1 + 2;
      traceLen = 0;
      for (Int_t i = 0; i < 8; i++) {
        traceLen |= static_cast<ULong64_t>(
                        static_cast<unsigned char>(headerBuf[tlOffset + i]))
                    << (static_cast<ULong64_t>(i) * 8);
      }

      if (traceLen > 0) {
        Long64_t traceBytes = traceLen * sizeof(Int_t);
        std::vector<char> traceBuf(traceBytes);
        inFile.read(traceBuf.data(), traceBytes);
        if (inFile.fail()) {
          break;
        }

        outFile.write(headerBuf.data(), 2 + headerSize);
        outFile.write(traceBuf.data(), traceBytes);

        if (firstBlock) {
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          firstBlock = kFALSE;
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        } else if (timestamp >= chunkEndTs) {
          outFile.close();
          chunkIndex++;
          currentOutPath = Form(outPath.Data(), chunkIndex);
          outFile.open(currentOutPath.Data(), std::ios::binary);
          if (!outFile.is_open()) {
            std::cerr << "ERROR: Cannot open output file " << currentOutPath
                      << std::endl;
            break;
          }
          outFile.write(headerBuf.data(), 2 + headerSize);
          outFile.write(traceBuf.data(), traceBytes);
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        }
      } else {
        outFile.write(headerBuf.data(), 2 + headerSize);

        if (firstBlock) {
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          firstBlock = kFALSE;
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        } else if (timestamp >= chunkEndTs) {
          outFile.close();
          chunkIndex++;
          currentOutPath = Form(outPath.Data(), chunkIndex);
          outFile.open(currentOutPath.Data(), std::ios::binary);
          if (!outFile.is_open()) {
            std::cerr << "ERROR: Cannot open output file " << currentOutPath
                      << std::endl;
            break;
          }
          outFile.write(headerBuf.data(), 2 + headerSize);
          chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
          outputFiles.push_back(currentOutPath);
          totalChunks++;
        }
      }

    } else if (dataType == SOLData::NoTrace) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2 + 1 + 2;
      inFile.read(headerBuf.data() + 2, headerSize);
      if (inFile.fail()) {
        break;
      }

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = 0;
      for (Int_t i = 0; i < 6; i++) {
        timestamp |= static_cast<Long64_t>(
                         static_cast<unsigned char>(headerBuf[tsOffset + i]))
                     << (static_cast<Long64_t>(i) * 8);
      }

      outFile.write(headerBuf.data(), 2 + headerSize);

      if (firstBlock) {
        chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
        firstBlock = kFALSE;
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      } else if (timestamp >= chunkEndTs) {
        outFile.close();
        chunkIndex++;
        currentOutPath = Form(outPath.Data(), chunkIndex);
        outFile.open(currentOutPath.Data(), std::ios::binary);
        if (!outFile.is_open()) {
          std::cerr << "ERROR: Cannot open output file " << currentOutPath
                    << std::endl;
          break;
        }
        outFile.write(headerBuf.data(), 2 + headerSize);
        chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      }

    } else if (dataType == SOLData::Minimum) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6;
      inFile.read(headerBuf.data() + 2, headerSize);
      if (inFile.fail()) {
        break;
      }

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = 0;
      for (Int_t i = 0; i < 6; i++) {
        timestamp |= static_cast<Long64_t>(
                         static_cast<unsigned char>(headerBuf[tsOffset + i]))
                     << (static_cast<Long64_t>(i) * 8);
      }

      outFile.write(headerBuf.data(), 2 + headerSize);

      if (firstBlock) {
        chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
        firstBlock = kFALSE;
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      } else if (timestamp >= chunkEndTs) {
        outFile.close();
        chunkIndex++;
        currentOutPath = Form(outPath.Data(), chunkIndex);
        outFile.open(currentOutPath.Data(), std::ios::binary);
        if (!outFile.is_open()) {
          std::cerr << "ERROR: Cannot open output file " << currentOutPath
                    << std::endl;
          break;
        }
        outFile.write(headerBuf.data(), 2 + headerSize);
        chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      }

    } else if (dataType == SOLData::MiniWithFineTime) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2;
      inFile.read(headerBuf.data() + 2, headerSize);
      if (inFile.fail()) {
        break;
      }

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = 0;
      for (Int_t i = 0; i < 6; i++) {
        timestamp |= static_cast<Long64_t>(
                         static_cast<unsigned char>(headerBuf[tsOffset + i]))
                     << (static_cast<Long64_t>(i) * 8);
      }

      outFile.write(headerBuf.data(), 2 + headerSize);

      if (firstBlock) {
        chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
        firstBlock = kFALSE;
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      } else if (timestamp >= chunkEndTs) {
        outFile.close();
        chunkIndex++;
        currentOutPath = Form(outPath.Data(), chunkIndex);
        outFile.open(currentOutPath.Data(), std::ios::binary);
        if (!outFile.is_open()) {
          std::cerr << "ERROR: Cannot open output file " << currentOutPath
                    << std::endl;
          break;
        }
        outFile.write(headerBuf.data(), 2 + headerSize);
        chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      }

    } else if (dataType == SOLData::Raw) {
      headerSize = 8;
      inFile.read(headerBuf.data() + 2, headerSize);
      if (inFile.fail()) {
        break;
      }

      ULong64_t rawDataSize = 0;
      for (Int_t i = 0; i < 8; i++) {
        rawDataSize |=
            static_cast<ULong64_t>(static_cast<unsigned char>(headerBuf[2 + i]))
            << (static_cast<ULong64_t>(i) * 8);
      }
      outFile.write(headerBuf.data(), 2 + headerSize);

      if (rawDataSize > 0) {
        std::vector<char> rawBuf(rawDataSize);
        inFile.read(rawBuf.data(), rawDataSize);
        if (inFile.fail()) {
          break;
        }
        outFile.write(rawBuf.data(), rawDataSize);
      }

      if (firstBlock) {
        chunkEndTs = static_cast<Long64_t>(chunkSeconds * 1e9);
        firstBlock = kFALSE;
        outputFiles.push_back(currentOutPath);
        totalChunks++;
      }

    } else {
      std::cerr << "ERROR: Unknown SOL data type 0x" << std::hex << dataType
                << std::dec << " at position " << inFile.tellg() - 2
                << std::endl;
      break;
    }

    totalBlocks++;
  }

  // Clean up
  outFile.close();
  inFile.close();

  std::cout << "Split " << totalBlocks << " blocks into " << totalChunks
            << " chunks in " << outputDir << std::endl;

  return outputFiles;
}

Bool_t SOLReader::ReadEvent() {
  if (IsEOF()) {
    return kFALSE;
  }

  // Reset current event
  current_event = SOLData();

  // Read 2-byte block header
  UShort_t block_header;
  file.read(reinterpret_cast<char *>(&block_header), sizeof(UShort_t));
  if (file.fail()) {
    return kFALSE;
  }
  bytes_read += sizeof(UShort_t);

  // Verify block start identifier
  if ((block_header & 0xAA00) != 0xAA00) {
    std::cerr << "WARNING: SOL block header mismatch at block " << block_id
              << " (got 0x" << std::hex << block_header << std::dec
              << ", expected 0xAAxx)" << std::endl;
    return kFALSE;
  }

  current_event.block_header = block_header;
  current_event.data_type = block_header & 0xF;
  current_event.is_psd = ((block_header >> 4) & 0xF) != 0;

  // Parse fields based on data type
  if (current_event.data_type == SOLData::ALL) {
    file.read(reinterpret_cast<char *>(&current_event.channel), 1);
    file.read(reinterpret_cast<char *>(&current_event.energy), 2);
    if (current_event.is_psd) {
      file.read(reinterpret_cast<char *>(&current_event.energy_short), 2);
    }
    file.read(reinterpret_cast<char *>(&current_event.timestamp), 6);
    file.read(reinterpret_cast<char *>(&current_event.fine_timestamp), 2);
    file.read(reinterpret_cast<char *>(&current_event.flags_high), 1);
    file.read(reinterpret_cast<char *>(&current_event.flags_low), 2);
    file.read(reinterpret_cast<char *>(&current_event.down_sampling), 1);
    file.read(reinterpret_cast<char *>(&current_event.board_fail), 1);
    file.read(reinterpret_cast<char *>(&current_event.flush), 1);
    file.read(reinterpret_cast<char *>(&current_event.trigger_thr), 2);
    file.read(reinterpret_cast<char *>(&current_event.event_size), 8);
    file.read(reinterpret_cast<char *>(&current_event.agg_counter), 4);
    file.read(reinterpret_cast<char *>(&current_event.trace_len), 8);
    file.read(reinterpret_cast<char *>(current_event.ana_probe_type), 2);
    file.read(reinterpret_cast<char *>(current_event.dig_probe_type), 4);

    if (file.fail()) {
      std::cerr << "WARNING: Incomplete ALL-format block at byte " << bytes_read
                << " (truncated file)" << std::endl;
      return kFALSE;
    }

    // Advance bytes_read for the fixed header portion
    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 +
                         2 + 1 + 1 + 1 + 2 + 8 + 4 + 8 + 2 + 4;
    bytes_read += header_bytes;

    // Read or skip trace data
    if (current_event.trace_len > 0) {
      ULong64_t n_samples = current_event.trace_len;

      // Cap to reasonable limit to avoid OOM on corrupt data
      const ULong64_t kMaxTraceLen = 100000;
      if (n_samples > kMaxTraceLen) {
        std::cerr << "WARNING: Implausible trace length " << n_samples
                  << " at block " << block_id << " (skipping trace data)"
                  << std::endl;
        // Skip the trace data (2 analog traces * 4 bytes + 4 digital * 1 byte =
        // 12 bytes/sample)
        Long64_t skip_bytes = n_samples * 12;
        file.seekg(skip_bytes, std::ios::cur);
        bytes_read += skip_bytes;
      } else if (skip_traces) {
        // Skip trace data without allocating (2 analog * 4B + 4 digital * 1B =
        // 12B/sample)
        Long64_t skip_bytes = n_samples * 12;
        file.seekg(skip_bytes, std::ios::cur);
        bytes_read += skip_bytes;
      } else {
        current_event.trace0.resize(n_samples);
        current_event.trace1.resize(n_samples);
        current_event.dig0.resize(n_samples);
        current_event.dig1.resize(n_samples);
        current_event.dig2.resize(n_samples);
        current_event.dig3.resize(n_samples);

        // Read analog traces (4 bytes per sample)
        file.read(reinterpret_cast<char *>(current_event.trace0.data()),
                  n_samples * sizeof(Int_t));
        if (file.fail() || file.gcount() != static_cast<std::streamsize>(
                                                n_samples * sizeof(Int_t))) {
          std::cerr << "WARNING: Incomplete trace0 at block " << block_id
                    << " (truncated file)" << std::endl;
          return kFALSE;
        }
        file.read(reinterpret_cast<char *>(current_event.trace1.data()),
                  n_samples * sizeof(Int_t));
        if (file.fail() || file.gcount() != static_cast<std::streamsize>(
                                                n_samples * sizeof(Int_t))) {
          std::cerr << "WARNING: Incomplete trace1 at block " << block_id
                    << " (truncated file)" << std::endl;
          return kFALSE;
        }

        // Read digital traces (1 byte per sample)
        file.read(reinterpret_cast<char *>(current_event.dig0.data()),
                  n_samples);
        file.read(reinterpret_cast<char *>(current_event.dig1.data()),
                  n_samples);
        file.read(reinterpret_cast<char *>(current_event.dig2.data()),
                  n_samples);
        file.read(reinterpret_cast<char *>(current_event.dig3.data()),
                  n_samples);
        if (file.fail()) {
          std::cerr << "WARNING: Incomplete digital traces at block "
                    << block_id << " (truncated file)" << std::endl;
          return kFALSE;
        }

        bytes_read += n_samples * (2 * sizeof(Int_t) + 4);
      }
    }

  } else if (current_event.data_type == SOLData::OneTrace) {
    file.read(reinterpret_cast<char *>(&current_event.channel), 1);
    file.read(reinterpret_cast<char *>(&current_event.energy), 2);
    if (current_event.is_psd) {
      file.read(reinterpret_cast<char *>(&current_event.energy_short), 2);
    }
    file.read(reinterpret_cast<char *>(&current_event.timestamp), 6);
    file.read(reinterpret_cast<char *>(&current_event.fine_timestamp), 2);
    file.read(reinterpret_cast<char *>(&current_event.flags_high), 1);
    file.read(reinterpret_cast<char *>(&current_event.flags_low), 2);
    file.read(reinterpret_cast<char *>(&current_event.trace_len), 8);
    file.read(reinterpret_cast<char *>(&current_event.ana_probe_type[0]), 1);

    if (file.fail()) {
      std::cerr << "WARNING: Incomplete OneTrace-format block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }

    Int_t header_bytes =
        1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 + 2 + 8 + 1;
    bytes_read += header_bytes;

    // Read or skip single analog trace
    if (current_event.trace_len > 0) {
      ULong64_t n_samples = current_event.trace_len;
      const ULong64_t kMaxTraceLen = 100000;

      if (n_samples > kMaxTraceLen) {
        std::cerr << "WARNING: Implausible trace length " << n_samples
                  << " at block " << block_id << " (skipping trace data)"
                  << std::endl;
        file.seekg(n_samples * sizeof(Int_t), std::ios::cur);
        bytes_read += n_samples * sizeof(Int_t);
      } else if (skip_traces) {
        file.seekg(n_samples * sizeof(Int_t), std::ios::cur);
        bytes_read += n_samples * sizeof(Int_t);
      } else {
        current_event.trace0.resize(n_samples);
        file.read(reinterpret_cast<char *>(current_event.trace0.data()),
                  n_samples * sizeof(Int_t));
        if (file.fail()) {
          std::cerr << "WARNING: Incomplete trace at block " << block_id
                    << " (truncated file)" << std::endl;
          return kFALSE;
        }
        bytes_read += n_samples * sizeof(Int_t);
      }
    }

  } else if (current_event.data_type == SOLData::NoTrace) {
    file.read(reinterpret_cast<char *>(&current_event.channel), 1);
    file.read(reinterpret_cast<char *>(&current_event.energy), 2);
    if (current_event.is_psd) {
      file.read(reinterpret_cast<char *>(&current_event.energy_short), 2);
    }
    file.read(reinterpret_cast<char *>(&current_event.timestamp), 6);
    file.read(reinterpret_cast<char *>(&current_event.fine_timestamp), 2);
    file.read(reinterpret_cast<char *>(&current_event.flags_high), 1);
    file.read(reinterpret_cast<char *>(&current_event.flags_low), 2);

    if (file.fail()) {
      std::cerr << "WARNING: Incomplete NoTrace-format block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }

    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 + 2;
    bytes_read += header_bytes;

    // Skip remaining bytes (down_sampling, board_fail, flush, trigger_thr,
    // event_size, agg_counter, etc.)
    Long64_t remaining =
        (current_event.is_psd ? 20 : 18) -
        (1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 + 2);
    if (remaining > 0) {
      file.seekg(remaining, std::ios::cur);
      bytes_read += remaining;
    }

  } else if (current_event.data_type == SOLData::MiniWithFineTime) {
    file.read(reinterpret_cast<char *>(&current_event.channel), 1);
    file.read(reinterpret_cast<char *>(&current_event.energy), 2);
    if (current_event.is_psd) {
      file.read(reinterpret_cast<char *>(&current_event.energy_short), 2);
    }
    file.read(reinterpret_cast<char *>(&current_event.timestamp), 6);
    file.read(reinterpret_cast<char *>(&current_event.fine_timestamp), 2);

    if (file.fail()) {
      std::cerr << "WARNING: Incomplete MiniWithFineTime block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }

    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2;
    bytes_read += header_bytes;

  } else if (current_event.data_type == SOLData::Minimum) {
    file.read(reinterpret_cast<char *>(&current_event.channel), 1);
    file.read(reinterpret_cast<char *>(&current_event.energy), 2);
    if (current_event.is_psd) {
      file.read(reinterpret_cast<char *>(&current_event.energy_short), 2);
    }
    file.read(reinterpret_cast<char *>(&current_event.timestamp), 6);

    if (file.fail()) {
      std::cerr << "WARNING: Incomplete Minimum-format block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }

    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6;
    bytes_read += header_bytes;

  } else if (current_event.data_type == SOLData::Raw) {
    ULong64_t data_size;
    file.read(reinterpret_cast<char *>(&data_size), 8);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete Raw-format block at byte " << bytes_read
                << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += 8;

    // Skip raw FPGA data
    if (data_size > 0) {
      file.seekg(data_size, std::ios::cur);
      bytes_read += data_size;
    }

  } else {
    std::cerr << "ERROR: Unknown SOL data type 0x" << std::hex
              << current_event.data_type << std::dec << " at block " << block_id
              << std::endl;
    return kFALSE;
  }

  current_event.block_id = block_id;
  block_id++;

  return kTRUE;
}

SOLHit SOLReader::ToHit() const {
  SOLHit hit;
  hit.channel = static_cast<Short_t>(current_event.channel);
  hit.energy = static_cast<Short_t>(current_event.energy);
  hit.energy_short = static_cast<Short_t>(current_event.energy_short);
  hit.timestamp = current_event.timestamp;
  hit.fine_timestamp = static_cast<Short_t>(current_event.fine_timestamp);
  hit.flags_high = static_cast<Short_t>(current_event.flags_high);
  hit.flags_low = static_cast<Short_t>(current_event.flags_low);
  hit.data_type = static_cast<Short_t>(current_event.data_type);
  hit.is_psd = static_cast<Char_t>(current_event.is_psd);
  hit.down_sampling = static_cast<Char_t>(current_event.down_sampling);
  hit.board_fail = static_cast<Char_t>(current_event.board_fail);
  hit.flush = static_cast<Char_t>(current_event.flush);
  hit.trigger_thr = static_cast<Short_t>(current_event.trigger_thr);
  hit.event_size = current_event.event_size;
  hit.agg_counter = static_cast<Int_t>(current_event.agg_counter);
  hit.block_id = current_event.block_id;
  hit.trace_len = static_cast<Int_t>(current_event.trace_len);
  hit.ana_probe_type[0] = static_cast<Char_t>(current_event.ana_probe_type[0]);
  hit.ana_probe_type[1] = static_cast<Char_t>(current_event.ana_probe_type[1]);
  hit.dig_probe_type[0] = static_cast<Char_t>(current_event.dig_probe_type[0]);
  hit.dig_probe_type[1] = static_cast<Char_t>(current_event.dig_probe_type[1]);
  hit.dig_probe_type[2] = static_cast<Char_t>(current_event.dig_probe_type[2]);
  hit.dig_probe_type[3] = static_cast<Char_t>(current_event.dig_probe_type[3]);
  return hit;
}
