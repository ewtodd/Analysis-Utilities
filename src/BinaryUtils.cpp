#include "BinaryUtils.hpp"
#include <cstring>

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
    const Int_t *a0 = getAnalog0();
    const Int_t *a1 = getAnalog1();
    if (data_type == ALL && a0) {
      std::cout << "Trace0 samples: " << getSamples() << std::endl;
      std::cout << "Trace1 samples: " << getSamples() << std::endl;
      UInt_t n = TMath::Min(5, static_cast<Int_t>(getSamples()));
      std::cout << "Trace0 first 5: ";
      for (UInt_t i = 0; i < n; i++) {
        std::cout << a0[i] << " ";
      }
      std::cout << std::endl;
      std::cout << "Trace1 first 5: ";
      for (UInt_t i = 0; i < n; i++) {
        std::cout << a1[i] << " ";
      }
      std::cout << std::endl;
    } else if (data_type == OneTrace) {
      const Int_t *ot = getOneTrace();
      if (ot) {
        std::cout << "Trace0 samples: " << getSamples() << std::endl;
        UInt_t n = TMath::Min(5, static_cast<Int_t>(getSamples()));
        std::cout << "Trace0 first 5: ";
        for (UInt_t i = 0; i < n; i++) {
          std::cout << ot[i] << " ";
        }
        std::cout << std::endl;
      }
    }
  }
}

// --- SOLReader ---

// Little-endian reads from a block buffer, mirroring the memcpy+shift
// pattern used elsewhere in this file.
static Long64_t ReadLe6(const char *p) {
  Long64_t v = 0;
  for (Int_t i = 0; i < 6; i++)
    v |= static_cast<Long64_t>(static_cast<unsigned char>(p[i]))
         << (static_cast<Long64_t>(i) * 8);
  return v;
}

static ULong64_t ReadLe8(const char *p) {
  ULong64_t v = 0;
  for (Int_t i = 0; i < 8; i++)
    v |= static_cast<ULong64_t>(static_cast<unsigned char>(p[i]))
         << (static_cast<ULong64_t>(i) * 8);
  return v;
}

// Keep at least needBytes valid bytes at the front of buf, compacting the
// unconsumed tail and reading from inFile as required. Returns kFALSE when
// the file is exhausted before the request can be satisfied.
static Bool_t RefillInput(std::ifstream &inFile, std::vector<char> &buf,
                          Int_t &pos, Int_t &len, Long64_t needBytes) {
  if (Long64_t(len - pos) >= needBytes)
    return kTRUE;
  if (pos > 0) {
    std::memmove(buf.data(), buf.data() + pos, size_t(len - pos));
    len -= pos;
    pos = 0;
  }
  if (Long64_t(buf.size()) < needBytes)
    buf.resize(size_t(needBytes));
  while (Long64_t(len) < needBytes) {
    inFile.read(buf.data() + len, Long64_t(buf.size()) - len);
    Long64_t got = inFile.gcount();
    len += Int_t(got);
    if (got <= 0)
      return kFALSE;
  }
  return kTRUE;
}

std::vector<TString> SOLReader::SplitSolFileByTime(const char *inputFile,
                                                   const char *outputDir,
                                                   Double_t chunkSeconds,
                                                   Int_t &totalBlocks,
                                                   Int_t &totalChunks) {
  std::vector<TString> outputFiles;
  totalBlocks = 0;
  totalChunks = 0;

  gSystem->mkdir(outputDir, kTRUE);

  TString baseName = TString(inputFile);
  Int_t lastSlash = baseName.Last('/');
  if (lastSlash >= 0) {
    baseName = baseName(lastSlash + 1, baseName.Length() - lastSlash - 1);
  }
  Int_t dotPos = baseName.Last('.sol');
  if (dotPos >= 3) {
    baseName = baseName(0, dotPos - 3);
  }

  std::ifstream inFile(inputFile, std::ios::binary);
  if (!inFile.is_open()) {
    std::cerr << "ERROR: Cannot open input file " << inputFile << std::endl;
    return outputFiles;
  }

  TString outPath = TString(outputDir) + "/" + baseName + "_chunk%03d.sol";
  TString currentOutPath = Form(outPath.Data(), 0);
  std::ofstream outFile(currentOutPath.Data(), std::ios::binary);

  if (!outFile.is_open()) {
    std::cerr << "ERROR: Cannot open output file " << currentOutPath
              << std::endl;
    inFile.close();
    return outputFiles;
  }

  Long64_t chunkEndTs = 0;
  Int_t chunkIndex = 0;
  Bool_t firstBlock = kTRUE;

  // Application-level I/O batching: blocks are parsed straight out of a
  // large input buffer and staged in a large output buffer, so per-block
  // parsing never costs individual stream calls. Minimum-format files are
  // ~11 bytes per block; without batching the per-block stream overhead
  // caps throughput around 15 MB/s instead of the disk's ~2 GB/s.
  const Int_t kBulkRead = 4 * 1024 * 1024;
  const Int_t kBulkWrite = 4 * 1024 * 1024;
  std::vector<char> ibuf(kBulkRead);
  std::vector<char> obuf;
  obuf.reserve(kBulkWrite);
  Int_t iPos = 0;       // next unconsumed byte in ibuf
  Int_t iLen = 0;       // valid bytes in ibuf
  Long64_t filePos = 0; // absolute file offset of ibuf[iPos]

  while (RefillInput(inFile, ibuf, iPos, iLen, 2)) {
    UShort_t blockHeader;
    std::memcpy(&blockHeader, ibuf.data() + iPos, 2);
    if ((blockHeader & 0xAA00) != 0xAA00) {
      std::cerr << "WARNING: Invalid block header 0x" << std::hex << blockHeader
                << std::dec << " at position " << filePos + 2 << std::endl;
      break;
    }

    UChar_t dataType = blockHeader & 0xF;
    Bool_t isPsd = ((blockHeader >> 4) & 0xF) != 0;

    Long64_t headerSize = 0;
    Long64_t timestamp = 0;
    Bool_t hasTimestamp = kTRUE;
    Long64_t dataBytes = 0;

    if (dataType == SOLData::ALL) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2 + 1 + 2 + 1 + 1 + 1 + 2 + 8 +
                   4 + 8 + 2 + 4;
      if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize))
        break;

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = ReadLe6(ibuf.data() + iPos + tsOffset);

      Int_t tlOffset = tsOffset + 6 + 2 + 1 + 2 + 1 + 1 + 1 + 2 + 8 + 4;
      ULong64_t traceLen = ReadLe8(ibuf.data() + iPos + tlOffset);
      if (traceLen > 0) {
        dataBytes = traceLen * 12;
      }

    } else if (dataType == SOLData::OneTrace) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2 + 1 + 2 + 8 + 1;
      if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize))
        break;

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = ReadLe6(ibuf.data() + iPos + tsOffset);

      Int_t tlOffset = tsOffset + 6 + 2 + 1 + 2;
      ULong64_t traceLen = ReadLe8(ibuf.data() + iPos + tlOffset);
      if (traceLen > 0) {
        dataBytes = traceLen * sizeof(Int_t);
      }

    } else if (dataType == SOLData::NoTrace) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2 + 1 + 2;
      if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize))
        break;

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = ReadLe6(ibuf.data() + iPos + tsOffset);

    } else if (dataType == SOLData::Minimum) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6;
      if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize))
        break;

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = ReadLe6(ibuf.data() + iPos + tsOffset);

    } else if (dataType == SOLData::MiniWithFineTime) {
      headerSize = 1 + 2 + (isPsd ? 2 : 0) + 6 + 2;
      if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize))
        break;

      Int_t tsOffset = 2 + 1 + 2 + (isPsd ? 2 : 0);
      timestamp = ReadLe6(ibuf.data() + iPos + tsOffset);

    } else if (dataType == SOLData::Raw) {
      headerSize = 8;
      if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize))
        break;

      hasTimestamp = kFALSE;
      dataBytes = ReadLe8(ibuf.data() + iPos + 2);

    } else {
      std::cerr << "ERROR: Unknown SOL data type 0x" << std::hex << dataType
                << std::dec << " at position " << filePos << std::endl;
      break;
    }

    // Ensure the full block (header + payload) is available before writing.
    if (!RefillInput(inFile, ibuf, iPos, iLen, 2 + headerSize + dataBytes))
      break;

    // Rotate to new chunk if timestamp crosses boundary.
    // Block that crosses boundary goes into the NEW chunk only.
    if (hasTimestamp && (firstBlock || timestamp >= chunkEndTs)) {
      if (!firstBlock) {
        if (!obuf.empty()) {
          outFile.write(obuf.data(), obuf.size());
          obuf.clear();
        }
        outFile.close();
        chunkIndex++;
        currentOutPath = Form(outPath.Data(), chunkIndex);
        outFile.open(currentOutPath.Data(), std::ios::binary);
        if (!outFile.is_open()) {
          std::cerr << "ERROR: Cannot open output file " << currentOutPath
                    << std::endl;
          break;
        }
      }
      chunkEndTs = timestamp + static_cast<Long64_t>(chunkSeconds * 1e9);
      firstBlock = kFALSE;
      outputFiles.push_back(currentOutPath);
      totalChunks++;
    } else if (!hasTimestamp && firstBlock) {
      firstBlock = kFALSE;
      outputFiles.push_back(currentOutPath);
      totalChunks++;
    }

    // Stage the block in the write buffer, flushing when it grows large.
    const Int_t blockBytes = Int_t(2 + headerSize + dataBytes);
    obuf.insert(obuf.end(), ibuf.begin() + iPos,
                ibuf.begin() + iPos + blockBytes);
    iPos += blockBytes;
    filePos += blockBytes;
    if (Long64_t(obuf.size()) >= kBulkWrite) {
      outFile.write(obuf.data(), obuf.size());
      obuf.clear();
    }

    totalBlocks++;
  }

  // Flush any staged output before closing.
  if (!obuf.empty()) {
    outFile.write(obuf.data(), obuf.size());
    obuf.clear();
  }
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

  current_event = SOLData();

  UShort_t block_header;
  file.read(reinterpret_cast<char *>(&block_header), sizeof(UShort_t));
  if (file.fail()) {
    return kFALSE;
  }
  bytes_read += sizeof(UShort_t);

  if ((block_header & 0xAA00) != 0xAA00) {
    std::cerr << "WARNING: SOL block header mismatch at block " << block_id
              << " (got 0x" << std::hex << block_header << std::dec
              << ", expected 0xAAxx)" << std::endl;
    return kFALSE;
  }

  current_event.block_header = block_header;
  current_event.data_type = block_header & 0xF;
  current_event.is_psd = ((block_header >> 4) & 0xF) != 0;

  if (current_event.data_type == SOLData::ALL) {
    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 +
                         2 + 1 + 1 + 1 + 2 + 8 + 4 + 8 + 2 + 4;
    std::vector<char> hdr(header_bytes);
    file.read(hdr.data(), header_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete ALL-format block at byte " << bytes_read
                << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += header_bytes;

    Int_t off = 0;
    current_event.channel = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.energy, hdr.data() + off, 2);
    off += 2;
    if (current_event.is_psd) {
      std::memcpy(&current_event.energy_short, hdr.data() + off, 2);
      off += 2;
    }
    std::memcpy(&current_event.timestamp, hdr.data() + off, 6);
    off += 6;
    std::memcpy(&current_event.fine_timestamp, hdr.data() + off, 2);
    off += 2;
    current_event.flags_high = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.flags_low, hdr.data() + off, 2);
    off += 2;
    current_event.down_sampling = static_cast<UChar_t>(hdr[off]);
    off += 1;
    current_event.board_fail = static_cast<UChar_t>(hdr[off]);
    off += 1;
    current_event.flush = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.trigger_thr, hdr.data() + off, 2);
    off += 2;
    std::memcpy(&current_event.event_size, hdr.data() + off, 8);
    off += 8;
    std::memcpy(&current_event.agg_counter, hdr.data() + off, 4);
    off += 4;
    std::memcpy(&current_event.trace_len, hdr.data() + off, 8);
    off += 8;
    std::memcpy(current_event.ana_probe_type, hdr.data() + off, 2);
    off += 2;
    std::memcpy(current_event.dig_probe_type, hdr.data() + off, 4);

    if (current_event.trace_len > 0) {
      ULong64_t n_samples = current_event.trace_len;
      const ULong64_t kMaxTraceLen = 100000;
      Long64_t trace_bytes = n_samples * 12;

      if (n_samples > kMaxTraceLen) {
        std::cerr << "WARNING: Implausible trace length " << n_samples
                  << " at block " << block_id << " (skipping trace data)"
                  << std::endl;
        file.seekg(trace_bytes, std::ios::cur);
        bytes_read += trace_bytes;
      } else if (skip_traces) {
        file.seekg(trace_bytes, std::ios::cur);
        bytes_read += trace_bytes;
      } else {
        current_event.trace_data.resize(trace_bytes);
        file.read(current_event.trace_data.data(), trace_bytes);
        if (file.fail() ||
            file.gcount() != static_cast<std::streamsize>(trace_bytes)) {
          std::cerr << "WARNING: Incomplete traces at block " << block_id
                    << " (truncated file)" << std::endl;
          return kFALSE;
        }
        bytes_read += trace_bytes;
      }
    }

  } else if (current_event.data_type == SOLData::OneTrace) {
    Int_t header_bytes =
        1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 + 2 + 8 + 1;
    std::vector<char> hdr(header_bytes);
    file.read(hdr.data(), header_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete OneTrace-format block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += header_bytes;

    Int_t off = 0;
    current_event.channel = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.energy, hdr.data() + off, 2);
    off += 2;
    if (current_event.is_psd) {
      std::memcpy(&current_event.energy_short, hdr.data() + off, 2);
      off += 2;
    }
    std::memcpy(&current_event.timestamp, hdr.data() + off, 6);
    off += 6;
    std::memcpy(&current_event.fine_timestamp, hdr.data() + off, 2);
    off += 2;
    current_event.flags_high = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.flags_low, hdr.data() + off, 2);
    off += 2;
    std::memcpy(&current_event.trace_len, hdr.data() + off, 8);
    off += 8;
    current_event.ana_probe_type[0] = static_cast<UChar_t>(hdr[off]);

    if (current_event.trace_len > 0) {
      ULong64_t n_samples = current_event.trace_len;
      const ULong64_t kMaxTraceLen = 100000;
      Long64_t trace_bytes = n_samples * sizeof(Int_t);

      if (n_samples > kMaxTraceLen) {
        std::cerr << "WARNING: Implausible trace length " << n_samples
                  << " at block " << block_id << " (skipping trace data)"
                  << std::endl;
        file.seekg(trace_bytes, std::ios::cur);
        bytes_read += trace_bytes;
      } else if (skip_traces) {
        file.seekg(trace_bytes, std::ios::cur);
        bytes_read += trace_bytes;
      } else {
        current_event.trace_data.resize(trace_bytes);
        file.read(current_event.trace_data.data(), trace_bytes);
        if (file.fail()) {
          std::cerr << "WARNING: Incomplete trace at block " << block_id
                    << " (truncated file)" << std::endl;
          return kFALSE;
        }
        bytes_read += trace_bytes;
      }
    }

  } else if (current_event.data_type == SOLData::NoTrace) {
    Int_t header_bytes =
        1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 + 2 +
        (current_event.is_psd ? 20 : 18) -
        (1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2 + 1 + 2);
    std::vector<char> hdr(header_bytes);
    file.read(hdr.data(), header_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete NoTrace-format block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += header_bytes;

    Int_t off = 0;
    current_event.channel = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.energy, hdr.data() + off, 2);
    off += 2;
    if (current_event.is_psd) {
      std::memcpy(&current_event.energy_short, hdr.data() + off, 2);
      off += 2;
    }
    std::memcpy(&current_event.timestamp, hdr.data() + off, 6);
    off += 6;
    std::memcpy(&current_event.fine_timestamp, hdr.data() + off, 2);
    off += 2;
    current_event.flags_high = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.flags_low, hdr.data() + off, 2);

  } else if (current_event.data_type == SOLData::MiniWithFineTime) {
    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6 + 2;
    std::vector<char> hdr(header_bytes);
    file.read(hdr.data(), header_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete MiniWithFineTime block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += header_bytes;

    Int_t off = 0;
    current_event.channel = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.energy, hdr.data() + off, 2);
    off += 2;
    if (current_event.is_psd) {
      std::memcpy(&current_event.energy_short, hdr.data() + off, 2);
      off += 2;
    }
    std::memcpy(&current_event.timestamp, hdr.data() + off, 6);
    off += 6;
    std::memcpy(&current_event.fine_timestamp, hdr.data() + off, 2);

  } else if (current_event.data_type == SOLData::Minimum) {
    Int_t header_bytes = 1 + 2 + (current_event.is_psd ? 2 : 0) + 6;
    std::vector<char> hdr(header_bytes);
    file.read(hdr.data(), header_bytes);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete Minimum-format block at byte "
                << bytes_read << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += header_bytes;

    Int_t off = 0;
    current_event.channel = static_cast<UChar_t>(hdr[off]);
    off += 1;
    std::memcpy(&current_event.energy, hdr.data() + off, 2);
    off += 2;
    if (current_event.is_psd) {
      std::memcpy(&current_event.energy_short, hdr.data() + off, 2);
      off += 2;
    }
    std::memcpy(&current_event.timestamp, hdr.data() + off, 6);

  } else if (current_event.data_type == SOLData::Raw) {
    ULong64_t data_size;
    file.read(reinterpret_cast<char *>(&data_size), 8);
    if (file.fail()) {
      std::cerr << "WARNING: Incomplete Raw-format block at byte " << bytes_read
                << " (truncated file)" << std::endl;
      return kFALSE;
    }
    bytes_read += 8;

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
  hit.channel = current_event.channel;
  hit.energy = current_event.energy;
  hit.energy_short = current_event.energy_short;
  hit.timestamp = current_event.timestamp;
  hit.fine_timestamp = current_event.fine_timestamp;
  hit.flags_high = current_event.flags_high;
  hit.flags_low = current_event.flags_low;
  hit.data_type = current_event.data_type;
  hit.is_psd = current_event.is_psd;
  hit.down_sampling = current_event.down_sampling;
  hit.board_fail = current_event.board_fail;
  hit.flush = current_event.flush;
  hit.trigger_thr = current_event.trigger_thr;
  hit.event_size = current_event.event_size;
  hit.agg_counter = current_event.agg_counter;
  hit.block_id = current_event.block_id;
  hit.trace_len = current_event.trace_len;
  hit.ana_probe_type[0] = current_event.ana_probe_type[0];
  hit.ana_probe_type[1] = current_event.ana_probe_type[1];
  hit.dig_probe_type[0] = current_event.dig_probe_type[0];
  hit.dig_probe_type[1] = current_event.dig_probe_type[1];
  hit.dig_probe_type[2] = current_event.dig_probe_type[2];
  hit.dig_probe_type[3] = current_event.dig_probe_type[3];
  return hit;
}
