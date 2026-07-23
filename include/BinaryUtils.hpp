#ifndef BINARYUTILS_H
#define BINARYUTILS_H

#include <TArrayS.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <bitset>
#include <fstream>
#include <iostream>
#include <vector>

struct RawHit {
  UShort_t board;
  UShort_t channel;
  UShort_t energy;
  ULong64_t timestamp;
  UInt_t flags;
};

class BinaryReader {
protected:
  std::ifstream file;
  TString filename;
  Bool_t file_open;
  Long64_t bytes_read;

public:
  BinaryReader() : file_open(kFALSE), bytes_read(0) {}

  virtual ~BinaryReader() {
    if (file.is_open()) {
      file.close();
    }
  }

  virtual Bool_t Open(const char *fname) {
    filename = fname;
    file.open(fname, std::ios::binary);

    if (!file.is_open()) {
      std::cerr << "ERROR: Cannot open file " << fname << std::endl;
      file_open = kFALSE;
      return kFALSE;
    }

    file_open = kTRUE;
    bytes_read = 0;
    return kTRUE;
  }

  virtual void Close() {
    if (file.is_open()) {
      file.close();
    }
    file_open = kFALSE;
  }

  Bool_t IsOpen() const { return file_open; }
  TString GetFileName() const { return filename; }
  Long64_t GetBytesRead() const { return bytes_read; }

  Bool_t IsEOF() const { return !file_open || file.eof(); }

  virtual Bool_t ReadEvent() = 0;
};

class CoMPASSData {
public:
  enum WaveformCode : UChar_t {
    INPUT = 1,
    RC_CR = 2,
    RC_CR2 = 3,
    TRAPEZOID = 4,
    BASELINE = 5,
    THRESHOLD = 6,
    CFD = 7,
    TRAPEZOID_BASELINE = 8,
    FAST_TRIANGLE = 33,
    SMOOTHED_INPUT = 41
  };

  static const UInt_t DEADTIME_OCCURRED = 0x1;
  static const UInt_t TIMESTAMP_ROLLOVER = 0x2;
  static const UInt_t TIMESTAMP_RESET_EXT = 0x4;
  static const UInt_t FAKE_EVENT = 0x8;
  static const UInt_t MEMORY_FULL = 0x10;
  static const UInt_t TRIGGER_LOST = 0x20;
  static const UInt_t N_TRIGGERS_LOST = 0x40;
  static const UInt_t SATURATION_IN_GATE = 0x80;
  static const UInt_t TRIGGERS_1024_COUNTED = 0x100;
  static const UInt_t FIRST_AFTER_BUSY = 0x200;
  static const UInt_t INPUT_SATURATING = 0x400;
  static const UInt_t N_TRIGGERS_COUNTED = 0x800;
  static const UInt_t NOT_MATCHED_TIMEFILTER = 0x1000;
  static const UInt_t FINE_TIMESTAMP = 0x4000;
  static const UInt_t PILEUP = 0x8000;
  static const UInt_t PLL_LOCK_LOSS = 0x80000;
  static const UInt_t OVER_TEMPERATURE = 0x100000;
  static const UInt_t ADC_SHUTDOWN = 0x200000;

  UShort_t header;
  UShort_t board;
  UShort_t channel;
  ULong64_t timestamp;
  UShort_t energy_ch;
  Double_t energy_cal;
  UShort_t energy_short_ch;
  UInt_t flags;
  UChar_t waveform_code;
  UInt_t num_samples;
  TArrayS samples;

  Bool_t hasEnergyCh() const { return (header & 0x0001); }
  Bool_t hasEnergyCal() const { return (header & 0x0002); }
  Bool_t hasEnergyShort() const { return (header & 0x0004); }
  Bool_t hasWaveform() const { return (header & 0x0008); }
  UChar_t getControlBits() const { return header & 0x000F; }

  Bool_t hasDeadtime() const { return flags & DEADTIME_OCCURRED; }
  Bool_t hasTimestampRollover() const { return flags & TIMESTAMP_ROLLOVER; }
  Bool_t hasTimestampResetExt() const { return flags & TIMESTAMP_RESET_EXT; }
  Bool_t isFakeEvent() const { return flags & FAKE_EVENT; }
  Bool_t hasMemoryFull() const { return flags & MEMORY_FULL; }
  Bool_t hasTriggerLost() const { return flags & TRIGGER_LOST; }
  Bool_t hasNTriggersLost() const { return flags & N_TRIGGERS_LOST; }
  Bool_t hasSaturation() const { return flags & SATURATION_IN_GATE; }
  Bool_t has1024Triggers() const { return flags & TRIGGERS_1024_COUNTED; }
  Bool_t isFirstAfterBusy() const { return flags & FIRST_AFTER_BUSY; }
  Bool_t isInputSaturating() const { return flags & INPUT_SATURATING; }
  Bool_t hasNTriggersCounted() const { return flags & N_TRIGGERS_COUNTED; }
  Bool_t isNotMatchedTimeFilter() const {
    return flags & NOT_MATCHED_TIMEFILTER;
  }
  Bool_t hasFineTimestamp() const { return flags & FINE_TIMESTAMP; }
  Bool_t isPileup() const { return flags & PILEUP; }
  Bool_t hasPLLLockLoss() const { return flags & PLL_LOCK_LOSS; }
  Bool_t isOverTemperature() const { return flags & OVER_TEMPERATURE; }
  Bool_t isADCShutdown() const { return flags & ADC_SHUTDOWN; }

  TString getWaveformCodeName() const;
  std::vector<TString> getActiveFlags() const;
  void Print() const;
  void PrintHeader() const;
  void PrintFlags() const;
  void PrintWaveform() const;

  CoMPASSData();
  virtual ~CoMPASSData() {}
};

class CoMPASSReader : public BinaryReader {
private:
  UShort_t global_header;
  CoMPASSData current_event;

public:
  CoMPASSReader() : BinaryReader(), global_header(0) {}

  virtual ~CoMPASSReader() {}

  Bool_t Open(const char *fname) override;
  Bool_t Open(const char *fname, UShort_t header_override);

  Bool_t ReadEvent() override;

  const CoMPASSData &GetCurrentEvent() const { return current_event; }
  CoMPASSData &GetCurrentEvent() { return current_event; }

  UShort_t GetGlobalHeader() const { return global_header; }
};

class WaveDump742Data {
public:
  UInt_t event_size;
  UInt_t board_id;
  UInt_t pattern;
  UInt_t channel;
  UInt_t event_counter;
  UInt_t group_trigger_time_tag;
  UInt_t dc_offset;
  UInt_t start_index_cell;
  TArrayS samples;

  void Print() const;

  WaveDump742Data();
  virtual ~WaveDump742Data() {}
};

class WaveDump742Reader : public BinaryReader {
private:
  WaveDump742Data current_event;
  Bool_t corrections_enabled;

public:
  WaveDump742Reader(Bool_t with_corrections = kTRUE)
      : BinaryReader(), corrections_enabled(with_corrections) {}

  virtual ~WaveDump742Reader() {}

  Bool_t ReadEvent() override;

  const WaveDump742Data &GetCurrentEvent() const { return current_event; }
  WaveDump742Data &GetCurrentEvent() { return current_event; }

  void SetCorrectionsEnabled(Bool_t enable) { corrections_enabled = enable; }
  Bool_t GetCorrectionsEnabled() const { return corrections_enabled; }
};

// Lightweight per-block hit struct for SOL format (no variable-length traces)
struct SOLHit {
  UChar_t channel;
  UShort_t energy;
  UShort_t energy_short;
  ULong64_t timestamp;
  UShort_t fine_timestamp;
  UShort_t flags_high;
  UShort_t flags_low;
  UChar_t data_type;
  Bool_t is_psd;
  UChar_t down_sampling;
  UChar_t board_fail;
  UChar_t flush;
  UShort_t trigger_thr;
  ULong64_t event_size;
  UInt_t agg_counter;
  ULong64_t block_id;
  ULong64_t trace_len;
  UChar_t ana_probe_type[2];
  UChar_t dig_probe_type[4];
};

// Per-block data for SOL format with optional trace arrays
class SOLData {
public:
  enum DataType : UChar_t {
    ALL = 0,
    OneTrace = 1,
    NoTrace = 2,
    Minimum = 3,
    MiniWithFineTime = 4,
    Raw = 10
  };

  UShort_t block_header;
  UChar_t channel;
  UShort_t energy;
  UShort_t energy_short;
  ULong64_t timestamp;
  UShort_t fine_timestamp;
  UShort_t flags_high;
  UShort_t flags_low;
  UChar_t data_type;
  Bool_t is_psd;
  UChar_t down_sampling;
  UChar_t board_fail;
  UChar_t flush;
  UShort_t trigger_thr;
  ULong64_t event_size;
  UInt_t agg_counter;
  ULong64_t trace_len;
  ULong64_t block_id;
  UChar_t ana_probe_type[2];
  UChar_t dig_probe_type[4];

  // Interleaved trace buffer: [trace0[0], trace1[0], dig0[0], dig1[0],
  // dig2[0], dig3[0], trace0[1], ...] = 12 bytes per sample for ALL format.
  // For OneTrace format: [trace0[0], trace0[1], ...] = 4 bytes per sample.
  std::vector<Char_t> trace_data;

  Bool_t hasTraces() const {
    return (data_type == ALL || data_type == OneTrace) && trace_len > 0;
  }

  UInt_t getSamples() const { return static_cast<UInt_t>(trace_len); }

  const Int_t *getAnalog0() const {
    if (data_type == ALL &&
        trace_data.size() >= static_cast<std::size_t>(trace_len) * 12) {
      return reinterpret_cast<const Int_t *>(trace_data.data());
    }
    return nullptr;
  }

  const Int_t *getAnalog1() const {
    if (data_type == ALL &&
        trace_data.size() >= static_cast<std::size_t>(trace_len) * 12) {
      return reinterpret_cast<const Int_t *>(trace_data.data() + 4);
    }
    return nullptr;
  }

  const UChar_t *getDigital(UInt_t ch) const {
    if (ch > 3 || data_type != ALL ||
        trace_data.size() < static_cast<std::size_t>(trace_len) * 12) {
      return nullptr;
    }
    return reinterpret_cast<const UChar_t *>(trace_data.data() + 8 + ch);
  }

  const Int_t *getOneTrace() const {
    if (data_type == OneTrace &&
        trace_data.size() >=
            static_cast<std::size_t>(trace_len) * sizeof(Int_t)) {
      return reinterpret_cast<const Int_t *>(trace_data.data());
    }
    return nullptr;
  }

  void clearTraces() { trace_data.clear(); }

  TString getDataTypeName() const;
  void Print() const;

  SOLData();
  ~SOLData() {}
};

class SOLReader : public BinaryReader {
private:
  SOLData current_event;
  Long64_t block_id;
  Bool_t skip_traces;

public:
  SOLReader() : BinaryReader(), block_id(0), skip_traces(kFALSE) {}

  virtual ~SOLReader() {}

  void SetSkipTraces(Bool_t skip) { skip_traces = skip; }
  Bool_t GetSkipTraces() const { return skip_traces; }

  Bool_t ReadEvent() override;

  const SOLData &GetCurrentEvent() const { return current_event; }
  SOLData &GetCurrentEvent() { return current_event; }

  Long64_t GetBlockID() const { return block_id; }

  // Convert current SOLData event to a SOLHit struct (strips traces)
  SOLHit ToHit() const;

  // Split a SOLARIS file into time-bounded chunks. Each chunk file contains
  // blocks whose timestamps fall within chunkSeconds of the chunk start.
  // Output files are named <basename>_chunk<NNN>.sol. Returns the list of
  // created files. totalBlocks and totalChunks are filled with summary stats.
  static std::vector<TString> SplitSolFileByTime(const char *inputFile,
                                                 const char *outputDir,
                                                 Double_t chunkSeconds,
                                                 Int_t &totalBlocks,
                                                 Int_t &totalChunks);
};

#endif
