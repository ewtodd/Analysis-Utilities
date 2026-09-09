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

/**
 * @file BinaryUtils.hpp
 * @brief Readers for CAEN and SOLARIS acquisition binary formats.
 *
 * Three readers share the BinaryReader base: CoMPASSReader for CAEN CoMPASS
 * `.bin`, WaveDump742Reader for CAEN WaveDump DT5742 output, and SOLReader for
 * SOLARIS DAQ `.sol`. Each owns one "current event" object that ReadEvent()
 * overwrites in place, so a reader is a cursor over the file rather than a
 * container.
 *
 * @warning The current-event object and any pointer obtained from it are
 *          invalidated by the next ReadEvent() call. Copy anything that must
 *          outlive the iteration — SOLReader::ToHit() exists for exactly this.
 *
 * @note This header declares the reader classes directly and defines no
 *       `BinaryUtils` symbol of its own. From PyROOT they are
 * `ROOT.CoMPASSReader`, `ROOT.WaveDump742Reader` and `ROOT.SOLReader`.
 */

/**
 * @brief Minimal CoMPASS hit: header fields only, no waveform.
 *
 * What InitUtils::ConvertCoMPASSBinToHits() returns, for callers that want
 * timing and energy over a whole run without paying for traces.
 */
struct RawHit {
  UShort_t board;      ///< Digitiser board id.
  UShort_t channel;    ///< Channel index on that board.
  UShort_t energy;     ///< Energy in ADC channel units.
  ULong64_t timestamp; ///< Acquisition timestamp, in picoseconds.
  UInt_t flags;        ///< CoMPASS status bits; see CoMPASSData's constants.
};

/**
 * @brief Base for the binary readers: file handling plus a ReadEvent()
 * contract.
 *
 * Derived classes implement ReadEvent() to decode one record into their own
 * current-event object.
 */
class BinaryReader {
protected:
  std::ifstream file;
  TString filename;
  Bool_t file_open;
  Long64_t bytes_read;

public:
  /// @brief Construct a reader with no file open.
  BinaryReader() : file_open(kFALSE), bytes_read(0) {}

  /// @brief Closes the file if it is still open.
  virtual ~BinaryReader() {
    if (file.is_open()) {
      file.close();
    }
  }

  /**
   * @brief Open a binary file for reading.
   * @param fname Path to the file.
   * @return `kTRUE` on success; `kFALSE` if the file cannot be opened, with the
   *         reason printed to stderr.
   */
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

  /// @brief Close the file. Safe to call when nothing is open.
  virtual void Close() {
    if (file.is_open()) {
      file.close();
    }
    file_open = kFALSE;
  }

  /// @brief Whether a file is currently open.
  Bool_t IsOpen() const { return file_open; }
  /// @brief Path most recently passed to Open().
  TString GetFileName() const { return filename; }
  /// @brief Bytes consumed so far, for progress reporting.
  Long64_t GetBytesRead() const { return bytes_read; }

  /// @brief Whether the file is exhausted or was never opened.
  Bool_t IsEOF() const { return !file_open || file.eof(); }

  /**
   * @brief Decode the next record into the derived reader's current event.
   * @return `kTRUE` if a record was read; `kFALSE` at end of file or on a
   *         malformed record.
   */
  virtual Bool_t ReadEvent() = 0;
};

/**
 * @brief One decoded CoMPASS event, with header-bit and status-flag accessors.
 *
 * Which fields carry meaningful values is dictated by the file's global header,
 * mirrored into #header. Test it with hasEnergyCh() and friends before reading
 * the corresponding member.
 */
class CoMPASSData {
public:
  /// @brief Which CoMPASS processing stage a stored waveform came from.
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

  /// @name CoMPASS status flag bits, as found in #flags.
  /// @{
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
  /// @}

  UShort_t header;     ///< Global header word; says which fields are valid.
  UShort_t board;      ///< Digitiser board id.
  UShort_t channel;    ///< Channel index on that board.
  ULong64_t timestamp; ///< Acquisition timestamp, in picoseconds.
  UShort_t energy_ch;  ///< Energy in ADC channel units; valid if hasEnergyCh().
  Double_t energy_cal; ///< Calibrated energy; valid if hasEnergyCal().
  UShort_t energy_short_ch; ///< Short-gate energy; valid if hasEnergyShort().
  UInt_t flags;             ///< Status bits; see the flag constants above.
  UChar_t waveform_code;    ///< Processing stage of #samples; see WaveformCode.
  UInt_t num_samples;       ///< Length of #samples.
  TArrayS samples;          ///< Waveform samples; valid if hasWaveform().

  /// @brief Whether the record carries an uncalibrated energy (header bit 0).
  Bool_t hasEnergyCh() const { return (header & 0x0001); }
  /// @brief Whether the record carries a calibrated energy (header bit 1).
  Bool_t hasEnergyCal() const { return (header & 0x0002); }
  /// @brief Whether the record carries a short-gate energy (header bit 2).
  Bool_t hasEnergyShort() const { return (header & 0x0004); }
  /// @brief Whether the record carries waveform samples (header bit 3).
  Bool_t hasWaveform() const { return (header & 0x0008); }
  /// @brief The four header control bits as a nibble.
  UChar_t getControlBits() const { return header & 0x000F; }

  /// @brief Deadtime occurred before this event.
  Bool_t hasDeadtime() const { return flags & DEADTIME_OCCURRED; }
  /// @brief The timestamp counter rolled over.
  Bool_t hasTimestampRollover() const { return flags & TIMESTAMP_ROLLOVER; }
  /// @brief An external timestamp reset was applied.
  Bool_t hasTimestampResetExt() const { return flags & TIMESTAMP_RESET_EXT; }
  /// @brief Synthetic event inserted by the acquisition.
  Bool_t isFakeEvent() const { return flags & FAKE_EVENT; }
  /// @brief Board memory was full.
  Bool_t hasMemoryFull() const { return flags & MEMORY_FULL; }
  /// @brief At least one trigger was lost.
  Bool_t hasTriggerLost() const { return flags & TRIGGER_LOST; }
  /// @brief The lost-trigger count field is present.
  Bool_t hasNTriggersLost() const { return flags & N_TRIGGERS_LOST; }
  /// @brief The input saturated within the integration gate.
  Bool_t hasSaturation() const { return flags & SATURATION_IN_GATE; }
  /// @brief 1024 triggers have been counted.
  Bool_t has1024Triggers() const { return flags & TRIGGERS_1024_COUNTED; }
  /// @brief First event after a busy period.
  Bool_t isFirstAfterBusy() const { return flags & FIRST_AFTER_BUSY; }
  /// @brief The input is currently saturating.
  Bool_t isInputSaturating() const { return flags & INPUT_SATURATING; }
  /// @brief The counted-trigger field is present.
  Bool_t hasNTriggersCounted() const { return flags & N_TRIGGERS_COUNTED; }
  /// @brief The event failed the coincidence time filter.
  Bool_t isNotMatchedTimeFilter() const {
    return flags & NOT_MATCHED_TIMEFILTER;
  }
  /// @brief A fine timestamp is available.
  Bool_t hasFineTimestamp() const { return flags & FINE_TIMESTAMP; }
  /// @brief Pileup was detected.
  Bool_t isPileup() const { return flags & PILEUP; }
  /// @brief The PLL lost lock.
  Bool_t hasPLLLockLoss() const { return flags & PLL_LOCK_LOSS; }
  /// @brief The board reported over-temperature.
  Bool_t isOverTemperature() const { return flags & OVER_TEMPERATURE; }
  /// @brief The ADC shut down.
  Bool_t isADCShutdown() const { return flags & ADC_SHUTDOWN; }

  /// @brief Human-readable name for #waveform_code.
  /// @return The WaveformCode name, or a placeholder if unrecognised.
  TString getWaveformCodeName() const;
  /// @brief Names of every status flag currently set in #flags.
  std::vector<TString> getActiveFlags() const;
  /// @brief Print the whole event to stdout.
  /// @brief Print the whole event to stdout.
  /// @brief Print the block to stdout.
  void Print() const;
  /// @brief Print the header fields to stdout.
  void PrintHeader() const;
  /// @brief Print the active status flags to stdout.
  void PrintFlags() const;
  /// @brief Print the waveform samples to stdout.
  void PrintWaveform() const;

  /// @brief Construct with all fields zeroed.
  CoMPASSData();
  virtual ~CoMPASSData() {}
};

/**
 * @brief Cursor over a CoMPASS binary file.
 *
 * Reads the two-byte global header on open, then decodes one event per
 * ReadEvent() into a reused CoMPASSData.
 */
class CoMPASSReader : public BinaryReader {
private:
  UShort_t global_header;
  CoMPASSData current_event;

public:
  /// @brief Construct a reader with no file open.
  CoMPASSReader() : BinaryReader(), global_header(0) {}

  virtual ~CoMPASSReader() {}

  /**
   * @brief Open a file and read its global header.
   * @param fname Path to the CoMPASS binary file.
   * @return `kTRUE` on success.
   */
  Bool_t Open(const char *fname) override;
  /**
   * @brief Open a file, assuming a given global header.
   *
   * For files whose two-byte header is missing or wrong. Passing an override
   * also suppresses the two-byte header skip, since there is nothing to skip.
   *
   * @param fname           Path to the CoMPASS binary file.
   * @param header_override Header word to assume. Must be non-zero.
   * @return `kTRUE` on success.
   */
  Bool_t Open(const char *fname, UShort_t header_override);

  /**
   * @brief Decode the next event into the current-event object.
   * @return `kTRUE` if an event was read; `kFALSE` at end of file or on a
   *         malformed record.
   */
  /**
   * @brief Decode the next event into the current-event object.
   * @return `kTRUE` if an event was read; `kFALSE` at end of file or on a
   *         malformed record.
   */
  /**
   * @brief Decode the next block into the current-event object.
   * @return `kTRUE` if a block was read; `kFALSE` at end of file.
   *
   * @note As a memory-safety measure, a block claiming more than 100 000
   *       samples is treated as corrupt and its trace payload is skipped.
   */
  Bool_t ReadEvent() override;

  /// @brief The event most recently read.
  /// @return A reference valid only until the next ReadEvent().
  const CoMPASSData &GetCurrentEvent() const { return current_event; }
  /// @brief Mutable access to the event most recently read.
  /// @return A reference valid only until the next ReadEvent().
  CoMPASSData &GetCurrentEvent() { return current_event; }

  /// @brief The global header in force, read from the file or overridden.
  UShort_t GetGlobalHeader() const { return global_header; }
};

/**
 * @brief One decoded WaveDump event from a DT5742-family digitiser.
 */
class WaveDump742Data {
public:
  UInt_t event_size;             ///< Record size in 32-bit words.
  UInt_t board_id;               ///< Digitiser board id.
  UInt_t pattern;                ///< Front-panel I/O pattern word.
  UInt_t channel;                ///< Channel index.
  UInt_t event_counter;          ///< Board's running event counter.
  UInt_t group_trigger_time_tag; ///< Group trigger time tag.
  UInt_t dc_offset;              ///< DC offset applied to this channel.
  /// Switched-capacitor start cell. Needed to unwrap the circular sampling
  /// array, and the basis of the timing corrections.
  UInt_t start_index_cell;
  TArrayS samples; ///< Waveform samples, in ADC counts.

  void Print() const;

  /// @brief Construct with all fields zeroed.
  WaveDump742Data();
  virtual ~WaveDump742Data() {}
};

/**
 * @brief Cursor over a WaveDump DT5742 binary file.
 */
class WaveDump742Reader : public BinaryReader {
private:
  WaveDump742Data current_event;
  Bool_t corrections_enabled;

public:
  /**
   * @brief Construct a reader with no file open.
   * @param with_corrections `kTRUE` (default) to apply the DT5742 timing
   *                         corrections as events are decoded.
   */
  WaveDump742Reader(Bool_t with_corrections = kTRUE)
      : BinaryReader(), corrections_enabled(with_corrections) {}

  virtual ~WaveDump742Reader() {}

  Bool_t ReadEvent() override;

  /// @brief The event most recently read.
  /// @return A reference valid only until the next ReadEvent().
  const WaveDump742Data &GetCurrentEvent() const { return current_event; }
  /// @brief Mutable access to the event most recently read.
  /// @return A reference valid only until the next ReadEvent().
  WaveDump742Data &GetCurrentEvent() { return current_event; }

  /// @brief Turn the DT5742 timing corrections on or off.
  /// @param enable `kTRUE` to apply them to subsequently read events.
  void SetCorrectionsEnabled(Bool_t enable) { corrections_enabled = enable; }
  /// @brief Whether the timing corrections are being applied.
  Bool_t GetCorrectionsEnabled() const { return corrections_enabled; }
};

/**
 * @brief Lightweight SOL block: header fields only, traces stripped.
 *
 * What SOLReader::ToHit() produces and what
 * InitUtils::ConvertSOLBinToHits() returns. Fixed size, so a whole run's worth
 * fits in a vector without the memory cost of the traces.
 */
struct SOLHit {
  UChar_t channel;           ///< Channel index.
  UShort_t energy;           ///< Long-gate energy, in ADC channel units.
  UShort_t energy_short;     ///< Short-gate energy; meaningful when #is_psd.
  ULong64_t timestamp;       ///< Coarse timestamp, in nanoseconds.
  UShort_t fine_timestamp;   ///< Sub-nanosecond refinement of #timestamp.
  UShort_t flags_high;       ///< Upper status word.
  UShort_t flags_low;        ///< Lower status word.
  UChar_t data_type;         ///< Block format; see SOLData::DataType.
  Bool_t is_psd;             ///< Whether this block came from a PSD firmware.
  UChar_t down_sampling;     ///< Trace downsampling factor applied on-board.
  UChar_t board_fail;        ///< Board failure indicator.
  UChar_t flush;             ///< Flush indicator.
  UShort_t trigger_thr;      ///< Trigger threshold in force.
  ULong64_t event_size;      ///< Size of the source block, in bytes.
  UInt_t agg_counter;        ///< Aggregate counter from the board.
  ULong64_t block_id;        ///< Zero-based index of this block in the file.
  ULong64_t trace_len;       ///< Samples the source block carried, though this
                             ///< struct stores none of them.
  UChar_t ana_probe_type[2]; ///< Analogue probe selection for traces 0 and 1.
  UChar_t dig_probe_type[4]; ///< Digital probe selection for traces 0 to 3.
};

/**
 * @brief One decoded SOL block, optionally carrying its traces.
 *
 * @warning The trace accessors hand back raw pointers into #trace_data, which
 *          the reader reuses. Every one of them is invalidated by
 *          clearTraces() and by the next SOLReader::ReadEvent(). Copy what you
 *          need before advancing the reader, or use SOLReader::ToHit() for the
 *          header fields.
 */
class SOLData {
public:
  /**
   * @brief SOL block layout, which decides both the fields and the traces
   *        present.
   *
   * Only #ALL and #OneTrace carry traces; the rest are header-only.
   */
  enum DataType : UChar_t {
    ALL = 0,
    OneTrace = 1,
    NoTrace = 2,
    Minimum = 3,
    MiniWithFineTime = 4,
    Raw = 10
  };

  /// Raw block header word.
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

  // Trace buffer, stored exactly as the SOLARIS DAQ writes it (SolReader.h in
  // SOLARIS_DAQ): contiguous blocks, NOT interleaved samples. ALL format:
  //   [trace0[0..len-1] as Int32][trace1[0..len-1] as Int32]
  //   [dig0[0..len-1]][dig1[0..len-1]][dig2[0..len-1]][dig3[0..len-1]]
  // i.e. 12 bytes per sample in total. OneTrace format: [trace0[0..len-1]] as
  // Int32, 4 bytes per sample. Each accessor returns a pointer to its block,
  // indexed by sample.
  /// Raw trace bytes, exactly as written by the DAQ. Prefer the typed
  /// accessors below to indexing this directly.
  std::vector<Char_t> trace_data;

  /// @brief Whether this block carries traces.
  /// @return `kTRUE` for #ALL and #OneTrace blocks with a non-zero length.
  Bool_t hasTraces() const {
    return (data_type == ALL || data_type == OneTrace) && trace_len > 0;
  }

  /// @brief Samples per trace in this block.
  UInt_t getSamples() const { return static_cast<UInt_t>(trace_len); }

  /**
   * @brief Pointer to the first analogue trace.
   * @return Pointer to @ref getSamples() samples, or null if this block format
   *         has no such trace or the buffer is short.
   *
   * @warning Invalidated by clearTraces() and by the next
   *          SOLReader::ReadEvent().
   */
  const Int_t *getAnalog0() const {
    const std::size_t len = static_cast<std::size_t>(trace_len);
    if ((data_type == ALL && trace_data.size() >= len * 12) ||
        (data_type == OneTrace && trace_data.size() >= len * 4)) {
      return reinterpret_cast<const Int_t *>(trace_data.data());
    }
    return nullptr;
  }

  /**
   * @brief Pointer to the second analogue trace.
   * @return Pointer to @ref getSamples() samples, or null for anything but an
   *         #ALL block.
   *
   * @warning Invalidated by clearTraces() and by the next
   *          SOLReader::ReadEvent().
   */
  const Int_t *getAnalog1() const {
    const std::size_t len = static_cast<std::size_t>(trace_len);
    if (data_type == ALL && trace_data.size() >= len * 12) {
      return reinterpret_cast<const Int_t *>(trace_data.data() + len * 4);
    }
    return nullptr;
  }

  /**
   * @brief Pointer to one digital probe trace.
   * @param ch Probe index, 0 to 3.
   * @return Pointer to @ref getSamples() samples, or null if @p ch is out of
   *         range, this is not an #ALL block, or the buffer is short.
   *
   * @warning Invalidated by clearTraces() and by the next
   *          SOLReader::ReadEvent().
   */
  const UChar_t *getDigital(UInt_t ch) const {
    const std::size_t len = static_cast<std::size_t>(trace_len);
    if (ch > 3 || data_type != ALL || trace_data.size() < len * 12) {
      return nullptr;
    }
    return reinterpret_cast<const UChar_t *>(trace_data.data() + len * 8 +
                                             len * ch);
  }

  /**
   * @brief Pointer to the single trace of a #OneTrace block.
   * @return Pointer to @ref getSamples() samples, or null for any other format.
   *
   * @warning Invalidated by clearTraces() and by the next
   *          SOLReader::ReadEvent().
   */
  const Int_t *getOneTrace() const {
    if (data_type == OneTrace &&
        trace_data.size() >=
            static_cast<std::size_t>(trace_len) * sizeof(Int_t)) {
      return reinterpret_cast<const Int_t *>(trace_data.data());
    }
    return nullptr;
  }

  /// @brief Release the trace buffer, invalidating every trace pointer.
  void clearTraces() { trace_data.clear(); }

  /// @brief Human-readable name for #data_type.
  TString getDataTypeName() const;
  void Print() const;

  /// @brief Construct with all fields zeroed and no traces.
  SOLData();
  ~SOLData() {}
};

/**
 * @brief Cursor over a SOLARIS DAQ `.sol` file.
 */
class SOLReader : public BinaryReader {
private:
  SOLData current_event;
  Long64_t block_id;
  Bool_t skip_traces;

public:
  /// @brief Construct a reader with no file open and traces enabled.
  SOLReader() : BinaryReader(), block_id(0), skip_traces(kFALSE) {}

  virtual ~SOLReader() {}

  /**
   * @brief Skip trace payloads while still parsing every block.
   * @param skip `kTRUE` to discard traces as blocks are read. Much faster and
   *             far cheaper in memory when only header fields are needed.
   */
  void SetSkipTraces(Bool_t skip) { skip_traces = skip; }
  /// @brief Whether trace payloads are being skipped.
  Bool_t GetSkipTraces() const { return skip_traces; }

  Bool_t ReadEvent() override;

  /// @brief The block most recently read.
  /// @return A reference valid only until the next ReadEvent().
  const SOLData &GetCurrentEvent() const { return current_event; }
  /// @brief Mutable access to the block most recently read.
  /// @return A reference valid only until the next ReadEvent().
  SOLData &GetCurrentEvent() { return current_event; }

  /// @brief Zero-based index of the block most recently read.
  Long64_t GetBlockID() const { return block_id; }

  /**
   * @brief Copy the current block's header fields into a standalone SOLHit.
   *
   * The way to retain a block past the next ReadEvent(): the result owns its
   * data and carries no traces.
   *
   * @return A SOLHit describing the current block.
   */
  SOLHit ToHit() const;

  /**
   * @brief Split a run file into time-bounded chunks.
   *
   * Each output file collects blocks whose timestamps fall within
   * @p chunkSeconds of that chunk's first block. Reads and writes are
   * bulk-buffered at 4 MiB, so splitting runs close to disk speed.
   *
   * @param inputFile   Path to the source `.sol` file.
   * @param outputDir   Directory for the chunk files, which are named
   *                    `<basename>_chunk<NNN>.sol`.
   * @param chunkSeconds Chunk duration in seconds. Compared against the
   *                    nanosecond block timestamps internally.
   * @param totalBlocks Set to the number of blocks processed.
   * @param totalChunks Set to the number of chunk files created.
   *
   * @return Paths of the chunk files created, in order.
   */
  static std::vector<TString> SplitSolFileByTime(const char *inputFile,
                                                 const char *outputDir,
                                                 Double_t chunkSeconds,
                                                 Int_t &totalBlocks,
                                                 Int_t &totalChunks);
};

#endif
