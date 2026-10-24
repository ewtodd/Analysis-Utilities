#ifndef WAVEFORMPROCESSINGUTILS_H
#define WAVEFORMPROCESSINGUTILS_H

/**
 * @file WaveformProcessingUtils.hpp
 * @brief Waveform sources and processors: digitiser records in, one tree of
 *        per-waveform features out.
 *
 * Two halves, each a small class family in the manner of BinaryReader:
 *
 * - A WaveformSource reads records from one input format and hands them over
 *   one at a time as WaveformRecord. There is one per ROOT layout (CoMPASS,
 *   WaveDump, SOLARIS) and per binary format (CoMPASSReader, WaveDump742Reader,
 *   SOLReader); CreateWaveformSource() picks one from an InputFormat, and a
 *   caller with another format derives its own.
 * - A WaveformProcessor runs every record through the shared steps (polarity,
 *   baseline, trigger, crop, clipping cut) and leaves the features to its
 *   variant: PulseIntegralProcessor integrates short and long gates (pulse
 *   shape discrimination), PulseHeightProcessor reads a trapezoid filter on
 *   the pole-zero-deconvolved step of a charge-sensitive preamplifier.
 *
 * From PyROOT the classes are `ROOT.PulseIntegralProcessor`,
 * `ROOT.PulseHeightProcessor`, `ROOT.FileProcessingConfig` and so on.
 */

#include "BinaryUtils.hpp"
#include "PlottingUtils.hpp"
#include <TArrayF.h>
#include <TArrayI.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <mutex>
#include <vector>

// ------------------------------------------------------------------- sources

/// @brief The input layouts CreateWaveformSource() knows.
enum class InputFormat {
  /// ROOT tree `Data_R` as CoMPASS (or InitUtils::ConvertCoMPASSBinToROOT)
  /// writes it: `Samples`, `Timestamp`, `Energy`, `Flags`, `Channel`, `Board`.
  kCOMPASS,
  /// ROOT tree `Data_R` from a DT5742 WaveDump conversion: `Samples`,
  /// `TriggerTimeTag`, `Channel`.
  kWAVEDUMP,
  /// ROOT tree `Data_R` from a SOLARIS conversion: `Trace0`, `Timestamp`,
  /// `Energy`, `Channel`, `FlagsHigh`, `FlagsLow`.
  kSOLARIS,
  /// CoMPASS `.bin`, read with CoMPASSReader.
  kCOMPASS_BINARY,
  /// DT5742 WaveDump binary, read with WaveDump742Reader.
  kWAVEDUMP742_BINARY,
  /// SOLARIS `.sol`, read with SOLReader; the first analogue trace.
  kSOLARIS_BINARY
};

/**
 * @brief One digitiser record: a waveform and what the digitiser said about
 *        it.
 *
 * Fields a format does not carry keep their defaults.
 */
struct WaveformRecord {
  TArrayI samples;         ///< Waveform samples, in ADC counts.
  ULong64_t timestamp = 0; ///< In the source's own units (CoMPASS ps, SOLARIS
                           ///< ns, WaveDump trigger time tag ticks).
  Int_t board = -1;        ///< Digitiser board, or -1 when not recorded.
  Int_t channel = -1;      ///< Channel, or -1 when not recorded.
  /// The digitiser's own energy (CoMPASS's channels or calibrated energy,
  /// SOLARIS's long-gate energy), or -1 when not recorded.
  Double_t energy = -1;
  /// Status flags (CoMPASS's flag word; SOLARIS's high word above its low
  /// word), or 0 when not recorded.
  UInt_t flags = 0;
};

/**
 * @brief Base for the waveform sources: a cursor over one file's records.
 *
 * Derived classes read one format; WaveformProcessor::ProcessSource() takes
 * any of them, so a format this header does not know needs only a new
 * derived class.
 */
class WaveformSource {
public:
  /// @brief Construct with nothing open.
  WaveformSource() {}
  /// @brief Derived sources close what they opened.
  virtual ~WaveformSource() {}

  /**
   * @brief Open a file.
   * @param path Path to the file.
   * @return `kTRUE` on success; `kFALSE` if it cannot be opened or is not in
   *         this source's format, with the reason printed.
   */
  virtual Bool_t Open(const TString &path) = 0;

  /**
   * @brief Read the next record that carries a waveform.
   * @param record Overwritten with the record; its samples are copied, so it
   *               stays valid after the next call.
   * @return `kTRUE` if a record was read; `kFALSE` at the end of the file.
   */
  virtual Bool_t Next(WaveformRecord &record) = 0;

  /// @brief Close the file. Safe to call when nothing is open.
  virtual void Close() = 0;
};

/**
 * @brief Shared part of the ROOT-tree sources: the file, the `Data_R` tree
 *        and the entry cursor.
 *
 * Scalar branches are read through their leaves, whatever their stored type,
 * because the writers differ (CoMPASS stores `Channel` as a short, the
 * InitUtils conversion as an unsigned short).
 */
class RootTreeWaveformSource : public WaveformSource {
public:
  /// @brief Construct with nothing open.
  RootTreeWaveformSource();
  /// @brief Closes the file if it is still open.
  virtual ~RootTreeWaveformSource();

  /**
   * @brief Open a ROOT file and its `Data_R` tree, and bind the branches.
   * @param path Path to the file.
   * @return `kTRUE` on success; `kFALSE` if the file, the tree or a required
   *         branch is missing.
   */
  Bool_t Open(const TString &path) override;

  /// @brief Close the file. Safe to call when nothing is open.
  void Close() override;

protected:
  /// The open file, or null.
  TFile *file_;
  /// The `Data_R` tree, owned by #file_.
  TTree *tree_;
  /// The next entry to read.
  Long64_t entry_;

  /**
   * @brief Bind the derived source's waveform branch.
   * @return `kTRUE` if the tree carries it.
   */
  virtual Bool_t BindBranches() = 0;

  /**
   * @brief A scalar leaf's value for the current entry.
   * @param name     Leaf name.
   * @param fallback Returned when the tree has no such leaf.
   * @return The value, converted to double.
   */
  Double_t LeafValue(const char *name, Double_t fallback) const;
};

/// @brief CoMPASS's ROOT layout (InputFormat::kCOMPASS).
class CoMPASSRootSource : public RootTreeWaveformSource {
public:
  /// @brief Construct with nothing open.
  CoMPASSRootSource();
  /// @brief Releases the sample buffer.
  ~CoMPASSRootSource();
  /// @copydoc WaveformSource::Next
  Bool_t Next(WaveformRecord &record) override;

protected:
  /// @copydoc RootTreeWaveformSource::BindBranches
  Bool_t BindBranches() override;

private:
  TArrayS *samples_;
};

/// @brief The WaveDump ROOT layout (InputFormat::kWAVEDUMP).
class WaveDumpRootSource : public RootTreeWaveformSource {
public:
  /// @brief Construct with nothing open.
  WaveDumpRootSource();
  /// @brief Releases the sample buffer.
  ~WaveDumpRootSource();
  /// @copydoc WaveformSource::Next
  Bool_t Next(WaveformRecord &record) override;

protected:
  /// @copydoc RootTreeWaveformSource::BindBranches
  Bool_t BindBranches() override;

private:
  TArrayS *samples_;
};

/// @brief The SOLARIS ROOT layout (InputFormat::kSOLARIS).
class SolarisRootSource : public RootTreeWaveformSource {
public:
  /// @brief Construct with nothing open.
  SolarisRootSource();
  /// @brief Releases the trace buffer.
  ~SolarisRootSource();
  /// @copydoc WaveformSource::Next
  Bool_t Next(WaveformRecord &record) override;

protected:
  /// @copydoc RootTreeWaveformSource::BindBranches
  Bool_t BindBranches() override;

private:
  TArrayI *trace_;
};

/// @brief CoMPASS `.bin` files through CoMPASSReader
///        (InputFormat::kCOMPASS_BINARY). Records without a waveform are
///        skipped.
class CoMPASSBinarySource : public WaveformSource {
public:
  /// @copydoc WaveformSource::Open
  Bool_t Open(const TString &path) override;
  /// @copydoc WaveformSource::Next
  Bool_t Next(WaveformRecord &record) override;
  /// @copydoc WaveformSource::Close
  void Close() override;

private:
  CoMPASSReader reader_;
};

/// @brief DT5742 WaveDump binaries through WaveDump742Reader
///        (InputFormat::kWAVEDUMP742_BINARY).
class WaveDump742BinarySource : public WaveformSource {
public:
  /// @copydoc WaveformSource::Open
  Bool_t Open(const TString &path) override;
  /// @copydoc WaveformSource::Next
  Bool_t Next(WaveformRecord &record) override;
  /// @copydoc WaveformSource::Close
  void Close() override;

private:
  WaveDump742Reader reader_;
};

/// @brief SOLARIS `.sol` files through SOLReader
///        (InputFormat::kSOLARIS_BINARY): the first analogue trace of every
///        block that carries traces.
class SolarisBinarySource : public WaveformSource {
public:
  /// @copydoc WaveformSource::Open
  Bool_t Open(const TString &path) override;
  /// @copydoc WaveformSource::Next
  Bool_t Next(WaveformRecord &record) override;
  /// @copydoc WaveformSource::Close
  void Close() override;

private:
  SOLReader reader_;
};

/**
 * @brief A new source for an input format.
 * @param format The format to read.
 * @return A source with nothing open, owned by the caller.
 */
WaveformSource *CreateWaveformSource(InputFormat format);

// ---------------------------------------------------------------- processors

/**
 * @brief What every processor variant takes: input, trigger, crop and output.
 *
 * All sample counts are in samples. The variants' own parameters are in
 * PulseIntegralSettings and PulseHeightSettings.
 */
struct FileProcessingConfig {
  /// Pulse polarity: `-1` for negative-going pulses (inverted during baseline
  /// subtraction), `+1` for positive-going.
  Int_t polarity = -1;
  /// Trigger level as a fraction of the waveform's peak, in (0, 1].
  Float_t trigger_threshold = 0.15;
  /// Leading samples averaged to estimate the baseline.
  Int_t num_samples_baseline = 10;
  /// Samples kept before the trigger when cropping.
  Int_t pre_samples = 17;
  /// Samples kept after the trigger when cropping.
  Int_t post_samples = 190;
  /// How many example waveforms to write out as figures, for eyeballing.
  Int_t sample_waveforms_to_save = 5;
  /// Cap on accepted waveforms per file; `-1` means no limit.
  Int_t max_events = -1;
  Bool_t verbose = kTRUE;         ///< Print progress and per-file summaries.
  Bool_t store_waveforms = kTRUE; ///< Write the cropped waveform to the tree.
  InputFormat input_format = InputFormat::kCOMPASS; ///< Input file format.
  /// ADC code that means full scale. A raw extremum equal to this (positive
  /// polarity) or to zero (negative polarity) is treated as clipped.
  Int_t adc_saturation_code = 16384;
  /// Process only this channel; `-1` processes every channel, which suits a
  /// file per channel. Records of other channels are skipped uncounted.
  Int_t channel = -1;
};

/**
 * @brief Running counters and accumulators for one processing run.
 *
 * The `rejected_*` counters partition the rejected waveforms by cause, so
 * `total_processed - accepted` equals their sum. A variant only ever uses the
 * counters of its own cuts.
 */
struct ProcessingStats {
  Int_t total_processed = 0;     ///< Waveforms examined.
  Int_t accepted = 0;            ///< Waveforms that passed every cut.
  Int_t rejected_no_trigger = 0; ///< No sample reached the trigger level.
  /// Trigger too close to either end to crop the requested window, or (pulse
  /// height) the trapezoid does not fit inside the crop.
  Int_t rejected_insufficient_samples = 0;
  /// Long-gate integral <= 0 (pulse integral).
  Int_t rejected_negative_integral = 0;
  /// More than 50% negative samples in the long gate (pulse integral).
  Int_t rejected_baseline = 0;
  Int_t rejected_clipped = 0;      ///< Raw extremum hit the ADC rail.
  Double_t sum_baseline_rms = 0.0; ///< Sum of per-waveform baseline RMS.
  Int_t baseline_rms_count = 0;    ///< Waveforms contributing to that sum.
  /// Sum of baseline RMS over accepted waveforms only.
  Double_t sum_baseline_rms_accepted = 0.0;
  Int_t baseline_rms_count_accepted = 0; ///< Accepted waveforms in that sum.
};

/**
 * @brief Base of the processor variants: the shared pipeline from a record to
 *        a row of the output tree.
 *
 * For every record: the raw extremum (for the clipping cut), the baseline
 * from the leading samples (subtracted, and inverted for negative polarity),
 * the trigger at a fraction of the peak, the crop around it, the clipping
 * cut, then the variant's features and cuts (ExtractFeatures()). Every
 * variant's tree `features` carries `timestamp`, `channel`, `board`,
 * `digitizer_energy`, `flags`, `trigger_position`, `pulse_height` (the
 * cropped waveform's maximum, ADC counts above baseline), `baseline_rms`, and
 * with FileProcessingConfig::store_waveforms the cropped waveform `Samples`;
 * the variant adds its own branches.
 *
 * @note One instance is not safe to share across threads.
 *       ProcessFilesParallel() gives every worker its own copy (Clone()).
 */
class WaveformProcessor {
public:
  /**
   * @brief Construct from the shared configuration.
   * @param config Copied into the instance.
   */
  WaveformProcessor(const FileProcessingConfig &config);

  /// @brief Closes any open output file and releases the buffers.
  virtual ~WaveformProcessor();

  /**
   * @brief A new processor of the same variant and settings.
   * @return A copy with fresh statistics, owned by the caller.
   */
  virtual WaveformProcessor *Clone() const = 0;

  /**
   * @brief Process one input file end to end.
   *
   * Opens it with the source FileProcessingConfig::input_format names, then
   * as ProcessSource().
   *
   * @param filepath    Path to the input file.
   * @param output_name Output basename, without extension, resolved against
   *                    the ROOT files base directory.
   * @return `kTRUE` on success; `kFALSE` if the input could not be opened or
   *         the output could not be written.
   */
  Bool_t ProcessFile(const TString &filepath, const TString &output_name);

  /**
   * @brief Process every record of an already-open source.
   *
   * Creates the output tree, runs every record of the configured channel
   * through the pipeline, then writes the tree, the sample figures and the
   * `.stats` summary. The source is left open.
   *
   * @param source      An open source; read to its end.
   * @param output_name Output basename, without extension, resolved against
   *                    the ROOT files base directory.
   * @return `kTRUE` on success; `kFALSE` if the output could not be written
   *         or the variant's settings do not fit the crop.
   */
  Bool_t ProcessSource(WaveformSource &source, const TString &output_name);

  /**
   * @brief Process many files concurrently, a copy of this processor per
   *        worker.
   *
   * Enables ROOT thread safety and IO::SetThreadSafe(), then dispatches the
   * files in batches of @p max_workers, waiting for each batch to finish
   * before starting the next.
   *
   * @param filepaths    Input files.
   * @param output_names Output basenames, positionally matched to
   *                     @p filepaths. Must be the same length.
   * @param max_workers  Concurrent workers. A value of `0` or less means
   *                     hardware concurrency. Capped at the number of files.
   *
   * @note Batches are synchronous, so one slow file holds up the start of the
   *       next batch. Per-file success or failure is printed as each
   *       completes; there is no aggregate return value.
   *
   * @warning Enables ROOT thread safety process-wide as a side effect, which
   *          cannot be undone.
   */
  void ProcessFilesParallel(const std::vector<TString> &filepaths,
                            const std::vector<TString> &output_names,
                            Int_t max_workers = 4) const;

  /**
   * @brief Print the run statistics, including mean baseline RMS.
   *
   * Reports the counters plus two means: baseline RMS over all processed
   * waveforms, and over accepted ones only. The same summary is written to
   * the `.stats` file alongside the output.
   */
  void PrintAllStatistics() const;

  /// @brief Snapshot of the current statistics.
  /// @return A copy of the counters accumulated so far.
  ProcessingStats GetStats() const { return stats_; }

  /// @brief The shared configuration.
  /// @return The configuration this processor was built with.
  const FileProcessingConfig &GetConfig() const { return config_; }

protected:
  /// The shared configuration.
  FileProcessingConfig config_;
  /// Counters of the current run.
  ProcessingStats stats_;
  /// Output basename of the current run, for figure and `.stats` names.
  TString current_output_name_;
  /// Sample figures written so far in the current run.
  Int_t sample_waveforms_saved_;
  /// Serialises figure drawing, which ROOT does not allow concurrently.
  static std::mutex canvas_mutex_;

  /**
   * @brief Add the variant's branches to the output tree.
   * @param tree The output tree; the branches point into the variant's own
   *             members.
   */
  virtual void BookBranches(TTree *tree) = 0;

  /**
   * @brief Whether the variant's settings fit the configured crop.
   * @return `kTRUE` if they do; otherwise the reason is printed and the run
   *         does not start.
   */
  virtual Bool_t CheckSettings() const { return kTRUE; }

  /**
   * @brief Compute the variant's features and apply its cuts.
   * @param cropped Cropped, baseline-subtracted waveform with a positive
   *                pulse, the trigger at FileProcessingConfig::pre_samples.
   * @return `kTRUE` if the waveform passes; on failure the matching
   *         ProcessingStats counter has been incremented.
   */
  virtual Bool_t ExtractFeatures(const TArrayF &cropped) = 0;

  /**
   * @brief Write one accepted waveform out as a figure, for inspection.
   * @param cropped The cropped, baseline-subtracted waveform.
   * @note Callers hold #canvas_mutex_.
   */
  virtual void DrawSampleWaveform(const TArrayF &cropped);

private:
  static Bool_t ProcessWithClone(const WaveformProcessor *prototype,
                                 const TString &filepath,
                                 const TString &output_name);
  Bool_t ProcessRecord(const WaveformRecord &record);
  void SubtractBaseline(const TArrayI &samples);
  Float_t FindTrigger(const TArrayF &waveform) const;
  void CropWaveform(Int_t trigger_pos);
  Bool_t OpenOutput(const TString &output_name);
  void CloseOutput();

  TFile *output_file_;
  TTree *output_tree_;
  TArrayF *baseline_subtracted_;
  TArrayF *cropped_;
  Float_t current_baseline_rms_;
  Bool_t current_baseline_rms_valid_;

  ULong64_t out_timestamp_;
  Int_t out_channel_;
  Int_t out_board_;
  Double_t out_digitizer_energy_;
  UInt_t out_flags_;
  Int_t out_trigger_position_;
  Float_t out_pulse_height_;
  Float_t out_baseline_rms_;
};

/// @brief The integration gates of PulseIntegralProcessor, in samples.
struct PulseIntegralSettings {
  /// Samples before FileProcessingConfig::pre_samples at which integration
  /// starts, so integration begins at `pre_samples - pre_gate` within the
  /// cropped waveform.
  Int_t pre_gate = 10;
  /// Width of the short integration gate.
  Int_t short_gate = 10;
  /// Width of the long integration gate. Should exceed #short_gate; the PSD
  /// ratio is meaningless otherwise.
  Int_t long_gate = 200;
};

/**
 * @brief Short- and long-gate integrals, for pulse shape discrimination.
 *
 * Adds `short_integral` and `long_integral` (ADC counts x samples) to the
 * tree. Rejects waveforms with more than 50% of the long gate below zero (a
 * bad baseline) and a non-positive long integral.
 */
class PulseIntegralProcessor : public WaveformProcessor {
public:
  /**
   * @brief Construct from the shared configuration and the gates.
   * @param config   Shared configuration, copied.
   * @param settings Integration gates, copied.
   */
  PulseIntegralProcessor(const FileProcessingConfig &config,
                         const PulseIntegralSettings &settings);

  /// @copydoc WaveformProcessor::Clone
  WaveformProcessor *Clone() const override;

  /// @brief The integration gates.
  /// @return The settings this processor was built with.
  const PulseIntegralSettings &GetSettings() const { return settings_; }

protected:
  /// @copydoc WaveformProcessor::BookBranches
  void BookBranches(TTree *tree) override;
  /// @copydoc WaveformProcessor::ExtractFeatures
  Bool_t ExtractFeatures(const TArrayF &cropped) override;

private:
  PulseIntegralSettings settings_;
  Float_t short_integral_;
  Float_t long_integral_;
};

/**
 * @brief The trapezoid filter of PulseHeightProcessor, in samples.
 *
 * The trapezoid is read relative to the trigger t (the crossing of
 * FileProcessingConfig::trigger_threshold): for a step at t it rises from t
 * to t + #rise_samples - 1 and stays flat for #flat_top_samples, and its
 * height is the mean of #peak_average_samples samples from
 * t + #rise_samples - 1 + #peaking_fraction x #flat_top_samples. A real step
 * starts a little before the trigger and takes its rise time to complete, so
 * the peaking point must leave room for both.
 */
struct PulseHeightSettings {
  /// Preamplifier decay constant for the pole-zero deconvolution. `0` or
  /// less skips the deconvolution, for inputs that are already steps.
  Double_t decay_time_samples = 0.0;
  /// Rise (and fall) time of the trapezoid, k.
  Int_t rise_samples = 100;
  /// Flat top of the trapezoid, m.
  Int_t flat_top_samples = 50;
  /// Where on the flat top the height is read, from 0 (its start) to 1 (its
  /// end).
  Double_t peaking_fraction = 0.5;
  /// Trapezoid samples averaged from the peaking point on; at least 1.
  Int_t peak_average_samples = 1;
  /// Start, after the trigger, of the window over which the deconvolved
  /// step's plateau is averaged for the rise time.
  Int_t rise_plateau_start = 150;
  /// Width of that window.
  Int_t rise_plateau_samples = 100;
  /// How far before the trigger the search for the 10% crossing starts.
  Int_t rise_search_back = 50;
};

/**
 * @brief Pulse height of a charge-sensitive preamplifier signal: pole-zero
 *        deconvolution, trapezoid filter, rise time.
 *
 * The cropped waveform is turned into a step by the exact discrete pole-zero
 * deconvolution y[n] = y[n-1] + x[n] - d x[n-1], d = exp(-1 / tau), which
 * maps a decaying exponential of amplitude A onto a step of height A. The
 * trapezoid T[n] is the mean of y over the last k samples minus its mean
 * over the k samples ending k + m earlier (Jordanov and Knoll's trapezoid as
 * a finite impulse response), so a step of height A gives a flat top of
 * height A.
 *
 * Adds `trapezoid_height` (ADC counts) and `rise_time` (the step's 10-90% rise,
 * in samples, linearly interpolated; -1 when not found) to the tree. Rejects
 * waveforms whose trapezoid does not fit inside the crop, which
 * CheckSettings() rules out for the trigger at its nominal place.
 */
class PulseHeightProcessor : public WaveformProcessor {
public:
  /**
   * @brief Construct from the shared configuration and the filter.
   * @param config   Shared configuration, copied.
   * @param settings Trapezoid and rise-time settings, copied.
   */
  PulseHeightProcessor(const FileProcessingConfig &config,
                       const PulseHeightSettings &settings);

  /// @copydoc WaveformProcessor::Clone
  WaveformProcessor *Clone() const override;

  /// @brief The trapezoid and rise-time settings.
  /// @return The settings this processor was built with.
  const PulseHeightSettings &GetSettings() const { return settings_; }

  /**
   * @brief Deconvolve a baseline-subtracted preamplifier signal into a step.
   * @param waveform           Baseline-subtracted samples, positive pulse.
   * @param decay_time_samples Decay constant, in samples; `0` or less returns
   *                           the input unchanged.
   * @return The step, same length as @p waveform.
   */
  static TArrayF Deconvolve(const TArrayF &waveform,
                            Double_t decay_time_samples);

  /**
   * @brief The trapezoid filter of a step.
   * @param step         Step, e.g. from Deconvolve().
   * @param rise_samples Rise time k, at least 1.
   * @param flat_samples Flat top m, at least 0.
   * @return T, same length as @p step. Samples before the start of @p step
   *         count as zero, its baseline.
   */
  static TArrayF Trapezoid(const TArrayF &step, Int_t rise_samples,
                           Int_t flat_samples);

protected:
  /// @copydoc WaveformProcessor::BookBranches
  void BookBranches(TTree *tree) override;
  /// @copydoc WaveformProcessor::CheckSettings
  Bool_t CheckSettings() const override;
  /// @copydoc WaveformProcessor::ExtractFeatures
  Bool_t ExtractFeatures(const TArrayF &cropped) override;
  /// @brief Draws the cropped waveform, its deconvolved step and its
  ///        trapezoid, with the peaking window marked.
  /// @param cropped The cropped, baseline-subtracted waveform.
  void DrawSampleWaveform(const TArrayF &cropped) override;

private:
  PulseHeightSettings settings_;
  Float_t trapezoid_height_;
  Float_t rise_time_;
  Int_t PeakingStart() const;
  Float_t RiseTime(const TArrayF &step) const;
};

#endif
