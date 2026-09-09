#ifndef WAVEFORMPROCESSOR_H
#define WAVEFORMPROCESSOR_H

#include "PlottingUtils.hpp"
#include <TArrayF.h>
#include <TArrayI.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <future>
#include <iostream>
#include <mutex>
#include <thread>
#include <vector>

/**
 * @brief Per-waveform quantities extracted by
 *        WaveformProcessingUtils::ExtractFeatures().
 *
 * Integrals and pulse heights are in baseline-subtracted ADC counts, positive
 * regardless of input polarity. Sample indices are relative to the *cropped*
 * waveform unless stated otherwise.
 */
struct WaveformFeatures {
  /// Extremum of the raw, pre-baseline-subtraction trace, as a magnitude.
  /// Compared against the saturation code to detect clipping.
  Int_t raw_pulse_height;
  /// Maximum of the baseline-subtracted cropped waveform, in ADC counts.
  Float_t pulse_height;
  /// Index of that maximum within the cropped waveform, in samples.
  Int_t peak_position;
  /// Trigger index re-evaluated on the cropped waveform, in samples.
  Int_t trigger_position;
  /// Sum over the short gate, in ADC counts x samples.
  Float_t short_integral;
  /// Sum over the long gate, in ADC counts x samples.
  Float_t long_integral;
  /// Fraction of samples in the long gate that are below zero, in [0, 1].
  /// A large value indicates a bad baseline; see
  /// WaveformProcessingUtils::ApplyQualityCuts().
  Float_t negative_fraction;
  /// Whether this waveform survived the quality cuts.
  Bool_t passes_cuts;
  /// Acquisition timestamp carried through from the input record.
  ULong64_t timestamp;
};

/**
 * @brief Running counters and accumulators for one processing run.
 *
 * The `rejected_*` counters partition the rejected waveforms by cause, so
 * `total_processed - accepted` equals their sum.
 */
struct ProcessingStats {
  Int_t total_processed = 0;     ///< Waveforms examined.
  Int_t accepted = 0;            ///< Waveforms that passed every cut.
  Int_t rejected_no_trigger = 0; ///< No sample reached the trigger level.
  /// Trigger too close to either end to crop the requested window.
  Int_t rejected_insufficient_samples = 0;
  Int_t rejected_negative_integral = 0; ///< Long-gate integral <= 0.
  Int_t rejected_baseline = 0;          ///< More than 50% negative samples.
  Int_t rejected_clipped = 0;           ///< Raw extremum hit the ADC rail.
  Double_t sum_baseline_rms = 0.0;      ///< Sum of per-waveform baseline RMS.
  Int_t baseline_rms_count = 0;         ///< Waveforms contributing to that sum.
  /// Sum of baseline RMS over accepted waveforms only.
  Double_t sum_baseline_rms_accepted = 0.0;
  Int_t baseline_rms_count_accepted = 0; ///< Accepted waveforms in that sum.
};

/// @brief Acquisition software that produced the input file.
enum class InputFormat {
  kCOMPASS,  ///< CAEN CoMPASS.
  kWAVEDUMP, ///< CAEN WaveDump (DT5742 family).
  kSOLARIS   ///< SOLARIS DAQ (SOL).
};

/**
 * @brief Everything needed to configure one processing run.
 *
 * Passed to the WaveformProcessingUtils constructor, and the only way to
 * configure the workers spawned by
 * WaveformProcessingUtils::ProcessFilesParallel(). All sample counts and gate
 * widths are in samples.
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
  /// Samples before `pre_samples` at which integration starts, so integration
  /// begins at `pre_samples - pre_gate` within the cropped waveform.
  Int_t pre_gate = 10;
  /// Width of the short integration gate.
  Int_t short_gate = 10;
  /// Width of the long integration gate.
  Int_t long_gate = 200;
  /// How many example waveforms to write out as figures, for eyeballing.
  Int_t sample_waveforms_to_save = 5;
  /// Cap on waveforms to read per file; `-1` means no limit.
  Int_t max_events = -1;
  Bool_t verbose = kTRUE;         ///< Print progress and per-file summaries.
  Bool_t store_waveforms = kTRUE; ///< Write the cropped waveform to the tree.
  InputFormat input_format = InputFormat::kCOMPASS; ///< Input file format.
  /// ADC code that means full scale. A raw extremum equal to this (positive
  /// polarity) or to zero (negative polarity) is treated as clipped.
  Int_t adc_saturation_code = 16384;
};

/**
 * @brief Waveform processing pipeline: baseline, trigger, crop, features, cuts.
 *
 * Typical use is to construct one from a FileProcessingConfig and call
 * ProcessFile(), or to hand a list of files to ProcessFilesParallel().
 *
 * @warning The pipeline stages — SubtractBaseline(), FindTrigger(),
 *          CropWaveform(), ExtractFeatures() — are public but are **not
 *          standalone utilities**. They communicate through an internal
 *          waveform buffer that each stage overwrites, and they mutate the
 *          statistics counters. Call them in pipeline order, or call
 *          ProcessWaveform(), which sequences them and enforces the
 *          preconditions each one assumes.
 *
 * @note One instance is not safe to share across threads.
 *       ProcessFilesParallel() gives every worker its own instance for exactly
 *       this reason.
 */
class WaveformProcessingUtils {
private:
  Int_t polarity_;
  Double_t trigger_threshold_;
  Int_t num_samples_baseline_;
  Int_t pre_samples_;
  Int_t post_samples_;
  Int_t pre_gate_;
  Int_t short_gate_;
  Int_t long_gate_;
  Int_t max_events_;
  Bool_t verbose_;
  Int_t adc_saturation_code_;

  static std::mutex canvas_mutex_;
  Int_t sample_waveforms_to_save_;
  Int_t sample_waveforms_saved_;
  TString current_output_name_;

  ProcessingStats stats_;

  Float_t current_baseline_rms_ = 0.0f;
  Bool_t current_baseline_rms_valid_ = kFALSE;

  TFile *output_file_;
  TTree *output_tree_;
  WaveformFeatures current_features_;
  Bool_t store_waveforms_;
  TArrayF *save_waveform_;
  ULong64_t current_timestamp_;
  InputFormat input_format_;

public:
  /// @brief Construct with the FileProcessingConfig defaults.
  WaveformProcessingUtils();

  /// @brief Construct from an explicit configuration.
  /// @param config Processing parameters; copied into the instance.
  WaveformProcessingUtils(const FileProcessingConfig &config);

  /// @brief Closes any open output file and releases the internal buffer.
  ~WaveformProcessingUtils();

  /// @brief Set the pulse polarity.
  /// @param polarity `-1` for negative-going pulses, `+1` for positive-going.
  void SetPolarity(const Int_t polarity) { polarity_ = polarity; }

  /// @brief Set the trigger level as a fraction of each waveform's peak.
  /// @param threshold Fraction in (0, 1]. Larger values trigger later on the
  ///                  rising edge.
  void SetTriggerThreshold(Double_t threshold) {
    trigger_threshold_ = threshold;
  }

  /// @brief Set how many leading samples are averaged for the baseline.
  /// @param num_samples_baseline Sample count. Must be at least 2 for a
  ///                             baseline RMS to be computed; with 1 the
  ///                             baseline is still subtracted but no RMS is
  ///                             recorded.
  void SetNumberOfSamplesForBaseline(Int_t num_samples_baseline) {
    num_samples_baseline_ = num_samples_baseline;
  }

  /// @brief Set the crop window around the trigger.
  /// @param pre_samples  Samples kept before the trigger.
  /// @param post_samples Samples kept after the trigger.
  void SetSampleWindows(Int_t pre_samples, Int_t post_samples) {
    pre_samples_ = pre_samples;
    post_samples_ = post_samples;
  }

  /// @brief Set the integration gates, in samples.
  /// @param pre_gate    Offset back from `pre_samples` at which integration
  ///                    starts.
  /// @param short_gate  Width of the short gate.
  /// @param long_gate   Width of the long gate. Should exceed @p short_gate;
  ///                    the PSD ratio is meaningless otherwise.
  void SetGates(Int_t pre_gate, Int_t short_gate, Int_t long_gate) {
    pre_gate_ = pre_gate;
    short_gate_ = short_gate;
    long_gate_ = long_gate;
  }

  /// @brief Cap the number of waveforms read per file.
  /// @param max_events Maximum count, or `-1` for no limit.
  void SetMaxEvents(Int_t max_events) { max_events_ = max_events; }

  /// @brief Enable or disable progress and summary output.
  /// @param verbose `kTRUE` to print.
  void SetVerbose(Bool_t verbose) { verbose_ = verbose; }

  /// @brief Choose whether cropped waveforms are written to the output tree.
  /// @param store `kTRUE` (default) to store them. Turning this off leaves only
  ///              the scalar features and shrinks the output substantially.
  void SetStoreWaveforms(Bool_t store = kTRUE) { store_waveforms_ = store; }

  /// @brief Set how many example waveforms are saved as figures.
  /// @param count Number of figures to write, for visual inspection.
  void SetSaveSampleWaveforms(Int_t count) {
    sample_waveforms_to_save_ = count;
  }

  /// @brief Set the ADC code that counts as full scale.
  /// @param adc_saturation_code Full-scale code, e.g. 16384 for a 14-bit
  ///                            digitiser. Used only for the clipping cut.
  void SetAdcSaturationCode(Int_t adc_saturation_code) {
    adc_saturation_code_ = adc_saturation_code;
  }

  /**
   * @brief Run one raw waveform through the whole pipeline.
   *
   * Subtracts the baseline, finds the trigger, crops, extracts features,
   * applies the quality cuts, and on success fills the output tree. Updates the
   * statistics counters in every case.
   *
   * @param samples Raw digitiser samples for one waveform, in ADC counts.
   *
   * @return `kTRUE` if the waveform was accepted and written; `kFALSE` if it
   *         was empty, failed to trigger, sat too close to an end of the record
   *         to crop, or failed a quality cut. The specific cause is recorded in
   *         the corresponding ProcessingStats counter.
   *
   * @note Requires an open output tree, so it is only meaningful between the
   *       setup and teardown performed by ProcessFile().
   */
  Bool_t ProcessWaveform(const TArrayI &samples);

  /**
   * @brief Estimate and remove the baseline, writing to the internal buffer.
   *
   * Averages the first `num_samples_baseline` samples, subtracts that level,
   * and inverts the result when polarity is `-1` so the pulse is positive
   * either way. Also computes the baseline RMS and folds it into the
   * statistics.
   *
   * @param samples Raw digitiser samples, in ADC counts.
   *
   * @note Writes the result into the internal waveform buffer, replacing
   *       whatever was there. The baseline RMS is only computed when at least
   *       two baseline samples are available.
   */
  void SubtractBaseline(const TArrayI &samples);

  /**
   * @brief Find the first sample crossing the fractional trigger level.
   *
   * The level is `trigger_threshold` times the waveform's own maximum, so it
   * adapts to pulse amplitude rather than being an absolute ADC threshold.
   *
   * @param waveform Baseline-subtracted waveform.
   *
   * @return Index of the first sample at or above the level, as a `Float_t`, or
   *         `-1.0` if no sample reaches it.
   */
  Float_t FindTrigger(const TArrayF &waveform);

  /**
   * @brief Cut the region of interest around the trigger.
   *
   * Keeps `pre_samples` before and `post_samples` after @p trigger_pos, and
   * writes the result into the internal waveform buffer.
   *
   * @param waveform    Baseline-subtracted waveform.
   * @param trigger_pos Trigger sample index, from FindTrigger().
   *
   * @warning @p trigger_pos must be at least `pre_samples`, and the waveform
   *          must extend more than `post_samples` beyond it. The start index is
   *          **not** clamped, so a smaller trigger position reads before the
   *          start of the array. ProcessWaveform() checks both bounds before
   *          calling this; a direct caller must do the same.
   */
  void CropWaveform(const TArrayF &waveform, Int_t trigger_pos);

  /**
   * @brief Compute the per-waveform features of a cropped waveform.
   *
   * Integration starts at `pre_samples - pre_gate` within the crop and runs for
   * `short_gate` and `long_gate` samples respectively, each clipped to the end
   * of the waveform.
   *
   * @param cropped_wf Cropped, baseline-subtracted waveform.
   *
   * @return The extracted features. `raw_pulse_height` is left unset here —
   *         ProcessWaveform() fills it from the raw trace, and the clipping cut
   *         depends on it.
   */
  WaveformFeatures ExtractFeatures(const TArrayF &cropped_wf);

  /**
   * @brief Decide whether a waveform survives the quality cuts.
   *
   * Rejects, in order: a raw extremum sitting on the ADC rail (the saturation
   * code for positive polarity, zero for negative); more than 50% of long-gate
   * samples below zero, which indicates a bad baseline; and a non-positive
   * long-gate integral.
   *
   * @param features Features to test. `raw_pulse_height` must already be set.
   *
   * @return `kTRUE` if the waveform passes every cut.
   *
   * @note Increments the matching ProcessingStats rejection counter on failure.
   */
  Bool_t ApplyQualityCuts(const WaveformFeatures &features);

  /**
   * @brief Write one waveform out as a figure for visual inspection.
   *
   * @param waveform Waveform to plot, normally the cropped one.
   *
   * @note Serialised across threads on an internal mutex, because ROOT canvas
   *       creation is not concurrency-safe. Stops having any effect once
   *       `sample_waveforms_to_save` figures have been written.
   */
  void SaveSampleWaveform(const TArrayF &waveform);

  /**
   * @brief Print the run statistics, including mean baseline RMS.
   *
   * Reports the counters plus two means: baseline RMS over all processed
   * waveforms, and over accepted ones only. The same summary is written to the
   * `.stats` file alongside the output.
   */
  void PrintAllStatistics() const;

  /// @brief Snapshot of the current statistics.
  /// @return A copy of the counters accumulated so far.
  ProcessingStats GetStats() const { return stats_; };

  /**
   * @brief Process one input file end to end.
   *
   * Opens the input in the configured InputFormat, creates the output tree,
   * runs every waveform through ProcessWaveform(), then writes the tree, the
   * sample figures, and the `.stats` summary.
   *
   * @param filepath    Path to the input file.
   * @param output_name Output basename, without extension, resolved against the
   *                    ROOT files base directory.
   *
   * @return `kTRUE` on success; `kFALSE` if the input could not be opened or
   *         the output could not be written.
   */
  Bool_t ProcessFile(const TString filepath, const TString output_name);

  /**
   * @brief Process many files concurrently, one instance per worker.
   *
   * Enables ROOT thread safety and IO::SetThreadSafe(), then dispatches the
   * files in batches of @p max_workers, waiting for each batch to finish before
   * starting the next. Every worker constructs its own
   * WaveformProcessingUtils from @p config, so no state is shared.
   *
   * @param filepaths    Input files.
   * @param output_names Output basenames, positionally matched to
   *                     @p filepaths. Must be the same length.
   * @param config       Configuration copied into every worker.
   * @param max_workers  Concurrent workers. `4` by default; a value of `0` or
   *                     less means hardware concurrency. Capped at the number
   *                     of files.
   *
   * @note Batches are synchronous, so one slow file holds up the start of the
   *       next batch. Per-file success or failure is printed as each completes;
   *       there is no aggregate return value.
   *
   * @warning Enables ROOT thread safety process-wide as a side effect, which
   *          cannot be undone.
   */
  static void ProcessFilesParallel(const std::vector<TString> &filepaths,
                                   const std::vector<TString> &output_names,
                                   const FileProcessingConfig &config,
                                   Int_t max_workers = 4);
};

#endif
