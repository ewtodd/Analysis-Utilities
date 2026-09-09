WaveformProcessingUtils {#waveform-processing}
=======================

Process raw waveforms and extract physical parameters. Outputs processed data to
ROOT TTrees with optional waveform storage.

- **Baseline subtraction** — configurable number of pre-trigger samples
- **Trigger finding** — fraction-of-peak threshold with configurable polarity
- **Waveform cropping** — extract region of interest around trigger
- **Feature extraction** — pulse height, peak position, short/long integrals,
  PSD ratio
- **Quality cuts** — reject clipped signals, baseline issues, negative integrals
- **Baseline RMS** — per-waveform baseline RMS is tracked, and the means (all
  processed, and accepted only) are reported to the console and the `.stats`
  file alongside the rejection counters
- **Parallel file processing** — process multiple files concurrently using
  `std::async` with configurable worker count (defaults to 4). Files are
  dispatched in batches, with each worker getting its own
  `WaveformProcessingUtils` instance. Requires ROOT thread safety
  (`ROOT::EnableThreadSafety()`), which is handled automatically.
