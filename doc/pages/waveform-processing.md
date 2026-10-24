Waveform processing {#waveform-processing}
===================

`WaveformProcessingUtils.hpp` turns digitiser waveforms into one ROOT tree of
per-waveform features. It has two halves, each a small class family in the
manner of the @ref binary-utils "binary readers": a **source** reads records
from one input format, and a **processor** variant computes the features.

Sources
-------

A `WaveformSource` hands over one `WaveformRecord` at a time: the samples, the
timestamp, board and channel, the digitiser's own energy and its status flags
(whatever the format records). `CreateWaveformSource(InputFormat)` picks one:

| `InputFormat` | Source | Reads |
|---|---|---|
| `kCOMPASS` | `CoMPASSRootSource` | CoMPASS's ROOT tree `Data_R` (or `InitUtils::ConvertCoMPASSBinToROOT`'s) |
| `kWAVEDUMP` | `WaveDumpRootSource` | `InitUtils::ConvertWavedumpBinToROOT`'s tree |
| `kSOLARIS` | `SolarisRootSource` | `InitUtils::ConvertSOLBinToROOT`'s tree, `Trace0` |
| `kCOMPASS_BINARY` | `CoMPASSBinarySource` | CoMPASS `.bin` through `CoMPASSReader` |
| `kWAVEDUMP742_BINARY` | `WaveDump742BinarySource` | DT5742 binaries through `WaveDump742Reader` |
| `kSOLARIS_BINARY` | `SolarisBinarySource` | SOLARIS `.sol` through `SOLReader`, the first analogue trace |

The binary sources need no conversion step. A format not listed needs only a
class deriving from `WaveformSource` (`Open`, `Next`, `Close`), passed to
`WaveformProcessor::ProcessSource()`. A file holding several channels is
processed one channel at a time with `FileProcessingConfig::channel`.

Processors
----------

`WaveformProcessor` runs every record through the shared steps and writes the
tree `features`:

- **Polarity and clipping** — the raw extremum is checked against the ADC rail
  (`adc_saturation_code` for positive pulses, zero for negative).
- **Baseline** — the mean of the leading `num_samples_baseline` samples is
  subtracted, and negative pulses are inverted; the baseline RMS is kept per
  waveform and summarised in the `.stats` file.
- **Trigger** — the first sample at `trigger_threshold` of the waveform's peak.
- **Crop** — `pre_samples` before the trigger and `post_samples` after it.
- **Common branches** — `timestamp`, `channel`, `board`, `digitizer_energy`,
  `flags`, `trigger_position`, `pulse_height` (the crop's maximum),
  `baseline_rms`, and with `store_waveforms` the cropped waveform `Samples`.

The variant adds its own features and cuts:

- **`PulseIntegralProcessor`** (`PulseIntegralSettings`): `short_integral` and
  `long_integral` over gates starting `pre_gate` before the trigger, for pulse
  shape discrimination. Rejects waveforms with more than half the long gate
  below zero, and non-positive long integrals.
- **`PulseHeightProcessor`** (`PulseHeightSettings`): for charge-sensitive
  preamplifiers. The exact discrete pole-zero deconvolution
  (`decay_time_samples`) turns the exponential into a step; a trapezoid filter
  (`rise_samples`, `flat_top_samples`) is read at `peaking_fraction` of its
  flat top, averaged over `peak_average_samples`, as `trapezoid_height`; and
  the step's 10–90% rise is `rise_time`, in samples. Its sample figures show
  the waveform, the step and the trapezoid together.

```cpp
FileProcessingConfig config;
config.polarity = 1;
config.num_samples_baseline = 500;
config.pre_samples = 500;
config.post_samples = 1000;
config.input_format = InputFormat::kSOLARIS_BINARY;
config.channel = 12;

PulseHeightSettings filter;
filter.decay_time_samples = 2400;
filter.rise_samples = 200;
filter.flat_top_samples = 100;
filter.peaking_fraction = 0.9;

PulseHeightProcessor processor(config, filter);
processor.ProcessFile("run_000.sol", "run_000_ch12");
```

**Parallel file processing** — `ProcessFilesParallel()` processes many files
concurrently, a copy of the processor (`Clone()`) per worker, in batches of
the requested worker count (4 by default). It turns on ROOT thread safety
(`ROOT::EnableThreadSafety()`) itself.
