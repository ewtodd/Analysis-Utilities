# Waveform processor variants — review checklist, 2026-10-24
<!---->
Branch `review/waveform-processor-variants`, on top of `review/simultaneous-fixed-lines`
(`6ceba20`).
**Merged into `main` on 2026-10-24 before ewt's review**, which is in
progress; the boxes below are still open. `WaveformProcessingUtils` was one class for
scintillator pulses read from CoMPASS-style ROOT trees. 73mGe's PSI CdTe and Si
detectors and MUSIC's SOLARIS traces are charge-sensitive preamplifier signals, which
need pole-zero deconvolution, a trapezoid and a rise time, and MUSIC's come as `.sol`
binaries.
The header now holds two class families in the manner of `BinaryReader`:
`WaveformSource` (one per ROOT layout and per binary format, `CreateWaveformSource()`
to pick one) and `WaveformProcessor`, the shared pipeline, with `PulseIntegralProcessor`
(the old features) and `PulseHeightProcessor` (the new ones) as variants.
The API
changes: `WaveformProcessingUtils` and `WaveformFeatures` are gone, the gates moved from
`FileProcessingConfig` into `PulseIntegralSettings`, and `ProcessFilesParallel` is a
member that clones the processor.
`PROMPT-waveform-processor-migration-2026-10-24.md`
is the brief for updating the repositories that use the old API.
Publishing this needs
a version bump: son-of-anton gates its runtime notes on the version.
<!---->
Line numbers are from `main` at `69098de`, where every open checklist points since the
review branches went onto `main` on 2026-10-24.
No
code is reproduced here; open the file at the line and review.
Checks indicate review
by ewt.
<!---->
## Core analysis logic — review these lines
<!---->
- [x] `src/WaveformProcessingUtils.cpp:311-335` — `ProcessSource`: the statistics and
      the figure count reset at the start of every run (before, a second `ProcessFile`
      on one instance added to the first's counters); records of other channels are
      skipped uncounted when `FileProcessingConfig::channel` is set; `max_events` still
      counts accepted waveforms. Check: resetting per run is what you want.
- [x] `src/WaveformProcessingUtils.cpp:379-442` — `ProcessRecord`, the old
      `ProcessWaveform` with the features split out: raw extremum, baseline, trigger,
      the two crop bounds, crop, `pulse_height` and `trigger_position` on the crop,
      then the clipping cut before the variant's `ExtractFeatures`, which keeps the
      old order of the cuts (clipped, bad baseline, negative integral).
      Check: the
      clipping cut moving into the base, ahead of the variant, rejects the same
      waveforms with the same counter.
- [x] `src/WaveformProcessingUtils.cpp:444-501` — `SubtractBaseline`, `FindTrigger`,
      `CropWaveform`, unchanged in substance.
      The trigger level is still formed in
      double and stored as float, so the trigger sample is the same to the last bit.
- [x] `src/WaveformProcessingUtils.cpp:661-691` — `PulseIntegralProcessor::
      ExtractFeatures`: the old gates and the old two cuts, verbatim.
- [x] `src/WaveformProcessingUtils.cpp:699-703` — `CrossingTime`: the linear crossing
      between samples `i-1` and `i`; `i` itself when `i-1` is already above the level
      (the search window starting on the edge), so it never divides by zero.
      ewt: i don't love where this is or how it is defined (random namespace)
- [x] `src/WaveformProcessingUtils.cpp:721-739` — `Deconvolve`: `y[n] = y[n-1] + x[n]
      - d x[n-1]`, `d = exp(-1/tau)`, the exact inverse of a discrete exponential
      (`x[n] = d x[n-1] + i[n]` gives `y[n] = sum i`); `tau <= 0` passes the input
      through. Check: this rather than the continuous `x + (1/tau) sum x`, which
      leaves a 1/(2 tau) gain error and a slope on the plateau.
- [x] `src/WaveformProcessingUtils.cpp:741-764` — `Trapezoid`: `T[n]` = mean of `y`
      over the last `k` samples minus its mean over the `k` samples ending `k + m`
      earlier, from prefix sums in double; samples before the crop count as zero
      (the baseline), so the trapezoid of a pulse near the crop's start is drawn
      whole.
      A step of height A gives a flat top of A from `s + k - 1` to
      `s + k + m - 1`.
      Check: zero-padding is right for the figures and leaves
      the height untouched, since `CheckSettings` keeps the peaking windows inside
      the crop.
      ewt: the trapezoid plots from this looked correct and so did the pulse height spectra. i did not verify the formula explicitly.
- [x] `src/WaveformProcessingUtils.cpp:766-770` — `PeakingStart`: the trigger (at
      `pre_samples` in the crop) plus `k - 1` plus `peaking_fraction x m`, rounded.
      Check: anchoring on the trigger (the `trigger_threshold` crossing) means a real
      step starts a few samples earlier and completes its rise later; the settings'
      documentation says the peaking point must leave room for both.
- [x] `src/WaveformProcessingUtils.cpp:772-802` — `CheckSettings`: the settings' ranges,
      and that the trapezoid at its peaking window and the rise-time windows fit the
      crop for the trigger at its nominal place; a run with settings that do not fit
      does not start.
- [x] `src/WaveformProcessingUtils.cpp:804-826` — `RiseTime`: the plateau is the mean
      step over `[trigger + rise_plateau_start, + rise_plateau_samples)`; the 10% and
      90% crossings are searched from `rise_search_back` before the trigger; -1 when
      the plateau is not positive or a crossing is not found.
- [x] `src/WaveformProcessingUtils.cpp:828-843` — `PulseHeightProcessor::
      ExtractFeatures`: `trapezoid_height` is the mean over `peak_average_samples` from
      the peaking point; a crop too short for it (a trigger late in a short record) is
      rejected as insufficient samples; no other cut, so noise triggers are kept with
      their small heights.
<!---->
## Tooling
<!---->
- [x] `include/WaveformProcessingUtils.hpp:42-269` — `InputFormat` (three new binary
      formats), `WaveformRecord`, `WaveformSource` and its six implementations, and
      `CreateWaveformSource`.
- [x] `src/WaveformProcessingUtils.cpp:17-63` — the ROOT sources' shared part: opens
      `Data_R`, binds the waveform branch through the derived class, reads scalars
      through their leaves whatever their stored type (CoMPASS writes `Channel` as a
      short, the InitUtils conversion as an unsigned short).
- [x] `src/WaveformProcessingUtils.cpp:65-174` — CoMPASS, WaveDump and SOLARIS ROOT
      sources. The 64-bit `Timestamp` is read through the leaf's value pointer, not
      as a double, which would lose picoseconds past 2.5 hours.
      Records with no
      samples are skipped (before, they counted as processed with no rejection
      counter, breaking the counters' partition).
- [x] `src/WaveformProcessingUtils.cpp:178-273` — the binary sources over
      `CoMPASSReader` (records without a waveform skipped; the channel energy, else
      the calibrated one), `WaveDump742Reader` and `SOLReader` (the first analogue
      trace; flags as high word over low word), and the factory.
- [x] `include/WaveformProcessingUtils.hpp:279-505` — `FileProcessingConfig` (gates
      moved out, `channel` added), `ProcessingStats`, and `WaveformProcessor`: public
      `Clone`, `ProcessFile`, `ProcessSource`, `ProcessFilesParallel`, statistics;
      protected virtual `BookBranches`, `CheckSettings`, `ExtractFeatures`,
      `DrawSampleWaveform`. The old public pipeline stages are private now; their
      documentation warned they were not standalone.
- [x] `include/WaveformProcessingUtils.hpp:508-662` — `PulseIntegralSettings`,
      `PulseIntegralProcessor`, `PulseHeightSettings`, `PulseHeightProcessor` (with
      `Deconvolve` and `Trapezoid` public and static, usable on their own).
- [x] `src/WaveformProcessingUtils.cpp:337-377` — `OpenOutput` and `CloseOutput`: the
      tree `features` with the common branches (`timestamp`, `channel`, `board`,
      `digitizer_energy`, `flags`, `trigger_position`, `pulse_height`,
      `baseline_rms`), the variant's, and `Samples`.
      The old five branches keep their
      names and types.
- [x] `src/WaveformProcessingUtils.cpp:603-643` — `ProcessWithClone` and
      `ProcessFilesParallel`: each worker clones the processor, so the settings travel
      with it; the batching is unchanged.
- [x] `src/WaveformProcessingUtils.cpp:845-900` — the pulse-height sample figure:
      waveform, deconvolved step and trapezoid, the peaking window marked.
- [x] `doc/pages/waveform-processing.md` — rewritten for the two families; `python/
      analysis_utilities/__init__.py:50-54`, `doc/pages/python.md:59`,
      `doc/pages/init-io.md:51`, `include/IOUtils.hpp:48`, `include/InitUtils.hpp:49`
      — the new names; `python/analysis_utilities/tests/test_load_cpp_library.py:34`,
      `:99-104` — the symbols the PyROOT test requires.
<!---->
## Checked and excluded
<!---->
- `PulseIntegralProcessor` reproduces the old class exactly: YAP-PSD `DataR_Am241.root`
  with YAP-PSD's configuration (20000 accepted of 29611, same rejection counters and
  baseline RMS) and a CeYAG DT5742 Am-241 run with YAG's (20000 of 64116), every
  branch of every event identical, the stored waveforms included.
  The DT5742 binary
  source on the `.dat` gives the same 20000 events as the ROOT source on its
  conversion.
- `PulseHeightProcessor` on synthetic preamplifier pulses (charge 5000 over 25 samples,
  tau 325, k 200, m 100, peaking 0.9): trapezoid height 5000.000 (rms 0.03 from integer
  sampling) and rise time 20.000 samples (truth 0.8 x 25) without noise; 4999.99 +-
  0.01 with 5 ADC of noise.
  The synthetic source was a class derived from
  `WaveformSource` inside a ROOT macro, so the extension point works.
- The SOLARIS binary source reads MUSIC 88Sr `.sol` traces (4000 samples, channel,
  energy, flags).
  The CoMPASS binary source is untested on waveforms: the only `.BIN`
  files at hand (MUSIC 87Rb) carry energies only (header 0xCAE1), and it returns no
  records from them, as it should.
- `make docs` leaves `doxygen-warnings.log` empty; the PyROOT tests pass (15).
