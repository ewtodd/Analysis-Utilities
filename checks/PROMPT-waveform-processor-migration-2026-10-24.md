# Prompt: move the repositories onto the waveform processor variants
<!---->
For an agent updating the code that uses Analysis-Utilities' old
`WaveformProcessingUtils` class. Hand it this file once the waveform processor
variants (on `main` since 2026-10-24, f97de2c and 69098de) have been reviewed, pushed
and given a version bump.
<!---->
---
<!---->
Analysis-Utilities replaced its `WaveformProcessingUtils` class with two class families
in `include/WaveformProcessingUtils.hpp` (the header keeps its name). Read that header and
`doc/pages/waveform-processing.md` in `~/Analysis-Utilities` first; this prompt lists
only what callers must change. Update each repository below, verify, and commit
following that repository's AGENTS.md (Conventional Commits with scope, author
`son-of-anton-bot <307402699+son-of-anton-bot@users.noreply.github.com>`, no
co-author trailer). Bump the Analysis-Utilities input with
`nix flake lock --update-input utils` (UM-ANSG) or the repository's equivalent; never
edit a lock file by hand.
<!---->
## What changed for a caller
<!---->
| Old | New |
|---|---|
| `WaveformProcessingUtils p(config); p.ProcessFile(in, out);` | `PulseIntegralProcessor p(config, gates); p.ProcessFile(in, out);` |
| `WaveformProcessingUtils::ProcessFilesParallel(files, names, config, n)` (static) | `PulseIntegralProcessor(config, gates).ProcessFilesParallel(files, names, n)` (a member; each worker clones the processor) |
| `FileProcessingConfig::pre_gate`, `short_gate`, `long_gate` | `PulseIntegralSettings::pre_gate`, `short_gate`, `long_gate` (same meaning, same defaults) |
| `WaveformFeatures` | gone; the features are the output tree's branches |
| `SetPolarity()`, `SetGates()`, other setters | gone; set the config and settings structs before constructing |
| `ProcessWaveform()`, `SubtractBaseline()`, `FindTrigger()`, `CropWaveform()`, `ExtractFeatures()`, `ApplyQualityCuts()` | private; use `ProcessFile()`, or `ProcessSource()` with a `WaveformSource` |
| `ROOT.WaveformProcessingUtils(cfg)` (PyROOT) | `ROOT.PulseIntegralProcessor(cfg, ROOT.PulseIntegralSettings())` |
<!---->
`FileProcessingConfig` keeps every other field, in the same order, and gains
`channel` (default `-1`, every channel) at the end. Designated initializers that
listed the gates must drop them: C++ requires the remaining designators in declaration
order, and the gates now go in a `PulseIntegralSettings`.
<!---->
The output tree `features` keeps `pulse_height`, `trigger_position`, `short_integral`,
`long_integral`, `timestamp` and `Samples` with the same names and types, and gains
`channel`, `board`, `digitizer_energy`, `flags` and `baseline_rms`. Readers need no
change. Two behaviours differ: the statistics reset at the start of every
`ProcessFile` call, and records with no samples are skipped instead of counted as
processed.
<!---->
The variant reproduces the old class exactly (YAP-PSD and CeYAG data, every branch of
every event). Verify that for each repository: before updating, process a subset (set
`max_events`, e.g. 20000, and an output name of your own) with the old library; after,
with the new one; compare every branch of every event. Report the comparison; a
difference is a bug to report, not to paper over.
<!---->
## Repositories
<!---->
1. **UM-ANSG** (`~/UM-ANSG`; enter shells with `nix develop .#<Project>` from the
   repository root, never `nix-shell`). Do not touch `73mGe`, which already uses the new
   API; do not `git add .` (the 73mGe shell decrypts secrets into the tree).
   - `YAP-PSD/macros/Constants.hpp:105-118`: split `DEFAULT_PROCESSING_CONFIG` into a
     `FileProcessingConfig` without the gates and a `PulseIntegralSettings`
     (`pre_gate 5, short_gate 10, long_gate 250`).
   - `YAP-PSD/macros/InitialProcessing.cpp:69`: the parallel call.
   - `YAP-PSD/macros/NoiseRMS.cpp:20` and `YAP-PSD/python/noise_study.py:99`, `:989`:
     comments naming `WaveformProcessingUtils::SubtractBaseline` and
     `ExtractFeatures`; point them at `WaveformProcessor`'s baseline step and
     `PulseIntegralProcessor::ExtractFeatures`.
   - `78mBr/macros/Constants.hpp:68-81` and `78mBr/macros/InitialProcessing.cpp:74`: as
     YAP-PSD.
   - `YAG-PreProcessing/macros/Constants.hpp:37-51` and
     `YAG-PreProcessing/macros/WaveformViewer.cpp:10`: as YAP-PSD; it reads WaveDump
     ROOT trees (`InputFormat::kWAVEDUMP`), and `kWAVEDUMP742_BINARY` would now read
     the `.dat` files directly.
2. **son-of-anton** (`~/Software/son-of-anton`):
   - `autophysicist/utils/runtime_notes.py:91-130`: the note telling the agent how to
     process raw waveforms. Its notes are selected by Analysis-Utilities version
     (`notes_for({"analysis_utilities": ...})`): keep the old text for versions before
     this change and add the new one (`PulseIntegralProcessor`, `PulseIntegralSettings`,
     the member `ProcessFilesParallel`, `ProcessSource` instead of the old
     "do not call ProcessWaveform" warning) from the bumped version on.
   - `autophysicist/utils/api_docs.py:366`, `:378`, `:407-408`: the example symbols.
   - `tests/test_physics_sandbox.py:313`: the assertion on the note's text.
   - `README.md:281`: the class name.
3. **check-source** (`~/Software/check-source`): `src/check_source/mining/draft.py:234`
   names the source file `WaveformProcessingUtils`, which keeps its name; confirm and
   leave it.
<!---->
Finish with a report per repository: what changed, the before/after comparison, the
commits, and anything that did not match.
