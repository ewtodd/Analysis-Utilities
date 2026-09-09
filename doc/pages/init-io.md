InitUtils and IOUtils {#init-io}
=====================

[TOC]

## InitUtils

Initialization and file conversion utilities.

- `SetROOTPreferences(save_format, plots_dir, root_files_dir, enable_mt = kTRUE)`
  — configure the ROOT environment, set the plot output base via
  `PlottingUtils::SetPlotsBaseDir`, and set the ROOT-files I/O base via
  `IO::SetRootFilesBaseDir`. With `enable_mt` (the default) it also calls
  `IO::SetThreadSafe()` and detaches new histograms from `gDirectory`
  (`TH1::AddDirectory(kFALSE)`), so parallel file processing is safe.

  @note Pass absolute paths so output is anchored to a project root regardless
  of CWD. If `plots_dir` or `root_files_dir` is omitted, a warning is printed
  and the CWD-relative defaults `"plots"` / `"root_files"` are used.

- `ConvertCoMPASSBinToROOT()` — convert CoMPASS binary files to ROOT format
- `ConvertCoMPASSBinToHits()` — read a CoMPASS binary file into an in-memory
  `std::vector<RawHit>` with no ROOT I/O
- `ConvertSOLBinToROOT()` — convert SOLARIS DAQ (SOL) binary files to ROOT;
  writes a `Data_R` tree with all per-block header fields plus `TArrayI` /
  `TArrayC` trace branches for trace-carrying blocks
- `ConvertSOLBinToHits()` — read a SOL file into an in-memory
  `std::vector<SOLHit>` (header fields only; traces are skipped)

## IOUtils

Path-aware ROOT file open helpers. Subpaths resolve against the base directory
configured via `IO::SetRootFilesBaseDir` (or `InitUtils::SetROOTPreferences`);
absolute paths pass through untouched. All filesystem operations use `gSystem`
(no `std::filesystem`).

- `IO::SetRootFilesBaseDir(dir)` — set the base directory (default
  `"root_files"`). Trailing slashes are stripped.
- `IO::GetRootFilesBaseDir()` — return the current base directory.
- `IO::OpenForReading(subpath)` — returns a `TFile*` opened in `"READ"` at
  `<base>/<subpath>` (or just `subpath` if absolute).
- `IO::OpenForWriting(subpath, mode = "RECREATE")` — same join semantics, plus
  creates parent directories via `gSystem->mkdir(..., kTRUE)` before opening.

### Thread safety

- `IO::SetThreadSafe(enabled = kTRUE)` / `IO::IsThreadSafe()` — enable ROOT
  thread safety (`ROOT::EnableThreadSafety()`) and route all file opens through
  a shared recursive mutex. Enabled by `InitUtils::SetROOTPreferences` (via
  `enable_mt`) and by the parallel file processing in
  @ref waveform-processing "WaveformProcessingUtils".
- `IO::ScopedRootLock` — RAII guard for code that opens `TFile`s directly
  instead of through `IO::OpenForReading` / `IO::OpenForWriting`; engages the
  same lock when thread-safe mode is on.
