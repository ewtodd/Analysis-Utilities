#ifndef INITUTILS_H
#define INITUTILS_H

#include "BinaryUtils.hpp"
#include "IOUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <cstdio>
#include <cstring>
#include <utility>

/**
 * @brief Environment setup and binary-to-ROOT file conversion.
 *
 * All-static. SetROOTPreferences() is the recommended single entry point for
 * configuring the ROOT environment at the top of an analysis; the remaining
 * methods convert acquisition-software binary files into ROOT trees or
 * in-memory hit vectors.
 */
class InitUtils {
public:
  /**
   * @brief Configure the ROOT environment and pin all output paths.
   *
   * Applies the house plotting style, forces it over any style carried by
   * objects read from file, and puts ROOT in batch mode so plots do not pop up.
   * Sets the plot output base via PlottingUtils::SetPlotsBaseDir() and the ROOT
   * file I/O base via IO::SetRootFilesBaseDir().
   *
   * @param save_format    Figure output format. PlotSaveFormat::kPNG by
   *                       default.
   * @param plots_dir      Base directory for saved figures. **Pass an absolute
   *                       path** so output is anchored to a project root
   *                       regardless of the current working directory. If
   *                       empty, a warning is printed and the CWD-relative
   *                       `"plots"` is used.
   * @param root_files_dir Base directory for ROOT file I/O. Same absolute-path
   *                       advice; if empty, a warning is printed and the
   *                       CWD-relative `"root_files"` is used.
   * @param enable_mt      When `kTRUE` (the default), enables ROOT thread
   *                       safety via IO::SetThreadSafe() and detaches new
   *                       histograms from `gDirectory` with
   *                       `TH1::AddDirectory(kFALSE)`, so concurrent `TFile`
   *                       opens on other threads cannot race on a directory's
   *                       child list. Required for the parallel file processing
   *                       in WaveformProcessingUtils.
   *
   * @note Mutates process-wide ROOT state: `gStyle`, `gROOT` batch mode and
   *       forced style, and — with @p enable_mt — histogram directory
   *       ownership. Call it once, early.
   *
   * @warning With @p enable_mt on, histograms are no longer owned by the
   *          current `TFile`. Anything you want written to disk must be
   *          `Write()`-ten explicitly, and every histogram you create becomes
   *          your responsibility to delete.
   */
  static void
  SetROOTPreferences(PlotSaveFormat save_format = PlotSaveFormat::kPNG,
                     const TString &plots_dir = "",
                     const TString &root_files_dir = "",
                     Bool_t enable_mt = kTRUE);

  /**
   * @brief Convert a WaveDump binary file (DT5742 family) to a ROOT tree.
   *
   * @param input_filename     Path to the WaveDump binary file.
   * @param output_name        Output basename, without extension. The file is
   *                           written to
   * `<root_files_base>/<output_name>.root`.
   * @param corrections_enabled When `kTRUE` (the default), applies the DT5742
   *                           timing corrections during conversion.
   *
   * @return `kTRUE` on success; `kFALSE` if the input is missing or cannot be
   *         parsed, with the reason printed to stdout.
   */
  static Bool_t ConvertWavedumpBinToROOT(const TString input_filename,
                                         const TString output_name,
                                         Bool_t corrections_enabled = kTRUE);

  /**
   * @brief Convert a CoMPASS binary file to a ROOT tree.
   *
   * Which branches are written is dictated by the file's global header, whose
   * low bits flag the fields present in each record: `0x0001` energy (channel
   * units), `0x0002` calibrated energy, `0x0004` short-gate energy, `0x0008`
   * waveform samples.
   *
   * @param input_filename         Path to the CoMPASS binary file.
   * @param output_name            Output basename, without extension. Written
   *                               to `<root_files_base>/<output_name>.root`;
   *                               the base directory is created if missing.
   * @param global_header_override Header word to assume instead of reading one
   *                               from the file. Pass `0` to read it from the
   *                               file, which is the normal case. Use a
   *                               non-zero value only for files whose two-byte
   *                               header is absent or wrong — passing one also
   *                               suppresses the two-byte header skip.
   * @param skip_bad_events        When `kTRUE`, records that fail validation
   *                               are skipped instead of aborting the
   *                               conversion. `kFALSE` by default.
   *
   * @return The global header actually used, which is non-zero on success, or
   *         `0` if the input does not exist, cannot be opened, or the
   *         conversion failed. Inspect the returned bits to learn which
   *         branches the output tree carries.
   */
  static UShort_t ConvertCoMPASSBinToROOT(const TString input_filename,
                                          const TString output_name,
                                          UShort_t global_header_override,
                                          Bool_t skip_bad_events = kFALSE);

  /**
   * @brief Read a CoMPASS binary file into memory, with no ROOT file I/O.
   *
   * The in-memory counterpart of ConvertCoMPASSBinToROOT(), for callers that
   * want the hits directly rather than a tree on disk.
   *
   * @param input_filename         Path to the CoMPASS binary file.
   * @param global_header_override Header word to assume instead of reading one
   *                               from the file; `0` (the default) reads it
   *                               from the file.
   * @param skip_bad_events        When `kTRUE`, records that fail validation
   *                               are skipped rather than aborting.
   *
   * @return A pair of the hits read and the global header used. On failure the
   *         vector is empty and the header is `0`.
   */
  static std::pair<std::vector<RawHit>, UShort_t>
  ConvertCoMPASSBinToHits(const TString input_filename,
                          UShort_t global_header_override = 0,
                          Bool_t skip_bad_events = kFALSE);

  /**
   * @brief Convert a SOLARIS DAQ (SOL) binary file to a ROOT tree.
   *
   * Writes a `Data_R` tree carrying every per-block header field, plus
   * `TArrayI` / `TArrayC` trace branches for the block formats that include
   * traces.
   *
   * @param input_filename Path to the SOL binary file.
   * @param output_name    Output basename, without extension. Written to
   *                       `<root_files_base>/<output_name>.root`.
   *
   * @return `kTRUE` on success; `kFALSE` if the input is missing or cannot be
   *         opened, with the reason printed to stdout.
   */
  static Bool_t ConvertSOLBinToROOT(const TString input_filename,
                                    const TString output_name);

  /**
   * @brief Read a SOL binary file into memory as lightweight hits.
   *
   * Traces are skipped entirely — each block is reduced to its header fields —
   * which makes this cheap enough to run over a whole run when only timing and
   * energy are needed.
   *
   * @param input_filename Path to the SOL binary file.
   *
   * @return A pair of the hits read and the total number of blocks processed.
   *         On failure the vector is empty and the count is `0`.
   */
  static std::pair<std::vector<SOLHit>, Long64_t>
  ConvertSOLBinToHits(const TString input_filename);
};

#endif
