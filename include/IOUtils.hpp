#ifndef IOUTILS_H
#define IOUTILS_H

#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

/**
 * @brief Path-aware ROOT file helpers with optional thread-safe opening.
 *
 * Every path passed to this namespace is interpreted relative to a single
 * configurable base directory, so a project can pin its ROOT I/O to an absolute
 * root once and have every downstream call land in the same tree regardless of
 * the current working directory. Absolute paths bypass the base entirely.
 *
 * All filesystem work goes through `gSystem` rather than `std::filesystem`.
 */
namespace IO {

/**
 * @brief Set the base directory that relative subpaths resolve against.
 *
 * @param dir Base directory. May be relative or absolute; pass an absolute path
 *            so output is anchored to a project root regardless of the current
 *            working directory. Trailing slashes are stripped.
 *
 * @note The default is the CWD-relative `"root_files"`.
 *       InitUtils::SetROOTPreferences() calls this for you.
 */
void SetRootFilesBaseDir(const TString &dir);

/// @brief Current base directory for relative subpaths, without trailing slash.
/// @return The directory last set by SetRootFilesBaseDir(), or `"root_files"`.
TString GetRootFilesBaseDir();

/**
 * @brief Enable ROOT thread safety and serialise file opening.
 *
 * When enabled this calls `ROOT::EnableThreadSafety()` and routes every
 * OpenForReading() / OpenForWriting() call through a shared recursive mutex.
 *
 * @param enabled `kTRUE` to enable, `kFALSE` to disable. Note that disabling
 *                only stops this namespace from taking the lock; ROOT's own
 *                thread-safety mode cannot be switched back off once enabled.
 *
 * @note Enabled by InitUtils::SetROOTPreferences() (via its `enable_mt`
 *       argument) and by the parallel file processing in
 *       WaveformProcessingUtils.
 */
void SetThreadSafe(Bool_t enabled = kTRUE);

/// @brief Whether thread-safe file opening is currently engaged.
/// @return `kTRUE` if SetThreadSafe(kTRUE) was called.
Bool_t IsThreadSafe();

/**
 * @brief Open a ROOT file for reading, resolved against the base directory.
 *
 * @param subpath Path to the file. Joined onto the base directory unless it is
 *                absolute, in which case it is used as-is.
 *
 * @return A newly allocated `TFile` opened in `"READ"` mode. **The caller owns
 *         this pointer** and is responsible for `Close()` and `delete`. Never
 *         null, but may be a zombie — check `IsZombie()` before use, as ROOT
 *         reports a missing or corrupt file that way rather than by returning
 *         null.
 *
 * @note In thread-safe mode the open is serialised and `gDirectory` is
 *       preserved across the call via `TDirectory::TContext`.
 */
TFile *OpenForReading(const TString &subpath);

/**
 * @brief Open a ROOT file for writing, creating parent directories first.
 *
 * Parent directories of the resolved path are created recursively before the
 * file is opened, so a caller need not pre-create the output tree.
 *
 * @param subpath Path to the file. Joined onto the base directory unless it is
 *                absolute, in which case it is used as-is.
 * @param mode    ROOT open mode, passed through to `TFile`. `"RECREATE"` by
 *                default; `"UPDATE"` and `"NEW"` behave as ROOT documents them.
 *
 * @return A newly allocated `TFile`. **The caller owns this pointer** and is
 *         responsible for `Close()` and `delete`. Check `IsZombie()` before
 *         use.
 *
 * @note The returned file is configured for ZSTD compression at level 5. This
 *       is set before any branches are created so that the setting propagates
 *       to the baskets of every `TTree` written into it.
 *
 * @warning Unlike OpenForReading(), this does **not** preserve `gDirectory`.
 *          Opening a file for writing makes it the current ROOT directory, so
 *          histograms created afterwards will be owned by it unless
 *          `TH1::AddDirectory(kFALSE)` is in effect (which
 *          InitUtils::SetROOTPreferences() sets when `enable_mt` is on).
 */
TFile *OpenForWriting(const TString &subpath, const TString mode = "RECREATE");

/**
 * @brief RAII guard engaging the same lock used by the open helpers.
 *
 * Use this around code that constructs `TFile` objects directly instead of
 * going through OpenForReading() / OpenForWriting(), so that hand-rolled opens
 * serialise against the helpers rather than racing them.
 *
 * The guard is a no-op when thread-safe mode is off, and the underlying mutex
 * is recursive, so nesting a guard inside a helper call is safe.
 *
 * @note Non-copyable. The engaged/idle decision is made at construction, so a
 *       guard created before SetThreadSafe(kTRUE) stays idle for its lifetime.
 */
class ScopedRootLock {
public:
  /// @brief Acquire the shared file-open lock if thread-safe mode is on.
  ScopedRootLock();
  /// @brief Release the lock if this guard acquired one.
  ~ScopedRootLock();

private:
  Bool_t engaged_;
  ScopedRootLock(const ScopedRootLock &);
  ScopedRootLock &operator=(const ScopedRootLock &);
};
} // namespace IO

#endif
