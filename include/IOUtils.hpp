#ifndef IOUTILS_H
#define IOUTILS_H

#include <TFile.h>
#include <TString.h>
#include <TSystem.h>

namespace IO {
void SetRootFilesBaseDir(const TString &dir);
TString GetRootFilesBaseDir();
TFile *OpenForReading(const TString &subpath);
TFile *OpenForWriting(const TString &subpath, const char *mode = "RECREATE");
} // namespace IO

#endif
