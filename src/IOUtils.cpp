#include "IOUtils.hpp"

namespace {
TString g_root_files_base = "root_files";

TString JoinPath(const TString &subpath) {
  if (gSystem->IsAbsoluteFileName(subpath))
    return subpath;
  char *tmp = gSystem->ConcatFileName(g_root_files_base.Data(), subpath.Data());
  TString full(tmp);
  delete[] tmp;
  return full;
}
} // namespace

void IO::SetRootFilesBaseDir(const TString &dir) {
  TString d = dir;
  while (d.Length() > 0 && d[d.Length() - 1] == '/')
    d.Chop();
  g_root_files_base = d;
}

TString IO::GetRootFilesBaseDir() { return g_root_files_base; }

TFile *IO::OpenForReading(const TString &subpath) {
  TString full = JoinPath(subpath);
  return new TFile(full, "READ");
}

TFile *IO::OpenForWriting(const TString &subpath, const TString mode) {
  TString full = JoinPath(subpath);
  TString parent = gSystem->DirName(full);
  gSystem->mkdir(parent, kTRUE);
  return new TFile(full, mode);
}
