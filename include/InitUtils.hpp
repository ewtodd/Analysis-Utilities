#ifndef INITUTILS_H
#define INITUTILS_H

#include "BinaryUtils.hpp"
#include "PlottingUtils.hpp"
#include <TFile.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>

class InitUtils {
public:
  static void SetROOTPreferences(Bool_t setupPlotting = kTRUE);
  static Bool_t ConvertWavedumpBinToROOT();
  static UShort_t ConvertCoMPASSBinToROOT(const TString input_filename,
                                          const TString output_name,
                                          UShort_t global_header_override,
                                          Bool_t skip_bad_events = kFALSE);
  static Bool_t ConvertCoMPASSCSVToROOT();
};

#endif
