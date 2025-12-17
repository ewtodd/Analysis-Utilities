#ifndef INITUTILS_H
#define INITUTILS_H

#include "PlottingUtils.hpp"
#include <TROOT.h>
#include <TSystem.h>

class InitUtils {
public:
  static void SetROOTPreferences(Bool_t setupPlotting = kTRUE);
  static Bool_t ConvertWavedumpBinToROOT();
  static Bool_t ConvertCoMPASSBinToROOT();
  static Bool_t ConvertCoMPASSCSVToROOT();
};

#endif
