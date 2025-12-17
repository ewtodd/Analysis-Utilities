#include "InitUtils.hpp"

void InitUtils::SetROOTPreferences(Bool_t setupPlotting) {
  if (setupPlotting) {
    PlottingUtils::SetStylePreferences();
  }
  gROOT->ForceStyle(kTRUE);
  gROOT->SetBatch(kTRUE);
  if (gSystem->AccessPathName("plots")) {
    gSystem->mkdir("plots", kTRUE);
  }
  if (gSystem->AccessPathName("root_files")) {
    gSystem->mkdir("root_files", kTRUE);
  }
}
