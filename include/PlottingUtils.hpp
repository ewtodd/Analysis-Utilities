#ifndef PLOTTINGUTILS_H
#define PLOTTINGUTILS_H

#include <TCanvas.h>
#include <TF1.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <vector>

class PlottingUtils {
public:
  static void SetROOTPreferences();
  static void ConfigureAndDrawGraph(TGraph *graph, Int_t color,
                                    const TString title);
  static void ConfigureAndDrawHistogram(TH1 *hist, Int_t color,
                                        const TString title = "");
  static void ConfigureAndDraw2DHistogram(TH2 *hist, TCanvas *canvas,
                                          Int_t color,
                                          const TString title = "");

  static void ConfigureGraph(TGraph *graph, Int_t color, const TString title);
  static void ConfigureHistogram(TH1 *hist, Int_t color,
                                 const TString title = "");
  static void Configure2DHistogram(TH2 *hist, TCanvas *canvas, Int_t color,
                                   const TString title = "");

  static void ConfigureCanvas(TCanvas *canvas, Bool_t logy = kFALSE);
  static void SaveFigure(TCanvas *canvas, TString output_name,
                         Bool_t log = kTRUE);

  static void AddLegend(Double_t x1 = 0.7, Double_t y1 = 0.7, Double_t x2 = 0.9,
                        Double_t y2 = 0.9);
  static void AddSubplotLabel(const TString label, Double_t x = 0.9,
                              Double_t y = 0.85);

  static std::vector<Int_t> GetDefaultColors();
  static Int_t GetSourceColor(Int_t source_id);
};

#endif
