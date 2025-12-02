#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <algorithm>

void PlottingUtils::SetROOTStyle() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadTopMargin(0.12);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetTitleSize(0.06, "XY");
  gStyle->SetLabelSize(0.06, "XY");
  gStyle->SetLegendFont(132);
  gStyle->SetTitleOffset(1.2, "X");
  gStyle->SetTitleOffset(1.2, "Y");
  gStyle->SetTextFont(42);
  gStyle->SetHistLineWidth(1);
  gStyle->SetLineWidth(1);

  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(1);
  gStyle->SetGridColor(kGray);

  gStyle->SetPadTickX(2);
  gStyle->SetPadTickY(2);

  gROOT->ForceStyle(kTRUE);
}

std::string PlottingUtils::CleanSourceName(const std::string &source_name) {
  std::string clean_name = source_name;

  // Replace problematic characters with single underscores
  std::replace(clean_name.begin(), clean_name.end(), ' ', '_');
  std::replace(clean_name.begin(), clean_name.end(), '-', '_');
  std::replace(clean_name.begin(), clean_name.end(), '&', '_');

  // Remove consecutive underscores
  std::string::iterator new_end =
      std::unique(clean_name.begin(), clean_name.end(),
                  [](char a, char b) { return a == '_' && b == '_'; });
  clean_name.erase(new_end, clean_name.end());

  return clean_name;
}

void PlottingUtils::ConfigureHistogram(TH1 *hist, Int_t color,
                                       const std::string &title) {
  if (!hist)
    return;

  hist->SetLineColor(color);
  hist->SetTitle(title.c_str());
  hist->SetFillColorAlpha(color, 0.2);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);

  hist->SetMinimum(10);
  hist->GetYaxis()->SetNdivisions(50109);
  hist->GetXaxis()->SetNdivisions(506);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.2);
}

void PlottingUtils::ConfigureCanvas(TCanvas *canvas, Bool_t logy) {
  if (!canvas)
    return;

  canvas->SetGridx(1);
  canvas->SetGridy(1);
  canvas->SetLogy(logy);

  canvas->SetTicks(1, 1);
  gPad->SetTicks(1, 1);
}

std::vector<Int_t> PlottingUtils::GetDefaultColors() {
  return {kRed + 1, kBlue + 1, kGreen + 2, kMagenta,
          kOrange,  kCyan + 1, kYellow + 2};
}

Int_t PlottingUtils::GetSourceColor(Int_t source_id) {
  auto colors = GetDefaultColors();
  return colors[source_id % colors.size()];
}

TLegend *PlottingUtils::CreateLegend(Double_t x1, Double_t y1, Double_t x2,
                                     Double_t y2) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.05);
  leg->SetTextFont(132);
  return leg;
}

void PlottingUtils::AddSubplotLabel(const std::string &label, Double_t x,
                                    Double_t y) {
  TText *latex = new TText(x, y, label.c_str());
  latex->SetNDC();
  latex->SetTextSize(0.06);
  latex->SetTextAlign(33);
  latex->Draw();
}
