#include "PlottingUtils.hpp"

void PlottingUtils::SetROOTPreferences() {
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.1);
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
  gROOT->SetBatch(kTRUE);
  if (gSystem->AccessPathName("plots")) {
    gSystem->mkdir("plots", kTRUE);
  }
}

void PlottingUtils::ConfigureGraph(TGraph *graph, Int_t color,
                                   const TString title) {
  graph->SetLineColor(color);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetLabelSize(0.06);
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.2);
}

void PlottingUtils::ConfigureHistogram(TH1 *hist, Int_t color,
                                       const TString title) {
  if (!hist)
    return;

  hist->SetLineColor(color);
  hist->SetTitle(title);
  hist->SetFillColorAlpha(color, 0.2);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->GetXaxis()->SetNoExponent(kTRUE);

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

void PlottingUtils::Configure2DHistogram(TH2 *hist, TCanvas *canvas,
                                         const TString title) {
  if (!hist)
    return;
  if (!canvas)
    return;

  hist->SetTitle(title);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->SetMinimum(1);
  canvas->SetLogz(kTRUE);
  canvas->SetRightMargin(0.15);
}

void PlottingUtils::ConfigureAndDrawGraph(TGraph *graph, Int_t color,
                                          const TString title) {
  if (!graph)
    return;

  ConfigureGraph(graph, color, title);
  graph->Draw();
}

void PlottingUtils::ConfigureAndDrawHistogram(TH1 *hist, Int_t color,
                                              const TString title) {
  if (!hist)
    return;

  ConfigureHistogram(hist, color, title);
  hist->Draw();
}

void PlottingUtils::ConfigureAndDraw2DHistogram(TH2 *hist, TCanvas *canvas,
                                                const TString title) {
  if (!hist)
    return;
  if (!canvas)
    return;

  ConfigureAndDraw2DHistogram(hist, canvas, title);
  hist->Draw();
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

void PlottingUtils::SaveFigure(TCanvas *canvas, TString output_filename,
                               Bool_t log) {
  canvas->Print("plots/" + output_filename);
  if (log) {
    canvas->SetLogy(kTRUE);
    canvas->Print("plots/log_" + output_filename);
  }
}

std::vector<Int_t> PlottingUtils::GetDefaultColors() {
  return {kRed + 1, kBlue + 1, kGreen + 2, kMagenta,
          kOrange,  kCyan + 1, kYellow + 2};
}

Int_t PlottingUtils::GetSourceColor(Int_t source_id) {
  auto colors = GetDefaultColors();
  return colors[source_id % colors.size()];
}

void PlottingUtils::AddLegend(Double_t x1, Double_t y1, Double_t x2,
                              Double_t y2) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.05);
  leg->SetTextFont(132);
  leg->Draw();
}

void PlottingUtils::AddSubplotLabel(const TString label, Double_t x,
                                    Double_t y) {
  TText *latex = new TText(x, y, label);
  latex->SetNDC();
  latex->SetTextSize(0.06);
  latex->SetTextAlign(33);
  latex->Draw();
}
