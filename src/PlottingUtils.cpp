#include "PlottingUtils.hpp"

PlotSaveFormat PlottingUtils::save_format_ = PlotSaveFormat::kPNG;
Bool_t PlottingUtils::preferences_set_ = kFALSE;
Width_t PlottingUtils::line_width_ = 2;

void PlottingUtils::WarnIfNotConfigured(const TString method_name) {
  if (!preferences_set_)
    std::cout << "WARNING: " << method_name
              << " called before SetStylePreferences()." << std::endl;
}

void PlottingUtils::SetStylePreferences(PlotSaveFormat save_format) {
  save_format_ = save_format;
  preferences_set_ = kTRUE;
  line_width_ = save_format_ == PlotSaveFormat::kPNG ? 2 : 1;

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
  gStyle->SetHistLineWidth(line_width_);
  gStyle->SetLineWidth(line_width_);
  gStyle->SetPadGridX(1);
  gStyle->SetPadGridY(1);
  gStyle->SetGridStyle(3);
  gStyle->SetGridWidth(line_width_);
  gStyle->SetGridColor(kGray);
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
}

void PlottingUtils::ConfigureGraph(TGraph *graph, Int_t color,
                                   const TString title) {
  WarnIfNotConfigured("ConfigureGraph");
  graph->SetLineColor(color);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetLabelSize(0.06);
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetLineWidth(line_width_);
}

void PlottingUtils::ConfigureGraph(TGraphErrors *graph, Int_t color,
                                   const TString title) {
  WarnIfNotConfigured("ConfigureGraph");
  graph->SetLineColor(color);
  graph->SetTitle(title);
  graph->GetXaxis()->SetTitleSize(0.06);
  graph->GetYaxis()->SetTitleSize(0.06);
  graph->GetXaxis()->SetLabelSize(0.06);
  graph->GetYaxis()->SetLabelSize(0.06);
  graph->GetXaxis()->SetTitleOffset(1.2);
  graph->GetYaxis()->SetTitleOffset(1.2);
  graph->GetXaxis()->SetNdivisions(506);
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(1.2);
  graph->SetMarkerColor(color);
  graph->SetLineWidth(line_width_);
}

void PlottingUtils::ConfigureHistogram(TH1 *hist, Int_t color,
                                       const TString title) {
  if (!hist)
    return;
  WarnIfNotConfigured("ConfigureHistogram");

  hist->SetLineColor(color);
  hist->SetTitle(title);
  hist->SetLineWidth(line_width_);
  hist->SetFillStyle(0);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->GetXaxis()->SetNoExponent(kTRUE);
  hist->GetYaxis()->SetNdivisions(50109);
  hist->GetXaxis()->SetNdivisions(505);
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
  WarnIfNotConfigured("Configure2DHistogram");

  hist->SetTitle(title);
  hist->GetYaxis()->SetMoreLogLabels(kFALSE);
  hist->GetYaxis()->SetNoExponent(kFALSE);
  hist->GetXaxis()->SetTitleSize(0.06);
  hist->GetYaxis()->SetTitleSize(0.06);
  hist->GetXaxis()->SetLabelSize(0.06);
  hist->GetYaxis()->SetLabelSize(0.06);
  hist->GetXaxis()->SetTitleOffset(1.2);
  hist->GetYaxis()->SetTitleOffset(1.2);
  hist->GetXaxis()->SetNdivisions(506);
  hist->GetYaxis()->SetNdivisions(506);

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
  hist->Draw("HIST");
}

void PlottingUtils::ConfigureAndDraw2DHistogram(TH2 *hist, TCanvas *canvas,
                                                const TString title) {
  if (!hist)
    return;
  if (!canvas)
    return;

  Configure2DHistogram(hist, canvas, title);
  hist->Draw("COLZ");
}

TCanvas *PlottingUtils::GetConfiguredCanvas(Bool_t logy) {
  WarnIfNotConfigured("GetConfiguredCanvas");
  TCanvas *canvas = new TCanvas(GetRandomName(), "", 1200, 800);

  canvas->SetGridx(1);
  canvas->SetGridy(1);
  canvas->SetLogy(logy);

  canvas->SetTicks(1, 1);
  gPad->SetTicks(1, 1);

  return canvas;
}

void PlottingUtils::SaveFigure(TCanvas *canvas, TString output_name,
                               PlotSaveOptions save_options) {
  WarnIfNotConfigured("SaveFigure");
  canvas->SetLogy(kFALSE);
  canvas->Modified();
  canvas->Update();

  TString extension = (save_format_ == PlotSaveFormat::kPNG) ? ".png" : ".pdf";
  TString output_filename = output_name + extension;

  if (save_options != PlotSaveOptions::kLOG)
    canvas->Print("plots/" + output_filename);

  if (save_options != PlotSaveOptions::kLINEAR) {
    canvas->SetLogy(kTRUE);
    canvas->Modified();
    canvas->Update();
    canvas->Print("plots/log_" + output_filename);

    canvas->SetLogy(kFALSE);
    canvas->Modified();
    canvas->Update();
  }
}
std::vector<Int_t> PlottingUtils::GetDefaultColors() {
  return {kRed + 1,   kBlue + 1,   kGreen + 2,  kOrange + 1,  kMagenta + 1,
          kCyan + 2,  kViolet + 1, kSpring - 1, kPink + 1,    kTeal + 2,
          kAzure + 2, kYellow + 1, kOrange - 3, kMagenta - 3, kCyan - 6,
          kRed - 4,   kBlue - 4,   kGreen - 6,  kViolet - 4,  kSpring + 5,
          kPink - 3,  kTeal - 5,   kAzure - 3,  kOrange + 7};
}

TLegend *PlottingUtils::AddLegend(Double_t x1, Double_t x2, Double_t y1,
                                  Double_t y2) {
  TLegend *leg = new TLegend(x1, y1, x2, y2);
  leg->SetBorderSize(1);
  leg->SetFillColor(kWhite);
  leg->SetTextSize(0.05);
  leg->SetTextFont(132);
  leg->Draw();

  return leg;
}

TText *PlottingUtils::AddSubplotLabel(const TString label, Double_t x,
                                      Double_t y) {
  TText *text = new TText(x, y, label);
  text->SetNDC();
  text->SetTextSize(0.06);
  text->SetTextAlign(33);
  text->Draw();

  return text;
}

TString PlottingUtils::GetRandomName() {
  static TRandom3 generator(0);
  Double_t number = generator.Rndm();
  TString name = Form("name%.7f", number);
  return name;
}
