#include "PlottingUtils.hpp"
#include <TAttLine.h>
#include <TAttMarker.h>
#include <TAxis.h>
#include <TList.h>
#include <TText.h>

PlotSaveFormat PlottingUtils::save_format_ = PlotSaveFormat::kPNG;
Bool_t PlottingUtils::preferences_set_ = kFALSE;
Width_t PlottingUtils::line_width_ = 2;
TString PlottingUtils::plots_base_dir_ = "plots";

void PlottingUtils::SetPlotsBaseDir(const TString &dir) {
  TString d = dir;
  while (d.Length() > 0 && d[d.Length() - 1] == '/')
    d.Chop();
  plots_base_dir_ = d;
}

TString PlottingUtils::GetPlotsBaseDir() { return plots_base_dir_; }

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

  SetTurboPalette();
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
  hist->GetYaxis()->SetTitleOffset(1);
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

TCanvas *PlottingUtils::GetConfiguredCanvasWithSideLegend(TPad *&plot,
                                                          TPad *&legend_pad,
                                                          Int_t legend_width_px,
                                                          Int_t plot_height_px,
                                                          Bool_t logy) {
  WarnIfNotConfigured("GetConfiguredCanvasWithSideLegend");
  const Int_t plot_width_px = 1200;
  const Int_t canvas_width_px = plot_width_px + TMath::Max(0, legend_width_px);
  TCanvas *canvas =
      new TCanvas(GetRandomName(), "", canvas_width_px, plot_height_px);
  canvas->SetTicks(1, 1);
  const Double_t split = Double_t(plot_width_px) / Double_t(canvas_width_px);
  plot = new TPad(GetRandomName(), "", 0.0, 0.0, split, 1.0);
  plot->SetGridx(1);
  plot->SetGridy(1);
  plot->SetTicks(1, 1);
  plot->SetLogy(logy);
  const Double_t frame_right = split * (1.0 - plot->GetRightMargin());
  legend_pad = new TPad(GetRandomName(), "", frame_right, 0.0, 1.0, 1.0);
  legend_pad->SetFillStyle(4000);
  legend_pad->SetMargin(0.0, 0.0, 0.0, 0.0);
  plot->Draw();
  legend_pad->Draw();
  plot->cd();
  return canvas;
}

void PlottingUtils::SplitPadForPanel(TVirtualPad *pad, Int_t panel_height_px,
                                     TPad *&upper, TPad *&lower) {
  upper = nullptr;
  lower = nullptr;
  if (!pad || !pad->GetCanvas())
    return;
  const Double_t pad_height_px = pad->GetCanvas()->GetWh() * pad->GetAbsHNDC();
  if (!(pad_height_px > 0.0) || panel_height_px <= 0 ||
      panel_height_px >= pad_height_px)
    return;
  const Double_t upper_height_px = pad_height_px - panel_height_px;
  const Double_t split = Double_t(panel_height_px) / pad_height_px;
  const Double_t lower_top =
      (panel_height_px + pad->GetBottomMargin() * upper_height_px) /
      pad_height_px;
  TVirtualPad *previous = gPad;
  pad->cd();
  upper = new TPad(GetRandomName(), "", 0.0, split, 1.0, 1.0);
  upper->SetLeftMargin(pad->GetLeftMargin());
  upper->SetRightMargin(pad->GetRightMargin());
  upper->SetTopMargin(pad->GetTopMargin());
  upper->SetBottomMargin(pad->GetBottomMargin());
  upper->SetGridx(pad->GetGridx());
  upper->SetGridy(pad->GetGridy());
  upper->SetTicks(pad->GetTickx(), pad->GetTicky());
  upper->SetLogy(pad->GetLogy());
  lower = new TPad(GetRandomName(), "", 0.0, 0.0, 1.0, lower_top);
  lower->SetFillStyle(4000);
  lower->SetLeftMargin(pad->GetLeftMargin());
  lower->SetRightMargin(pad->GetRightMargin());
  lower->SetTopMargin(0.03);
  lower->SetBottomMargin(pad->GetBottomMargin());
  lower->SetGridx(pad->GetGridx());
  lower->SetGridy(pad->GetGridy());
  lower->SetTicks(pad->GetTickx(), pad->GetTicky());
  upper->Draw();
  lower->Draw();
  if (previous)
    previous->cd();
}

void PlottingUtils::SaveFigure(TCanvas *canvas, TString output_name,
                               TString output_subdirectory,
                               PlotSaveOptions save_options) {
  WarnIfNotConfigured("SaveFigure");
  canvas->SetLogy(kFALSE);
  canvas->Modified();
  canvas->Update();

  TString extension = (save_format_ == PlotSaveFormat::kPNG) ? ".png" : ".pdf";
  TString output_filename = output_name + extension;

  TString base = GetPlotsBaseDir();
  TString full_dir =
      output_subdirectory == "" ? base : base + "/" + output_subdirectory;

  if (gSystem->AccessPathName(full_dir)) {
    gSystem->mkdir(full_dir, kTRUE);
  }

  TString prefix = full_dir + "/";

  if (save_options != PlotSaveOptions::kLOG)
    canvas->Print(prefix + output_filename);

  if (save_options != PlotSaveOptions::kLINEAR) {
    TList *primitives = canvas->GetListOfPrimitives();
    Int_t size = primitives->GetSize();

    for (Int_t i = 0; i < size; i++) {
      TObject *object = primitives->At(i);
      if (object->InheritsFrom(TH2::Class())) {
        std::cout << std::endl;
        std::cerr << "ERROR: Used PlotSaveOptions::kLOG for 2D histogram."
                  << std::endl;
        std::cout
            << "This option is exclusive to 1D histograms/graphs because it refers to the y axis."
            << std::endl;
        std::cout << "Use PlotSaveOptions::kLINEAR." << std::endl;
        std::cout
            << "The z axis is already log by default if you used PlottingUtils::Configure2DHistogram()."
            << std::endl;
        std::cout << "Plot " << prefix + "log_" + output_filename
                  << " was not saved." << std::endl;
        std::exit(1);
      };
    };

    canvas->SetLogy(kTRUE);
    canvas->Modified();
    canvas->Update();
    canvas->Print(prefix + "log_" + output_filename);

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
  leg->SetTextSize(30);
  leg->SetTextFont(43);
  leg->Draw();

  return leg;
}

Double_t PlottingUtils::TextWidthPx(const TString &text, Double_t size_px,
                                    Int_t font) {
  if (text.Length() == 0)
    return 0.0;
  TVirtualPad *previous = gPad;
  const Int_t scratch_px = 1000;
  TCanvas *scratch = new TCanvas(GetRandomName(), "", scratch_px, scratch_px);
  TLatex latex(0.1, 0.5, text);
  latex.SetNDC();
  latex.SetTextFont(font);
  latex.SetTextSize(size_px);
  const Double_t width_px = latex.GetXsize() * scratch_px;
  delete scratch;
  if (previous)
    previous->cd();
  return width_px;
}

Double_t PlottingUtils::LegendWidthPx(const std::vector<TString> &labels,
                                      const TString &header, Double_t size_px,
                                      Double_t margin, Double_t padding_px) {
  Double_t width_px = TextWidthPx(header, size_px);
  for (std::size_t i = 0; i < labels.size(); i++)
    width_px = TMath::Max(width_px, TextWidthPx(labels[i], size_px) /
                                        TMath::Max(0.05, 1.0 - margin));
  return width_px + 2.0 * padding_px;
}

Double_t PlottingUtils::FigureScale(TVirtualPad *pad,
                                    Int_t reference_width_px) {
  if (!pad || !pad->GetCanvas() || reference_width_px <= 0)
    return 1.0;
  const Double_t pad_width_px = pad->GetCanvas()->GetWw() * pad->GetAbsWNDC();
  return pad_width_px / Double_t(reference_width_px);
}

void PlottingUtils::ScaleFigure(TVirtualPad *pad, Int_t reference_width_px,
                                Int_t reference_height_px) {
  WarnIfNotConfigured("ScaleFigure");
  if (!pad)
    return;
  ScalePad(pad, FigureScale(pad, reference_width_px), reference_width_px,
           reference_height_px);
}

void PlottingUtils::ScalePad(TVirtualPad *pad, Double_t scale,
                             Int_t reference_width_px,
                             Int_t reference_height_px) {
  if (!pad || !pad->GetCanvas())
    return;
  const Double_t pad_width_px = pad->GetCanvas()->GetWw() * pad->GetAbsWNDC();
  const Double_t pad_height_px = pad->GetCanvas()->GetWh() * pad->GetAbsHNDC();
  if (!(pad_width_px > 0.0) || !(pad_height_px > 0.0))
    return;
  const Double_t x_factor = reference_width_px * scale / pad_width_px;
  const Double_t y_factor = reference_height_px * scale / pad_height_px;
  pad->SetLeftMargin(TMath::Min(0.9, pad->GetLeftMargin() * x_factor));
  pad->SetRightMargin(TMath::Min(0.9, pad->GetRightMargin() * x_factor));
  pad->SetTopMargin(TMath::Min(0.9, pad->GetTopMargin() * y_factor));
  pad->SetBottomMargin(TMath::Min(0.9, pad->GetBottomMargin() * y_factor));

  TList *primitives = pad->GetListOfPrimitives();
  if (!primitives)
    return;
  TIter next(primitives);
  for (TObject *object = next(); object; object = next()) {
    if (object->InheritsFrom(TVirtualPad::Class())) {
      ScalePad(static_cast<TVirtualPad *>(object), scale, reference_width_px,
               reference_height_px);
      continue;
    }
    if (object->InheritsFrom(TLegend::Class())) {
      ScaleLegend(static_cast<TLegend *>(object), scale);
      continue;
    }
    if (object->InheritsFrom(TText::Class())) {
      ScaleText(static_cast<TText *>(object), pad, scale, reference_height_px);
      continue;
    }
    if (object->InheritsFrom(TH1::Class()))
      ScaleAxisFontsToPad(static_cast<TH1 *>(object), pad, scale,
                          reference_width_px, reference_height_px);
    ScaleLine(dynamic_cast<TAttLine *>(object), scale);
    ScaleMarker(dynamic_cast<TAttMarker *>(object), scale);
  }
}

void PlottingUtils::ScaleAxisFontsToPad(TH1 *frame, TVirtualPad *pad,
                                        Double_t scale,
                                        Int_t reference_width_px,
                                        Int_t reference_height_px) {
  WarnIfNotConfigured("ScaleAxisFontsToPad");
  if (!frame || !pad || !pad->GetCanvas() || reference_width_px <= 0 ||
      reference_height_px <= 0)
    return;
  const Double_t pad_width_px = pad->GetCanvas()->GetWw() * pad->GetAbsWNDC();
  const Double_t pad_height_px = pad->GetCanvas()->GetWh() * pad->GetAbsHNDC();
  const Double_t pad_px = TMath::Min(pad_width_px, pad_height_px);
  if (!(pad_px > 0.0))
    return;
  const Double_t size = 0.06 * reference_height_px * scale / pad_px;
  const Double_t reference_aspect =
      Double_t(TMath::Min(reference_width_px, reference_height_px)) /
      Double_t(reference_width_px);
  TAxis *axes[2] = {frame->GetXaxis(), frame->GetYaxis()};
  for (Int_t i = 0; i < 2; i++) {
    if (axes[i]->GetTitleSize() > 0.0)
      axes[i]->SetTitleSize(size);
    if (axes[i]->GetLabelSize() > 0.0)
      axes[i]->SetLabelSize(size);
  }
  frame->GetXaxis()->SetTitleOffset(1.2 * pad_px / pad_height_px);
  frame->GetYaxis()->SetTitleOffset(1.2 * (pad_px / pad_width_px) /
                                    reference_aspect);
}

void PlottingUtils::ScaleLegend(TLegend *legend, Double_t scale) {
  if (!legend)
    return;
  legend->SetTextSize(30.0 * scale);
}

void PlottingUtils::ScaleText(TText *text, TVirtualPad *pad, Double_t scale,
                              Int_t reference_height_px) {
  if (!text || !pad || !pad->GetCanvas())
    return;
  const Int_t precision = text->GetTextFont() % 10;
  if (precision == 3) {
    text->SetTextSize(text->GetTextSize() * scale);
    return;
  }
  const Double_t pad_width_px = pad->GetCanvas()->GetWw() * pad->GetAbsWNDC();
  const Double_t pad_height_px = pad->GetCanvas()->GetWh() * pad->GetAbsHNDC();
  const Double_t pad_px = TMath::Min(pad_width_px, pad_height_px);
  if (!(pad_px > 0.0))
    return;
  text->SetTextSize(text->GetTextSize() * reference_height_px * scale / pad_px);
}

TLatex *PlottingUtils::DrawTitle(TVirtualPad *pad, const TString &title,
                                 Double_t scale) {
  if (!pad || title.Length() == 0)
    return nullptr;
  TVirtualPad *previous = gPad;
  pad->cd();
  TLatex *text = new TLatex(0.5, 1.0 - 0.5 * pad->GetTopMargin(), title);
  text->SetNDC();
  text->SetTextFont(43);
  text->SetTextSize(36.0 * scale);
  text->SetTextAlign(22);
  text->Draw();
  if (previous)
    previous->cd();
  return text;
}

void PlottingUtils::ScaleLine(TAttLine *line, Double_t scale) {
  if (!line || line->GetLineWidth() == 0)
    return;
  line->SetLineWidth(
      Width_t(TMath::Max(1, TMath::Nint(line->GetLineWidth() * scale))));
}

void PlottingUtils::ScaleMarker(TAttMarker *marker, Double_t scale) {
  if (!marker)
    return;
  marker->SetMarkerSize(marker->GetMarkerSize() * scale);
}

TLatex *PlottingUtils::AddText(const TString label, Double_t x, Double_t y,
                               Double_t angle) {
  TLatex *text = new TLatex(x, y, label);
  text->SetNDC();
  text->SetTextSize(30);
  text->SetTextAlign(33);
  text->SetTextFont(43);
  text->SetTextAngle(angle);
  text->Draw();

  return text;
}

TString PlottingUtils::GetRandomName() {
  static TRandom3 generator(0);
  Double_t number = generator.Rndm();
  TString name = Form("name%.7f", number);
  return name;
}

void PlottingUtils::PlotFitWithResiduals(
    TH1 *hist, TGraph *total_graph,
    const std::vector<TGraph *> &component_graphs, Float_t fit_range_low,
    Float_t fit_range_high, const TString &output_name,
    const TString &output_subdirectory, const TString &label, Bool_t logy) {
  TCanvas *canvas = GetConfiguredCanvas(kFALSE);

  TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
  TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.3);
  pad1->SetBottomMargin(0.04);
  pad1->SetGridx(1);
  pad1->SetGridy(1);
  pad1->SetTopMargin(0.12);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.35);
  pad2->SetGridx(1);
  pad2->SetGridy(1);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();

  Float_t min_hist_value = 0.9 * fit_range_low;
  Float_t max_hist_value = 1.1 * fit_range_high;

  hist->GetXaxis()->SetRangeUser(min_hist_value, max_hist_value);
  hist->GetXaxis()->SetLabelSize(0);
  hist->GetXaxis()->SetTitleSize(0);
  hist->SetLineColor(kViolet);
  hist->GetYaxis()->SetTitleOffset(1);
  hist->SetLineWidth(line_width_);
  hist->Draw();

  pad1->SetTickx(0);
  if (total_graph)
    total_graph->Draw("L same");
  for (std::size_t i = 0; i < component_graphs.size(); i++) {
    if (component_graphs[i])
      component_graphs[i]->Draw("L same");
  }

  pad2->cd();

  Int_t nbins = hist->GetNbinsX();
  TGraph *residuals = new TGraph();
  Int_t point_counter = 0;
  for (Int_t i = 1; i <= nbins; i++) {
    Double_t x = hist->GetBinCenter(i);
    if (x < fit_range_low || x > fit_range_high)
      continue;
    Double_t data = hist->GetBinContent(i);
    Double_t fit_val = total_graph ? total_graph->Eval(x) : 0.0;
    Double_t error = hist->GetBinError(i);

    if (error > 0 && data > 0) {
      Double_t pull = (data - fit_val) / error;
      residuals->SetPoint(point_counter, x, pull);
      point_counter++;
    }
  }

  residuals->SetMarkerStyle(20);
  residuals->SetMarkerSize(0.8);
  residuals->SetMarkerColor(kAzure);
  residuals->SetLineColor(kAzure);
  residuals->SetTitle("");
  Double_t actual_min =
      hist->GetXaxis()->GetBinLowEdge(hist->GetXaxis()->GetFirst());
  Double_t actual_max =
      hist->GetXaxis()->GetBinUpEdge(hist->GetXaxis()->GetLast());
  residuals->GetXaxis()->SetLimits(actual_min, actual_max);
  residuals->GetYaxis()->SetTitle("#delta/#sigma");
  residuals->GetXaxis()->SetTitle(hist->GetXaxis()->GetTitle());
  residuals->GetXaxis()->SetTitleSize(0.13);
  residuals->GetYaxis()->SetTitleSize(0.13);
  residuals->GetXaxis()->SetLabelSize(0.12);
  residuals->GetYaxis()->SetLabelSize(0.12);
  residuals->GetXaxis()->SetTitleOffset(1.0);
  residuals->GetYaxis()->SetTitleOffset(0.3);
  residuals->GetYaxis()->SetNdivisions(505);
  residuals->GetXaxis()->SetNdivisions(510);
  residuals->GetYaxis()->CenterTitle(kTRUE);
  residuals->GetYaxis()->SetRangeUser(-5.5, 5.5);
  residuals->Draw("AP");

  TF1 *zero_line = new TF1("zero_line", "0", actual_min, actual_max);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineStyle(2);
  zero_line->SetLineWidth(line_width_);
  zero_line->Draw("same");

  TF1 *plus3_line = new TF1("plus3_line", "3", actual_min, actual_max);
  plus3_line->SetLineColor(kGray + 2);
  plus3_line->SetLineStyle(3);
  plus3_line->SetLineWidth(line_width_);
  plus3_line->Draw("same");

  TF1 *minus3_line = new TF1("minus3_line", "-3", actual_min, actual_max);
  minus3_line->SetLineColor(kGray + 2);
  minus3_line->SetLineStyle(3);
  minus3_line->SetLineWidth(line_width_);
  minus3_line->Draw("same");

  pad1->cd();
  pad1->SetLogy(logy);
  if (label.Length() > 0) {
    AddText(label, 0.85, 0.85);
  }

  PlotSaveOptions save_opts =
      logy ? PlotSaveOptions::kLOG : PlotSaveOptions::kLINEAR;
  SaveFigure(canvas, output_name, output_subdirectory, save_opts);

  PlotPullHistogram(residuals, output_name, output_subdirectory);

  gROOT->GetListOfCanvases()->Remove(canvas);
  canvas->SetBatch(kTRUE);
  delete canvas;
  delete residuals;
  delete zero_line;
  delete plus3_line;
  delete minus3_line;
}

void PlottingUtils::PlotPullHistogram(TGraph *residuals,
                                      const TString &output_name,
                                      const TString &output_subdirectory) {
  Int_t npoints = residuals->GetN();
  if (npoints == 0)
    return;

  TH1D *pull_hist =
      new TH1D("pull_hist", ";#delta/#sigma;Counts", 82, -5.5, 5.5);

  Double_t *y = residuals->GetY();
  for (Int_t i = 0; i < npoints; i++) {
    pull_hist->Fill(y[i]);
  }

  TH1D *gauss_ref = new TH1D("gauss_ref", "", 82, -5.5, 5.5);
  for (Int_t i = 1; i <= gauss_ref->GetNbinsX(); i++) {
    Double_t lo = gauss_ref->GetBinLowEdge(i);
    Double_t hi = lo + gauss_ref->GetBinWidth(i);
    gauss_ref->SetBinContent(i, npoints * (TMath::Freq(hi) - TMath::Freq(lo)));
  }

  Double_t ks_pvalue = pull_hist->KolmogorovTest(gauss_ref);

  TCanvas *hist_canvas = GetConfiguredCanvas(kFALSE);
  ConfigureAndDrawHistogram(pull_hist, kAzure);
  AddText(TString::Format("KS p = %.3f", ks_pvalue), 0.85, 0.85);

  SaveFigure(hist_canvas, "residuals_" + output_name,
             output_subdirectory + "/residual_hists", PlotSaveOptions::kLINEAR);

  delete hist_canvas;
  delete gauss_ref;
  delete pull_hist;
}

void PlottingUtils::SetTurboPalette() {
  Int_t N_COLORS = 255;

  std::vector<Double_t> stops(N_COLORS);
  std::vector<Double_t> reds(N_COLORS);
  std::vector<Double_t> greens(N_COLORS);
  std::vector<Double_t> blues(N_COLORS);

  for (Int_t i = 0; i < N_COLORS; ++i) {
    Double_t x = Double_t(i) / Double_t(N_COLORS - 1);

    Double_t r = 0.13572138 + 4.61539260 * x - 42.66032258 * x * x +
                 132.13108234 * x * x * x - 152.94239396 * x * x * x * x +
                 59.28637943 * x * x * x * x * x;

    Double_t g = 0.09140261 + 2.19418839 * x + 4.84296658 * x * x -
                 14.18503333 * x * x * x + 4.27729857 * x * x * x * x +
                 2.82956604 * x * x * x * x * x;

    Double_t b = 0.10667330 + 12.64194608 * x - 60.58204836 * x * x +
                 110.36276771 * x * x * x - 89.90310912 * x * x * x * x +
                 27.34824973 * x * x * x * x * x;

    stops[i] = x;
    reds[i] = std::max(0.0, std::min(1.0, r));
    greens[i] = std::max(0.0, std::min(1.0, g));
    blues[i] = std::max(0.0, std::min(1.0, b));
  }

  Int_t firstColor =
      TColor::CreateGradientColorTable(N_COLORS, stops.data(), reds.data(),
                                       greens.data(), blues.data(), N_COLORS);

  std::vector<Int_t> palette(N_COLORS);
  for (Int_t i = 0; i < N_COLORS; ++i) {
    palette[i] = firstColor + i;
  }

  gStyle->SetNumberContours(N_COLORS);
  gStyle->SetPalette(N_COLORS, palette.data());
}
