#include "PlottingUtils.h"
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
#include <iostream>

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

void PlottingUtils::PlotIndividualSpectra(
    const std::map<Int_t, TH1F *> &spectra,
    const std::map<Int_t, std::string> &source_names,
    const std::string &output_prefix) {
  for (const auto &pair : spectra) {
    Int_t src_id = pair.first;
    TH1F *hist = pair.second;
    if (!hist)
      continue;

    auto name_it = source_names.find(src_id);
    std::string source_name = (name_it != source_names.end())
                                  ? name_it->second
                                  : "Source_" + std::to_string(src_id);

    TCanvas *c = new TCanvas(("c_" + source_name).c_str(), source_name.c_str(),
                             1200, 800);
    ConfigureCanvas(c, kTRUE);

    Int_t color = GetSourceColor(src_id);
    ConfigureHistogram(hist, color);
    hist->Draw("HIST");

    AddSubplotLabel("(" + std::string(1, 'a' + src_id) + ") " + source_name,
                    0.9, 0.83);

    c->Modified();
    c->Update();

    std::string clean_name = CleanSourceName(source_name);
    std::string pdf_name = output_prefix + clean_name + ".pdf";
    c->SaveAs(pdf_name.c_str());

    delete c;
  }
}

void PlottingUtils::PlotCalibrationCurve(TGraphErrors *calibration_curve,
                                         TF1 *linear_function,
                                         const std::string &output_name) {

  if (!calibration_curve) {
    std::cout << "PlotCalibrationCurve: No calibration curve provided."
              << std::endl;
    return;
  }

  TCanvas *c = new TCanvas("c_cal_curve", "Calibration Curve", 1200, 800);
  ConfigureCanvas(c, kFALSE);

  calibration_curve->SetMarkerStyle(20);
  calibration_curve->SetMarkerSize(1.2);
  calibration_curve->SetMarkerColor(kBlack);
  calibration_curve->SetLineColor(kBlack);

  TF1 *quad_func = calibration_curve->GetFunction("calibration_function");
  if (quad_func) {
    Double_t x_min, x_max, y_dummy;
    calibration_curve->ComputeRange(x_min, y_dummy, x_max, y_dummy);
    quad_func->SetRange(
        0, x_max + 1000); // Start from 0 to show zero point behavior
    quad_func->SetLineColor(kRed);
    quad_func->SetLineStyle(1);
  }

  if (linear_function) {
    linear_function->SetLineColor(kBlue);
    linear_function->SetLineStyle(2);
  }

  // Set up axes - make sure to start from 0 to show origin
  calibration_curve->SetTitle("");
  calibration_curve->GetXaxis()->SetTitle("Pulse Integral [a.u.]");
  calibration_curve->GetYaxis()->SetTitle("Deposited Energy [keV]");
  calibration_curve->GetXaxis()->SetTitleSize(0.06);
  calibration_curve->GetXaxis()->SetNdivisions(506);
  calibration_curve->GetYaxis()->SetTitleSize(0.06);
  calibration_curve->GetXaxis()->SetLabelSize(0.06);
  calibration_curve->GetYaxis()->SetLabelSize(0.06);
  calibration_curve->GetXaxis()->SetTitleOffset(1.2);
  calibration_curve->GetYaxis()->SetTitleOffset(1.2);

  // Ensure the plot shows from origin
  Double_t x_min, x_max, y_min, y_max;
  calibration_curve->ComputeRange(x_min, y_min, x_max, y_max);
  calibration_curve->GetHistogram()->SetMaximum(y_max * 1.1);

  // Draw all components
  calibration_curve->Draw("APE");

  if (linear_function) {
    linear_function->Draw("SAME");
  }

  // Check if (0,0) is in the calibration points
  Double_t *x_vals = calibration_curve->GetX();
  Double_t *y_vals = calibration_curve->GetY();
  Int_t n_points = calibration_curve->GetN();

  for (Int_t i = 0; i < n_points; ++i) {
    if (TMath::Abs(x_vals[i]) < 1e-6 && TMath::Abs(y_vals[i]) < 1e-6) {
      break;
    }
  }

  // Create legend
  TLegend *leg = CreateLegend(0.2, 0.65, 0.6, 0.84);
  leg->AddEntry(calibration_curve, "Calibration points", "pe");
  leg->SetMargin(0.2);
  if (quad_func) {
    leg->AddEntry(quad_func, "Quadratic fit", "l");
  }
  if (linear_function) {
    leg->AddEntry(linear_function, "Linear fit (high energy)", "l");
  }
  leg->Draw();

  c->Modified();
  c->Update();

  std::string pdf_name = output_name + ".pdf";
  c->SaveAs(pdf_name.c_str());
  delete c;
}
void PlottingUtils::PlotFittedPeaks(const std::vector<CalibrationPeak> &peaks,
                                    const std::map<Int_t, TH1F *> &spectra,
                                    const std::string &output_prefix) {
  if (peaks.empty() || spectra.empty()) {
    std::cout << "Error: No peaks or spectra available for plotting!"
              << std::endl;
    return;
  }

  for (size_t i = 0; i < peaks.size(); ++i) {
    const CalibrationPeak &peak = peaks[i];

    // FIX: Use the peak's source_id to find correct spectrum
    auto spectrum_it = spectra.find(peak.source_id);
    if (spectrum_it == spectra.end()) {
      std::cout << "Error: No spectrum found for source_id " << peak.source_id
                << " for peak " << peak.isotope << std::endl;
      continue;
    }

    TH1F *hist = spectrum_it->second;
    std::string filename = output_prefix + "_" + peak.isotope + ".pdf";
    PlotSinglePeakFit(peak, hist, filename);
  }
}

void PlottingUtils::PlotSinglePeakFit(const CalibrationPeak &peak,
                                      TH1F *histogram,
                                      const std::string &output_filename) {
  if (gPad) {
    TList *prims = gPad->GetListOfPrimitives();
    TIter next(prims);
    TObject *obj;
    std::vector<TObject *> to_remove;

    while ((obj = next())) {
      if (obj->InheritsFrom("TLegend")) {
        to_remove.push_back(obj);
      }
    }

    for (auto *legend_obj : to_remove) {
      prims->Remove(legend_obj);
    }
  }
  TH1F *hist = histogram;
  ConfigureHistogram(hist, kBlue, "");

  // Create canvas for this peak
  TCanvas *canvas =
      new TCanvas(("canvas_" + peak.isotope).c_str(),
                  ("Peak Fit: " + peak.isotope).c_str(), 1200, 800);
  ConfigureCanvas(canvas, kTRUE);

  // Configure and draw the fine-binned histogram
  hist->SetTitle(("Peak Fit: " + peak.isotope + " (" +
                  std::to_string(peak.deposited_energy_kev) + " keV)")
                     .c_str());
  hist->Draw("HIST");

  if (peak.fit_successful) {
    // Create fit functions with power law background
    TF1 *complete_fit = nullptr;
    TF1 *gauss_only = nullptr;

    // Complete fit: Gaussian + linear background + power law background
    complete_fit = new TF1(
        "complete_fit",
        "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x + [5]*pow(x, -[6])",
        peak.fit_range_low, peak.fit_range_high);
    complete_fit->SetParameter(0, peak.fitted_amplitude);
    complete_fit->SetParameter(1, peak.fitted_position);
    complete_fit->SetParameter(2, peak.fitted_sigma);
    complete_fit->SetParameter(3, peak.fitted_bkg_const);
    complete_fit->SetParameter(4, peak.fitted_bkg_slope);
    complete_fit->SetParameter(5, peak.fitted_powerlaw_amp);
    complete_fit->SetParameter(6, peak.fitted_powerlaw_exp);

    // Gaussian component only
    gauss_only = new TF1("gauss_only", "[0]*exp(-0.5*((x-[1])/[2])^2)",
                         peak.fit_range_low, peak.fit_range_high);
    gauss_only->SetParameter(0, peak.fitted_amplitude);
    gauss_only->SetParameter(1, peak.fitted_position);
    gauss_only->SetParameter(2, peak.fitted_sigma);

    // Background component (linear + power law)
    TF1 *bkg_only = new TF1("bkg_only", "[0] + [1]*x + [2]*pow(x, -[3])",
                            peak.fit_range_low, peak.fit_range_high);
    bkg_only->SetParameter(0, peak.fitted_bkg_const);
    bkg_only->SetParameter(1, peak.fitted_bkg_slope);
    bkg_only->SetParameter(2, peak.fitted_powerlaw_amp);
    bkg_only->SetParameter(3, peak.fitted_powerlaw_exp);

    // Set line styles and colors
    complete_fit->SetLineColor(kRed);
    complete_fit->Draw("same");

    gauss_only->SetLineColor(kGray + 2);
    gauss_only->SetLineStyle(2);
    gauss_only->Draw("same");

    bkg_only->SetLineColor(kGreen + 2);
    bkg_only->SetLineStyle(2);
    bkg_only->Draw("same");

    // LEGEND - positioned at bottom right
    TLegend *legend = new TLegend(0.55, 0.75, 0.88, 0.85);
    legend->SetBorderSize(1);
    legend->SetFillColor(kWhite);
    legend->SetTextSize(0.04);
    legend->SetTextFont(42);

    legend->AddEntry(complete_fit, "Total fit", "l");
    legend->AddEntry(gauss_only, "Gaussian", "l");
    legend->AddEntry(bkg_only, "Background", "l");
    legend->Draw();

    // Text positioned at bottom left to avoid legend
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);

    Double_t y_pos = 0.45;
    Double_t y_step = 0.05;
    Double_t x_pos = 0.15;
    char buf[256];

    // Position text with proper formatting
    sprintf(buf, "Position: %.1f +/- %.1f", peak.fitted_position,
            peak.fitted_position_error);
    text->DrawLatex(x_pos, y_pos, buf);
    y_pos -= y_step;

    // FWHM text
    sprintf(buf, "FWHM: %.2f", peak.fitted_fwhm);
    text->DrawLatex(x_pos, y_pos, buf);
    y_pos -= y_step;

    // Energy text
    sprintf(buf, "Energy: %.1f keV", peak.deposited_energy_kev);
    text->DrawLatex(x_pos, y_pos, buf);
    y_pos -= y_step;

    // Background parameters
    sprintf(buf, "Bkg: %.1f + %.2e*x + %.1e*x^{-%.2f}", peak.fitted_bkg_const,
            peak.fitted_bkg_slope, peak.fitted_powerlaw_amp,
            peak.fitted_powerlaw_exp);
    text->DrawLatex(x_pos, y_pos, buf);

  } else {
    // Mark failed fits
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.05);
    text->SetTextColor(kRed);
    text->DrawLatex(0.15, 0.85, "FIT FAILED");
  }

  // Save standard scale version
  canvas->SaveAs(output_filename.c_str());

  delete canvas;
}
void PlottingUtils::PlotAverageWaveforms(
    const std::map<Int_t, std::vector<Float_t>> &waveforms,
    const std::map<Int_t, std::string> &source_names,
    const std::string &output_name) {
  if (waveforms.empty()) {
    std::cout << "PlotAverageWaveforms: No waveforms provided." << std::endl;
    return;
  }

  TCanvas *c = new TCanvas("c_avg_waveforms", "Average Waveforms", 1200, 800);
  ConfigureCanvas(c, kFALSE);

  Bool_t first = kTRUE;
  TLegend *leg = CreateLegend(0.6, 0.7, 0.9, 0.85);
  std::vector<TGraph *> graphs;

  for (const auto &pair : waveforms) {
    Int_t src_id = pair.first;
    const std::vector<Float_t> &waveform = pair.second;
    if (waveform.empty())
      continue;

    Int_t color = GetSourceColor(src_id);

    // Create TGraph
    Int_t n_points = waveform.size();
    Double_t *time_array = new Double_t[n_points];
    Double_t *amp_array = new Double_t[n_points];

    for (Int_t i = 0; i < n_points; ++i) {
      time_array[i] = 2.0 * i; // 2 ns per sample
      amp_array[i] = waveform[i];
    }

    TGraph *gr = new TGraph(n_points, time_array, amp_array);
    gr->SetLineColor(color);
    gr->SetLineWidth(2);
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0);
    gr->SetMarkerColor(color);

    if (first) {
      gr->SetTitle("");
      gr->GetXaxis()->SetTitle("Time [ns]");
      gr->GetYaxis()->SetTitle("Normalized Amplitude [a.u.]");
      gr->GetXaxis()->SetTitleSize(0.06);
      gr->GetYaxis()->SetTitleSize(0.06);
      gr->GetXaxis()->SetLabelSize(0.06);
      gr->GetYaxis()->SetLabelSize(0.06);
      gr->Draw("APL");
      first = kFALSE;
    } else {
      gr->Draw("PL SAME");
    }

    auto it = source_names.find(src_id);
    std::string label = (it != source_names.end())
                            ? it->second
                            : "Source " + std::to_string(src_id);

    if (label == "Am-241") {
      label += " (#alpha)";
    } else if (label == "Na-22") {
      label += " (#gamma)";
    } else if (label == "Cs-137") {
      label += " (#gamma)";
    }

    leg->AddEntry(gr, label.c_str(), "lp");
    graphs.push_back(gr);

    delete[] time_array;
    delete[] amp_array;
  }

  leg->Draw();

  c->Modified();
  c->Update();

  std::string pdf_name = output_name + ".pdf";
  c->SaveAs(pdf_name.c_str());

  // Log scale version
  c->SetLogy();
  if (!graphs.empty()) {
    graphs[0]->SetMinimum(5e-4);
    graphs[0]->SetMaximum(2);
  }
  c->Modified();
  c->Update();
  c->SaveAs((output_name + "_log.pdf").c_str());

  delete c;
}

void PlottingUtils::PlotAverageWaveforms(const std::vector<Float_t> &f_alpha,
                                         const std::vector<Float_t> &f_gamma,
                                         const std::string &output_name) {

  if (f_alpha.empty() || f_gamma.empty()) {
    std::cout << "PlotAverageWaveforms: No waveforms provided." << std::endl;
    return;
  }

  TCanvas *c = new TCanvas("c_avg_waveforms", "Average Waveforms", 1200, 800);
  ConfigureCanvas(c, kFALSE);

  TLegend *leg = CreateLegend(0.72, 0.7, 0.92, 0.85);
  std::vector<TGraph *> graphs;

  // Determine common length
  Int_t n_points = std::min(f_alpha.size(), f_gamma.size());
  Double_t *time_array = new Double_t[n_points];

  for (Int_t i = 0; i < n_points; ++i) {
    time_array[i] = 2.0 * i; // 2 ns per sample
  }

  // 1. Alpha waveform
  Double_t *alpha_array = new Double_t[n_points];
  for (Int_t i = 0; i < n_points; ++i) {
    alpha_array[i] = f_alpha[i];
  }
  TGraph *gr_alpha = new TGraph(n_points, time_array, alpha_array);
  gr_alpha->SetLineColor(kRed + 1);
  gr_alpha->SetLineWidth(2);
  gr_alpha->SetMarkerStyle(20);
  gr_alpha->SetMarkerSize(0);
  gr_alpha->SetMarkerColor(kRed + 1);

  // Set up axes on first graph
  gr_alpha->SetTitle("");
  gr_alpha->GetXaxis()->SetTitle("Time [ns]");
  gr_alpha->GetYaxis()->SetTitle("Normalized Amplitude [a.u.]");
  gr_alpha->GetXaxis()->SetTitleSize(0.06);
  gr_alpha->GetYaxis()->SetTitleSize(0.06);
  gr_alpha->GetXaxis()->SetLabelSize(0.06);
  gr_alpha->GetYaxis()->SetLabelSize(0.06);
  gr_alpha->Draw("APL");
  graphs.push_back(gr_alpha);

  // 2. Gamma waveform
  Double_t *gamma_array = new Double_t[n_points];
  for (Int_t i = 0; i < n_points; ++i) {
    gamma_array[i] = f_gamma[i];
  }
  TGraph *gr_gamma = new TGraph(n_points, time_array, gamma_array);
  gr_gamma->SetLineColor(kGreen + 2);
  gr_gamma->SetLineWidth(2);
  gr_gamma->SetMarkerStyle(20);
  gr_gamma->SetMarkerSize(0);
  gr_gamma->SetMarkerColor(kGreen + 2);
  gr_gamma->Draw("PL SAME");
  graphs.push_back(gr_gamma);

  // Legend
  leg->AddEntry(gr_alpha, "Am-241 (#alpha)", "lp");
  leg->AddEntry(gr_gamma, "Na-22 (#gamma)", "lp");
  leg->Draw();

  c->Modified();
  c->Update();

  std::string pdf_name = output_name + ".pdf";
  c->SaveAs(pdf_name.c_str());

  // Log scale version
  c->SetLogy();
  gr_alpha->SetMinimum(5e-4);
  gr_alpha->SetMaximum(2);
  c->Modified();
  c->Update();
  c->SaveAs((output_name + "_log.pdf").c_str());

  // Cleanup
  delete[] time_array;
  delete[] alpha_array;
  delete[] gamma_array;
  delete c;
}

void PlottingUtils::PlotSIWeightingFactor(
    const std::vector<Float_t> &f_alpha, const std::vector<Float_t> &f_gamma,
    const std::vector<Float_t> &weighting_factor,
    const std::string &output_name) {

  if (f_alpha.empty() || f_gamma.empty() || weighting_factor.empty()) {
    std::cout << "PlotSIWeightingFactor: Invalid input data." << std::endl;
    return;
  }

  // Create canvas
  TCanvas *c = new TCanvas("c_si_weighting", "SI Weighting Factor", 1200, 800);
  ConfigureCanvas(c, kFALSE);

  // Create time array
  Int_t n_points = weighting_factor.size();
  Double_t *time_array = new Double_t[n_points];
  for (Int_t i = 0; i < n_points; ++i) {
    time_array[i] = 2.0 * i; // 2 ns per sample
  }

  // Create graphs
  std::vector<TGraph *> graphs;
  TLegend *leg = CreateLegend(0.82, 0.6, 0.9, 0.85);

  // 1. f_alpha(t) graph (normalized)
  Double_t *alpha_array = new Double_t[n_points];
  Float_t alpha_peak = *std::max_element(f_alpha.begin(), f_alpha.end());
  for (Int_t i = 0; i < n_points; ++i) {
    alpha_array[i] = f_alpha[i] / alpha_peak;
  }
  TGraph *gr_alpha = new TGraph(n_points, time_array, alpha_array);
  gr_alpha->SetLineColor(kRed + 1);
  gr_alpha->SetLineStyle(kSolid);
  graphs.push_back(gr_alpha);

  // 2. f_gamma(t) graph (normalized)
  Double_t *gamma_array = new Double_t[n_points];
  Float_t gamma_peak = *std::max_element(f_gamma.begin(), f_gamma.end());
  for (Int_t i = 0; i < n_points; ++i) {
    gamma_array[i] = f_gamma[i] / gamma_peak;
  }
  TGraph *gr_gamma = new TGraph(n_points, time_array, gamma_array);
  gr_gamma->SetLineColor(kGreen + 2);
  gr_gamma->SetLineStyle(kSolid);
  graphs.push_back(gr_gamma);

  // 3. p(t) weighting factor graph
  Double_t *pt_array = new Double_t[n_points];
  for (Int_t i = 0; i < n_points; ++i) {
    pt_array[i] = weighting_factor[i];
  }
  TGraph *gr_pt = new TGraph(n_points, time_array, pt_array);
  gr_pt->SetLineColor(kBlue);
  gr_pt->SetLineStyle(kSolid);
  graphs.push_back(gr_pt);

  // Draw first graph to set up axes
  gr_alpha->SetTitle("");
  gr_alpha->GetXaxis()->SetTitle("Time [ns]");
  gr_alpha->GetYaxis()->SetTitle("Normalized Amplitude [a.u.]");
  gr_alpha->GetXaxis()->SetTitleSize(0.05);
  gr_alpha->GetYaxis()->SetTitleSize(0.05);
  gr_alpha->GetXaxis()->SetLabelSize(0.04);
  gr_alpha->GetYaxis()->SetLabelSize(0.04);
  gr_alpha->SetMinimum(-1.2);
  gr_alpha->SetMaximum(1.2);
  gr_alpha->Draw("AL");

  // Draw other graphs
  gr_gamma->Draw("L SAME");
  gr_pt->Draw("L SAME");

  // Add horizontal line at y=0
  TLine *zero_line = new TLine(time_array[0], 0, time_array[n_points - 1], 0);
  zero_line->SetLineStyle(kDotted);
  zero_line->SetLineColor(kBlack);
  zero_line->SetLineWidth(1);
  zero_line->Draw("same");

  // Create legend
  leg->AddEntry(gr_alpha, "f_{#alpha}(t)", "l");
  leg->AddEntry(gr_gamma, "f_{#gamma}(t)", "l");
  leg->AddEntry(gr_pt, "p(t)", "l");
  leg->Draw();

  c->Modified();
  c->Update();

  // Save plots
  std::string pdf_name = output_name + "_with_weighting.pdf";
  c->SaveAs(pdf_name.c_str());

  std::cout << "SI weighting factor plot saved: " << pdf_name << std::endl;

  // Cleanup
  delete[] time_array;
  delete[] alpha_array;
  delete[] gamma_array;
  delete[] pt_array;
  delete c;
  for (auto *gr : graphs) {
    delete gr;
  }
  delete zero_line;
}

void PlottingUtils::PlotPSDSpectra(
    const std::map<Int_t, TH1F *> &psd_spectra,
    const std::map<Int_t, std::string> &source_names,
    const std::string &output_prefix, const std::string &psd_type) {

  for (const auto &pair : psd_spectra) {
    Int_t src_id = pair.first;
    TH1F *hist = pair.second;
    if (!hist)
      continue;

    auto name_it = source_names.find(src_id);
    std::string source_name = (name_it != source_names.end())
                                  ? name_it->second
                                  : "Source_" + std::to_string(src_id);

    TCanvas *c =
        new TCanvas(("c_psd_" + CleanSourceName(source_name)).c_str(),
                    (psd_type + " - " + source_name).c_str(), 1200, 800);
    ConfigureCanvas(c, kTRUE);

    Int_t color = GetSourceColor(src_id);

    // Configure histogram using existing function
    ConfigureHistogram(hist, color);
    hist->GetXaxis()->SetTitle((psd_type).c_str());
    hist->GetYaxis()->SetTitle("Counts");

    hist->Draw("HIST");

    AddSubplotLabel("(" + std::string(1, 'a' + src_id) + ") " + source_name,
                    0.9, 0.83);

    c->Modified();
    c->Update();

    std::string clean_name = CleanSourceName(source_name);
    std::string pdf_name = output_prefix + "_" + clean_name + ".pdf";
    c->SaveAs(pdf_name.c_str());

    delete c;
  }
}

void PlottingUtils::PlotPSD2D(const std::string &waveforms_file,
                              const std::string &psd_branch,
                              const std::string &output_name,
                              const std::string &psd_label) {

  TFile *file = TFile::Open(waveforms_file.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << waveforms_file
              << " for 2D PSD plot!" << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("features"));
  if (!tree) {
    std::cout << "Error: Could not find features tree!" << std::endl;
    file->Close();
    return;
  }

  // Check if required branches exist
  if (!tree->GetBranch(psd_branch.c_str())) {
    std::cout << "Error: Branch " << psd_branch << " not found!" << std::endl;
    file->Close();
    return;
  }
  if (!tree->GetBranch("light_output_keVee")) {
    std::cout << "Error: Branch light_output_keVee not found!" << std::endl;
    file->Close();
    return;
  }

  // Set up variables
  Float_t psd_value;
  Float_t light_output_keVee;
  Int_t source_id;

  tree->SetBranchAddress(psd_branch.c_str(), &psd_value);
  tree->SetBranchAddress("light_output_keVee", &light_output_keVee);
  tree->SetBranchAddress("source_id", &source_id);

  // First pass: determine ranges and collect source IDs
  Long64_t n_entries = tree->GetEntries();
  std::set<Int_t> source_ids;

  for (Long64_t i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);
    source_ids.insert(source_id);
  }

  // Add padding to ranges
  Float_t min_psd = psd_label == "PSP_{CC}" ? 0 : -0.5;
  Int_t max_psd = psd_label == "PSP_{CC}" ? 1 : 0;
  Int_t min_light = 0;
  Int_t max_light = 2000;

  // Create 2D histograms for each source
  std::map<Int_t, TH2F *> hist2d_map;
  std::map<Int_t, std::string> source_names = {
      {0, "Am-241"}, {1, "Cs-137"}, {2, "Na-22"}, {3, "Am-241 & Cs-137"}};

  for (Int_t sid : source_ids) {
    auto name_it = source_names.find(sid);
    std::string source_name = (name_it != source_names.end())
                                  ? name_it->second
                                  : "Source_" + std::to_string(sid);

    TH2F *hist2d =
        new TH2F(("h2d_psd_" + std::to_string(sid)).c_str(),
                 (psd_label + " vs Light Output - " + source_name).c_str(), 100,
                 min_light, max_light,   // x-axis: light output
                 100, min_psd, max_psd); // y-axis: PSD
    hist2d_map[sid] = hist2d;
  }

  // Second pass: fill histograms
  for (Long64_t i = 0; i < n_entries; ++i) {
    tree->GetEntry(i);
    if (hist2d_map.count(source_id)) {
      hist2d_map[source_id]->Fill(light_output_keVee, psd_value);
    }
  }

  for (auto &pair : hist2d_map) {
    Int_t sid = pair.first;
    TH2F *hist2d = pair.second;

    if (hist2d->GetEntries() == 0)
      continue;

    auto name_it = source_names.find(sid);
    std::string source_name = (name_it != source_names.end())
                                  ? name_it->second
                                  : "Source_" + std::to_string(sid);

    TCanvas *c =
        new TCanvas(("c_2d_psd_" + std::to_string(sid)).c_str(),
                    (psd_label + " 2D - " + source_name).c_str(), 1200, 800);

    c->SetLogz(1);           // Log scale for z-axis (counts)
    c->SetRightMargin(0.15); // Make room for color palette
    ConfigureCanvas(c, kFALSE);

    // Configure histogram
    hist2d->SetTitle("");
    hist2d->GetXaxis()->SetTitle("Light Output [keVee]");
    hist2d->GetYaxis()->SetTitle((psd_label).c_str());
    hist2d->GetXaxis()->SetTitleSize(0.05);
    hist2d->GetYaxis()->SetTitleSize(0.05);
    hist2d->GetXaxis()->SetLabelSize(0.045);
    hist2d->GetYaxis()->SetLabelSize(0.045);
    hist2d->GetXaxis()->SetTitleOffset(1.2);
    hist2d->GetYaxis()->SetTitleOffset(1.3);
    hist2d->SetMinimum(1); // Minimum for log scale

    hist2d->Draw("COLZ");

    AddSubplotLabel("(" + std::string(1, 'a' + sid) + ") " + source_name, 0.8,
                    0.83);

    c->Modified();
    c->Update();

    std::string clean_name = CleanSourceName(source_name);
    std::string pdf_name = output_name + "_2d_" + clean_name + ".pdf";
    c->SaveAs(pdf_name.c_str());

    delete c;
  }

  file->Close();
  std::cout << "2D PSD plots saved with prefix: " << output_name << std::endl;
}

Double_t PlottingUtils::CalculateFigureOfMerit(TH1F *hist_alpha,
                                               TH1F *hist_gamma,
                                               Double_t &separation,
                                               Double_t &fwhm_alpha,
                                               Double_t &fwhm_gamma) {
  if (!hist_alpha || !hist_gamma) {
    std::cout << "Error: Invalid histograms for FOM calculation" << std::endl;
    return -1;
  }

  if (hist_alpha->GetEntries() == 0 || hist_gamma->GetEntries() == 0) {
    std::cout << "Error: Empty histograms for FOM calculation" << std::endl;
    return -1;
  }

  // Find meaningful fit ranges based on histogram content
  Int_t first_bin_alpha =
      hist_alpha->FindFirstBinAbove(hist_alpha->GetMaximum() * 0.01);
  Int_t last_bin_alpha =
      hist_alpha->FindLastBinAbove(hist_alpha->GetMaximum() * 0.01);
  Int_t first_bin_gamma =
      hist_gamma->FindFirstBinAbove(hist_gamma->GetMaximum() * 0.01);
  Int_t last_bin_gamma =
      hist_gamma->FindLastBinAbove(hist_gamma->GetMaximum() * 0.01);

  Double_t fit_min_alpha = hist_alpha->GetBinCenter(first_bin_alpha);
  Double_t fit_max_alpha = hist_alpha->GetBinCenter(last_bin_alpha);
  Double_t fit_min_gamma = hist_gamma->GetBinCenter(first_bin_gamma);
  Double_t fit_max_gamma = hist_gamma->GetBinCenter(last_bin_gamma);

  // Create fit functions with reasonable ranges
  TF1 *fit_alpha = new TF1("fit_alpha", "gaus", fit_min_alpha, fit_max_alpha);
  TF1 *fit_gamma = new TF1("fit_gamma", "gaus", fit_min_gamma, fit_max_gamma);

  // Set initial parameters based on histogram statistics
  Double_t alpha_max = hist_alpha->GetMaximum();
  Double_t alpha_mean = hist_alpha->GetMean();
  Double_t alpha_rms = hist_alpha->GetRMS();

  Double_t gamma_max = hist_gamma->GetMaximum();
  Double_t gamma_mean = hist_gamma->GetMean();
  Double_t gamma_rms = hist_gamma->GetRMS();

  // Set reasonable parameter limits
  fit_alpha->SetParameters(alpha_max, alpha_mean, alpha_rms);
  fit_alpha->SetParLimits(0, alpha_max * 0.1, alpha_max * 2.0);
  fit_alpha->SetParLimits(1, alpha_mean - 3 * alpha_rms,
                          alpha_mean + 3 * alpha_rms);
  fit_alpha->SetParLimits(2, alpha_rms * 0.1, alpha_rms * 10.0);

  fit_gamma->SetParameters(gamma_max, gamma_mean, gamma_rms);
  fit_gamma->SetParLimits(0, gamma_max * 0.1, gamma_max * 2.0);
  fit_gamma->SetParLimits(1, gamma_mean - 3 * gamma_rms,
                          gamma_mean + 3 * gamma_rms);
  fit_gamma->SetParLimits(2, gamma_rms * 0.1, gamma_rms * 10.0);

  std::cout << "Fitting alpha distribution..." << std::endl;
  std::cout << "  Range: " << fit_min_alpha << " to " << fit_max_alpha
            << std::endl;
  std::cout << "  Initial params: amp=" << alpha_max << ", mean=" << alpha_mean
            << ", sigma=" << alpha_rms << std::endl;

  // Perform fits with error handling
  Int_t fit_status_alpha = hist_alpha->Fit(fit_alpha, "RQS");
  if (fit_status_alpha != 0) {
    std::cout << "Warning: Alpha fit failed with status " << fit_status_alpha
              << std::endl;
  }

  std::cout << "Fitting gamma distribution..." << std::endl;
  std::cout << "  Range: " << fit_min_gamma << " to " << fit_max_gamma
            << std::endl;
  std::cout << "  Initial params: amp=" << gamma_max << ", mean=" << gamma_mean
            << ", sigma=" << gamma_rms << std::endl;

  Int_t fit_status_gamma = hist_gamma->Fit(fit_gamma, "RQS");
  if (fit_status_gamma != 0) {
    std::cout << "Warning: Gamma fit failed with status " << fit_status_gamma
              << std::endl;
  }

  // Extract parameters with error checking
  Double_t mean_alpha = fit_alpha->GetParameter(1);
  Double_t sigma_alpha = std::abs(fit_alpha->GetParameter(2));
  Double_t mean_gamma = fit_gamma->GetParameter(1);
  Double_t sigma_gamma = std::abs(fit_gamma->GetParameter(2));

  // Sanity check on fit results
  if (sigma_alpha <= 0 || sigma_gamma <= 0) {
    std::cout << "Error: Invalid sigma values from fits" << std::endl;
    return -1;
  }

  // Calculate FWHM = 2.355 * sigma
  fwhm_alpha = 2.355 * sigma_alpha;
  fwhm_gamma = 2.355 * sigma_gamma;

  // Calculate separation
  separation = std::abs(mean_alpha - mean_gamma);

  // Calculate Figure of Merit = separation / (FWHM_alpha + FWHM_gamma)
  Double_t fom = separation / (fwhm_alpha + fwhm_gamma);

  return fom;
}

void PlottingUtils::PlotPSDWithFOM(TH1F *hist_alpha, TH1F *hist_gamma,
                                   const std::string &psd_type,
                                   const std::string &output_name) {
  if (!hist_alpha || !hist_gamma) {
    std::cout << "Error: Invalid histograms for PSD+FOM plot" << std::endl;
    return;
  }

  if (hist_alpha->GetEntries() == 0 || hist_gamma->GetEntries() == 0) {
    std::cout << "Error: Empty histograms for PSD+FOM plot" << std::endl;
    return;
  }

  TCanvas *c = new TCanvas("c_psd_fom", ("PSD with FOM - " + psd_type).c_str(),
                           1200, 800);
  ConfigureCanvas(c, kTRUE);

  // Configure histograms
  ConfigureHistogram(hist_alpha, kRed + 1, "");
  ConfigureHistogram(hist_gamma, kGreen + 2, "");

  hist_alpha->SetFillColorAlpha(kRed + 1, 0.3);
  hist_gamma->SetFillColorAlpha(kGreen + 2, 0.3);

  // Set axes labels
  hist_alpha->GetXaxis()->SetTitle(psd_type.c_str());
  hist_alpha->GetYaxis()->SetTitle("Counts");

  Double_t max_alpha = hist_alpha->GetMaximum();
  Double_t max_gamma = hist_gamma->GetMaximum();
  Double_t y_max = std::max(max_alpha, max_gamma) * 1.2;

  hist_alpha->SetMaximum(y_max);
  hist_alpha->Draw("HIST");
  hist_gamma->Draw("HIST SAME");

  Double_t separation, fwhm_alpha, fwhm_gamma;
  Double_t fom = CalculateFigureOfMerit(hist_alpha, hist_gamma, separation,
                                        fwhm_alpha, fwhm_gamma);

  if (fom < 0) {
    std::cout << "Warning: FOM calculation failed, creating plot without fits"
              << std::endl;
  } else {
    // Draw fits if successful
    TF1 *fit_alpha = hist_alpha->GetFunction("fit_alpha");
    TF1 *fit_gamma = hist_gamma->GetFunction("fit_gamma");

    if (fit_alpha) {
      fit_alpha->SetLineColor(kRed + 2);
      fit_alpha->SetLineWidth(2);
      fit_alpha->SetLineStyle(1);
      fit_alpha->Draw("SAME");
    }

    if (fit_gamma) {
      fit_gamma->SetLineColor(kGreen + 3);
      fit_gamma->SetLineWidth(2);
      fit_gamma->SetLineStyle(1);
      fit_gamma->Draw("SAME");
    }
  }

  // Create legend on the left side
  TLegend *leg = CreateLegend(0.645, 0.55, 0.92, 0.83);
  leg->AddEntry(hist_alpha, "Am-241 (#alpha)", "f");
  leg->AddEntry(hist_gamma, "Na-22 (#gamma)", "f");

  if (fom >= 0) {
    TF1 *fit_alpha = hist_alpha->GetFunction("fit_alpha");
    TF1 *fit_gamma = hist_gamma->GetFunction("fit_gamma");
    if (fit_alpha)
      leg->AddEntry(fit_alpha, "#alpha fit", "l");
    if (fit_gamma)
      leg->AddEntry(fit_gamma, "#gamma fit", "l");

    // Add FOM to legend
    char fom_text[100];
    sprintf(fom_text, "FOM: %.3f", fom);
    leg->AddEntry((TObject *)0, fom_text, "");
  }
  leg->Draw();

  c->Modified();
  c->Update();

  std::string pdf_name = output_name + ".pdf";
  c->SaveAs(pdf_name.c_str());

  std::cout << "PSD+FOM plot saved: " << pdf_name << std::endl;

  delete c;
}

void PlottingUtils::PlotAverageWaveformsByLightOutput(
    const std::string &waveforms_file, const std::string &output_prefix,
    Float_t bin_width, Float_t min_light, Float_t max_light,
    Int_t min_events_per_bin) {

  TFile *file = TFile::Open(waveforms_file.c_str(), "READ");
  if (!file || file->IsZombie()) {
    std::cout << "Error: Could not open " << waveforms_file << std::endl;
    return;
  }

  TTree *tree = static_cast<TTree *>(file->Get("features"));
  if (!tree) {
    std::cout << "Error: Could not find features tree" << std::endl;
    file->Close();
    return;
  }

  // Set up branch addresses
  TArrayS *samples = nullptr;
  Int_t source_id;
  Float_t light_output_keVee;

  tree->SetBranchAddress("Samples", &samples);
  tree->SetBranchAddress("source_id", &source_id);
  tree->SetBranchAddress("light_output_keVee", &light_output_keVee);

  // Create light output bins
  std::vector<Float_t> bin_edges;
  for (Float_t light = min_light; light <= max_light; light += bin_width) {
    bin_edges.push_back(light);
  }
  Int_t n_bins = bin_edges.size() - 1;

  // Source information
  std::map<Int_t, std::string> source_names = {{0, "Am-241"}, {2, "Na-22"}};
  std::vector<Int_t> sources_to_plot = {0, 2}; // Am-241 and Na-22

  // Get waveform length
  tree->GetEntry(0);
  Int_t waveform_length = samples->GetSize();

  for (Int_t src_id : sources_to_plot) {
    std::cout << "Processing " << source_names[src_id] << "..." << std::endl;

    // Storage for waveforms by light output bin
    std::map<Int_t, std::vector<std::vector<Float_t>>> waveforms_by_bin;

    // Collect waveforms
    Long64_t n_entries = tree->GetEntries();
    for (Long64_t i = 0; i < n_entries; ++i) {
      tree->GetEntry(i);

      if (source_id != src_id)
        continue;

      // Find which bin this event belongs to
      Int_t bin_idx = -1;
      for (Int_t b = 0; b < n_bins; ++b) {
        if (light_output_keVee >= bin_edges[b] &&
            light_output_keVee < bin_edges[b + 1]) {
          bin_idx = b;
          break;
        }
      }

      if (bin_idx < 0)
        continue; // Outside range

      // Convert TArrayS to vector
      std::vector<Float_t> waveform;
      for (Int_t j = 0; j < waveform_length; ++j) {
        waveform.push_back(Float_t(samples->At(j)));
      }

      waveforms_by_bin[bin_idx].push_back(waveform);
    }

    // Calculate averages and create plot
    TCanvas *c = new TCanvas(
        ("c_avg_wf_" + source_names[src_id]).c_str(),
        ("Average Waveforms - " + source_names[src_id]).c_str(), 1200, 800);
    ConfigureCanvas(c, kFALSE);

    TLegend *leg = CreateLegend(0.6, 0.5, 0.9, 0.85);
    std::vector<TGraph *> graphs;
    Bool_t first_graph = kTRUE;

    // Time array (2 ns per sample)
    Double_t *time_array = new Double_t[waveform_length];
    for (Int_t i = 0; i < waveform_length; ++i) {
      time_array[i] = 2.0 * i;
    }

    // Get default colors
    auto colors = GetDefaultColors();
    Int_t color_idx = 0;

    for (const auto &bin_pair : waveforms_by_bin) {
      Int_t bin_idx = bin_pair.first;
      const std::vector<std::vector<Float_t>> &waveforms = bin_pair.second;

      if (waveforms.size() < min_events_per_bin) {
        std::cout << "  Skipping bin " << bin_idx << " (only "
                  << waveforms.size() << " events)" << std::endl;
        continue;
      }

      // Calculate average waveform
      std::vector<Float_t> avg_waveform(waveform_length, 0.0);
      for (const auto &wf : waveforms) {
        for (Int_t i = 0; i < waveform_length; ++i) {
          avg_waveform[i] += wf[i];
        }
      }

      // Normalize by number of waveforms
      for (Int_t i = 0; i < waveform_length; ++i) {
        avg_waveform[i] /= Float_t(waveforms.size());
      }

      // Normalize to peak = 1
      Float_t peak =
          *std::max_element(avg_waveform.begin(), avg_waveform.end());
      if (peak > 0) {
        for (Float_t &val : avg_waveform) {
          val /= peak;
        }
      }

      // Create TGraph
      Double_t *amp_array = new Double_t[waveform_length];
      for (Int_t i = 0; i < waveform_length; ++i) {
        amp_array[i] = avg_waveform[i];
      }

      TGraph *gr = new TGraph(waveform_length, time_array, amp_array);
      Int_t color = colors[color_idx % colors.size()];
      gr->SetLineColor(color);
      gr->SetLineWidth(2);
      gr->SetMarkerColor(color);

      if (first_graph) {
        gr->SetTitle("");
        gr->GetXaxis()->SetTitle("Time [ns]");
        gr->GetYaxis()->SetTitle("Normalized Amplitude [a.u.]");
        gr->GetXaxis()->SetTitleSize(0.06);
        gr->GetYaxis()->SetTitleSize(0.06);
        gr->GetXaxis()->SetLabelSize(0.06);
        gr->GetYaxis()->SetLabelSize(0.06);
        gr->GetXaxis()->SetTitleOffset(1.2);
        gr->GetYaxis()->SetTitleOffset(1.2);
        gr->SetMinimum(0.0);
        gr->SetMaximum(1.1);
        gr->Draw("AL");
        first_graph = kFALSE;
      } else {
        gr->Draw("L SAME");
      }

      graphs.push_back(gr);

      // Create legend entry
      Float_t bin_start = bin_edges[bin_idx];
      Float_t bin_end = bin_edges[bin_idx + 1];
      char legend_text[100];
      sprintf(legend_text, "%.0f-%.0f keVee", bin_start, bin_end,
              (Int_t)waveforms.size());
      leg->AddEntry(gr, legend_text, "l");

      delete[] amp_array;
      color_idx++;

      std::cout << "  Added bin " << bin_idx << ": " << bin_start << "-"
                << bin_end << " keVee (" << waveforms.size() << " events)"
                << std::endl;
    }

    // Add subplot label
    std::string subplot_label = src_id == 0 ? "(a) " + source_names[src_id]
                                            : "(b) " + source_names[src_id];
    AddSubplotLabel(subplot_label, 0.5, 0.8);

    // Draw legend
    leg->Draw();

    c->Modified();
    c->Update();

    // Save plot
    std::string clean_name = source_names[src_id];
    std::replace(clean_name.begin(), clean_name.end(), '-', '_');
    std::string pdf_name = output_prefix + "_" + clean_name + ".pdf";
    c->SaveAs(pdf_name.c_str());

    std::cout << "Saved: " << pdf_name << std::endl;

    // Cleanup
    delete[] time_array;
    for (auto *gr : graphs) {
      delete gr;
    }
    delete c;
  }

  file->Close();
  std::cout << "Average waveforms by light output plots completed."
            << std::endl;
}
