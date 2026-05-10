#include "RooFitUtils.hpp"

#include "InteractiveRooFitEditor.hpp"

#include <RooMsgService.h>
#include <TGraph.h>
#include <TMath.h>
#include <TSystem.h>
#include <fstream>
#include <iomanip>
#include <iostream>

RooAbsPdf *RooFitFunctions::MakeGaussian(const TString &name, RooRealVar &x,
                                         RooRealVar &mu, RooRealVar &sigma) {
  return new RooGaussian(name.Data(), name.Data(), x, mu, sigma);
}

RooAbsPdf *RooFitFunctions::MakeStepShelf(const TString &name, RooRealVar &x,
                                          RooRealVar &mu, RooRealVar &sigma) {
  return new RooGenericPdf(name.Data(), name.Data(),
                           "1.0/TMath::Power(1.0+exp((@0-@1)/@2),2)",
                           RooArgList(x, mu, sigma));
}

RooAbsPdf *RooFitFunctions::MakeLowExpTail(const TString &name, RooRealVar &x,
                                           RooRealVar &mu, RooRealVar &sigma,
                                           RooRealVar &tau) {
  return new RooGenericPdf(
      name.Data(), name.Data(),
      "exp((@0-@1)/@3)*TMath::Erfc((@0-@1)/(TMath::Sqrt(2.0)*@2))",
      RooArgList(x, mu, sigma, tau));
}

RooAbsPdf *RooFitFunctions::MakeLowLinTail(const TString &name, RooRealVar &x,
                                           RooRealVar &mu, RooRealVar &sigma,
                                           RooRealVar &slope) {
  return new RooGenericPdf(
      name.Data(), name.Data(),
      "TMath::Max(0.0,1.0+@3*(@0-@1))*TMath::Erfc((@0-@1)/(TMath::Sqrt(2.0)*@2))",
      RooArgList(x, mu, sigma, slope));
}

RooAbsPdf *RooFitFunctions::MakeHighExpTail(const TString &name, RooRealVar &x,
                                            RooRealVar &mu, RooRealVar &sigma,
                                            RooRealVar &tau) {
  return new RooGenericPdf(
      name.Data(), name.Data(),
      "exp((@1-@0)/@3)*TMath::Erfc((@1-@0)/(TMath::Sqrt(2.0)*@2))",
      RooArgList(x, mu, sigma, tau));
}

RooAbsPdf *RooFitFunctions::MakeLinearBackground(const TString &name,
                                                 RooRealVar &x,
                                                 RooRealVar &slope) {
  return new RooPolynomial(name.Data(), name.Data(), x, RooArgList(slope));
}

RooAbsPdf *RooFitFunctions::MakeFlatBackground(const TString &name,
                                               RooRealVar &x) {
  return new RooPolynomial(name.Data(), name.Data(), x, RooArgList());
}

void RooFitUtils::RegisterOwned(RooAbsArg *arg) { owned_args_.push_back(arg); }

RooFitUtils::RooFitUtils(TH1 *working_hist, Float_t fit_range_low,
                         Float_t fit_range_high, Bool_t use_flat_background,
                         Bool_t use_step, Bool_t use_low_exp_tail,
                         Bool_t use_low_lin_tail, Bool_t use_high_exp_tail) {
  working_hist_ = static_cast<TH1 *>(working_hist->Clone());
  fit_range_low_ = fit_range_low;
  fit_range_high_ = fit_range_high;
  use_flat_background_ = use_flat_background;
  use_step_ = use_step;
  use_low_exp_tail_ = use_low_exp_tail;
  use_low_lin_tail_ = use_low_lin_tail;
  use_high_exp_tail_ = use_high_exp_tail;
  use_manual_init_ = kFALSE;
  interactive_ = kFALSE;

  x_ = nullptr;
  data_hist_ = nullptr;
  total_pdf_ = nullptr;
  num_peaks_ = 0;

  RooMsgService::instance().setGlobalKillBelow(RooFit::WARNING);
  RooRealVar::enableSilentClipping();

  std::cout << "Fit configuration:" << std::endl;
  std::cout << std::endl;
  if (use_flat_background_) {
    std::cout << "Background: FLAT" << std::endl;
  } else {
    std::cout << "Background: LINEAR" << std::endl;
  }
  std::cout << "Step function: " << (use_step_ ? "ENABLED" : "DISABLED")
            << std::endl;
  std::cout << "Low exponential tail: "
            << (use_low_exp_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "Low linear tail: "
            << (use_low_lin_tail_ ? "ENABLED" : "DISABLED") << std::endl;
  std::cout << "High exponential tail: "
            << (use_high_exp_tail_ ? "ENABLED" : "DISABLED") << std::endl;
}

RooFitUtils::~RooFitUtils() {
  for (size_t i = 0; i < owned_args_.size(); i++) {
    delete owned_args_[i];
  }
  owned_args_.clear();
  delete data_hist_;
  delete working_hist_;
}

void RooFitUtils::SetManualParameters(const std::vector<Double_t> &params) {
  Int_t expected = num_peaks_ * 10 + 2;
  if (num_peaks_ == 0 || params.size() != (size_t)expected) {
    std::cerr << "ERROR: Manual parameters size (" << params.size()
              << ") doesn't match expected (" << expected << ")" << std::endl;
    return;
  }

  manual_params_ = params;
  use_manual_init_ = kTRUE;

  std::vector<RooRealVar *> all = CollectAllParams();
  for (size_t i = 0; i < all.size(); i++) {
    all[i]->setVal(params[i]);
    all[i]->setConstant(kTRUE);
  }

  std::cout << "Manual parameters set:" << std::endl;
  for (size_t i = 0; i < all.size(); i++) {
    std::cout << "  Par[" << i << "] " << all[i]->GetName() << " = "
              << params[i] << std::endl;
  }
}

void RooFitUtils::SetManualParameter(Int_t index, Double_t value) {
  std::vector<RooRealVar *> all = CollectAllParams();
  if (index < 0 || index >= (Int_t)all.size()) {
    std::cerr << "ERROR: Parameter index " << index << " out of range [0, "
              << (Int_t)all.size() - 1 << "]" << std::endl;
    return;
  }

  if (!use_manual_init_) {
    manual_params_.resize(all.size(), 0.0);
    use_manual_init_ = kTRUE;
  }

  manual_params_[index] = value;
  all[index]->setVal(value);
  all[index]->setConstant(kTRUE);

  std::cout << "Set Par[" << index << "] " << all[index]->GetName() << " = "
            << value << std::endl;
}

Double_t RooFitUtils::EstimateBackground() {
  Int_t left_bin = working_hist_->FindBin(fit_range_low_);
  Int_t right_bin = working_hist_->FindBin(fit_range_high_);

  Int_t n_sideband = (right_bin - left_bin) / 10;
  Double_t left_avg = 0;
  Double_t right_avg = 0;

  for (Int_t i = 0; i < n_sideband; i++) {
    left_avg += working_hist_->GetBinContent(left_bin + i);
    right_avg += working_hist_->GetBinContent(right_bin - i);
  }

  return (left_avg + right_avg) / (2.0 * n_sideband);
}

void RooFitUtils::BuildPeak(Int_t peak_idx, Double_t mu_init,
                            Double_t sigma_init, Double_t peak_height,
                            Double_t range_width) {
  RooFitPeakModel p;
  TString suffix = TString::Format("%d", peak_idx + 1);

  p.mu = new RooRealVar("Mu" + suffix, "Mu" + suffix, mu_init, fit_range_low_,
                        fit_range_high_);
  p.sigma = new RooRealVar("Sigma" + suffix, "Sigma" + suffix, sigma_init,
                           range_width * 0.001, range_width * 0.5);

  Int_t mu_bin = working_hist_->FindBin(mu_init);
  Double_t local_height = working_hist_->GetBinContent(mu_bin);
  Double_t bkg_floor = EstimateBackground();
  Double_t net_height = local_height - bkg_floor;
  if (net_height < 0.1 * local_height)
    net_height = 0.1 * local_height;
  Double_t total_init = net_height * sigma_init *
                         TMath::Sqrt(2.0 * TMath::Pi());
  p.gaus_yield = new RooRealVar("GausAmplitude" + suffix,
                                "GausAmplitude" + suffix, total_init, 0,
                                peak_height * range_width * 10.0);

  p.ratio_step = new RooRealVar("StepAmplitude" + suffix,
                                "StepAmplitude" + suffix, 0.0, 0.0, 0.5);
  p.ratio_low_exp = new RooRealVar("LowExpTailAmplitude" + suffix,
                                   "LowExpTailAmplitude" + suffix, 0.0, 0.0,
                                   0.5);
  p.tau_low_exp = new RooRealVar("LowExpTailDecay" + suffix,
                                 "LowExpTailDecay" + suffix, 1.0, 0.9, 100.0);
  p.ratio_low_lin = new RooRealVar("LowLinTailAmplitude" + suffix,
                                   "LowLinTailAmplitude" + suffix, 0.0, 0.0,
                                   0.5);
  p.slope_low_lin = new RooRealVar("LowLinTailSlope" + suffix,
                                   "LowLinTailSlope" + suffix, 0.0, -0.1, 0.1);
  p.ratio_high_exp = new RooRealVar("HighExpTailAmplitude" + suffix,
                                    "HighExpTailAmplitude" + suffix, 0.0, 0.0,
                                    0.5);
  p.tau_high_exp = new RooRealVar("HighExpTailDecay" + suffix,
                                  "HighExpTailDecay" + suffix, 1.0, 0.9, 100.0);

  RegisterOwned(p.mu);
  RegisterOwned(p.sigma);
  RegisterOwned(p.gaus_yield);
  RegisterOwned(p.ratio_step);
  RegisterOwned(p.ratio_low_exp);
  RegisterOwned(p.tau_low_exp);
  RegisterOwned(p.ratio_low_lin);
  RegisterOwned(p.slope_low_lin);
  RegisterOwned(p.ratio_high_exp);
  RegisterOwned(p.tau_high_exp);

  p.gauss_pdf = RooFitFunctions::MakeGaussian("gauss_pdf" + suffix, *x_, *p.mu,
                                              *p.sigma);
  p.step_pdf = RooFitFunctions::MakeStepShelf("step_pdf" + suffix, *x_, *p.mu,
                                              *p.sigma);
  p.low_exp_pdf = RooFitFunctions::MakeLowExpTail("low_exp_pdf" + suffix, *x_,
                                                  *p.mu, *p.sigma,
                                                  *p.tau_low_exp);
  p.low_lin_pdf = RooFitFunctions::MakeLowLinTail("low_lin_pdf" + suffix, *x_,
                                                  *p.mu, *p.sigma,
                                                  *p.slope_low_lin);
  p.high_exp_pdf = RooFitFunctions::MakeHighExpTail("high_exp_pdf" + suffix,
                                                    *x_, *p.mu, *p.sigma,
                                                    *p.tau_high_exp);

  RegisterOwned(p.gauss_pdf);
  RegisterOwned(p.step_pdf);
  RegisterOwned(p.low_exp_pdf);
  RegisterOwned(p.low_lin_pdf);
  RegisterOwned(p.high_exp_pdf);

  p.step_yield = new RooFormulaVar("step_yield" + suffix, "@0*@1",
                                   RooArgList(*p.gaus_yield, *p.ratio_step));
  p.low_exp_yield = new RooFormulaVar(
      "low_exp_yield" + suffix, "@0*@1",
      RooArgList(*p.gaus_yield, *p.ratio_low_exp));
  p.low_lin_yield = new RooFormulaVar(
      "low_lin_yield" + suffix, "@0*@1",
      RooArgList(*p.gaus_yield, *p.ratio_low_lin));
  p.high_exp_yield = new RooFormulaVar(
      "high_exp_yield" + suffix, "@0*@1",
      RooArgList(*p.gaus_yield, *p.ratio_high_exp));

  RegisterOwned(p.step_yield);
  RegisterOwned(p.low_exp_yield);
  RegisterOwned(p.low_lin_yield);
  RegisterOwned(p.high_exp_yield);

  peaks_.push_back(p);
}

void RooFitUtils::BuildBackground(Double_t bkg_estimate, Double_t peak_height,
                                  Double_t range_width) {
  bkg_.bkg_yield = new RooRealVar("BkgConstant", "BkgConstant",
                                  bkg_estimate * range_width, 0,
                                  peak_height * range_width * 10.0);
  bkg_.bkg_slope = new RooRealVar("BkgSlope", "BkgSlope", 0.0, -1000.0, 1000.0);

  RegisterOwned(bkg_.bkg_yield);
  RegisterOwned(bkg_.bkg_slope);

  if (use_flat_background_) {
    bkg_.bkg_pdf = RooFitFunctions::MakeFlatBackground("bkg_pdf", *x_);
    bkg_.bkg_slope->setVal(0.0);
    bkg_.bkg_slope->setConstant(kTRUE);
  } else {
    bkg_.bkg_pdf = RooFitFunctions::MakeLinearBackground("bkg_pdf", *x_,
                                                          *bkg_.bkg_slope);
  }
  RegisterOwned(bkg_.bkg_pdf);
}

void RooFitUtils::BuildTotalModel() {
  RooArgList pdf_list;
  RooArgList coef_list;
  for (size_t pi = 0; pi < peaks_.size(); pi++) {
    pdf_list.add(*peaks_[pi].gauss_pdf);
    coef_list.add(*peaks_[pi].gaus_yield);
    pdf_list.add(*peaks_[pi].step_pdf);
    coef_list.add(*peaks_[pi].step_yield);
    pdf_list.add(*peaks_[pi].low_exp_pdf);
    coef_list.add(*peaks_[pi].low_exp_yield);
    pdf_list.add(*peaks_[pi].low_lin_pdf);
    coef_list.add(*peaks_[pi].low_lin_yield);
    pdf_list.add(*peaks_[pi].high_exp_pdf);
    coef_list.add(*peaks_[pi].high_exp_yield);
  }
  pdf_list.add(*bkg_.bkg_pdf);
  coef_list.add(*bkg_.bkg_yield);

  total_pdf_ = new RooAddPdf("total_pdf", "total_pdf", pdf_list, coef_list);
  RegisterOwned(total_pdf_);
}

void RooFitUtils::ConfigureComponentFlagsForPeak(Int_t peak_idx) {
  RooFitPeakModel &p = peaks_[peak_idx];

  if (use_step_) {
    p.ratio_step->setVal(0.0);
    p.ratio_step->setConstant(kFALSE);
  } else {
    p.ratio_step->setVal(0.0);
    p.ratio_step->setConstant(kTRUE);
  }

  if (use_low_exp_tail_) {
    p.ratio_low_exp->setVal(0.1);
    p.ratio_low_exp->setConstant(kFALSE);
    p.tau_low_exp->setVal(1.0);
    p.tau_low_exp->setConstant(kFALSE);
  } else {
    p.ratio_low_exp->setVal(0.0);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(1.0);
    p.tau_low_exp->setConstant(kTRUE);
  }

  if (use_low_lin_tail_) {
    p.ratio_low_lin->setVal(0.1);
    p.ratio_low_lin->setConstant(kFALSE);
    p.slope_low_lin->setVal(0.0);
    p.slope_low_lin->setConstant(kFALSE);
  } else {
    p.ratio_low_lin->setVal(0.0);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(0.0);
    p.slope_low_lin->setConstant(kTRUE);
  }

  if (use_high_exp_tail_) {
    p.ratio_high_exp->setVal(0.1);
    p.ratio_high_exp->setConstant(kFALSE);
    p.tau_high_exp->setVal(1.0);
    p.tau_high_exp->setConstant(kFALSE);
  } else {
    p.ratio_high_exp->setVal(0.0);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
}

void RooFitUtils::FixComponent(Int_t peak_idx, const TString &component) {
  RooFitPeakModel &p = peaks_[peak_idx];
  if (component == "step") {
    p.ratio_step->setVal(0.0);
    p.ratio_step->setConstant(kTRUE);
  } else if (component == "low_exp") {
    p.ratio_low_exp->setVal(0.0);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(1.0);
    p.tau_low_exp->setConstant(kTRUE);
  } else if (component == "low_lin") {
    p.ratio_low_lin->setVal(0.0);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(0.0);
    p.slope_low_lin->setConstant(kTRUE);
  } else if (component == "high_exp") {
    p.ratio_high_exp->setVal(0.0);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
}

void RooFitUtils::ReleaseComponent(Int_t peak_idx, const TString &component) {
  RooFitPeakModel &p = peaks_[peak_idx];
  if (component == "step") {
    p.ratio_step->setConstant(kFALSE);
    p.ratio_step->setVal(0.15);
  } else if (component == "low_exp") {
    p.ratio_low_exp->setConstant(kFALSE);
    p.tau_low_exp->setConstant(kFALSE);
    p.ratio_low_exp->setVal(0.15);
    p.tau_low_exp->setVal(1.0);
  } else if (component == "low_lin") {
    p.ratio_low_lin->setConstant(kFALSE);
    p.slope_low_lin->setConstant(kFALSE);
    p.ratio_low_lin->setVal(0.15);
    p.slope_low_lin->setVal(0.0);
  } else if (component == "high_exp") {
    p.ratio_high_exp->setConstant(kFALSE);
    p.tau_high_exp->setConstant(kFALSE);
    p.ratio_high_exp->setVal(0.15);
    p.tau_high_exp->setVal(1.0);
  }
}

std::vector<RooRealVar *> RooFitUtils::CollectAllParams() {
  std::vector<RooRealVar *> out;
  for (size_t pi = 0; pi < peaks_.size(); pi++) {
    out.push_back(peaks_[pi].mu);
    out.push_back(peaks_[pi].sigma);
    out.push_back(peaks_[pi].gaus_yield);
    out.push_back(peaks_[pi].ratio_step);
    out.push_back(peaks_[pi].ratio_low_exp);
    out.push_back(peaks_[pi].tau_low_exp);
    out.push_back(peaks_[pi].ratio_low_lin);
    out.push_back(peaks_[pi].slope_low_lin);
    out.push_back(peaks_[pi].ratio_high_exp);
    out.push_back(peaks_[pi].tau_high_exp);
  }
  out.push_back(bkg_.bkg_yield);
  out.push_back(bkg_.bkg_slope);
  return out;
}

std::vector<RooRealVar *> RooFitUtils::CollectFloatingParams() {
  std::vector<RooRealVar *> all = CollectAllParams();
  std::vector<RooRealVar *> out;
  for (size_t i = 0; i < all.size(); i++) {
    if (!all[i]->isConstant()) {
      out.push_back(all[i]);
    }
  }
  return out;
}

RooFitResult *RooFitUtils::RunFit(Bool_t quiet) {
  Int_t print_level = quiet ? -1 : 0;
  RooFitResult *result =
      total_pdf_->fitTo(*data_hist_, RooFit::Save(kTRUE),
                         RooFit::Extended(kTRUE),
                         RooFit::Range(kFitRangeName),
                         RooFit::SumW2Error(kFALSE),
                         RooFit::PrintLevel(print_level),
                         RooFit::Strategy(2),
                         RooFit::Minimizer("Minuit2", "migrad"),
                         RooFit::EvalBackend::Cpu());
  return result;
}

Double_t RooFitUtils::ComputeReducedChi2(RooFitResult *fit_result,
                                          Int_t &ndof) {
  RooArgSet nset(*x_);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  Double_t bin_width = working_hist_->GetBinWidth(1);
  Double_t saved = x_->getVal();

  Double_t chi2 = 0;
  Int_t nbins_in_range = 0;
  Int_t nbins_hist = working_hist_->GetNbinsX();
  for (Int_t i = 1; i <= nbins_hist; i++) {
    Double_t xv = working_hist_->GetBinCenter(i);
    if (xv < fit_range_low_ || xv > fit_range_high_)
      continue;
    Double_t data = working_hist_->GetBinContent(i);
    Double_t error = working_hist_->GetBinError(i);
    if (error <= 0 || data <= 0)
      continue;
    x_->setVal(xv);
    Double_t fit_val = total_exp * total_pdf_->getVal(&nset) * bin_width;
    Double_t residual = (data - fit_val) / error;
    chi2 += residual * residual;
    nbins_in_range++;
  }
  x_->setVal(saved);

  Int_t npars = fit_result ? fit_result->floatParsFinal().size()
                            : (Int_t)CollectFloatingParams().size();
  ndof = nbins_in_range - npars;
  if (ndof <= 0)
    return -1;
  return chi2 / ndof;
}

void RooFitUtils::SnapshotParams(std::vector<Double_t> &vals,
                                  std::vector<Double_t> &errs,
                                  std::vector<Bool_t> &consts) {
  std::vector<RooRealVar *> all = CollectAllParams();
  vals.resize(all.size());
  errs.resize(all.size());
  consts.resize(all.size());
  for (size_t i = 0; i < all.size(); i++) {
    vals[i] = all[i]->getVal();
    errs[i] = all[i]->getError();
    consts[i] = all[i]->isConstant();
  }
}

void RooFitUtils::RestoreParams(const std::vector<Double_t> &vals,
                                 const std::vector<Double_t> &errs,
                                 const std::vector<Bool_t> &consts) {
  std::vector<RooRealVar *> all = CollectAllParams();
  for (size_t i = 0; i < all.size() && i < vals.size(); i++) {
    all[i]->setVal(vals[i]);
    all[i]->setError(errs[i]);
    all[i]->setConstant(consts[i]);
  }
}

void RooFitUtils::TestLowSideGroup(Int_t peak_idx, Double_t &best_chi2,
                                    std::vector<Double_t> &best_vals,
                                    std::vector<Double_t> &best_errs,
                                    std::vector<Bool_t> &best_const) {
  Bool_t any_low_side = use_step_ || use_low_exp_tail_ || use_low_lin_tail_;
  if (!any_low_side)
    return;

  std::cout << "Testing low-side component group for peak " << peak_idx + 1
            << "..." << std::endl;
  if (use_step_)
    ReleaseComponent(peak_idx, "step");
  if (use_low_exp_tail_)
    ReleaseComponent(peak_idx, "low_exp");
  if (use_low_lin_tail_)
    ReleaseComponent(peak_idx, "low_lin");

  RooFitResult *group_fit = RunFit(kTRUE);
  Bool_t group_ok = group_fit && group_fit->status() == 0;
  Int_t tmp_ndof = 0;
  Double_t chi2_group = group_ok ? ComputeReducedChi2(group_fit, tmp_ndof) : -1;
  delete group_fit;

  if (group_ok && chi2_group < best_chi2) {
    std::cout << "Low-side group peak " << peak_idx + 1
              << " ACCEPTED, pruning..." << std::endl;
    best_chi2 = chi2_group;
    SnapshotParams(best_vals, best_errs, best_const);

    const TString comps[3] = {"step", "low_exp", "low_lin"};
    Bool_t enabled[3] = {use_step_, use_low_exp_tail_, use_low_lin_tail_};
    for (Int_t ci = 0; ci < 3; ci++) {
      if (!enabled[ci])
        continue;
      FixComponent(peak_idx, comps[ci]);
      RooFitResult *pf = RunFit(kTRUE);
      Bool_t ok = pf && pf->status() == 0;
      Int_t nd = 0;
      Double_t c2 = ok ? ComputeReducedChi2(pf, nd) : -1;
      delete pf;
      if (ok && c2 <= best_chi2) {
        std::cout << "  " << comps[ci] << " peak " << peak_idx + 1 << " pruned"
                  << std::endl;
        best_chi2 = c2;
        SnapshotParams(best_vals, best_errs, best_const);
      } else {
        std::cout << "  " << comps[ci] << " peak " << peak_idx + 1 << " retained"
                  << std::endl;
        ReleaseComponent(peak_idx, comps[ci]);
        RestoreParams(best_vals, best_errs, best_const);
      }
    }
  } else {
    std::cout << "Low-side group peak " << peak_idx + 1 << " REJECTED"
              << std::endl;
    if (use_step_)
      FixComponent(peak_idx, "step");
    if (use_low_exp_tail_)
      FixComponent(peak_idx, "low_exp");
    if (use_low_lin_tail_)
      FixComponent(peak_idx, "low_lin");
    RestoreParams(best_vals, best_errs, best_const);
  }
}

void RooFitUtils::TestHighTailIndependent(Int_t peak_idx, Double_t &best_chi2,
                                           std::vector<Double_t> &best_vals,
                                           std::vector<Double_t> &best_errs,
                                           std::vector<Bool_t> &best_const) {
  if (!use_high_exp_tail_)
    return;

  std::cout << "Testing high exponential tail for peak " << peak_idx + 1
            << "..." << std::endl;
  ReleaseComponent(peak_idx, "high_exp");
  RooFitResult *htail_fit = RunFit(kTRUE);
  Bool_t htail_ok = htail_fit && htail_fit->status() == 0;
  Int_t tmp_ndof = 0;
  Double_t chi2_htail = htail_ok ? ComputeReducedChi2(htail_fit, tmp_ndof) : -1;
  delete htail_fit;
  if (htail_ok && chi2_htail < best_chi2) {
    std::cout << "High exp tail peak " << peak_idx + 1 << " ACCEPTED"
              << std::endl;
    best_chi2 = chi2_htail;
    SnapshotParams(best_vals, best_errs, best_const);
  } else {
    std::cout << "High exp tail peak " << peak_idx + 1 << " REJECTED"
              << std::endl;
    FixComponent(peak_idx, "high_exp");
    RestoreParams(best_vals, best_errs, best_const);
  }
}

PeakFitResult RooFitUtils::ExtractPeakResult(Int_t peak_idx) {
  RooFitPeakModel &p = peaks_[peak_idx];
  PeakFitResult result;
  result.mu = p.mu->getVal();
  result.mu_error = p.mu->getError();
  result.sigma = p.sigma->getVal();
  result.sigma_error = p.sigma->getError();
  result.gaus_amplitude = p.gaus_yield->getVal();
  result.gaus_amplitude_error = p.gaus_yield->getError();

  Double_t ga = result.gaus_amplitude;
  result.step_amplitude = p.ratio_step->getVal() * ga;
  result.step_amplitude_error = p.ratio_step->getError() * ga;
  result.low_exp_tail_amplitude = p.ratio_low_exp->getVal() * ga;
  result.low_exp_tail_amplitude_error = p.ratio_low_exp->getError() * ga;
  result.low_exp_tail_decay = p.tau_low_exp->getVal();
  result.low_exp_tail_decay_error = p.tau_low_exp->getError();
  result.low_lin_tail_amplitude = p.ratio_low_lin->getVal() * ga;
  result.low_lin_tail_amplitude_error = p.ratio_low_lin->getError() * ga;
  result.low_lin_tail_slope = p.slope_low_lin->getVal();
  result.low_lin_tail_slope_error = p.slope_low_lin->getError();
  result.high_exp_tail_amplitude = p.ratio_high_exp->getVal() * ga;
  result.high_exp_tail_amplitude_error = p.ratio_high_exp->getError() * ga;
  result.high_exp_tail_decay = p.tau_high_exp->getVal();
  result.high_exp_tail_decay_error = p.tau_high_exp->getError();
  return result;
}

void RooFitUtils::SortPeaksByMu(Int_t num_peaks) {
  for (Int_t i = 0; i < num_peaks - 1; i++) {
    for (Int_t j = 0; j < num_peaks - i - 1; j++) {
      Double_t mu_j = peaks_[j].mu->getVal();
      Double_t mu_next = peaks_[j + 1].mu->getVal();
      if (mu_j > mu_next) {
        std::cout << "Sorting peaks: swapping peak " << j + 1 << " (mu=" << mu_j
                  << ") and peak " << j + 2 << " (mu=" << mu_next << ")"
                  << std::endl;
        RooFitPeakModel tmp = peaks_[j];
        peaks_[j] = peaks_[j + 1];
        peaks_[j + 1] = tmp;
      }
    }
  }
}

void RooFitUtils::SaveInteractiveParams(const TString &input_name,
                                        const TString &peak_name) {
  TString fits_dir = PlottingUtils::GetPlotsBaseDir() + "/fits";
  gSystem->mkdir(fits_dir, kTRUE);
  TString filename = fits_dir + "/" + peak_name + "_" + input_name + ".roofits";
  std::ofstream out(filename.Data());
  if (!out.is_open()) {
    std::cerr << "WARNING: Could not save interactive params to " << filename
              << std::endl;
    return;
  }
  out << std::setprecision(15);
  out << "RANGE " << fit_range_low_ << " " << fit_range_high_ << "\n";
  std::vector<RooRealVar *> all = CollectAllParams();
  for (size_t i = 0; i < all.size(); i++) {
    out << all[i]->GetName() << " " << all[i]->getVal() << " "
        << (all[i]->isConstant() ? 1 : 0) << "\n";
  }
  out.close();
  std::cout << "Saved interactive params to " << filename << std::endl;
}

Bool_t RooFitUtils::LoadInteractiveParams(const TString &input_name,
                                          const TString &peak_name) {
  TString filename = PlottingUtils::GetPlotsBaseDir() + "/fits/" + peak_name +
                     "_" + input_name + ".roofits";
  std::ifstream in(filename.Data());
  if (!in.is_open())
    return kFALSE;

  std::vector<RooRealVar *> all = CollectAllParams();
  std::string token;
  Int_t idx = 0;

  in >> token;
  if (token == "RANGE") {
    Double_t rlo, rhi;
    in >> rlo >> rhi;
    fit_range_low_ = rlo;
    fit_range_high_ = rhi;
    x_->setRange(rlo, rhi);
    x_->setRange(kFitRangeName, rlo, rhi);
  }

  Double_t value;
  Int_t fixed;
  while (in >> token >> value >> fixed && idx < (Int_t)all.size()) {
    all[idx]->setVal(value);
    all[idx]->setConstant(fixed ? kTRUE : kFALSE);
    idx++;
  }
  in.close();

  if (idx != (Int_t)all.size()) {
    std::cerr << "WARNING: Parameter count mismatch in " << filename
              << " (expected " << all.size() << ", got " << idx << ")"
              << std::endl;
    return kFALSE;
  }

  std::cout << "Loaded interactive params from " << filename << std::endl;
  return kTRUE;
}

void RooFitUtils::AppendPeakGraphs(std::vector<TGraph *> &components,
                                   Int_t peak_idx, Style_t line_style,
                                   RooAbsPdf *background_pdf,
                                   Double_t bkg_yield_val, Int_t npts,
                                   Double_t x_step, Double_t bin_width) {
  RooFitPeakModel &p = peaks_[peak_idx];
  Width_t line_width = PlottingUtils::GetLineWidth();
  RooArgSet nset(*x_);

  TGraph *peak_graph = new TGraph(npts);
  Double_t gy = p.gaus_yield->getVal();
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = gy * p.gauss_pdf->getVal(&nset) * bin_width;
    Double_t bkg_v =
        bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
    peak_graph->SetPoint(i, xv, y + bkg_v);
  }
  peak_graph->SetLineColor(kBlack);
  peak_graph->SetLineStyle(line_style);
  peak_graph->SetLineWidth(line_width);
  components.push_back(peak_graph);

  if (TMath::Abs(p.ratio_step->getVal()) > 1e-6) {
    TGraph *step_graph = new TGraph(npts);
    Double_t sy = p.step_yield->getVal();
    for (Int_t i = 0; i < npts; i++) {
      Double_t xv = fit_range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t y = sy * p.step_pdf->getVal(&nset) * bin_width;
      Double_t bkg_v =
          bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
      step_graph->SetPoint(i, xv, y + bkg_v);
    }
    step_graph->SetLineColor(kGray);
    step_graph->SetLineStyle(line_style);
    step_graph->SetLineWidth(line_width);
    components.push_back(step_graph);
  }

  if (TMath::Abs(p.ratio_low_exp->getVal()) > 1e-6 ||
      TMath::Abs(p.ratio_low_lin->getVal()) > 1e-6) {
    TGraph *low_tail_graph = new TGraph(npts);
    Double_t lexp_y = p.low_exp_yield->getVal();
    Double_t llin_y = p.low_lin_yield->getVal();
    for (Int_t i = 0; i < npts; i++) {
      Double_t xv = fit_range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t y_exp = lexp_y * p.low_exp_pdf->getVal(&nset) * bin_width;
      Double_t y_lin = llin_y * p.low_lin_pdf->getVal(&nset) * bin_width;
      Double_t bkg_v =
          bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
      low_tail_graph->SetPoint(i, xv, y_exp + y_lin + bkg_v);
    }
    low_tail_graph->SetLineColor(kRed);
    low_tail_graph->SetLineStyle(line_style);
    low_tail_graph->SetLineWidth(line_width);
    components.push_back(low_tail_graph);
  }

  if (TMath::Abs(p.ratio_high_exp->getVal()) > 1e-6) {
    TGraph *high_tail_graph = new TGraph(npts);
    Double_t hexp_y = p.high_exp_yield->getVal();
    for (Int_t i = 0; i < npts; i++) {
      Double_t xv = fit_range_low_ + i * x_step;
      x_->setVal(xv);
      Double_t y = hexp_y * p.high_exp_pdf->getVal(&nset) * bin_width;
      Double_t bkg_v =
          bkg_yield_val * background_pdf->getVal(&nset) * bin_width;
      high_tail_graph->SetPoint(i, xv, y + bkg_v);
    }
    high_tail_graph->SetLineColor(kOrange);
    high_tail_graph->SetLineStyle(line_style);
    high_tail_graph->SetLineWidth(line_width);
    components.push_back(high_tail_graph);
  }
}

void RooFitUtils::PlotFitSinglePeak(const TString input_name,
                                    const TString peak_name,
                                    const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();
  Double_t bin_width = working_hist_->GetBinWidth(1);
  RooArgSet nset(*x_);

  TGraph *total_graph = new TGraph(npts);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = total_exp * total_pdf_->getVal(&nset) * bin_width;
    total_graph->SetPoint(i, xv, y);
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  Double_t bkg_yield_val = bkg_.bkg_yield->getVal();
  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = bkg_yield_val * bkg_.bkg_pdf->getVal(&nset) * bin_width;
    background_graph->SetPoint(i, xv, y);
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);
  AppendPeakGraphs(components, 0, 1, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                    bin_width);

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);
}

void RooFitUtils::PlotFitDoublePeak(const TString input_name,
                                    const TString peak_name,
                                    const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();
  Double_t bin_width = working_hist_->GetBinWidth(1);
  RooArgSet nset(*x_);

  TGraph *total_graph = new TGraph(npts);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = total_exp * total_pdf_->getVal(&nset) * bin_width;
    total_graph->SetPoint(i, xv, y);
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  Double_t bkg_yield_val = bkg_.bkg_yield->getVal();
  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = bkg_yield_val * bkg_.bkg_pdf->getVal(&nset) * bin_width;
    background_graph->SetPoint(i, xv, y);
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);
  AppendPeakGraphs(components, 0, 1, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                    bin_width);
  AppendPeakGraphs(components, 1, 3, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                    bin_width);

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);
}

void RooFitUtils::PlotFitTriplePeak(const TString input_name,
                                    const TString peak_name,
                                    const TString label) {
  Int_t npts = 1000;
  Double_t x_step = (fit_range_high_ - fit_range_low_) / (npts - 1);
  Width_t line_width = PlottingUtils::GetLineWidth();
  Double_t bin_width = working_hist_->GetBinWidth(1);
  RooArgSet nset(*x_);

  TGraph *total_graph = new TGraph(npts);
  Double_t total_exp = total_pdf_->expectedEvents(&nset);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = total_exp * total_pdf_->getVal(&nset) * bin_width;
    total_graph->SetPoint(i, xv, y);
  }
  total_graph->SetLineColor(kAzure);
  total_graph->SetLineWidth(line_width);

  Double_t bkg_yield_val = bkg_.bkg_yield->getVal();
  TGraph *background_graph = new TGraph(npts);
  for (Int_t i = 0; i < npts; i++) {
    Double_t xv = fit_range_low_ + i * x_step;
    x_->setVal(xv);
    Double_t y = bkg_yield_val * bkg_.bkg_pdf->getVal(&nset) * bin_width;
    background_graph->SetPoint(i, xv, y);
  }
  background_graph->SetLineColor(kGreen);
  background_graph->SetLineWidth(line_width);

  std::vector<TGraph *> components;
  components.push_back(background_graph);
  AppendPeakGraphs(components, 0, 1, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                    bin_width);
  AppendPeakGraphs(components, 1, 3, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                    bin_width);
  AppendPeakGraphs(components, 2, 4, bkg_.bkg_pdf, bkg_yield_val, npts, x_step,
                    bin_width);

  PlottingUtils::PlotFitWithResiduals(
      working_hist_, total_graph, components, fit_range_low_, fit_range_high_,
      peak_name + "_" + input_name, "fits", label, kTRUE);
}

FitResult RooFitUtils::FitSinglePeak(const TString input_name,
                                     const TString peak_name) {
  FitResult results;
  results.peaks.emplace_back();

  num_peaks_ = 1;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t mu_init = (fit_range_low_ + fit_range_high_) / 2;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  BuildPeak(0, mu_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  data_hist_ = new RooDataHist("data_hist", "data_hist", RooArgList(*x_),
                                working_hist_);
  x_->setRange(fit_range_low_, fit_range_high_);

  ConfigureComponentFlagsForPeak(0);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(working_hist_, total_pdf_, x_,
                                         data_hist_, &peaks_, &bkg_,
                                         fit_range_low_, fit_range_high_,
                                         peak_name + " / " + input_name)) {
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(0, "step");
    FixComponent(0, "low_exp");
    FixComponent(0, "low_lin");
    FixComponent(0, "high_exp");

    if (use_manual_init_) {
      std::cout << "Using manually initialized parameters" << std::endl;
      std::vector<RooRealVar *> all = CollectAllParams();
      for (size_t i = 0; i < manual_params_.size() && i < all.size(); i++) {
        all[i]->setVal(manual_params_[i]);
      }
    }

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial fit failed" << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(0, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(0, best_chi2, best_vals, best_errs, best_const);

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitSinglePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Fit did not converge" << std::endl;
  }

  return results;
}

FitResult RooFitUtils::FitDoublePeak(const TString input_name,
                                     const TString peak_name,
                                     Double_t mu1_init, Double_t mu2_init) {
  FitResult results;
  results.peaks.emplace_back();
  results.peaks.emplace_back();

  if (mu1_init > mu2_init) {
    std::cout << "Warning: mu1_init > mu2_init, swapping initial values"
              << std::endl;
    Double_t tmp = mu1_init;
    mu1_init = mu2_init;
    mu2_init = tmp;
  }

  num_peaks_ = 2;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  BuildPeak(0, mu1_init, sigma_init, peak_height, range_width);
  BuildPeak(1, mu2_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  data_hist_ = new RooDataHist("data_hist", "data_hist", RooArgList(*x_),
                                working_hist_);
  x_->setRange(fit_range_low_, fit_range_high_);

  ConfigureComponentFlagsForPeak(0);
  ConfigureComponentFlagsForPeak(1);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      std::cout << "Refit from saved params chi2/ndf = " << final_chi2
                << std::endl;
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(working_hist_, total_pdf_, x_,
                                         data_hist_, &peaks_, &bkg_,
                                         fit_range_low_, fit_range_high_,
                                         peak_name + " / " + input_name)) {
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        std::cout << "Interactive chi2/ndf = " << final_chi2 << std::endl;
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(0, "step");
    FixComponent(0, "low_exp");
    FixComponent(0, "low_lin");
    FixComponent(0, "high_exp");
    FixComponent(1, "step");
    FixComponent(1, "low_exp");
    FixComponent(1, "low_lin");
    FixComponent(1, "high_exp");

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial double peak fit failed" << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(0, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(1, best_chi2, best_vals, best_errs, best_const);

    {
      std::cout
          << "Testing inter-peak group (peak1 high tail + peak2 low-side)..."
          << std::endl;
      if (use_high_exp_tail_)
        ReleaseComponent(0, "high_exp");
      if (use_step_)
        ReleaseComponent(1, "step");
      if (use_low_exp_tail_)
        ReleaseComponent(1, "low_exp");
      if (use_low_lin_tail_)
        ReleaseComponent(1, "low_lin");

      RooFitResult *group_fit = RunFit(kTRUE);
      Bool_t ok = group_fit && group_fit->status() == 0;
      Int_t nd = 0;
      Double_t c2 = ok ? ComputeReducedChi2(group_fit, nd) : -1;
      delete group_fit;

      if (ok && c2 < best_chi2) {
        std::cout << "Inter-peak group ACCEPTED, pruning..." << std::endl;
        best_chi2 = c2;
        SnapshotParams(best_vals, best_errs, best_const);

        struct CompRef { Int_t peak_idx; TString comp; Bool_t enabled; };
        CompRef refs[4] = {
            {0, "high_exp", use_high_exp_tail_},
            {1, "step", use_step_},
            {1, "low_exp", use_low_exp_tail_},
            {1, "low_lin", use_low_lin_tail_},
        };
        for (Int_t ri = 0; ri < 4; ri++) {
          if (!refs[ri].enabled)
            continue;
          FixComponent(refs[ri].peak_idx, refs[ri].comp);
          RooFitResult *pf = RunFit(kTRUE);
          Bool_t pok = pf && pf->status() == 0;
          Int_t pnd = 0;
          Double_t pc2 = pok ? ComputeReducedChi2(pf, pnd) : -1;
          delete pf;
          if (pok && pc2 <= best_chi2) {
            std::cout << "  " << refs[ri].comp << " peak "
                      << refs[ri].peak_idx + 1 << " pruned" << std::endl;
            best_chi2 = pc2;
            SnapshotParams(best_vals, best_errs, best_const);
          } else {
            std::cout << "  " << refs[ri].comp << " peak "
                      << refs[ri].peak_idx + 1 << " retained" << std::endl;
            ReleaseComponent(refs[ri].peak_idx, refs[ri].comp);
            RestoreParams(best_vals, best_errs, best_const);
          }
        }
      } else {
        std::cout << "Inter-peak group REJECTED" << std::endl;
        if (use_high_exp_tail_)
          FixComponent(0, "high_exp");
        if (use_step_)
          FixComponent(1, "step");
        if (use_low_exp_tail_)
          FixComponent(1, "low_exp");
        if (use_low_lin_tail_)
          FixComponent(1, "low_lin");
        RestoreParams(best_vals, best_errs, best_const);
      }
    }

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Double peak fit converged successfully" << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    } else {
      std::cout << "ERROR: Double peak fit failed to converge" << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    SortPeaksByMu(2);
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitDoublePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.peaks[1] = ExtractPeakResult(1);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  } else {
    std::cout << "ERROR: Double peak fit failed" << std::endl;
  }

  return results;
}

FitResult RooFitUtils::FitDoublePeak(const TString input_name,
                                     const TString peak_name,
                                     const PeakFitResult &constrained_peak,
                                     Double_t mu2_init) {
  FitResult results;
  results.peaks.emplace_back();
  results.peaks.emplace_back();

  num_peaks_ = 2;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  BuildPeak(0, constrained_peak.mu, constrained_peak.sigma, peak_height,
            range_width);
  BuildPeak(1, mu2_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  data_hist_ = new RooDataHist("data_hist", "data_hist", RooArgList(*x_),
                                working_hist_);
  x_->setRange(fit_range_low_, fit_range_high_);

  {
    RooFitPeakModel &p = peaks_[0];
    Double_t cga = constrained_peak.gaus_amplitude;
    p.mu->setVal(constrained_peak.mu);
    p.mu->setConstant(kTRUE);
    p.sigma->setVal(constrained_peak.sigma);
    p.sigma->setConstant(kTRUE);
    p.gaus_yield->setVal(cga);
    p.gaus_yield->setRange(0, peak_height * range_width * 10.0);
    p.gaus_yield->setConstant(kFALSE);
    p.ratio_step->setVal(constrained_peak.step_amplitude / cga);
    p.ratio_step->setConstant(kTRUE);
    p.ratio_low_exp->setVal(constrained_peak.low_exp_tail_amplitude / cga);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(constrained_peak.low_exp_tail_decay > 0
                              ? constrained_peak.low_exp_tail_decay
                              : 1.0);
    p.tau_low_exp->setConstant(kTRUE);
    p.ratio_low_lin->setVal(constrained_peak.low_lin_tail_amplitude / cga);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(constrained_peak.low_lin_tail_slope);
    p.slope_low_lin->setConstant(kTRUE);
    p.ratio_high_exp->setVal(constrained_peak.high_exp_tail_amplitude / cga);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(constrained_peak.high_exp_tail_decay > 0
                                ? constrained_peak.high_exp_tail_decay
                                : 1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
  ConfigureComponentFlagsForPeak(1);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(working_hist_, total_pdf_, x_,
                                         data_hist_, &peaks_, &bkg_,
                                         fit_range_low_, fit_range_high_,
                                         peak_name + " / " + input_name)) {
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(1, "step");
    FixComponent(1, "low_exp");
    FixComponent(1, "low_lin");
    FixComponent(1, "high_exp");

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial constrained double peak fit failed"
                << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(1, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(1, best_chi2, best_vals, best_errs, best_const);

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    } else {
      std::cout << "ERROR: Constrained double peak fit failed" << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    SortPeaksByMu(2);
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitDoublePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.peaks[1] = ExtractPeakResult(1);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  }

  return results;
}

FitResult RooFitUtils::FitTriplePeak(const TString input_name,
                                     const TString peak_name,
                                     const FitResult &constrained_peaks,
                                     Double_t mu3_init) {
  FitResult results;
  results.peaks.emplace_back();
  results.peaks.emplace_back();
  results.peaks.emplace_back();

  num_peaks_ = 3;
  Double_t range_width = fit_range_high_ - fit_range_low_;
  Double_t sigma_init = range_width * 0.01;
  Double_t peak_height =
      working_hist_->GetBinContent(working_hist_->GetMaximumBin());
  Double_t bkg_estimate = EstimateBackground();

  Double_t hist_xmin = working_hist_->GetXaxis()->GetXmin();
  Double_t hist_xmax = working_hist_->GetXaxis()->GetXmax();
  x_ = new RooRealVar("x", "x", hist_xmin, hist_xmax);
  x_->setRange(kFitRangeName, fit_range_low_, fit_range_high_);
  RegisterOwned(x_);

  for (Int_t pi = 0; pi < 2; pi++) {
    const PeakFitResult &cp = constrained_peaks.peaks[pi];
    BuildPeak(pi, cp.mu, cp.sigma, peak_height, range_width);
  }
  BuildPeak(2, mu3_init, sigma_init, peak_height, range_width);
  BuildBackground(bkg_estimate, peak_height, range_width);
  BuildTotalModel();

  data_hist_ = new RooDataHist("data_hist", "data_hist", RooArgList(*x_),
                                working_hist_);
  x_->setRange(fit_range_low_, fit_range_high_);

  for (Int_t pi = 0; pi < 2; pi++) {
    const PeakFitResult &cp = constrained_peaks.peaks[pi];
    RooFitPeakModel &p = peaks_[pi];
    Double_t cga = cp.gaus_amplitude;
    p.mu->setVal(cp.mu);
    p.mu->setConstant(kTRUE);
    p.sigma->setVal(cp.sigma);
    p.sigma->setConstant(kTRUE);
    p.gaus_yield->setVal(cga);
    p.gaus_yield->setRange(0, peak_height * range_width * 10.0);
    p.gaus_yield->setConstant(kFALSE);
    p.ratio_step->setVal(cp.step_amplitude / cga);
    p.ratio_step->setConstant(kTRUE);
    p.ratio_low_exp->setVal(cp.low_exp_tail_amplitude / cga);
    p.ratio_low_exp->setConstant(kTRUE);
    p.tau_low_exp->setVal(cp.low_exp_tail_decay > 0 ? cp.low_exp_tail_decay
                                                      : 1.0);
    p.tau_low_exp->setConstant(kTRUE);
    p.ratio_low_lin->setVal(cp.low_lin_tail_amplitude / cga);
    p.ratio_low_lin->setConstant(kTRUE);
    p.slope_low_lin->setVal(cp.low_lin_tail_slope);
    p.slope_low_lin->setConstant(kTRUE);
    p.ratio_high_exp->setVal(cp.high_exp_tail_amplitude / cga);
    p.ratio_high_exp->setConstant(kTRUE);
    p.tau_high_exp->setVal(cp.high_exp_tail_decay > 0 ? cp.high_exp_tail_decay
                                                         : 1.0);
    p.tau_high_exp->setConstant(kTRUE);
  }
  ConfigureComponentFlagsForPeak(2);

  Bool_t fit_valid = kFALSE;
  Double_t final_chi2 = 0;
  Int_t final_ndof = 0;

  if (interactive_) {
    if (LoadInteractiveParams(input_name, peak_name)) {
      RooFitResult *refit = RunFit(kTRUE);
      final_chi2 = ComputeReducedChi2(refit, final_ndof);
      fit_valid = kTRUE;
      delete refit;
    } else {
      Bool_t was_batch = gROOT->IsBatch();
      gROOT->SetBatch(kFALSE);
      if (LaunchInteractiveRooFitEditor(working_hist_, total_pdf_, x_,
                                         data_hist_, &peaks_, &bkg_,
                                         fit_range_low_, fit_range_high_,
                                         peak_name + " / " + input_name)) {
        final_chi2 = ComputeReducedChi2(nullptr, final_ndof);
        SaveInteractiveParams(input_name, peak_name);
        fit_valid = kTRUE;
      }
      gROOT->SetBatch(was_batch);
    }
  } else {
    FixComponent(2, "step");
    FixComponent(2, "low_exp");
    FixComponent(2, "low_lin");
    FixComponent(2, "high_exp");

    RooFitResult *initial_fit = RunFit(kTRUE);
    if (!initial_fit || initial_fit->status() != 0) {
      std::cout << "ERROR: Initial triple peak fit failed" << std::endl;
      delete initial_fit;
      return results;
    }

    Int_t tmp_ndof = 0;
    Double_t best_chi2 = ComputeReducedChi2(initial_fit, tmp_ndof);
    std::cout << "Initial chi2/ndf = " << best_chi2 << std::endl;
    delete initial_fit;

    std::vector<Double_t> best_vals;
    std::vector<Double_t> best_errs;
    std::vector<Bool_t> best_const;
    SnapshotParams(best_vals, best_errs, best_const);

    TestLowSideGroup(2, best_chi2, best_vals, best_errs, best_const);
    TestHighTailIndependent(2, best_chi2, best_vals, best_errs, best_const);

    std::cout << "Final fit with selected components..." << std::endl;
    RestoreParams(best_vals, best_errs, best_const);
    RooFitResult *final_fit = RunFit(kFALSE);
    if (final_fit && final_fit->status() == 0) {
      final_chi2 = ComputeReducedChi2(final_fit, final_ndof);
      fit_valid = kTRUE;
      std::cout << "Triple peak fit converged successfully" << std::endl;
      std::cout << "Final chi2/ndf = " << final_chi2 << std::endl;
    } else {
      std::cout << "ERROR: Triple peak fit failed to converge" << std::endl;
    }
    delete final_fit;
  }

  if (fit_valid) {
    SortPeaksByMu(3);
    TString chi2label = Form("#chi^{2}/ndf = %.3f", final_chi2);
    PlotFitTriplePeak(input_name, peak_name, chi2label);

    results.peaks[0] = ExtractPeakResult(0);
    results.peaks[1] = ExtractPeakResult(1);
    results.peaks[2] = ExtractPeakResult(2);
    results.bkg_constant = bkg_.bkg_yield->getVal();
    results.bkg_constant_error = bkg_.bkg_yield->getError();
    results.lin_bkg_slope = bkg_.bkg_slope->getVal();
    results.lin_bkg_slope_error = bkg_.bkg_slope->getError();
    results.reduced_chi2 = final_chi2;
    results.valid = kTRUE;
  }

  return results;
}
