#ifndef ROOFITUTILS_H
#define ROOFITUTILS_H

#include "FittingUtils.hpp"
#include "PlottingUtils.hpp"

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooDataHist.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>

#include <TH1.h>
#include <TString.h>
#include <map>
#include <vector>

namespace RooFitFunctions {
RooAbsPdf *MakeGaussian(const TString &name, RooRealVar &x, RooRealVar &mu,
                        RooRealVar &sigma);
RooAbsPdf *MakeStepShelf(const TString &name, RooRealVar &x, RooRealVar &mu,
                         RooRealVar &sigma);
RooAbsPdf *MakeLowExpTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                          RooRealVar &sigma, RooRealVar &tau);
RooAbsPdf *MakeLowLinTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                          RooRealVar &sigma, RooRealVar &slope);
RooAbsPdf *MakeHighExpTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                           RooRealVar &sigma, RooRealVar &tau);
RooAbsPdf *MakeLinearBackground(const TString &name, RooRealVar &x,
                                RooRealVar &slope);
RooAbsPdf *MakeFlatBackground(const TString &name, RooRealVar &x);
} // namespace RooFitFunctions

struct RooFitPeakModel {
  RooRealVar *mu = nullptr;
  RooRealVar *sigma = nullptr;
  RooRealVar *gaus_yield = nullptr;
  RooRealVar *ratio_step = nullptr;
  RooRealVar *ratio_low_exp = nullptr;
  RooRealVar *tau_low_exp = nullptr;
  RooRealVar *ratio_low_lin = nullptr;
  RooRealVar *slope_low_lin = nullptr;
  RooRealVar *ratio_high_exp = nullptr;
  RooRealVar *tau_high_exp = nullptr;

  RooAbsPdf *gauss_pdf = nullptr;
  RooAbsPdf *step_pdf = nullptr;
  RooAbsPdf *low_exp_pdf = nullptr;
  RooAbsPdf *low_lin_pdf = nullptr;
  RooAbsPdf *high_exp_pdf = nullptr;

  RooFormulaVar *step_yield = nullptr;
  RooFormulaVar *low_exp_yield = nullptr;
  RooFormulaVar *low_lin_yield = nullptr;
  RooFormulaVar *high_exp_yield = nullptr;
};

struct RooFitBackgroundModel {
  RooRealVar *bkg_yield = nullptr;
  RooRealVar *bkg_slope = nullptr;
  RooAbsPdf *bkg_pdf = nullptr;
};

class RooFitUtils {
private:
  TH1 *working_hist_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;

  Bool_t use_flat_background_;
  Bool_t use_step_;
  Bool_t use_low_exp_tail_;
  Bool_t use_low_lin_tail_;
  Bool_t use_high_exp_tail_;

  Bool_t use_manual_init_;
  Bool_t interactive_;
  std::vector<Double_t> manual_params_;

  RooRealVar *x_;
  RooDataHist *data_hist_;
  RooAddPdf *total_pdf_;
  Int_t num_peaks_;

  std::vector<RooFitPeakModel> peaks_;
  RooFitBackgroundModel bkg_;
  std::vector<RooAbsArg *> owned_args_;

  static constexpr const char *kFitRangeName = "fitrange";

  Double_t EstimateBackground();

  void BuildPeak(Int_t peak_idx, Double_t mu_init, Double_t sigma_init,
                 Double_t peak_height, Double_t range_width);
  void BuildBackground(Double_t bkg_estimate, Double_t peak_height,
                       Double_t range_width);
  void BuildTotalModel();
  void ConfigureComponentFlagsForPeak(Int_t peak_idx);

  void FixComponent(Int_t peak_idx, const TString &component);
  void ReleaseComponent(Int_t peak_idx, const TString &component);
  Bool_t ComponentIsActive(Int_t peak_idx, const TString &component);

  std::vector<RooRealVar *> CollectFloatingParams();
  std::vector<RooRealVar *> CollectAllParams();
  RooFitResult *RunFit(Bool_t quiet);
  Double_t ComputeReducedChi2(RooFitResult *fit_result, Int_t &ndof);

  void SnapshotParams(std::vector<Double_t> &vals, std::vector<Double_t> &errs,
                       std::vector<Bool_t> &consts);
  void RestoreParams(const std::vector<Double_t> &vals,
                      const std::vector<Double_t> &errs,
                      const std::vector<Bool_t> &consts);
  void TestLowSideGroup(Int_t peak_idx, Double_t &best_chi2,
                         std::vector<Double_t> &best_vals,
                         std::vector<Double_t> &best_errs,
                         std::vector<Bool_t> &best_const);
  void TestHighTailIndependent(Int_t peak_idx, Double_t &best_chi2,
                                 std::vector<Double_t> &best_vals,
                                 std::vector<Double_t> &best_errs,
                                 std::vector<Bool_t> &best_const);

  PeakFitResult ExtractPeakResult(Int_t peak_idx);

  void SaveInteractiveParams(const TString &input_name,
                             const TString &peak_name);
  Bool_t LoadInteractiveParams(const TString &input_name,
                               const TString &peak_name);

  void SortPeaksByMu(Int_t num_peaks);
  void AppendPeakGraphs(std::vector<TGraph *> &components, Int_t peak_idx,
                        Style_t line_style, RooAbsPdf *background_pdf,
                        Double_t bkg_yield_val, Int_t npts, Double_t x_step,
                        Double_t bin_width);

  void RegisterOwned(RooAbsArg *arg);

public:
  RooFitUtils(TH1 *working_hist, Float_t fit_range_low, Float_t fit_range_high,
              Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
              Bool_t use_low_exp_tail = kFALSE,
              Bool_t use_low_lin_tail = kFALSE,
              Bool_t use_high_exp_tail = kFALSE);
  ~RooFitUtils();

  void SetBackgroundModel(Bool_t use_flat_background) {
    use_flat_background_ = use_flat_background;
  }
  void SetStep(Bool_t use_step = kTRUE) { use_step_ = use_step; }
  void SetLowExpTail(Bool_t use_low_exp_tail = kTRUE) {
    use_low_exp_tail_ = use_low_exp_tail;
  }
  void SetLowLinTail(Bool_t use_low_lin_tail = kTRUE) {
    use_low_lin_tail_ = use_low_lin_tail;
  }
  void SetHighExpTail(Bool_t use_high_exp_tail = kTRUE) {
    use_high_exp_tail_ = use_high_exp_tail;
  }
  void SetInteractive(Bool_t interactive = kTRUE) {
    interactive_ = interactive;
  }

  void SetManualParameters(const std::vector<Double_t> &params);
  void SetManualParameter(Int_t index, Double_t value);
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  void PlotFitSinglePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitDoublePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  void PlotFitTriplePeak(const TString input_name, const TString peak_name,
                         const TString label = "");

  FitResult FitSinglePeak(const TString input_name, const TString peak_name);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init);
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakFitResult &constrained_peak,
                          Double_t mu2_init);
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks,
                          Double_t mu3_init);
};

#endif
