#ifndef ROOFITUTILS_H
#define ROOFITUTILS_H

#include "FittingUtils.hpp"
#include "PlottingUtils.hpp"
#include "RooFitPhotopeakPdfs.hpp"

#include <RooAbsData.h>
#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooAddPdf.h>
#include <RooArgList.h>
#include <RooArgSet.h>
#include <RooCategory.h>
#include <RooCmdArg.h>
#include <RooDataSet.h>
#include <RooFitResult.h>
#include <RooFormulaVar.h>
#include <RooGaussian.h>
#include <RooGenericPdf.h>
#include <RooGlobalFunc.h>
#include <RooMsgService.h>
#include <RooPolynomial.h>
#include <RooRealVar.h>
#include <RooSimultaneous.h>

#include <TBranch.h>
#include <TGraph.h>
#include <TH1.h>
#include <TLeaf.h>
#include <TMath.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <vector>

/**
 * @brief The best RooFit evaluation backend this build supports.
 *
 * Defaults to CPU, where evaluation is still batched through `doEval()`.
 * Compiling with `-DAU_ROOFIT_BACKEND_CUDA=1`, against a ROOT built with CUDA
 * support, switches this to the GPU backend.
 *
 * @return A `RooCmdArg` to hand to `fitTo()`.
 *
 * @warning A downstream project must define the same macro as the prebuilt
 *          library, or its call sites will disagree with the library about
 *          which backend is in use.
 */
inline RooCmdArg BestAvailableBackend() {
#if defined(AU_ROOFIT_BACKEND_CUDA) && AU_ROOFIT_BACKEND_CUDA
  return RooFit::EvalBackend::Cuda();
#else
  return RooFit::EvalBackend::Cpu();
#endif
}

/**
 * @brief Factories for the custom photopeak PDF components.
 *
 * Thin wrappers over the classes in RooFitPhotopeakPdfs.hpp, for building a
 * model by hand instead of going through RooFitUtils.
 *
 * @note Every factory returns a heap-allocated PDF that the caller owns. The
 *       PDF holds references to the variables passed in, so those must outlive
 *       it.
 */
namespace RooFitFunctions {
/**
 * @brief Build a Gaussian component.
 * @param name  RooFit object name, unique within the workspace.
 * @param x     Observable.
 * @param mu    Peak centroid.
 * @param sigma Gaussian resolution.
 * @return A newly allocated PDF; the caller owns it.
 */
RooAbsPdf *MakeGaussian(const TString &name, RooRealVar &x, RooRealVar &mu,
                        RooRealVar &sigma);
/**
 * @brief Build a StepShelf component.
 * @param name  RooFit object name, unique within the workspace.
 * @param x     Observable.
 * @param mu    Peak centroid.
 * @param sigma Gaussian resolution.
 * @return A newly allocated PDF; the caller owns it.
 */
RooAbsPdf *MakeStepShelf(const TString &name, RooRealVar &x, RooRealVar &mu,
                         RooRealVar &sigma);
/**
 * @brief Build a LowExpTail component.
 * @param name  RooFit object name, unique within the workspace.
 * @param x     Observable.
 * @param mu    Peak centroid.
 * @param sigma Gaussian resolution.
 * @param tau_ratio Decay constant in units of @p sigma.
 * @return A newly allocated PDF; the caller owns it.
 */
RooAbsPdf *MakeLowExpTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                          RooRealVar &sigma, RooRealVar &tau_ratio);
/**
 * @brief Build a LowLinTail component.
 * @param name  RooFit object name, unique within the workspace.
 * @param x     Observable.
 * @param mu    Peak centroid.
 * @param sigma Gaussian resolution.
 * @param slope Slope of the linear tail factor.
 * @return A newly allocated PDF; the caller owns it.
 */
RooAbsPdf *MakeLowLinTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                          RooRealVar &sigma, RooRealVar &slope);
/**
 * @brief Build a HighExpTail component.
 * @param name  RooFit object name, unique within the workspace.
 * @param x     Observable.
 * @param mu    Peak centroid.
 * @param sigma Gaussian resolution.
 * @param tau_ratio Decay constant in units of @p sigma.
 * @return A newly allocated PDF; the caller owns it.
 */
RooAbsPdf *MakeHighExpTail(const TString &name, RooRealVar &x, RooRealVar &mu,
                           RooRealVar &sigma, RooRealVar &tau_ratio);
/**
 * @brief Build a linear background component.
 * @param name  RooFit object name, unique within the workspace.
 * @param x     Observable.
 * @param slope Background slope; fix it at zero for a flat background.
 * @return A newly allocated PDF; the caller owns it.
 */
RooAbsPdf *MakeLinearBackground(const TString &name, RooRealVar &x,
                                RooRealVar &slope);
} // namespace RooFitFunctions

/**
 * @brief Every RooFit object making up one peak.
 *
 * Held so a fit can reach in and fix, free or inspect an individual parameter.
 * Component amplitudes are `RooFormulaVar` yields expressed as ratios of the
 * Gaussian yield, which is what lets a constrained fit carry a shape between
 * peaks whose scales differ.
 *
 * @note Non-owning pointers. The RooFitUtils instance that built them owns
 *       them and destroys them with itself.
 */
struct RooFitPeakModel {
  RooRealVar *mu = nullptr;
  RooRealVar *sigma = nullptr;
  RooRealVar *gaus_yield = nullptr;
  RooRealVar *ratio_step = nullptr;
  RooRealVar *ratio_low_exp = nullptr;
  RooRealVar *tau_ratio_low_exp = nullptr;
  RooRealVar *ratio_low_lin = nullptr;
  RooRealVar *slope_low_lin = nullptr;
  RooRealVar *ratio_high_exp = nullptr;
  RooRealVar *tau_ratio_high_exp = nullptr;

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

/**
 * @brief The RooFit objects making up one channel's background.
 * @note Non-owning pointers, as with RooFitPeakModel.
 */
struct RooFitBackgroundModel {
  RooRealVar *bkg_yield = nullptr;
  RooRealVar *bkg_slope = nullptr;
  RooAbsPdf *bkg_pdf = nullptr;
};

/**
 * @brief One channel of a simultaneous fit: its data, model and locks.
 *
 * A channel is one spectrum with its own events, fit range, peaks and
 * background, sharing a single observable with the others. Populated by
 * RooFitUtils::AddChannel().
 */
struct RooFitChannelConfig {
  TString name;
  TH1 *hist;
  std::vector<Double_t> events;
  Float_t fit_range_low;
  Float_t fit_range_high;
  Float_t display_bin_width_kev;
  Int_t num_peaks;
  std::vector<Double_t> mu_inits;
  std::vector<Bool_t> mu_fixed;
  std::vector<Bool_t> use_step_per_peak;
  Bool_t bkg_yield_fixed = kFALSE;
  Bool_t bkg_slope_fixed = kFALSE;
  Bool_t lock_shape_after_seed = kFALSE;
  // Per-peak shape-lock override. When non-empty and lock_shape_after_seed is
  // true, each entry controls whether THAT peak's shape is locked. When empty,
  // all peaks are locked (backward-compatible with lock_shape_after_seed).
  std::vector<Bool_t> shape_lock_per_peak;
  Bool_t use_flat_background;
  Bool_t use_step;
  Bool_t use_low_exp_tail;
  Bool_t use_low_lin_tail;
  Bool_t use_high_exp_tail;
};

/**
 * @brief A tie between two parameters, fitted as one degree of freedom.
 *
 * Created by RooFitUtils::LinkParameter() and
 * RooFitUtils::LinkPeakShape().
 */
struct RooFitParamLink {
  TString target_channel;
  TString target_param;
  TString source_channel;
  TString source_param;
};

/**
 * @brief A known spacing between two peaks, imposed as a Gaussian penalty.
 *
 * Applies a penalty on `mu_hi - mu_lo` instead of fixing either centroid.
 *
 * @note Locking two centroids and constraining their separation are different
 *       statements. Locking pins the pair to an absolute energy scale as well;
 *       this leaves the pair free to slide together and forbids only the
 *       *stretch* between them.
 *
 * @note #sigma is the uncertainty on #delta, which keeps the constraint a
 *       measurement rather than a hard equality. Set it to the literature error
 *       on the spacing.
 */
struct RooFitSeparationConstraint {
  TString channel;    ///< Channel the two peaks belong to.
  Int_t peak_hi = 1;  ///< Index of the upper peak.
  Int_t peak_lo = 0;  ///< Index of the lower peak.
  Double_t delta = 0; ///< Known separation, in the observable's units.
  Double_t sigma = 0; ///< Uncertainty on @ref delta; must be positive.
};

/**
 * @brief Unbinned extended maximum-likelihood photopeak fitting on RooFit.
 *
 * The counterpart to FittingUtils, which fits the same model by binned
 * chi-squared. Both take the same constructor arguments, expose the same `Set*`
 * component flags, and return the same FitResult.
 *
 * Two modes, selected by which constructor is used:
 *
 * - **Single channel** — construct with events and a range, then call a `Fit*`
 *   method.
 * - **Simultaneous** — default-construct, AddChannel() each spectrum,
 * optionally LinkParameter() / LinkPeakShape() / ConstrainPeakSeparation() /
 *   SeedChannel(), then FitSimultaneous() for one joint fit across all of them.
 *
 * Input is an event-level `std::vector<Double_t>`. The bin width is used only
 * for display and for the post-hoc reduced chi-squared that drives component
 * pruning — the fit itself never bins.
 *
 * @warning Yields, not heights. The amplitude fields of PeakFitResult carry
 *          RooFit event counts here, where FittingUtils puts peak heights.
 *          Nothing in the type distinguishes them, so a result moved between
 *          backends is silently wrong by the normalisation. Only the ratio
 *          `component_amplitude / gaus_amplitude` is the same quantity in
 *          both.
 *
 * @note `RooRealVar::enableSilentClipping()` is on, so observable values
 *       outside the current range are clipped rather than throwing.
 */
class RooFitUtils {
private:
  TH1 *working_hist_;
  std::vector<Double_t> events_;
  Float_t fit_range_low_;
  Float_t fit_range_high_;
  Float_t display_bin_width_kev_;

  Bool_t use_flat_background_;
  Bool_t use_step_;
  Bool_t use_low_exp_tail_;
  Bool_t use_low_lin_tail_;
  Bool_t use_high_exp_tail_;

  Bool_t use_manual_init_;
  Bool_t interactive_;
  Bool_t fit_debug_;
  // See SetRefitAfterLoad(). Default kFALSE preserves the historic behaviour.
  Bool_t refit_after_load_ = kFALSE;
  // Upper bound on the exponential tail decay ratios tau/sigma. See
  // SetTailRatioMax().
  Double_t tail_ratio_max_ = 100.0;
  std::vector<Double_t> manual_params_;

  RooRealVar *x_;
  RooDataSet *unbinned_data_;
  RooAddPdf *total_pdf_;
  Int_t num_peaks_;

  std::vector<RooFitPeakModel> peaks_;
  RooFitBackgroundModel bkg_;
  std::vector<RooAbsArg *> owned_args_;

  std::vector<RooFitChannelConfig> sim_channels_;
  std::vector<RooFitParamLink> sim_links_;
  std::vector<RooFitSeparationConstraint> sim_sep_constraints_;
  RooArgSet sim_constraint_set_;
  std::map<TString, FitResult> sim_seeds_;
  std::map<TString, std::vector<RooFitPeakModel>> sim_channel_peaks_;
  std::map<TString, RooFitBackgroundModel> sim_channel_bkg_;
  std::map<TString, RooAbsPdf *> sim_channel_pdfs_;
  std::map<TString, RooDataSet *> sim_channel_data_;
  std::map<TString, TString> sim_channel_range_names_;
  RooCategory *sim_category_;
  RooSimultaneous *sim_pdf_;
  RooDataSet *sim_combined_data_;
  Bool_t sim_mode_;

  void InitState();
  void BuildDisplayHistogram();
  void BuildUnbinnedData();
  static RooDataSet *BuildUnbinnedDataFrom(const std::vector<Double_t> &events,
                                           RooRealVar *x);
  RooRealVar *ResolveOrCreate(const TString &channel, const TString &param_name,
                              std::map<TString, RooRealVar *> &registry,
                              Double_t init_val, Double_t lo, Double_t hi);
  Bool_t BuildChannelModel(const RooFitChannelConfig &cfg,
                           std::map<TString, RooRealVar *> &registry);
  void ApplySeedToChannel(const TString &channel);
  void ApplyChannelMuLocks();
  void ApplyChannelBkgLocks();
  void ApplyChannelShapeLocks();
  void SaveSimInteractiveParams(const TString &input_name,
                                const TString &base_label);
  Bool_t LoadSimInteractiveParams(const TString &input_name,
                                  const TString &base_label);
  TString ParamFullName(const TString &channel, const TString &param);
  TString SourceForTarget(const TString &target);
  void BuildSeparationConstraints();
  Double_t ComputeChannelChi2(const TString &channel,
                              const std::vector<RooFitPeakModel> &peaks,
                              const RooFitBackgroundModel &bkg, Int_t &ndof);
  void PlotChannel(const TString &channel, Int_t num_peaks,
                   const std::vector<RooFitPeakModel> &peaks,
                   const RooFitBackgroundModel &bkg, const TString &input_name,
                   const TString &base_label, const TString &chi2_label);
  PeakFitResult ExtractPeakResultFor(const RooFitPeakModel &p);

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

  // Extract per-parameter diagnostics (value, error, limits, near-limit flags)
  // from all RooRealVar objects in a channel's peak and background models.
  std::vector<FitParameterDiagnostic>
  ExtractParameterDiagnostics(const TString &channel);

  // Single-channel variant for non-simultaneous fits.
  std::vector<FitParameterDiagnostic> ExtractParameterDiagnosticsSingle();

  void SaveInteractiveParams(const TString &input_name,
                             const TString &peak_name);
  Bool_t LoadInteractiveParams(const TString &input_name,
                               const TString &peak_name);
  void AdoptSavedRange(const TString &input_name, const TString &peak_name);
  void AdoptSavedSimRange(const TString &input_name, const TString &base_label);

  void SortPeaksByMu(Int_t num_peaks);
  void AppendPeakGraphs(std::vector<TGraph *> &components, Int_t peak_idx,
                        Style_t line_style, RooAbsPdf *background_pdf,
                        Double_t bkg_yield_val, Int_t npts, Double_t x_step,
                        Double_t bin_width);

  void RegisterOwned(RooAbsArg *arg);

public:
  /// @brief Name of the RooFit range this class fits over.
  static constexpr const char *kFitRangeName = "fitrange";

  /**
   * @brief Construct in simultaneous mode.
   *
   * Register spectra with AddChannel(), then call FitSimultaneous(). The
   * single-channel `Fit*` methods do not apply to an instance built this way.
   */
  RooFitUtils();
  /**
   * @brief Construct in single-channel mode.
   *
   * @param events                Event-level observable values. Copied.
   * @param fit_range_low         Lower fit bound.
   * @param fit_range_high        Upper fit bound.
   * @param display_bin_width_kev Bin width for the display histogram and the
   *                              post-hoc reduced chi-squared. Does not affect
   *                              the likelihood, which is unbinned.
   * @param use_flat_background   `kTRUE` for a constant background, `kFALSE`
   *                              (default) for a linear one.
   * @param use_step              Enable the step shelf.
   * @param use_low_exp_tail      Enable the low-energy exponential tail.
   * @param use_low_lin_tail      Enable the low-energy linear tail.
   * @param use_high_exp_tail     Enable the high-energy exponential tail.
   */
  RooFitUtils(const std::vector<Double_t> &events, Float_t fit_range_low,
              Float_t fit_range_high, Float_t display_bin_width_kev,
              Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
              Bool_t use_low_exp_tail = kFALSE,
              Bool_t use_low_lin_tail = kFALSE,
              Bool_t use_high_exp_tail = kFALSE);
  /// @brief Destroys every RooFit object this instance created.
  ~RooFitUtils();

  /**
   * @brief Read a branch into the event vector the constructor expects.
   * @param tree        Tree to read. Must not be null.
   * @param branch_name Branch holding one observable value per entry.
   * @return The values, in entry order; empty if the branch is absent.
   */
  static std::vector<Double_t> LoadEventsFromTree(TTree *tree,
                                                  const TString &branch_name);
  /**
   * @brief Bin events into a display histogram.
   * @param events                Event-level values.
   * @param fit_range_low         Lower bound.
   * @param fit_range_high        Upper bound.
   * @param display_bin_width_kev Bin width.
   * @return A newly allocated histogram; the caller owns it.
   */
  static TH1F *BuildDisplayHistogramFrom(const std::vector<Double_t> &events,
                                         Float_t fit_range_low,
                                         Float_t fit_range_high,
                                         Float_t display_bin_width_kev);
  /**
   * @brief Rebin an existing display histogram in place.
   *
   * What the interactive editor uses when the range slider moves, so binning
   * follows the zoom without reallocating.
   *
   * @param hist                  Histogram to refill. Must not be null.
   * @param events                Event-level values.
   * @param fit_range_low         New lower bound.
   * @param fit_range_high        New upper bound.
   * @param display_bin_width_kev Bin width.
   */
  static void RefillDisplayHistogram(TH1 *hist,
                                     const std::vector<Double_t> &events,
                                     Float_t fit_range_low,
                                     Float_t fit_range_high,
                                     Float_t display_bin_width_kev);

  /// @brief Choose the background shape.
  /// @param use_flat_background `kTRUE` for constant, `kFALSE` for linear.
  void SetBackgroundModel(Bool_t use_flat_background) {
    use_flat_background_ = use_flat_background;
  }
  /// @brief Offer the step shelf to the component search.
  /// @param use_step `kTRUE` to enable.
  void SetStep(Bool_t use_step = kTRUE) { use_step_ = use_step; }
  /// @brief Offer the low-energy exponential tail.
  /// @param use_low_exp_tail `kTRUE` to enable.
  void SetLowExpTail(Bool_t use_low_exp_tail = kTRUE) {
    use_low_exp_tail_ = use_low_exp_tail;
  }
  /// @brief Offer the low-energy linear tail.
  /// @param use_low_lin_tail `kTRUE` to enable.
  void SetLowLinTail(Bool_t use_low_lin_tail = kTRUE) {
    use_low_lin_tail_ = use_low_lin_tail;
  }
  /// @brief Offer the high-energy exponential tail.
  /// @param use_high_exp_tail `kTRUE` to enable.
  void SetHighExpTail(Bool_t use_high_exp_tail = kTRUE) {
    use_high_exp_tail_ = use_high_exp_tail;
  }
  /**
   * @brief Open the GUI editor after the automated fit.
   *
   * Accepted parameters are saved to
   * `<plots_base>/fits/<peak_name>_<input_name>.roofits` and reloaded on later
   * runs, so the manual step happens once. Delete the file to redo it.
   *
   * @param interactive `kTRUE` to enable.
   *
   * @note The `.roofits` extension keeps these distinct from the `.fits` files
   *       written by FittingUtils, so both backends can be used on one input.
   */
  void SetInteractive(Bool_t interactive = kTRUE) {
    interactive_ = interactive;
  }
  /**
   * @brief Un-suppress RooFit evaluation errors during a simultaneous fit.
   *
   * Also prints the seed NLL, so an invalid-NLL failure names the offending
   * PDF instead of failing anonymously.
   *
   * @param fit_debug `kTRUE` to enable. Off by default; set it before calling
   *                  FitSimultaneous().
   */
  void SetFitDebug(Bool_t fit_debug = kTRUE) { fit_debug_ = fit_debug; }

  /**
   * @brief Treat saved interactive parameters as a seed, not as the answer.
   *
   * Interactive simultaneous fits have two historic modes. With a saved
   * `.simroofits` file present, the parameters were loaded and **adopted as the
   * result**, with no minimisation at all; without one, the GUI editor opened.
   *
   * The first mode makes a re-run silently replay a previous session's answer:
   * it still prints a chi-squared, redraws every plot and writes fresh output,
   * while the fitted values cannot respond to anything upstream that has
   * changed. Enabling this makes the saved parameters seed a real
   * minimisation instead.
   *
   * @param refit `kTRUE` to re-minimise after loading.
   *
   * @warning Defaults to `kFALSE` so existing callers keep the replay
   *          behaviour. Turn it on for any fit whose inputs may change.
   */
  void SetRefitAfterLoad(Bool_t refit = kTRUE) { refit_after_load_ = refit; }

  /**
   * @brief Cap the dimensionless tail decay ratio `tau / sigma`.
   *
   * Applies to both `LowExpTailRatio` and `HighExpTailRatio`.
   *
   * @param ratio_max Upper bound, in units of sigma.
   *
   * @warning The default of 100 is rarely what you want. On a typical fit
   *          window it makes `tau` longer than the window itself, over which
   *          the exponential is flat — so the "tail" becomes a constant
   *          pedestal, degenerate with the background yield, and absorbs
   *          background rather than describing a tail. Set a physically
   *          motivated value so the component either describes a real tail or
   *          fits to nothing.
   */
  void SetTailRatioMax(Double_t ratio_max) { tail_ratio_max_ = ratio_max; }

  /**
   * @brief Supply explicit starting values for every parameter.
   * @param params Initial values, in the model's parameter order.
   */
  void SetManualParameters(const std::vector<Double_t> &params);
  /**
   * @brief Override one starting value.
   * @param index Parameter index in the active model.
   * @param value Initial value.
   */
  void SetManualParameter(Int_t index, Double_t value);
  /// @brief Discard manual starting values and return to automatic guesses.
  void ClearManualParameters() {
    use_manual_init_ = kFALSE;
    manual_params_.clear();
  }

  /**
   * @brief Draw and save the single-peak fit with its residual panel.
   * @param input_name Input identifier, used in the output filename.
   * @param peak_name  Peak identifier, used in the output filename.
   * @param label      Optional annotation drawn on the plot.
   */
  void PlotFitSinglePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  /**
   * @brief Draw and save the double-peak fit with its residual panel.
   * @param input_name Input identifier, used in the output filename.
   * @param peak_name  Peak identifier, used in the output filename.
   * @param label      Optional annotation drawn on the plot.
   */
  void PlotFitDoublePeak(const TString input_name, const TString peak_name,
                         const TString label = "");
  /**
   * @brief Draw and save the triple-peak fit with its residual panel.
   * @param input_name Input identifier, used in the output filename.
   * @param peak_name  Peak identifier, used in the output filename.
   * @param label      Optional annotation drawn on the plot.
   */
  void PlotFitTriplePeak(const TString input_name, const TString peak_name,
                         const TString label = "");

  /**
   * @brief Fit one peak, pruning components that do not earn their place.
   *
   * Low-side components are enabled as a group and pruned individually; the
   * high-energy tail is tested separately. A component is kept only if it
   * improves the reduced chi-squared.
   *
   * @param input_name Input identifier, used for saved parameters and output
   *                   filenames.
   * @param peak_name  Peak identifier, used the same way.
   * @return The fitted result; check FitResult::valid.
   */
  FitResult FitSinglePeak(const TString input_name, const TString peak_name);
  /**
   * @brief Fit two peaks from centroid guesses.
   *
   * @param input_name Input identifier.
   * @param peak_name  Peak identifier.
   * @param mu1_init   Starting centroid for the first peak.
   * @param mu2_init   Starting centroid for the second.
   * @param link_sigma When `kTRUE`, ties both peaks to a single Gaussian width:
   *                   the second peak's PDFs are rebuilt to reference the
   *                   first's sigma and `Sigma2` is frozen. Right for close,
   *                   same-detector lines such as the Pb K-beta group, where a
   *                   free second width would otherwise trade against the
   *                   background.
   * @return The fitted result, with peaks sorted by centroid.
   */
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init,
                          Bool_t link_sigma = kFALSE);
  /**
   * @brief Fit two peaks with the first constrained by an earlier fit.
   * @param input_name       Input identifier.
   * @param peak_name        Peak identifier.
   * @param constrained_peak Previously fitted peak supplying the constraint.
   * @param mu2_init         Starting centroid for the free peak.
   * @return The fitted result, with peaks sorted by centroid.
   */
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          const PeakFitResult &constrained_peak,
                          Double_t mu2_init);
  /**
   * @brief Fit three peaks with the first two constrained by an earlier fit.
   * @param input_name        Input identifier.
   * @param peak_name         Peak identifier.
   * @param constrained_peaks Previous result supplying the two constraints.
   * @param mu3_init          Starting centroid for the free peak.
   * @return The fitted result, with peaks sorted by centroid.
   */
  FitResult FitTriplePeak(const TString input_name, const TString peak_name,
                          const FitResult &constrained_peaks,
                          Double_t mu3_init);

  /**
   * @brief Register one spectrum as a channel of the simultaneous fit.
   *
   * @param name                  Channel name, unique within the fit; used to
   *                              address it from the linking methods.
   * @param events                Event-level observable values. Copied.
   * @param fit_range_low         Lower fit bound for this channel.
   * @param fit_range_high        Upper fit bound for this channel.
   * @param display_bin_width_kev Display bin width; does not affect the
   *                              likelihood.
   * @param num_peaks             Peaks in this channel.
   * @param mu_inits              Starting centroid per peak; @p num_peaks
   *                              entries.
   * @param use_flat_background   `kTRUE` for constant, `kFALSE` for linear.
   * @param use_step              Enable the step shelf for this channel.
   * @param use_low_exp_tail      Enable the low-energy exponential tail.
   * @param use_low_lin_tail      Enable the low-energy linear tail.
   * @param use_high_exp_tail     Enable the high-energy exponential tail.
   * @param mu_fixed              Per-peak; where `kTRUE`, that centroid is held
   *                              at its @p mu_inits value instead of floating.
   *                              Empty means all float.
   * @param bkg_yield_fixed       Hold this channel's background yield constant.
   * @param bkg_slope_fixed       Hold this channel's background slope constant.
   * @param lock_shape_after_seed After seeding, fix sigma, the Gaussian yield
   *                              and every tail and step shape parameter, so
   *                              only peak positions and normalisations float.
   * @param use_step_per_peak     Per-peak override of @p use_step, for a step
   *                              on some peaks but not others sharing a
   *                              channel. Empty means @p use_step applies to
   *                              all.
   * @param shape_lock_per_peak   Per-peak refinement of
   *                              @p lock_shape_after_seed: entry `i` controls
   *                              whether peak `i`'s shape locks. Empty locks
   *                              every peak.
   *
   * @note Simultaneous mode only — the instance must have been
   *       default-constructed.
   */
  void AddChannel(
      const TString &name, const std::vector<Double_t> &events,
      Float_t fit_range_low, Float_t fit_range_high,
      Float_t display_bin_width_kev, Int_t num_peaks,
      const std::vector<Double_t> &mu_inits,
      Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
      Bool_t use_low_exp_tail = kFALSE, Bool_t use_low_lin_tail = kFALSE,
      Bool_t use_high_exp_tail = kFALSE,
      const std::vector<Bool_t> &mu_fixed = std::vector<Bool_t>(),
      Bool_t bkg_yield_fixed = kFALSE, Bool_t bkg_slope_fixed = kFALSE,
      Bool_t lock_shape_after_seed = kFALSE,
      const std::vector<Bool_t> &use_step_per_peak = std::vector<Bool_t>(),
      const std::vector<Bool_t> &shape_lock_per_peak = std::vector<Bool_t>());
  /**
   * @brief Tie one parameter to another, fitting them as one degree of freedom.
   * @param target Parameter to tie, as `"<channel>:<parameter>"`.
   * @param source Parameter it follows, in the same form.
   */
  void LinkParameter(const TString &target, const TString &source);
  /**
   * @brief Constrain the spacing between two peaks to a known value.
   *
   * Imposes a Gaussian penalty on `mu[peak_hi] - mu[peak_lo]`.
   *
   * @param channel Channel holding both peaks.
   * @param peak_hi Index of the upper peak.
   * @param peak_lo Index of the lower peak.
   * @param delta   Known separation, in the observable's units.
   * @param sigma   Uncertainty on @p delta; set it to the literature error, and
   *                keep it positive.
   *
   * @note This is not the same as locking both centroids — see
   *       RooFitSeparationConstraint. Call it before FitSimultaneous().
   */
  void ConstrainPeakSeparation(const TString &channel, Int_t peak_hi,
                               Int_t peak_lo, Double_t delta, Double_t sigma);
  /**
   * @brief Tie a peak's whole shape to another peak's.
   *
   * Links sigma and every tail and step shape parameter at once, leaving
   * position and normalisation free.
   *
   * @param target_channel Channel holding the peak to tie.
   * @param target_peak    Index of that peak.
   * @param source_channel Channel holding the peak it follows.
   * @param source_peak    Index of that peak.
   */
  void LinkPeakShape(const TString &target_channel, Int_t target_peak,
                     const TString &source_channel, Int_t source_peak);
  /**
   * @brief Use a prior single-channel fit as starting values for a channel.
   * @param channel_name Channel to seed.
   * @param result       Previous result supplying the starting values.
   */
  void SeedChannel(const TString &channel_name, const FitResult &result);
  /**
   * @brief Run one joint extended-likelihood fit across every channel.
   *
   * Builds the `RooSimultaneous`, applies the links and separation constraints,
   * and minimises. With SetInteractive() the simultaneous editor opens
   * afterwards, and accepted parameters are saved and reloaded as for the
   * single-channel editors.
   *
   * @param input_name Input identifier, used for saved parameters and output
   *                   filenames.
   * @param base_label Base label for the per-channel outputs.
   *
   * @return One FitResult per channel, in registration order. These are the
   *         only results that carry the Minuit diagnostics — `fit_status`,
   *         `cov_qual`, `edm`, `min_nll` and `parameter_diagnostics`.
   *
   * @note Check `fit_status == 0` and `cov_qual == 0` before trusting the
   *       errors, and inspect `parameter_diagnostics` for parameters resting
   *       against a limit.
   */
  std::vector<FitResult> FitSimultaneous(const TString &input_name,
                                         const TString &base_label);

  /**
   * @brief Export a converged channel to CSV for plotting elsewhere.
   *
   * Writes two tidy files that `pandas.read_csv` reads directly, using the same
   * `expectedEvents * getVal * bin_width` convention as the plotted curve:
   *
   * - `<base>_spectrum.csv` — one row per histogram bin: `energy_keV`,
   *   `data_counts`, `fit_total`, `fit_background`, `residual_pull`.
   * - `<base>_fitcurve.csv` — a smooth overlay of @p npts samples:
   * `energy_keV`, `fit_total`, `fit_background`, in counts per bin.
   *
   * @param channel  Channel to export.
   * @param csv_path Output path **base**; a trailing `.csv` is stripped before
   *                 the two suffixes are appended.
   * @param npts     Samples along the smooth overlay curve.
   *
   * @note Read-only: evaluates the current parameter state and changes no fit
   *       state. Call it after FitSimultaneous().
   */
  void DumpChannelCSV(const TString &channel, const TString &csv_path,
                      Int_t npts = 1000);
};

#endif
