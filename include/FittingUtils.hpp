#ifndef FITTINGUTILS_H
#define FITTINGUTILS_H

#include "PlottingUtils.hpp"
#include <TCanvas.h>
#include <TF1.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TH1.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TSystem.h>
#include <TTree.h>
#include <fstream>
#include <iomanip>

/**
 * @file FittingUtils.hpp
 * @brief Binned chi-squared photopeak fitting on a `TF1` + Minuit2 backend.
 *
 * The counterpart to RooFitUtils, which fits the same model by unbinned
 * extended maximum likelihood. Both take the same constructor arguments, expose
 * the same `Set*` component flags, and return the same FitResult, so a call
 * site can swap one for the other. This backend is faster, because the model is
 * evaluated pointwise and never normalised; the RooFit backend is the one to
 * reach for when you need correct likelihood errors or simultaneous fits.
 *
 * Saved interactive parameters use the `.fits` extension here and `.roofits`
 * in the RooFit backend, so both can coexist for the same input.
 */

// Forward-declared so FittingUtils can launch the editor without pulling in the
// GUI headers. Documented at its real declaration in InteractiveFitEditor.hpp —
// a doxygen block here too would be merged with that one, doubling the @param
// list for a single function.
Bool_t LaunchInteractiveFitEditor(TH1 *hist, TF1 *fit_func, Double_t range_low,
                                  Double_t range_high, Int_t num_peaks = 1,
                                  const TString &info_label = "");

/**
 * @brief Raw model functions in `TF1` form: `f(Double_t *x, Double_t *par)`.
 *
 * Exposed so a caller can build a custom `TF1` from the same components the
 * class uses. Every function is un-normalised — amplitudes are peak heights in
 * counts, not yields.
 *
 * In the composite functions the step and tail amplitudes are expressed as
 * **ratios of the Gaussian amplitude** rather than absolute heights, which is
 * what keeps a constrained fit meaningful when the peak height changes. Tail
 * `ratio` parameters are decay constants in units of `sigma`, so
 * `tau = ratio * sigma`.
 */
namespace FittingFunctions {
/**
 * @brief Gaussian peak.
 * @param x   Observable; only `x[0]` is used.
 * @param par `[0]` centroid, `[1]` sigma, `[2]` peak height.
 * @return The un-normalised Gaussian value.
 */
Double_t Gaussian(Double_t *x, Double_t *par);
/**
 * @brief Straight-line background.
 * @param x   Observable; only `x[0]` is used.
 * @param par `[0]` constant, `[1]` slope.
 * @return `slope * x + constant`.
 */
Double_t LinearBackground(Double_t *x, Double_t *par);
/**
 * @brief Resolution-smeared step shelf below the peak.
 * @param x   Observable; only `x[0]` is used.
 * @param par `[0]` centroid, `[1]` sigma, `[2]` step amplitude.
 * @return The shelf value; zero when `sigma <= 0` or the denominator
 * underflows.
 */
Double_t Step(Double_t *x, Double_t *par);
/**
 * @brief Combined exponential and linear tail below the peak.
 * @param x   Observable; only `x[0]` is used.
 * @param par `[0]` centroid, `[1]` sigma, `[2]` exponential amplitude,
 *            `[3]` exponential decay in units of sigma, `[4]` linear
 *            amplitude, `[5]` linear slope.
 * @return The summed tail, cut off above the centroid by a complementary error
 *         function. Zero when `sigma <= 0` or both amplitudes are zero.
 */
Double_t LowTail(Double_t *x, Double_t *par);
/**
 * @brief Exponential tail above the peak, from pileup.
 * @param x   Observable; only `x[0]` is used.
 * @param par `[0]` centroid, `[1]` sigma, `[2]` amplitude, `[3]` decay constant
 *            in units of sigma.
 * @return The tail value; zero when `sigma <= 0`, the amplitude is zero, or the
 *         decay constant is non-positive.
 */
Double_t HighTail(Double_t *x, Double_t *par);
/**
 * @brief One peak with every optional component, plus a linear background.
 *
 * Twelve parameters: `[0]` centroid, `[1]` sigma, `[2]` Gaussian height,
 * `[3]` step amplitude, `[4]` low exponential tail amplitude, `[5]` its decay
 * in units of sigma, `[6]` low linear tail amplitude, `[7]` its slope,
 * `[8]` high exponential tail amplitude, `[9]` its decay in units of sigma,
 * `[10]` background constant, `[11]` background slope.
 *
 * @param x   Observable; only `x[0]` is used.
 * @param par The twelve parameters above.
 * @return The summed model value.
 *
 * @note Parameters 3, 4, 6 and 8 are **ratios of the Gaussian amplitude**, not
 *       absolute heights. Setting a ratio to zero disables that component.
 */
Double_t PeakFunction(Double_t *x, Double_t *par);
/**
 * @brief Two peaks sharing one linear background.
 *
 * Each peak occupies a ten-parameter block laid out as parameters 0 to 9 of
 * PeakFunction(): peak one at `[0..9]`, peak two at `[10..19]`, with the
 * background constant and slope appended at `[20]` and `[21]`.
 *
 * @param x   Observable; only `x[0]` is used.
 * @param par The 22 parameters above.
 * @return The summed model value.
 */
Double_t DoublePeakFunction(Double_t *x, Double_t *par);
/**
 * @brief Three peaks sharing one linear background.
 *
 * As DoublePeakFunction(), with a third ten-parameter block at `[20..29]` and
 * the background appended at `[30]` and `[31]`.
 *
 * @param x   Observable; only `x[0]` is used.
 * @param par The 32 parameters above.
 * @return The summed model value.
 */
Double_t TriplePeakFunction(Double_t *x, Double_t *par);
} // namespace FittingFunctions

/**
 * @brief Fitted parameters and errors for one peak.
 *
 * Every field defaults to `-1`, which is also what a failed fit leaves behind —
 * check FitResult::valid rather than testing individual members.
 *
 * @warning Amplitude semantics differ between the two backends, and nothing in
 *          the type says which you are holding. FittingUtils fills these with
 *          peak heights in counts; RooFitUtils fills them with RooFit yields
 *          (event counts). Mixing them up produces numbers that look plausible
 *          and are wrong by the normalisation. The ratio
 *          `component_amplitude / gaus_amplitude` is the one quantity that
 *          means the same thing in both, which is what makes constrained fits
 *          portable between them.
 */
struct PeakFitResult {
  Float_t mu = -1, mu_error = -1; ///< Centroid, in the histogram's x units.
  Float_t sigma = -1, sigma_error = -1; ///< Gaussian resolution.
  Float_t gaus_amplitude = -1, gaus_amplitude_error = -1; ///< Gaussian scale.
  Float_t step_amplitude = -1, step_amplitude_error = -1; ///< Step shelf scale.
  /// Low-energy exponential tail scale.
  Float_t low_exp_tail_amplitude = -1, low_exp_tail_amplitude_error = -1;
  /// Low-energy exponential decay constant, in units of #sigma.
  Float_t low_exp_tail_ratio = -1, low_exp_tail_ratio_error = -1;
  /// Low-energy linear tail scale.
  Float_t low_lin_tail_amplitude = -1, low_lin_tail_amplitude_error = -1;
  /// Slope of the linear tail factor.
  Float_t low_lin_tail_slope = -1, low_lin_tail_slope_error = -1;
  /// High-energy exponential tail scale.
  Float_t high_exp_tail_amplitude = -1, high_exp_tail_amplitude_error = -1;
  /// High-energy exponential decay constant, in units of #sigma.
  Float_t high_exp_tail_ratio = -1, high_exp_tail_ratio_error = -1;
};

/**
 * @brief One fitted parameter's value, error and proximity to its limits.
 *
 * A parameter resting against a bound usually means the model is being forced
 * somewhere it does not want to go, so the near-limit flags are worth checking
 * before trusting a converged fit.
 *
 * @note Only populated by RooFitUtils simultaneous fits.
 */
struct FitParameterDiagnostic {
  std::string name;           ///< Parameter name as RooFit knows it.
  Double_t value = 0;         ///< Fitted value.
  Double_t error = 0;         ///< Fitted uncertainty.
  Double_t lo = 0;            ///< Lower bound, if #has_limits.
  Double_t hi = 0;            ///< Upper bound, if #has_limits.
  Bool_t has_limits = kFALSE; ///< Whether the parameter was bounded.
  Bool_t near_lower = kFALSE; ///< Sits within a small margin of #lo.
  Bool_t near_upper = kFALSE; ///< Sits within a small margin of #hi.
  Bool_t near_limit = kFALSE; ///< Either of the two above.
};

/**
 * @brief Result of one fit: peaks, background, quality, and diagnostics.
 *
 * Returned by both backends. The diagnostic fields below the background are
 * only populated by RooFitUtils simultaneous fits; everything else leaves them
 * at their sentinel defaults, flagged by #has_fit_diagnostics.
 */
struct FitResult {
  std::vector<PeakFitResult> peaks; ///< One entry per peak; 1 to 3 supported.
  Float_t bkg_constant = -1, bkg_constant_error = -1;   ///< Background offset.
  Float_t lin_bkg_slope = -1, lin_bkg_slope_error = -1; ///< Background slope.
  Float_t reduced_chi2 = -1; ///< Chi-squared per degree of freedom.
  /// In the RooFit backend this is computed post hoc against the display
  /// binning, and drives component pruning; the fit itself is unbinned.
  Bool_t valid = kFALSE; ///< Whether the fit succeeded. Check this first.

  /// @name Minuit and RooFit diagnostics
  /// Populated only by RooFitUtils simultaneous fits, so that downstream code
  /// can tell a true minimum from a local or failed one. Every other fit leaves
  /// these at their sentinel defaults.
  /// @{
  Int_t fit_status = -999; ///< Minuit migrad status; `0` means success.
  /// Covariance quality: `0` Exact, `1` NotPositiveDefinite, `2` Approximate,
  /// `3` External. Anything but `0` means the errors are unreliable.
  Int_t cov_qual = -999;
  Double_t edm = -1;     ///< Estimated distance to the minimum.
  Double_t min_nll = -1; ///< Minimised negative log-likelihood.
  /// Whether a `RooFitResult` was available, and hence whether the fields in
  /// this group mean anything.
  Bool_t has_fit_diagnostics = kFALSE;
  Int_t n_invalid_nll = -1; ///< Invalid NLL evaluations; `-1` if unknown.
  /// One entry per fitted `RooRealVar`, peak and background alike.
  std::vector<FitParameterDiagnostic> parameter_diagnostics;
  /// @}
};

/**
 * @brief Binned chi-squared photopeak fitting over a `TH1`.
 *
 * Construct with the histogram, the fit range, and which optional components
 * the model should offer, then call one of the `Fit*` methods. Component
 * selection is not final: the fitter tests the low-side components as a group
 * and prunes those that do not improve the reduced chi-squared, and tests the
 * high-energy tail separately.
 *
 * @note Does not take ownership of the histogram passed to the constructor, but
 *       does own the `TF1` it builds. See GetFitFunction() and SetFitFunction()
 *       for the ownership caveat when substituting your own.
 */
class FittingUtils {
private:
  TF1 *fit_function_;
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
  Double_t tail_ratio_max_ = 100.0;
  std::vector<Double_t> manual_params_;

  Double_t EstimateBackground();
  Double_t ClampToBounds(Int_t param_index, Double_t value);

  void SaveInteractiveParams(const TString &input_name,
                             const TString &peak_name);
  Bool_t LoadInteractiveParams(const TString &input_name,
                               const TString &peak_name);

  void SortPeaksByMu(Int_t num_peaks);
  void AppendPeakGraphs(std::vector<TGraph *> &components, Int_t param_offset,
                        Style_t line_style, TF1 *background, Int_t npts,
                        Double_t x_step);

public:
  /**
   * @brief Build a fitter for one histogram and range.
   *
   * @param working_hist        Histogram to fit. Borrowed, not owned; it must
   *                            outlive this object.
   * @param fit_range_low       Lower fit bound, in the histogram's x units.
   * @param fit_range_high      Upper fit bound, in the histogram's x units.
   * @param use_flat_background `kTRUE` for a constant background, `kFALSE`
   *                            (default) for a linear one.
   * @param use_step            Enable the step shelf.
   * @param use_low_exp_tail    Enable the low-energy exponential tail.
   * @param use_low_lin_tail    Enable the low-energy linear tail.
   * @param use_high_exp_tail   Enable the high-energy exponential tail.
   */
  FittingUtils(TH1 *working_hist, Float_t fit_range_low, Float_t fit_range_high,
               Bool_t use_flat_background = kFALSE, Bool_t use_step = kFALSE,
               Bool_t use_low_exp_tail = kFALSE,
               Bool_t use_low_lin_tail = kFALSE,
               Bool_t use_high_exp_tail = kFALSE);
  /// @brief Releases the fit function this object owns.
  ~FittingUtils();

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
   * On the first run the editor opens once the automated fit finishes, and the
   * accepted parameters and range are written to
   * `<plots_base>/fits/<peak_name>_<input_name>.fits`. On later runs that file
   * is loaded and the editor is skipped, so the manual step happens once.
   * Delete the file to redo it.
   *
   * @param interactive `kTRUE` to enable.
   *
   * @note Re-running the macro after an interactive session is worthwhile: the
   *       saved values are much better starting points and often give a lower
   *       chi-squared on the second pass.
   */
  void SetInteractive(Bool_t interactive = kTRUE) {
    interactive_ = interactive;
  }
  /// @brief Cap the tail decay constants during fitting.
  /// @param ratio_max Upper bound on any tail ratio, in units of sigma.
  ///                  Defaults to 100.
  void SetTailRatioMax(Double_t ratio_max) { tail_ratio_max_ = ratio_max; }

  /**
   * @brief Supply explicit starting values for every parameter.
   * @param params Initial values, in the parameter order documented on the
   *               matching FittingFunctions model. Switches the fitter off its
   *               automatic initial-guess path.
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

  /// @brief The fit function, for inspection or manual adjustment.
  /// @return The owned `TF1`, or null before a fit has built one.
  TF1 *GetFitFunction() { return fit_function_; }
  /**
   * @brief Substitute a fit function of your own.
   * @param func Replacement function.
   * @warning Overwrites the pointer without deleting the previous function, so
   *          the old one leaks unless you retrieved and freed it first. The
   *          replacement is destroyed with this object.
   */
  void SetFitFunction(TF1 *func) { fit_function_ = func; }

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
   * @param input_name Input identifier, used for saved parameters and output
   *                   filenames.
   * @param peak_name  Peak identifier, used the same way.
   * @return The fitted result; check FitResult::valid, since a failed fit
   *         returns `-1` in every parameter.
   */
  FitResult FitSinglePeak(const TString input_name, const TString peak_name);
  /**
   * @brief Fit two peaks from centroid guesses.
   * @param input_name Input identifier.
   * @param peak_name  Peak identifier.
   * @param mu1_init   Starting centroid for the first peak.
   * @param mu2_init   Starting centroid for the second.
   * @return The fitted result, with peaks sorted by centroid.
   */
  FitResult FitDoublePeak(const TString input_name, const TString peak_name,
                          Double_t mu1_init, Double_t mu2_init);
  /**
   * @brief Fit two peaks with the first constrained by an earlier fit.
   *
   * Useful when one peak is well measured in isolation and the other is not:
   * the known shape anchors the fit instead of floating freely.
   *
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
};

#endif
