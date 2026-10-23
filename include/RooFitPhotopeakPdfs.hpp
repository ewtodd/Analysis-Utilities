#ifndef ROOFITPHOTOPEAKPDFS_H
#define ROOFITPHOTOPEAKPDFS_H

#include <RooAbsPdf.h>
#include <RooAbsReal.h>
#include <RooFit/EvalContext.h>
#include <RooRealProxy.h>

/**
 * @file RooFitPhotopeakPdfs.hpp
 * @brief Custom RooFit PDF components for gamma-ray photopeak fitting.
 *
 * Four `RooAbsPdf` components modelling the departures from a pure Gaussian
 * seen in semiconductor detectors: a step shelf from energy escaping the active
 * volume, low-energy tails from incomplete charge collection, and a high-energy
 * tail from pileup. They are combined with a Gaussian and a background into a
 * `RooAddPdf` by RooFitUtils. A fifth class, RooExpTailWindowFactor, is not a
 * PDF but the scalar that turns an exponential tail's full-range yield ratio
 * into its yield inside the fit window.
 *
 * Conventions shared by all four:
 *
 * - Every component is written in terms of `y = x - mu`, and every shape
 *   parameter is expressed relative to the resolution `sigma`, so the shapes
 *   stay meaningful as the resolution floats.
 * - The exponential tails take a dimensionless `tau_ratio`; the decay constant
 *   is `tau = tau_ratio * sigma`.
 * - Densities are floored at `1e-300` rather than being allowed to reach zero.
 *   A component that evaluates to exactly zero across the fit range integrates
 *   to zero, and RooFit's extended `Range()` projection then forms a 0/0 ratio
 *   and yields NaN. The floor is negligible wherever the real density is not.
 * - Non-physical parameters (`sigma <= 0`, `tau_ratio <= 0`) yield zero from
 *   the scalar path and the density floor from the batch path, rather than
 *   raising.
 * - Each implements `doEval()` so a likelihood pass is vectorised over events
 *   instead of dispatching a scalar `evaluate()` per event, and each provides
 *   an analytic integral over `x` so normalisation is not numeric.
 * - When built with `-DAU_USE_CUDA=ON`, `canComputeBatchWithCuda()` reports
 *   `true` and `doEval()` dispatches to the kernels in
 *   RooFitPhotopeakKernels.hpp whenever RooFit supplies device-resident
 *   buffers. In a CPU build it reports `false` and the OpenMP host loop runs.
 *
 * @note These classes carry ROOT dictionaries (see `include/LinkDef.h`), which
 *       is what makes them usable from Cling and PyROOT. The default
 *       constructors exist for I/O and leave the proxies unbound; do not
 *       evaluate a default-constructed instance.
 */

/**
 * @brief Resolution-smeared step on the low side of the photopeak.
 *
 * Models events where part of the photon energy escapes the active volume,
 * producing a step-like shelf below the peak, smeared by the detector
 * resolution. Built from a sigmoid in `z = (x - mu) / sigma`, squared, which
 * gives a shelf that falls smoothly to zero above the centroid.
 */
class RooStepShelf : public RooAbsPdf {
public:
  /// @brief Default constructor for ROOT I/O; leaves proxies unbound.
  RooStepShelf() {}
  /**
   * @brief Construct the component and bind its parameters.
   * @param name  RooFit object name, unique within the workspace.
   * @param title Human-readable title.
   * @param x     Observable (energy).
   * @param mu    Peak centroid, in the observable's units.
   * @param sigma Gaussian resolution; must be positive.
   */
  RooStepShelf(const char *name, const char *title, RooAbsReal &x,
               RooAbsReal &mu, RooAbsReal &sigma);
  /// @brief Copy constructor used by RooFit when cloning a model.
  /// @param other Instance to copy.
  /// @param name  New object name, or null to keep @p other's.
  RooStepShelf(const RooStepShelf &other, const char *name = nullptr);
  /// @brief RooFit clone hook.
  /// @param newname Name for the clone, or null to reuse this one's.
  /// @return A heap-allocated copy; ownership passes to the caller.
  TObject *clone(const char *newname = nullptr) const override {
    return new RooStepShelf(*this, newname);
  }

  /**
   * @brief Advertise the analytic integral over the observable.
   * @param allVars  Variables RooFit offers to integrate over.
   * @param analVars Filled with those handled analytically.
   * @param rangeName Named range, or null for the full range.
   * @return `1` if the integral over `x` is available, `0` otherwise.
   */
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  /**
   * @brief Evaluate the analytic integral advertised above.
   * @param code      Code returned by getAnalyticalIntegral(); only `1` is
   *                  supported.
   * @param rangeName Named range, or null for the full range.
   * @return The integral over the requested range.
   */
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  /**
   * @brief Batched evaluation over every event in one pass.
   *
   * Hoists the parameter-dependent constants out of the loop and writes one
   * density per event. Dispatches to a CUDA kernel when built with CUDA support
   * and handed device-resident buffers; otherwise runs an OpenMP host loop.
   *
   * @param ctx RooFit evaluation context supplying inputs and the output span.
   */
  void doEval(RooFit::EvalContext &ctx) const override;
  /// @brief Whether doEval() can consume device-resident buffers.
  /// @return `true` in a CUDA build, `false` in a CPU build.
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  /// @brief Scalar evaluation for a single observable value.
  /// @return The unnormalised density, floored at `1e-300`; zero if the
  ///         parameters are non-physical.
  Double_t evaluate() const override;

  ClassDefOverride(RooStepShelf, 1)
};

/**
 * @brief Exponential tail below the photopeak, convolved with the resolution.
 *
 * Models incomplete charge collection and charge trapping, which move counts
 * from the peak to lower energies. The density is the analytic convolution of a
 * one-sided exponential with a Gaussian, evaluated as `exp(a) * erfc(b)` in a
 * form that avoids intermediate overflow by folding the two exponentials
 * together and switching to the large-argument expansion of the scaled
 * complementary error function.
 */
class RooLowExpTail : public RooAbsPdf {
public:
  /// @brief Default constructor for ROOT I/O; leaves proxies unbound.
  RooLowExpTail() {}
  /**
   * @brief Construct the component and bind its parameters.
   * @param name  RooFit object name, unique within the workspace.
   * @param title Human-readable title.
   * @param x     Observable (energy).
   * @param mu    Peak centroid, in the observable's units.
   * @param sigma Gaussian resolution; must be positive.
   * @param tau_ratio Tail decay constant in units of @p sigma, so
   *                  `tau = tau_ratio * sigma`. Dimensionless and positive.
   */
  RooLowExpTail(const char *name, const char *title, RooAbsReal &x,
                RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &tau_ratio);
  /// @brief Copy constructor used by RooFit when cloning a model.
  /// @param other Instance to copy.
  /// @param name  New object name, or null to keep @p other's.
  RooLowExpTail(const RooLowExpTail &other, const char *name = nullptr);
  /// @brief RooFit clone hook.
  /// @param newname Name for the clone, or null to reuse this one's.
  /// @return A heap-allocated copy; ownership passes to the caller.
  TObject *clone(const char *newname = nullptr) const override {
    return new RooLowExpTail(*this, newname);
  }

  /**
   * @brief Advertise the analytic integral over the observable.
   * @param allVars  Variables RooFit offers to integrate over.
   * @param analVars Filled with those handled analytically.
   * @param rangeName Named range, or null for the full range.
   * @return `1` if the integral over `x` is available, `0` otherwise.
   */
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  /**
   * @brief Evaluate the analytic integral advertised above.
   * @param code      Code returned by getAnalyticalIntegral(); only `1` is
   *                  supported.
   * @param rangeName Named range, or null for the full range.
   * @return The integral over the requested range.
   */
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  /**
   * @brief Batched evaluation over every event in one pass.
   *
   * Hoists the parameter-dependent constants out of the loop and writes one
   * density per event. Dispatches to a CUDA kernel when built with CUDA support
   * and handed device-resident buffers; otherwise runs an OpenMP host loop.
   *
   * @param ctx RooFit evaluation context supplying inputs and the output span.
   */
  void doEval(RooFit::EvalContext &ctx) const override;
  /// @brief Whether doEval() can consume device-resident buffers.
  /// @return `true` in a CUDA build, `false` in a CPU build.
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy tau_ratio_;
  /// @brief Scalar evaluation for a single observable value.
  /// @return The unnormalised density, floored at `1e-300`; zero if the
  ///         parameters are non-physical.
  Double_t evaluate() const override;

  ClassDefOverride(RooLowExpTail, 1)
};

/**
 * @brief Linear tail below the photopeak, convolved with the resolution.
 *
 * Captures asymmetric low-side tailing that a single exponential describes
 * poorly. The linear factor `(1 + slope * (x - mu))` is floored at zero to keep
 * the density non-negative, as RooFit normalisation requires; the `TF1` backend
 * in FittingUtils applies no such constraint, so a fit that hits the floor will
 * not agree between the two backends.
 */
class RooLowLinTail : public RooAbsPdf {
public:
  /// @brief Default constructor for ROOT I/O; leaves proxies unbound.
  RooLowLinTail() {}
  /**
   * @brief Construct the component and bind its parameters.
   * @param name  RooFit object name, unique within the workspace.
   * @param title Human-readable title.
   * @param x     Observable (energy).
   * @param mu    Peak centroid, in the observable's units.
   * @param sigma Gaussian resolution; must be positive.
   * @param slope Slope of the linear factor `(1 + slope * (x - mu))`,
   *              in inverse units of the observable.
   */
  RooLowLinTail(const char *name, const char *title, RooAbsReal &x,
                RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &slope);
  /// @brief Copy constructor used by RooFit when cloning a model.
  /// @param other Instance to copy.
  /// @param name  New object name, or null to keep @p other's.
  RooLowLinTail(const RooLowLinTail &other, const char *name = nullptr);
  /// @brief RooFit clone hook.
  /// @param newname Name for the clone, or null to reuse this one's.
  /// @return A heap-allocated copy; ownership passes to the caller.
  TObject *clone(const char *newname = nullptr) const override {
    return new RooLowLinTail(*this, newname);
  }

  /**
   * @brief Advertise the analytic integral over the observable.
   * @param allVars  Variables RooFit offers to integrate over.
   * @param analVars Filled with those handled analytically.
   * @param rangeName Named range, or null for the full range.
   * @return `1` if the integral over `x` is available, `0` otherwise.
   */
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  /**
   * @brief Evaluate the analytic integral advertised above.
   * @param code      Code returned by getAnalyticalIntegral(); only `1` is
   *                  supported.
   * @param rangeName Named range, or null for the full range.
   * @return The integral over the requested range.
   */
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  /**
   * @brief Batched evaluation over every event in one pass.
   *
   * Hoists the parameter-dependent constants out of the loop and writes one
   * density per event. Dispatches to a CUDA kernel when built with CUDA support
   * and handed device-resident buffers; otherwise runs an OpenMP host loop.
   *
   * @param ctx RooFit evaluation context supplying inputs and the output span.
   */
  void doEval(RooFit::EvalContext &ctx) const override;
  /// @brief Whether doEval() can consume device-resident buffers.
  /// @return `true` in a CUDA build, `false` in a CPU build.
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy slope_;
  /// @brief Scalar evaluation for a single observable value.
  /// @return The unnormalised density, floored at `1e-300`; zero if the
  ///         parameters are non-physical.
  Double_t evaluate() const override;

  ClassDefOverride(RooLowLinTail, 1)
};

/**
 * @brief Exponential tail above the photopeak, convolved with the resolution.
 *
 * Models pileup, which moves counts to higher energies. Identical in form to
 * RooLowExpTail but mirrored about the centroid.
 */
class RooHighExpTail : public RooAbsPdf {
public:
  /// @brief Default constructor for ROOT I/O; leaves proxies unbound.
  RooHighExpTail() {}
  /**
   * @brief Construct the component and bind its parameters.
   * @param name  RooFit object name, unique within the workspace.
   * @param title Human-readable title.
   * @param x     Observable (energy).
   * @param mu    Peak centroid, in the observable's units.
   * @param sigma Gaussian resolution; must be positive.
   * @param tau_ratio Tail decay constant in units of @p sigma, so
   *                  `tau = tau_ratio * sigma`. Dimensionless and positive.
   */
  RooHighExpTail(const char *name, const char *title, RooAbsReal &x,
                 RooAbsReal &mu, RooAbsReal &sigma, RooAbsReal &tau_ratio);
  /// @brief Copy constructor used by RooFit when cloning a model.
  /// @param other Instance to copy.
  /// @param name  New object name, or null to keep @p other's.
  RooHighExpTail(const RooHighExpTail &other, const char *name = nullptr);
  /// @brief RooFit clone hook.
  /// @param newname Name for the clone, or null to reuse this one's.
  /// @return A heap-allocated copy; ownership passes to the caller.
  TObject *clone(const char *newname = nullptr) const override {
    return new RooHighExpTail(*this, newname);
  }

  /**
   * @brief Advertise the analytic integral over the observable.
   * @param allVars  Variables RooFit offers to integrate over.
   * @param analVars Filled with those handled analytically.
   * @param rangeName Named range, or null for the full range.
   * @return `1` if the integral over `x` is available, `0` otherwise.
   */
  Int_t getAnalyticalIntegral(RooArgSet &allVars, RooArgSet &analVars,
                              const char *rangeName = nullptr) const override;
  /**
   * @brief Evaluate the analytic integral advertised above.
   * @param code      Code returned by getAnalyticalIntegral(); only `1` is
   *                  supported.
   * @param rangeName Named range, or null for the full range.
   * @return The integral over the requested range.
   */
  Double_t analyticalIntegral(Int_t code,
                              const char *rangeName = nullptr) const override;

  /**
   * @brief Batched evaluation over every event in one pass.
   *
   * Hoists the parameter-dependent constants out of the loop and writes one
   * density per event. Dispatches to a CUDA kernel when built with CUDA support
   * and handed device-resident buffers; otherwise runs an OpenMP host loop.
   *
   * @param ctx RooFit evaluation context supplying inputs and the output span.
   */
  void doEval(RooFit::EvalContext &ctx) const override;
  /// @brief Whether doEval() can consume device-resident buffers.
  /// @return `true` in a CUDA build, `false` in a CPU build.
#ifdef AU_ROOFIT_BACKEND_CUDA
  inline bool canComputeBatchWithCuda() const override { return true; }
#else
  inline bool canComputeBatchWithCuda() const override { return false; }
#endif

protected:
  RooRealProxy x_;
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy tau_ratio_;
  /// @brief Scalar evaluation for a single observable value.
  /// @return The unnormalised density, floored at `1e-300`; zero if the
  ///         parameters are non-physical.
  Double_t evaluate() const override;

  ClassDefOverride(RooHighExpTail, 1)
};

/**
 * @brief Converts a full-range tail/Gaussian yield ratio into the window.
 *
 * Evaluates `f_T / f_G`, where `f_T` is the fraction of an exponential tail's
 * integral over the whole line that falls inside the window
 * `[window_low, window_high]`, and `f_G` is the same fraction for the Gaussian
 * of the same peak. RooFitUtils normalises every component over the window and
 * parametrises each exponential tail by its full-range yield ratio `R` (tail
 * integral over the whole line divided by the Gaussian's), so the tail's yield
 * inside the window is
 *
 *     GausAmplitude * R * f_T / f_G
 *
 * with `GausAmplitude` the Gaussian's yield inside the window. `R` depends only
 * on the tail's shape, not on the window or on where the peak sits in it, so a
 * `LinkParameter()` on it gives the linked peaks the same physical tail.
 *
 * With `y = x - mu` for the low-side tail, `y = mu - x` for the high-side one
 * (RooHighExpTail is RooLowExpTail mirrored), `tau = tau_ratio * sigma` and
 * `[y_lo, y_hi]` the window in `y`:
 *
 *     f_G = [erf((window_high - mu)/(sqrt2 sigma))
 *            - erf((window_low - mu)/(sqrt2 sigma))] / 2
 *     f_T = [A(y_hi) - A(y_lo)] / (2 tau exp(sigma^2 / (2 tau^2)))
 *
 * where `A` is the antiderivative of the tail shape
 * `exp(y/tau) erfc(y/(sqrt2 sigma))` that the tails' analyticalIntegral() uses,
 * and the denominator of `f_T` is that shape's integral over the whole line.
 *
 * @note A scalar of the peak's parameters, not of the observable. RooFit
 *       evaluates it once per likelihood evaluation on the CPU, in CUDA fits
 *       too, as it does the `RooFormulaVar` yields it feeds.
 * @note The window edges are ordinary `RooAbsReal` servers, so the factor
 *       follows a change of window only when they are changed. RooFitUtils
 *       holds them as constant `RooRealVar`s and moves them with the range its
 *       yields count events in; see RooFitPeakModel::SetWindow().
 */
class RooExpTailWindowFactor : public RooAbsReal {
public:
  /// @brief Default constructor for ROOT I/O; leaves proxies unbound.
  RooExpTailWindowFactor() {}
  /**
   * @brief Construct the factor and bind its parameters.
   * @param name        RooFit object name, unique within the workspace.
   * @param title       Human-readable title.
   * @param mu          Peak centroid, in the observable's units.
   * @param sigma       Gaussian resolution; must be positive.
   * @param tau_ratio   Tail decay constant in units of @p sigma, the same
   *                    variable the tail PDF takes. Dimensionless and positive.
   * @param window_low  Lower edge of the window, in the observable's units.
   * @param window_high Upper edge of the window; above @p window_low.
   * @param high_side   `kTRUE` for a RooHighExpTail, `kFALSE` for a
   *                    RooLowExpTail.
   */
  RooExpTailWindowFactor(const char *name, const char *title, RooAbsReal &mu,
                         RooAbsReal &sigma, RooAbsReal &tau_ratio,
                         RooAbsReal &window_low, RooAbsReal &window_high,
                         Bool_t high_side);
  /// @brief Copy constructor used by RooFit when cloning a model.
  /// @param other Instance to copy.
  /// @param name  New object name, or null to keep @p other's.
  RooExpTailWindowFactor(const RooExpTailWindowFactor &other,
                         const char *name = nullptr);
  /// @brief RooFit clone hook.
  /// @param newname Name for the clone, or null to reuse this one's.
  /// @return A heap-allocated copy; ownership passes to the caller.
  TObject *clone(const char *newname = nullptr) const override {
    return new RooExpTailWindowFactor(*this, newname);
  }

protected:
  RooRealProxy mu_;
  RooRealProxy sigma_;
  RooRealProxy tau_ratio_;
  RooRealProxy window_low_;
  RooRealProxy window_high_;
  Bool_t high_side_ = kFALSE;
  /// @brief Evaluate `f_T / f_G` at the current parameter values.
  /// @return The factor; zero if `sigma <= 0` or `tau_ratio <= 0`. `f_G` is
  ///         floored at `1e-300` so a window that misses the Gaussian gives a
  ///         large finite factor instead of a division by zero.
  Double_t evaluate() const override;

  ClassDefOverride(RooExpTailWindowFactor, 1)
};

#endif
