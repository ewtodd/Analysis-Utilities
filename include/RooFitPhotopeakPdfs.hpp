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
 * `RooAddPdf` by RooFitUtils.
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

#endif
