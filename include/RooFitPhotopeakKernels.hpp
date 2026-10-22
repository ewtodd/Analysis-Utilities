#ifndef ROOFITPHOTOPEAKKERNELS_HPP
#define ROOFITPHOTOPEAKKERNELS_HPP

#include <cstddef>

/**
 * @file RooFitPhotopeakKernels.hpp
 * @brief Host-callable launch wrappers for the photopeak PDF CUDA kernels.
 *
 * Each launcher takes device pointers plus the constants hoisted once per
 * Minuit evaluation — reciprocals are passed in already inverted so the kernel
 * multiplies rather than divides per element — dispatches the kernel on the
 * CUDA stream RooFit's evaluator assigned to the calling node, and returns as
 * soon as that stream has accepted the launch.
 *
 * These are declared unconditionally, but only defined when the library is
 * built with `-DAU_USE_CUDA=ON`. The CPU build never compiles `gpu/*.cu`, so
 * calling one of these in a CPU build is a link error. The custom PDFs reach
 * them only from a `doEval` path already guarded by
 * `canComputeBatchWithCuda()`.
 *
 * @warning Every pointer argument must be device-resident. The callers obtain
 *          them from RooFit's evaluation context and check with
 *          `cudaPointerGetAttributes` before dispatching.
 *
 * @note No implicit synchronisation: the launch is asynchronous. Ordering
 *       comes from the stream. Before calling `doEval`, RooFit's evaluator
 *       makes the node's stream wait for the node's inputs, and afterwards it
 *       records the node's completion event on that same stream. A kernel
 *       launched on that stream is therefore ordered against the rest of the
 *       computation graph without any synchronisation here.
 */

/// @cond
// Only forward-declared: CudaStream's conversion to cudaStream_t exists only
// under __CUDACC__, so the full definition (ROOT's private CudaInterface.h,
// supplied by flake.nix as AU_ROOT_BATCHCOMPUTE_INCLUDE_DIR) is needed in
// gpu/*.cu but must not leak into host code or into this installed header.
namespace RooBatchCompute {
namespace CudaInterface {
class CudaStream;
} // namespace CudaInterface
} // namespace RooBatchCompute
/// @endcond

/**
 * @brief Launch the low-energy exponential tail kernel.
 *
 * Computes `output[i] = exp((x[i] - mu) / tau) * erfc((x[i] - mu) /
 * (sqrt(2) * sigma))` for all `i < n`, using an overflow-safe `exp * erfc`
 * and flooring the result at `1e-300`.
 *
 * @param output          Device buffer for @p n PDF values.
 * @param x_vals          Device buffer of @p n observable values.
 * @param mu              Peak centroid, in the observable's units.
 * @param inv_tau         Reciprocal of the tail decay constant.
 * @param inv_sqrt2_sigma `1 / (sqrt(2) * sigma)`, the resolution factor.
 * @param n               Number of elements.
 * @param stream          Stream RooFit's evaluator assigned to the calling
 *                        node, read in `doEval` as
 *                        `ctx.config(this).cudaStream()`. Owned by the
 *                        evaluator. Null launches on the legacy default
 *                        stream, which is still correctly ordered but
 *                        serialises against every other node.
 */
void RooLowExpTail_launchKernel(
    double *output, const double *x_vals, double mu, double inv_tau,
    double inv_sqrt2_sigma, size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

/**
 * @brief Launch the resolution-smeared step-shelf kernel.
 *
 * Computes `output[i] = sigmoid(-(x[i] - mu) / sigma)^2`, a smoothed step that
 * is ~1 below @p mu and ~0 above it, for all `i < n`, flooring the result at
 * `1e-300`.
 *
 * @param output    Device buffer for @p n PDF values.
 * @param x_vals    Device buffer of @p n observable values.
 * @param mu        Peak centroid, in the observable's units.
 * @param inv_sigma Reciprocal of the Gaussian resolution.
 * @param n         Number of elements.
 * @param stream    Stream RooFit's evaluator assigned to the calling node, read
 *                  in `doEval` as `ctx.config(this).cudaStream()`. Owned by the
 *                  evaluator. Null launches on the legacy default stream, which
 *                  is still correctly ordered but serialises against every
 *                  other node.
 */
void RooStepShelf_launchKernel(
    double *output, const double *x_vals, double mu, double inv_sigma, size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

/**
 * @brief Launch the low-energy linear tail kernel.
 *
 * Computes `output[i] = (1 + slope * (x[i] - mu)) * erfc((x[i] - mu) /
 * (sqrt(2) * sigma))` for all `i < n`, flooring the linear factor at zero and
 * the result at `1e-300`.
 *
 * @param output          Device buffer for @p n PDF values.
 * @param x_vals          Device buffer of @p n observable values.
 * @param mu              Peak centroid, in the observable's units.
 * @param slope           Linear tail slope `s` in the `(1 + s*(x - mu))`
 *                        factor, which is floored at zero to keep the PDF
 *                        non-negative.
 * @param inv_sqrt2_sigma `1 / (sqrt(2) * sigma)`, the resolution factor.
 * @param n               Number of elements.
 * @param stream          Stream RooFit's evaluator assigned to the calling
 *                        node, read in `doEval` as
 *                        `ctx.config(this).cudaStream()`. Owned by the
 *                        evaluator. Null launches on the legacy default
 *                        stream, which is still correctly ordered but
 *                        serialises against every other node.
 */
void RooLowLinTail_launchKernel(
    double *output, const double *x_vals, double mu, double slope,
    double inv_sqrt2_sigma, size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

/**
 * @brief Launch the high-energy exponential tail kernel.
 *
 * Computes `output[i] = exp(-(x[i] - mu) / tau) * erfc(-(x[i] - mu) /
 * (sqrt(2) * sigma))` for all `i < n`, using an overflow-safe `exp * erfc`
 * and flooring the result at `1e-300`.
 *
 * @param output          Device buffer for @p n PDF values.
 * @param x_vals          Device buffer of @p n observable values.
 * @param mu              Peak centroid, in the observable's units.
 * @param inv_tau         Reciprocal of the tail decay constant.
 * @param inv_sqrt2_sigma `1 / (sqrt(2) * sigma)`, the resolution factor.
 * @param n               Number of elements.
 * @param stream          Stream RooFit's evaluator assigned to the calling
 *                        node, read in `doEval` as
 *                        `ctx.config(this).cudaStream()`. Owned by the
 *                        evaluator. Null launches on the legacy default
 *                        stream, which is still correctly ordered but
 *                        serialises against every other node.
 */
void RooHighExpTail_launchKernel(
    double *output, const double *x_vals, double mu, double inv_tau,
    double inv_sqrt2_sigma, size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

#endif
