#ifndef ROOFITPHOTOPEAKKERNELS_H
#define ROOFITPHOTOPEAKKERNELS_H

#include <cstddef>

/**
 * @file RooFitPhotopeakKernels.hpp
 * @brief Host-callable launch wrappers for the photopeak PDF CUDA kernels.
 *
 * Each launcher takes device pointers plus the constants hoisted once per
 * Minuit evaluation — reciprocals are passed in already inverted so the kernel
 * multiplies rather than divides per element — dispatches the kernel, and
 * returns as soon as the stream has accepted the launch.
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
 * @note No implicit synchronisation: the launch is asynchronous and the caller
 *       is responsible for any synchronisation the surrounding evaluation
 *       needs.
 */

/**
 * @brief Launch the low-energy exponential tail kernel.
 * @param output          Device buffer for @p n PDF values.
 * @param x_vals          Device buffer of @p n observable values.
 * @param mu              Peak centroid, in the observable's units.
 * @param inv_tau         Reciprocal of the tail decay constant.
 * @param inv_sqrt2_sigma `1 / (sqrt(2) * sigma)`, the resolution factor.
 * @param n               Number of elements.
 */
void RooLowExpTail_launchKernel(double *output, const double *x_vals, double mu,
                                double inv_tau, double inv_sqrt2_sigma,
                                size_t n);

/**
 * @brief Launch the resolution-smeared step-shelf kernel.
 * @param output    Device buffer for @p n PDF values.
 * @param x_vals    Device buffer of @p n observable values.
 * @param mu        Peak centroid, in the observable's units.
 * @param inv_sigma Reciprocal of the Gaussian resolution.
 * @param n         Number of elements.
 */
void RooStepShelf_launchKernel(double *output, const double *x_vals, double mu,
                               double inv_sigma, size_t n);

/**
 * @brief Launch the low-energy linear tail kernel.
 * @param output          Device buffer for @p n PDF values.
 * @param x_vals          Device buffer of @p n observable values.
 * @param mu              Peak centroid, in the observable's units.
 * @param slope           Linear tail slope `s` in the `(1 + s*(x - mu))`
 *                        factor, which is floored at zero to keep the PDF
 *                        non-negative.
 * @param inv_sqrt2_sigma `1 / (sqrt(2) * sigma)`, the resolution factor.
 * @param n               Number of elements.
 */
void RooLowLinTail_launchKernel(double *output, const double *x_vals, double mu,
                                double slope, double inv_sqrt2_sigma, size_t n);

/**
 * @brief Launch the high-energy exponential tail kernel.
 * @param output          Device buffer for @p n PDF values.
 * @param x_vals          Device buffer of @p n observable values.
 * @param mu              Peak centroid, in the observable's units.
 * @param inv_tau         Reciprocal of the tail decay constant.
 * @param inv_sqrt2_sigma `1 / (sqrt(2) * sigma)`, the resolution factor.
 * @param n               Number of elements.
 */
void RooHighExpTail_launchKernel(double *output, const double *x_vals,
                                 double mu, double inv_tau,
                                 double inv_sqrt2_sigma, size_t n);

#endif
