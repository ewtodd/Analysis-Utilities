#ifndef ROOFITPHOTOPEAKKERNELS_HPP
#define ROOFITPHOTOPEAKKERNELS_HPP

// Host-side launchers for the photopeak CUDA kernels in gpu/*.cu. They are
// compiled only for the CUDA build (see CMakeLists.txt), so a CPU build of the
// library never references them. The PDF classes in
// include/RooFitPhotopeakPdfs.hpp call these from their doEval() methods when
// the output buffer is a device pointer.
//
// Each launcher takes the per-node CUDA stream from ROOT's batch-compute
// Config (RooBatchCompute::CudaInterface::CudaStream) and launches its kernel
// on that stream, so simultaneous multi-node fits overlap on the GPU instead of
// serializing on the legacy default stream. CudaStream is only forward-declared
// here: the CudaStream -> cudaStream_t conversion operator is declared
// exclusively under __CUDACC__ in CudaInterface.h, so the conversion is
// performed inside the .cu files (which include that header), not in the host
// .cpp. That header and RooBatchCompute.h are vendored verbatim from ROOT
// 6.40.00's roofit/batchcompute tree (include/roofit_batchcompute/), because
// ROOT's nix package does not install them. A null stream selects the default
// stream.

#include <cstddef>

namespace RooBatchCompute {
namespace CudaInterface {
class CudaStream;
} // namespace CudaInterface
} // namespace RooBatchCompute

/// \brief Launch the low-energy exponential-tail kernel.
///
/// Computes output[i] = exp((x[i]-mu)/tau) * erfc((x[i]-mu)/(sigma*sqrt(2)))
/// with an overflow-safe exp*erfc (clamped to a 1e-300 floor) for all i < n.
///
/// \param output Device buffer of n doubles to fill (owned by the caller).
/// \param x_vals Device buffer of n doubles, the observable values.
/// \param mu Peak position, in the same units as x_vals.
/// \param inv_tau Inverse low-energy exponential tail width (1/tau).
/// \param inv_sqrt2_sigma Inverse (sigma*sqrt(2)).
/// \param n Number of elements to process.
/// \param stream Per-node CUDA stream from the batch-compute Config to launch
///        on, so this kernel overlaps with the other graph nodes' work; nullptr
///        selects the default stream. Not owned; must outlive the launch.
void RooLowExpTail_launchKernel(
    double *output, const double *x_vals, double mu, double inv_tau,
    double inv_sqrt2_sigma, std::size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

/// \brief Launch the low-energy linear-tail kernel.
///
/// Computes output[i] = (1 + slope*(x[i]-mu)) * erfc((x[i]-mu)/(sigma*sqrt(2)))
/// for all i < n, clamping the linear factor at zero and the result to a 1e-300
/// floor.
///
/// \param output Device buffer of n doubles to fill (owned by the caller).
/// \param x_vals Device buffer of n doubles, the observable values.
/// \param mu Peak position, in the same units as x_vals.
/// \param slope Low-energy linear tail slope (dimensionless per unit x).
/// \param inv_sqrt2_sigma Inverse (sigma*sqrt(2)).
/// \param n Number of elements to process.
/// \param stream Per-node CUDA stream from the batch-compute Config to launch
///        on, so this kernel overlaps with the other graph nodes' work; nullptr
///        selects the default stream. Not owned; must outlive the launch.
void RooLowLinTail_launchKernel(
    double *output, const double *x_vals, double mu, double slope,
    double inv_sqrt2_sigma, std::size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

/// \brief Launch the step-function shelf kernel.
///
/// Computes output[i] = sigmoid(-((x[i]-mu)/sigma))^2, i.e. a smoothed step
/// that is ~1 below the shelf and ~0 above it, clamped to a 1e-300 floor, for
/// all i < n.
///
/// \param output Device buffer of n doubles to fill (owned by the caller).
/// \param x_vals Device buffer of n doubles, the observable values.
/// \param mu Shelf position, in the same units as x_vals.
/// \param inv_sigma Inverse shelf transition width (1/sigma).
/// \param n Number of elements to process.
/// \param stream Per-node CUDA stream from the batch-compute Config to launch
///        on, so this kernel overlaps with the other graph nodes' work; nullptr
///        selects the default stream. Not owned; must outlive the launch.
void RooStepShelf_launchKernel(
    double *output, const double *x_vals, double mu, double inv_sigma,
    std::size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

/// \brief Launch the high-energy exponential-tail kernel.
///
/// Computes output[i] = exp(-(x[i]-mu)/tau) * erfc((x[i]-mu)/(sigma*sqrt(2)))
/// with an overflow-safe exp*erfc (clamped to a 1e-300 floor) for all i < n.
///
/// \param output Device buffer of n doubles to fill (owned by the caller).
/// \param x_vals Device buffer of n doubles, the observable values.
/// \param mu Peak position, in the same units as x_vals.
/// \param inv_tau Inverse high-energy exponential tail width (1/tau).
/// \param inv_sqrt2_sigma Inverse (sigma*sqrt(2)).
/// \param n Number of elements to process.
/// \param stream Per-node CUDA stream from the batch-compute Config to launch
///        on, so this kernel overlaps with the other graph nodes' work; nullptr
///        selects the default stream. Not owned; must outlive the launch.
void RooHighExpTail_launchKernel(
    double *output, const double *x_vals, double mu, double inv_tau,
    double inv_sqrt2_sigma, std::size_t n,
    RooBatchCompute::CudaInterface::CudaStream *stream = nullptr);

#endif
