#include "RooFitPhotopeakKernels.hpp"

#include <cuda_runtime.h>

#include "roofit_batchcompute/CudaInterface.h"

namespace {

__global__ void RooLowLinTailKernel(double *output, const double *x_vals,
                                    double mu, double slope,
                                    double inv_sqrt2_sigma, size_t n) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double y = x_vals[i] - mu;
    double lin = 1.0 + slope * y;
    if (lin < 0.0) {
      lin = 0.0;
    }
    double v = lin * erfc(y * inv_sqrt2_sigma);
    output[i] = v > 1e-300 ? v : 1e-300;
  }
}

} // namespace

void RooLowLinTail_launchKernel(double *output, const double *x_vals, double mu,
                                double slope, double inv_sqrt2_sigma,
                                size_t n,
                                RooBatchCompute::CudaInterface::CudaStream *stream) {
  const int threads = 256;
  const int blocks = (n + threads - 1) / threads;
  // Launch on the per-node stream supplied by ROOT's batch-compute Config so
  // this kernel overlaps with the other nodes' work instead of serializing on
  // the legacy default stream. nullptr (CPU-built or stream-less context)
  // falls back to the default stream. The CudaStream -> cudaStream_t
  // conversion only exists under __CUDACC__, hence it lives here in the .cu.
  cudaStream_t s = stream ? static_cast<cudaStream_t>(*stream) : cudaStream_t(0);
  RooLowLinTailKernel<<<blocks, threads, 0, s>>>(output, x_vals, mu, slope,
                                                 inv_sqrt2_sigma, n);
}
