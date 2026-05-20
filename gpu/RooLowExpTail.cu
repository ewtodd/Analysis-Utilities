#include "RooFitPhotopeakKernels.hpp"

#include <cuda_runtime.h>

namespace {

__global__ void RooLowExpTailKernel(double *output, const double *x_vals,
                                    double mu, double inv_tau,
                                    double inv_sqrt2_sigma, size_t n) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double y = x_vals[i] - mu;
    output[i] = exp(y * inv_tau) * erfc(y * inv_sqrt2_sigma);
  }
}

} // namespace

void RooLowExpTail_launchKernel(double *output, const double *x_vals, double mu,
                                double inv_tau, double inv_sqrt2_sigma,
                                size_t n) {
  const int threads = 256;
  const int blocks = (n + threads - 1) / threads;
  RooLowExpTailKernel<<<blocks, threads>>>(output, x_vals, mu, inv_tau,
                                           inv_sqrt2_sigma, n);
}
