#include "RooFitPhotopeakKernels.hpp"

#include <cuda_runtime.h>

namespace {

__global__ void RooHighExpTailKernel(double *output, const double *x_vals,
                                     double mu, double inv_tau,
                                     double inv_sqrt2_sigma, size_t n) {
  size_t i = blockIdx.x * blockDim.x + threadIdx.x;
  if (i < n) {
    double z = mu - x_vals[i];
    output[i] = exp(z * inv_tau) * erfc(z * inv_sqrt2_sigma);
  }
}

} // namespace

void RooHighExpTail_launchKernel(double *output, const double *x_vals,
                                 double mu, double inv_tau,
                                 double inv_sqrt2_sigma, size_t n) {
  const int threads = 256;
  const int blocks = (n + threads - 1) / threads;
  RooHighExpTailKernel<<<blocks, threads>>>(output, x_vals, mu, inv_tau,
                                            inv_sqrt2_sigma, n);
}
