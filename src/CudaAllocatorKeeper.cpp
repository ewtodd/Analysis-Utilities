#include "RooFitUtils.hpp"

#ifdef AU_ROOFIT_BACKEND_CUDA
#include <cuda_runtime.h>

namespace {

void *AllocateChunkKeeper() {
  void *keeper = nullptr;
  if (cudaMalloc(&keeper, 256) != cudaSuccess) {
    // Reset the runtime's last-error state so this failure is not reported
    // later by an unrelated CUDA call.
    cudaGetLastError();
    return nullptr;
  }
  return keeper;
}

} // namespace
#endif

void KeepCudaAllocatorChunkMapped() {
#ifdef AU_ROOFIT_BACKEND_CUDA
  // A function-local static is initialised exactly once, even with concurrent
  // first calls. The allocation is never freed on purpose; see the header.
  static void *const keeper = AllocateChunkKeeper();
  (void)keeper;
#endif
}
