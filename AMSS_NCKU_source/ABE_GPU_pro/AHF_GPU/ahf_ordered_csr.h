#ifndef AHF_ORDERED_CSR_H
#define AHF_ORDERED_CSR_H
#include <vector>
#include <cuda_runtime.h>
namespace AHFinderDirect {
class patch_system;
namespace gpu {
enum class OrderedTermKind : int { rho, sigma, mixed };
// One original sum_into_element call, including structural zero contributions.
// Nominal terms bypass ghost multiplication; first contributions assign directly.
struct OrderedInstruction {
    int slot;
    OrderedTermKind kind;
    int coefficient0, coefficient1;
    double molecule0, molecule1, ghost_factor;
    int ghost, first;
};
struct OrderedPlan {
    std::vector<int> rows{0}, columns, instruction_rows{0};
    std::vector<OrderedInstruction> instructions;
};
struct DeviceOrderedPlan {
    int n = 0;
    const int *instruction_rows = nullptr;
    const OrderedInstruction *instructions = nullptr;
};
// Requires ps.synchronize() to have initialized interpolation tables. Call only
// at topology setup, never per Newton iteration. No runtime field arithmetic.
OrderedPlan make_ordered_plan(patch_system &ps);
// No host coefficients, CPU injection, atomics, or per-iteration allocation.
void assemble_ordered(cudaStream_t stream, DeviceOrderedPlan plan,
                      const double *coefficients, double *values);
}
}
#endif
