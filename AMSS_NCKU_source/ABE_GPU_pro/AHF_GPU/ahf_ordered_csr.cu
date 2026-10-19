#include "ahf_ordered_csr.h"
#include <stdexcept>
namespace AHFinderDirect { namespace gpu {
namespace {
__global__ void ordered_kernel(DeviceOrderedPlan plan, const double *coefficients,
                               double *values) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= plan.n) return;
    const double *c = coefficients + 5 * row;
    // One thread per row: stable instruction order, disjoint destinations.
    for (int j = plan.instruction_rows[row]; j < plan.instruction_rows[row + 1]; ++j) {
        const OrderedInstruction op = plan.instructions[j];
        double mol = c[op.coefficient0] * op.molecule0;
        if (op.kind != OrderedTermKind::mixed) {
            const double second = c[op.coefficient1] * op.molecule1;
            mol = mol + second; // Group Jac_sum before multiplying ghost factor.
        }
        const double contribution = op.ghost ? mol * op.ghost_factor : mol;
        if (op.first) values[op.slot] = contribution;
        else values[op.slot] += contribution;
    }
}
}
void assemble_ordered(cudaStream_t stream, DeviceOrderedPlan plan,
                      const double *coefficients, double *values) {
    if (plan.n <= 0 || !plan.instruction_rows || !plan.instructions ||
        !coefficients || !values)
        throw std::invalid_argument("AHF invalid device ordered plan");
    ordered_kernel<<<(plan.n + 127) / 128, 128, 0, stream>>>(plan, coefficients, values);
    const cudaError_t result = cudaGetLastError();
    if (result != cudaSuccess) throw std::runtime_error(cudaGetErrorString(result));
}
}}
