#ifndef AHF_GPU_H
#define AHF_GPU_H
#include <cuda_runtime.h>
#include <cstddef>
#include <vector>
#include "ahf_native_types.h"
#include "ahf_ordered_csr.h"
#include "ahf_newton_gpu.h"
#include "ahf_diagnostics_gpu.h"

struct GpuAhfLinearSolver;
class Block;
class bssn_class;
class var;
template<class T> class MyList;
namespace AHFinderDirect {
class patch_system;
struct BH_diagnostics;
namespace gpu {

struct SparsePlan {
    int *rows = nullptr, *columns = nullptr;
};
struct DeviceSurface {
    int n = 0, ng = 0, nnz = 0, symmetry = 0;
    // Static angles only; all radius-dependent geometry is evaluated on device.
    ahf_native::StaticMetadata *metadata = nullptr;
    ahf_native::RuntimeGeometry *runtime_geometry = nullptr;
    ahf_diagnostics::DevicePlan diagnostic_plan;
    SparsePlan angular;
    int *instruction_rows = nullptr, *nominal_ghost = nullptr;
    OrderedInstruction *instructions = nullptr;
    AhfNewtonState *newton = nullptr;
    double *initial_direction = nullptr;
    double *h = nullptr, *oldh = nullptr, *ghost_h = nullptr;
    double *xyz = nullptr, *raw = nullptr, *geometry = nullptr, *gauge = nullptr;
    double *theta = nullptr, *coefficients = nullptr, *values = nullptr;
    double *delta = nullptr, *saved_theta = nullptr, *scratch6 = nullptr;
    double *scalars = nullptr, *diagnostics = nullptr, *signature = nullptr;
    double *origin = nullptr;
    int *status = nullptr;
};
struct Surface {
    DeviceSurface d;
    cudaStream_t stream = nullptr;
    GpuAhfLinearSolver *linear_solver = nullptr;
    patch_system *reference = nullptr;
    int exact_pole_signature_nans = 0;
    Surface(patch_system &ps, int symmetry);
    ~Surface();
    Surface(const Surface &) = delete;
    Surface &operator=(const Surface &) = delete;
};
struct TransferStats {
    unsigned long long metadata_h2d = 0, mpi_h2d = 0, mpi_d2h = 0;
    unsigned long long output_d2h = 0, scalar_d2h = 0;
};
TransferStats &transfer_stats();
void check(cudaError_t error, const char *operation);

// Inputs and outputs are logical-origin padded device pointers, not allocation bases.
void derivatives(cudaStream_t stream, const int shape[3], const int padded[3],
                 const double *const coordinates[3], const double *const fields[7],
                 double *const outputs[21], int symmetry, int shell_type = -1,
                 const double *const transform[9] = nullptr);
void initialize(Surface &surface, const double center[3], const double radii[3]);
void save_old_shape(Surface &surface);
void synchronize(Surface &surface);
void coordinates(Surface &surface);
void transform_geometry(cudaStream_t stream, int n, const double *raw, double *physical);
void evaluate(Surface &surface, double add_to_expansion, bool coefficients = true);
void norms(Surface &surface, double result[3]); // Theta inf, mean h, status
void assemble(Surface &surface);
void perturb(Surface &surface);
void restore_and_radial(Surface &surface, bool success);
int linear_step(Surface &surface);
void update(Surface &surface);
void diagnostics(Surface &surface);
void signature(Surface &surface, double dt);
void recenter(Surface &surface);
void download_diagnostics(Surface &surface, BH_diagnostics &result);
void download_output(Surface &surface, bool with_signature);

// Collective, including dummy/nonowner ranks. Refreshes local grid addresses each call.
// Only horizon-sized communication buffers are staged; no host numerical field work.
int sample_collective(Surface *surface, int owner, bssn_class &adm,
                      MyList<var> *fields, int symmetry, bool gauge = false);
void setup();
void cleanup();
void cleanup_sampling();
void prepare_derivatives(bssn_class &adm, MyList<var> *outputs);
void find_horizons(int hn, int *dumpid, double *xc, double *yc, double *zc,
                   double *xr, double *yr, double *zr, bool *trigger, double *dt);
}
}
#endif
