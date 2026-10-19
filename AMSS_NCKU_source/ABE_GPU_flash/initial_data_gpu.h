#ifndef INITIAL_DATA_GPU_H
#define INITIAL_DATA_GPU_H

#include <cuda_runtime.h>

enum GpuInitialDataMethod {
    AMSS_GPU_INIT_CAO = 1,
    AMSS_GPU_INIT_LOUSTO = 2
};

void gpu_initial_data_puncture_launch(
    cudaStream_t stream,
    int method,
    const int shape[3],
    const double* x,
    const double* y,
    const double* z,
    double* chi,
    double* trK,
    double* gxx,
    double* gxy,
    double* gxz,
    double* gyy,
    double* gyz,
    double* gzz,
    double* Axx,
    double* Axy,
    double* Axz,
    double* Ayy,
    double* Ayz,
    double* Azz,
    double* Gmx,
    double* Gmy,
    double* Gmz,
    double* Lap,
    double* Sfx,
    double* Sfy,
    double* Sfz,
    double* dtSfx,
    double* dtSfy,
    double* dtSfz,
    const double* mass,
    const double* porg,
    const double* pmom,
    const double* spin,
    int bh_count);

#endif /* INITIAL_DATA_GPU_H */
