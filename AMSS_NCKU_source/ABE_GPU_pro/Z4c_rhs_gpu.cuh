#pragma once

#include <math.h>

__device__ __forceinline__ int z4c_idx3d(int i, int j, int k, int pex0, int pex1)
{
    return i + pex0 * (j + pex1 * k);
}
