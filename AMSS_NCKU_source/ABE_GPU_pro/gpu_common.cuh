#ifndef GPU_COMMON_CUH
#define GPU_COMMON_CUH

// Decompose a flat index into (i, j, k) given dimensions (ex0, ex1).
// Uses mul-sub instead of modulo to avoid expensive integer division where possible.
__device__ __forceinline__
void linear_to_ijk(int tid, int ex0, int ex1, int &i, int &j, int &k) {
    int jk = tid / ex0;
    i = tid - jk * ex0;   // tid % ex0
    k = jk / ex1;
    j = jk - k * ex1;     // jk % ex1
}

#endif /* GPU_COMMON_CUH */
