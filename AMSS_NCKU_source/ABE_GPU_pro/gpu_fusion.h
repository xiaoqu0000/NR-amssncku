#ifndef GPU_FUSION_H
#define GPU_FUSION_H

#define MAX_FUSED_VARS 32

struct FusedVars {
    int num_vars;
    double* ptrs[MAX_FUSED_VARS];
    double SoA0[MAX_FUSED_VARS];
    double SoA1[MAX_FUSED_VARS];
    double SoA2[MAX_FUSED_VARS];
};

// 用于 t+dt, t -> t+a*dt 的三变量插值
struct FusedTimeVars3 {
    int num_vars;
    double* ptrs1[MAX_FUSED_VARS];
    double* ptrs2[MAX_FUSED_VARS];
    double* ptrs3[MAX_FUSED_VARS];
};

// 用于 t+dt, t, t-dt -> t+a*dt 的四变量插值
struct FusedTimeVars4 {
    int num_vars;
    double* ptrs1[MAX_FUSED_VARS];
    double* ptrs2[MAX_FUSED_VARS];
    double* ptrs3[MAX_FUSED_VARS];
    double* ptrs4[MAX_FUSED_VARS];
};

struct FusedRK4Vars {
    int num_vars;
    double* f0[MAX_FUSED_VARS];
    double* f1[MAX_FUSED_VARS];
    double* f_rhs[MAX_FUSED_VARS];
};

struct FusedSommerfeldVars {
    int num_vars;
    double* f0[MAX_FUSED_VARS];
    double* f_out[MAX_FUSED_VARS]; // 用于 f 或 f_rhs
    double propspeed[MAX_FUSED_VARS];
    double SoA0[MAX_FUSED_VARS];
    double SoA1[MAX_FUSED_VARS];
    double SoA2[MAX_FUSED_VARS];
};

#endif /* GPU_FUSION_H */