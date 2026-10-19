#ifndef LOPISIDEDIFF_GPU_CUH
#define LOPISIDEDIFF_GPU_CUH

#include "fmisc.h"
#include "fmisc_gpu.cuh"

#include "macrodef.fh"
#include <cmath>

__device__ __forceinline__ int flatten_grid_index(int i, int j, int k, const int ex[3]) {
    return i + ex[0] * (j + ex[1] * k);
}

__device__ __forceinline__ int reflected_lower_bound(int symmetry, double axis_origin, double axis_step, int threshold_symmetry) {
#if ghost_width == 2
    return (symmetry > threshold_symmetry && fabs(axis_origin) < axis_step) ? -2 : 0;
#elif ghost_width == 3
    return (symmetry > threshold_symmetry && fabs(axis_origin) < axis_step) ? -3 : 0;
#elif ghost_width == 4
    return (symmetry > threshold_symmetry && fabs(axis_origin) < axis_step) ? -4 : 0;
#elif ghost_width == 5
    return (symmetry > threshold_symmetry && fabs(axis_origin) < axis_step) ? -5 : 0;
#else
#error Unsupported ghost_width in reflected_lower_bound; expected ghost_width == 2, 3, 4, or 5.
#endif
}

__device__ __forceinline__ double read_lopsided_sample(
    const double* field,
    const int ex[3],
    const double soa[3],
    int i,
    int j,
    int k
) {
#if ghost_width == 2
    return d_symmetry_bd_0b(2, ex[0], ex[1], ex[2], field, i, j, k, soa[0], soa[1], soa[2]);
#elif ghost_width == 3
    return d_symmetry_bd_0b(3, ex[0], ex[1], ex[2], field, i, j, k, soa[0], soa[1], soa[2]);
#elif ghost_width == 4
    return d_symmetry_bd_0b(4, ex[0], ex[1], ex[2], field, i, j, k, soa[0], soa[1], soa[2]);
#elif ghost_width == 5
    return d_symmetry_bd_0b(5, ex[0], ex[1], ex[2], field, i, j, k, soa[0], soa[1], soa[2]);
#else
#error Unsupported ghost_width in read_lopsided_sample; expected ghost_width == 2, 3, 4, or 5.
#endif
}

#if ghost_width == 2 || ghost_width == 3
__device__ __forceinline__ double upwind_bias_contribution(
    double velocity,
    double inverse_step,
    int center,
    int lower_bound,
    int upper_bound,
    double far_backward,
    double backward,
    double center_value,
    double forward,
    double far_forward,
    double extra_backward,
    double extra_forward
) {
    const double zero_value = 0.0;

#if ghost_width == 2
    if (velocity > zero_value) {
        if (center + 2 <= upper_bound) {
            return velocity * inverse_step * (-3.0 * center_value + 4.0 * forward - far_forward);
        }
        if (center + 1 <= upper_bound) {
            return velocity * inverse_step * (-center_value + forward);
        }
        return zero_value;
    }

    if (velocity < zero_value) {
        if (center - 2 >= lower_bound) {
            return -velocity * inverse_step * (-3.0 * center_value + 4.0 * backward - far_backward);
        }
        if (center - 1 >= lower_bound) {
            return -velocity * inverse_step * (-center_value + backward);
        }
    }

    (void)extra_backward;
    (void)extra_forward;
    return zero_value;
#elif ghost_width == 3
    const double factor_three = 3.0;
    const double factor_six = 6.0;
    const double factor_eight = 8.0;
    const double factor_ten = 10.0;
    const double factor_eighteen = 18.0;

    if (velocity > zero_value) {
        if (center + 3 <= upper_bound) {
            return velocity * inverse_step * (
                -factor_three * backward - factor_ten * center_value + factor_eighteen * forward
                -factor_six * far_forward + extra_forward
            );
        }
        if (center + 2 <= upper_bound) {
            return velocity * inverse_step * (
                far_backward - factor_eight * backward + factor_eight * forward - far_forward
            );
        }
        if (center + 1 <= upper_bound) {
            return -velocity * inverse_step * (
                -factor_three * forward - factor_ten * center_value + factor_eighteen * backward
                -factor_six * far_backward + extra_backward
            );
        }
        return zero_value;
    }

    if (velocity < zero_value) {
        if (center - 3 >= lower_bound) {
            return -velocity * inverse_step * (
                -factor_three * forward - factor_ten * center_value + factor_eighteen * backward
                -factor_six * far_backward + extra_backward
            );
        }
        if (center - 2 >= lower_bound) {
            return velocity * inverse_step * (
                far_backward - factor_eight * backward + factor_eight * forward - far_forward
            );
        }
        if (center - 1 >= lower_bound) {
            return velocity * inverse_step * (
                -factor_three * backward - factor_ten * center_value + factor_eighteen * forward
                -factor_six * far_forward + extra_forward
            );
        }
    }

    return zero_value;
#else
#error Unsupported ghost_width in upwind_bias_contribution; expected ghost_width == 2 or 3.
#endif
}
#endif

#if ghost_width == 5
__device__ __forceinline__ double upwind_bias_contribution_fd8(
    double velocity,
    double d840,
    double d60,
    double d12,
    double d2,
    int center,
    int lower_bound,
    int upper_bound,
    double m5,
    double m4,
    double m3,
    double m2,
    double m1,
    double center_value,
    double p1,
    double p2,
    double p3,
    double p4,
    double p5
) {
    const double zero_value = 0.0;

    if (velocity > zero_value) {
        if (center + 5 <= upper_bound && center - 3 >= lower_bound) {
            return velocity * d840 * (-5.0*m3 + 60.0*m2 - 420.0*m1 - 378.0*center_value
                                      + 1050.0*p1 - 420.0*p2 + 140.0*p3 - 30.0*p4 + 3.0*p5);
        }
    } else if (velocity < zero_value) {
        if (center - 5 >= lower_bound && center + 3 <= upper_bound) {
            return -velocity * d840 * (-5.0*p3 + 60.0*p2 - 420.0*p1 - 378.0*center_value
                                       + 1050.0*m1 - 420.0*m2 + 140.0*m3 - 30.0*m4 + 3.0*m5);
        }
    }

    if (center + 4 <= upper_bound && center - 4 >= lower_bound) {
        return velocity * d840 * (3.0*m4 - 32.0*m3 + 168.0*m2 - 672.0*m1
                                  + 672.0*p1 - 168.0*p2 + 32.0*p3 - 3.0*p4);
    }
    if (center + 3 <= upper_bound && center - 3 >= lower_bound) {
        return velocity * d60 * (-m3 + 9.0*m2 - 45.0*m1 + 45.0*p1 - 9.0*p2 + p3);
    }
    if (center + 2 <= upper_bound && center - 2 >= lower_bound) {
        return velocity * d12 * (m2 - 8.0*m1 + 8.0*p1 - p2);
    }
    if (center + 1 <= upper_bound && center - 1 >= lower_bound) {
        return velocity * d2 * (-m1 + p1);
    }
    (void)m5;
    return zero_value;
}
#endif

#if ghost_width == 4
__device__ __forceinline__ double upwind_bias_contribution_fd6(
    double velocity,
    double d60,
    double d12,
    double d2,
    int center,
    int lower_bound,
    int upper_bound,
    double m5,
    double m4,
    double m3,
    double m2,
    double m1,
    double center_value,
    double p1,
    double p2,
    double p3,
    double p4,
    double p5
) {
    const double zero_value = 0.0;

    if (velocity > zero_value) {
        if (center + 4 <= upper_bound && center - 2 >= lower_bound) {
            return velocity * d60 * (2.0*m2 - 24.0*m1 - 35.0*center_value + 80.0*p1 - 30.0*p2 + 8.0*p3 - p4);
        }
        if (center + 5 <= upper_bound && center - 1 >= lower_bound) {
            return velocity * d60 * (-10.0*m1 - 77.0*center_value + 150.0*p1 - 100.0*p2 + 50.0*p3 - 15.0*p4 + 2.0*p5);
        }
    } else if (velocity < zero_value) {
        if (center - 4 >= lower_bound && center + 2 <= upper_bound) {
            return -velocity * d60 * (2.0*p2 - 24.0*p1 - 35.0*center_value + 80.0*m1 - 30.0*m2 + 8.0*m3 - m4);
        }
        if (center - 5 >= lower_bound && center + 1 <= upper_bound) {
            return -velocity * d60 * (-10.0*p1 - 77.0*center_value + 150.0*m1 - 100.0*m2 + 50.0*m3 - 15.0*m4 + 2.0*m5);
        }
    }

    if (center + 3 <= upper_bound && center - 3 >= lower_bound) {
        return velocity * d60 * (-m3 + 9.0*m2 - 45.0*m1 + 45.0*p1 - 9.0*p2 + p3);
    }
    if (center + 2 <= upper_bound && center - 2 >= lower_bound) {
        return velocity * d12 * (m2 - 8.0*m1 + 8.0*p1 - p2);
    }
    if (center + 1 <= upper_bound && center - 1 >= lower_bound) {
        return velocity * d2 * (-m1 + p1);
    }
    return zero_value;
}
#endif

__device__ __forceinline__ double d_lopsided_point(
    const int ex[3], const double* f,
    const double* f_rhs, const double* Sfx, const double* Sfy, const double* Sfz,
    const double* X, const double* Y, const double* Z,
    int symmetry, double SYM1, double SYM2, double SYM3,
    int i, int j, int k // i, j, k is 0-based
) {
    const double ONE = 1.0;
    const int NO_SYMM = 0, EQ_SYMM = 1;

    const double dX = X[1] - X[0];
    const double dY = Y[1] - Y[0];
    const double dZ = Z[1] - Z[0];

#if ghost_width == 2
    const double upwind_dx = ONE / 2.0 / dX;
    const double upwind_dy = ONE / 2.0 / dY;
    const double upwind_dz = ONE / 2.0 / dZ;
#elif ghost_width == 3
    const double upwind_dx = ONE / 12.0 / dX;
    const double upwind_dy = ONE / 12.0 / dY;
    const double upwind_dz = ONE / 12.0 / dZ;
#elif ghost_width == 4
    const double d60dx = ONE / 60.0 / dX;
    const double d60dy = ONE / 60.0 / dY;
    const double d60dz = ONE / 60.0 / dZ;
    const double d12dx = ONE / 12.0 / dX;
    const double d12dy = ONE / 12.0 / dY;
    const double d12dz = ONE / 12.0 / dZ;
    const double d2dx = ONE / 2.0 / dX;
    const double d2dy = ONE / 2.0 / dY;
    const double d2dz = ONE / 2.0 / dZ;
#elif ghost_width == 5
    const double d840dx = ONE / 840.0 / dX;
    const double d840dy = ONE / 840.0 / dY;
    const double d840dz = ONE / 840.0 / dZ;
    const double d60dx = ONE / 60.0 / dX;
    const double d60dy = ONE / 60.0 / dY;
    const double d60dz = ONE / 60.0 / dZ;
    const double d12dx = ONE / 12.0 / dX;
    const double d12dy = ONE / 12.0 / dY;
    const double d12dz = ONE / 12.0 / dZ;
    const double d2dx = ONE / 2.0 / dX;
    const double d2dy = ONE / 2.0 / dY;
    const double d2dz = ONE / 2.0 / dZ;
#else
#error Unsupported ghost_width in d_lopsided_point; expected ghost_width == 2, 3, 4, or 5.
#endif

    const int imax = ex[0] - 1;
    const int jmax = ex[1] - 1;
    const int kmax = ex[2] - 1;

    if (i >= imax || j >= jmax || k >= kmax) return 0.0;

    const int imin = reflected_lower_bound(symmetry, X[0], dX, EQ_SYMM);
    const int jmin = reflected_lower_bound(symmetry, Y[0], dY, EQ_SYMM);
    const int kmin = reflected_lower_bound(symmetry, Z[0], dZ, NO_SYMM);

    const double soa[3] = {SYM1, SYM2, SYM3};

    const int idx = flatten_grid_index(i, j, k, ex);

    const double vx = Sfx[idx];
    const double vy = Sfy[idx];
    const double vz = Sfz[idx];

    const double center_value = read_lopsided_sample(f, ex, soa, i, j, k);

#if ghost_width == 2 || ghost_width == 3
    const double rhs_add_x = upwind_bias_contribution(
        vx,
        upwind_dx,
        i,
        imin,
        imax,
        read_lopsided_sample(f, ex, soa, i - 2, j, k),
        read_lopsided_sample(f, ex, soa, i - 1, j, k),
        center_value,
        read_lopsided_sample(f, ex, soa, i + 1, j, k),
        read_lopsided_sample(f, ex, soa, i + 2, j, k),
        read_lopsided_sample(f, ex, soa, i - 3, j, k),
        read_lopsided_sample(f, ex, soa, i + 3, j, k)
    );

    const double rhs_add_y = upwind_bias_contribution(
        vy,
        upwind_dy,
        j,
        jmin,
        jmax,
        read_lopsided_sample(f, ex, soa, i, j - 2, k),
        read_lopsided_sample(f, ex, soa, i, j - 1, k),
        center_value,
        read_lopsided_sample(f, ex, soa, i, j + 1, k),
        read_lopsided_sample(f, ex, soa, i, j + 2, k),
        read_lopsided_sample(f, ex, soa, i, j - 3, k),
        read_lopsided_sample(f, ex, soa, i, j + 3, k)
    );

    const double rhs_add_z = upwind_bias_contribution(
        vz,
        upwind_dz,
        k,
        kmin,
        kmax,
        read_lopsided_sample(f, ex, soa, i, j, k - 2),
        read_lopsided_sample(f, ex, soa, i, j, k - 1),
        center_value,
        read_lopsided_sample(f, ex, soa, i, j, k + 1),
        read_lopsided_sample(f, ex, soa, i, j, k + 2),
        read_lopsided_sample(f, ex, soa, i, j, k - 3),
        read_lopsided_sample(f, ex, soa, i, j, k + 3)
    );
#elif ghost_width == 4
    const double rhs_add_x = upwind_bias_contribution_fd6(
        vx, d60dx, d12dx, d2dx, i, imin, imax,
        read_lopsided_sample(f, ex, soa, i - 5, j, k),
        read_lopsided_sample(f, ex, soa, i - 4, j, k),
        read_lopsided_sample(f, ex, soa, i - 3, j, k),
        read_lopsided_sample(f, ex, soa, i - 2, j, k),
        read_lopsided_sample(f, ex, soa, i - 1, j, k),
        center_value,
        read_lopsided_sample(f, ex, soa, i + 1, j, k),
        read_lopsided_sample(f, ex, soa, i + 2, j, k),
        read_lopsided_sample(f, ex, soa, i + 3, j, k),
        read_lopsided_sample(f, ex, soa, i + 4, j, k),
        read_lopsided_sample(f, ex, soa, i + 5, j, k)
    );

    const double rhs_add_y = upwind_bias_contribution_fd6(
        vy, d60dy, d12dy, d2dy, j, jmin, jmax,
        read_lopsided_sample(f, ex, soa, i, j - 5, k),
        read_lopsided_sample(f, ex, soa, i, j - 4, k),
        read_lopsided_sample(f, ex, soa, i, j - 3, k),
        read_lopsided_sample(f, ex, soa, i, j - 2, k),
        read_lopsided_sample(f, ex, soa, i, j - 1, k),
        center_value,
        read_lopsided_sample(f, ex, soa, i, j + 1, k),
        read_lopsided_sample(f, ex, soa, i, j + 2, k),
        read_lopsided_sample(f, ex, soa, i, j + 3, k),
        read_lopsided_sample(f, ex, soa, i, j + 4, k),
        read_lopsided_sample(f, ex, soa, i, j + 5, k)
    );

    const double rhs_add_z = upwind_bias_contribution_fd6(
        vz, d60dz, d12dz, d2dz, k, kmin, kmax,
        read_lopsided_sample(f, ex, soa, i, j, k - 5),
        read_lopsided_sample(f, ex, soa, i, j, k - 4),
        read_lopsided_sample(f, ex, soa, i, j, k - 3),
        read_lopsided_sample(f, ex, soa, i, j, k - 2),
        read_lopsided_sample(f, ex, soa, i, j, k - 1),
        center_value,
        read_lopsided_sample(f, ex, soa, i, j, k + 1),
        read_lopsided_sample(f, ex, soa, i, j, k + 2),
        read_lopsided_sample(f, ex, soa, i, j, k + 3),
        read_lopsided_sample(f, ex, soa, i, j, k + 4),
        read_lopsided_sample(f, ex, soa, i, j, k + 5)
    );
#elif ghost_width == 5
    const double rhs_add_x = upwind_bias_contribution_fd8(
        vx, d840dx, d60dx, d12dx, d2dx, i, imin, imax,
        read_lopsided_sample(f, ex, soa, i - 5, j, k),
        read_lopsided_sample(f, ex, soa, i - 4, j, k),
        read_lopsided_sample(f, ex, soa, i - 3, j, k),
        read_lopsided_sample(f, ex, soa, i - 2, j, k),
        read_lopsided_sample(f, ex, soa, i - 1, j, k),
        center_value,
        read_lopsided_sample(f, ex, soa, i + 1, j, k),
        read_lopsided_sample(f, ex, soa, i + 2, j, k),
        read_lopsided_sample(f, ex, soa, i + 3, j, k),
        read_lopsided_sample(f, ex, soa, i + 4, j, k),
        read_lopsided_sample(f, ex, soa, i + 5, j, k)
    );

    const double rhs_add_y = upwind_bias_contribution_fd8(
        vy, d840dy, d60dy, d12dy, d2dy, j, jmin, jmax,
        read_lopsided_sample(f, ex, soa, i, j - 5, k),
        read_lopsided_sample(f, ex, soa, i, j - 4, k),
        read_lopsided_sample(f, ex, soa, i, j - 3, k),
        read_lopsided_sample(f, ex, soa, i, j - 2, k),
        read_lopsided_sample(f, ex, soa, i, j - 1, k),
        center_value,
        read_lopsided_sample(f, ex, soa, i, j + 1, k),
        read_lopsided_sample(f, ex, soa, i, j + 2, k),
        read_lopsided_sample(f, ex, soa, i, j + 3, k),
        read_lopsided_sample(f, ex, soa, i, j + 4, k),
        read_lopsided_sample(f, ex, soa, i, j + 5, k)
    );

    const double rhs_add_z = upwind_bias_contribution_fd8(
        vz, d840dz, d60dz, d12dz, d2dz, k, kmin, kmax,
        read_lopsided_sample(f, ex, soa, i, j, k - 5),
        read_lopsided_sample(f, ex, soa, i, j, k - 4),
        read_lopsided_sample(f, ex, soa, i, j, k - 3),
        read_lopsided_sample(f, ex, soa, i, j, k - 2),
        read_lopsided_sample(f, ex, soa, i, j, k - 1),
        center_value,
        read_lopsided_sample(f, ex, soa, i, j, k + 1),
        read_lopsided_sample(f, ex, soa, i, j, k + 2),
        read_lopsided_sample(f, ex, soa, i, j, k + 3),
        read_lopsided_sample(f, ex, soa, i, j, k + 4),
        read_lopsided_sample(f, ex, soa, i, j, k + 5)
    );
#else
#error Unsupported ghost_width in d_lopsided_point rhs assembly; expected ghost_width == 2, 3, 4, or 5.
#endif

    (void)f_rhs;
    return rhs_add_x + rhs_add_y + rhs_add_z;
}

#endif // LOPISIDEDIFF_GPU_CUH
