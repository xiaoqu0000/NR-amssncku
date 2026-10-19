#include "bssn_rhs.h"

#include <cuda_runtime.h>

#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

void check_cuda(cudaError_t status, const char* operation) {
    if (status != cudaSuccess) {
        std::cerr << operation << " failed: " << cudaGetErrorString(status) << std::endl;
        std::exit(EXIT_FAILURE);
    }
}

int parse_int_arg(int argc, char** argv, int index, int default_value) {
    if (argc <= index) return default_value;
    return std::atoi(argv[index]);
}

double parse_double_arg(int argc, char** argv, int index, double default_value) {
    if (argc <= index) return default_value;
    return std::atof(argv[index]);
}

void allocate_zero_field(double** field, size_t byte_count, std::vector<double*>& allocations) {
    check_cuda(cudaMalloc(reinterpret_cast<void**>(field), byte_count), "cudaMalloc(field)");
    check_cuda(cudaMemset(*field, 0, byte_count), "cudaMemset(field)");
    allocations.push_back(*field);
}

void initialize_axis(double** field, int point_count, double spacing, std::vector<double*>& allocations) {
    std::vector<double> host_values(point_count, 0.0);
    for (int index = 0; index < point_count; ++index) host_values[index] = spacing * static_cast<double>(index);
    check_cuda(cudaMalloc(reinterpret_cast<void**>(field), point_count * sizeof(double)), "cudaMalloc(axis)");
    check_cuda(cudaMemcpy(*field, host_values.data(), point_count * sizeof(double), cudaMemcpyHostToDevice), "cudaMemcpy(axis)");
    allocations.push_back(*field);
}

}  

int main(int argc, char** argv) {
    int nx = parse_int_arg(argc, argv, 1, 64);
    int ny = parse_int_arg(argc, argv, 2, nx);
    int nz = parse_int_arg(argc, argv, 3, nx);
    int warmup_count = parse_int_arg(argc, argv, 4, 10);
    int repeat_count = parse_int_arg(argc, argv, 5, 50);
    int co = parse_int_arg(argc, argv, 6, 1);
    int symmetry = parse_int_arg(argc, argv, 7, 0);
    int lev = parse_int_arg(argc, argv, 8, 0);
    double eps = parse_double_arg(argc, argv, 9, 0.0);

    if (nx < 4 || ny < 4 || nz < 4 || warmup_count < 0 || repeat_count <= 0) {
        std::cerr << "Usage: ./rhs_gpu_bench [nx] [ny] [nz] [warmup] [repeat] [co] [symmetry] [lev] [eps]" << std::endl;
        return EXIT_FAILURE;
    }

    check_cuda(cudaSetDevice(0), "cudaSetDevice");

    cudaDeviceProp device_properties;
    check_cuda(cudaGetDeviceProperties(&device_properties, 0), "cudaGetDeviceProperties");

    int extents[3] = {nx, ny, nz};
    int cell_count = nx * ny * nz;
    size_t field_bytes = static_cast<size_t>(cell_count) * sizeof(double);

    std::vector<double*> allocations;

    double *X = nullptr, *Y = nullptr, *Z = nullptr;
    initialize_axis(&X, nx, 1.0, allocations);
    initialize_axis(&Y, ny, 1.0, allocations);
    initialize_axis(&Z, nz, 1.0, allocations);

    double *chi = nullptr, *trK = nullptr;
    double *dxx = nullptr, *gxy = nullptr, *gxz = nullptr, *dyy = nullptr, *gyz = nullptr, *dzz = nullptr;
    double *Axx = nullptr, *Axy = nullptr, *Axz = nullptr, *Ayy = nullptr, *Ayz = nullptr, *Azz = nullptr;
    double *Gamx = nullptr, *Gamy = nullptr, *Gamz = nullptr, *Lap = nullptr;
    double *betax = nullptr, *betay = nullptr, *betaz = nullptr;
    double *dtSfx = nullptr, *dtSfy = nullptr, *dtSfz = nullptr;

    double *chi_rhs = nullptr, *trK_rhs = nullptr;
    double *gxx_rhs = nullptr, *gxy_rhs = nullptr, *gxz_rhs = nullptr, *gyy_rhs = nullptr, *gyz_rhs = nullptr, *gzz_rhs = nullptr;
    double *Axx_rhs = nullptr, *Axy_rhs = nullptr, *Axz_rhs = nullptr, *Ayy_rhs = nullptr, *Ayz_rhs = nullptr, *Azz_rhs = nullptr;
    double *Gamx_rhs = nullptr, *Gamy_rhs = nullptr, *Gamz_rhs = nullptr, *Lap_rhs = nullptr;
    double *betax_rhs = nullptr, *betay_rhs = nullptr, *betaz_rhs = nullptr;
    double *dtSfx_rhs = nullptr, *dtSfy_rhs = nullptr, *dtSfz_rhs = nullptr;

    double *rho = nullptr, *Sx = nullptr, *Sy = nullptr, *Sz = nullptr;
    double *Sxx = nullptr, *Sxy = nullptr, *Sxz = nullptr, *Syy = nullptr, *Syz = nullptr, *Szz = nullptr;

    double *Gamxxx = nullptr, *Gamxxy = nullptr, *Gamxxz = nullptr, *Gamxyy = nullptr, *Gamxyz = nullptr, *Gamxzz = nullptr;
    double *Gamyxx = nullptr, *Gamyxy = nullptr, *Gamyxz = nullptr, *Gamyyy = nullptr, *Gamyyz = nullptr, *Gamyzz = nullptr;
    double *Gamzxx = nullptr, *Gamzxy = nullptr, *Gamzxz = nullptr, *Gamzyy = nullptr, *Gamzyz = nullptr, *Gamzzz = nullptr;
    double *Rxx = nullptr, *Rxy = nullptr, *Rxz = nullptr, *Ryy = nullptr, *Ryz = nullptr, *Rzz = nullptr;
    double *ham_Res = nullptr, *movx_Res = nullptr, *movy_Res = nullptr, *movz_Res = nullptr, *Gmx_Res = nullptr, *Gmy_Res = nullptr, *Gmz_Res = nullptr;

    allocate_zero_field(&chi, field_bytes, allocations);
    allocate_zero_field(&trK, field_bytes, allocations);
    allocate_zero_field(&dxx, field_bytes, allocations);
    allocate_zero_field(&gxy, field_bytes, allocations);
    allocate_zero_field(&gxz, field_bytes, allocations);
    allocate_zero_field(&dyy, field_bytes, allocations);
    allocate_zero_field(&gyz, field_bytes, allocations);
    allocate_zero_field(&dzz, field_bytes, allocations);
    allocate_zero_field(&Axx, field_bytes, allocations);
    allocate_zero_field(&Axy, field_bytes, allocations);
    allocate_zero_field(&Axz, field_bytes, allocations);
    allocate_zero_field(&Ayy, field_bytes, allocations);
    allocate_zero_field(&Ayz, field_bytes, allocations);
    allocate_zero_field(&Azz, field_bytes, allocations);
    allocate_zero_field(&Gamx, field_bytes, allocations);
    allocate_zero_field(&Gamy, field_bytes, allocations);
    allocate_zero_field(&Gamz, field_bytes, allocations);
    allocate_zero_field(&Lap, field_bytes, allocations);
    allocate_zero_field(&betax, field_bytes, allocations);
    allocate_zero_field(&betay, field_bytes, allocations);
    allocate_zero_field(&betaz, field_bytes, allocations);
    allocate_zero_field(&dtSfx, field_bytes, allocations);
    allocate_zero_field(&dtSfy, field_bytes, allocations);
    allocate_zero_field(&dtSfz, field_bytes, allocations);
    allocate_zero_field(&chi_rhs, field_bytes, allocations);
    allocate_zero_field(&trK_rhs, field_bytes, allocations);
    allocate_zero_field(&gxx_rhs, field_bytes, allocations);
    allocate_zero_field(&gxy_rhs, field_bytes, allocations);
    allocate_zero_field(&gxz_rhs, field_bytes, allocations);
    allocate_zero_field(&gyy_rhs, field_bytes, allocations);
    allocate_zero_field(&gyz_rhs, field_bytes, allocations);
    allocate_zero_field(&gzz_rhs, field_bytes, allocations);
    allocate_zero_field(&Axx_rhs, field_bytes, allocations);
    allocate_zero_field(&Axy_rhs, field_bytes, allocations);
    allocate_zero_field(&Axz_rhs, field_bytes, allocations);
    allocate_zero_field(&Ayy_rhs, field_bytes, allocations);
    allocate_zero_field(&Ayz_rhs, field_bytes, allocations);
    allocate_zero_field(&Azz_rhs, field_bytes, allocations);
    allocate_zero_field(&Gamx_rhs, field_bytes, allocations);
    allocate_zero_field(&Gamy_rhs, field_bytes, allocations);
    allocate_zero_field(&Gamz_rhs, field_bytes, allocations);
    allocate_zero_field(&Lap_rhs, field_bytes, allocations);
    allocate_zero_field(&betax_rhs, field_bytes, allocations);
    allocate_zero_field(&betay_rhs, field_bytes, allocations);
    allocate_zero_field(&betaz_rhs, field_bytes, allocations);
    allocate_zero_field(&dtSfx_rhs, field_bytes, allocations);
    allocate_zero_field(&dtSfy_rhs, field_bytes, allocations);
    allocate_zero_field(&dtSfz_rhs, field_bytes, allocations);
    allocate_zero_field(&rho, field_bytes, allocations);
    allocate_zero_field(&Sx, field_bytes, allocations);
    allocate_zero_field(&Sy, field_bytes, allocations);
    allocate_zero_field(&Sz, field_bytes, allocations);
    allocate_zero_field(&Sxx, field_bytes, allocations);
    allocate_zero_field(&Sxy, field_bytes, allocations);
    allocate_zero_field(&Sxz, field_bytes, allocations);
    allocate_zero_field(&Syy, field_bytes, allocations);
    allocate_zero_field(&Syz, field_bytes, allocations);
    allocate_zero_field(&Szz, field_bytes, allocations);
    allocate_zero_field(&Gamxxx, field_bytes, allocations);
    allocate_zero_field(&Gamxxy, field_bytes, allocations);
    allocate_zero_field(&Gamxxz, field_bytes, allocations);
    allocate_zero_field(&Gamxyy, field_bytes, allocations);
    allocate_zero_field(&Gamxyz, field_bytes, allocations);
    allocate_zero_field(&Gamxzz, field_bytes, allocations);
    allocate_zero_field(&Gamyxx, field_bytes, allocations);
    allocate_zero_field(&Gamyxy, field_bytes, allocations);
    allocate_zero_field(&Gamyxz, field_bytes, allocations);
    allocate_zero_field(&Gamyyy, field_bytes, allocations);
    allocate_zero_field(&Gamyyz, field_bytes, allocations);
    allocate_zero_field(&Gamyzz, field_bytes, allocations);
    allocate_zero_field(&Gamzxx, field_bytes, allocations);
    allocate_zero_field(&Gamzxy, field_bytes, allocations);
    allocate_zero_field(&Gamzxz, field_bytes, allocations);
    allocate_zero_field(&Gamzyy, field_bytes, allocations);
    allocate_zero_field(&Gamzyz, field_bytes, allocations);
    allocate_zero_field(&Gamzzz, field_bytes, allocations);
    allocate_zero_field(&Rxx, field_bytes, allocations);
    allocate_zero_field(&Rxy, field_bytes, allocations);
    allocate_zero_field(&Rxz, field_bytes, allocations);
    allocate_zero_field(&Ryy, field_bytes, allocations);
    allocate_zero_field(&Ryz, field_bytes, allocations);
    allocate_zero_field(&Rzz, field_bytes, allocations);
    allocate_zero_field(&ham_Res, field_bytes, allocations);
    allocate_zero_field(&movx_Res, field_bytes, allocations);
    allocate_zero_field(&movy_Res, field_bytes, allocations);
    allocate_zero_field(&movz_Res, field_bytes, allocations);
    allocate_zero_field(&Gmx_Res, field_bytes, allocations);
    allocate_zero_field(&Gmy_Res, field_bytes, allocations);
    allocate_zero_field(&Gmz_Res, field_bytes, allocations);

    cudaStream_t stream = nullptr;
    cudaEvent_t start_event = nullptr;
    cudaEvent_t stop_event = nullptr;

    check_cuda(cudaStreamCreate(&stream), "cudaStreamCreate");
    check_cuda(cudaEventCreate(&start_event), "cudaEventCreate(start)");
    check_cuda(cudaEventCreate(&stop_event), "cudaEventCreate(stop)");

    std::cout << std::fixed << std::setprecision(3);
    std::cout << "GPU            : " << device_properties.name << std::endl;
    std::cout << "Grid           : " << nx << " x " << ny << " x " << nz << std::endl;
    std::cout << "Cells          : " << cell_count << std::endl;
    std::cout << "Warmup         : " << warmup_count << std::endl;
    std::cout << "Repeat         : " << repeat_count << std::endl;
    std::cout << "Constraints    : " << (co == 0 ? "on" : "off") << std::endl;
    std::cout << "Symmetry       : " << symmetry << std::endl;
    std::cout << "Level          : " << lev << std::endl;
    std::cout << "Dissipation    : " << eps << std::endl;

    for (int iteration = 0; iteration < warmup_count; ++iteration) {
        gpu_compute_rhs_bssn_launch(
            stream, extents, 0.0, X, Y, Z,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gamx, Gamy, Gamz, Lap,
            betax, betay, betaz,
            dtSfx, dtSfy, dtSfz,
            chi_rhs, trK_rhs,
            gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
            Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
            Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
            betax_rhs, betay_rhs, betaz_rhs,
            dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
            rho, Sx, Sy, Sz,
            Sxx, Sxy, Sxz, Syy, Syz, Szz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res,
            symmetry, lev, eps, co
        );
    }
    check_cuda(cudaStreamSynchronize(stream), "cudaStreamSynchronize(warmup)");

    check_cuda(cudaEventRecord(start_event, stream), "cudaEventRecord(start)");
    for (int iteration = 0; iteration < repeat_count; ++iteration) {
        gpu_compute_rhs_bssn_launch(
            stream, extents, 0.0, X, Y, Z,
            chi, trK, dxx, gxy, gxz, dyy, gyz, dzz,
            Axx, Axy, Axz, Ayy, Ayz, Azz,
            Gamx, Gamy, Gamz, Lap,
            betax, betay, betaz,
            dtSfx, dtSfy, dtSfz,
            chi_rhs, trK_rhs,
            gxx_rhs, gxy_rhs, gxz_rhs, gyy_rhs, gyz_rhs, gzz_rhs,
            Axx_rhs, Axy_rhs, Axz_rhs, Ayy_rhs, Ayz_rhs, Azz_rhs,
            Gamx_rhs, Gamy_rhs, Gamz_rhs, Lap_rhs,
            betax_rhs, betay_rhs, betaz_rhs,
            dtSfx_rhs, dtSfy_rhs, dtSfz_rhs,
            rho, Sx, Sy, Sz,
            Sxx, Sxy, Sxz, Syy, Syz, Szz,
            Gamxxx, Gamxxy, Gamxxz, Gamxyy, Gamxyz, Gamxzz,
            Gamyxx, Gamyxy, Gamyxz, Gamyyy, Gamyyz, Gamyzz,
            Gamzxx, Gamzxy, Gamzxz, Gamzyy, Gamzyz, Gamzzz,
            Rxx, Rxy, Rxz, Ryy, Ryz, Rzz,
            ham_Res, movx_Res, movy_Res, movz_Res,
            Gmx_Res, Gmy_Res, Gmz_Res,
            symmetry, lev, eps, co
        );
    }
    check_cuda(cudaEventRecord(stop_event, stream), "cudaEventRecord(stop)");
    check_cuda(cudaEventSynchronize(stop_event), "cudaEventSynchronize(stop)");

    float total_milliseconds = 0.0f;
    check_cuda(cudaEventElapsedTime(&total_milliseconds, start_event, stop_event), "cudaEventElapsedTime");

    double chi_rhs_sample = 0.0;
    double gamx_rhs_sample = 0.0;
    check_cuda(cudaMemcpy(&chi_rhs_sample, chi_rhs + cell_count / 2, sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy(chi_rhs_sample)");
    check_cuda(cudaMemcpy(&gamx_rhs_sample, Gamx_rhs + cell_count / 2, sizeof(double), cudaMemcpyDeviceToHost), "cudaMemcpy(gamx_rhs_sample)");

    double average_milliseconds = total_milliseconds / static_cast<double>(repeat_count);
    double cell_updates_per_second = static_cast<double>(cell_count) * static_cast<double>(repeat_count) / (static_cast<double>(total_milliseconds) * 1.0e-3);

    std::cout << "Average RHS call : " << average_milliseconds << " ms" << std::endl;
    std::cout << "Total elapsed    : " << total_milliseconds << " ms" << std::endl;
    std::cout << "Throughput       : " << cell_updates_per_second / 1.0e6 << " Mcells/s" << std::endl;
    std::cout << "Sample chi_rhs   : " << chi_rhs_sample << std::endl;
    std::cout << "Sample Gamx_rhs  : " << gamx_rhs_sample << std::endl;

    check_cuda(cudaEventDestroy(start_event), "cudaEventDestroy(start)");
    check_cuda(cudaEventDestroy(stop_event), "cudaEventDestroy(stop)");
    check_cuda(cudaStreamDestroy(stream), "cudaStreamDestroy");

    for (double* allocation : allocations) check_cuda(cudaFree(allocation), "cudaFree");

    return EXIT_SUCCESS;
}