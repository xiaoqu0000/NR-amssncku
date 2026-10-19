#include "macrodef.h"
#include "surface_integral.h"

#include <cuda_runtime.h>
#include <device_functions.h>
#include <device_launch_parameters.h>

#include "helper.h"
#include "fadmquantites_bssn.h"
#include "misc_gpu.cuh"
#include "getnpem2_gpu.cuh"
#include "getnpem2.h"
#include "reduce_helper.h"

#define PI M_PI

__global__ void scale_normals_kernel(
    int n_tot, double rex, 
    const double* d_nx, const double* d_ny, const double* d_nz,
    double* d_pox_0, double* d_pox_1, double* d_pox_2
) {
    int n = blockIdx.x * blockDim.x + threadIdx.x;
    if (n < n_tot) {
        d_pox_0[n] = rex * d_nx[n];
        d_pox_1[n] = rex * d_ny[n];
        d_pox_2[n] = rex * d_nz[n];
    }
}

void gpu_scale_normals_launch(
    cudaStream_t stream, int n_tot, double rex, 
    const double* d_nx, const double* d_ny, const double* d_nz,
    double* d_pox_0, double* d_pox_1, double* d_pox_2
) {
    int blockSize = 256;
    int gridSize = (n_tot + blockSize - 1) / blockSize;
    
    scale_normals_kernel<<<gridSize, blockSize, 0, stream>>>(
        n_tot, rex, d_nx, d_ny, d_nz, d_pox_0, d_pox_1, d_pox_2
    );
}

__global__ void surf_MassPAng_kernel(
    int Nmin, int Nmax, int N_phi, int InList, int Symmetry,
    const double* d_shellf, 
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_nx, const double* d_ny, const double* d_nz,
    const double* d_wtcostheta,
    double* d_reductions
) {
    int n = Nmin + blockIdx.x * blockDim.x + threadIdx.x;
    
    // 线程局部的累加器
    double mass = 0, px = 0, py = 0, pz = 0, sx = 0, sy = 0, sz = 0;

    if (n <= Nmax) {
        int i = n / N_phi; 
        
        // 提取物理量
        double Chi = d_shellf[InList * n + 3];
        double TRK = d_shellf[InList * n + 4];
        double Gxx = d_shellf[InList * n + 5] + 1.0;
        double Gxy = d_shellf[InList * n + 6];
        double Gxz = d_shellf[InList * n + 7];
        double Gyy = d_shellf[InList * n + 8] + 1.0;
        double Gyz = d_shellf[InList * n + 9];
        double Gzz = d_shellf[InList * n + 10] + 1.0;
        double axx = d_shellf[InList * n + 11];
        double axy = d_shellf[InList * n + 12];
        double axz = d_shellf[InList * n + 13];
        double ayy = d_shellf[InList * n + 14];
        double ayz = d_shellf[InList * n + 15];
        double azz = d_shellf[InList * n + 16];

        Chi = 1.0 / (1.0 + Chi); // exp(4*phi)
        double Psi = Chi * sqrt(Chi); // Psi^6

        double wt = d_wtcostheta[i];
        double nx = d_nx[n];
        double ny = d_ny[n];
        double nz = d_nz[n];
        
        // 1. 计算 Mass
#ifdef GaussInt
        mass = (d_shellf[InList * n] * nx + d_shellf[InList * n + 1] * ny + d_shellf[InList * n + 2] * nz) * wt;
#else
        mass = (d_shellf[InList * n] * nx + d_shellf[InList * n + 1] * ny + d_shellf[InList * n + 2] * nz);
#endif

        // 2. 逆度规运算
        double gupzz_denom = Gxx * Gyy * Gzz + Gxy * Gyz * Gxz + Gxz * Gxy * Gyz -
                             Gxz * Gyy * Gxz - Gxy * Gxy * Gzz - Gxx * Gyz * Gyz;
        double gupxx = (Gyy * Gzz - Gyz * Gyz) / gupzz_denom;
        double gupxy = -(Gxy * Gzz - Gyz * Gxz) / gupzz_denom;
        double gupxz = (Gxy * Gyz - Gyy * Gxz) / gupzz_denom;
        double gupyy = (Gxx * Gzz - Gxz * Gxz) / gupzz_denom;
        double gupyz = -(Gxx * Gyz - Gxy * Gxz) / gupzz_denom;
        double gupzz = (Gxx * Gyy - Gxy * Gxy) / gupzz_denom;

        double aupxx = gupxx * axx + gupxy * axy + gupxz * axz;
        double aupxy = gupxx * axy + gupxy * ayy + gupxz * ayz;
        double aupxz = gupxx * axz + gupxy * ayz + gupxz * azz;
        double aupyx = gupxy * axx + gupyy * axy + gupyz * axz;
        double aupyy = gupxy * axy + gupyy * ayy + gupyz * ayz;
        double aupyz = gupxy * axz + gupyy * ayz + gupyz * azz;
        double aupzx = gupxz * axx + gupyz * axy + gupzz * axz;
        double aupzy = gupxz * axy + gupyz * ayy + gupzz * ayz;
        double aupzz = gupxz * axz + gupyz * ayz + gupzz * azz;

        double px_coord = d_pox_0[n];
        double py_coord = d_pox_1[n];
        double pz_coord = d_pox_2[n];
        const double f1o8 = 0.125;

        // 3. 计算 Angular Momentum
        if (Symmetry == 0) {
#ifdef GaussInt
            sx = f1o8 * Psi * (nx * (py_coord * aupxz - pz_coord * aupxy) + ny * (py_coord * aupyz - pz_coord * aupyy) + nz * (py_coord * aupzz - pz_coord * aupzy)) * wt;
            sy = f1o8 * Psi * (nx * (pz_coord * aupxx - px_coord * aupxz) + ny * (pz_coord * aupyx - px_coord * aupyz) + nz * (pz_coord * aupzx - px_coord * aupzz)) * wt;
            sz = f1o8 * Psi * (nx * (px_coord * aupxy - py_coord * aupxx) + ny * (px_coord * aupyy - py_coord * aupyx) + nz * (px_coord * aupzy - py_coord * aupzx)) * wt;
#else
            sx = f1o8 * Psi * (nx * (py_coord * aupxz - pz_coord * aupxy) + ny * (py_coord * aupyz - pz_coord * aupyy) + nz * (py_coord * aupzz - pz_coord * aupzy));
            sy = f1o8 * Psi * (nx * (pz_coord * aupxx - px_coord * aupxz) + ny * (pz_coord * aupyx - px_coord * aupyz) + nz * (pz_coord * aupzx - px_coord * aupzz));
            sz = f1o8 * Psi * (nx * (px_coord * aupxy - py_coord * aupxx) + ny * (px_coord * aupyy - py_coord * aupyx) + nz * (px_coord * aupzy - py_coord * aupzx));
#endif
        } else if (Symmetry == 1) {
#ifdef GaussInt
            sz = f1o8 * Psi * (nx * (px_coord * aupxy - py_coord * aupxx) + ny * (px_coord * aupyy - py_coord * aupyx) + nz * (px_coord * aupzy - py_coord * aupzx)) * wt;
#else
            sz = f1o8 * Psi * (nx * (px_coord * aupxy - py_coord * aupxx) + ny * (px_coord * aupyy - py_coord * aupyx) + nz * (px_coord * aupzy - py_coord * aupzx));
#endif
        }

        // 4. 计算 Linear Momentum
        axx = Chi * (axx + Gxx * TRK / 3.0) - TRK;
        axy = Chi * (axy + Gxy * TRK / 3.0);
        axz = Chi * (axz + Gxz * TRK / 3.0);
        ayy = Chi * (ayy + Gyy * TRK / 3.0) - TRK;
        ayz = Chi * (ayz + Gyz * TRK / 3.0);
        azz = Chi * (azz + Gzz * TRK / 3.0) - TRK;

        if (Symmetry == 0) {
#ifdef GaussInt
            px = f1o8 * Psi * (nx * axx + ny * axy + nz * axz) * wt;
            py = f1o8 * Psi * (nx * axy + ny * ayy + nz * ayz) * wt;
            pz = f1o8 * Psi * (nx * axz + ny * ayz + nz * azz) * wt;
#else
            px = f1o8 * Psi * (nx * axx + ny * axy + nz * axz);
            py = f1o8 * Psi * (nx * axy + ny * ayy + nz * ayz);
            pz = f1o8 * Psi * (nx * axz + ny * ayz + nz * azz);
#endif
        } else if (Symmetry == 1) {
#ifdef GaussInt
            px = f1o8 * Psi * (nx * axx + ny * axy + nz * axz) * wt;
            py = f1o8 * Psi * (nx * axy + ny * ayy + nz * ayz) * wt;
#else
            px = f1o8 * Psi * (nx * axx + ny * axy + nz * axz);
            py = f1o8 * Psi * (nx * axy + ny * ayy + nz * ayz);
#endif
        }
    }

    __shared__ double s_red[7];
    if (threadIdx.x == 0) {
        for (int k = 0; k < 7; k++) s_red[k] = 0.0;
    }
    __syncthreads();

    atomicAdd(&s_red[0], mass);
    atomicAdd(&s_red[1], px);
    atomicAdd(&s_red[2], py);
    atomicAdd(&s_red[3], pz);
    atomicAdd(&s_red[4], sx);
    atomicAdd(&s_red[5], sy);
    atomicAdd(&s_red[6], sz);
    __syncthreads();

    if (threadIdx.x == 0) {
        atomicAdd(&d_reductions[0], s_red[0]);
        atomicAdd(&d_reductions[1], s_red[1]);
        atomicAdd(&d_reductions[2], s_red[2]);
        atomicAdd(&d_reductions[3], s_red[3]);
        atomicAdd(&d_reductions[4], s_red[4]);
        atomicAdd(&d_reductions[5], s_red[5]);
        atomicAdd(&d_reductions[6], s_red[6]);
    }
}

void gpu_surf_MassPAng_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int InList, int Symmetry,
    const double* d_shellf, 
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_nx, const double* d_ny, const double* d_nz,
    const double* d_wtcostheta, double* d_reductions
) {
    int num_elements = Nmax - Nmin + 1;
    if (num_elements <= 0) return;

    int blockSize = 256;
    int gridSize = (num_elements + blockSize - 1) / blockSize;

    surf_MassPAng_kernel<<<gridSize, blockSize, 0, stream>>>(
        Nmin, Nmax, N_phi, InList, Symmetry,
        d_shellf, 
        d_pox_0, d_pox_1, d_pox_2,
        d_nx, d_ny, d_nz, d_wtcostheta,
        d_reductions
    );
}

void surface_integral::gpu_surf_MassPAng(
    double rex, int lev, cgh *GH, var *chi, var *trK,
    var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
    var *Gmx, var *Gmy, var *Gmz,
    var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs, // temparay memory for mass^i
    double *Rout, monitor *Monitor
) {
    if (myrank == 0 && GH->grids[lev] != 1)
        if (Monitor && Monitor->outfile)
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        else
            cout << "WARNING: surface integral on multipatches" << endl;

    const int InList = 17;

    MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
    DG_List->insert(Sfy_rhs); DG_List->insert(Sfz_rhs);
    DG_List->insert(chi); DG_List->insert(trK);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);
    DG_List->insert(Axx); DG_List->insert(Axy); DG_List->insert(Axz);
    DG_List->insert(Ayy); DG_List->insert(Ayz); DG_List->insert(Azz);

    MyList<Patch> *Pp = GH->PatL[lev];
    while (Pp) {
        MyList<Block> *BP = Pp->data->blb;
        while (BP) {
            Block *cg = BP->data;
            if (myrank == cg->rank) {
                gpu_admmass_bssn_launch(
                    cg->stream,
                    cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                    cg->d_fgfs[chi->sgfn], cg->d_fgfs[trK->sgfn],
                    cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn], cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                    cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn], cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn],
                    cg->d_fgfs[Gmx->sgfn], cg->d_fgfs[Gmy->sgfn], cg->d_fgfs[Gmz->sgfn],
                    cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
                    Symmetry
                );
            }
            if (BP == Pp->data->ble)
                break;
            BP = BP->next;
        }
        Pp = Pp->next;
    }
    GPUManager::getInstance().synchronize_all();

    double *d_pox[3];
    for (int i = 0; i < 3; i++) d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp, Lp, Nmin, Nmax;
    mp = n_tot / cpusize;
    Lp = n_tot - cpusize * mp;

    if (Lp > myrank) {
        Nmin = myrank * mp + myrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = myrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double* d_reductions = GPUManager::getInstance().allocate_comm_buffer<double>(7, this->stream);
    GPUManager::getInstance().memset_device_memory(d_reductions, 7, this->stream);
    gpu_surf_MassPAng_launch(
        this->stream,
        Nmin, Nmax, N_phi, InList, Symmetry,
        d_shellf, d_pox[0], d_pox[1], d_pox[2], d_nx_g, d_ny_g, d_nz_g, d_wtcostheta, d_reductions
    );

    // 7-double reduction. Only rank 0 reads Rout (MAPMonitor->writefile is
    // rank-0-only), so manual gather + local sum on root is fine and noticeably
    // faster than MPI_Allreduce/Reduce on small N + small msg.
    double h_reductions[7] = {0};
    GPUManager::getInstance().sync_to_cpu(h_reductions, d_reductions, 7, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(h_reductions, 7, 0, MPI_COMM_WORLD);

    // 最后乘上积分面元权重及常数系数
#ifdef GaussInt
    double scale_mass = rex * rex * dphi * factor;
    double scale_ang  = scale_mass * (1.0 / PI);
#else
    double scale_mass = rex * rex * dphi * dcostheta * factor;
    double scale_ang  = scale_mass * (1.0 / PI);
#endif

    int myrank_; MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    if (myrank_ == 0) {
        Rout[0] = h_reductions[0] * scale_mass; // mass
        Rout[1] = h_reductions[1] * scale_ang;  // px
        Rout[2] = h_reductions[2] * scale_ang;  // py
        Rout[3] = h_reductions[3] * scale_ang;  // pz
        Rout[4] = h_reductions[4] * scale_ang;  // sx
        Rout[5] = h_reductions[5] * scale_ang;  // sy
        Rout[6] = h_reductions[6] * scale_ang;  // sz
    }

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_reductions, 7, this->stream);

    DG_List->clearList();
    GPUManager::synchronize_stream(this->stream);
}

__global__ void surf_Wave_kernel(
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const double* d_shellf, 
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double Rpsi4_SoA_0, double Rpsi4_SoA_1, double Rpsi4_SoA_2, // 通过值传递小数组
    double Ipsi4_SoA_0, double Ipsi4_SoA_1, double Ipsi4_SoA_2,
    double* d_RP_out, double* d_IP_out
) {
    int n = Nmin + blockIdx.x * blockDim.x + threadIdx.x;
    if (n > Nmax) return;

    int i = n / N_phi;
    int j = n - i * N_phi;

    int lpsy = 0;
    if (Symmetry == 0) lpsy = 1;
    else if (Symmetry == 1) lpsy = 2;
    else if (Symmetry == 2) lpsy = 8;

    double wt = d_wtcostheta[i];

    int countlm = 0;
    for (int pl = spinw; pl <= maxl; pl++) {
        for (int pm = -pl; pm <= pl; pm++) {
            double costheta_base = d_arcostheta[i];
            double phi_base = (j + 0.5) * dphi;
            double psi4RR_base = d_shellf[2 * n];
            double psi4II_base = d_shellf[2 * n + 1];

            double local_RP = 0.0;
            double local_IP = 0.0;

            for (int lp = 0; lp < lpsy; lp++) {
                double costheta = costheta_base;
                double phi_angle = phi_base;
                double psi4RR = psi4RR_base;
                double psi4II = psi4II_base;

                switch (lp) {
                    case 0: // +++ (theta, phi)
                        break;
                    case 1: // ++- (pi-theta, phi)
                        costheta = -costheta_base;
                        psi4RR *= Rpsi4_SoA_2;
                        psi4II *= Ipsi4_SoA_2;
                        break;
                    case 2: // +-+ (theta, 2*pi-phi) => 等价于 -phi
                        phi_angle = -phi_base; 
                        psi4RR *= Rpsi4_SoA_1;
                        psi4II *= Ipsi4_SoA_1;
                        break;
                    case 3: // +-- (pi-theta, 2*pi-phi)
                        costheta = -costheta_base;
                        phi_angle = -phi_base;
                        psi4RR *= (Rpsi4_SoA_2 * Rpsi4_SoA_1);
                        psi4II *= (Ipsi4_SoA_2 * Ipsi4_SoA_1);
                        break;
                    case 4: // -++ (theta, pi-phi)
                        phi_angle = M_PI - phi_base;
                        psi4RR *= Rpsi4_SoA_0;
                        psi4II *= Ipsi4_SoA_0;
                        break;
                    case 5: // -+- (pi-theta, pi-phi)
                        costheta = -costheta_base;
                        phi_angle = M_PI - phi_base;
                        psi4RR *= (Rpsi4_SoA_2 * Rpsi4_SoA_0);
                        psi4II *= (Ipsi4_SoA_2 * Ipsi4_SoA_0);
                        break;
                    case 6: // --+ (theta, pi+phi)
                        phi_angle = M_PI + phi_base;
                        psi4RR *= (Rpsi4_SoA_1 * Rpsi4_SoA_0);
                        psi4II *= (Ipsi4_SoA_1 * Ipsi4_SoA_0);
                        break;
                    case 7: // --- (pi-theta, pi+phi)
                        costheta = -costheta_base;
                        phi_angle = M_PI + phi_base;
                        psi4RR *= (Rpsi4_SoA_2 * Rpsi4_SoA_1 * Rpsi4_SoA_0);
                        psi4II *= (Ipsi4_SoA_2 * Ipsi4_SoA_1 * Ipsi4_SoA_0);
                        break;
                }

                double cosmphi = cos(pm * phi_angle);
                double sinmphi = sin(pm * phi_angle);

                double thetap = sqrt((2.0 * pl + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(pl, pm, spinw, costheta);

#ifdef GaussInt
                local_RP += thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wt;
                local_IP += thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wt;
#else
                local_RP += thetap * (psi4RR * cosmphi + psi4II * sinmphi);
                local_IP += thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
            }

            atomicAdd(&d_RP_out[countlm], local_RP);
            atomicAdd(&d_IP_out[countlm], local_IP);
            
            countlm++;
        }
    }
}

void gpu_surf_Wave_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const double* d_shellf,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double Rpsi4_SoA_0, double Rpsi4_SoA_1, double Rpsi4_SoA_2,
    double Ipsi4_SoA_0, double Ipsi4_SoA_1, double Ipsi4_SoA_2,
    double* d_RP_out, double* d_IP_out
) {
    int num_elements = Nmax - Nmin + 1;
    if (num_elements <= 0) return;

    int blockSize = 256;
    int gridSize = (num_elements + blockSize - 1) / blockSize;

    surf_Wave_kernel<<<gridSize, blockSize, 0, stream>>>(
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf,
        d_wtcostheta, d_arcostheta,
        dphi,
        Rpsi4_SoA_0, Rpsi4_SoA_1, Rpsi4_SoA_2,
        Ipsi4_SoA_0, Ipsi4_SoA_1, Ipsi4_SoA_2,
        d_RP_out, d_IP_out
    );
}

void surface_integral::gpu_surf_Wave(
    double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    if (myrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile)
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        else
            cout << "WARNING: surface integral on multipatches" << endl;
    }

    const int InList = 2; // Rpsi4 和 Ipsi4

    MyList<var> *DG_List = new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    double *d_pox[3];
    for (int i = 0; i < 3; i++) {
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    }
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp, Lp, Nmin, Nmax;
    mp = n_tot / cpusize;
    Lp = n_tot - cpusize * mp;

    if (Lp > myrank) {
        Nmin = myrank * mp + myrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = myrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    double R_SoA[3] = {Rpsi4->SoA[0], Rpsi4->SoA[1], Rpsi4->SoA[2]};
    double I_SoA[3] = {Ipsi4->SoA[0], Ipsi4->SoA[1], Ipsi4->SoA[2]};

    gpu_surf_Wave_launch(
        this->stream,
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf, d_wtcostheta, d_arcostheta, dphi,
        R_SoA[0], R_SoA[1], R_SoA[2],
        I_SoA[0], I_SoA[1], I_SoA[2],
        d_RP_out, d_IP_out
    );

    // Manual gather to root: only rank 0 reads RP/IP (Psi4Monitor->writefile is
    // rank-0-only), so a point-to-point gather + local sum is cheaper than
    // either MPI_Allreduce or MPI_Reduce on small N + small NN.
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(RP, NN, 0, MPI_COMM_WORLD);
    amss_reduce::gather_sum_to_root(IP, NN, 0, MPI_COMM_WORLD);

    // 根据分配律，最后在目标 CPU 数组上统一乘上常数项即可，避免在规约前浪费 N 次运算开销
#ifdef GaussInt
    double scale = rex * dphi;
#else
    double scale = rex * dphi * dcostheta;
#endif
    int myrank_; MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    if (myrank_ == 0) {
        for (int ii = 0; ii < NN; ii++) {
            RP[ii] *= scale;
            IP[ii] *= scale;
        }
    }

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    DG_List->clearList();
}

// ============================================================
// Tier 1: MPI_Comm variants
// ============================================================

void surface_integral::gpu_surf_Wave(
    double rex, int lev, cgh *GH, var *Rpsi4, var *Ipsi4,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor, MPI_Comm Comm_here
) {
    int lmyrank;
    MPI_Comm_rank(Comm_here, &lmyrank);
    if (lmyrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile)
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        else
            cout << "WARNING: surface integral on multipatches" << endl;
    }

    const int InList = 2;

    MyList<var> *DG_List = new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    double *d_pox[3];
    for (int i = 0; i < 3; i++) {
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    }
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int cpusize_here;
    MPI_Comm_size(Comm_here, &cpusize_here);

    int mp = n_tot / cpusize_here;
    int Lp = n_tot - cpusize_here * mp;

    int Nmin, Nmax;
    if (Lp > lmyrank) {
        Nmin = lmyrank * mp + lmyrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = lmyrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    double R_SoA[3] = {Rpsi4->SoA[0], Rpsi4->SoA[1], Rpsi4->SoA[2]};
    double I_SoA[3] = {Ipsi4->SoA[0], Ipsi4->SoA[1], Ipsi4->SoA[2]};

    gpu_surf_Wave_launch(
        this->stream,
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf, d_wtcostheta, d_arcostheta, dphi,
        R_SoA[0], R_SoA[1], R_SoA[2],
        I_SoA[0], I_SoA[1], I_SoA[2],
        d_RP_out, d_IP_out
    );

#if MPI_CUDA_AWARE
    GPUManager::synchronize_stream(this->stream);
    MPI_Allreduce(MPI_IN_PLACE, d_RP_out, NN, MPI_DOUBLE, MPI_SUM, Comm_here);
    MPI_Allreduce(MPI_IN_PLACE, d_IP_out, NN, MPI_DOUBLE, MPI_SUM, Comm_here);
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
#else
    double *h_RP_local = new double[NN];
    double *h_IP_local = new double[NN];
    GPUManager::getInstance().sync_to_cpu(h_RP_local, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(h_IP_local, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
    MPI_Allreduce(h_RP_local, RP, NN, MPI_DOUBLE, MPI_SUM, Comm_here);
    MPI_Allreduce(h_IP_local, IP, NN, MPI_DOUBLE, MPI_SUM, Comm_here);
    delete[] h_RP_local;
    delete[] h_IP_local;
#endif

#ifdef GaussInt
    double scale = rex * dphi;
#else
    double scale = rex * dphi * dcostheta;
#endif
    for (int ii = 0; ii < NN; ii++) {
        RP[ii] *= scale;
        IP[ii] *= scale;
    }

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    DG_List->clearList();
}

// ============================================================
// Tier 1: ShellPatch surf_Wave
// ============================================================

void surface_integral::gpu_surf_Wave(
    double rex, int lev, ShellPatch *GH, var *Rpsi4, var *Ipsi4,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    const int InList = 2;

    MyList<var> *DG_List = new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    double *d_pox[3];
    for (int i = 0; i < 3; i++) {
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    }
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;

    int Nmin, Nmax;
    if (Lp > myrank) {
        Nmin = myrank * mp + myrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = myrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    double R_SoA[3] = {Rpsi4->SoA[0], Rpsi4->SoA[1], Rpsi4->SoA[2]};
    double I_SoA[3] = {Ipsi4->SoA[0], Ipsi4->SoA[1], Ipsi4->SoA[2]};

    gpu_surf_Wave_launch(
        this->stream,
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf, d_wtcostheta, d_arcostheta, dphi,
        R_SoA[0], R_SoA[1], R_SoA[2],
        I_SoA[0], I_SoA[1], I_SoA[2],
        d_RP_out, d_IP_out
    );

    // Manual gather to root: only rank 0 reads RP/IP (Psi4Monitor->writefile
    // is rank-0-only). Cheaper than MPI_Allreduce/Reduce on small N+small NN.
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(RP, NN, 0, MPI_COMM_WORLD);
    amss_reduce::gather_sum_to_root(IP, NN, 0, MPI_COMM_WORLD);

#ifdef GaussInt
    double scale = rex * dphi;
#else
    double scale = rex * dphi * dcostheta;
#endif
    for (int ii = 0; ii < NN; ii++) {
        RP[ii] *= scale;
        IP[ii] *= scale;
    }

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    DG_List->clearList();
}

// ============================================================
// NullShellPatch kernel — different symmetry: no SoA on psi4RR,
// only sign flips on psi4II (cases 1,2,4,7 flip sign)
// ============================================================

__global__ void surf_Wave_NullShell_kernel(
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const double* d_shellf,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int n = Nmin + blockIdx.x * blockDim.x + threadIdx.x;
    if (n > Nmax) return;

    int i = n / N_phi;
    int j = n - i * N_phi;

    int lpsy = 0;
    if (Symmetry == 0) lpsy = 1;
    else if (Symmetry == 1) lpsy = 2;
    else if (Symmetry == 2) lpsy = 8;

    double wt = d_wtcostheta[i];
    double costheta_base = d_arcostheta[i];

    int countlm = 0;
    for (int pl = spinw; pl <= maxl; pl++) {
        for (int pm = -pl; pm <= pl; pm++) {
            double phi_base = (j + 0.5) * dphi;
            double psi4RR_base = d_shellf[2 * n];
            double psi4II_base = d_shellf[2 * n + 1];

            double local_RP = 0.0;
            double local_IP = 0.0;

            for (int lp = 0; lp < lpsy; lp++) {
                double costheta = costheta_base;
                double phi_angle = phi_base;
                double cosmphi, sinmphi;
                double psi4RR = psi4RR_base; // no SoA multiplication
                double psi4II = psi4II_base;

                switch (lp) {
                    case 0:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        break;
                    case 1:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        psi4II = -psi4II_base;
                        break;
                    case 2:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        psi4II = -psi4II_base;
                        break;
                    case 3:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        // psi4II = +psi4II_base (no flip)
                        break;
                    case 4:
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        psi4II = -psi4II_base;
                        break;
                    case 5:
                        costheta = -costheta_base;
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        // psi4II = +psi4II_base (no flip)
                        break;
                    case 6:
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        // psi4II = +psi4II_base (no flip)
                        break;
                    case 7:
                        costheta = -costheta_base;
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        psi4II = -psi4II_base;
                        break;
                }

                double thetap = sqrt((2.0 * pl + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(pl, pm, spinw, costheta);

#ifdef GaussInt
                local_RP += thetap * (psi4RR * cosmphi + psi4II * sinmphi) * wt;
                local_IP += thetap * (psi4II * cosmphi - psi4RR * sinmphi) * wt;
#else
                local_RP += thetap * (psi4RR * cosmphi + psi4II * sinmphi);
                local_IP += thetap * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
            }

            atomicAdd(&d_RP_out[countlm], local_RP);
            atomicAdd(&d_IP_out[countlm], local_IP);

            countlm++;
        }
    }
}

void gpu_surf_Wave_NullShell_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const double* d_shellf,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int num_elements = Nmax - Nmin + 1;
    if (num_elements <= 0) return;

    int blockSize = 256;
    int gridSize = (num_elements + blockSize - 1) / blockSize;

    surf_Wave_NullShell_kernel<<<gridSize, blockSize, 0, stream>>>(
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf,
        d_wtcostheta, d_arcostheta,
        dphi,
        d_RP_out, d_IP_out
    );
}

// ============================================================
// Tier 1: NullShellPatch2 surf_Wave (CPU interp, GPU decompose)
// ============================================================

void surface_integral::gpu_surf_Wave(
    double rex, int lev, NullShellPatch2 *GH, var *Rpsi4, var *Ipsi4,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    const int InList = 2;

    MyList<var> *DG_List = new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    int n;
    double *pox[3];
    for (int i = 0; i < 3; i++)
        pox[i] = new double[n_tot];
    for (n = 0; n < n_tot; n++) {
        pox[0][n] = rex * nx_g[n];
        pox[1][n] = rex * ny_g[n];
        pox[2][n] = rex * nz_g[n];
    }

    double *shellf = new double[n_tot * InList];
    GH->Interp_Points_2D(DG_List, n_tot, pox, shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;

    int Nmin, Nmax;
    if (Lp > myrank) {
        Nmin = myrank * mp + myrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = myrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double *d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GPUManager::getInstance().sync_to_gpu(shellf, d_shellf, n_tot * InList, this->stream);

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    gpu_surf_Wave_NullShell_launch(
        this->stream,
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf, d_wtcostheta, d_arcostheta, dphi,
        d_RP_out, d_IP_out
    );

    // Manual gather to root: only rank 0 reads RP/IP (Psi4Monitor->writefile
    // is rank-0-only). Cheaper than MPI_Allreduce/Reduce on small N+small NN.
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(RP, NN, 0, MPI_COMM_WORLD);
    amss_reduce::gather_sum_to_root(IP, NN, 0, MPI_COMM_WORLD);

    // No rex multiplication for null shell
#ifdef GaussInt
    double scale = dphi;
#else
    double scale = dphi * dcostheta;
#endif
    for (int ii = 0; ii < NN; ii++) {
        RP[ii] *= scale;
        IP[ii] *= scale;
    }

    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    delete[] pox[0]; delete[] pox[1]; delete[] pox[2];
    delete[] shellf;
    DG_List->clearList();
}

// ============================================================
// Tier 1: NullShellPatch surf_Wave (identical to NullShellPatch2)
// ============================================================

void surface_integral::gpu_surf_Wave(
    double rex, int lev, NullShellPatch *GH, var *Rpsi4, var *Ipsi4,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    const int InList = 2;

    MyList<var> *DG_List = new MyList<var>(Rpsi4);
    DG_List->insert(Ipsi4);

    int n;
    double *pox[3];
    for (int i = 0; i < 3; i++)
        pox[i] = new double[n_tot];
    for (n = 0; n < n_tot; n++) {
        pox[0][n] = rex * nx_g[n];
        pox[1][n] = rex * ny_g[n];
        pox[2][n] = rex * nz_g[n];
    }

    double *shellf = new double[n_tot * InList];
    GH->Interp_Points_2D(DG_List, n_tot, pox, shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;

    int Nmin, Nmax;
    if (Lp > myrank) {
        Nmin = myrank * mp + myrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = myrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double *d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GPUManager::getInstance().sync_to_gpu(shellf, d_shellf, n_tot * InList, this->stream);

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    gpu_surf_Wave_NullShell_launch(
        this->stream,
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN,
        d_shellf, d_wtcostheta, d_arcostheta, dphi,
        d_RP_out, d_IP_out
    );

    // Manual gather to root: only rank 0 reads RP/IP (Psi4Monitor->writefile
    // is rank-0-only). Cheaper than MPI_Allreduce/Reduce on small N+small NN.
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(RP, NN, 0, MPI_COMM_WORLD);
    amss_reduce::gather_sum_to_root(IP, NN, 0, MPI_COMM_WORLD);

#ifdef GaussInt
    double scale = dphi;
#else
    double scale = dphi * dcostheta;
#endif
    for (int ii = 0; ii < NN; ii++) {
        RP[ii] *= scale;
        IP[ii] *= scale;
    }

    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    delete[] pox[0]; delete[] pox[1]; delete[] pox[2];
    delete[] shellf;
    DG_List->clearList();
}

// ============================================================
// Tier 1: MPI_Comm MassPAng variant
// ============================================================

void surface_integral::gpu_surf_MassPAng(
    double rex, int lev, cgh *GH, var *chi, var *trK,
    var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
    var *Gmx, var *Gmy, var *Gmz,
    var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
    double *Rout, monitor *Monitor, MPI_Comm Comm_here
) {
    int lmyrank;
    MPI_Comm_rank(Comm_here, &lmyrank);
    if (lmyrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile)
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        else
            cout << "WARNING: surface integral on multipatches" << endl;
    }

    const int InList = 17;

    MyList<Patch> *Pp = GH->PatL[lev];
    while (Pp) {
        MyList<Block> *BP = Pp->data->blb;
        while (BP) {
            Block *cg = BP->data;
            if (myrank == cg->rank) {
                gpu_admmass_bssn_launch(
                    cg->stream,
                    cg->shape, cg->shape_pad, cg->d_X[0], cg->d_X[1], cg->d_X[2],
                    cg->d_fgfs[chi->sgfn], cg->d_fgfs[trK->sgfn],
                    cg->d_fgfs[gxx->sgfn], cg->d_fgfs[gxy->sgfn], cg->d_fgfs[gxz->sgfn], cg->d_fgfs[gyy->sgfn], cg->d_fgfs[gyz->sgfn], cg->d_fgfs[gzz->sgfn],
                    cg->d_fgfs[Axx->sgfn], cg->d_fgfs[Axy->sgfn], cg->d_fgfs[Axz->sgfn], cg->d_fgfs[Ayy->sgfn], cg->d_fgfs[Ayz->sgfn], cg->d_fgfs[Azz->sgfn],
                    cg->d_fgfs[Gmx->sgfn], cg->d_fgfs[Gmy->sgfn], cg->d_fgfs[Gmz->sgfn],
                    cg->d_fgfs[Sfx_rhs->sgfn], cg->d_fgfs[Sfy_rhs->sgfn], cg->d_fgfs[Sfz_rhs->sgfn],
                    Symmetry
                );
            }
            if (BP == Pp->data->ble)
                break;
            BP = BP->next;
        }
        Pp = Pp->next;
    }
    GPUManager::getInstance().synchronize_all();

    MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
    DG_List->insert(Sfy_rhs); DG_List->insert(Sfz_rhs);
    DG_List->insert(chi); DG_List->insert(trK);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);
    DG_List->insert(Axx); DG_List->insert(Axy); DG_List->insert(Axz);
    DG_List->insert(Ayy); DG_List->insert(Ayz); DG_List->insert(Azz);

    double *d_pox[3];
    for (int i = 0; i < 3; i++)
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int cpusize_here;
    MPI_Comm_size(Comm_here, &cpusize_here);

    int mp = n_tot / cpusize_here;
    int Lp = n_tot - cpusize_here * mp;

    int Nmin, Nmax;
    if (Lp > lmyrank) {
        Nmin = lmyrank * mp + lmyrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = lmyrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double* d_reductions = GPUManager::getInstance().allocate_comm_buffer<double>(7, this->stream);
    GPUManager::getInstance().memset_device_memory(d_reductions, 7, this->stream);
    gpu_surf_MassPAng_launch(
        this->stream,
        Nmin, Nmax, N_phi, InList, Symmetry,
        d_shellf, d_pox[0], d_pox[1], d_pox[2], d_nx_g, d_ny_g, d_nz_g, d_wtcostheta, d_reductions
    );

    // 7-double Allreduce — result is consumed on host (Rout scales below).
    // Host-stage unconditionally; cuda-aware Allreduce here adds ~30 µs of
    // rendezvous setup for what's a 56-byte payload.
    double h_reductions[7] = {0};
    double h_reductions_local[7] = {0};
    GPUManager::getInstance().sync_to_cpu(h_reductions_local, d_reductions, 7, this->stream);
    GPUManager::synchronize_stream(this->stream);
    MPI_Allreduce(h_reductions_local, h_reductions, 7, MPI_DOUBLE, MPI_SUM, Comm_here);

#ifdef GaussInt
    double scale_mass = rex * rex * dphi * factor;
    double scale_ang  = scale_mass * (1.0 / PI);
#else
    double scale_mass = rex * rex * dphi * dcostheta * factor;
    double scale_ang  = scale_mass * (1.0 / PI);
#endif

    Rout[0] = h_reductions[0] * scale_mass;
    Rout[1] = h_reductions[1] * scale_ang;
    Rout[2] = h_reductions[2] * scale_ang;
    Rout[3] = h_reductions[3] * scale_ang;
    Rout[4] = h_reductions[4] * scale_ang;
    Rout[5] = h_reductions[5] * scale_ang;
    Rout[6] = h_reductions[6] * scale_ang;

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_reductions, 7, this->stream);

    DG_List->clearList();
    GPUManager::synchronize_stream(this->stream);
}

// ============================================================
// Tier 1: ShellPatch surf_MassPAng
// Uses CPU preprocessing for f_admmass_bssn_ss, GPU for rest
// ============================================================

void surface_integral::gpu_surf_MassPAng(
    double rex, int lev, ShellPatch *GH, var *chi, var *trK,
    var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
    var *Gmx, var *Gmy, var *Gmz,
    var *Sfx_rhs, var *Sfy_rhs, var *Sfz_rhs,
    double *Rout, monitor *Monitor
) {
    if (lev != 0) {
        if (myrank == 0) {
            if (Monitor && Monitor->outfile)
                Monitor->outfile << "WARNING: shell surface integral not on level 0" << endl;
            else
                cout << "WARNING: shell surface integral not on level 0" << endl;
        }
        return;
    }

    // CPU preprocessing: f_admmass_bssn_ss
    MyList<ss_patch> *Pp = GH->PatL;
    while (Pp) {
        MyList<Block> *BL = Pp->data->blb;
        int fngfs = Pp->data->fngfs;
        while (BL) {
            Block *cg = BL->data;
            if (myrank == cg->rank) {
                f_admmass_bssn_ss(cg->shape, cg->X[0], cg->X[1], cg->X[2],
                    cg->fgfs[fngfs + ShellPatch::gx], cg->fgfs[fngfs + ShellPatch::gy], cg->fgfs[fngfs + ShellPatch::gz],
                    cg->fgfs[fngfs + ShellPatch::drhodx], cg->fgfs[fngfs + ShellPatch::drhody], cg->fgfs[fngfs + ShellPatch::drhodz],
                    cg->fgfs[fngfs + ShellPatch::dsigmadx], cg->fgfs[fngfs + ShellPatch::dsigmady], cg->fgfs[fngfs + ShellPatch::dsigmadz],
                    cg->fgfs[fngfs + ShellPatch::dRdx], cg->fgfs[fngfs + ShellPatch::dRdy], cg->fgfs[fngfs + ShellPatch::dRdz],
                    cg->fgfs[fngfs + ShellPatch::drhodxx], cg->fgfs[fngfs + ShellPatch::drhodxy], cg->fgfs[fngfs + ShellPatch::drhodxz],
                    cg->fgfs[fngfs + ShellPatch::drhodyy], cg->fgfs[fngfs + ShellPatch::drhodyz], cg->fgfs[fngfs + ShellPatch::drhodzz],
                    cg->fgfs[fngfs + ShellPatch::dsigmadxx], cg->fgfs[fngfs + ShellPatch::dsigmadxy], cg->fgfs[fngfs + ShellPatch::dsigmadxz],
                    cg->fgfs[fngfs + ShellPatch::dsigmadyy], cg->fgfs[fngfs + ShellPatch::dsigmadyz], cg->fgfs[fngfs + ShellPatch::dsigmadzz],
                    cg->fgfs[fngfs + ShellPatch::dRdxx], cg->fgfs[fngfs + ShellPatch::dRdxy], cg->fgfs[fngfs + ShellPatch::dRdxz],
                    cg->fgfs[fngfs + ShellPatch::dRdyy], cg->fgfs[fngfs + ShellPatch::dRdyz], cg->fgfs[fngfs + ShellPatch::dRdzz],
                    cg->fgfs[chi->sgfn], cg->fgfs[trK->sgfn],
                    cg->fgfs[gxx->sgfn], cg->fgfs[gxy->sgfn], cg->fgfs[gxz->sgfn], cg->fgfs[gyy->sgfn], cg->fgfs[gyz->sgfn], cg->fgfs[gzz->sgfn],
                    cg->fgfs[Axx->sgfn], cg->fgfs[Axy->sgfn], cg->fgfs[Axz->sgfn], cg->fgfs[Ayy->sgfn], cg->fgfs[Ayz->sgfn], cg->fgfs[Azz->sgfn],
                    cg->fgfs[Gmx->sgfn], cg->fgfs[Gmy->sgfn], cg->fgfs[Gmz->sgfn],
                    cg->fgfs[Sfx_rhs->sgfn], cg->fgfs[Sfy_rhs->sgfn], cg->fgfs[Sfz_rhs->sgfn],
                    Symmetry, Pp->data->sst);
            }
            if (BL == Pp->data->ble)
                break;
            BL = BL->next;
        }
        Pp = Pp->next;
    }

    const int InList = 17;

    MyList<var> *DG_List = new MyList<var>(Sfx_rhs);
    DG_List->insert(Sfy_rhs); DG_List->insert(Sfz_rhs);
    DG_List->insert(chi); DG_List->insert(trK);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);
    DG_List->insert(Axx); DG_List->insert(Axy); DG_List->insert(Axz);
    DG_List->insert(Ayy); DG_List->insert(Ayz); DG_List->insert(Azz);

    double *d_pox[3];
    for (int i = 0; i < 3; i++)
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;

    int Nmin, Nmax;
    if (Lp > myrank) {
        Nmin = myrank * mp + myrank;
        Nmax = Nmin + mp;
    } else {
        Nmin = myrank * mp + Lp;
        Nmax = Nmin + mp - 1;
    }

    double* d_reductions = GPUManager::getInstance().allocate_comm_buffer<double>(7, this->stream);
    GPUManager::getInstance().memset_device_memory(d_reductions, 7, this->stream);
    gpu_surf_MassPAng_launch(
        this->stream,
        Nmin, Nmax, N_phi, InList, Symmetry,
        d_shellf, d_pox[0], d_pox[1], d_pox[2], d_nx_g, d_ny_g, d_nz_g, d_wtcostheta, d_reductions
    );

    // 7-double reduction — only rank 0 reads Rout (MAPMonitor->writefile is
    // rank-0-only). Manual gather is cheaper than MPI_Allreduce on small N+msg.
    double h_reductions[7] = {0};
    GPUManager::getInstance().sync_to_cpu(h_reductions, d_reductions, 7, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(h_reductions, 7, 0, MPI_COMM_WORLD);

#ifdef GaussInt
    double scale_mass = rex * rex * dphi * factor;
    double scale_ang  = scale_mass * (1.0 / PI);
#else
    double scale_mass = rex * rex * dphi * dcostheta * factor;
    double scale_ang  = scale_mass * (1.0 / PI);
#endif

    int myrank_; MPI_Comm_rank(MPI_COMM_WORLD, &myrank_);
    if (myrank_ == 0) {
        Rout[0] = h_reductions[0] * scale_mass;
        Rout[1] = h_reductions[1] * scale_ang;
        Rout[2] = h_reductions[2] * scale_ang;
        Rout[3] = h_reductions[3] * scale_ang;
        Rout[4] = h_reductions[4] * scale_ang;
        Rout[5] = h_reductions[5] * scale_ang;
        Rout[6] = h_reductions[6] * scale_ang;
    }

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_reductions, 7, this->stream);

    DG_List->clearList();
    GPUManager::synchronize_stream(this->stream);
}
// ============================================================
// Tier 2: EM Psi4 kernel — 13 variables with getnpem2_point
// ============================================================

__global__ void surf_Wave_EM_kernel(
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const int InList,
    const double* d_shellf,
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int n = Nmin + blockIdx.x * blockDim.x + threadIdx.x;
    if (n > Nmax) return;

    int i = n / N_phi;
    int j = n - i * N_phi;

    int lpsy = 0;
    if (Symmetry == 0) lpsy = 1;
    else if (Symmetry == 1) lpsy = 2;
    else if (Symmetry == 2) lpsy = 8;

    double wt = d_wtcostheta[i];
    double costheta_base = d_arcostheta[i];

    int countlm = 0;
    for (int pl = spinw; pl <= maxl; pl++) {
        for (int pm = -pl; pm <= pl; pm++) {
            double phi_base = (j + 0.5) * dphi;

            double px  = d_pox_0[n];
            double py  = d_pox_1[n];
            double pz  = d_pox_2[n];
            double pEx = d_shellf[InList * n];
            double pEy = d_shellf[InList * n + 1];
            double pEz = d_shellf[InList * n + 2];
            double pBx = d_shellf[InList * n + 3];
            double pBy = d_shellf[InList * n + 4];
            double pBz = d_shellf[InList * n + 5];
            double pchi = d_shellf[InList * n + 6];
            double pgxx = d_shellf[InList * n + 7];
            double pgxy = d_shellf[InList * n + 8];
            double pgxz = d_shellf[InList * n + 9];
            double pgyy = d_shellf[InList * n + 10];
            double pgyz = d_shellf[InList * n + 11];
            double pgzz = d_shellf[InList * n + 12];

            double local_RP = 0.0;
            double local_IP = 0.0;

            for (int lp = 0; lp < lpsy; lp++) {
                double costheta = costheta_base;
                double phi_angle = phi_base;
                double cosmphi, sinmphi;

                double lpx = px, lpy = py, lpz = pz;
                double lEx = pEx, lEy = pEy, lEz = pEz;
                double lBx = pBx, lBy = pBy, lBz = pBz;
                double lgxx = pgxx, lgxy = pgxy, lgxz = pgxz;
                double lgyy = pgyy, lgyz = pgyz, lgzz = pgzz;

                switch (lp) {
                    case 0:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        break;
                    case 1:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpz = -pz; lEz = -pEz;
                        lBx = -pBx; lBy = -pBy;
                        lgxz = -pgxz; lgyz = -pgyz;
                        break;
                    case 2:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        lpy = -py; lEy = -pEy;
                        lBx = -pBx; lBz = -pBz;
                        lgxy = -pgxy; lgyz = -pgyz;
                        break;
                    case 3:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        lpy = -py; lpz = -pz;
                        lEz = -pEz; lBz = -pBz;
                        lgxz = -pgxz; lEy = -pEy;
                        lBy = -pBy; lgxy = -pgxy;
                        break;
                    case 4:
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpx = -px; lEx = -pEx;
                        lBy = -pBy; lBz = -pBz;
                        lgxy = -pgxy; lgxz = -pgxz;
                        break;
                    case 5:
                        costheta = -costheta_base;
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpz = -pz; lpx = -px;
                        lEz = -pEz; lBz = -pBz;
                        lgyz = -pgyz; lEx = -pEx;
                        lBx = -pBx; lgxy = -pgxy;
                        break;
                    case 6:
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpx = -px; lpy = -py;
                        lEx = -pEx; lBx = -pBx;
                        lgxz = -pgxz; lEy = -pEy;
                        lBy = -pBy; lgyz = -pgyz;
                        break;
                    case 7:
                        costheta = -costheta_base;
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpx = -px; lpy = -py; lpz = -pz;
                        lEx = -pEx; lEy = -pEy; lEz = -pEz;
                        break;
                }

                double psi4RR, psi4II;
                getnpem2_point_device(lpx, lpy, lpz, pchi, lgxx, lgxy, lgxz, lgyy, lgyz, lgzz,
                                      lEx, lEy, lEz, lBx, lBy, lBz, psi4RR, psi4II);

                double thetap = sqrt((2.0 * pl + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(pl, pm, spinw, costheta);
                double chi_scale = 1.0 / ((pchi + 1.0) * (pchi + 1.0));

#ifdef GaussInt
                local_RP += thetap * chi_scale * (psi4RR * cosmphi + psi4II * sinmphi) * wt;
                local_IP += thetap * chi_scale * (psi4II * cosmphi - psi4RR * sinmphi) * wt;
#else
                local_RP += thetap * chi_scale * (psi4RR * cosmphi + psi4II * sinmphi);
                local_IP += thetap * chi_scale * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
            }

            atomicAdd(&d_RP_out[countlm], local_RP);
            atomicAdd(&d_IP_out[countlm], local_IP);
            countlm++;
        }
    }
}

// Tier 2-EM1: Phi1 variant of surf_Wave_EM_kernel, using getnpem1_point_device
__global__ void surf_Wave_EM1_kernel(
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const int InList,
    const double* d_shellf,
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int n = Nmin + blockIdx.x * blockDim.x + threadIdx.x;
    if (n > Nmax) return;

    int i = n / N_phi;
    int j = n - i * N_phi;

    int lpsy = 0;
    if (Symmetry == 0) lpsy = 1;
    else if (Symmetry == 1) lpsy = 2;
    else if (Symmetry == 2) lpsy = 8;

    double wt = d_wtcostheta[i];
    double costheta_base = d_arcostheta[i];

    int countlm = 0;
    for (int pl = spinw; pl <= maxl; pl++) {
        for (int pm = -pl; pm <= pl; pm++) {
            double phi_base = (j + 0.5) * dphi;

            double px  = d_pox_0[n];
            double py  = d_pox_1[n];
            double pz  = d_pox_2[n];
            double pEx = d_shellf[InList * n];
            double pEy = d_shellf[InList * n + 1];
            double pEz = d_shellf[InList * n + 2];
            double pBx = d_shellf[InList * n + 3];
            double pBy = d_shellf[InList * n + 4];
            double pBz = d_shellf[InList * n + 5];
            double pchi = d_shellf[InList * n + 6];
            double pgxx = d_shellf[InList * n + 7];
            double pgxy = d_shellf[InList * n + 8];
            double pgxz = d_shellf[InList * n + 9];
            double pgyy = d_shellf[InList * n + 10];
            double pgyz = d_shellf[InList * n + 11];
            double pgzz = d_shellf[InList * n + 12];

            double local_RP = 0.0;
            double local_IP = 0.0;

            for (int lp = 0; lp < lpsy; lp++) {
                double costheta = costheta_base;
                double phi_angle = phi_base;
                double cosmphi, sinmphi;

                double lpx = px, lpy = py, lpz = pz;
                double lEx = pEx, lEy = pEy, lEz = pEz;
                double lBx = pBx, lBy = pBy, lBz = pBz;
                double lgxx = pgxx, lgxy = pgxy, lgxz = pgxz;
                double lgyy = pgyy, lgyz = pgyz, lgzz = pgzz;

                switch (lp) {
                    case 0:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        break;
                    case 1:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpz = -pz; lEz = -pEz;
                        lBx = -pBx; lBy = -pBy;
                        lgxz = -pgxz; lgyz = -pgyz;
                        break;
                    case 2:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        lpy = -py; lEy = -pEy;
                        lBx = -pBx; lBz = -pBz;
                        lgxy = -pgxy; lgyz = -pgyz;
                        break;
                    case 3:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        lpy = -py; lpz = -pz;
                        lEz = -pEz; lBz = -pBz;
                        lgxz = -pgxz; lEy = -pEy;
                        lBy = -pBy; lgxy = -pgxy;
                        break;
                    case 4:
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpx = -px; lEx = -pEx;
                        lBy = -pBy; lBz = -pBz;
                        lgxy = -pgxy; lgxz = -pgxz;
                        break;
                    case 5:
                        costheta = -costheta_base;
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpz = -pz; lpx = -px;
                        lEz = -pEz; lBz = -pBz;
                        lgyz = -pgyz; lEx = -pEx;
                        lBx = -pBx; lgxy = -pgxy;
                        break;
                    case 6:
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpx = -px; lpy = -py;
                        lEx = -pEx; lBx = -pBx;
                        lgxz = -pgxz; lEy = -pEy;
                        lBy = -pBy; lgyz = -pgyz;
                        break;
                    case 7:
                        costheta = -costheta_base;
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lpx = -px; lpy = -py; lpz = -pz;
                        lEx = -pEx; lEy = -pEy; lEz = -pEz;
                        break;
                }

                double psi1RR, psi1II;
                getnpem1_point_device(lpx, lpy, lpz, pchi, lgxx, lgxy, lgxz, lgyy, lgyz, lgzz,
                                      lEx, lEy, lEz, lBx, lBy, lBz, psi1RR, psi1II);

                double thetap = sqrt((2.0 * pl + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(pl, pm, spinw, costheta);
                double chi_scale = 1.0 / ((pchi + 1.0) * (pchi + 1.0));

#ifdef GaussInt
                local_RP += thetap * chi_scale * (psi1RR * cosmphi + psi1II * sinmphi) * wt;
                local_IP += thetap * chi_scale * (psi1II * cosmphi - psi1RR * sinmphi) * wt;
#else
                local_RP += thetap * chi_scale * (psi1RR * cosmphi + psi1II * sinmphi);
                local_IP += thetap * chi_scale * (psi1II * cosmphi - psi1RR * sinmphi);
#endif
            }

            atomicAdd(&d_RP_out[countlm], local_RP);
            atomicAdd(&d_IP_out[countlm], local_IP);
            countlm++;
        }
    }
}

static void gpu_surf_Wave_EM1_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN, int InList,
    const double* d_shellf,
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int num_elements = Nmax - Nmin + 1;
    if (num_elements <= 0) return;
    int blockSize = 256;
    int gridSize = (num_elements + blockSize - 1) / blockSize;
    surf_Wave_EM1_kernel<<<gridSize, blockSize, 0, stream>>>(
        Nmin, Nmax, N_phi, Symmetry, spinw, maxl, NN, InList,
        d_shellf, d_pox_0, d_pox_1, d_pox_2,
        d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
    );
}

static void gpu_surf_Wave_EM_launch(
    cudaStream_t stream,
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN, int InList,
    const double* d_shellf,
    const double* d_pox_0, const double* d_pox_1, const double* d_pox_2,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int num_elements = Nmax - Nmin + 1;
    if (num_elements <= 0) return;
    int blockSize = 256;
    int gridSize = (num_elements + blockSize - 1) / blockSize;
    surf_Wave_EM_kernel<<<gridSize, blockSize, 0, stream>>>(
        Nmin, Nmax, N_phi, Symmetry, spinw, maxl, NN, InList,
        d_shellf, d_pox_0, d_pox_1, d_pox_2,
        d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
    );
}

static void gpu_em_reduce_and_scale(
    cudaStream_t stream, int NN,
    double* d_RP_out, double* d_IP_out,
    double* RP, double* IP, double rex,
    double dphi_val, double dcostheta_val
) {
#if MPI_CUDA_AWARE
    GPUManager::synchronize_stream(stream);
    MPI_Allreduce(MPI_IN_PLACE, d_RP_out, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, d_IP_out, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, stream);
    GPUManager::synchronize_stream(stream);
#else
    double *h_RP_local = new double[NN];
    double *h_IP_local = new double[NN];
    GPUManager::getInstance().sync_to_cpu(h_RP_local, d_RP_out, NN, stream);
    GPUManager::getInstance().sync_to_cpu(h_IP_local, d_IP_out, NN, stream);
    GPUManager::synchronize_stream(stream);
    MPI_Allreduce(h_RP_local, RP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(h_IP_local, IP, NN, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    delete[] h_RP_local;
    delete[] h_IP_local;
#endif

#ifdef GaussInt
    double scale = rex * dphi_val;
#else
    double scale = rex * dphi_val * dcostheta_val;
#endif
    for (int ii = 0; ii < NN; ii++) {
        RP[ii] *= scale;
        IP[ii] *= scale;
    }
}

// Tier 2.3: gpu_surf_Wave(ShellPatch*, Ex..Bz, chi, gxx..gzz)
void surface_integral::gpu_surf_Wave(
    double rex, int lev, ShellPatch *GH,
    var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
    var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    const int InList = 13;

    MyList<var> *DG_List = new MyList<var>(Ex);
    DG_List->insert(Ey); DG_List->insert(Ez);
    DG_List->insert(Bx); DG_List->insert(By); DG_List->insert(Bz);
    DG_List->insert(chi);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);

    double *d_pox[3];
    for (int i = 0; i < 3; i++)
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;
    int Nmin, Nmax;
    if (Lp > myrank) { Nmin = myrank * mp + myrank; Nmax = Nmin + mp; }
    else { Nmin = myrank * mp + Lp; Nmax = Nmin + mp - 1; }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    gpu_surf_Wave_EM_launch(
        this->stream, Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN, InList,
        d_shellf, d_pox[0], d_pox[1], d_pox[2],
        d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
    );

    gpu_em_reduce_and_scale(this->stream, NN, d_RP_out, d_IP_out, RP, IP, rex, dphi, dcostheta);

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    DG_List->clearList();
}

// Tier 2.4: gpu_surf_Wave(ShellPatch*, Ex..Bz, chi, gxx..gzz, funcs)
void surface_integral::gpu_surf_Wave(
    double rex, int lev, ShellPatch *GH,
    var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
    var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor,
    void (*funcs)(double &, double &, double &,
                  double &, double &, double &, double &, double &, double &, double &,
                  double &, double &, double &, double &, double &, double &,
                  double &, double &)
) {
    const int InList = 13;

    MyList<var> *DG_List = new MyList<var>(Ex);
    DG_List->insert(Ey); DG_List->insert(Ez);
    DG_List->insert(Bx); DG_List->insert(By); DG_List->insert(Bz);
    DG_List->insert(chi);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);

    double *d_pox[3];
    for (int i = 0; i < 3; i++)
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;
    int Nmin, Nmax;
    if (Lp > myrank) { Nmin = myrank * mp + myrank; Nmax = Nmin + mp; }
    else { Nmin = myrank * mp + Lp; Nmax = Nmin + mp - 1; }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    if (funcs == f_getnpem1_point) {
        gpu_surf_Wave_EM1_launch(
            this->stream, Nmin, Nmax, N_phi, Symmetry,
            spinw, maxl, NN, InList,
            d_shellf, d_pox[0], d_pox[1], d_pox[2],
            d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
        );
    } else {
        gpu_surf_Wave_EM_launch(
            this->stream, Nmin, Nmax, N_phi, Symmetry,
            spinw, maxl, NN, InList,
            d_shellf, d_pox[0], d_pox[1], d_pox[2],
            d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
        );
    }

    gpu_em_reduce_and_scale(this->stream, NN, d_RP_out, d_IP_out, RP, IP, rex, dphi, dcostheta);

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    DG_List->clearList();
}

// Tier 2.5: gpu_surf_Wave(cgh*, Ex..Bz, chi, gxx..gzz, funcs)
void surface_integral::gpu_surf_Wave(
    double rex, int lev, cgh *GH,
    var *Ex, var *Ey, var *Ez, var *Bx, var *By, var *Bz,
    var *chi, var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor,
    void (*funcs)(double &, double &, double &,
                  double &, double &, double &, double &, double &, double &, double &,
                  double &, double &, double &, double &, double &, double &,
                  double &, double &)
) {
    if (myrank == 0 && GH->grids[lev] != 1) {
        if (Monitor && Monitor->outfile)
            Monitor->outfile << "WARNING: surface integral on multipatches" << endl;
        else
            cout << "WARNING: surface integral on multipatches" << endl;
    }

    const int InList = 13;

    MyList<var> *DG_List = new MyList<var>(Ex);
    DG_List->insert(Ey); DG_List->insert(Ez);
    DG_List->insert(Bx); DG_List->insert(By); DG_List->insert(Bz);
    DG_List->insert(chi);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);

    double *d_pox[3];
    for (int i = 0; i < 3; i++)
        d_pox[i] = GPUManager::getInstance().allocate_device_memory<double>(n_tot, this->stream);
    gpu_scale_normals_launch(this->stream, n_tot, rex, d_nx_g, d_ny_g, d_nz_g, d_pox[0], d_pox[1], d_pox[2]);

    double* d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GH->PatL[lev]->data->Interp_Points_GPU(this->stream, DG_List, n_tot, d_pox, d_shellf, Symmetry);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;
    int Nmin, Nmax;
    if (Lp > myrank) { Nmin = myrank * mp + myrank; Nmax = Nmin + mp; }
    else { Nmin = myrank * mp + Lp; Nmax = Nmin + mp - 1; }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    if (funcs == f_getnpem1_point) {
        gpu_surf_Wave_EM1_launch(
            this->stream, Nmin, Nmax, N_phi, Symmetry,
            spinw, maxl, NN, InList,
            d_shellf, d_pox[0], d_pox[1], d_pox[2],
            d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
        );
    } else {
        gpu_surf_Wave_EM_launch(
            this->stream, Nmin, Nmax, N_phi, Symmetry,
            spinw, maxl, NN, InList,
            d_shellf, d_pox[0], d_pox[1], d_pox[2],
            d_wtcostheta, d_arcostheta, dphi, d_RP_out, d_IP_out
        );
    }

    gpu_em_reduce_and_scale(this->stream, NN, d_RP_out, d_IP_out, RP, IP, rex, dphi, dcostheta);

    GPUManager::getInstance().free_device_memory(d_pox[0], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[1], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_pox[2], n_tot, this->stream);
    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    DG_List->clearList();
}

// ============================================================
// Tier 3: Point_Psi4 — CPU interp, GPU Wigner-d decomposition
// ============================================================

#include "getnp4_gpu.cuh"

// Kernel: applies symmetry flips + getnp4_point + Wigner-d decomposition
// InList=62, psi4 at offset 0 (R) and 1 (I) after getnp4_point
__global__ void surf_Wave_PointPsi4_kernel(
    int Nmin, int Nmax, int N_phi, int Symmetry,
    int spinw, int maxl, int NN,
    const int InList,
    const double* d_shellf,
    const double* d_wtcostheta, const double* d_arcostheta,
    double dphi,
    double* d_RP_out, double* d_IP_out
) {
    int n = Nmin + blockIdx.x * blockDim.x + threadIdx.x;
    if (n > Nmax) return;

    int i = n / N_phi;
    int j = n - i * N_phi;

    int lpsy = 0;
    if (Symmetry == 0) lpsy = 1;
    else if (Symmetry == 1) lpsy = 2;
    else if (Symmetry == 2) lpsy = 8;

    double wt = d_wtcostheta[i];
    double costheta_base = d_arcostheta[i];

    int countlm = 0;
    for (int pl = spinw; pl <= maxl; pl++) {
        for (int pm = -pl; pm <= pl; pm++) {
            double phi_base = (j + 0.5) * dphi;

            // Read all 62 vars at once
            double pchi = d_shellf[InList * n];
            double ptrK = d_shellf[InList * n + 1];
            double pgxx = d_shellf[InList * n + 2];
            double pgxy = d_shellf[InList * n + 3];
            double pgxz = d_shellf[InList * n + 4];
            double pgyy = d_shellf[InList * n + 5];
            double pgyz = d_shellf[InList * n + 6];
            double pgzz = d_shellf[InList * n + 7];
            double pAxx = d_shellf[InList * n + 8];
            double pAxy = d_shellf[InList * n + 9];
            double pAxz = d_shellf[InList * n + 10];
            double pAyy = d_shellf[InList * n + 11];
            double pAyz = d_shellf[InList * n + 12];
            double pAzz = d_shellf[InList * n + 13];
            double pchix = d_shellf[InList * n + 14];
            double pchiy = d_shellf[InList * n + 15];
            double pchiz = d_shellf[InList * n + 16];
            double ptrKx = d_shellf[InList * n + 17];
            double ptrKy = d_shellf[InList * n + 18];
            double ptrKz = d_shellf[InList * n + 19];
            double pAxxx = d_shellf[InList * n + 20];
            double pAxxy = d_shellf[InList * n + 21];
            double pAxxz = d_shellf[InList * n + 22];
            double pAxyx = d_shellf[InList * n + 23];
            double pAxyy = d_shellf[InList * n + 24];
            double pAxyz = d_shellf[InList * n + 25];
            double pAxzx = d_shellf[InList * n + 26];
            double pAxzy = d_shellf[InList * n + 27];
            double pAxzz = d_shellf[InList * n + 28];
            double pAyyx = d_shellf[InList * n + 29];
            double pAyyy = d_shellf[InList * n + 30];
            double pAyyz = d_shellf[InList * n + 31];
            double pAyzx = d_shellf[InList * n + 32];
            double pAyzy = d_shellf[InList * n + 33];
            double pAyzz = d_shellf[InList * n + 34];
            double pAzzx = d_shellf[InList * n + 35];
            double pAzzy = d_shellf[InList * n + 36];
            double pAzzz = d_shellf[InList * n + 37];
            double pGamxxx = d_shellf[InList * n + 38];
            double pGamxxy = d_shellf[InList * n + 39];
            double pGamxxz = d_shellf[InList * n + 40];
            double pGamxyy = d_shellf[InList * n + 41];
            double pGamxyz = d_shellf[InList * n + 42];
            double pGamxzz = d_shellf[InList * n + 43];
            double pGamyxx = d_shellf[InList * n + 44];
            double pGamyxy = d_shellf[InList * n + 45];
            double pGamyxz = d_shellf[InList * n + 46];
            double pGamyyy = d_shellf[InList * n + 47];
            double pGamyyz = d_shellf[InList * n + 48];
            double pGamyzz = d_shellf[InList * n + 49];
            double pGamzxx = d_shellf[InList * n + 50];
            double pGamzxy = d_shellf[InList * n + 51];
            double pGamzxz = d_shellf[InList * n + 52];
            double pGamzyy = d_shellf[InList * n + 53];
            double pGamzyz = d_shellf[InList * n + 54];
            double pGamzzz = d_shellf[InList * n + 55];
            double pRxx = d_shellf[InList * n + 56];
            double pRxy = d_shellf[InList * n + 57];
            double pRxz = d_shellf[InList * n + 58];
            double pRyy = d_shellf[InList * n + 59];
            double pRyz = d_shellf[InList * n + 60];
            double pRzz = d_shellf[InList * n + 61];

            double local_RP = 0.0;
            double local_IP = 0.0;

            for (int lp = 0; lp < lpsy; lp++) {
                double costheta = costheta_base;
                double phi_angle = phi_base;
                double cosmphi, sinmphi;

                double lchi = pchi, ltrK = ptrK;
                double lgxx = pgxx, lgxy = pgxy, lgxz = pgxz;
                double lgyy = pgyy, lgyz = pgyz, lgzz = pgzz;
                double lAxx = pAxx, lAxy = pAxy, lAxz = pAxz;
                double lAyy = pAyy, lAyz = pAyz, lAzz = pAzz;
                double lchix = pchix, lchiy = pchiy, lchiz = pchiz;
                double ltrKx = ptrKx, ltrKy = ptrKy, ltrKz = ptrKz;
                double lAxxx = pAxxx, lAxxy = pAxxy, lAxxz = pAxxz;
                double lAxyx = pAxyx, lAxyy = pAxyy, lAxyz = pAxyz;
                double lAxzx = pAxzx, lAxzy = pAxzy, lAxzz = pAxzz;
                double lAyyx = pAyyx, lAyyy = pAyyy, lAyyz = pAyyz;
                double lAyzx = pAyzx, lAyzy = pAyzy, lAyzz = pAyzz;
                double lAzzx = pAzzx, lAzzy = pAzzy, lAzzz = pAzzz;
                double lGamxxx = pGamxxx, lGamxxy = pGamxxy, lGamxxz = pGamxxz;
                double lGamxyy = pGamxyy, lGamxyz = pGamxyz, lGamxzz = pGamxzz;
                double lGamyxx = pGamyxx, lGamyxy = pGamyxy, lGamyxz = pGamyxz;
                double lGamyyy = pGamyyy, lGamyyz = pGamyyz, lGamyzz = pGamyzz;
                double lGamzxx = pGamzxx, lGamzxy = pGamzxy, lGamzxz = pGamzxz;
                double lGamzyy = pGamzyy, lGamzyz = pGamzyz, lGamzzz = pGamzzz;
                double lRxx = pRxx, lRxy = pRxy, lRxz = pRxz;
                double lRyy = pRyy, lRyz = pRyz, lRzz = pRzz;

                // NOTE: px,py,pz are not flipped in this variant (see CPU code)
                switch (lp) {
                    case 0:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        break;
                    case 1:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lgxz = -lgxz; lgyz = -lgyz;
                        lAxz = -lAxz; lAyz = -lAyz;
                        lchiz = -lchiz; ltrKz = -ltrKz;
                        lAxxz = -lAxxz; lAxyz = -lAxyz;
                        lAxzx = -lAxzx; lAxzy = -lAxzy;
                        lAyyz = -lAyyz; lAyzx = -lAyzx;
                        lAyzy = -lAyzy; lAzzz = -lAzzz;
                        lGamxxz = -lGamxxz; lGamxyz = -lGamxyz;
                        lGamyxz = -lGamyxz; lGamyyz = -lGamyyz;
                        lGamzxx = -lGamzxx; lGamzxy = -lGamzxy;
                        lGamzyy = -lGamzyy; lGamzzz = -lGamzzz;
                        lRxz = -lRxz; lRyz = -lRyz;
                        break;
                    case 2:
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        lgxy = -lgxy; lgyz = -lgyz;
                        lAxy = -lAxy; lAyz = -lAyz;
                        lchiy = -lchiy; ltrKy = -ltrKy;
                        lAxxy = -lAxxy; lAxyx = -lAxyx;
                        lAxyz = -lAxyz; lAxzy = -lAxzy;
                        lAyyy = -lAyyy; lAyzx = -lAyzx;
                        lAyzz = -lAyzz; lAzzy = -lAzzy;
                        lGamxxy = -lGamxxy; lGamxyz = -lGamxyz;
                        lGamyxx = -lGamyxx; lGamyxz = -lGamyxz;
                        lGamyyy = -lGamyyy; lGamyyz = -lGamyyz;
                        lGamzxy = -lGamzxy; lGamzyz = -lGamzyz;
                        lRxy = -lRxy; lRyz = -lRyz;
                        break;
                    case 3:
                        costheta = -costheta_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = -sin(pm * phi_angle);
                        lgxy = -lgxy; lgxz = -lgxz;
                        lAxy = -lAxy; lAxz = -lAxz;
                        lchiy = -lchiy; lchiz = -lchiz;
                        ltrKy = -ltrKy; ltrKz = -ltrKz;
                        lAxxy = -lAxxy; lAxxz = -lAxxz;
                        lAxyx = -lAxyx; lAxzx = -lAxzx;
                        lAyyy = -lAyyy; lAyyz = -lAyyz;
                        lAyzy = -lAyzy; lAyzz = -lAyzz;
                        lAzzy = -lAzzy; lAzzz = -lAzzz;
                        lGamxxy = -lGamxxy; lGamxxz = -lGamxxz;
                        lGamyxx = -lGamyxx; lGamyyy = -lGamyyy;
                        lGamyyz = -lGamyyz; lGamyzz = -lGamyzz;
                        lGamzxx = -lGamzxx; lGamzyy = -lGamzyy;
                        lGamzyz = -lGamzyz; lGamzzz = -lGamzzz;
                        lRxy = -lRxy; lRxz = -lRxz;
                        break;
                    case 4:
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lgxy = -lgxy; lgxz = -lgxz;
                        lAxy = -lAxy; lAxz = -lAxz;
                        lchix = -lchix; ltrKx = -ltrKx;
                        lAxxx = -lAxxx; lAxyy = -lAxyy;
                        lAxyz = -lAxyz; lAxzy = -lAxzy;
                        lAxzz = -lAxzz; lAyyx = -lAyyx;
                        lAyzx = -lAyzx; lAzzx = -lAzzx;
                        lGamxxx = -lGamxxx; lGamxyy = -lGamxyy;
                        lGamxyz = -lGamxyz; lGamxzz = -lGamxzz;
                        lGamyxy = -lGamyxy; lGamyxz = -lGamyxz;
                        lGamzxy = -lGamzxy; lGamzxz = -lGamzxz;
                        lRxy = -lRxy; lRxz = -lRxz;
                        break;
                    case 5:
                        costheta = -costheta_base;
                        phi_angle = M_PI - phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lgxy = -lgxy; lgyz = -lgyz;
                        lAxy = -lAxy; lAyz = -lAyz;
                        lchix = -lchix; lchiz = -lchiz;
                        ltrKx = -ltrKx; ltrKz = -ltrKz;
                        lAxxx = -lAxxx; lAxxz = -lAxxz;
                        lAxyy = -lAxyy; lAxzx = -lAxzx;
                        lAxzz = -lAxzz; lAyyx = -lAyyx;
                        lAyyz = -lAyyz; lAyzy = -lAyzy;
                        lAzzx = -lAzzx; lAzzz = -lAzzz;
                        lGamxxx = -lGamxxx; lGamxxz = -lGamxxz;
                        lGamxyy = -lGamxyy; lGamxzz = -lGamxzz;
                        lGamyxy = -lGamyxy; lGamyyz = -lGamyyz;
                        lGamzxx = -lGamzxx; lGamzxz = -lGamzxz;
                        lGamzyy = -lGamzyy; lGamzzz = -lGamzzz;
                        lRxy = -lRxy; lRyz = -lRyz;
                        break;
                    case 6:
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lgxz = -lgxz; lgyz = -lgyz;
                        lAxz = -lAxz; lAyz = -lAyz;
                        lchix = -lchix; lchiy = -lchiy;
                        ltrKx = -ltrKx; ltrKy = -ltrKy;
                        lAxxx = -lAxxx; lAxxy = -lAxxy;
                        lAxyx = -lAxyx; lAxyy = -lAxyy;
                        lAxzz = -lAxzz; lAyyx = -lAyyx;
                        lAyyy = -lAyyy; lAyzz = -lAyzz;
                        lAzzx = -lAzzx; lAzzy = -lAzzy;
                        lGamxxx = -lGamxxx; lGamxxy = -lGamxxy;
                        lGamxyy = -lGamxyy; lGamxzz = -lGamxzz;
                        lGamyxx = -lGamyxx; lGamyxy = -lGamyxy;
                        lGamyyy = -lGamyyy; lGamyzz = -lGamyzz;
                        lGamzxz = -lGamzxz; lGamzyz = -lGamzyz;
                        lRxz = -lRxz; lRyz = -lRyz;
                        break;
                    case 7:
                        costheta = -costheta_base;
                        phi_angle = M_PI + phi_base;
                        cosmphi = cos(pm * phi_angle);
                        sinmphi = sin(pm * phi_angle);
                        lchix = -lchix; lchiy = -lchiy; lchiz = -lchiz;
                        ltrKx = -ltrKx; ltrKy = -ltrKy; ltrKz = -ltrKz;
                        lAxxx = -lAxxx; lAxxy = -lAxxy; lAxxz = -lAxxz;
                        lAxyx = -lAxyx; lAxyy = -lAxyy; lAxyz = -lAxyz;
                        lAxzx = -lAxzx; lAxzy = -lAxzy; lAxzz = -lAxzz;
                        lAyyx = -lAyyx; lAyyy = -lAyyy; lAyyz = -lAyyz;
                        lAyzx = -lAyzx; lAyzy = -lAyzy; lAyzz = -lAyzz;
                        lAzzx = -lAzzx; lAzzy = -lAzzy; lAzzz = -lAzzz;
                        lGamxxx = -lGamxxx; lGamxxy = -lGamxxy;
                        lGamxxz = -lGamxxz; lGamxyy = -lGamxyy;
                        lGamxyz = -lGamxyz; lGamxzz = -lGamxzz;
                        lGamyxx = -lGamyxx; lGamyxy = -lGamyxy;
                        lGamyxz = -lGamyxz; lGamyyy = -lGamyyy;
                        lGamyyz = -lGamyyz; lGamyzz = -lGamyzz;
                        lGamzxx = -lGamzxx; lGamzxy = -lGamzxy;
                        lGamzxz = -lGamzxz; lGamzyy = -lGamzyy;
                        lGamzyz = -lGamzyz; lGamzzz = -lGamzzz;
                        break;
                }

                double psi4RR, psi4II;
                d_getnp4_point<2>(0, 0, 0, lchi, ltrK,
                    lgxx, lgxy, lgxz, lgyy, lgyz, lgzz,
                    lAxx, lAxy, lAxz, lAyy, lAyz, lAzz,
                    lchix, lchiy, lchiz, ltrKx, ltrKy, ltrKz,
                    lAxxx, lAxxy, lAxxz, lAxyx, lAxyy, lAxyz,
                    lAxzx, lAxzy, lAxzz, lAyyx, lAyyy, lAyyz,
                    lAyzx, lAyzy, lAyzz, lAzzx, lAzzy, lAzzz,
                    lGamxxx, lGamxxy, lGamxxz, lGamxyy, lGamxyz, lGamxzz,
                    lGamyxx, lGamyxy, lGamyxz, lGamyyy, lGamyyz, lGamyzz,
                    lGamzxx, lGamzxy, lGamzxz, lGamzyy, lGamzyz, lGamzzz,
                    lRxx, lRxy, lRxz, lRyy, lRyz, lRzz,
                    psi4RR, psi4II);

                double thetap = sqrt((2.0 * pl + 1.0) / (4.0 * M_PI)) * misc::wigner_d_device(pl, pm, spinw, costheta);
                double chi_scale = 1.0 / ((lchi + 1.0) * (lchi + 1.0));

#ifdef GaussInt
                local_RP += thetap * chi_scale * (psi4RR * cosmphi + psi4II * sinmphi) * wt;
                local_IP += thetap * chi_scale * (psi4II * cosmphi - psi4RR * sinmphi) * wt;
#else
                local_RP += thetap * chi_scale * (psi4RR * cosmphi + psi4II * sinmphi);
                local_IP += thetap * chi_scale * (psi4II * cosmphi - psi4RR * sinmphi);
#endif
            }

            atomicAdd(&d_RP_out[countlm], local_RP);
            atomicAdd(&d_IP_out[countlm], local_IP);
            countlm++;
        }
    }
}

void surface_integral::gpu_surf_Wave(
    double rex, cgh *GH, ShellPatch *SH,
    var *chi, var *trK,
    var *gxx, var *gxy, var *gxz, var *gyy, var *gyz, var *gzz,
    var *Axx, var *Axy, var *Axz, var *Ayy, var *Ayz, var *Azz,
    var *chix, var *chiy, var *chiz,
    var *trKx, var *trKy, var *trKz,
    var *Axxx, var *Axxy, var *Axxz,
    var *Axyx, var *Axyy, var *Axyz,
    var *Axzx, var *Axzy, var *Axzz,
    var *Ayyx, var *Ayyy, var *Ayyz,
    var *Ayzx, var *Ayzy, var *Ayzz,
    var *Azzx, var *Azzy, var *Azzz,
    var *Gamxxx, var *Gamxxy, var *Gamxxz, var *Gamxyy, var *Gamxyz, var *Gamxzz,
    var *Gamyxx, var *Gamyxy, var *Gamyxz, var *Gamyyy, var *Gamyyz, var *Gamyzz,
    var *Gamzxx, var *Gamzxy, var *Gamzxz, var *Gamzyy, var *Gamzyz, var *Gamzzz,
    var *Rxx, var *Rxy, var *Rxz, var *Ryy, var *Ryz, var *Rzz,
    int spinw, int maxl, int NN, double *RP, double *IP,
    monitor *Monitor
) {
    const int InList = 62;

    MyList<var> *DG_List = new MyList<var>(chi);
    DG_List->insert(trK);
    DG_List->insert(gxx); DG_List->insert(gxy); DG_List->insert(gxz);
    DG_List->insert(gyy); DG_List->insert(gyz); DG_List->insert(gzz);
    DG_List->insert(Axx); DG_List->insert(Axy); DG_List->insert(Axz);
    DG_List->insert(Ayy); DG_List->insert(Ayz); DG_List->insert(Azz);
    DG_List->insert(chix); DG_List->insert(chiy); DG_List->insert(chiz);
    DG_List->insert(trKx); DG_List->insert(trKy); DG_List->insert(trKz);
    DG_List->insert(Axxx); DG_List->insert(Axxy); DG_List->insert(Axxz);
    DG_List->insert(Axyx); DG_List->insert(Axyy); DG_List->insert(Axyz);
    DG_List->insert(Axzx); DG_List->insert(Axzy); DG_List->insert(Axzz);
    DG_List->insert(Ayyx); DG_List->insert(Ayyy); DG_List->insert(Ayyz);
    DG_List->insert(Ayzx); DG_List->insert(Ayzy); DG_List->insert(Ayzz);
    DG_List->insert(Azzx); DG_List->insert(Azzy); DG_List->insert(Azzz);
    DG_List->insert(Gamxxx); DG_List->insert(Gamxxy); DG_List->insert(Gamxxz);
    DG_List->insert(Gamxyy); DG_List->insert(Gamxyz); DG_List->insert(Gamxzz);
    DG_List->insert(Gamyxx); DG_List->insert(Gamyxy); DG_List->insert(Gamyxz);
    DG_List->insert(Gamyyy); DG_List->insert(Gamyyz); DG_List->insert(Gamyzz);
    DG_List->insert(Gamzxx); DG_List->insert(Gamzxy); DG_List->insert(Gamzxz);
    DG_List->insert(Gamzyy); DG_List->insert(Gamzyz); DG_List->insert(Gamzzz);
    DG_List->insert(Rxx); DG_List->insert(Rxy); DG_List->insert(Rxz);
    DG_List->insert(Ryy); DG_List->insert(Ryz); DG_List->insert(Rzz);

    int n;
    double *pox[3];
    for (int i = 0; i < 3; i++)
        pox[i] = new double[n_tot];
    for (n = 0; n < n_tot; n++) {
        pox[0][n] = rex * nx_g[n];
        pox[1][n] = rex * ny_g[n];
        pox[2][n] = rex * nz_g[n];
    }

    double *shellf = new double[n_tot * InList];
    SR_Interp_Points(DG_List, GH, SH, n_tot, pox, shellf);

    // Copy shellf to GPU
    double *d_shellf = GPUManager::getInstance().allocate_device_memory<double>(n_tot * InList, this->stream);
    GPUManager::getInstance().sync_to_gpu(shellf, d_shellf, n_tot * InList, this->stream);

    int mp = n_tot / cpusize;
    int Lp = n_tot - cpusize * mp;
    int Nmin, Nmax;
    if (Lp > myrank) { Nmin = myrank * mp + myrank; Nmax = Nmin + mp; }
    else { Nmin = myrank * mp + Lp; Nmax = Nmin + mp - 1; }

    double *d_RP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    double *d_IP_out = GPUManager::getInstance().allocate_comm_buffer<double>(NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_RP_out, NN, this->stream);
    GPUManager::getInstance().memset_device_memory(d_IP_out, NN, this->stream);

    int blockSize = 256;
    int gridSize = (n_tot + blockSize - 1) / blockSize;
    surf_Wave_PointPsi4_kernel<<<gridSize, blockSize, 0, this->stream>>>(
        Nmin, Nmax, N_phi, Symmetry,
        spinw, maxl, NN, InList,
        d_shellf, d_wtcostheta, d_arcostheta, dphi,
        d_RP_out, d_IP_out
    );

    // Manual gather to root: only rank 0 reads RP/IP (Psi4Monitor->writefile
    // is rank-0-only). Cheaper than MPI_Allreduce/Reduce on small N+small NN.
    GPUManager::getInstance().sync_to_cpu(RP, d_RP_out, NN, this->stream);
    GPUManager::getInstance().sync_to_cpu(IP, d_IP_out, NN, this->stream);
    GPUManager::synchronize_stream(this->stream);
    amss_reduce::gather_sum_to_root(RP, NN, 0, MPI_COMM_WORLD);
    amss_reduce::gather_sum_to_root(IP, NN, 0, MPI_COMM_WORLD);

#ifdef GaussInt
    double scale = rex * dphi;
#else
    double scale = rex * dphi * dcostheta;
#endif
    for (int ii = 0; ii < NN; ii++) {
        RP[ii] *= scale;
        IP[ii] *= scale;
    }

    GPUManager::getInstance().free_device_memory(d_shellf, n_tot * InList, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_RP_out, NN, this->stream);
    GPUManager::getInstance().free_comm_buffer(d_IP_out, NN, this->stream);

    delete[] pox[0]; delete[] pox[1]; delete[] pox[2];
    delete[] shellf;
    DG_List->clearList();
}
