
#ifndef SHELLPATCH_H
#define SHELLPATCH_H

#include <mpi.h>
#include <cuda_runtime.h>
#include <vector>
#include "MyList.h"
#include "Block.h"
#include "Parallel.h"
#include "var.h"
#include "monitor.h"
#include "macrodef.h" //need dim here; Vertex or Cell; ghost_width

#if (dim != 3)
#error shellpatch only supports 3 dimensional stuff yet
#endif

class ss_patch
{

public:
   int sst; // ss_patch type: 0:zp, 1:zm, 2:xp, 3:xm, 4:yp, 5:ym
   int myrank;
   int shape[dim];
   double bbox[2 * dim]; // this bbox includes nominal points and overlap points
   MyList<Block> *blb, *ble;
   int ingfs, fngfs;

   ss_patch() {};
   ss_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki);

   ~ss_patch();

   virtual void setupcordtrans() {};
   void Sync(MyList<var> *VarList, int Symmetry);
#ifdef USE_GPU
   void Sync_GPU(MyList<var> *VarList, int Symmetry);
#endif
   MyList<Parallel::gridseg> *build_bulk_gsl(Block *bp);
   MyList<Parallel::gridseg> *build_ghost_gsl();
   MyList<Parallel::gridseg> *build_owned_gsl0(int rank_in);
};

class xp_patch : public ss_patch
{
public:
   xp_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 2; };
   void setupcordtrans();
};

class xm_patch : public ss_patch
{
public:
   xm_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 3; };
   void setupcordtrans();
};
class yp_patch : public ss_patch
{
public:
   yp_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 4; };
   void setupcordtrans();
};

class ym_patch : public ss_patch
{
public:
   ym_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 5; };
   void setupcordtrans();
};
class zp_patch : public ss_patch
{
public:
   zp_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 0; };
   void setupcordtrans();
};

class zm_patch : public ss_patch
{
public:
   zm_patch(int ingfsi, int fngfsi, int *shapei, double *bboxi, int myranki) : ss_patch(ingfsi, fngfsi, shapei, bboxi, myranki) { sst = 1; };
   void setupcordtrans();
};
// Shell Patch system
// for derivatives usage we ask 27 more double type grid functions
// here we use **sngfs corresponding to fngfs to store them:
//                    drho/dx, drho/dy, drho/dz
//                    dsigma/dx, dsigma/dy, dsigma/dz
//                    dR/dx, dR/dy, dR/dz
//                    drho/dxdx, drho/dxdy, drho/dxdz, drho/dydy, drho/dydz, drho/dzdz
//                    dsigma/dxdx, dsigma/dxdy, dsigma/dxdz, dsigma/dydy, dsigma/dydz, dsigma/dzdz
//                    dR/dxdx, dR/dxdy, dR/dxdz, dR/dydy, dR/dydz, dR/dzdz
class ShellPatch
{

public:
   struct pointstru
   {
      double gpox[dim]; // global cordinate
      double lpox[dim]; // local cordinate
      Block *Bg;
      int ssst; //-1: cardisian, others as sst of ss_patch source sst
      int tsst; //-1: cardisian, others as sst of ss_patch target sst
      double *coef;
      int *sind;
      int dumyd; // the dimension which has common lines, only useful in interdata_packer
                 //-1: means no dumy dimension at all; 0: means rho; 1: means sigma
   };

   int myrank;
   int shape[dim];   // for (rho, sigma, R), for rho and sigma means number of points for every pi/2
   double Rrange[2]; // for Rmin and Rmax
   int Symmetry;
   int ingfs, fngfs;

   MyList<ss_patch> *PatL;

   // we use fngfs+v to reference the variable
   enum
   {
      gx = 0,
      gy,
      gz,
      drhodx,
      drhody,
      drhodz,
      dsigmadx,
      dsigmady,
      dsigmadz,
      dRdx,
      dRdy,
      dRdz,
      drhodxx,
      drhodxy,
      drhodxz,
      drhodyy,
      drhodyz,
      drhodzz,
      dsigmadxx,
      dsigmadxy,
      dsigmadxz,
      dsigmadyy,
      dsigmadyz,
      dsigmadzz,
      dRdxx,
      dRdxy,
      dRdxz,
      dRdyy,
      dRdyz,
      dRdzz
   };

   MyList<pointstru> **ss_src, **ss_dst;
   // at means target
   MyList<pointstru> **csatc_src, **csatc_dst;
   MyList<pointstru> **csats_src, **csats_dst;

#ifdef USE_GPU
   // GPU block metadata table — maps block index to GPU-accessible metadata.
   // Stored as flat arrays for direct GPU kernel access.
   // Built once during init, shared by all pointstru caches.
   int gpu_num_blocks;
   std::vector<Block *> gpu_blocks;
   std::vector<unsigned long long> h_gpu_block_fgfs;
   int    *d_gpu_block_shape;      // [num_blocks * 3]
   int    *d_gpu_block_shape_pad;  // [num_blocks * 3]
   double *d_gpu_block_bbox;       // [num_blocks * 6]
   unsigned long long *d_gpu_block_fgfs; // [num_blocks * fngfs], device ptrs as uint64

   void build_gpu_block_table(MyList<Patch> *CPatL);
   void free_gpu_block_table();

   // Persistent GPU cache for pointstru interpolation data.
   // Built once during init, used by interdata_packer_GPU every RK4 substep.
   // Eliminates per-point H2D transfers and enables batched kernels.
   struct gpu_ps_cache {
      int num_points;       // total points across all ranks
      int max_coef;         // max coef entries per point (ordn * 3 for 3D)
      int ordn;

      double *d_coef;       // [num_points * max_coef], padded, row-major per point
      int    *d_sind;       // [num_points * 3], 0-based starting indices
      int    *d_dumyd;      // [num_points], -1=3D, 0=rho, 1=sigma
      int    *d_ssst;       // [num_points], shell patch type or -1 for Cartesian

      // For UNPACK: destination point local coordinates
      double *d_dst_lpox;   // [num_points * 3]

      // Block indices into gpu_block_table (src for PACK, dst for UNPACK)
      int    *d_src_block_idx; // [num_points]
      int    *d_dst_block_idx; // [num_points]

      // Per-rank offset arrays: point range [offsets[r], offsets[r+1]) belongs to rank r
      int    *d_rank_offsets; // [cpusize+1], GPU copy
      int    *h_rank_offsets; // [cpusize+1], host copy (for index calculation in interdata_packer_GPU)

      gpu_ps_cache() : num_points(0), max_coef(0), ordn(0),
                       d_coef(nullptr), d_sind(nullptr), d_dumyd(nullptr), d_ssst(nullptr),
                       d_dst_lpox(nullptr), d_src_block_idx(nullptr), d_dst_block_idx(nullptr),
                       d_rank_offsets(nullptr), h_rank_offsets(nullptr) {}
   };

   gpu_ps_cache gpu_ss_cache;    // ss_src / ss_dst
   gpu_ps_cache gpu_csatc_cache; // csatc_src / csatc_dst
   gpu_ps_cache gpu_csats_cache; // csats_src / csats_dst

   void build_gpu_pointstru_cache(gpu_ps_cache &cache,
                                  MyList<pointstru> **src_lists,
                                  MyList<pointstru> **dst_lists,
                                  int cpusize,
                                  MyList<Patch> *CPatL);
   void free_gpu_pointstru_cache(gpu_ps_cache &cache);
#endif

public:
   ShellPatch(int ingfsi, int fngfsi, char *filename, int Symmetry, int myranki, monitor *ErrorMonitor);

   ~ShellPatch();

   MyList<Block> *compose_sh(int cpusize, int nodes = 0);
   MyList<Block> *compose_shr(int cpusize, int nodes = 0);
   void setupcordtrans();
   double getR(double r);
   double getsr(double R);
   void checkPatch();
   void checkBlock(int sst);
   void check_pointstrul(MyList<pointstru> *pp, bool first_only);
   void check_pointstrul2(MyList<pointstru> *pp, int first_last_only);
   double getdX(int dir); //(rho, sigma, R)
   void Dump_xyz(char *tag, double time, double dT);
   void Dump_Data(MyList<var> *DumpList, char *tag, double time, double dT);
   double *Collect_Data(ss_patch *PP, var *VP);
   void getlocalpoxsst(double gx, double gy, double gz, int sst, double &lx, double &ly, double &lz);
   void getlocalpox(double gx, double gy, double gz, int &sst, double &lx, double &ly, double &lz);
   void getglobalpox(double &x, double &y, double &z, int sst, double lx, double ly, double lz);
   void prolongpointstru(MyList<pointstru> *&psul, MyList<ss_patch> *sPp, double DH[dim],
                         MyList<Patch> *Pp, double CDH[dim], MyList<pointstru> *pss);
   bool prolongpointstru(MyList<pointstru> *&psul, bool ssyn, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                         MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in);
   void setupintintstuff(int cpusize, MyList<Patch> *CPatL, int Symmetry);
   void intertransfer(MyList<pointstru> **src, MyList<pointstru> **dst,
                      MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                      int Symmetry);
   int interdata_packer(double *data, MyList<pointstru> *src, MyList<pointstru> *dst,
                        int rank_in, int dir,
                        MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */,
                        int Symmetry);
   void Synch(MyList<var> *VarList, int Symmetry);
   void CS_Inter(MyList<var> *VarList, int Symmetry);
#ifdef USE_GPU
   void intertransfer_GPU(MyList<pointstru> **src, MyList<pointstru> **dst,
                           MyList<var> *VarList1, MyList<var> *VarList2,
                           int Symmetry);
   int interdata_packer_GPU(gpu_ps_cache &cache,
                            double *d_data, int buf_size,
                            MyList<pointstru> *src, MyList<pointstru> *dst,
                            int list_rank, int rank_in, int dir,
                            MyList<var> *VarLists, MyList<var> *VarListd,
                            int Symmetry);
   void Synch_GPU(MyList<var> *VarList, int Symmetry);
   void CS_Inter_GPU(MyList<var> *VarList, int Symmetry);
   // GPU-accelerated Phase 2: Cartesian boundary point processing
   void gpu_setup_phase2(int cpusize, MyList<Patch> *CPatL, int Symmetry,
                         double Rrange[2], double DH_radial);
#endif
   void destroypsuList(MyList<pointstru> *ct);
   int getdumydimension(int acsst, int posst); // -1 means no dumy dimension
   void matchcheck(MyList<Patch> *CPatL);
   void shellname(char *sn, int i);
   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX, /*input global Cartesian coordinate*/
                      double *Shellf, int Symmetry);
   void Interp_Points_GPU(
       cudaStream_t stream,
       MyList<var> *VarList,
       int NN, double **d_XX,
       double *d_Shellf, int Symmetry
   );
   bool Interp_One_Point(MyList<var> *VarList,
                         double *XX, /*input global Cartesian coordinate*/
                         double *Shellf, int Symmetry);
   void write_Pablo_file_ss(int *ext, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                            char *filename, int sst);
   double L2Norm(var *vf);
   void Find_Maximum(MyList<var> *VarList, double *XX, double *Shellf);
};

#endif /* SHELLPATCH_H */
