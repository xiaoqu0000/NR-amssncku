
#ifndef SHELLPATCH_H
#define SHELLPATCH_H

#include <mpi.h>
#include "MyList.h"
#include "Block.h"
#include "Parallel.h"
#include "var.h"
#include "monitor.h"
#include "macrodef.h" //need dim here; Vertex or Cell; ghost_width

#if (dim != 3)
#error shellpatch only supports 3 dimensional stuff yet
#endif

struct ShellPatchDevicePoint
{
   double gpox[dim];
   double lpox[dim];
   int block_index;
   int block_rank;
   int ssst;
   int tsst;
   int dumyd;
   int shape[dim];
   double bbox[2 * dim];
   int sind[dim];
   int coef_offset;
   int coef_count;
};

struct ShellPatchDeviceTransferList
{
   int count;
   int coef_count;
   ShellPatchDevicePoint *d_src;
   ShellPatchDevicePoint *d_dst;
   double *d_coef;
};

struct ShellPatchDeviceMirror
{
   int node_count;
   int shell_block_count;
   int cartesian_block_count;
   ShellPatchDeviceTransferList *ss;
   ShellPatchDeviceTransferList *csats;
   ShellPatchDeviceTransferList *csatc;
};

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
   void Sync_GPU(MyList<var> *VarList, int Symmetry);
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
   MyList<Patch> *CPatL_for_gpu;

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
   ShellPatchDeviceMirror device_mirror;
   cudaStream_t device_stream;

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
   void prolongpointstru(MyList<pointstru> *&psul, MyList<pointstru> *&psul_tail,
                         MyList<ss_patch> *sPp, double DH[dim],
                         MyList<Patch> *Pp, double CDH[dim], MyList<pointstru> *pss);
   bool prolongpointstru(MyList<pointstru> *&psul, bool ssyn, int tsst, MyList<ss_patch> *sPp, double DH[dim],
                         MyList<Patch> *Pp, double CDH[dim], double x, double y, double z, int Symmetry, int rank_in);
   void setupintintstuff(int cpusize, MyList<Patch> *CPatL, int Symmetry);
   void intertransfer(MyList<pointstru> **src, MyList<pointstru> **dst,
                      MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                      int Symmetry);
   void intertransfer_gpu(MyList<pointstru> **src, MyList<pointstru> **dst,
                          MyList<var> *VarList1 /* source */, MyList<var> *VarList2 /*target */,
                          int Symmetry);
   int interdata_packer(double *data, MyList<pointstru> *src, MyList<pointstru> *dst,
                        int rank_in, int dir,
                        MyList<var> *VarLists /* source */, MyList<var> *VarListd /* target */,
                        int Symmetry);
   void setup_device_mirror(int cpusize);
   void release_device_mirror();
   void move_to_gpu(MyList<var> *VarList);
   void move_to_cpu(MyList<var> *VarList);
   int device_transfer_count(const char *family, int rank_in) const;
   int local_shell_block_count() const;
   void Synch(MyList<var> *VarList, int Symmetry);
   void CS_Inter(MyList<var> *VarList, int Symmetry);
   void Synch_GPU(MyList<var> *VarList, int Symmetry);
   void CS_Inter_GPU(MyList<var> *VarList, int Symmetry);
   void destroypsuList(MyList<pointstru> *ct);
   int getdumydimension(int acsst, int posst); // -1 means no dumy dimension
   void matchcheck(MyList<Patch> *CPatL);
   void shellname(char *sn, int i);
   void Interp_Points(MyList<var> *VarList,
                      int NN, double **XX, /*input global Cartesian coordinate*/
                      double *Shellf, int Symmetry);
   bool Interp_One_Point(MyList<var> *VarList,
                         double *XX, /*input global Cartesian coordinate*/
                         double *Shellf, int Symmetry);
   void write_Pablo_file_ss(int *ext, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax,
                            char *filename, int sst);
   double L2Norm(var *vf);
   void Find_Maximum(MyList<var> *VarList, double *XX, double *Shellf);
};

#endif /* SHELLPATCH_H */
