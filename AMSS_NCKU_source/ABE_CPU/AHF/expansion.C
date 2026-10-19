

#include "macrodef.h"
#ifdef With_AHF

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include "util_Table.h"
#include "cctk.h"

#include "config.h"
#include "stdc.h"
#include "myglobal.h"
#include "util.h"
#include "array.h"
#include "cpm_map.h"
#include "linear_map.h"

#include "coords.h"
#include "tgrid.h"
#include "fd_grid.h"
#include "patch.h"
#include "patch_edge.h"
#include "patch_interp.h"
#include "ghost_zone.h"
#include "patch_system.h"

#include "Jacobian.h"

#include "gfns.h"
#include "gr.h"

// all the code in this file is inside this namespace
namespace AHFinderDirect
{
      using jtutil::error_exit;
      using jtutil::pow2;
      using jtutil::pow4;

      namespace
      {

            void setup_xyz_posns(patch_system &ps, bool print_msg_flag);
            enum expansion_status
            interpolate_geometry(patch_system *ps_ptr,
                                 bool initial_flag,
                                 bool print_msg_flag);
            void convert_conformal_to_physical(patch_system &ps,
                                               bool print_msg_flag);

            bool h_is_finite(patch_system &ps, bool initial_flag,
                             bool print_msg_flag);
            bool geometry_is_finite(patch_system &ps, bool initial_flag,
                                    bool print_msg_flag);

            bool compute_Theta(patch_system &ps, fp add_to_expansion,
                               bool Jacobian_flag, jtutil::norm<fp> *Theta_norms_ptr,
                               bool initial_flag,
                               bool print_msg_flag);
      }

      extern struct state state;
      //******************************************************************************
      enum expansion_status
      expansion(patch_system *ps_ptr, fp add_to_expansion,
                bool initial_flag,
                bool Jacobian_flag /* = false */,
                jtutil::norm<fp> *Theta_norms_ptr /* = NULL */)
      {
            const bool active_flag = (ps_ptr != NULL);

            if (active_flag)
                  then
                  {
                        //
                        // normal computation
                        //

                        // fill in values of all ghosted gridfns in ghost zones
                        ps_ptr->synchronize();

                        if (!h_is_finite(*ps_ptr, initial_flag, false))
                              then return expansion_failure__surface_nonfinite;

                        // set up xyz positions of grid points
                        setup_xyz_posns(*ps_ptr, false);
                  }

            {
                  // this is the only function we call unconditionally; it looks at
                  // ps_ptr (non-NULL vs NULL) to choose a normal vs dummy computation
                  const enum expansion_status status = interpolate_geometry(ps_ptr,
                                                                            initial_flag,
                                                                            false);

                  if (status != expansion_success)
                        then return status; // *** ERROR RETURN ***
                  if (active_flag)
                        convert_conformal_to_physical(*ps_ptr, false);
            }

            if (active_flag)
                  then
                  {
                        if (!geometry_is_finite(*ps_ptr, initial_flag, false))
                              then return expansion_failure__geometry_nonfinite;

                        // compute remaining gridfns --> $\Theta$
                        // and optionally also the Jacobian coefficients
                        // by algebraic ops and angular finite differencing
                        if (!compute_Theta(*ps_ptr, add_to_expansion,
                                           Jacobian_flag, Theta_norms_ptr,
                                           initial_flag,
                                           false))
                              then return expansion_failure__gij_not_positive_definite;
                        // *** ERROR RETURN ***
                  }

            return expansion_success; // *** NORMAL RETURN ***
      }

      //******************************************************************************
      namespace
      {
            void setup_xyz_posns(patch_system &ps, bool print_msg_flag)
            {
                  if (print_msg_flag)
                        then CCTK_VInfo(CCTK_THORNSTRING,
                                        "      xyz positions and derivative coefficients");

                  for (int pn = 0; pn < ps.N_patches(); ++pn)
                  {
                        patch &p = ps.ith_patch(pn);

                        for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
                        {
                              for (int isigma = p.min_isigma();
                                   isigma <= p.max_isigma();
                                   ++isigma)
                              {
                                    const fp r = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
                                    const fp rho = p.rho_of_irho(irho);
                                    const fp sigma = p.sigma_of_isigma(isigma);

                                    fp local_x, local_y, local_z;
                                    p.xyz_of_r_rho_sigma(r, rho, sigma, local_x, local_y, local_z);

                                    const fp global_x = ps.origin_x() + local_x;
                                    const fp global_y = ps.origin_y() + local_y;
                                    const fp global_z = ps.origin_z() + local_z;

                                    p.gridfn(gfns::gfn__global_x, irho, isigma) = global_x;
                                    p.gridfn(gfns::gfn__global_y, irho, isigma) = global_y;
                                    p.gridfn(gfns::gfn__global_z, irho, isigma) = global_z;

                                    const fp global_xx = global_x * global_x;
                                    const fp global_xy = global_x * global_y;
                                    const fp global_xz = global_x * global_z;
                                    const fp global_yy = global_y * global_y;
                                    const fp global_yz = global_y * global_z;
                                    const fp global_zz = global_z * global_z;

                                    p.gridfn(gfns::gfn__global_xx, irho, isigma) = global_xx;
                                    p.gridfn(gfns::gfn__global_xy, irho, isigma) = global_xy;
                                    p.gridfn(gfns::gfn__global_xz, irho, isigma) = global_xz;
                                    p.gridfn(gfns::gfn__global_yy, irho, isigma) = global_yy;
                                    p.gridfn(gfns::gfn__global_yz, irho, isigma) = global_yz;
                                    p.gridfn(gfns::gfn__global_zz, irho, isigma) = global_zz;
                              }
                        }
                  }
            }
      }

      //******************************************************************************
      namespace
      {
            enum expansion_status
            interpolate_geometry(patch_system *ps_ptr,
                                 bool initial_flag,
                                 bool print_msg_flag)
            {
                  int status = 1;

#define CAST_PTR_OR_NULL(type_, ptr_) \
      (ps_ptr == NULL) ? NULL : static_cast<type_>(ptr_)

                  //
                  // ***** interpolation points *****
                  //
                  const int N_interp_points = (ps_ptr == NULL) ? 0 : ps_ptr->N_grid_points();
                  double *interp_coords[3] = {
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_x)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_y)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__global_z)),
                  };

                  double *const output_arrays[] = {
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__g_dd_11)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_111)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_211)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_311)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__g_dd_12)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_112)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_212)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_312)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__g_dd_13)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_113)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_213)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_313)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__g_dd_22)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_122)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_222)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_322)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__g_dd_23)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_123)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_223)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_323)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__g_dd_33)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_133)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_233)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_g_dd_333)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__psi)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_1)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_2)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__partial_d_psi_3)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__K_dd_11)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__K_dd_12)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__K_dd_13)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__K_dd_22)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__K_dd_23)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__K_dd_33)),
                      CAST_PTR_OR_NULL(double *, ps_ptr->gridfn_data(gfns::gfn__trK)),
                  };

                  const int N_output_arrays_dim = sizeof(output_arrays) / sizeof(output_arrays[0]);
                  const int N_output_arrays_use = N_output_arrays_dim;

                  int s;
                  int Npts = 0;
                  for (int ncpu = 0; ncpu < state.N_procs; ncpu++)
                  {

                        if (state.my_proc == ncpu)
                              Npts = N_interp_points;

                        MPI_Bcast(&Npts, 1, MPI_INT, ncpu, MPI_COMM_WORLD);

                        if (Npts != 0)
                        {
                              if (state.my_proc == ncpu)
                              {
                                    memcpy(state.oX, interp_coords[0], Npts * sizeof(double));
                                    memcpy(state.oY, interp_coords[1], Npts * sizeof(double));
                                    memcpy(state.oZ, interp_coords[2], Npts * sizeof(double));
                              }
                              MPI_Bcast(state.oX, Npts, MPI_DOUBLE, ncpu, MPI_COMM_WORLD);
                              MPI_Bcast(state.oY, Npts, MPI_DOUBLE, ncpu, MPI_COMM_WORLD);
                              MPI_Bcast(state.oZ, Npts, MPI_DOUBLE, ncpu, MPI_COMM_WORLD);

                              // each cpu calls interpolator
                              s = globalInterpGFL(state.oX, state.oY, state.oZ, Npts, state.Data); // 1 succuss; 0 fail

                              if (state.my_proc == ncpu)
                              {
                                    status = s;

                                    if (status == 1)
                                    {
                                          for (int ngf = 0; ngf < N_output_arrays_use; ngf++)
                                          {
                                                memcpy(output_arrays[ngf], state.Data + ngf * N_interp_points,
                                                       sizeof(double) * N_interp_points);
                                          }
                                    }
                                    else
                                    {
                                          char filename[100];
                                          sprintf(filename, "check%05d.dat", state.my_proc);
                                          if (ps_ptr)
                                                ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_11, true, gfns::gfn__h, filename);
                                          //           MPI_Abort(MPI_COMM_WORLD,1);
                                          return expansion_failure__surface_outside_grid;
                                    }
                              }
                        }
                  }

#if 0
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_11,true,gfns::gfn__h,"check.dat");
           char filename[100];
           sprintf(filename,"g311%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_311,true,gfns::gfn__h,filename);
           sprintf(filename,"g12%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_12,true,gfns::gfn__h,filename);
           sprintf(filename,"g112%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_112,true,gfns::gfn__h,filename);
           sprintf(filename,"g212%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_212,true,gfns::gfn__h,filename);
           sprintf(filename,"g312%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_312,true,gfns::gfn__h,filename);
           sprintf(filename,"g13%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_13,true,gfns::gfn__h,filename);
           sprintf(filename,"g113%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_113,true,gfns::gfn__h,filename);
           sprintf(filename,"g213%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_213,true,gfns::gfn__h,filename);
           sprintf(filename,"g313%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_313,true,gfns::gfn__h,filename);
           sprintf(filename,"g22%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_22,true,gfns::gfn__h,filename);
           sprintf(filename,"g122%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_122,true,gfns::gfn__h,filename);
           sprintf(filename,"g222%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_222,true,gfns::gfn__h,filename);
           sprintf(filename,"g322%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_322,true,gfns::gfn__h,filename);
           sprintf(filename,"g23%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_23,true,gfns::gfn__h,filename);
           sprintf(filename,"g123%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_123,true,gfns::gfn__h,filename);
           sprintf(filename,"g223%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_223,true,gfns::gfn__h,filename);
           sprintf(filename,"g323%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_323,true,gfns::gfn__h,filename);
           sprintf(filename,"g33%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__g_dd_33,true,gfns::gfn__h,filename);
           sprintf(filename,"g133%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_133,true,gfns::gfn__h,filename);
           sprintf(filename,"g233%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_233,true,gfns::gfn__h,filename);
           sprintf(filename,"g333%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_g_dd_333,true,gfns::gfn__h,filename);
           sprintf(filename,"psi%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__psi,true,gfns::gfn__h,filename);
           sprintf(filename,"psi1%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_psi_1,true,gfns::gfn__h,filename);
           sprintf(filename,"psi2%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_psi_2,true,gfns::gfn__h,filename);
           sprintf(filename,"psi3%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__partial_d_psi_3,true,gfns::gfn__h,filename);
           sprintf(filename,"K11%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__K_dd_11,true,gfns::gfn__h,filename);
           sprintf(filename,"K12%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__K_dd_12,true,gfns::gfn__h,filename);
           sprintf(filename,"K13%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__K_dd_13,true,gfns::gfn__h,filename);
           sprintf(filename,"K22%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__K_dd_22,true,gfns::gfn__h,filename);
           sprintf(filename,"K23%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__K_dd_23,true,gfns::gfn__h,filename);
           sprintf(filename,"K33%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__K_dd_33,true,gfns::gfn__h,filename);
           sprintf(filename,"trK%02d.dat",state.my_proc);
  if(ps_ptr) ps_ptr->print_gridfn_with_xyz(gfns::gfn__trK,true,gfns::gfn__h,filename);

  MPI_Abort(MPI_COMM_WORLD,1);
#endif

                  if (status == 0)
                        then error_exit(ERROR_EXIT,
                                        "***** interpolate_geometry(): error return %d from interpolator!\n",
                                        status); /*NOTREACHED*/

                  return expansion_success; // *** NORMAL RETURN ***
            }
      }

      //******************************************************************************
      namespace
      {
            void convert_conformal_to_physical(patch_system &ps, bool print_msg_flag)
            {
                  for (int pn = 0; pn < ps.N_patches(); ++pn)
                  {
                        patch &p = ps.ith_patch(pn);

                        for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
                              for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma)
                              {

                                    const fp psi = (p.gridfn(gfns::gfn__psi, irho, isigma));
                                    const fp psi3 = jtutil::pow3(psi);
                                    const fp psi4 = jtutil::pow4(psi);

                                    const fp partial_d_psi_1 = p.gridfn(gfns::gfn__partial_d_psi_1, irho, isigma);
                                    const fp partial_d_psi_2 = p.gridfn(gfns::gfn__partial_d_psi_2, irho, isigma);
                                    const fp partial_d_psi_3 = p.gridfn(gfns::gfn__partial_d_psi_3, irho, isigma);

                                    const fp stored_g_dd_11 = p.gridfn(gfns::gfn__g_dd_11, irho, isigma);
                                    const fp stored_g_dd_12 = p.gridfn(gfns::gfn__g_dd_12, irho, isigma);
                                    const fp stored_g_dd_13 = p.gridfn(gfns::gfn__g_dd_13, irho, isigma);
                                    const fp stored_g_dd_22 = p.gridfn(gfns::gfn__g_dd_22, irho, isigma);
                                    const fp stored_g_dd_23 = p.gridfn(gfns::gfn__g_dd_23, irho, isigma);
                                    const fp stored_g_dd_33 = p.gridfn(gfns::gfn__g_dd_33, irho, isigma);

                                    p.gridfn(gfns::gfn__g_dd_11, irho, isigma) *= psi4;
                                    p.gridfn(gfns::gfn__g_dd_12, irho, isigma) *= psi4;
                                    p.gridfn(gfns::gfn__g_dd_13, irho, isigma) *= psi4;
                                    p.gridfn(gfns::gfn__g_dd_22, irho, isigma) *= psi4;
                                    p.gridfn(gfns::gfn__g_dd_23, irho, isigma) *= psi4;
                                    p.gridfn(gfns::gfn__g_dd_33, irho, isigma) *= psi4;

                                    p.gridfn(gfns::gfn__partial_d_g_dd_111, irho, isigma) = 4.0 * psi3 * partial_d_psi_1 * stored_g_dd_11 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_111, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_112, irho, isigma) = 4.0 * psi3 * partial_d_psi_1 * stored_g_dd_12 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_112, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_113, irho, isigma) = 4.0 * psi3 * partial_d_psi_1 * stored_g_dd_13 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_113, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_122, irho, isigma) = 4.0 * psi3 * partial_d_psi_1 * stored_g_dd_22 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_122, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_123, irho, isigma) = 4.0 * psi3 * partial_d_psi_1 * stored_g_dd_23 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_123, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_133, irho, isigma) = 4.0 * psi3 * partial_d_psi_1 * stored_g_dd_33 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_133, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_211, irho, isigma) = 4.0 * psi3 * partial_d_psi_2 * stored_g_dd_11 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_211, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_212, irho, isigma) = 4.0 * psi3 * partial_d_psi_2 * stored_g_dd_12 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_212, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_213, irho, isigma) = 4.0 * psi3 * partial_d_psi_2 * stored_g_dd_13 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_213, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_222, irho, isigma) = 4.0 * psi3 * partial_d_psi_2 * stored_g_dd_22 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_222, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_223, irho, isigma) = 4.0 * psi3 * partial_d_psi_2 * stored_g_dd_23 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_223, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_233, irho, isigma) = 4.0 * psi3 * partial_d_psi_2 * stored_g_dd_33 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_233, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_311, irho, isigma) = 4.0 * psi3 * partial_d_psi_3 * stored_g_dd_11 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_311, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_312, irho, isigma) = 4.0 * psi3 * partial_d_psi_3 * stored_g_dd_12 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_312, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_313, irho, isigma) = 4.0 * psi3 * partial_d_psi_3 * stored_g_dd_13 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_313, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_322, irho, isigma) = 4.0 * psi3 * partial_d_psi_3 * stored_g_dd_22 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_322, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_323, irho, isigma) = 4.0 * psi3 * partial_d_psi_3 * stored_g_dd_23 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_323, irho, isigma);
                                    p.gridfn(gfns::gfn__partial_d_g_dd_333, irho, isigma) = 4.0 * psi3 * partial_d_psi_3 * stored_g_dd_33 + psi4 * p.gridfn(gfns::gfn__partial_d_g_dd_333, irho, isigma);

                                    // K_ij = psi4 \tilde{A}_ij + (1/3) g_ij TrK,    g_ij = psi4 \tilde{g}_ij
                                    const fp stored_trKo3 = p.gridfn(gfns::gfn__trK, irho, isigma) / 3.0;
                                    const fp stored_K_dd_11 = p.gridfn(gfns::gfn__K_dd_11, irho, isigma);
                                    const fp stored_K_dd_12 = p.gridfn(gfns::gfn__K_dd_12, irho, isigma);
                                    const fp stored_K_dd_13 = p.gridfn(gfns::gfn__K_dd_13, irho, isigma);
                                    const fp stored_K_dd_22 = p.gridfn(gfns::gfn__K_dd_22, irho, isigma);
                                    const fp stored_K_dd_23 = p.gridfn(gfns::gfn__K_dd_23, irho, isigma);
                                    const fp stored_K_dd_33 = p.gridfn(gfns::gfn__K_dd_33, irho, isigma);

                                    p.gridfn(gfns::gfn__K_dd_11, irho, isigma) = psi4 *
                                                                                 (stored_K_dd_11 + stored_g_dd_11 * stored_trKo3);
                                    p.gridfn(gfns::gfn__K_dd_12, irho, isigma) = psi4 *
                                                                                 (stored_K_dd_12 + stored_g_dd_12 * stored_trKo3);
                                    p.gridfn(gfns::gfn__K_dd_13, irho, isigma) = psi4 *
                                                                                 (stored_K_dd_13 + stored_g_dd_13 * stored_trKo3);
                                    p.gridfn(gfns::gfn__K_dd_22, irho, isigma) = psi4 *
                                                                                 (stored_K_dd_22 + stored_g_dd_22 * stored_trKo3);
                                    p.gridfn(gfns::gfn__K_dd_23, irho, isigma) = psi4 *
                                                                                 (stored_K_dd_23 + stored_g_dd_23 * stored_trKo3);
                                    p.gridfn(gfns::gfn__K_dd_33, irho, isigma) = psi4 *
                                                                                 (stored_K_dd_33 + stored_g_dd_33 * stored_trKo3);

                              } // end for irho isigma
                  }
            }
      }

      namespace
      {
            bool h_is_finite(patch_system &ps, bool initial_flag,
                             bool print_msg_flag)
            {
                  if (print_msg_flag)
                        then CCTK_VInfo(CCTK_THORNSTRING, "      checking that h is finite");

                  for (int pn = 0; pn < ps.N_patches(); ++pn)
                  {
                        patch &p = ps.ith_patch(pn);

                        for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
                        {
                              for (int isigma = p.min_isigma();
                                   isigma <= p.max_isigma();
                                   ++isigma)
                              {
                                    const fp h = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
                                    if (!finite(h))
                                          then
                                          {
                                                const fp rho = p.rho_of_irho(irho);
                                                const fp sigma = p.sigma_of_isigma(isigma);
                                                const fp drho = jtutil::degrees_of_radians(rho);
                                                const fp dsigma = jtutil::degrees_of_radians(sigma);
                                                CCTK_VWarn(1,
                                                           __LINE__, __FILE__, CCTK_THORNSTRING,
                                                           "\n"
                                                           "   h=%g isn't finite!\n"
                                                           "   %s patch (rho,sigma)=(%g,%g) (drho,dsigma)=(%g,%g)\n",
                                                           double(h),
                                                           p.name(), double(rho), double(sigma),
                                                           double(drho), double(dsigma));
                                                return false; // *** found a NaN ***
                                          }
                              }
                        }
                  }
                  return true; // *** all values finite ***
            }
      }

      //******************************************************************************
      namespace
      {
            bool geometry_is_finite(patch_system &ps, bool initial_flag,
                                    bool print_msg_flag)
            {
                  if (print_msg_flag)
                        then CCTK_VInfo(CCTK_THORNSTRING, "      checking that geometry is finite");

                  for (int pn = 0; pn < ps.N_patches(); ++pn)
                  {
                        patch &p = ps.ith_patch(pn);

                        for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
                        {
                              for (int isigma = p.min_isigma();
                                   isigma <= p.max_isigma();
                                   ++isigma)
                              {
                                    const fp g_dd_11 = p.gridfn(gfns::gfn__g_dd_11, irho, isigma);
                                    const fp g_dd_12 = p.gridfn(gfns::gfn__g_dd_12, irho, isigma);
                                    const fp g_dd_13 = p.gridfn(gfns::gfn__g_dd_13, irho, isigma);
                                    const fp g_dd_22 = p.gridfn(gfns::gfn__g_dd_22, irho, isigma);
                                    const fp g_dd_23 = p.gridfn(gfns::gfn__g_dd_23, irho, isigma);
                                    const fp g_dd_33 = p.gridfn(gfns::gfn__g_dd_33, irho, isigma);

                                    const fp K_dd_11 = p.gridfn(gfns::gfn__K_dd_11, irho, isigma);
                                    const fp K_dd_12 = p.gridfn(gfns::gfn__K_dd_12, irho, isigma);
                                    const fp K_dd_13 = p.gridfn(gfns::gfn__K_dd_13, irho, isigma);
                                    const fp K_dd_22 = p.gridfn(gfns::gfn__K_dd_22, irho, isigma);
                                    const fp K_dd_23 = p.gridfn(gfns::gfn__K_dd_23, irho, isigma);
                                    const fp K_dd_33 = p.gridfn(gfns::gfn__K_dd_33, irho, isigma);

                                    const fp partial_d_g_dd_111 = p.gridfn(gfns::gfn__partial_d_g_dd_111, irho, isigma);
                                    const fp partial_d_g_dd_112 = p.gridfn(gfns::gfn__partial_d_g_dd_112, irho, isigma);
                                    const fp partial_d_g_dd_113 = p.gridfn(gfns::gfn__partial_d_g_dd_113, irho, isigma);
                                    const fp partial_d_g_dd_122 = p.gridfn(gfns::gfn__partial_d_g_dd_122, irho, isigma);
                                    const fp partial_d_g_dd_123 = p.gridfn(gfns::gfn__partial_d_g_dd_123, irho, isigma);
                                    const fp partial_d_g_dd_133 = p.gridfn(gfns::gfn__partial_d_g_dd_133, irho, isigma);
                                    const fp partial_d_g_dd_211 = p.gridfn(gfns::gfn__partial_d_g_dd_211, irho, isigma);
                                    const fp partial_d_g_dd_212 = p.gridfn(gfns::gfn__partial_d_g_dd_212, irho, isigma);
                                    const fp partial_d_g_dd_213 = p.gridfn(gfns::gfn__partial_d_g_dd_213, irho, isigma);
                                    const fp partial_d_g_dd_222 = p.gridfn(gfns::gfn__partial_d_g_dd_222, irho, isigma);
                                    const fp partial_d_g_dd_223 = p.gridfn(gfns::gfn__partial_d_g_dd_223, irho, isigma);
                                    const fp partial_d_g_dd_233 = p.gridfn(gfns::gfn__partial_d_g_dd_233, irho, isigma);
                                    const fp partial_d_g_dd_311 = p.gridfn(gfns::gfn__partial_d_g_dd_311, irho, isigma);
                                    const fp partial_d_g_dd_312 = p.gridfn(gfns::gfn__partial_d_g_dd_312, irho, isigma);
                                    const fp partial_d_g_dd_313 = p.gridfn(gfns::gfn__partial_d_g_dd_313, irho, isigma);
                                    const fp partial_d_g_dd_322 = p.gridfn(gfns::gfn__partial_d_g_dd_322, irho, isigma);
                                    const fp partial_d_g_dd_323 = p.gridfn(gfns::gfn__partial_d_g_dd_323, irho, isigma);
                                    const fp partial_d_g_dd_333 = p.gridfn(gfns::gfn__partial_d_g_dd_333, irho, isigma);

                                    if (!finite(g_dd_11) || !finite(g_dd_12) || !finite(g_dd_13) || !finite(g_dd_22) || !finite(g_dd_23) || !finite(g_dd_33) || !finite(K_dd_11) || !finite(K_dd_12) || !finite(K_dd_13) || !finite(K_dd_22) || !finite(K_dd_23) || !finite(K_dd_33) || !finite(partial_d_g_dd_111) || !finite(partial_d_g_dd_112) || !finite(partial_d_g_dd_113) || !finite(partial_d_g_dd_122) || !finite(partial_d_g_dd_123) || !finite(partial_d_g_dd_133) || !finite(partial_d_g_dd_211) || !finite(partial_d_g_dd_212) || !finite(partial_d_g_dd_213) || !finite(partial_d_g_dd_222) || !finite(partial_d_g_dd_223) || !finite(partial_d_g_dd_233) || !finite(partial_d_g_dd_311) || !finite(partial_d_g_dd_312) || !finite(partial_d_g_dd_313) || !finite(partial_d_g_dd_322) || !finite(partial_d_g_dd_323) || !finite(partial_d_g_dd_333))
                                          then
                                          {
                                                const fp h = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
                                                const fp rho = p.rho_of_irho(irho);
                                                const fp sigma = p.sigma_of_isigma(isigma);
                                                const fp drho = jtutil::degrees_of_radians(rho);
                                                const fp dsigma = jtutil::degrees_of_radians(sigma);
                                                fp local_x, local_y, local_z;
                                                p.xyz_of_r_rho_sigma(h, rho, sigma, local_x, local_y, local_z);
                                                const fp global_x = ps.origin_x() + local_x;
                                                const fp global_y = ps.origin_y() + local_y;
                                                const fp global_z = ps.origin_z() + local_z;
                                                CCTK_VWarn(1,
                                                           __LINE__, __FILE__, CCTK_THORNSTRING,
                                                           "\n"
                                                           "   geometry isn't finite at %s patch\n"
                                                           "   h=%g (rho,sigma)=(%g,%g) (drho,dsigma)=(%g,%g)\n"
                                                           "   local_(x,y,z)=(%g,%g,%g)\n"
                                                           "   global_(x,y,z)=(%g,%g,%g)\n"
                                                           "   g_dd_11=%g   _12=%g   _13=%g\n"
                                                           "       _22=%g   _23=%g   _33=%g\n"
                                                           "   K_dd_11=%g   _12=%g   _13=%g\n"
                                                           "       _22=%g   _23=%g   _33=%g\n"
                                                           "   partial_d_g_dd_111=%g   _112=%g   _113=%g\n"
                                                           "                 _122=%g   _123=%g   _133=%g\n"
                                                           "   partial_d_g_dd_211=%g   _212=%g   _213=%g\n"
                                                           "                 _222=%g   _223=%g   _233=%g\n"
                                                           "   partial_d_g_dd_311=%g   _312=%g   _313=%g\n"
                                                           "                 _322=%g   _323=%g   _333=%g\n",
                                                           p.name(),
                                                           double(h), double(rho), double(sigma),
                                                           double(drho), double(dsigma),
                                                           double(local_x), double(local_y), double(local_z),
                                                           double(global_x), double(global_y), double(global_z),
                                                           double(g_dd_11), double(g_dd_12), double(g_dd_13),
                                                           double(g_dd_22), double(g_dd_23), double(g_dd_33),
                                                           double(K_dd_11), double(K_dd_12), double(K_dd_13),
                                                           double(K_dd_22), double(K_dd_23), double(K_dd_33),
                                                           double(partial_d_g_dd_111),
                                                           double(partial_d_g_dd_112),
                                                           double(partial_d_g_dd_113),
                                                           double(partial_d_g_dd_122),
                                                           double(partial_d_g_dd_123),
                                                           double(partial_d_g_dd_133),
                                                           double(partial_d_g_dd_211),
                                                           double(partial_d_g_dd_212),
                                                           double(partial_d_g_dd_213),
                                                           double(partial_d_g_dd_222),
                                                           double(partial_d_g_dd_223),
                                                           double(partial_d_g_dd_233),
                                                           double(partial_d_g_dd_311),
                                                           double(partial_d_g_dd_312),
                                                           double(partial_d_g_dd_313),
                                                           double(partial_d_g_dd_322),
                                                           double(partial_d_g_dd_323),
                                                           double(partial_d_g_dd_333));
                                                return false; // *** found a NaN ***
                                          }
                              }
                        }
                  }
                  return true; // *** no NaNs found ***
            }
      }

      //******************************************************************************
      //******************************************************************************
      //******************************************************************************

      //
      // This function computes the expansion Theta(h), and optionally also
      // its Jacobian coefficients, (from which the Jacobian matrix may be
      // computed later).  This function uses a mixture of algebraic operations
      // and (rho,sigma) finite differencing.  The computation is done entirely
      // on the nominal angular grid.
      //
      // N.b. This function #includes "cg.hh", which defines "dangerous" macros
      //      which will stay in effect for the rest of this compilation unit!
      //
      // Arguments:
      // Jacobian_flag = true to compute the Jacobian coefficients,
      //		   false to skip this.
      //
      // Results:
      // This function returns true for a successful computation, or false
      // if the computation failed because Theta_D <= 0 (this means the interpolated
      // g_ij isn't positive definite).
      //
      namespace
      {
            bool compute_Theta(patch_system &ps, fp add_to_expansion,
                               bool Jacobian_flag, jtutil::norm<fp> *Theta_norms_ptr,
                               bool initial_flag,
                               bool print_msg_flag)
            {
                  if (print_msg_flag)
                        then CCTK_VInfo(CCTK_THORNSTRING, "      computing Theta(h)");

                  for (int pn = 0; pn < ps.N_patches(); ++pn)
                  {
                        patch &p = ps.ith_patch(pn);

                        for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
                        {
                              for (int isigma = p.min_isigma();
                                   isigma <= p.max_isigma();
                                   ++isigma)
                              {
                                    //
                                    // compute the X_ud and X_udd derivative coefficients
                                    // ... n.b. this uses the *local* (x,y,z) coordinates
                                    //
                                    const fp r = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
                                    const fp rho = p.rho_of_irho(irho);
                                    const fp sigma = p.sigma_of_isigma(isigma);
                                    fp xx, yy, zz;
                                    p.xyz_of_r_rho_sigma(r, rho, sigma, xx, yy, zz);

                                    // 1st derivative coefficients X_ud
                                    const fp X_ud_11 = p.partial_rho_wrt_x(xx, yy, zz);
                                    const fp X_ud_12 = p.partial_rho_wrt_y(xx, yy, zz);
                                    const fp X_ud_13 = p.partial_rho_wrt_z(xx, yy, zz);
                                    const fp X_ud_21 = p.partial_sigma_wrt_x(xx, yy, zz);
                                    const fp X_ud_22 = p.partial_sigma_wrt_y(xx, yy, zz);
                                    const fp X_ud_23 = p.partial_sigma_wrt_z(xx, yy, zz);

                                    // 2nd derivative coefficient gridfns X_udd
                                    const fp X_udd_111 = p.partial2_rho_wrt_xx(xx, yy, zz);
                                    const fp X_udd_112 = p.partial2_rho_wrt_xy(xx, yy, zz);
                                    const fp X_udd_113 = p.partial2_rho_wrt_xz(xx, yy, zz);
                                    const fp X_udd_122 = p.partial2_rho_wrt_yy(xx, yy, zz);
                                    const fp X_udd_123 = p.partial2_rho_wrt_yz(xx, yy, zz);
                                    const fp X_udd_133 = p.partial2_rho_wrt_zz(xx, yy, zz);
                                    const fp X_udd_211 = p.partial2_sigma_wrt_xx(xx, yy, zz);
                                    const fp X_udd_212 = p.partial2_sigma_wrt_xy(xx, yy, zz);
                                    const fp X_udd_213 = p.partial2_sigma_wrt_xz(xx, yy, zz);
                                    const fp X_udd_222 = p.partial2_sigma_wrt_yy(xx, yy, zz);
                                    const fp X_udd_223 = p.partial2_sigma_wrt_yz(xx, yy, zz);
                                    const fp X_udd_233 = p.partial2_sigma_wrt_zz(xx, yy, zz);

#define RATIONAL(num, den) (num / den)

#define PARTIAL_RHO(ghosted_gridfn_name) \
      p.partial_rho(gfns::gfn__##ghosted_gridfn_name, irho, isigma)
#define PARTIAL_SIGMA(ghosted_gridfn_name) \
      p.partial_sigma(gfns::gfn__##ghosted_gridfn_name, irho, isigma)
#define PARTIAL_RHO_RHO(ghosted_gridfn_name) \
      p.partial_rho_rho(gfns::gfn__##ghosted_gridfn_name, irho, isigma)
#define PARTIAL_RHO_SIGMA(ghosted_gridfn_name) \
      p.partial_rho_sigma(gfns::gfn__##ghosted_gridfn_name, irho, isigma)
#define PARTIAL_SIGMA_SIGMA(ghosted_gridfn_name) \
      p.partial_sigma_sigma(gfns::gfn__##ghosted_gridfn_name, irho, isigma)

#define h p.ghosted_gridfn(gfns::gfn__h, irho, isigma)
#define r h

#define g_dd_11 p.gridfn(gfns::gfn__g_dd_11, irho, isigma)
#define g_dd_12 p.gridfn(gfns::gfn__g_dd_12, irho, isigma)
#define g_dd_13 p.gridfn(gfns::gfn__g_dd_13, irho, isigma)
#define g_dd_22 p.gridfn(gfns::gfn__g_dd_22, irho, isigma)
#define g_dd_23 p.gridfn(gfns::gfn__g_dd_23, irho, isigma)
#define g_dd_33 p.gridfn(gfns::gfn__g_dd_33, irho, isigma)
#define K_dd_11 p.gridfn(gfns::gfn__K_dd_11, irho, isigma)
#define K_dd_12 p.gridfn(gfns::gfn__K_dd_12, irho, isigma)
#define K_dd_13 p.gridfn(gfns::gfn__K_dd_13, irho, isigma)
#define K_dd_22 p.gridfn(gfns::gfn__K_dd_22, irho, isigma)
#define K_dd_23 p.gridfn(gfns::gfn__K_dd_23, irho, isigma)
#define K_dd_33 p.gridfn(gfns::gfn__K_dd_33, irho, isigma)

#define partial_d_g_dd_111 p.gridfn(gfns::gfn__partial_d_g_dd_111, irho, isigma)
#define partial_d_g_dd_112 p.gridfn(gfns::gfn__partial_d_g_dd_112, irho, isigma)
#define partial_d_g_dd_113 p.gridfn(gfns::gfn__partial_d_g_dd_113, irho, isigma)
#define partial_d_g_dd_122 p.gridfn(gfns::gfn__partial_d_g_dd_122, irho, isigma)
#define partial_d_g_dd_123 p.gridfn(gfns::gfn__partial_d_g_dd_123, irho, isigma)
#define partial_d_g_dd_133 p.gridfn(gfns::gfn__partial_d_g_dd_133, irho, isigma)
#define partial_d_g_dd_211 p.gridfn(gfns::gfn__partial_d_g_dd_211, irho, isigma)
#define partial_d_g_dd_212 p.gridfn(gfns::gfn__partial_d_g_dd_212, irho, isigma)
#define partial_d_g_dd_213 p.gridfn(gfns::gfn__partial_d_g_dd_213, irho, isigma)
#define partial_d_g_dd_222 p.gridfn(gfns::gfn__partial_d_g_dd_222, irho, isigma)
#define partial_d_g_dd_223 p.gridfn(gfns::gfn__partial_d_g_dd_223, irho, isigma)
#define partial_d_g_dd_233 p.gridfn(gfns::gfn__partial_d_g_dd_233, irho, isigma)
#define partial_d_g_dd_311 p.gridfn(gfns::gfn__partial_d_g_dd_311, irho, isigma)
#define partial_d_g_dd_312 p.gridfn(gfns::gfn__partial_d_g_dd_312, irho, isigma)
#define partial_d_g_dd_313 p.gridfn(gfns::gfn__partial_d_g_dd_313, irho, isigma)
#define partial_d_g_dd_322 p.gridfn(gfns::gfn__partial_d_g_dd_322, irho, isigma)
#define partial_d_g_dd_323 p.gridfn(gfns::gfn__partial_d_g_dd_323, irho, isigma)
#define partial_d_g_dd_333 p.gridfn(gfns::gfn__partial_d_g_dd_333, irho, isigma)

#define Theta p.gridfn(gfns::gfn__Theta, irho, isigma)

#define partial_Theta_wrt_partial_d_h_1 \
      p.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_1, irho, isigma)
#define partial_Theta_wrt_partial_d_h_2 \
      p.gridfn(gfns::gfn__partial_Theta_wrt_partial_d_h_2, irho, isigma)
#define partial_Theta_wrt_partial_dd_h_11 \
      p.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_11, irho, isigma)
#define partial_Theta_wrt_partial_dd_h_12 \
      p.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_12, irho, isigma)
#define partial_Theta_wrt_partial_dd_h_22 \
      p.gridfn(gfns::gfn__partial_Theta_wrt_partial_dd_h_22, irho, isigma)

#define save_Theta p.gridfn(gfns::gfn__save_Theta, irho, isigma)
#define Delta_h p.gridfn(gfns::gfn__Delta_h, irho, isigma)

                                    fp g_uu_11;
                                    fp g_uu_12;
                                    fp g_uu_13;
                                    fp g_uu_22;
                                    fp g_uu_23;
                                    fp g_uu_33;
                                    fp K;
                                    fp K_uu_11;
                                    fp K_uu_12;
                                    fp K_uu_13;
                                    fp K_uu_22;
                                    fp K_uu_23;
                                    fp K_uu_33;

                                    fp partial_d_ln_sqrt_g_1;
                                    fp partial_d_ln_sqrt_g_2;
                                    fp partial_d_ln_sqrt_g_3;

                                    fp partial_d_g_uu_111;
                                    fp partial_d_g_uu_112;
                                    fp partial_d_g_uu_113;
                                    fp partial_d_g_uu_122;
                                    fp partial_d_g_uu_123;
                                    fp partial_d_g_uu_133;
                                    fp partial_d_g_uu_211;
                                    fp partial_d_g_uu_212;
                                    fp partial_d_g_uu_213;
                                    fp partial_d_g_uu_222;
                                    fp partial_d_g_uu_223;
                                    fp partial_d_g_uu_233;
                                    fp partial_d_g_uu_311;
                                    fp partial_d_g_uu_312;
                                    fp partial_d_g_uu_313;
                                    fp partial_d_g_uu_322;
                                    fp partial_d_g_uu_323;
                                    fp partial_d_g_uu_333;

                                    fp Theta_A;
                                    fp Theta_B;
                                    fp Theta_C;
                                    fp Theta_D;

                                    {
                                          // g_uu
                                          fp t1, t2, t4, t5, t7, t8, t11, t12, t14, t15;
                                          fp t18, t21;
                                          t1 = g_dd_22;
                                          t2 = g_dd_33;
                                          t4 = g_dd_23;
                                          t5 = t4 * t4;
                                          t7 = g_dd_11;
                                          t8 = t7 * t1;
                                          t11 = g_dd_12;
                                          t12 = t11 * t11;
                                          t14 = g_dd_13;
                                          t15 = t11 * t14;
                                          t18 = t14 * t14;
                                          t21 = 1 / (t8 * t2 - t7 * t5 - t12 * t2 + 2.0 * t15 * t4 - t18 * t1);
                                          g_uu_11 = (t1 * t2 - t5) * t21;
                                          g_uu_12 = -(t11 * t2 - t14 * t4) * t21;
                                          g_uu_13 = -(-t11 * t4 + t14 * t1) * t21;
                                          g_uu_22 = (t7 * t2 - t18) * t21;
                                          g_uu_23 = -(t7 * t4 - t15) * t21;
                                          g_uu_33 = (t8 - t12) * t21;
                                    }

                                    {
                                          // K, K_uu
                                          fp t1, t2, t4, t5, t8, t9, t12, t13, t15, t16;
                                          fp t19, t20, t22, t24, t27, t30, t32, t35, t42, t44;
                                          fp t46, t48, t50, t60, t62, t69, t71, t74, t85, t95;
                                          t1 = g_uu_11;
                                          t2 = K_dd_11;
                                          t4 = g_uu_12;
                                          t5 = K_dd_12;
                                          t8 = g_uu_13;
                                          t9 = K_dd_13;
                                          t12 = g_uu_22;
                                          t13 = K_dd_22;
                                          t15 = g_uu_23;
                                          t16 = K_dd_23;
                                          t19 = g_uu_33;
                                          t20 = K_dd_33;
                                          K = t1 * t2 + 2.0 * t4 * t5 + 2.0 * t8 * t9 + t12 * t13 + 2.0 * t15 * t16 + t19 * t20;
                                          t22 = t1 * t1;
                                          t24 = t4 * t1;
                                          t27 = t8 * t1;
                                          t30 = t4 * t4;
                                          t32 = t8 * t4;
                                          t35 = t8 * t8;
                                          K_uu_11 = t22 * t2 + 2.0 * t24 * t5 + 2.0 * t27 * t9 + t30 * t13 + 2.0 * t32 * t16 + t35 * t20;
                                          t42 = t4 * t12;
                                          t44 = t8 * t12;
                                          t46 = t1 * t15;
                                          t48 = t15 * t4;
                                          t50 = t8 * t15;
                                          K_uu_12 = t24 * t2 + t30 * t5 + t32 * t9 + t1 * t12 * t5 + t42 * t13 + t44 * t16 + t46 * t9 + t48 * t16 +
                                                    t50 * t20;
                                          t60 = t4 * t19;
                                          t62 = t8 * t19;
                                          K_uu_13 = t27 * t2 + t32 * t5 + t35 * t9 + t46 * t5 + t48 * t13 + t50 * t16 + t1 * t19 * t9 + t60 * t16 +
                                                    t62 * t20;
                                          t69 = t12 * t12;
                                          t71 = t15 * t12;
                                          t74 = t15 * t15;
                                          K_uu_22 = t30 * t2 + 2.0 * t42 * t5 + 2.0 * t48 * t9 + t69 * t13 + 2.0 * t71 * t16 + t74 * t20;
                                          t85 = t15 * t19;
                                          K_uu_23 = t32 * t2 + t44 * t5 + t50 * t9 + t48 * t5 + t71 * t13 + t74 * t16 + t60 * t9 + t12 * t19 * t16 +
                                                    t85 * t20;
                                          t95 = t19 * t19;
                                          K_uu_33 = t35 * t2 + 2.0 * t50 * t5 + 2.0 * t62 * t9 + t74 * t13 + 2.0 * t85 * t16 + t95 * t20;
                                    }

                                    {
                                          // partial_d_g_uu
                                          fp t1, t2, t3, t5, t6, t7, t10, t11, t12, t15;
                                          fp t16, t18, t19, t22, t23, t28, t29, t31, t33, t35;
                                          fp t36, t38, t40, t48, t49, t51, t53, t60, t62, t65;
                                          fp t74, t76, t86, t88, t90, t93, t96, t98, t101, t148;
                                          fp t150, t153, t156, t158, t161;
                                          t1 = g_uu_11;
                                          t2 = t1 * t1;
                                          t3 = partial_d_g_dd_111;
                                          t5 = g_uu_12;
                                          t6 = t5 * t1;
                                          t7 = partial_d_g_dd_112;
                                          t10 = g_uu_13;
                                          t11 = t10 * t1;
                                          t12 = partial_d_g_dd_113;
                                          t15 = t5 * t5;
                                          t16 = partial_d_g_dd_122;
                                          t18 = t10 * t5;
                                          t19 = partial_d_g_dd_123;
                                          t22 = t10 * t10;
                                          t23 = partial_d_g_dd_133;
                                          partial_d_g_uu_111 = -t2 * t3 - 2.0 * t6 * t7 - 2.0 * t11 * t12 - t15 * t16 - 2.0 * t18 * t19 - t22 * t23;
                                          t28 = g_uu_22;
                                          t29 = t1 * t28;
                                          t31 = t5 * t28;
                                          t33 = t10 * t28;
                                          t35 = g_uu_23;
                                          t36 = t1 * t35;
                                          t38 = t5 * t35;
                                          t40 = t10 * t35;
                                          partial_d_g_uu_112 = -t6 * t3 - t15 * t7 - t18 * t12 - t29 * t7 - t31 * t16 - t33 * t19 - t36 * t12 - t38 * t19 - t40 * t23;
                                          t48 = g_uu_33;
                                          t49 = t1 * t48;
                                          t51 = t48 * t5;
                                          t53 = t10 * t48;
                                          partial_d_g_uu_113 = -t11 * t3 - t18 * t7 - t22 * t12 - t36 * t7 - t38 * t16 - t40 * t19 - t49 * t12 - t51 * t19 - t53 * t23;
                                          t60 = t28 * t28;
                                          t62 = t35 * t28;
                                          t65 = t35 * t35;
                                          partial_d_g_uu_122 = -t15 * t3 - 2.0 * t31 * t7 - 2.0 * t38 * t12 - t60 * t16 - 2.0 * t62 * t19 -
                                                               t65 * t23;
                                          t74 = t28 * t48;
                                          t76 = t35 * t48;
                                          partial_d_g_uu_123 = -t18 * t3 - t33 * t7 - t40 * t12 - t38 * t7 - t62 * t16 - t65 * t19 - t51 * t12 - t74 * t19 - t76 * t23;
                                          t86 = t48 * t48;
                                          partial_d_g_uu_133 = -t22 * t3 - 2.0 * t40 * t7 - 2.0 * t53 * t12 - t65 * t16 - 2.0 * t76 * t19 -
                                                               t86 * t23;
                                          t88 = partial_d_g_dd_211;
                                          t90 = partial_d_g_dd_212;
                                          t93 = partial_d_g_dd_213;
                                          t96 = partial_d_g_dd_222;
                                          t98 = partial_d_g_dd_223;
                                          t101 = partial_d_g_dd_233;
                                          partial_d_g_uu_211 = -t2 * t88 - 2.0 * t6 * t90 - 2.0 * t11 * t93 - t15 * t96 - 2.0 * t18 * t98 -
                                                               t22 * t101;
                                          partial_d_g_uu_212 = -t6 * t88 - t15 * t90 - t18 * t93 - t29 * t90 - t31 * t96 - t33 * t98 - t36 * t93 - t38 * t98 - t40 * t101;
                                          partial_d_g_uu_213 = -t11 * t88 - t18 * t90 - t22 * t93 - t36 * t90 - t38 * t96 - t40 * t98 - t49 * t93 - t51 * t98 - t53 * t101;
                                          partial_d_g_uu_222 = -t15 * t88 - 2.0 * t31 * t90 - 2.0 * t38 * t93 - t60 * t96 - 2.0 * t62 * t98 - t65 * t101;
                                          partial_d_g_uu_223 = -t18 * t88 - t33 * t90 - t40 * t93 - t38 * t90 - t62 * t96 - t65 * t98 - t51 * t93 - t74 * t98 - t76 * t101;
                                          partial_d_g_uu_233 = -t22 * t88 - 2.0 * t40 * t90 - 2.0 * t53 * t93 - t65 * t96 - 2.0 * t76 * t98 - t86 * t101;
                                          t148 = partial_d_g_dd_311;
                                          t150 = partial_d_g_dd_312;
                                          t153 = partial_d_g_dd_313;
                                          t156 = partial_d_g_dd_322;
                                          t158 = partial_d_g_dd_323;
                                          t161 = partial_d_g_dd_333;
                                          partial_d_g_uu_311 = -t2 * t148 - 2.0 * t6 * t150 - 2.0 * t11 * t153 - t15 * t156 - 2.0 * t18 * t158 - t22 * t161;
                                          partial_d_g_uu_312 = -t6 * t148 - t15 * t150 - t18 * t153 - t29 * t150 - t31 * t156 - t33 * t158 - t36 * t153 - t38 * t158 - t40 * t161;
                                          partial_d_g_uu_313 = -t11 * t148 - t18 * t150 - t22 * t153 - t36 * t150 - t38 * t156 - t40 * t158 - t49 * t153 - t51 * t158 - t53 * t161;
                                          partial_d_g_uu_322 = -t15 * t148 - 2.0 * t31 * t150 - 2.0 * t38 * t153 - t60 * t156 - 2.0 * t62 * t158 - t65 * t161;
                                          partial_d_g_uu_323 = -t18 * t148 - t33 * t150 - t40 * t153 - t38 * t150 - t62 * t156 - t65 * t158 - t51 * t153 - t74 * t158 - t76 * t161;
                                          partial_d_g_uu_333 = -t22 * t148 - 2.0 * t40 * t150 - 2.0 * t53 * t153 - t65 * t156 - 2.0 * t76 * t158 - t86 * t161;
                                    }

                                    {
                                          // partial_d_ln_sqrt_g
                                          fp t1, t5, t8, t11, t15, t18;
                                          t1 = g_uu_11;
                                          t5 = g_uu_12;
                                          t8 = g_uu_13;
                                          t11 = g_uu_22;
                                          t15 = g_uu_23;
                                          t18 = g_uu_33;
                                          partial_d_ln_sqrt_g_1 = RATIONAL(1.0, 2.0) * t1 * partial_d_g_dd_111 + t5 * partial_d_g_dd_112 + t8 * partial_d_g_dd_113 + RATIONAL(1.0, 2.0) * t11 * partial_d_g_dd_122 + t15 * partial_d_g_dd_123 + RATIONAL(1.0, 2.0) * t18 * partial_d_g_dd_133;
                                          partial_d_ln_sqrt_g_2 = RATIONAL(1.0, 2.0) * t1 * partial_d_g_dd_211 + t5 * partial_d_g_dd_212 + t8 * partial_d_g_dd_213 + RATIONAL(1.0, 2.0) * t11 * partial_d_g_dd_222 + t15 * partial_d_g_dd_223 + RATIONAL(1.0, 2.0) * t18 * partial_d_g_dd_233;
                                          partial_d_ln_sqrt_g_3 = RATIONAL(1.0, 2.0) * t1 * partial_d_g_dd_311 + t5 * partial_d_g_dd_312 + t8 * partial_d_g_dd_313 + RATIONAL(1.0, 2.0) * t11 * partial_d_g_dd_322 + t15 * partial_d_g_dd_323 + RATIONAL(1.0, 2.0) * t18 * partial_d_g_dd_333;
                                    }

                                    {
                                          // Theta_A, Theta_B, Theta_C, Theta_D
                                          fp t1, t2, t3, t5, t6, t8, t9, t11, t12, t14;
                                          fp t15, t17, t19, t25, t26, t27, t29, t31, t34, t35;
                                          fp t37, t39, t40, t42, t44, t46, t47, t49, t56, t61;
                                          fp t63, t65, t66, t67, t82, t93, t98, t100, t102, t106;
                                          fp t107, t110, t111, t112, t116, t119, t120, t121, t123, t124;
                                          fp t127, t128, t129, t130, t131, t133, t134, t135, t137, t138;
                                          fp t139, t141, t142, t143, t148, t149, t150, t153, t154, t155;
                                          fp t158, t159, t160, t163, t164, t167, t168, t171, t172, t177;
                                          fp t181, t182, t185, t186, t189, t191, t197, t198, t200, t205;
                                          fp t220, t224, t232, t239, t266, t273, t276, t280, t283, t289;
                                          fp t292, t302, t303, t306, t307, t310, t311, t314, t317, t326;
                                          fp t330, t334, t337, t340, t343, t353, t355, t356, t360, t362;
                                          fp t366, t382, t387, t394, t431, t440, t444, t447, t450, t465;
                                          t1 = g_uu_13;
                                          t2 = t1 * t1;
                                          t3 = 1 / r;
                                          t5 = X_ud_13;
                                          t6 = PARTIAL_RHO(h);
                                          t8 = X_ud_23;
                                          t9 = PARTIAL_SIGMA(h);
                                          t11 = zz * t3 - t5 * t6 - t8 * t9;
                                          t12 = t11 * t11;
                                          t14 = yy * yy;
                                          t15 = zz * zz;
                                          t17 = r * r;
                                          t19 = 1 / t17 / r;
                                          t25 = X_ud_11;
                                          t26 = t25 * t25;
                                          t27 = PARTIAL_RHO_RHO(h);
                                          t29 = X_ud_21;
                                          t31 = PARTIAL_RHO_SIGMA(h);
                                          t34 = t29 * t29;
                                          t35 = PARTIAL_SIGMA_SIGMA(h);
                                          t37 = (t14 + t15) * t19 - X_udd_111 * t6 - X_udd_211 * t9 - t26 * t27 - 2.0 * t29 * t25 * t31 - t34 * t35;
                                          t39 = g_uu_23;
                                          t40 = t39 * t39;
                                          t42 = X_ud_12;
                                          t44 = X_ud_22;
                                          t46 = yy * t3 - t42 * t6 - t44 * t9;
                                          t47 = t46 * t46;
                                          t49 = xx * xx;
                                          t56 = t5 * t5;
                                          t61 = t8 * t8;
                                          t63 = (t49 + t14) * t19 - X_udd_133 * t6 - X_udd_233 * t9 - t56 * t27 - 2.0 * t8 * t5 * t31 - t61 * t35;
                                          t65 = t1 * t11;
                                          t66 = g_uu_22;
                                          t67 = t66 * t46;
                                          t82 = -xx * yy * t19 - X_udd_112 * t6 - X_udd_212 * t9 - t25 * t42 * t27 - t29 * t42 * t31 - t25 * t44 * t31 - t29 * t44 * t35;
                                          t93 = t42 * t42;
                                          t98 = t44 * t44;
                                          t100 = (t49 + t15) * t19 - X_udd_122 * t6 - X_udd_222 * t9 - t93 * t27 - 2.0 * t44 * t42 * t31 -
                                                 t98 * t35;
                                          t102 = t39 * t11;
                                          t106 = t1 * t12;
                                          t107 = partial_d_g_uu_123;
                                          t110 = g_uu_12;
                                          t111 = t110 * t47;
                                          t112 = partial_d_g_uu_112;
                                          t116 = xx * t3 - t25 * t6 - t29 * t9;
                                          t119 = t66 * t47;
                                          t120 = partial_d_g_uu_212;
                                          t121 = t120 * t116;
                                          t123 = t39 * t47;
                                          t124 = partial_d_g_uu_312;
                                          t127 = g_uu_11;
                                          t128 = t116 * t116;
                                          t129 = t127 * t128;
                                          t130 = partial_d_g_uu_113;
                                          t131 = t130 * t11;
                                          t133 = t1 * t128;
                                          t134 = partial_d_g_uu_313;
                                          t135 = t134 * t11;
                                          t137 = g_uu_33;
                                          t138 = t137 * t12;
                                          t139 = t134 * t116;
                                          t141 = -t2 * t12 * t37 - t40 * t47 * t63 - 2.0 * t65 * t67 * t82 - t40 * t12 * t100 - 2.0 * t102 * t67 * t100 - t106 * t107 * t46 - t111 * t112 * t116 - t119 * t121 - t123 * t124 * t116 - t129 * t131 - t133 * t135 -
                                                 t138 * t139;
                                          t142 = t39 * t12;
                                          t143 = partial_d_g_uu_213;
                                          t148 = t1 * t116;
                                          t149 = partial_d_g_uu_322;
                                          t150 = t149 * t47;
                                          t153 = t110 * t116;
                                          t154 = partial_d_g_uu_222;
                                          t155 = t154 * t47;
                                          t158 = t127 * t116;
                                          t159 = partial_d_g_uu_122;
                                          t160 = t159 * t47;
                                          t163 = partial_d_g_uu_333;
                                          t164 = t163 * t12;
                                          t167 = partial_d_g_uu_133;
                                          t168 = t167 * t12;
                                          t171 = partial_d_g_uu_233;
                                          t172 = t171 * t12;
                                          t177 = t110 * t46;
                                          t181 = partial_d_g_uu_323;
                                          t182 = t181 * t11;
                                          t185 = t137 * t11;
                                          t186 = t124 * t46;
                                          t189 = -t142 * t143 * t116 - t106 * t130 * t116 + RATIONAL(-1.0, 2.0) * t148 * t150 +
                                                 RATIONAL(-1.0, 2.0) * t153 * t155 + RATIONAL(-1.0, 2.0) * t158 * t160 + RATIONAL(-1.0, 2.0) * t148 * t164 + RATIONAL(-1.0, 2.0) * t158 * t168 + RATIONAL(-1.0, 2.0) * t153 * t172 + RATIONAL(-1.0, 2.0) * t65 * t160 - 2.0 * t65 * t177 * t37 - t148 * t182 * t46 - t185 * t186 * t116;
                                          t191 = t127 * t127;
                                          t197 = t110 * t128;
                                          t198 = t143 * t11;
                                          t200 = t137 * t137;
                                          t205 = t39 * t46;
                                          t220 = -xx * zz * t19 - X_udd_113 * t6 - X_udd_213 * t9 - t25 * t5 * t27 - t29 * t5 * t31 - t25 * t8 * t31 - t29 * t8 * t35;
                                          t224 = t12 * t11;
                                          t232 = t1 * t220;
                                          t239 = -t191 * t128 * t37 - 2.0 * t142 * t1 * t82 - t197 * t198 - t200 * t12 * t63 - t177 * t131 * t116 - 2.0 * t65 * t205 * t220 + RATIONAL(-1.0, 2.0) * t39 * t224 * t171 - t67 * t198 * t116 - t205 * t135 * t116 - 2.0 * t138 * t232 + RATIONAL(-1.0, 2.0) * t205 * t164 + RATIONAL(-1.0, 2.0) * t177 * t168;
                                          t266 = -yy * zz * t19 - X_udd_123 * t6 - X_udd_223 * t9 - t42 * t5 * t27 - t44 * t5 * t31 - t42 * t8 * t31 - t44 * t8 * t35;
                                          t273 = t110 * t110;
                                          t276 = t47 * t46;
                                          t280 = t39 * t266;
                                          t283 = t158 * t37;
                                          t289 = t148 * t266;
                                          t292 = RATIONAL(-1.0, 2.0) * t67 * t172 + RATIONAL(-1.0, 2.0) * t185 * t150 + RATIONAL(-1.0, 2.0) * t102 * t155 - 2.0 * t197 * t127 * t82 - 2.0 * t133 * t127 * t220 - 2.0 * t133 * t110 * t266 +
                                                 RATIONAL(-1.0, 2.0) * t1 * t224 * t167 - t273 * t128 * t100 + RATIONAL(-1.0, 2.0) * t39 * t276 * t149 - 2.0 * t138 * t280 - 2.0 * t65 * t283 + RATIONAL(-1.0, 2.0) * t110 * t276 * t159 - 2.0 * t67 * t289;
                                          t302 = partial_d_g_uu_311;
                                          t303 = t302 * t128;
                                          t306 = partial_d_g_uu_211;
                                          t307 = t306 * t128;
                                          t310 = partial_d_g_uu_111;
                                          t311 = t310 * t128;
                                          t314 = t148 * t63;
                                          t317 = t153 * t266;
                                          t326 = t107 * t11;
                                          t330 = RATIONAL(-1.0, 2.0) * t66 * t276 * t154 - 2.0 * t273 * t46 * t116 * t82 + RATIONAL(-1.0, 2.0) * t205 * t303 + RATIONAL(-1.0, 2.0) * t67 * t307 + RATIONAL(-1.0, 2.0) * t177 * t311 - 2.0 * t205 * t314 - 2.0 * t205 * t317 + RATIONAL(-1.0, 2.0) * t185 * t303 + RATIONAL(-1.0, 2.0) * t102 * t307 + RATIONAL(-1.0, 2.0) * t65 * t311 - t111 * t326 - t158 * t326 * t46;
                                          t334 = t158 * t82;
                                          t337 = t110 * t82;
                                          t340 = t158 * t220;
                                          t343 = t153 * t100;
                                          t353 = t112 * t46;
                                          t355 = partial_d_g_uu_223;
                                          t356 = t355 * t11;
                                          t360 = t120 * t46;
                                          t362 = -2.0 * t177 * t148 * t220 - 2.0 * t67 * t334 - 2.0 * t119 * t337 - 2.0 * t205 * t340 - 2.0 * t67 * t343 + RATIONAL(-1.0, 2.0) * t137 * t224 * t163 - t2 * t128 * t63 - t273 * t47 * t37 - t129 * t353 -
                                                 t119 * t356 - t123 * t182 - t133 * t186 - t197 * t360;
                                          t366 = t181 * t46;
                                          t382 = t66 * t66;
                                          t387 = t128 * t116;
                                          t394 = -t142 * t355 * t46 - t138 * t366 - 2.0 * t177 * t283 - 2.0 * t123 * t110 * t220 - 2.0 * t123 * t66 * t266 - t153 * t356 * t46 - t65 * t353 * t116 - t102 * t360 * t116 - t382 * t47 * t100 - 2.0 * t185 * t317 + RATIONAL(-1.0, 2.0) * t127 * t387 * t310 + RATIONAL(-1.0, 2.0) * t110 * t387 * t306;
                                          t431 = RATIONAL(-1.0, 2.0) * t1 * t387 * t302 - 2.0 * t2 * t11 * t116 * t220 - 2.0 * t185 * t314 - 2.0 * t102 * t289 - 2.0 * t65 * t153 * t82 - 2.0 * t185 * t205 * t63 - 2.0 * t40 * t11 * t46 * t266 - 2.0 * t102 * t343 - 2.0 * t102 * t334 - 2.0 * t185 * t340 - 2.0 * t102 * t177 * t82 - 2.0 * t185 * t67 * t266 - 2.0 * t185 * t177 * t220;
                                          Theta_A = t141 + t189 + t239 + t292 + t330 + t362 + t394 + t431;
                                          t440 = t310 * t116 + t121 + t139 + t353 + t154 * t46 + t366 + t131 + t356 + t163 * t11 + t127 * t37 + 2.0 * t337 + 2.0 * t232;
                                          t444 = partial_d_ln_sqrt_g_1;
                                          t447 = partial_d_ln_sqrt_g_2;
                                          t450 = partial_d_ln_sqrt_g_3;
                                          t465 = t66 * t100 + 2.0 * t280 + t137 * t63 + t127 * t444 * t116 + t110 * t447 * t116 + t1 * t450 * t116 + t110 * t444 * t46 + t66 * t447 * t46 + t39 * t450 * t46 + t1 * t444 * t11 + t39 * t447 * t11 + t137 * t450 * t11;
                                          Theta_B = t440 + t465;
                                          Theta_C = K_uu_11 * t128 + 2.0 * K_uu_12 * t46 * t116 + 2.0 * K_uu_13 * t11 * t116 + K_uu_22 * t47 + 2.0 * K_uu_23 * t11 * t46 + K_uu_33 * t12;
                                          Theta_D = t129 + 2.0 * t177 * t116 + 2.0 * t65 * t116 + t119 + 2.0 * t102 * t46 + t138;
                                    }

                                    if (Theta_D <= 0)
                                          then
                                          {
                                                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                                                           "\n"
                                                           "   compute_Theta(): Theta_D = $g^{ij} s_i s_j$ = %g <= 0\n"
                                                           "                    at %s patch rho=%g sigma=%g!\n"
                                                           "                    (i.e. the interpolated g_ij isn't positive definite)",
                                                           double(Theta_D),
                                                           p.name(), double(rho), double(sigma));

                                                cout << g_dd_11 << "," << g_dd_12 << "," << g_dd_13 << "," << g_dd_22 << "," << g_dd_23 << "," << g_dd_33 << endl;
                                                cout << xx << "," << yy << "," << zz << endl;
                                                return false; // *** ERROR RETURN ***
                                          }

                                    // compute H via equation (14) of my 1996 horizon finding paper
                                    const fp sqrt_Theta_D = sqrt(Theta_D);
                                    Theta = +Theta_A / (Theta_D * sqrt_Theta_D) + Theta_B / sqrt_Theta_D + Theta_C / Theta_D - K + add_to_expansion;

                                    // update running norms of Theta(h) function
                                    if (Theta_norms_ptr != NULL)
                                          then Theta_norms_ptr->data(Theta);

                                    if (Jacobian_flag)
                                          then
                                          {
                                                // partial_Theta_wrt_partial_d_h,
                                                // partial_Theta_wrt_partial_dd_h
                                                fp t1, t2, t3, t4, t5, t7, t8, t10, t11, t13;
                                                fp t14, t16, t18, t20, t22, t24, t26, t28, t29, t31;
                                                fp t32, t35, t37, t38, t41, t42, t43, t46, t48, t52;
                                                fp t54, t55, t59, t60, t63, t67, t68, t69, t70, t71;
                                                fp t74, t76, t78, t80, t83, t85, t86, t92, t93, t94;
                                                fp t98, t99, t102, t103, t104, t107, t108, t112, t113, t114;
                                                fp t115, t116, t118, t119, t120, t122, t123, t126, t127, t128;
                                                fp t133, t136, t140, t141, t142, t143, t153, t156, t158, t160;
                                                fp t162, t165, t167, t168, t171, t172, t173, t174, t179, t183;
                                                fp t185, t189, t190, t193, t194, t195, t197, t198, t202, t205;
                                                fp t208, t209, t212, t216, t217, t218, t220, t222, t223, t224;
                                                fp t226, t227, t232, t235, t236, t237, t238, t240, t247, t248;
                                                fp t249, t254, t259, t263, t266, t267, t275, t278, t281, t284;
                                                fp t287, t288, t291, t296, t297, t298, t300, t307, t309, t311;
                                                fp t314, t316, t317, t322, t325, t326, t329, t334, t335, t336;
                                                fp t340, t346, t350, t351, t352, t354, t357, t358, t359, t361;
                                                fp t364, t365, t366, t368, t370, t373, t374, t376, t381, t385;
                                                fp t386, t392, t398, t401, t404, t405, t407, t408, t411, t414;
                                                fp t416, t417, t419, t421, t422, t424, t428, t431, t432, t434;
                                                fp t437, t440, t442, t449, t454, t458, t461, t467, t470, t471;
                                                fp t474, t475, t481, t485, t489, t494, t498, t503, t504, t505;
                                                fp t507, t514, t518, t534, t536, t542, t545, t548, t551, t552;
                                                fp t559, t561, t562, t565, t569, t571, t572, t573, t575, t576;
                                                fp t588, t589, t590, t593, t594, t599, t601, t605, t608, t609;
                                                fp t612, t613, t627, t632, t633, t640, t644, t652, t656, t664;
                                                fp t669, t672, t677, t678, t680, t694, t704, t707, t712, t716;
                                                fp t723, t738, t741, t746, t748, t750, t774, t776, t780, t785;
                                                fp t787, t792, t796, t797, t799, t800, t802, t803, t805, t807;
                                                fp t809, t811, t813, t815, t817, t819, t822, t824, t827, t829;
                                                fp t832, t835, t837, t840, t843, t847, t860, t869, t871, t876;
                                                fp t882, t886, t890, t891, t897, t899, t900, t902, t904, t905;
                                                fp t907, t913, t920, t929, t930, t933, t938, t944, t947, t949;
                                                fp t962, t970, t971, t976, t979, t983, t996, t997, t1000, t1001;
                                                fp t1004, t1010, t1012, t1015, t1033, t1036, t1039, t1047, t1048, t1050;
                                                fp t1062, t1065, t1070, t1074, t1075, t1078, t1080, t1082, t1087, t1093;
                                                fp t1095, t1097, t1103, t1107, t1112, t1114, t1138, t1139, t1141, t1145;
                                                fp t1150, t1163, t1166, t1169, t1174, t1186, t1189, t1192, t1200, t1214;
                                                fp t1234, t1266, t1281, t1289, t1300, t1301, t1308, t1335, t1342, t1345;
                                                fp t1364, t1370, t1405, t1414, t1427, t1457, t1460, t1463, t1465, t1469;
                                                fp t1475, t1476, t1477, t1483, t1486, t1487, t1491, t1492, t1493, t1497;
                                                fp t1505, t1508, t1510, t1513, t1516, t1517, t1520, t1526, t1536, t1547;
                                                fp t1552, t1555, t1558, t1561, t1572, t1580, t1594, t1600, t1606, t1610;
                                                fp t1622, t1629, t1639, t1641, t1643, t1645, t1648, t1655, t1659, t1660;
                                                fp t1666, t1667, t1684, t1697, t1704, t1718, t1721, t1739, t1748, t1751;
                                                fp t1757, t1760, t1761, t1768, t1771, t1783, t1785, t1788, t1791, t1803;
                                                fp t1809, t1812, t1825;
                                                t1 = g_uu_13;
                                                t2 = X_ud_13;
                                                t3 = t1 * t2;
                                                t4 = g_uu_12;
                                                t5 = 1 / r;
                                                t7 = X_ud_11;
                                                t8 = PARTIAL_RHO(h);
                                                t10 = X_ud_21;
                                                t11 = PARTIAL_SIGMA(h);
                                                t13 = xx * t5 - t7 * t8 - t10 * t11;
                                                t14 = t4 * t13;
                                                t16 = r * r;
                                                t18 = 1 / t16 / r;
                                                t20 = X_udd_112;
                                                t22 = X_udd_212;
                                                t24 = X_ud_12;
                                                t26 = PARTIAL_RHO_RHO(h);
                                                t28 = t10 * t24;
                                                t29 = PARTIAL_RHO_SIGMA(h);
                                                t31 = X_ud_22;
                                                t32 = t7 * t31;
                                                t35 = PARTIAL_SIGMA_SIGMA(h);
                                                t37 = -xx * yy * t18 - t20 * t8 - t22 * t11 - t7 * t24 * t26 - t28 * t29 - t32 * t29 - t10 * t31 * t35;
                                                t38 = t14 * t37;
                                                t41 = g_uu_22;
                                                t42 = t41 * t24;
                                                t43 = t1 * t13;
                                                t46 = X_udd_123;
                                                t48 = X_udd_223;
                                                t52 = t31 * t2;
                                                t54 = X_ud_23;
                                                t55 = t24 * t54;
                                                t59 = -yy * zz * t18 - t46 * t8 - t48 * t11 - t24 * t2 * t26 - t52 * t29 - t55 * t29 - t31 * t54 * t35;
                                                t60 = t43 * t59;
                                                t63 = g_uu_23;
                                                t67 = yy * t5 - t24 * t8 - t31 * t11;
                                                t68 = t63 * t67;
                                                t69 = t1 * t7;
                                                t70 = xx * xx;
                                                t71 = yy * yy;
                                                t74 = X_udd_133;
                                                t76 = X_udd_233;
                                                t78 = t2 * t2;
                                                t80 = t54 * t2;
                                                t83 = t54 * t54;
                                                t85 = (t70 + t71) * t18 - t74 * t8 - t76 * t11 - t78 * t26 - 2.0 * t80 * t29 - t83 * t35;
                                                t86 = t69 * t85;
                                                t92 = zz * t5 - t2 * t8 - t54 * t11;
                                                t93 = t63 * t92;
                                                t94 = t4 * t67;
                                                t98 = t41 * t67;
                                                t99 = t69 * t59;
                                                t102 = g_uu_33;
                                                t103 = t102 * t92;
                                                t104 = t43 * t74;
                                                t107 = t1 * t92;
                                                t108 = t4 * t7;
                                                t112 = g_uu_11;
                                                t113 = t112 * t13;
                                                t114 = partial_d_g_uu_123;
                                                t115 = t114 * t2;
                                                t116 = t115 * t67;
                                                t118 = partial_d_g_uu_211;
                                                t119 = t118 * t13;
                                                t120 = t119 * t7;
                                                t122 = t63 * t2;
                                                t123 = t94 * t37;
                                                t126 = partial_d_g_uu_122;
                                                t127 = t126 * t67;
                                                t128 = t127 * t24;
                                                t133 = t98 * t37;
                                                t136 = X_udd_113;
                                                t140 = 2.0 * t3 * t38 + 2.0 * t42 * t60 + 2.0 * t68 * t86 + 2.0 * t93 * t94 * t20 + 2.0 * t98 * t99 + 2.0 * t103 * t104 + 2.0 * t107 * t108 * t37 + t113 * t116 + t93 * t120 + 2.0 * t122 * t123 + t113 * t128 + 2.0 * t107 * t14 * t20 + 2.0 * t3 * t133 + 2.0 * t107 * t68 * t136;
                                                t141 = partial_d_g_uu_311;
                                                t142 = t141 * t13;
                                                t143 = t142 * t7;
                                                t153 = zz * zz;
                                                t156 = X_udd_122;
                                                t158 = X_udd_222;
                                                t160 = t24 * t24;
                                                t162 = t31 * t24;
                                                t165 = t31 * t31;
                                                t167 = (t70 + t153) * t18 - t156 * t8 - t158 * t11 - t160 * t26 - 2.0 * t162 * t29 - t165 * t35;
                                                t168 = t108 * t167;
                                                t171 = t13 * t13;
                                                t172 = t112 * t171;
                                                t173 = partial_d_g_uu_112;
                                                t174 = t173 * t24;
                                                t179 = X_udd_213;
                                                t183 = t10 * t2;
                                                t185 = t7 * t54;
                                                t189 = -xx * zz * t18 - t136 * t8 - t179 * t11 - t7 * t2 * t26 - t183 * t29 - t185 * t29 - t10 * t54 * t35;
                                                t190 = t68 * t189;
                                                t193 = t112 * t7;
                                                t194 = t114 * t92;
                                                t195 = t194 * t67;
                                                t197 = t4 * t4;
                                                t198 = t197 * t67;
                                                t202 = t108 * t59;
                                                t205 = t193 * t37;
                                                t208 = t102 * t2;
                                                t209 = t14 * t59;
                                                t212 = t63 * t24;
                                                t216 = t63 * t63;
                                                t217 = t92 * t92;
                                                t218 = t216 * t217;
                                                t220 = t103 * t143 + 2.0 * t94 * t43 * t136 + 2.0 * t107 * t98 * t20 + 2.0 * t68 * t104 + 2.0 * t93 * t168 + t172 * t174 + 2.0 * t3 * t190 + t193 * t195 + 2.0 * t198 * t7 * t37 + 2.0 * t103 * t202 + 2.0 * t93 * t205 + 2.0 * t208 * t209 + 2.0 * t107 * t212 * t189 + t218 * t156;
                                                t222 = t1 * t1;
                                                t223 = t222 * t217;
                                                t224 = X_udd_111;
                                                t226 = t102 * t102;
                                                t227 = t226 * t217;
                                                t232 = t113 * t189;
                                                t235 = t67 * t67;
                                                t236 = t41 * t235;
                                                t237 = partial_d_g_uu_223;
                                                t238 = t237 * t2;
                                                t240 = t194 * t24;
                                                t247 = partial_d_g_uu_333;
                                                t248 = t247 * t92;
                                                t249 = t248 * t2;
                                                t254 = t113 * t136;
                                                t259 = t1 * t171;
                                                t263 = t193 * t189;
                                                t266 = t223 * t224 + t227 * t74 + 2.0 * t107 * t42 * t37 + 2.0 * t208 * t232 + t236 * t238 + t113 * t240 + 2.0 * t93 * t98 * t156 + 2.0 * t68 * t202 + t43 * t249 + 2.0 * t93 * t42 * t167 + 2.0 * t103 * t254 + 2.0 * t212 * t209 + 2.0 * t259 * t4 * t46 + 2.0 * t103 * t263;
                                                t267 = t98 * t167;
                                                t275 = t14 * t46;
                                                t278 = t43 * t46;
                                                t281 = t113 * t224;
                                                t284 = t113 * t37;
                                                t287 = t102 * t217;
                                                t288 = t63 * t46;
                                                t291 = t113 * t20;
                                                t296 = partial_d_g_uu_312;
                                                t297 = t296 * t67;
                                                t298 = t297 * t13;
                                                t300 = t222 * t92;
                                                t307 = X_udd_211;
                                                t309 = t7 * t7;
                                                t311 = t10 * t7;
                                                t314 = t10 * t10;
                                                t316 = (t71 + t153) * t18 - t224 * t8 - t307 * t11 - t309 * t26 - 2.0 * t311 * t29 - t314 * t35;
                                                t317 = t113 * t316;
                                                t322 = 2.0 * t122 * t267 + 2.0 * t94 * t69 * t189 + 4.0 * t43 * t263 + 2.0 * t103 * t275 + 2.0 * t98 * t278 + 2.0 * t107 * t281 + 2.0 * t122 * t284 + 2.0 * t287 * t288 + 2.0 * t93 * t291 + 2.0 * t68 * t275 + t208 * t298 + 2.0 * t300 * t7 * t189 + 2.0 * t3 * t317 + 2.0 * t103 * t86;
                                                t325 = t4 * t24;
                                                t326 = t325 * t189;
                                                t329 = t43 * t85;
                                                t334 = partial_d_g_uu_313;
                                                t335 = t334 * t92;
                                                t336 = t335 * t13;
                                                t340 = t335 * t7;
                                                t346 = t63 * t59;
                                                t350 = partial_d_g_uu_111;
                                                t351 = t350 * t13;
                                                t352 = t351 * t7;
                                                t354 = t193 * t316;
                                                t357 = partial_d_g_uu_113;
                                                t358 = t357 * t2;
                                                t359 = t358 * t13;
                                                t361 = t94 * t189;
                                                t364 = partial_d_g_uu_323;
                                                t365 = t364 * t2;
                                                t366 = t365 * t67;
                                                t368 = 2.0 * t103 * t326 + 2.0 * t208 * t329 + 2.0 * t212 * t329 + t212 * t336 + 4.0 * t68 * t326 +
                                                       t68 * t340 + 2.0 * t93 * t278 + 4.0 * t43 * t202 + 4.0 * t103 * t346 * t2 + t94 * t352 + 2.0 * t107 * t354 + t94 * t359 + 2.0 * t208 * t361 + t43 * t366;
                                                t370 = t41 * t59 * t24;
                                                t373 = t357 * t92;
                                                t374 = t373 * t13;
                                                t376 = t1 * t189;
                                                t381 = t63 * t235;
                                                t385 = partial_d_g_uu_133;
                                                t386 = t385 * t217;
                                                t392 = t4 * t20;
                                                t398 = t350 * t171;
                                                t401 = t118 * t171;
                                                t404 = t334 * t2;
                                                t405 = t404 * t13;
                                                t407 = t4 * t37;
                                                t408 = t407 * t24;
                                                t411 = t43 * t189;
                                                t414 = 4.0 * t68 * t370 + t325 * t374 + 4.0 * t103 * t376 * t2 + t98 * t120 + 2.0 * t381 * t41 * t46 +
                                                       RATIONAL(1.0, 2.0) * t193 * t386 + 2.0 * t381 * t4 * t136 + 2.0 * t236 * t392 + 2.0 * t259 * t112 * t136 +
                                                       RATIONAL(1.0, 2.0) * t3 * t398 + RATIONAL(1.0, 2.0) * t122 * t401 + t68 * t405 + 4.0 * t98 * t408 + 2.0 * t325 * t411;
                                                t416 = t364 * t92;
                                                t417 = t416 * t67;
                                                t419 = t297 * t7;
                                                t421 = t296 * t24;
                                                t422 = t421 * t13;
                                                t424 = t1 * t37;
                                                t428 = t94 * t316;
                                                t431 = t41 * t41;
                                                t432 = t431 * t235;
                                                t434 = t126 * t235;
                                                t437 = t247 * t217;
                                                t440 = t416 * t24;
                                                t442 = t373 * t7;
                                                t449 = t431 * t67;
                                                t454 = t69 * t417 + t103 * t419 + t103 * t422 + 4.0 * t93 * t424 * t2 + 2.0 * t3 * t428 + t432 * t156 + RATIONAL(1.0, 2.0) * t193 * t434 + RATIONAL(1.0, 2.0) * t69 * t437 + t43 * t440 + t94 * t442 + 2.0 * t300 * t13 * t136 + t381 * t296 * t7 + 2.0 * t449 * t167 * t24 + t259 * t421;
                                                t458 = t350 * t7;
                                                t461 = t4 * t235;
                                                t467 = t13 * t189;
                                                t470 = t237 * t92;
                                                t471 = t470 * t24;
                                                t474 = t385 * t92;
                                                t475 = t474 * t2;
                                                t481 = t13 * t37;
                                                t485 = t67 * t59;
                                                t489 = t238 * t67;
                                                t494 = RATIONAL(3.0, 2.0) * t259 * t141 * t7 + RATIONAL(3.0, 2.0) * t172 * t458 + t461 * t115 + 2.0 * t198 * t13 * t20 + 2.0 * t222 * t2 * t467 + 2.0 * t98 * t471 + t113 * t475 + 2.0 * t107 * t94 * t224 + 2.0 * t197 * t24 * t481 + 2.0 * t216 * t2 * t485 + t68 * t249 + t14 * t489 + t107 * t128 + 2.0 * t93 * t99;
                                                t498 = t470 * t67;
                                                t503 = partial_d_g_uu_233;
                                                t504 = t503 * t92;
                                                t505 = t504 * t2;
                                                t507 = t4 * t171;
                                                t514 = t216 * t92;
                                                t518 = t334 * t7;
                                                t534 = t108 * t498 + 2.0 * t103 * t94 * t136 + t14 * t505 + RATIONAL(3.0, 2.0) * t507 * t118 * t7 + 2.0 * t107 * t325 * t316 + 2.0 * t514 * t24 * t59 + t287 * t518 + t259 * t404 + RATIONAL(3.0, 2.0) * t461 * t126 * t24 + 2.0 * t514 * t67 * t46 + RATIONAL(1.0, 2.0) * t3 * t434 + 2.0 * t68 * t440 + t172 * t358 + 2.0 * t68 * t422;
                                                t536 = partial_d_g_uu_213;
                                                t542 = t98 * t59;
                                                t545 = t68 * t85;
                                                t548 = t216 * t235;
                                                t551 = t536 * t13;
                                                t552 = t551 * t2;
                                                t559 = t174 * t13;
                                                t561 = t536 * t92;
                                                t562 = t561 * t7;
                                                t565 = t226 * t92;
                                                t569 = t94 * t475 + t507 * t536 * t2 + 2.0 * t43 * t340 + t14 * t471 + 2.0 * t208 * t542 + 2.0 * t208 * t545 + t548 * t74 + t98 * t505 + 2.0 * t93 * t552 + 2.0 * t94 * t240 + 2.0 * t113 * t442 + t107 * t559 + 2.0 * t14 * t562 + 2.0 * t565 * t85 * t2;
                                                t571 = partial_d_g_uu_322;
                                                t572 = t571 * t67;
                                                t573 = t572 * t24;
                                                t575 = t173 * t67;
                                                t576 = t575 * t13;
                                                t588 = partial_d_g_uu_212;
                                                t589 = t588 * t24;
                                                t590 = t589 * t13;
                                                t593 = t588 * t67;
                                                t594 = t593 * t13;
                                                t599 = t575 * t7;
                                                t601 = t63 * t217;
                                                t605 = t141 * t171;
                                                t608 = t43 * t573 + t3 * t576 + 2.0 * t103 * t405 + 2.0 * t43 * t419 + t103 * t573 + 2.0 * t107 * t359 + 2.0 * t514 * t167 * t2 + t93 * t590 + t381 * t365 + t122 * t594 + 2.0 * t103 * t98 * t46 + t107 * t599 +
                                                       2.0 * t601 * t1 * t20 + RATIONAL(1.0, 2.0) * t208 * t605;
                                                t609 = t593 * t7;
                                                t612 = partial_d_g_uu_222;
                                                t613 = t612 * t24;
                                                t627 = t588 * t7;
                                                t632 = t612 * t67;
                                                t633 = t632 * t24;
                                                t640 = t216 * t67;
                                                t644 = 2.0 * t14 * t609 + RATIONAL(3.0, 2.0) * t236 * t613 + t93 * t609 + 2.0 * t113 * t599 +
                                                       RATIONAL(1.0, 2.0) * t42 * t401 + 2.0 * t107 * t116 + RATIONAL(1.0, 2.0) * t325 * t398 + 2.0 * t103 * t366 + t236 * t627 + 2.0 * t103 * t212 * t85 + t14 * t633 + 2.0 * t93 * t489 + RATIONAL(3.0, 2.0) * t381 * t571 * t24 + 2.0 * t640 * t85 * t24;
                                                t652 = t364 * t24;
                                                t656 = t1 * t217;
                                                t664 = t247 * t2;
                                                t669 = t1 * t136;
                                                t672 = t503 * t217;
                                                t677 = t112 * t112;
                                                t678 = t677 * t171;
                                                t680 = 4.0 * t14 * t205 + 2.0 * t103 * t68 * t74 + t287 * t652 + t461 * t173 * t7 + t656 * t114 * t24 + t601 * t237 * t24 + t507 * t589 + t601 * t536 * t7 + RATIONAL(3.0, 2.0) * t287 * t664 + RATIONAL(1.0, 2.0) * t212 * t437 + 2.0 * t287 * t669 + RATIONAL(1.0, 2.0) * t108 * t672 + RATIONAL(1.0, 2.0) * t42 * t672 + t678 * t224;
                                                t694 = t677 * t13;
                                                t704 = t571 * t235;
                                                t707 = t612 * t235;
                                                t712 = t222 * t13;
                                                t716 = 2.0 * t98 * t590 + 2.0 * t300 * t316 * t2 + 2.0 * t94 * t559 + t98 * t562 + 2.0 * t122 * t60 +
                                                       t93 * t633 + 2.0 * t103 * t370 + 2.0 * t694 * t316 * t7 + RATIONAL(3.0, 2.0) * t656 * t385 * t2 + RATIONAL(3.0, 2.0) * t601 * t503 * t2 + RATIONAL(1.0, 2.0) * t208 * t704 + RATIONAL(1.0, 2.0) * t122 * t707 +
                                                       RATIONAL(1.0, 2.0) * t69 * t704 + 2.0 * t712 * t85 * t7;
                                                t723 = t197 * t13;
                                                t738 = t14 * t167;
                                                t741 = t14 * t156;
                                                t746 = t561 * t13;
                                                t748 = t197 * t235;
                                                t750 = 2.0 * t198 * t316 * t24 + t656 * t357 * t7 + 2.0 * t723 * t167 * t7 + t68 * t143 + 2.0 * t507 * t112 * t20 + 2.0 * t94 * t354 + t98 * t552 + RATIONAL(1.0, 2.0) * t108 * t707 + RATIONAL(1.0, 2.0) * t212 * t605 + 2.0 * t122 * t738 + 2.0 * t98 * t741 + 2.0 * t93 * t408 + t42 * t746 + t748 * t224;
                                                t774 = t197 * t171;
                                                t776 = t222 * t171;
                                                t780 = 2.0 * t94 * t281 + 2.0 * t42 * t284 + 2.0 * t98 * t168 + t107 * t352 + 2.0 * t212 * t232 + 2.0 * t93 * t741 + RATIONAL(1.0, 2.0) * t325 * t386 + 2.0 * t42 * t738 + 2.0 * t98 * t205 + 2.0 * t98 * t291 +
                                                       2.0 * t325 * t317 + 2.0 * t68 * t254 + t774 * t156 + t776 * t74 + 2.0 * t68 * t263;
                                                t785 = pow(Theta_D, 1.0 * RATIONAL(1.0, 2.0));
                                                t787 = 1 / t785 / Theta_D;
                                                t792 = -t458 - t627 - t518 - t174 - t613 - t652 - t358 - t238 - t664 - t112 * t224 - 2.0 * t392 - 2.0 * t669;
                                                t796 = partial_d_ln_sqrt_g_1;
                                                t797 = t112 * t796;
                                                t799 = partial_d_ln_sqrt_g_2;
                                                t800 = t4 * t799;
                                                t802 = partial_d_ln_sqrt_g_3;
                                                t803 = t1 * t802;
                                                t805 = t4 * t796;
                                                t807 = t41 * t799;
                                                t809 = t63 * t802;
                                                t811 = t1 * t796;
                                                t813 = t63 * t799;
                                                t815 = t102 * t802;
                                                t817 = -t41 * t156 - 2.0 * t288 - t102 * t74 - t797 * t7 - t800 * t7 - t803 * t7 - t805 * t24 - t807 * t24 - t809 * t24 - t811 * t2 - t813 * t2 - t815 * t2;
                                                t819 = 1 / t785;
                                                t822 = K_uu_11 * t13;
                                                t824 = K_uu_12;
                                                t827 = t824 * t67;
                                                t829 = K_uu_13;
                                                t832 = t829 * t92;
                                                t835 = K_uu_22 * t67;
                                                t837 = K_uu_23;
                                                t840 = t837 * t92;
                                                t843 = K_uu_33 * t92;
                                                t847 = 1 / Theta_D;
                                                t860 = Theta_D * Theta_D;
                                                t869 = RATIONAL(3.0, 2.0) * Theta_A / t785 / t860 + RATIONAL(1.0, 2.0) * Theta_B * t787 + Theta_C / t860;
                                                partial_Theta_wrt_partial_d_h_1 = (t140 + t220 + t266 + t322 + t368 + t414 + t454 +
                                                                                   t494 + t534 + t569 + t608 + t644 + t680 + t716 + t750 + t780) *
                                                                                      t787 +
                                                                                  (t792 + t817) * t819 + (-2.0 * t822 * t7 - 2.0 * t824 * t24 * t13 - 2.0 * t827 * t7 - 2.0 * t829 * t2 * t13 - 2.0 * t832 * t7 - 2.0 * t835 * t24 - 2.0 * t837 * t2 * t67 - 2.0 * t840 * t24 - 2.0 * t843 * t2) * t847 - (-2.0 * t113 * t7 - 2.0 * t325 * t13 - 2.0 * t94 * t7 - 2.0 * t3 * t13 - 2.0 * t107 * t7 - 2.0 * t98 * t24 - 2.0 * t122 * t67 - 2.0 * t93 * t24 - 2.0 * t103 * t2) * t869;
                                                t871 = t113 * t22;
                                                t876 = t63 * t54;
                                                t882 = t551 * t54;
                                                t886 = t561 * t10;
                                                t890 = t112 * t10;
                                                t891 = t890 * t316;
                                                t897 = t334 * t10;
                                                t899 = 2.0 * t93 * t871 + t381 * t296 * t10 + t876 * t594 + 2.0 * t93 * t94 * t22 + t432 * t158 + 2.0 * t93 * t882 + t218 * t158 + 2.0 * t14 * t886 + t748 * t307 + 2.0 * t94 * t891 + t890 * t195 + t548 * t76 + t223 * t307 + t287 * t897;
                                                t900 = t194 * t31;
                                                t902 = t334 * t54;
                                                t904 = t114 * t54;
                                                t905 = t904 * t67;
                                                t907 = t63 * t31;
                                                t913 = t102 * t54;
                                                t920 = t14 * t48;
                                                t929 = t4 * t10;
                                                t930 = t929 * t59;
                                                t933 = t335 * t10;
                                                t938 = t113 * t900 + t259 * t902 + t113 * t905 + 2.0 * t907 * t209 + 2.0 * t300 * t13 * t179 + 2.0 * t913 * t545 + t507 * t536 * t54 + t601 * t536 * t10 + 2.0 * t68 * t920 + 2.0 * t712 * t85 * t10 + 2.0 * t449 * t167 * t31 + 2.0 * t68 * t930 + t68 * t933 + 2.0 * t197 * t31 * t481;
                                                t944 = t1 * t54;
                                                t947 = t588 * t31;
                                                t949 = t113 * t307;
                                                t962 = t364 * t54;
                                                t970 = t4 * t31;
                                                t971 = t970 * t189;
                                                t976 = t913 * t298 + 2.0 * t103 * t907 * t85 + 2.0 * t944 * t317 + t507 * t947 + 2.0 * t107 * t949 +
                                                       RATIONAL(3.0, 2.0) * t507 * t118 * t10 + 2.0 * t259 * t4 * t48 + t259 * t296 * t31 + 2.0 * t107 * t891 +
                                                       t381 * t962 + 2.0 * t198 * t13 * t22 + 2.0 * t103 * t68 * t76 + 2.0 * t103 * t971 + 2.0 * t876 * t60;
                                                t979 = t416 * t31;
                                                t983 = t351 * t10;
                                                t996 = t1 * t10;
                                                t997 = t996 * t85;
                                                t1000 = t41 * t31;
                                                t1001 = t1000 * t59;
                                                t1004 = t996 * t59;
                                                t1010 = t142 * t10;
                                                t1012 = 2.0 * t907 * t329 + t43 * t979 + 2.0 * t913 * t361 + t107 * t983 + 2.0 * t944 * t38 + 4.0 * t93 * t424 * t54 + 2.0 * t107 * t929 * t37 + 2.0 * t103 * t94 * t179 + 2.0 * t68 * t997 + 2.0 * t103 * t1001 +
                                                        2.0 * t98 * t1004 + t970 * t374 + 2.0 * t913 * t542 + t103 * t1010;
                                                t1015 = t119 * t10;
                                                t1033 = t43 * t48;
                                                t1036 = t297 * t10;
                                                t1039 = t373 * t10;
                                                t1047 = t357 * t54;
                                                t1048 = t1047 * t13;
                                                t1050 = t93 * t1015 + 2.0 * t107 * t1000 * t37 + 2.0 * t259 * t112 * t179 + 2.0 * t970 * t411 + 2.0 * t944 * t133 + 2.0 * t93 * t1004 + 2.0 * t103 * t98 * t48 + t774 * t158 + 2.0 * t98 * t1033 + 2.0 * t43 * t1036 + 2.0 * t113 * t1039 + 2.0 * t1000 * t60 + 2.0 * t94 * t43 * t179 + t94 * t1048;
                                                t1062 = t43 * t76;
                                                t1065 = t113 * t179;
                                                t1070 = t470 * t31;
                                                t1074 = t237 * t54;
                                                t1075 = t1074 * t67;
                                                t1078 = t504 * t54;
                                                t1080 = t474 * t54;
                                                t1082 = 2.0 * t107 * t98 * t22 + 2.0 * t94 * t996 * t189 + t98 * t1015 + 2.0 * t970 * t317 + 2.0 * t876 * t123 + 2.0 * t68 * t1062 + 2.0 * t68 * t1065 + 2.0 * t944 * t428 + t14 * t1070 + 2.0 * t94 * t949 + t14 * t1075 + t94 * t1039 + t14 * t1078 + t113 * t1080;
                                                t1087 = t112 * t189 * t10;
                                                t1093 = t248 * t54;
                                                t1095 = t127 * t31;
                                                t1097 = t572 * t31;
                                                t1103 = t296 * t13 * t31;
                                                t1107 = t407 * t31;
                                                t1112 = t632 * t31;
                                                t1114 = 4.0 * t43 * t930 + 4.0 * t43 * t1087 + t227 * t76 + 4.0 * t68 * t971 + t68 * t1093 + t107 * t1095 + t103 * t1097 + 4.0 * t68 * t1001 + t98 * t1078 + 2.0 * t68 * t1103 + t678 * t307 + 4.0 * t98 * t1107 +
                                                        RATIONAL(1.0, 2.0) * t1000 * t672 + t93 * t1112;
                                                t1138 = t173 * t31;
                                                t1139 = t1138 * t13;
                                                t1141 = t4 * t22;
                                                t1145 = 2.0 * t381 * t41 * t48 + 2.0 * t216 * t54 * t485 + t103 * t1103 + RATIONAL(1.0, 2.0) * t996 * t704 + t103 * t1036 + t601 * t237 * t31 + t172 * t1047 + RATIONAL(3.0, 2.0) * t601 * t503 * t54 +
                                                        2.0 * t514 * t31 * t59 + 2.0 * t514 * t67 * t48 + t776 * t76 + t107 * t1139 + 2.0 * t236 * t1141 + t996 * t417;
                                                t1150 = t593 * t10;
                                                t1163 = t947 * t13;
                                                t1166 = t962 * t67;
                                                t1169 = t575 * t10;
                                                t1174 = t944 * t576 + t93 * t1150 + 2.0 * t723 * t167 * t10 + RATIONAL(1.0, 2.0) * t890 * t386 + RATIONAL(1.0, 2.0) * t929 * t672 + RATIONAL(1.0, 2.0) * t944 * t434 + RATIONAL(1.0, 2.0) * t907 * t437 + t93 * t1163 + t98 * t882 + t43 * t1166 + t1000 * t746 + t107 * t1169 + t43 * t1097 + 2.0 * t107 * t1048;
                                                t1186 = t112 * t37 * t10;
                                                t1189 = t929 * t167;
                                                t1192 = t14 * t158;
                                                t1200 = t43 * t1093 + t907 * t336 + 2.0 * t98 * t1070 + t113 * t1095 + 2.0 * t113 * t1169 + t68 * t1010 + t98 * t886 + t94 * t1080 + 4.0 * t14 * t1186 + 2.0 * t93 * t1189 + 2.0 * t93 * t1192 + 2.0 * t913 * t232 + 2.0 * t876 * t738 + t14 * t1112;
                                                t1214 = t902 * t13;
                                                t1234 = 2.0 * t107 * t14 * t22 + 2.0 * t1000 * t738 + 2.0 * t876 * t267 + 2.0 * t107 * t68 * t179 +
                                                        2.0 * t94 * t900 + t68 * t1214 + 2.0 * t103 * t920 + 2.0 * t944 * t190 + 2.0 * t68 * t1087 + 2.0 * t93 * t98 * t158 + 2.0 * t907 * t232 + 2.0 * t93 * t1000 * t167 + 2.0 * t98 * t1192 + 2.0 * t98 * t1189;
                                                t1266 = 2.0 * t103 * t1065 + t94 * t983 + 2.0 * t1000 * t284 + 2.0 * t198 * t316 * t31 + 2.0 * t107 * t907 * t189 + RATIONAL(1.0, 2.0) * t890 * t434 + 2.0 * t103 * t1166 + 2.0 * t43 * t933 + 2.0 * t103 * t930 + 2.0 * t94 * t1139 + 2.0 * t98 * t1163 + 2.0 * t98 * t1186 + RATIONAL(3.0, 2.0) * t259 * t141 * t10 +
                                                        2.0 * t222 * t54 * t467;
                                                t1281 = t588 * t10;
                                                t1289 = t247 * t54;
                                                t1300 = 2.0 * t300 * t10 * t189 + RATIONAL(3.0, 2.0) * t656 * t385 * t54 + t461 * t904 + t172 * t1138 + t236 * t1074 + 2.0 * t694 * t316 * t10 + t236 * t1281 + RATIONAL(3.0, 2.0) * t381 * t571 * t31 +
                                                        RATIONAL(3.0, 2.0) * t461 * t126 * t31 + RATIONAL(3.0, 2.0) * t287 * t1289 + 2.0 * t103 * t1087 + 2.0 * t107 * t905 + 2.0 * t68 * t979 + 2.0 * t98 * t871;
                                                t1301 = t612 * t31;
                                                t1308 = t364 * t31;
                                                t1335 = RATIONAL(3.0, 2.0) * t236 * t1301 + 2.0 * t93 * t1075 + 2.0 * t14 * t1150 + t287 * t1308 + t656 * t114 * t31 + t461 * t173 * t10 + RATIONAL(1.0, 2.0) * t1000 * t401 + t656 * t357 * t10 +
                                                        2.0 * t300 * t316 * t54 + 2.0 * t507 * t112 * t22 + 2.0 * t640 * t85 * t31 + 2.0 * t601 * t1 * t22 + RATIONAL(1.0, 2.0) * t876 * t707 + 2.0 * t565 * t85 * t54;
                                                t1342 = t63 * t48;
                                                t1345 = t1 * t179;
                                                t1364 = t350 * t10;
                                                t1370 = RATIONAL(1.0, 2.0) * t913 * t704 + 2.0 * t514 * t167 * t54 + 2.0 * t287 * t1342 + 2.0 * t287 * t1345 + RATIONAL(1.0, 2.0) * t970 * t398 + RATIONAL(1.0, 2.0) * t996 * t437 + RATIONAL(1.0, 2.0) * t907 * t605 + RATIONAL(1.0, 2.0) * t944 * t398 + RATIONAL(1.0, 2.0) * t876 * t401 +
                                                        RATIONAL(1.0, 2.0) * t913 * t605 + RATIONAL(1.0, 2.0) * t929 * t707 + RATIONAL(1.0, 2.0) * t970 * t386 + RATIONAL(3.0, 2.0) * t172 * t1364 + 2.0 * t198 * t10 * t37;
                                                t1405 = 4.0 * t103 * t376 * t54 + 2.0 * t103 * t1214 + 4.0 * t103 * t346 * t54 + 2.0 * t93 * t1107 +
                                                        2.0 * t107 * t94 * t307 + 2.0 * t107 * t970 * t316 + 2.0 * t913 * t209 + 2.0 * t381 * t4 * t179 + t929 * t498 +
                                                        2.0 * t93 * t1033 + 2.0 * t103 * t997 + 2.0 * t103 * t1062 + 2.0 * t913 * t329 + 2.0 * t876 * t284 + 2.0 * t93 * t1186;
                                                t1414 = -t1364 - t1281 - t897 - t1138 - t1301 - t1308 - t1047 - t1074 - t1289 - t112 * t307 - 2.0 * t1141 - 2.0 * t1345;
                                                t1427 = -t41 * t158 - 2.0 * t1342 - t102 * t76 - t797 * t10 - t800 * t10 - t803 * t10 - t805 * t31 -
                                                        t807 * t31 - t809 * t31 - t811 * t54 - t813 * t54 - t815 * t54;
                                                partial_Theta_wrt_partial_d_h_2 = (t899 + t938 + t976 + t1012 + t1050 + t1082 + t1114 + t1145 + t1174 + t1200 + t1234 + t1266 + t1300 + t1335 + t1370 + t1405) * t787 + (t1414 + t1427) * t819 + (-2.0 * t822 * t10 - 2.0 * t824 * t31 * t13 - 2.0 * t827 * t10 - 2.0 * t829 * t54 * t13 - 2.0 * t832 * t10 - 2.0 * t835 * t31 - 2.0 * t837 * t54 * t67 - 2.0 * t840 * t31 - 2.0 * t843 * t54) * t847 - (-2.0 * t113 * t10 - 2.0 * t970 * t13 - 2.0 * t94 * t10 - 2.0 * t944 * t13 - 2.0 * t107 * t10 - 2.0 * t98 * t31 - 2.0 * t876 * t67 - 2.0 * t93 * t31 - 2.0 * t103 * t54) * t869;
                                                t1457 = t14 * t160;
                                                t1460 = t69 * t2;
                                                t1463 = t68 * t4;
                                                t1465 = t13 * t24 * t2;
                                                t1469 = t43 * t78;
                                                t1475 = t103 * t4;
                                                t1476 = t67 * t7;
                                                t1477 = t1476 * t2;
                                                t1483 = t212 * t2;
                                                t1486 = t107 * t41;
                                                t1487 = t1476 * t24;
                                                t1491 = 2.0 * t98 * t1457 + 2.0 * t287 * t1460 + 2.0 * t1463 * t1465 + t774 * t160 + 2.0 * t68 * t1469 + 2.0 * t107 * t94 * t309 + 2.0 * t1475 * t1477 + 2.0 * t381 * t108 * t2 + 2.0 * t287 * t1483 + 2.0 * t1486 * t1487 + t218 * t160;
                                                t1492 = t13 * t7;
                                                t1493 = t1492 * t24;
                                                t1497 = t113 * t309;
                                                t1505 = t98 * t1;
                                                t1508 = t103 * t41;
                                                t1510 = t67 * t24 * t2;
                                                t1513 = t93 * t4;
                                                t1516 = t94 * t1;
                                                t1517 = t1492 * t2;
                                                t1520 = t107 * t4;
                                                t1526 = 2.0 * t198 * t1493 + t223 * t309 + 2.0 * t107 * t1497 + 2.0 * t259 * t193 * t2 + 2.0 * t93 * t1457 + 2.0 * t1505 * t1465 + 2.0 * t1508 * t1510 + 2.0 * t1513 * t1487 + 2.0 * t1516 * t1517 + 2.0 * t1520 * t1493 + 2.0 * t103 * t68 * t78;
                                                t1536 = t93 * t112;
                                                t1547 = t93 * t1;
                                                t1552 = t432 * t160 + 2.0 * t93 * t98 * t160 + t548 * t78 + 2.0 * t514 * t1510 + t748 * t309 + 2.0 * t1536 * t1493 + 2.0 * t300 * t1517 + 2.0 * t381 * t42 * t2 + 2.0 * t507 * t193 * t24 + 2.0 * t1547 * t1465 +
                                                        2.0 * t1475 * t1465;
                                                t1555 = t103 * t112;
                                                t1558 = t98 * t112;
                                                t1561 = t108 * t24;
                                                t1572 = t68 * t112;
                                                t1580 = 2.0 * t1547 * t1477 + 2.0 * t1555 * t1517 + 2.0 * t1558 * t1493 + 2.0 * t236 * t1561 +
                                                        2.0 * t259 * t325 * t2 + t776 * t78 + t678 * t309 + t227 * t78 + 2.0 * t94 * t1497 + 2.0 * t1572 * t1517 + 2.0 * t103 * t1469 + 2.0 * t601 * t69 * t24;
                                                partial_Theta_wrt_partial_dd_h_11 = (t1491 + t1526 + t1552 + t1580) * t787 + (-t112 * t309 - 2.0 * t1561 - 2.0 * t1460 - t41 * t160 - 2.0 * t1483 - t102 * t78) * t819;
                                                t1594 = -t183 - t185;
                                                t1600 = t67 * t10;
                                                t1606 = -t28 - t32;
                                                t1610 = t1 * t1594;
                                                t1622 = 2.0 * t218 * t162 - 2.0 * t107 * t68 * t1594 + 2.0 * t432 * t162 + 4.0 * t1520 * t1600 * t7 + 2.0 * t776 * t80 - 2.0 * t601 * t1 * t1606 - 2.0 * t287 * t1610 + 2.0 * t223 * t311 + 2.0 * t748 * t311 - 2.0 * t93 * t94 * t1606 + 2.0 * t774 * t162;
                                                t1629 = -t52 - t55;
                                                t1639 = t113 * t1606;
                                                t1641 = t113 * t1594;
                                                t1643 = t14 * t1629;
                                                t1645 = -t381 * t4 * t1594 - t300 * t13 * t1594 - t107 * t98 * t1606 - t259 * t4 * t1629 - t103 * t94 * t1594 - t107 * t14 * t1606 + t678 * t311 - t507 * t112 * t1606 - t93 * t1639 - t103 * t1641 - t103 * t1643;
                                                t1648 = t43 * t1629;
                                                t1655 = t13 * t54 * t2;
                                                t1659 = t13 * t10;
                                                t1660 = t1659 * t7;
                                                t1666 = t13 * t31;
                                                t1667 = t1666 * t24;
                                                t1684 = -2.0 * t93 * t1648 + 2.0 * t227 * t80 + 4.0 * t103 * t1 * t1655 + 4.0 * t94 * t112 * t1660 - 2.0 * t198 * t13 * t1606 + 4.0 * t1513 * t1667 - 2.0 * t514 * t67 * t1629 - 2.0 * t94 * t43 * t1594 + 4.0 * t68 * t1 * t1655 + 4.0 * t98 * t4 * t1667 - 2.0 * t98 * t1639;
                                                t1697 = t4 * t1606;
                                                t1704 = t67 * t31;
                                                t1718 = t63 * t1629;
                                                t1721 = -2.0 * t259 * t112 * t1594 - 2.0 * t68 * t1643 - 2.0 * t68 * t1641 - 2.0 * t98 * t1648 +
                                                        4.0 * t107 * t112 * t1660 - 2.0 * t236 * t1697 - 2.0 * t381 * t41 * t1629 + 4.0 * t93 * t41 * t1704 * t24 - 2.0 * t103 * t98 * t1629 + 2.0 * t548 * t80 + 4.0 * t103 * t63 * t67 * t54 * t2 - 2.0 * t287 * t1718;
                                                partial_Theta_wrt_partial_dd_h_12 = (t1622 + 2.0 * t1645 + t1684 + t1721) * t787 + (-2.0 * t890 * t7 + 2.0 * t1697 + 2.0 * t1610 - 2.0 * t1000 * t24 + 2.0 * t1718 - 2.0 * t913 * t2) * t819;
                                                t1739 = t996 * t54;
                                                t1748 = t1704 * t54;
                                                t1751 = t1600 * t54;
                                                t1757 = t1600 * t31;
                                                t1760 = 2.0 * t507 * t890 * t31 + 2.0 * t93 * t98 * t165 + t227 * t83 + t548 * t83 + 2.0 * t287 * t1739 + 2.0 * t259 * t970 * t54 + 2.0 * t601 * t996 * t31 + 2.0 * t514 * t1748 + 2.0 * t1547 * t1751 + 2.0 * t107 * t94 * t314 + 2.0 * t1486 * t1757;
                                                t1761 = t907 * t54;
                                                t1768 = t1659 * t31;
                                                t1771 = t113 * t314;
                                                t1783 = 2.0 * t287 * t1761 + t748 * t314 + t774 * t165 + t678 * t314 + t223 * t314 + 2.0 * t198 * t1768 + 2.0 * t94 * t1771 + 2.0 * t1513 * t1757 + 2.0 * t1520 * t1768 + 2.0 * t1475 * t1751 + 2.0 * t103 * t68 * t83;
                                                t1785 = t1666 * t54;
                                                t1788 = t1659 * t54;
                                                t1791 = t43 * t83;
                                                t1803 = t14 * t165;
                                                t1809 = 2.0 * t1547 * t1785 + 2.0 * t1516 * t1788 + 2.0 * t68 * t1791 + 2.0 * t1558 * t1768 + 2.0 * t259 * t890 * t54 + 2.0 * t107 * t1771 + 2.0 * t1463 * t1785 + 2.0 * t98 * t1803 + t218 * t165 + t776 * t83 +
                                                        t432 * t165;
                                                t1812 = t929 * t31;
                                                t1825 = t1572 * t1788 + t1505 * t1785 + t236 * t1812 + t103 * t1791 + t300 * t1788 + t381 * t1000 * t54 + t381 * t929 * t54 + t93 * t1803 + t1555 * t1788 + t1536 * t1768 + t1475 * t1785 + t1508 * t1748;
                                                partial_Theta_wrt_partial_dd_h_22 = (t1760 + t1783 + t1809 + 2.0 * t1825) * t787 + (-t112 * t314 - 2.0 * t1812 - 2.0 * t1739 - t41 * t165 - 2.0 * t1761 - t102 * t83) * t819;
                                          }
                              }
                        }
                  }

                  return true; // *** NORMAL RETURN ***
            }
      }

} // namespace AHFinderDirect
#endif
