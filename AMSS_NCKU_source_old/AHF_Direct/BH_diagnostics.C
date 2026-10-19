#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "util_Table.h"
#include "cctk.h"

#include "config.h"
#include "stdc.h"
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
#include "myglobal.h"

#include "horizon_sequence.h"
#include "BH_diagnostics.h"
#include "driver.h"

namespace AHFinderDirect
{
	using jtutil::error_exit;

	BH_diagnostics::BH_diagnostics()
		: centroid_x(0.0), centroid_y(0.0), centroid_z(0.0),
		  quadrupole_xx(0.0), quadrupole_xy(0.0), quadrupole_xz(0.0),
		  quadrupole_yy(0.0), quadrupole_yz(0.0),
		  quadrupole_zz(0.0),
		  min_radius(0.0), max_radius(0.0),
		  mean_radius(0.0),
		  min_x(0.0), max_x(0.0),
		  min_y(0.0), max_y(0.0),
		  min_z(0.0), max_z(0.0),
		  circumference_xy(0.0), circumference_xz(0.0), circumference_yz(0.0),
		  area(0.0), irreducible_mass(0.0), areal_radius(0.0) // no comma
	{
	}

	void BH_diagnostics::copy_to_buffer(double buffer[N_buffer])
		const
	{
		buffer[posn__centroid_x] = centroid_x;
		buffer[posn__centroid_y] = centroid_y;
		buffer[posn__centroid_z] = centroid_z;

		buffer[posn__quadrupole_xx] = quadrupole_xx;
		buffer[posn__quadrupole_xy] = quadrupole_xy;
		buffer[posn__quadrupole_xz] = quadrupole_xz;
		buffer[posn__quadrupole_yy] = quadrupole_yy;
		buffer[posn__quadrupole_xz] = quadrupole_yz;
		buffer[posn__quadrupole_zz] = quadrupole_zz;

		buffer[posn__min_radius] = min_radius;
		buffer[posn__max_radius] = max_radius;
		buffer[posn__mean_radius] = mean_radius;

		buffer[posn__min_x] = min_x;
		buffer[posn__max_x] = max_x;
		buffer[posn__min_y] = min_y;
		buffer[posn__max_y] = max_y;
		buffer[posn__min_z] = min_z;
		buffer[posn__max_z] = max_z;

		buffer[posn__circumference_xy] = circumference_xy;
		buffer[posn__circumference_xz] = circumference_xz;
		buffer[posn__circumference_yz] = circumference_yz;

		buffer[posn__area] = area;
		buffer[posn__irreducible_mass] = irreducible_mass;
		buffer[posn__areal_radius] = areal_radius;
	}

	void BH_diagnostics::copy_from_buffer(const double buffer[N_buffer])
	{
		centroid_x = buffer[posn__centroid_x];
		centroid_y = buffer[posn__centroid_y];
		centroid_z = buffer[posn__centroid_z];

		quadrupole_xx = buffer[posn__quadrupole_xx];
		quadrupole_xy = buffer[posn__quadrupole_xy];
		quadrupole_xz = buffer[posn__quadrupole_xz];
		quadrupole_yy = buffer[posn__quadrupole_yy];
		quadrupole_yz = buffer[posn__quadrupole_yz];
		quadrupole_zz = buffer[posn__quadrupole_zz];

		min_radius = buffer[posn__min_radius];
		max_radius = buffer[posn__max_radius];
		mean_radius = buffer[posn__mean_radius];

		min_x = buffer[posn__min_x];
		max_x = buffer[posn__max_x];
		min_y = buffer[posn__min_y];
		max_y = buffer[posn__max_y];
		min_z = buffer[posn__min_z];
		max_z = buffer[posn__max_z];

		circumference_xy = buffer[posn__circumference_xy];
		circumference_xz = buffer[posn__circumference_xz];
		circumference_yz = buffer[posn__circumference_yz];

		area = buffer[posn__area];
		irreducible_mass = buffer[posn__irreducible_mass];
		areal_radius = buffer[posn__areal_radius];
	}
	void BH_diagnostics::compute(patch_system &ps)
	{
		jtutil::norm<fp> h_norms;
		ps.ghosted_gridfn_norms(gfns::gfn__h, h_norms);
		min_radius = h_norms.min_abs_value();
		max_radius = h_norms.max_abs_value();

		jtutil::norm<fp> x_norms;
		jtutil::norm<fp> y_norms;
		jtutil::norm<fp> z_norms;

		ps.gridfn_norms(gfns::gfn__global_x, x_norms);
		ps.gridfn_norms(gfns::gfn__global_y, y_norms);
		ps.gridfn_norms(gfns::gfn__global_z, z_norms);

		min_x = x_norms.min_value();
		max_x = x_norms.max_value();
		min_y = y_norms.min_value();
		max_y = y_norms.max_value();
		min_z = z_norms.min_value();
		max_z = z_norms.max_value();

// adjust the bounding box for the symmetries
#define REFLECT(origin_, max_) (origin_ - (max_ - origin_))
		switch (ps.type())
		{
		case patch_system::patch_system__full_sphere:
			break;
		case patch_system::patch_system__plus_z_hemisphere:
			min_z = REFLECT(ps.origin_z(), max_z);
			break;
		case patch_system::patch_system__plus_xy_quadrant_mirrored:
		case patch_system::patch_system__plus_xy_quadrant_rotating:
			min_x = REFLECT(ps.origin_x(), max_x);
			min_y = REFLECT(ps.origin_y(), max_y);
			break;
		case patch_system::patch_system__plus_xz_quadrant_mirrored:
		case patch_system::patch_system__plus_xz_quadrant_rotating:
			min_x = REFLECT(ps.origin_x(), max_x);
			min_z = REFLECT(ps.origin_z(), max_z);
			break;
		case patch_system::patch_system__plus_xyz_octant_mirrored:
		case patch_system::patch_system__plus_xyz_octant_rotating:
			min_x = REFLECT(ps.origin_x(), max_x);
			min_y = REFLECT(ps.origin_y(), max_y);
			min_z = REFLECT(ps.origin_z(), max_z);
			break;
		default:
			error_exit(PANIC_EXIT,
					   "***** BH_diagnostics::compute(): unknown patch system type()=(int)%d!\n"
					   "                                 (this should never happen!)\n",
					   int(ps.type())); /*NOTREACHED*/
		}

		//
		// surface integrals
		//
		const fp integral_one = surface_integral(ps,
												 gfns::gfn__one, true, true, true,
												 patch::integration_method__automatic_choice);
		const fp integral_h = surface_integral(ps,
											   gfns::gfn__h, true, true, true,
											   patch::integration_method__automatic_choice);
		const fp integral_x = surface_integral(ps,
											   gfns::gfn__global_x, true, true, false,
											   patch::integration_method__automatic_choice);
		const fp integral_y = surface_integral(ps,
											   gfns::gfn__global_y, true, false, true,
											   patch::integration_method__automatic_choice);
		const fp integral_z = surface_integral(ps,
											   gfns::gfn__global_z, false, true, true,
											   patch::integration_method__automatic_choice);
		const fp integral_xx = surface_integral(ps,
												gfns::gfn__global_xx, true, true, true,
												patch::integration_method__automatic_choice);
		const fp integral_xy = surface_integral(ps,
												gfns::gfn__global_xy, true, false, false,
												patch::integration_method__automatic_choice);
		const fp integral_xz = surface_integral(ps,
												gfns::gfn__global_xz, false, true, false,
												patch::integration_method__automatic_choice);
		const fp integral_yy = surface_integral(ps,
												gfns::gfn__global_yy, true, true, true,
												patch::integration_method__automatic_choice);
		const fp integral_yz = surface_integral(ps,
												gfns::gfn__global_yz, false, false, true,
												patch::integration_method__automatic_choice);
		const fp integral_zz = surface_integral(ps,
												gfns::gfn__global_zz, true, true, true,
												patch::integration_method__automatic_choice);

		//
		// centroids
		//
		centroid_x = integral_x / integral_one;
		centroid_y = integral_y / integral_one;
		centroid_z = integral_z / integral_one;

		//
		// quadrupoles (taken about centroid position)
		//
		quadrupole_xx = integral_xx / integral_one - centroid_x * centroid_x;
		quadrupole_xy = integral_xy / integral_one - centroid_x * centroid_y;
		quadrupole_xz = integral_xz / integral_one - centroid_x * centroid_z;
		quadrupole_yy = integral_yy / integral_one - centroid_y * centroid_y;
		quadrupole_yz = integral_yz / integral_one - centroid_y * centroid_z;
		quadrupole_zz = integral_zz / integral_one - centroid_z * centroid_z;

		//
		// mean radius of horizon
		//
		mean_radius = integral_h / integral_one;

		//
		// surface area and quantities derived from it
		//
		area = integral_one;
		irreducible_mass = sqrt(area / (16.0 * PI));
		areal_radius = sqrt(area / (4.0 * PI));

		//
		// proper circumferences
		//
		circumference_xy = ps.circumference("xy", gfns::gfn__h,
											gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
											gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
											gfns::gfn__g_dd_33,
											patch::integration_method__automatic_choice);
		circumference_xz = ps.circumference("xz", gfns::gfn__h,
											gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
											gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
											gfns::gfn__g_dd_33,
											patch::integration_method__automatic_choice);
		circumference_yz = ps.circumference("yz", gfns::gfn__h,
											gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
											gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
											gfns::gfn__g_dd_33,
											patch::integration_method__automatic_choice);

		// prepare P^i,S^i in xx,xy,xz and yy,yz,zz
		{
			for (int pn = 0; pn < ps.N_patches(); ++pn)
			{
				patch &p = ps.ith_patch(pn);

				for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
				{
					for (int isigma = p.min_isigma();
						 isigma <= p.max_isigma();
						 ++isigma)
					{
						const fp g_xx = p.gridfn(gfns::gfn__g_dd_11, irho, isigma);
						const fp g_xy = p.gridfn(gfns::gfn__g_dd_12, irho, isigma);
						const fp g_xz = p.gridfn(gfns::gfn__g_dd_13, irho, isigma);
						const fp g_yy = p.gridfn(gfns::gfn__g_dd_22, irho, isigma);
						const fp g_yz = p.gridfn(gfns::gfn__g_dd_23, irho, isigma);
						const fp g_zz = p.gridfn(gfns::gfn__g_dd_33, irho, isigma);

						const fp k_xx = p.gridfn(gfns::gfn__K_dd_11, irho, isigma);
						const fp k_xy = p.gridfn(gfns::gfn__K_dd_12, irho, isigma);
						const fp k_xz = p.gridfn(gfns::gfn__K_dd_13, irho, isigma);
						const fp k_yy = p.gridfn(gfns::gfn__K_dd_22, irho, isigma);
						const fp k_yz = p.gridfn(gfns::gfn__K_dd_23, irho, isigma);
						const fp k_zz = p.gridfn(gfns::gfn__K_dd_33, irho, isigma);
						const fp trk = p.gridfn(gfns::gfn__trK, irho, isigma);

						const fp r = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
						const fp rho = p.rho_of_irho(irho);
						const fp sigma = p.sigma_of_isigma(isigma);
						fp xx, yy, zz; // local Cardesian coordinate
						p.xyz_of_r_rho_sigma(r, rho, sigma, xx, yy, zz);
						const fp X_ud_11 = p.partial_rho_wrt_x(xx, yy, zz);
						const fp X_ud_12 = p.partial_rho_wrt_y(xx, yy, zz);
						const fp X_ud_13 = p.partial_rho_wrt_z(xx, yy, zz);
						const fp X_ud_21 = p.partial_sigma_wrt_x(xx, yy, zz);
						const fp X_ud_22 = p.partial_sigma_wrt_y(xx, yy, zz);
						const fp X_ud_23 = p.partial_sigma_wrt_z(xx, yy, zz);
#if 0 // for P^i and S^i
	  // F,i = x^i/r-X_ud_1i(dh/drho)-X_ud_2i(dh/dsigma)
		  double nx,ny,nz;
		  nx = xx/r-X_ud_11*p.partial_rho(gfns::gfn__h, irho,isigma)-X_ud_21*p.partial_sigma(gfns::gfn__h, irho,isigma);
		  ny = yy/r-X_ud_12*p.partial_rho(gfns::gfn__h, irho,isigma)-X_ud_22*p.partial_sigma(gfns::gfn__h, irho,isigma);
		  nz = zz/r-X_ud_13*p.partial_rho(gfns::gfn__h, irho,isigma)-X_ud_23*p.partial_sigma(gfns::gfn__h, irho,isigma);
		  double eps; // volume element
		  fp g_uu_11, g_uu_12, g_uu_13, g_uu_22, g_uu_23, g_uu_33;
		  double pxx,pxy,pxz,pyy,pyz,pzz;
  		    {
		    fp t1, t2, t4, t5, t7, t8, t11, t12, t14, t15;
		    fp t18, t21;
	      	    t1 = g_yy;
	      	    t2 = g_zz;
	      	    t4 = g_yz;
	      	    t5 = t4*t4;
	      	    t7 = g_xx;
	      	    t8 = t7*t1;
	      	    t11 = g_xy;
	      	    t12 = t11*t11;
	      	    t14 = g_xz;
	      	    t15 = t11*t14;
	      	    t18 = t14*t14;
		    eps = t8*t2-t7*t5-t12*t2+2.0*t15*t4-t18*t1;
	      	    t21 = 1/eps;
		    eps = sqrt(eps);
	      	    g_uu_11 = (t1*t2-t5)*t21;
	      	    g_uu_12 = -(t11*t2-t14*t4)*t21;
	      	    g_uu_13 = -(-t11*t4+t14*t1)*t21;
	      	    g_uu_22 = (t7*t2-t18)*t21;
	      	    g_uu_23 = -(t7*t4-t15)*t21;
	      	    g_uu_33 = (t8-t12)*t21;

		    t5 = g_uu_11*nx*nx+g_uu_22*ny*ny+g_uu_33*nz*nz+2*(g_uu_12*nx*ny+g_uu_13*nx*nz+g_uu_23*ny*nz);
		    t5 = sqrt(t5);
		    nx = nx/t5;  // lower index
		    ny = ny/t5;
		    nz = nz/t5;

		    pxx= g_uu_11*(g_uu_11*k_xx+g_uu_12*k_xy+g_uu_13*k_xz)
			+g_uu_12*(g_uu_11*k_xy+g_uu_12*k_yy+g_uu_13*k_yz)
			+g_uu_13*(g_uu_11*k_xz+g_uu_12*k_yz+g_uu_13*k_zz); //k^xx
		    pxy= g_uu_11*(g_uu_12*k_xx+g_uu_22*k_xy+g_uu_23*k_xz)
			+g_uu_12*(g_uu_12*k_xy+g_uu_22*k_yy+g_uu_23*k_yz)
			+g_uu_13*(g_uu_12*k_xz+g_uu_22*k_yz+g_uu_23*k_zz); //k^xy
		    pxz= g_uu_11*(g_uu_13*k_xx+g_uu_23*k_xy+g_uu_33*k_xz)
			+g_uu_12*(g_uu_13*k_xy+g_uu_23*k_yy+g_uu_33*k_yz)
			+g_uu_13*(g_uu_13*k_xz+g_uu_23*k_yz+g_uu_33*k_zz); //k^xz
		    pyy= g_uu_12*(g_uu_12*k_xx+g_uu_22*k_xy+g_uu_23*k_xz)
			+g_uu_22*(g_uu_12*k_xy+g_uu_22*k_yy+g_uu_23*k_yz)
			+g_uu_23*(g_uu_12*k_xz+g_uu_22*k_yz+g_uu_23*k_zz); //k^yy
		    pyz= g_uu_12*(g_uu_13*k_xx+g_uu_23*k_xy+g_uu_33*k_xz)
			+g_uu_22*(g_uu_13*k_xy+g_uu_23*k_yy+g_uu_33*k_yz)
			+g_uu_23*(g_uu_13*k_xz+g_uu_23*k_yz+g_uu_33*k_zz); //k^yz
		    pzz= g_uu_13*(g_uu_13*k_xx+g_uu_23*k_xy+g_uu_33*k_xz)
			+g_uu_23*(g_uu_13*k_xy+g_uu_23*k_yy+g_uu_33*k_yz)
			+g_uu_33*(g_uu_13*k_xz+g_uu_23*k_yz+g_uu_33*k_zz); //k^zz
		  }

		  pxx = pxx-g_uu_11*trk; // tracefree
		  pyy = pyy-g_uu_22*trk;
		  pzz = pzz-g_uu_33*trk;
		  double tx,ty,tz;
		  double sxx,sxy,sxz,syx,syy,syz,szx,szy,szz;
		  tx = nx*pxx + ny*pxy + nz*pxz;
		  ty = nx*pxy + ny*pyy + nz*pyz;
		  tz = nx*pxz + ny*pyz + nz*pzz;
		  sxx = xx*tx;
		  sxy = xx*ty;
		  sxz = xx*tz;
		  syx = yy*tx;
		  syy = yy*ty;
		  syz = yy*tz;
		  szx = zz*tx;
	          szy = zz*ty;
		  szz = zz*tz;
                  p.gridfn(gfns::gfn__global_xx, irho,isigma) = tx; //p^x
                  p.gridfn(gfns::gfn__global_xy, irho,isigma) = ty; //p^y
                  p.gridfn(gfns::gfn__global_xz, irho,isigma) = tz; //p^z
		  tx = eps*(syz-szy); //s_x
		  ty = eps*(szx-sxz);
		  tz = eps*(sxy-syx);
                  p.gridfn(gfns::gfn__global_yy, irho,isigma) = g_uu_11*tx+g_uu_12*ty+g_uu_13*tz; //s^x
                  p.gridfn(gfns::gfn__global_yz, irho,isigma) = g_uu_12*tx+g_uu_22*ty+g_uu_23*tz; //s^y
                  p.gridfn(gfns::gfn__global_zz, irho,isigma) = g_uu_13*tx+g_uu_23*ty+g_uu_33*tz; //s^z
#endif
#if 1 // for P_i and S_i
	  // F,i = x^i/r-X_ud_1i(dh/drho)-X_ud_2i(dh/dsigma)
						double nx, ny, nz;
						nx = xx / r - X_ud_11 * p.partial_rho(gfns::gfn__h, irho, isigma) - X_ud_21 * p.partial_sigma(gfns::gfn__h, irho, isigma);
						ny = yy / r - X_ud_12 * p.partial_rho(gfns::gfn__h, irho, isigma) - X_ud_22 * p.partial_sigma(gfns::gfn__h, irho, isigma);
						nz = zz / r - X_ud_13 * p.partial_rho(gfns::gfn__h, irho, isigma) - X_ud_23 * p.partial_sigma(gfns::gfn__h, irho, isigma);
						{
							fp g_uu_11, g_uu_12, g_uu_13, g_uu_22, g_uu_23, g_uu_33;
							fp t1, t2, t4, t5, t7, t8, t11, t12, t14, t15;
							fp t18, t21;
							t1 = g_yy;
							t2 = g_zz;
							t4 = g_yz;
							t5 = t4 * t4;
							t7 = g_xx;
							t8 = t7 * t1;
							t11 = g_xy;
							t12 = t11 * t11;
							t14 = g_xz;
							t15 = t11 * t14;
							t18 = t14 * t14;
							t21 = 1 / (t8 * t2 - t7 * t5 - t12 * t2 + 2.0 * t15 * t4 - t18 * t1);
							g_uu_11 = (t1 * t2 - t5) * t21;
							g_uu_12 = -(t11 * t2 - t14 * t4) * t21;
							g_uu_13 = -(-t11 * t4 + t14 * t1) * t21;
							g_uu_22 = (t7 * t2 - t18) * t21;
							g_uu_23 = -(t7 * t4 - t15) * t21;
							g_uu_33 = (t8 - t12) * t21;

							t1 = g_uu_11 * nx + g_uu_12 * ny + g_uu_13 * nz;
							t2 = g_uu_12 * nx + g_uu_22 * ny + g_uu_23 * nz;
							t4 = g_uu_13 * nx + g_uu_23 * ny + g_uu_33 * nz;
							t5 = g_uu_11 * nx * nx + g_uu_22 * ny * ny + g_uu_33 * nz * nz + 2 * (g_uu_12 * nx * ny + g_uu_13 * nx * nz + g_uu_23 * ny * nz);
							t5 = sqrt(t5);
							nx = t1 / t5; // uper index
							ny = t2 / t5;
							nz = t4 / t5;
						}

						double pxx, pxy, pxz, pyy, pyz, pzz;
						double sxx, sxy, sxz, syx, syy, syz, szx, szy, szz;
						// these tensor components are same for local Cardisean and global Cardisean
						pxx = k_xx - g_xx * trk; // lower index
						pxy = k_xy;
						pxz = k_xz;
						pyy = k_yy - g_yy * trk;
						pyz = k_yz;
						pzz = k_zz - g_zz * trk;
						/*
								  sxx = yy*pxy - zz*pxz;
								  sxy = yy*pyy - zz*pyz;
								  sxz = yy*pyz - zz*pzz;
								  syx = zz*pxy - yy*pxz;
								  syy = zz*pyy - yy*pyz;
								  syz = zz*pyz - yy*pzz;
								  szx = xx*pxy - yy*pxx;
								  szy = xx*pyy - yy*pxy;
								  szz = xx*pyz - yy*pxz;
						*/
						// we need Cardisean coordinate whose original point coincide with centroid_x^i
						xx = p.gridfn(gfns::gfn__global_x, irho, isigma) - centroid_x;
						yy = p.gridfn(gfns::gfn__global_y, irho, isigma) - centroid_y;
						zz = p.gridfn(gfns::gfn__global_z, irho, isigma) - centroid_z;
						sxx = yy * pxz - zz * pxy;
						sxy = zz * pxx - xx * pxz;
						sxz = xx * pxy - yy * pxx;
						syx = yy * pyz - zz * pyy;
						syy = zz * pxy - xx * pyz;
						syz = xx * pyy - yy * pxy;
						szx = yy * pzz - zz * pyz;
						szy = zz * pxz - xx * pzz;
						szz = xx * pyz - yy * pxz;

						p.gridfn(gfns::gfn__global_xx, irho, isigma) = nx * pxx + ny * pxy + nz * pxz; // p_x
						p.gridfn(gfns::gfn__global_xy, irho, isigma) = nx * pxy + ny * pyy + nz * pyz; // p_y
						p.gridfn(gfns::gfn__global_xz, irho, isigma) = nx * pxz + ny * pyz + nz * pzz; // p_z
						p.gridfn(gfns::gfn__global_yy, irho, isigma) = nx * sxx + ny * syx + nz * szx; // s_x
						p.gridfn(gfns::gfn__global_yz, irho, isigma) = nx * sxy + ny * syy + nz * szy; // s_y
						p.gridfn(gfns::gfn__global_zz, irho, isigma) = nx * sxz + ny * syz + nz * szz; // s_z
#endif
					}
				}
			}
		}

		Px = surface_integral(ps,
							  gfns::gfn__global_xx, true, true, false, // z,y,x direction, even or odd function
							  patch::integration_method__automatic_choice);
		Py = surface_integral(ps,
							  gfns::gfn__global_xy, true, false, true,
							  patch::integration_method__automatic_choice);
		Pz = surface_integral(ps,
							  gfns::gfn__global_xz, false, true, true,
							  patch::integration_method__automatic_choice);
		Sx = surface_integral(ps,
							  gfns::gfn__global_yy, false, false, true,
							  patch::integration_method__automatic_choice);
		Sy = surface_integral(ps,
							  gfns::gfn__global_yz, false, true, false,
							  patch::integration_method__automatic_choice);
		Sz = surface_integral(ps,
							  gfns::gfn__global_zz, true, false, false,
							  patch::integration_method__automatic_choice);
		const double F1o8pi = 1.0 / 8 / PI;
		Px = Px * F1o8pi;
		Py = Py * F1o8pi;
		Pz = Pz * F1o8pi;
		Sx = Sx * F1o8pi;
		Sy = Sy * F1o8pi;
		Sz = Sz * F1o8pi;
	}

	//******************************************************************************

	//
	// This function computes the surface integral of a gridfn over the
	// horizon.
	//
	fp BH_diagnostics::surface_integral(const patch_system &ps,
										int src_gfn, bool src_gfn_is_even_across_xy_plane,
										bool src_gfn_is_even_across_xz_plane,
										bool src_gfn_is_even_across_yz_plane,
										enum patch::integration_method method)
	{
		return ps.integrate_gridfn(src_gfn, src_gfn_is_even_across_xy_plane,
								   src_gfn_is_even_across_xz_plane,
								   src_gfn_is_even_across_yz_plane,
								   gfns::gfn__h,
								   gfns::gfn__g_dd_11, gfns::gfn__g_dd_12, gfns::gfn__g_dd_13,
								   gfns::gfn__g_dd_22, gfns::gfn__g_dd_23,
								   gfns::gfn__g_dd_33,
								   method);
	}
	// with triad theta and phi
	// since Thornburg uses vertex center, we will meet nan at pole points
	void BH_diagnostics::compute_signature(patch_system &ps, const double dT)
	{
		for (int pn = 0; pn < ps.N_patches(); ++pn)
		{
			patch &p = ps.ith_patch(pn);

			for (int irho = p.min_irho(); irho <= p.max_irho(); ++irho)
				for (int isigma = p.min_isigma(); isigma <= p.max_isigma(); ++isigma)
				{
					const fp r = p.ghosted_gridfn(gfns::gfn__h, irho, isigma);
					const fp rho = p.rho_of_irho(irho);
					const fp sigma = p.sigma_of_isigma(isigma);
					fp xx, yy, zz;
					p.xyz_of_r_rho_sigma(r, rho, sigma, xx, yy, zz);

					const fp sintheta = sqrt(1 - zz * zz / r / r);

					const fp X_ud_11 = xx * zz / r / r / sqrt(xx * xx + yy * yy);
					const fp X_ud_12 = yy * zz / r / r / sqrt(xx * xx + yy * yy);
					const fp X_ud_13 = -sqrt(xx * xx + yy * yy) / r / r;
					const fp X_ud_21 = -yy / (xx * xx + yy * yy);
					const fp X_ud_22 = xx / (xx * xx + yy * yy);
					const fp X_ud_23 = 0;

					const fp g_dd_11 = p.gridfn(gfns::gfn__g_dd_11, irho, isigma);
					const fp g_dd_12 = p.gridfn(gfns::gfn__g_dd_12, irho, isigma);
					const fp g_dd_13 = p.gridfn(gfns::gfn__g_dd_13, irho, isigma);
					const fp g_dd_22 = p.gridfn(gfns::gfn__g_dd_22, irho, isigma);
					const fp g_dd_23 = p.gridfn(gfns::gfn__g_dd_23, irho, isigma);
					const fp g_dd_33 = p.gridfn(gfns::gfn__g_dd_33, irho, isigma);

					const fp Lap = 1.0 + p.gridfn(gfns::gfn__global_xx, irho, isigma);
					const fp Sfx = p.gridfn(gfns::gfn__global_xy, irho, isigma);
					const fp Sfy = p.gridfn(gfns::gfn__global_xz, irho, isigma);
					const fp Sfz = p.gridfn(gfns::gfn__global_yy, irho, isigma);

					const fp dfdt = (r - p.gridfn(gfns::gfn__oldh, irho, isigma)) / dT;

					double Br = Sfx * xx / r + Sfy * yy / r + Sfz * zz / r;
					double Brho = Sfx * X_ud_11 + Sfy * X_ud_12 + Sfz * X_ud_13;
					double Bsigma = Sfx * X_ud_21 + Sfy * X_ud_22 + Sfz * X_ud_23;

					double g_uu_11, g_uu_12, g_uu_13, g_uu_22, g_uu_23, g_uu_33;
					double g11, g12, g13, g22, g23, g33;
					{
						// g^uu
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
						g11 = (t1 * t2 - t5) * t21;
						g12 = -(t11 * t2 - t14 * t4) * t21;
						g13 = -(-t11 * t4 + t14 * t1) * t21;
						g22 = (t7 * t2 - t18) * t21;
						g23 = -(t7 * t4 - t15) * t21;
						g33 = (t8 - t12) * t21;
					}
					// 1 r;2 rho; 3 sigma
					g_uu_22 = (g11 * X_ud_11 + g12 * X_ud_12 + g13 * X_ud_13) * X_ud_11 + (g12 * X_ud_11 + g22 * X_ud_12 + g23 * X_ud_13) * X_ud_12 + (g13 * X_ud_11 + g23 * X_ud_12 + g33 * X_ud_13) * X_ud_13;
					g_uu_23 = (g11 * X_ud_11 + g12 * X_ud_12 + g13 * X_ud_13) * X_ud_21 + (g12 * X_ud_11 + g22 * X_ud_12 + g23 * X_ud_13) * X_ud_22 + (g13 * X_ud_11 + g23 * X_ud_12 + g33 * X_ud_13) * X_ud_23;
					g_uu_12 = (g11 * X_ud_11 + g12 * X_ud_12 + g13 * X_ud_13) * xx / r + (g12 * X_ud_11 + g22 * X_ud_12 + g23 * X_ud_13) * yy / r + (g13 * X_ud_11 + g23 * X_ud_12 + g33 * X_ud_13) * zz / r;
					g_uu_33 = (g11 * X_ud_21 + g12 * X_ud_22 + g13 * X_ud_23) * X_ud_21 + (g12 * X_ud_21 + g22 * X_ud_22 + g23 * X_ud_23) * X_ud_22 + (g13 * X_ud_21 + g23 * X_ud_22 + g33 * X_ud_23) * X_ud_23;
					g_uu_13 = (g11 * X_ud_21 + g12 * X_ud_22 + g13 * X_ud_23) * xx / r + (g12 * X_ud_21 + g22 * X_ud_22 + g23 * X_ud_23) * yy / r + (g13 * X_ud_21 + g23 * X_ud_22 + g33 * X_ud_23) * zz / r;
					g_uu_11 = (g11 * xx / r + g12 * yy / r + g13 * zz / r) * xx / r + (g12 * xx / r + g22 * yy / r + g23 * zz / r) * yy / r + (g13 * xx / r + g23 * yy / r + g33 * zz / r) * zz / r;
					{
						// g_uu
						fp t1, t2, t4, t5, t7, t8, t11, t12, t14, t15;
						fp t18, t21;
						t1 = g_uu_22;
						t2 = g_uu_33;
						t4 = g_uu_23;
						t5 = t4 * t4;
						t7 = g_uu_11;
						t8 = t7 * t1;
						t11 = g_uu_12;
						t12 = t11 * t11;
						t14 = g_uu_13;
						t15 = t11 * t14;
						t18 = t14 * t14;
						t21 = 1 / (t8 * t2 - t7 * t5 - t12 * t2 + 2.0 * t15 * t4 - t18 * t1);
						g11 = (t1 * t2 - t5) * t21;
						g12 = -(t11 * t2 - t14 * t4) * t21;
						g13 = -(-t11 * t4 + t14 * t1) * t21;
						g22 = (t7 * t2 - t18) * t21;
						g23 = -(t7 * t4 - t15) * t21;
						g33 = (t8 - t12) * t21;
					}

					double q11 = g22, q12 = g23, q13 = Br + dfdt * g12;
					double q22 = g33, q23 = Bsigma + dfdt * g13;
					double q33 = (-Lap * Lap + g11 * Br * Br + g22 * Brho * Brho + g33 * Bsigma * Bsigma +
								  2 * (g12 * Br * Brho + g13 * Br * Bsigma + g23 * Brho * Bsigma)) +
								 2 * dfdt * Br + dfdt * dfdt * g11;
					q12 = q12 / sintheta;
					q22 = q22 / sintheta / sintheta;
					q23 = q23 / sintheta;
					// we use gfns::gfn__global_zz to store determinant
					p.gridfn(gfns::gfn__global_zz, irho, isigma) = q11 * q22 * q33 + q12 * q23 * q13 + q13 * q12 * q23 - q13 * q22 * q13 - q12 * q12 * q33 - q11 * q23 * q23;
				} // end for irho isigma
		}
	}
	FILE *BH_diagnostics::setup_output_file(int N_horizons, int hn)
		const
	{
		char file_name_buffer[50];
		sprintf(file_name_buffer, "infoah%02d.dat", hn);
		const char *const file_open_mode = "w";

		FILE *fileptr = fopen(file_name_buffer, file_open_mode);
		if (fileptr == NULL)
			printf("\n"
				   "   BH_diagnostics::setup_output_file():\n"
				   "        can't open BH-diagnostics output file\n"
				   "        \"%s\"!",
				   file_name_buffer);
		/*
		fprintf(fileptr, "# apparent horizon %d/%d\n", hn, N_horizons);
		fprintf(fileptr, "#\n");
		fprintf(fileptr, "# column  1 = cctk_time\n");
		fprintf(fileptr, "# column  2 = centroid_x\n");
		fprintf(fileptr, "# column  3 = centroid_y\n");
		fprintf(fileptr, "# column  4 = centroid_z\n");
		fprintf(fileptr, "# column  5 = min radius\n");
		fprintf(fileptr, "# column  6 = max radius\n");
		fprintf(fileptr, "# column  7 = mean radius\n");
		fprintf(fileptr, "# column  8 = quadrupole_xx\n");
		fprintf(fileptr, "# column  9 = quadrupole_xy\n");
		fprintf(fileptr, "# column 10 = quadrupole_xz\n");
		fprintf(fileptr, "# column 11 = quadrupole_yy\n");
		fprintf(fileptr, "# column 12 = quadrupole_yz\n");
		fprintf(fileptr, "# column 13 = quadrupole_zz\n");
		fprintf(fileptr, "# column 14 = min x\n");
		fprintf(fileptr, "# column 15 = max x\n");
		fprintf(fileptr, "# column 16 = min y\n");
		fprintf(fileptr, "# column 17 = max y\n");
		fprintf(fileptr, "# column 18 = min z\n");
		fprintf(fileptr, "# column 19 = max z\n");
		fprintf(fileptr, "# column 20 = xy-plane circumference\n");
		fprintf(fileptr, "# column 21 = xz-plane circumference\n");
		fprintf(fileptr, "# column 22 = yz-plane circumference\n");
		fprintf(fileptr, "# column 23 = ratio of xz/xy-plane circumferences\n");
		fprintf(fileptr, "# column 24 = ratio of yz/xy-plane circumferences\n");
		fprintf(fileptr, "# column 25 = area\n");
		fprintf(fileptr, "# column 26 = irreducible mass\n");
		fprintf(fileptr, "# column 27 = areal radius\n");
		*/

		fprintf(fileptr, "#time Mass x y z Px Py Pz Sx Sy Sz\n");
		fflush(fileptr);

		return fileptr;
	}
	void BH_diagnostics::output(FILE *fileptr, double time)
		const
	{
		assert(fileptr != NULL);
		/*
		fprintf(fileptr,
			"%f\t%f\t%f\t%f\t%#.10g\t%#.10g\t%#.10g\t",
			double(time),
			double(centroid_x), double(centroid_y), double(centroid_z),
			double(min_radius), double(max_radius), double(mean_radius));

		fprintf(fileptr,
			"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
			double(quadrupole_xx), double(quadrupole_xy), double(quadrupole_xz),
						   double(quadrupole_yy), double(quadrupole_yz),
									  double(quadrupole_zz));

		fprintf(fileptr,
			"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
			double(min_x), double(max_x),
			double(min_y), double(max_y),
			double(min_z), double(max_z));

		fprintf(fileptr,
			"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t",
			double(circumference_xy),
			double(circumference_xz),
			double(circumference_yz),
			double(circumference_xz / circumference_xy),
			double(circumference_yz / circumference_xy));

		fprintf(fileptr,
			"%#.10g\t%#.10g\t%#.10g\n",
			double(area), double(irreducible_mass), double(areal_radius));
		*/

		fprintf(fileptr,
				"%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\t%#.10g\n",
				double(time), double(irreducible_mass),
				double(centroid_x), double(centroid_y), double(centroid_z),
				double(Px), double(Py), double(Pz), double(Sx), double(Sy), double(Sz));

		fflush(fileptr);
	}

} // namespace AHFinderDirect
