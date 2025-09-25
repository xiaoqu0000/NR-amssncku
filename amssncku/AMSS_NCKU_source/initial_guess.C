#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>

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

#include "horizon_sequence.h"
#include "BH_diagnostics.h"
#include "myglobal.h"

namespace AHFinderDirect
{
	extern struct state state;
	//******************************************************************************

	// ellipsoid has global-coordinates center (A,B,C), radius (a,b,c)
	// angular coordinate system has center (U,V,W)
	//
	// direction cosines wrt angular coordinate center are (xcos,ycos,zcos)
	// i.e. a point has coordinates (U+xcos*r, V+ycos*r, W+zcos*r)
	//
	// then the equation of the ellipsoid is
	//	(U+xcos*r - A)^2     (V+ycos*r - B)^2     (W+zcos*r - C)^2
	//	-----------------  +  ----------------  +  -----------------  =  1
	//	        a^2                  b^2                   c^2
	//
	// to solve this, we introduce intermediate variables
	//	AU = A - U
	//	BV = B - V
	//	CW = C - W
	//
	void setup_initial_guess(patch_system &ps,
							 fp x_center, fp y_center, fp z_center,
							 fp x_radius, fp y_radius, fp z_radius)
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
					const fp rho = p.rho_of_irho(irho);
					const fp sigma = p.sigma_of_isigma(isigma);
					fp xcos, ycos, zcos;
					p.xyzcos_of_rho_sigma(rho, sigma, xcos, ycos, zcos);

					// set up variables used by Maple-generated code
					const fp AU = x_center - ps.origin_x();
					const fp BV = y_center - ps.origin_y();
					const fp CW = z_center - ps.origin_z();
					const fp a = x_radius;
					const fp b = y_radius;
					const fp c = z_radius;

					// compute the solutions r_plus and r_minus
					fp r_plus, r_minus;
					{
						fp t1, t2, t3, t5, t6, t7, t9, t10, t12, t28;
						fp t30, t33, t35, t36, t40, t42, t43, t48, t49, t52;
						fp t55;
						t1 = a * a;
						t2 = b * b;
						t3 = t1 * t2;
						t5 = t3 * zcos * CW;
						t6 = c * c;
						t7 = t1 * t6;
						t9 = t7 * ycos * BV;
						t10 = t2 * t6;
						t12 = t10 * xcos * AU;
						t28 = xcos * xcos;
						t30 = CW * CW;
						t33 = BV * BV;
						t35 = t10 * t28;
						t36 = ycos * ycos;
						t40 = AU * AU;
						t42 = t7 * t36;
						t43 = zcos * zcos;
						t48 = t3 * t43;
						t49 = -2.0 * t1 * zcos * CW * ycos * BV - 2.0 * t2 * zcos * CW * xcos * AU - 2.0 * t6 * ycos * BV * xcos * AU + t2 * t28 * t30 + t6 * t28 * t33 - t35 + t1 * t36 * t30 + t6 * t36 * t40 - t42 + t1 * t43 * t33 + t2 * t43 * t40 -
							  t48;
						t52 = sqrt(-t3 * t6 * t49);
						t55 = 1 / (t35 + t42 + t48);
						r_plus = (t5 + t9 + t12 + t52) * t55;
						r_minus = (t5 + t9 + t12 - t52) * t55;
					}

					// exactly one of the solutions (call it r) should be positive
					fp r;
					if ((r_plus > 0.0) && (r_minus < 0.0))
						then r = r_plus;
					else if ((r_plus < 0.0) && (r_minus > 0.0))
						then r = r_minus;
					else if (state.my_proc == 0)
						printf("\nsetup_coord_ellipsoid():\nexpected exactly one r>0 solution to quadratic, got 0 or 2!\n%s patch (irho,isigma)=(%d,%d) ==> (rho,sigma)=(%g,%g)\ndirection cosines (xcos,ycos,zcos)=(%g,%g,%g)\nr_plus=%g r_minus=%g\n==> this probably means the initial guess surface doesn't contain\nthe local origin point, or more generally that the initial\nguess surface isn't a Strahlkoerper (\"star-shaped region\")\nwith respect to the local origin point\n", p.name(), irho, isigma, double(rho), double(sigma), double(xcos), double(ycos), double(zcos), double(r_plus), double(r_minus));

					// r = horizon radius at this grid point
					p.ghosted_gridfn(gfns::gfn__h, irho, isigma) = r;
				}
			}
		}
	}

	//******************************************************************************

} // namespace AHFinderDirect
