#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "cctk.h"

#include "config.h"
#include "stdc.h"
#include "util.h"
#include "array.h"
#include "linear_map.h"

#include "coords.h"
#include "tgrid.h"
#include "fd_grid.h"

namespace AHFinderDirect
{
  using jtutil::error_exit;

  //*****************************************************************************

  //
  // This function computes a single coefficient of a 1st derivative
  // molecule, for unit grid spacing.
  //
  // static
  fp fd_grid::dx_coeff(int m)
  {
    switch (m)
    {
    case -2:
      return FD_GRID__ORDER4__DX__COEFF_M2;
    case -1:
      return FD_GRID__ORDER4__DX__COEFF_M1;
    case 0:
      return FD_GRID__ORDER4__DX__COEFF_0;
    case +1:
      return FD_GRID__ORDER4__DX__COEFF_P1;
    case +2:
      return FD_GRID__ORDER4__DX__COEFF_P2;

    default:
      cout << "***** fd_grid::dx_coeff(): m=" << m << " is outside order=4 molecule radius=" << FD_GRID__MOL_RADIUS << endl;
      abort();
    }
  }

  //*****************************************************************************

  //
  // This function computes a single coefficient of a 2nd derivative
  // molecule, for unit grid spacing.
  //
  // static
  fp fd_grid::dxx_coeff(int m)
  {
    switch (m)
    {
    case -2:
      return FD_GRID__ORDER4__DXX__COEFF_M2;
    case -1:
      return FD_GRID__ORDER4__DXX__COEFF_M1;
    case 0:
      return FD_GRID__ORDER4__DXX__COEFF_0;
    case +1:
      return FD_GRID__ORDER4__DXX__COEFF_P1;
    case +2:
      return FD_GRID__ORDER4__DXX__COEFF_P2;

    default:
      cout << "***** fd_grid::dx_coeff(): m=" << m << " is outside order=4 molecule radius=" << FD_GRID__MOL_RADIUS << endl;
      abort();
    }
  }

  //******************************************************************************

} // namespace AHFinderDirect
