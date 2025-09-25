#ifndef GR_H
#define GR_H
namespace AHFinderDirect
{

	enum expansion_status
	{
		expansion_success,

		expansion_failure__surface_nonfinite,

		expansion_failure__surface_too_large,

		expansion_failure__surface_outside_grid,

		expansion_failure__surface_in_excised_region,

		expansion_failure__geometry_nonfinite,

		expansion_failure__gij_not_positive_definite // no comma
	};

	// expansion.cc
	enum expansion_status
	expansion(patch_system *ps_ptr, fp add_to_expansion,
			  bool initial_flag,
			  bool Jacobian_flag = false,
			  jtutil::norm<fp> *H_norms_ptr = NULL);

	// expansion_Jacobian.cc
	enum expansion_status
	expansion_Jacobian(patch_system *ps_ptr, Jacobian *Jac_ptr,
					   fp add_to_expansion,
					   bool initial_flag,
					   bool print_msg_flag = false);

	//******************************************************************************

} // namespace AHFinderDirect
#endif /* GR_H  */
