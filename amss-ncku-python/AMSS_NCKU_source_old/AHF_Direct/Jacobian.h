#ifndef AHFINDERDIRECT__JACOBIAN_HH
#define AHFINDERDIRECT__JACOBIAN_HH

namespace AHFinderDirect
{
	class Jacobian
	{
	public:
		// basic meta-info
		patch_system &my_patch_system() const { return ps_; }
		int N_rows() const { return N_rows_; }

		// convert (patch,irho,isigma) <--> row/column index
		int II_of_patch_irho_isigma(const patch &p, int irho, int isigma)
			const
		{
			return ps_.gpn_of_patch_irho_isigma(p, irho, isigma);
		}
		const patch &patch_irho_isigma_of_II(int II, int &irho, int &isigma)
			const
		{
			return ps_.patch_irho_isigma_of_gpn(II, irho, isigma);
		}

		double element(int II, int JJ) const;

		// is the matrix element (II,JJ) stored explicitly?
		bool is_explicitly_stored(int II, int JJ) const
		{
			return find_element(II, JJ) > 0;
		}

		int IO() const { return IO_; }
		enum
		{
			C_index_origin = 0,
			Fortran_index_origin = 1
		};

		void zero_matrix();

		void set_element(int II, int JJ, fp value);

		void sum_into_element(int II, int JJ, fp value);

		int find_element(int II, int JJ) const;

		int insert_element(int II, int JJ, fp value);

		void grow_arrays();

		enum
		{
			base_growth_amount = 1000
		};

		void sort_each_row_into_column_order();

		double solve_linear_system(int rhs_gfn, int x_gfn,
								   bool print_msg_flag);

	public:
		Jacobian(patch_system &ps);
		~Jacobian();

	protected:
		patch_system &ps_;
		int N_rows_;

		int IO_;

		int N_nonzeros_;
		int current_N_rows_;

		int N_nonzeros_allocated_;

		int *IA_;

		int *JA_;

		double *A_;

		int *itemp_;
		double *rtemp_;
	};

	//******************************************************************************

} // namespace AHFinderDirect
#endif /* AHFINDERDIRECT__JACOBIAN_HH */
