#ifndef BH_DIAGNOSTICS_H
#define BH_DIAGNOSTICS_H
namespace AHFinderDirect
{

	struct BH_diagnostics
	{
	public:
		// mean x,y,z
		fp centroid_x, centroid_y, centroid_z;

		// these are quadrupole moments about the centroid, i.e.
		// mean(xi*xj) - centroid_i*centroid_j
		fp quadrupole_xx, quadrupole_xy, quadrupole_xz,
			quadrupole_yy, quadrupole_yz,
			quadrupole_zz;

		// min,max,mean surface radius about local coordinate origin
		fp min_radius, max_radius, mean_radius;

		// xyz bounding box
		fp min_x, max_x,
			min_y, max_y,
			min_z, max_z;

		// proper circumference
		// (computed using induced metric along these local-coordinate planes)
		fp circumference_xy,
			circumference_xz,
			circumference_yz;

		// surface area (computed using induced metric)
		// and quantities derived from it
		fp area, irreducible_mass, areal_radius;

		double Px, Py, Pz, Sx, Sy, Sz;

	public:
		// position of diagnostics in buffer and number of diagnostics
		enum
		{
			posn__centroid_x = 0,
			posn__centroid_y,
			posn__centroid_z,
			posn__quadrupole_xx,
			posn__quadrupole_xy,
			posn__quadrupole_xz,
			posn__quadrupole_yy,
			posn__quadrupole_yz,
			posn__quadrupole_zz,
			posn__min_radius,
			posn__max_radius,
			posn__mean_radius,

			posn__min_x,
			posn__max_x,
			posn__min_y,
			posn__max_y,
			posn__min_z,
			posn__max_z,

			posn__circumference_xy,
			posn__circumference_xz,
			posn__circumference_yz,

			posn__area,
			posn__irreducible_mass,
			posn__areal_radius,

			N_buffer // no comma	// size of buffer
		};

		// copy diagnostics to/from buffer
		void copy_to_buffer(double buffer[N_buffer]) const;
		void copy_from_buffer(const double buffer[N_buffer]);

	public:
		void compute(patch_system &ps);

		void compute_signature(patch_system &ps, const double dT);

		FILE *setup_output_file(int N_horizons, int hn)
			const;

		void output(FILE *fileptr, double time)
			const;

		BH_diagnostics();

	private:
		static double surface_integral(const patch_system &ps,
									   int src_gfn, bool src_gfn_is_even_across_xy_plane,
									   bool src_gfn_is_even_across_xz_plane,
									   bool src_gfn_is_even_across_yz_plane,
									   enum patch::integration_method method);
	};

	//******************************************************************************

} // namespace AHFinderDirect
#endif /* BH_DIAGNOSTICS_H */
