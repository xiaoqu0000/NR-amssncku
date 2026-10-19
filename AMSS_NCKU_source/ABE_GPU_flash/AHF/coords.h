#ifndef COORDS_H
#define COORDS_H
namespace AHFinderDirect
{
	namespace local_coords
	{

		// compare if two angles are fuzzily equal mod 2*pi radians (360 degrees)
		bool fuzzy_EQ_ang(fp ang1, fp ang2);	// radians
		bool fuzzy_EQ_dang(fp dang1, fp dang2); // degrees

		// modulo-reduce  {ang,dang}  to be (fuzzily) within the range
		// [min,max]_{ang,dang}, or error_exit() if no such value exists
		fp modulo_reduce_ang(fp ang, fp min_ang, fp max_ang);
		fp modulo_reduce_dang(fp dang, fp min_dang, fp max_dang);

	} // close namespace local_coords::

	namespace local_coords
	{
		// (r,(mu,nu,phi)) <--> (x,y,z)
		void xyz_of_r_mu_nu(fp r, fp mu, fp nu, fp &x, fp &y, fp &z);
		void xyz_of_r_mu_phi(fp r, fp mu, fp phi, fp &x, fp &y, fp &z);
		void xyz_of_r_nu_phi(fp r, fp nu, fp phi, fp &x, fp &y, fp &z);
		fp r_of_xyz(fp x, fp y, fp z);
		fp mu_of_yz(fp y, fp z);
		fp nu_of_xz(fp x, fp z);
		fp phi_of_xy(fp x, fp y);

		// ((mu,nu,phi)) --> the 3rd
		fp phi_of_mu_nu(fp mu, fp nu);
		fp nu_of_mu_phi(fp mu, fp phi);
		fp mu_of_nu_phi(fp nu, fp phi);

		// partial {x,y,z} / partial {mu,nu,phi}
		void partial_xyz_wrt_r_mu_nu(fp r, fp mu, fp nu,
									 fp &partial_x_wrt_r, fp &partial_x_wrt_mu, fp &partial_x_wrt_nu,
									 fp &partial_y_wrt_r, fp &partial_y_wrt_mu, fp &partial_y_wrt_nu,
									 fp &partial_z_wrt_r, fp &partial_z_wrt_mu, fp &partial_z_wrt_nu);
		void partial_xyz_wrt_r_mu_phi(fp r, fp mu, fp phi,
									  fp &partial_x_wrt_r, fp &partial_x_wrt_mu, fp &partial_x_wrt_phi,
									  fp &partial_y_wrt_r, fp &partial_y_wrt_mu, fp &partial_y_wrt_phi,
									  fp &partial_z_wrt_r, fp &partial_z_wrt_mu, fp &partial_z_wrt_phi);
		void partial_xyz_wrt_r_nu_phi(fp r, fp nu, fp phi,
									  fp &partial_x_wrt_r, fp &partial_x_wrt_nu, fp &partial_x_wrt_phi,
									  fp &partial_y_wrt_r, fp &partial_y_wrt_nu, fp &partial_y_wrt_phi,
									  fp &partial_z_wrt_r, fp &partial_z_wrt_nu, fp &partial_z_wrt_phi);

		// partial {mu,nu,phi} / partial {x,y,z}
		fp partial_mu_wrt_y(fp y, fp z);
		fp partial_mu_wrt_z(fp y, fp z);
		fp partial_nu_wrt_x(fp x, fp z);
		fp partial_nu_wrt_z(fp x, fp z);
		fp partial_phi_wrt_x(fp x, fp y);
		fp partial_phi_wrt_y(fp x, fp y);

		// partial^2 {mu,nu,phi} / partial {x,y,z}{x,y,z}
		fp partial2_mu_wrt_yy(fp y, fp z);
		fp partial2_mu_wrt_yz(fp y, fp z);
		fp partial2_mu_wrt_zz(fp y, fp z);
		fp partial2_nu_wrt_xx(fp x, fp z);
		fp partial2_nu_wrt_xz(fp x, fp z);
		fp partial2_nu_wrt_zz(fp x, fp z);
		fp partial2_phi_wrt_xx(fp x, fp y);
		fp partial2_phi_wrt_xy(fp x, fp y);
		fp partial2_phi_wrt_yy(fp x, fp y);

		// usual polar spherical (r,theta,phi) <--> (x,y,z)
		void xyz_of_r_theta_phi(fp r, fp theta, fp phi, fp &x, fp &y, fp &z);
		void r_theta_phi_of_xyz(fp x, fp y, fp z, fp &r, fp &theta, fp &phi);
		// ... already have r_of_xyz()
		// ... already have phi_of_xy()
		fp theta_of_xyz(fp x, fp y, fp z);

		// ((mu,nu,phi)) <--> usual polar spherical (theta,phi)
		// ... note phi is the same coordinate in both systems
		void theta_phi_of_mu_nu(fp mu, fp nu, fp &ps_theta, fp &ps_phi);
		void theta_phi_of_mu_phi(fp mu, fp phi, fp &ps_theta, fp &ps_phi);
		void theta_phi_of_nu_phi(fp nu, fp phi, fp &ps_theta, fp &ps_phi);
		void mu_nu_of_theta_phi(fp ps_theta, fp ps_phi, fp &mu, fp &nu);
		void mu_phi_of_theta_phi(fp ps_theta, fp ps_phi, fp &mu, fp &phi);
		void nu_phi_of_theta_phi(fp ps_theta, fp ps_phi, fp &nu, fp &phi);

		// ((mu,nu,phi)) --> direction cosines (xcos,ycos,zcos)
		void xyzcos_of_mu_nu(fp mu, fp nu, fp &xcos, fp &ycos, fp &zcos);
		void xyzcos_of_mu_phi(fp mu, fp phi, fp &xcos, fp &ycos, fp &zcos);
		void xyzcos_of_nu_phi(fp nu, fp phi, fp &xcos, fp &ycos, fp &zcos);
	} // close namespace local_coords::

	//*****************************************************************************

	//
	// ***** bit masks for coordinates ****
	//

	//
	// We need to manipulate coordinates to do calculations like "which
	// coordinate do these two patches have in common".  We do these by
	// Boolean operations on integers using the following bit masks:
	//

	namespace local_coords
	{

		typedef int coords_set;

		enum
		{
			coords_set_mu = 0x1,
			coords_set_nu = 0x2,
			coords_set_phi = 0x4,

			coords_set_empty = 0x0,
			coords_set_all = coords_set_mu | coords_set_nu | coords_set_phi // no comma
		};

		// human-readable coordinate names for debugging etc
		const char *name_of_coords_set(coords_set S);

		// set complement of coordinates
		inline coords_set coords_set_not(coords_set S)
		{
			return coords_set_all & ~S;
		}

	} // close namespace local_coords::

	//******************************************************************************

	//
	// This class stores the origin point of our local coordinates, and
	// provides conversions between local and global coordinates.
	//
	class global_coords
	{
	public:
		// get global (x,y,z) coordinates of local origin point
		fp origin_x() const { return origin_x_; }
		fp origin_y() const { return origin_y_; }
		fp origin_z() const { return origin_z_; }

		// constructor: specify global (x,y,z) coordinates of local origin point
		global_coords(fp origin_x_in, fp origin_y_in, fp origin_z_in)
			: origin_x_(origin_x_in),
			  origin_y_(origin_y_in),
			  origin_z_(origin_z_in)
		{
		}
		// destructor: compiler-generated no-op is ok

		void recentering(fp x, fp y, fp z)
		{
			origin_x_ = x;
			origin_y_ = y;
			origin_z_ = z;
		}

	private:
		// we forbid copying and passing by value
		// by declaring the copy constructor and assignment operator
		// private, but never defining them
		global_coords(const global_coords &rhs);
		global_coords &operator=(const global_coords &rhs);

	private:
		// global (x,y,z) coordinates of local origin point
		fp origin_x_, origin_y_, origin_z_;
	};

	//******************************************************************************

} // namespace AHFinderDirect
#endif /*  COORDS_H  */
