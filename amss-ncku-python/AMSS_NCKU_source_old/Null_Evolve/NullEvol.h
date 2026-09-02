
#ifndef NULLEVOL_H
#define NULLEVOL_H

#ifdef fortran1
#define f_setup_dyad setup_dyad
#define f_eth_derivs eth_derivs
#define f_eth_dderivs eth_dderivs
#define f_fill_symmetric_boundarybuffer fill_symmetric_boundarybuffer
#define f_fill_symmetric_boundarybuffer2 fill_symmetric_boundarybuffer2
#define f_calculate_K calculate_k
#define f_NullEvol_beta nullevol_beta
#define f_NullEvol_Q nullevol_q
#define f_NullEvol_U nullevol_u
#define f_NullEvol_W nullevol_w
#define f_NullEvol_Theta nullevol_theta
#define f_NullEvol_Theta_givenx nullevol_theta_givenx
#define f_Eq_Theta eq_theta
#define f_Eq_Theta_2 eq_theta_2
#define f_NullEvol_g01 nullevol_g01
#define f_NullEvol_pg0A nullevol_pg0a
#define f_NullEvol_Theta2 nullevol_theta2
#define f_NullEvol_Thetag00 nullevol_thetag00
#endif
#ifdef fortran2
#define f_setup_dyad SETUP_DYAD
#define f_eth_derivs ETH_DERIVS
#define f_eth_dderivs ETH_DDERIVS
#define f_fill_symmetric_boundarybuffer FILL_SYMMETRIC_BOUNDARYBUFFER
#define f_fill_symmetric_boundarybuffer2 FILL_SYMMETRIC_BOUNDARYBUFFER2
#define f_calculate_K CALCULATE_K
#define f_NullEvol_beta NULLEVOL_BETA
#define f_NullEvol_Q NULLEVOL_Q
#define f_NullEvol_U NULLEVOL_U
#define f_NullEvol_W NULLEVOL_W
#define f_NullEvol_Theta NULLEVOL_THETA
#define f_NullEvol_Theta_givenx NULLEVOL_THETA_GIVENX
#define f_Eq_Theta EQ_THETA
#define f_Eq_Theta_2 EQ_THETA_2
#define f_NullEvol_g01 NULLEVOL_G01
#define f_NullEvol_pg0A NULLEVOL_PG0A
#define f_NullEvol_Theta2 NULLEVOL_THETA2
#define f_NullEvol_Thetag00 NULLEVOL_THETAG00
#endif
#ifdef fortran3
#define f_setup_dyad setup_dyad_
#define f_eth_derivs eth_derivs_
#define f_eth_dderivs eth_dderivs_
#define f_fill_symmetric_boundarybuffer fill_symmetric_boundarybuffer_
#define f_fill_symmetric_boundarybuffer2 fill_symmetric_boundarybuffer2_
#define f_calculate_K calculate_k_
#define f_NullEvol_beta nullevol_beta_
#define f_NullEvol_Q nullevol_q_
#define f_NullEvol_U nullevol_u_
#define f_NullEvol_W nullevol_w_
#define f_NullEvol_Theta nullevol_theta_
#define f_NullEvol_Theta_givenx nullevol_theta_givenx_
#define f_Eq_Theta eq_theta_
#define f_Eq_Theta_2 eq_theta_2_
#define f_NullEvol_g01 nullevol_g01_
#define f_NullEvol_pg0A nullevol_pg0a_
#define f_NullEvol_Theta2 nullevol_theta2_
#define f_NullEvol_Thetag00 nullevol_thetag00_
#endif

extern "C"
{
	void f_setup_dyad(int *, double *, double *, double *,
					  double *, double *, double *, double *,
					  double *, double *, double *, double *,
					  double *, double *,
					  double *, double *, double *, double *,
					  double *, double *, double *, double *,
					  double *, double *, double *, double *,
					  double *, double *, double *,
					  int &, double &);
}

extern "C"
{
	void f_eth_derivs(int *, double *, double *,
					  double *, double *,
					  double *, double *,
					  int &, int &,
					  double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	void f_eth_dderivs(int *, double *, double *,
					   double *, double *,
					   double *, double *,
					   int &, int &, int &,
					   double *, double *, double *, double *, double *, double *,
					   double *, double *, double *, double *,
					   double *, double *, double *, double *,
					   double *, double *, double *, double *);
}

extern "C"
{
	void f_fill_symmetric_boundarybuffer(int *, double *, double *, double *,
										 double &, double &,
										 double *, double *, double *, double *, double *, double *, double *, double *,
										 double *, double *, int &, int &, int &);
}

extern "C"
{
	void f_fill_symmetric_boundarybuffer2(int *, double *, double *, double *,
										  double &, double &,
										  double *, int &, int &, double *);
}

extern "C"
{
	void f_calculate_K(int *, double *, double *, double *,
					   double *, double *,
					   double *, double *, double *, double *);
}

extern "C"
{
	int f_NullEvol_beta(int *, double *, double *, double *,
						double *, double *, double *, double *, double *);
}

extern "C"
{
	int f_NullEvol_Q(int *, double *, double *, double *,
					 double *, double *, double *, double *, double *, double *,
					 double *, double *, double *, double *,
					 double *, double *, double *, double *,
					 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	int f_NullEvol_U(int *, double *, double *, double *,
					 double *, double *, double *, double *,
					 double *, double *, double *,
					 double *, double *, double &);
}

extern "C"
{
	int f_NullEvol_W(int *, double *, double *, double *,
					 double *, double *, double *, double *, double *, double *, double *, double *,
					 double *, double *, double *, double *,
					 double *, double *, double *, double *, double &,
					 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	int f_NullEvol_Theta(int *, double *, double *, double *,
						 double *, double *, double *, double *, double *, double *, double *,
						 double *, double *, double *, double *, double *, double *, double *,
						 double *, double *, double *, double *,
						 double &,
						 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	int f_NullEvol_Theta_givenx(int *, double *, double *, double *,
								double *, double *, double *, double *, double *, double *, double *,
								double *, double *, double *, double *, double *, double *, double *,
								double *, double *, double *, double *,
								double &,
								double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
								double *, double *, double *, double *,
								double *, double *, double *, double *,
								double *, double *, double *, double *,
								double &, int &);
}

extern "C"
{
	int f_Eq_Theta(int *, double *, double *, double *,
				   double *, double *, double *, double *, double *, double *, double *,
				   double *, double *, double *, double *, double *, double *, double *, double &,
				   double *, double *, double *, double *, double *, double *, double *, double *, double *, double *);
}

extern "C"
{
	int f_Eq_Theta_2(int *, double *, double *, double *,
					 double *, double *, double *, double *, double *, double *, double *,
					 double *, double *, double *, double *, double *, double *, double *, double &,
					 double *, double *, double *, double *, double *, double *, double *, double *, double *, double *,
					 double &, int &);
}

extern "C"
{
	int f_NullEvol_g01(int *, double *, double *, double *,
					   double *, double *, double *, double *,
					   double &);
}

extern "C"
{
	int f_NullEvol_pg0A(int *, double *, double *, double *,
						double *, double *, double *, double *,
						double *, double *, double *, double *,
						double &);
}

extern "C"
{
	int f_NullEvol_Theta2(int *, double *, double *, double *,
						  double *, double *, double *, double *, double *, double *, double *, double *, double *,
						  double *, double *, double *,
						  double &);
}

extern "C"
{
	int f_NullEvol_Thetag00(int *, double *, double *, double *,
							double *, double *, double *, double *, double *, double *, double *, double *, double *,
							double *, double *, double *,
							double &);
}
#endif /* NULLEVOL_H */
