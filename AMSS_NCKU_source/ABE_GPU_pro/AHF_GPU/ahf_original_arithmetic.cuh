#pragma once
// Generated original arithmetic; see extraction.json. No CPU source changes.
using fp=double;
__device__ double induced_metric(double partial_x_wrt_r,double partial_x_wrt_rho,double partial_x_wrt_sigma,double partial_y_wrt_r,double partial_y_wrt_rho,double partial_y_wrt_sigma,double partial_z_wrt_r,double partial_z_wrt_rho,double partial_z_wrt_sigma,double partial_surface_r_wrt_rho,double partial_surface_r_wrt_sigma,double g_xx,double g_xy,double g_xz,double g_yy,double g_yz,double g_zz,double &g_rho_rho,double &g_rho_sigma,double &g_sigma_sigma){
		const fp dx_wrt_rho = partial_x_wrt_rho + partial_x_wrt_r * partial_surface_r_wrt_rho;
		const fp dx_wrt_sigma = partial_x_wrt_sigma + partial_x_wrt_r * partial_surface_r_wrt_sigma;
		const fp dy_wrt_rho = partial_y_wrt_rho + partial_y_wrt_r * partial_surface_r_wrt_rho;
		const fp dy_wrt_sigma = partial_y_wrt_sigma + partial_y_wrt_r * partial_surface_r_wrt_sigma;
		const fp dz_wrt_rho = partial_z_wrt_rho + partial_z_wrt_r * partial_surface_r_wrt_rho;
		const fp dz_wrt_sigma = partial_z_wrt_sigma + partial_z_wrt_r * partial_surface_r_wrt_sigma;

		g_rho_rho = +dx_wrt_rho * dx_wrt_rho * g_xx + 2.0 * dx_wrt_rho * dy_wrt_rho * g_xy + 2.0 * dx_wrt_rho * dz_wrt_rho * g_xz + dy_wrt_rho * dy_wrt_rho * g_yy + 2.0 * dy_wrt_rho * dz_wrt_rho * g_yz + dz_wrt_rho * dz_wrt_rho * g_zz;
		g_rho_sigma = +dx_wrt_rho * dx_wrt_sigma * g_xx + (dx_wrt_rho * dy_wrt_sigma + dy_wrt_rho * dx_wrt_sigma) * g_xy + (dx_wrt_rho * dz_wrt_sigma + dz_wrt_rho * dx_wrt_sigma) * g_xz + dy_wrt_rho * dy_wrt_sigma * g_yy + (dy_wrt_rho * dz_wrt_sigma + dz_wrt_rho * dy_wrt_sigma) * g_yz + dz_wrt_rho * dz_wrt_sigma * g_zz;
		g_sigma_sigma = +dx_wrt_sigma * dx_wrt_sigma * g_xx + 2.0 * dx_wrt_sigma * dy_wrt_sigma * g_xy + 2.0 * dx_wrt_sigma * dz_wrt_sigma * g_xz + dy_wrt_sigma * dy_wrt_sigma * g_yy + 2.0 * dy_wrt_sigma * dz_wrt_sigma * g_yz + dz_wrt_sigma * dz_wrt_sigma * g_zz;

		return g_rho_rho * g_sigma_sigma - (g_rho_sigma*g_rho_sigma);
}
__device__ void original_ps(double r,double xx,double yy,double zz,const double *global,const double *centroid,const double *a,const double *physical,const double ud[2][3],double *out){
const double centroid_x=centroid[0],centroid_y=centroid[1],centroid_z=centroid[2];
const double g_xx=physical[0], k_xx=physical[28];
const double g_xy=physical[4], k_xy=physical[29];
const double g_xz=physical[8], k_xz=physical[30];
const double g_yy=physical[12], k_yy=physical[31];
const double g_yz=physical[16], k_yz=physical[32];
const double g_zz=physical[20], k_zz=physical[33];
const double trk=physical[34];
const double X_ud_11=ud[0][0];
const double X_ud_12=ud[0][1];
const double X_ud_13=ud[0][2];
const double X_ud_21=ud[1][0];
const double X_ud_22=ud[1][1];
const double X_ud_23=ud[1][2];
						double nx, ny, nz;
						nx = xx / r - X_ud_11 * a[0] - X_ud_21 * a[1];
						ny = yy / r - X_ud_12 * a[0] - X_ud_22 * a[1];
						nz = zz / r - X_ud_13 * a[0] - X_ud_23 * a[1];
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
						// we need Cardisean coordinate whose original point coincide with centroid_x^i
						xx = global[0] - centroid_x;
						yy = global[1] - centroid_y;
						zz = global[2] - centroid_z;
						sxx = yy * pxz - zz * pxy;
						sxy = zz * pxx - xx * pxz;
						sxz = xx * pxy - yy * pxx;
						syx = yy * pyz - zz * pyy;
						syy = zz * pxy - xx * pyz;
						syz = xx * pyy - yy * pxy;
						szx = yy * pzz - zz * pyz;
						szy = zz * pxz - xx * pzz;
						szz = xx * pyz - yy * pxz;

						out[0] = nx * pxx + ny * pxy + nz * pxz; // p_x
						out[1] = nx * pxy + ny * pyy + nz * pyz; // p_y
						out[2] = nx * pxz + ny * pyz + nz * pzz; // p_z
						out[3] = nx * sxx + ny * syx + nz * szx; // s_x
						out[4] = nx * sxy + ny * syy + nz * szy; // s_y
						out[5] = nx * sxz + ny * syz + nz * szz; // s_z
}
__device__ double original_signature(double r,double xx,double yy,double zz,const double *physical,const double *gauge,double oldh,double dT){
					const fp sintheta = sqrt(1 - zz * zz / r / r);

					const fp X_ud_11 = xx * zz / r / r / sqrt(xx * xx + yy * yy);
					const fp X_ud_12 = yy * zz / r / r / sqrt(xx * xx + yy * yy);
					const fp X_ud_13 = -sqrt(xx * xx + yy * yy) / r / r;
					const fp X_ud_21 = -yy / (xx * xx + yy * yy);
					const fp X_ud_22 = xx / (xx * xx + yy * yy);
					const fp X_ud_23 = 0;

					const fp g_dd_11 = physical[0];
					const fp g_dd_12 = physical[4];
					const fp g_dd_13 = physical[8];
					const fp g_dd_22 = physical[12];
					const fp g_dd_23 = physical[16];
					const fp g_dd_33 = physical[20];

					const fp Lap = 1.0 + gauge[0];
					const fp Sfx = gauge[1];
					const fp Sfy = gauge[2];
					const fp Sfz = gauge[3];

					const fp dfdt = (r - oldh) / dT;

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
					return q11 * q22 * q33 + q12 * q23 * q13 + q13 * q12 * q23 - q13 * q22 * q13 - q12 * q12 * q33 - q11 * q23 * q23;
}
__device__ void tangent_0(const Point &m,double r,double x,double y,double z,double *v){
const double rinv=1.0/r;
double partial_x_wrt_r=x*rinv,partial_x_wrt_rho,partial_x_wrt_sigma;
double partial_y_wrt_r=y*rinv,partial_y_wrt_rho,partial_y_wrt_sigma;
double partial_z_wrt_r=z*rinv,partial_z_wrt_rho,partial_z_wrt_sigma;
			const fp t = 1 + (m.first*m.first) + (m.second*m.second); // = $r^2/z^2$
			const fp partial_t_wrt_mu = 2.0 * m.first * (1.0 + (m.first*m.first));
			const fp partial_t_wrt_nu = 2.0 * m.second * (1.0 + (m.second*m.second));

			const fp r2_over_zt2 = (r * r) / (z * t * t);
			partial_z_wrt_rho = -0.5 * r2_over_zt2 * partial_t_wrt_mu;
			partial_z_wrt_sigma = -0.5 * r2_over_zt2 * partial_t_wrt_nu;

			partial_x_wrt_rho = m.second * partial_z_wrt_rho;
			partial_x_wrt_sigma = m.second * partial_z_wrt_sigma + z * (1.0 + (m.second*m.second));
			partial_y_wrt_rho = m.first * partial_z_wrt_rho + z * (1.0 + (m.first*m.first));
			partial_y_wrt_sigma = m.first * partial_z_wrt_sigma;
v[0]=partial_x_wrt_r;v[1]=partial_x_wrt_rho;v[2]=partial_x_wrt_sigma;
v[3]=partial_y_wrt_r;v[4]=partial_y_wrt_rho;v[5]=partial_y_wrt_sigma;
v[6]=partial_z_wrt_r;v[7]=partial_z_wrt_rho;v[8]=partial_z_wrt_sigma;
}
__device__ void tangent_1(const Point &m,double r,double x,double y,double z,double *v){
const double rinv=1.0/r;
double partial_x_wrt_r=x*rinv,partial_x_wrt_rho,partial_x_wrt_sigma;
double partial_y_wrt_r=y*rinv,partial_y_wrt_rho,partial_y_wrt_sigma;
double partial_z_wrt_r=z*rinv,partial_z_wrt_rho,partial_z_wrt_sigma;
			const fp t = 1 + (m.first*m.first) + (m.second*m.second); // = $r^2/x^2$
			const fp partial_t_wrt_nu_bar = 2.0 * m.first * (1.0 + (m.first*m.first));
			const fp partial_t_wrt_phi = 2.0 * m.second * (1.0 + (m.second*m.second));

			const fp r2_over_xt2 = (r * r) / (x * t * t);
			partial_x_wrt_rho = 0.5 * r2_over_xt2 * partial_t_wrt_nu_bar;
			partial_x_wrt_sigma = -0.5 * r2_over_xt2 * partial_t_wrt_phi;

			partial_y_wrt_rho = m.second * partial_x_wrt_rho;
			partial_y_wrt_sigma = m.second * partial_x_wrt_sigma + x * (1.0 + (m.second*m.second));
			partial_z_wrt_rho = m.first * partial_x_wrt_rho - x * (1.0 + (m.first*m.first));
			partial_z_wrt_sigma = m.first * partial_x_wrt_sigma;
v[0]=partial_x_wrt_r;v[1]=partial_x_wrt_rho;v[2]=partial_x_wrt_sigma;
v[3]=partial_y_wrt_r;v[4]=partial_y_wrt_rho;v[5]=partial_y_wrt_sigma;
v[6]=partial_z_wrt_r;v[7]=partial_z_wrt_rho;v[8]=partial_z_wrt_sigma;
}
__device__ void tangent_2(const Point &m,double r,double x,double y,double z,double *v){
const double rinv=1.0/r;
double partial_x_wrt_r=x*rinv,partial_x_wrt_rho,partial_x_wrt_sigma;
double partial_y_wrt_r=y*rinv,partial_y_wrt_rho,partial_y_wrt_sigma;
double partial_z_wrt_r=z*rinv,partial_z_wrt_rho,partial_z_wrt_sigma;
			const fp t = 1 + (m.first*m.first) + (m.second*m.second); // = $r^2/y^2$
			const fp partial_t_wrt_mu_bar = 2.0 * m.first * (1.0 + (m.first*m.first));
			const fp partial_t_wrt_phi_bar = 2.0 * m.second * (1.0 + (m.second*m.second));

			const fp r2_over_yt2 = (r * r) / (y * t * t);
			partial_y_wrt_rho = 0.5 * r2_over_yt2 * partial_t_wrt_mu_bar;
			partial_y_wrt_sigma = 0.5 * r2_over_yt2 * partial_t_wrt_phi_bar;

			partial_x_wrt_rho = m.second * partial_y_wrt_rho;
			partial_x_wrt_sigma = m.second * partial_y_wrt_sigma - y * (1.0 + (m.second*m.second));
			partial_z_wrt_rho = m.first * partial_y_wrt_rho - y * (1.0 + (m.first*m.first));
			partial_z_wrt_sigma = m.first * partial_y_wrt_sigma;
v[0]=partial_x_wrt_r;v[1]=partial_x_wrt_rho;v[2]=partial_x_wrt_sigma;
v[3]=partial_y_wrt_r;v[4]=partial_y_wrt_rho;v[5]=partial_y_wrt_sigma;
v[6]=partial_z_wrt_r;v[7]=partial_z_wrt_rho;v[8]=partial_z_wrt_sigma;
}
