#ifndef AHF_NATIVE_GEOMETRY_CUH
#define AHF_NATIVE_GEOMETRY_CUH
#include "ahf_original_scalar.cuh"
#ifdef __CUDACC__
#define AHF_HD __host__ __device__ __forceinline__
#else
#define AHF_HD inline
#endif
namespace ahf_native {
AHF_HD double pow2(double x){return x*x;}
// Exact original coords.C derivative expressions, evaluated at actual local XYZ.
		AHF_HD double partial_mu_wrt_y(double y, double z) { return z / (y * y + z * z); }
		AHF_HD double partial_mu_wrt_z(double y, double z) { return -y / (y * y + z * z); }

		AHF_HD double partial_nu_wrt_x(double x, double z) { return z / (x * x + z * z); }
		AHF_HD double partial_nu_wrt_z(double x, double z) { return -x / (x * x + z * z); }

		AHF_HD double partial_phi_wrt_x(double x, double y) { return -y / (x * x + y * y); }
		AHF_HD double partial_phi_wrt_y(double x, double y) { return x / (x * x + y * y); }
		AHF_HD double partial2_mu_wrt_yy(double y, double z) { return -2.0 * y * z / pow2(y * y + z * z); }
		AHF_HD double partial2_mu_wrt_yz(double y, double z) { return (y * y - z * z) / pow2(y * y + z * z); }
		AHF_HD double partial2_mu_wrt_zz(double y, double z) { return 2.0 * y * z / pow2(y * y + z * z); }

		AHF_HD double partial2_nu_wrt_xx(double x, double z) { return -2.0 * x * z / pow2(x * x + z * z); }
		AHF_HD double partial2_nu_wrt_xz(double x, double z) { return (x * x - z * z) / pow2(x * x + z * z); }
		AHF_HD double partial2_nu_wrt_zz(double x, double z) { return 2.0 * x * z / pow2(x * x + z * z); }

		AHF_HD double partial2_phi_wrt_xx(double x, double y) { return 2.0 * x * y / pow2(x * x + y * y); }
		AHF_HD double partial2_phi_wrt_xy(double x, double y) { return (y * y - x * x) / pow2(x * x + y * y); }
		AHF_HD double partial2_phi_wrt_yy(double x, double y) { return -2.0 * x * y / pow2(x * x + y * y); }

AHF_HD void angular_geometry(int angle,double x,double y,double z,double ud[3],double udd[6]){
    for(int i=0;i<3;++i)ud[i]=0.0;
    for(int i=0;i<6;++i)udd[i]=0.0;
    if(angle==0){
        ud[1]=partial_mu_wrt_y(y,z); ud[2]=partial_mu_wrt_z(y,z);
        udd[3]=partial2_mu_wrt_yy(y,z); udd[4]=partial2_mu_wrt_yz(y,z); udd[5]=partial2_mu_wrt_zz(y,z);
    }else if(angle==1){
        ud[0]=partial_nu_wrt_x(x,z); ud[2]=partial_nu_wrt_z(x,z);
        udd[0]=partial2_nu_wrt_xx(x,z); udd[2]=partial2_nu_wrt_xz(x,z); udd[5]=partial2_nu_wrt_zz(x,z);
    }else{
        ud[0]=partial_phi_wrt_x(x,y); ud[1]=partial_phi_wrt_y(x,y);
        udd[0]=partial2_phi_wrt_xx(x,y); udd[1]=partial2_phi_wrt_xy(x,y); udd[3]=partial2_phi_wrt_yy(x,y);
    }
}
// Family z=0, x=1, y=2. Local coordinates are never translated for derivatives.
AHF_HD bool native_geometry(const StaticMetadata &m,double radius,RuntimeGeometry &g,double global[3]){
    if(!(radius>0.) || !finite_value(radius) || m.family<0 || m.family>2 ||
       !(m.sign==1. || m.sign==-1.) || !(m.normalizer>0.) || !finite_value(m.normalizer) ||
       !finite_value(m.ratio_first) || !finite_value(m.ratio_second))return false;
    for(double v:m.origin)if(!finite_value(v))return false;
    g.r=radius;
    double &x=g.xyz[0],&y=g.xyz[1],&z=g.xyz[2];
    if(m.family==0){
        z=m.sign*radius*m.normalizer;
        x=m.ratio_second*z;
        y=m.ratio_first*z;
    }else if(m.family==1){
        x=m.sign*radius*m.normalizer;
        z=m.ratio_first*x;
        y=m.ratio_second*x;
    }else{
        y=m.sign*radius*m.normalizer;
        z=m.ratio_first*y;
        x=m.ratio_second*y;
    }
    const int rho=m.family==1?1:0;
    const int sigma=m.family==0?1:2;
    // Reject singular/overflowing denominator paths without changing legal expressions.
    const double d0=rho==0?y*y+z*z:x*x+z*z;
    const double d1=sigma==1?x*x+z*z:x*x+y*y;
    if(!(d0>0.) || !(d1>0.) || !finite_value(d0) || !finite_value(d1) ||
       !(d0*d0>0.) || !(d1*d1>0.) || !finite_value(d0*d0) || !finite_value(d1*d1))return false;
    angular_geometry(rho,x,y,z,g.xud[0],g.xudd[0]);
    angular_geometry(sigma,x,y,z,g.xud[1],g.xudd[1]);
    for(int i=0;i<3;++i){
        if(!finite_value(g.xyz[i]))return false;
        global[i]=m.origin[i]+g.xyz[i];
        if(!finite_value(global[i]))return false;
    }
    return true;
}
}
#undef AHF_HD
#endif
