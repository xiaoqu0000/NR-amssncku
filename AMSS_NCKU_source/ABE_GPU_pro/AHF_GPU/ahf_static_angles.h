#pragma once
#include <cmath>
inline bool make_static_metadata(int family,double rho,double sigma,const double origin[3],
                                Point &out,int patch_number=0){
    if(family<0 || family>2 || !std::isfinite(rho) || !std::isfinite(sigma) || !origin)return false;
    for(int i=0;i<3;++i)if(!std::isfinite(origin[i]))return false;
    const auto signum=[](double v){return v==0.0?0.0:(v>0.0?1.0:-1.0);};
    const auto pow2=[](double v){return v*v;};
    constexpr double pi_value=3.14159265358979323846; // Original pi_value == M_pi_value, same double bits.
    Point m{};m.family=family;
    if(family==0){
        const double sign_z_via_mu=signum(std::cos(rho));
        const double y_over_z=std::tan(rho);
        const double sign_z_via_nu=signum(std::cos(sigma));
        const double x_over_z=std::tan(sigma);
        if(sign_z_via_mu==0.0 || sign_z_via_mu!=sign_z_via_nu)return false;
        const double sign_z=sign_z_via_mu;
        const double temp=1.0/std::sqrt(1.0+pow2(y_over_z)+pow2(x_over_z));
        m.first=y_over_z;m.second=x_over_z;m.sign=sign_z;m.normalizer=temp;
    }else if(family==1){
        const double nu_bar=0.5*pi_value-rho;
        const double sign_x_via_nu_bar=signum(std::cos(nu_bar));
        const double z_over_x=std::tan(nu_bar);
        const double sign_x_via_phi=signum(std::cos(sigma));
        const double y_over_x=std::tan(sigma);
        if(sign_x_via_nu_bar==0.0 || sign_x_via_nu_bar!=sign_x_via_phi)return false;
        const double sign_x=sign_x_via_nu_bar;
        const double temp=1.0/std::sqrt(1.0+pow2(z_over_x)+pow2(y_over_x));
        m.first=z_over_x;m.second=y_over_x;m.sign=sign_x;m.normalizer=temp;
    }else{
        const double mu_bar=0.5*pi_value-rho;
        const double phi_bar=0.5*pi_value-sigma;
        const double sign_y_via_mu_bar=signum(std::cos(mu_bar));
        const double z_over_y=std::tan(mu_bar);
        const double sign_y_via_phi_bar=signum(std::cos(phi_bar));
        const double x_over_y=std::tan(phi_bar);
        if(sign_y_via_mu_bar==0.0 || sign_y_via_mu_bar!=sign_y_via_phi_bar)return false;
        const double sign_y=sign_y_via_mu_bar;
        const double temp=1.0/std::sqrt(1.0+pow2(z_over_y)+pow2(x_over_y));
        m.first=z_over_y;m.second=x_over_y;m.sign=sign_y;m.normalizer=temp;
    }
    if(!std::isfinite(m.first) || !std::isfinite(m.second) || !(m.normalizer>0.) || !std::isfinite(m.normalizer))return false;
    out=m;return true;
}
