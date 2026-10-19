#ifndef AHF_EXPANSION_POINT_CUH
#define AHF_EXPANSION_POINT_CUH
#include "ahf_gpu.h"
#include <cmath>

namespace AHFinderDirect { namespace gpu { namespace point {
#define AHF_HD __host__ __device__ inline
// Five-component forward differentiation is analytic chain-rule propagation,
// not an angular finite difference. The host entry is for independent fixtures;
// production invokes this same operator only from a device kernel.
struct Dual {
    double v,d[5];
    AHF_HD Dual(double x=0):v(x),d{0,0,0,0,0}{}
    AHF_HD static Dual variable(double x,int i){Dual a(x);a.d[i]=1;return a;}
};
AHF_HD Dual operator+(Dual a,Dual b){Dual c(a.v+b.v);for(int i=0;i<5;++i)c.d[i]=a.d[i]+b.d[i];return c;}
AHF_HD Dual operator-(Dual a,Dual b){Dual c(a.v-b.v);for(int i=0;i<5;++i)c.d[i]=a.d[i]-b.d[i];return c;}
AHF_HD Dual operator*(Dual a,Dual b){Dual c(a.v*b.v);for(int i=0;i<5;++i)c.d[i]=a.d[i]*b.v+a.v*b.d[i];return c;}
AHF_HD Dual operator/(Dual a,Dual b){Dual c(a.v/b.v);for(int i=0;i<5;++i)c.d[i]=(a.d[i]-c.v*b.d[i])/b.v;return c;}
AHF_HD Dual root(Dual a){Dual c(sqrt(a.v));for(int i=0;i<5;++i)c.d[i]=a.d[i]/(2*c.v);return c;}
AHF_HD void symmetric(const double *a,double m[3][3]){
    m[0][0]=a[0];m[0][1]=m[1][0]=a[1];m[0][2]=m[2][0]=a[2];m[1][1]=a[3];m[1][2]=m[2][1]=a[4];m[2][2]=a[5];
}
AHF_HD bool inverse(const double g[3][3],double u[3][3]){
    double det=g[0][0]*(g[1][1]*g[2][2]-g[1][2]*g[1][2])-g[0][1]*(g[0][1]*g[2][2]-g[0][2]*g[1][2])+g[0][2]*(g[0][1]*g[1][2]-g[0][2]*g[1][1]);
    if(!(g[0][0]>0&&g[0][0]*g[1][1]-g[0][1]*g[0][1]>0&&det>0))return false;
    u[0][0]=(g[1][1]*g[2][2]-g[1][2]*g[1][2])/det;
    u[0][1]=u[1][0]=(g[0][2]*g[1][2]-g[0][1]*g[2][2])/det;
    u[0][2]=u[2][0]=(g[0][1]*g[1][2]-g[0][2]*g[1][1])/det;
    u[1][1]=(g[0][0]*g[2][2]-g[0][2]*g[0][2])/det;
    u[1][2]=u[2][1]=(g[0][1]*g[0][2]-g[0][0]*g[1][2])/det;
    u[2][2]=(g[0][0]*g[1][1]-g[0][1]*g[0][1])/det;return true;
}
AHF_HD bool expansion(const Point &p,double r,const double a[5],const double physical[35],double add,double &theta,double coefficients[5]){
    double g[3][3],u[3][3],k[3][3],dg[3][3][3],b[6];
    for(int j=0;j<6;++j)b[j]=physical[4*j];symmetric(b,g);
    if(!(r>0)||!inverse(g,u))return false;
    for(int j=0;j<6;++j)b[j]=physical[28+j];symmetric(b,k);
    for(int axis=0;axis<3;++axis){for(int j=0;j<6;++j)b[j]=physical[4*j+1+axis];symmetric(b,dg[axis]);}
    Dual d[5];for(int j=0;j<5;++j)d[j]=Dual::variable(a[j],j);
    Dual normal[3],H[3][3],v[3],D,B,A,C;
    for(int i=0;i<3;++i){
        normal[i]=Dual(p.direction[i])-d[0]*(p.angular_gradient[0][i]/r)-d[1]*(p.angular_gradient[1][i]/r);
        for(int j=0;j<3;++j)H[i][j]=Dual(((i==j?1.:0)-p.direction[i]*p.direction[j])/r)
            -d[0]*(p.angular_hessian[0][i][j]/(r*r))-d[1]*(p.angular_hessian[1][i][j]/(r*r))
            -d[2]*(p.angular_gradient[0][i]*p.angular_gradient[0][j]/(r*r))
            -d[3]*((p.angular_gradient[0][i]*p.angular_gradient[1][j]+p.angular_gradient[1][i]*p.angular_gradient[0][j])/(r*r))
            -d[4]*(p.angular_gradient[1][i]*p.angular_gradient[1][j]/(r*r));
    }
    double du[3][3][3]={},lng[3]={},trk=0;
    for(int axis=0;axis<3;++axis)for(int i=0;i<3;++i)for(int j=0;j<3;++j){
        lng[axis]+=.5*u[i][j]*dg[axis][i][j];
        for(int l=0;l<3;++l)for(int m=0;m<3;++m)du[axis][i][j]-=u[i][l]*dg[axis][l][m]*u[m][j];
    }
    for(int i=0;i<3;++i)for(int j=0;j<3;++j)v[i]=v[i]+normal[j]*u[i][j];
    for(int i=0;i<3;++i){
        D=D+normal[i]*v[i];B=B+v[i]*lng[i];
        for(int j=0;j<3;++j){
            B=B+normal[j]*du[i][i][j]+H[i][j]*u[i][j];
            C=C+v[i]*v[j]*k[i][j];trk+=u[i][j]*k[i][j];A=A-v[i]*v[j]*H[i][j];
            for(int l=0;l<3;++l)A=A-v[i]*normal[j]*normal[l]*(.5*du[i][j][l]);
        }
    }
    if(!(D.v>0))return false;
    Dual result=A/(D*root(D))+B/root(D)+C/D-Dual(trk)+Dual(add);
    theta=result.v;for(int j=0;j<5;++j)coefficients[j]=result.d[j];return true;
}
#undef AHF_HD
}}}
#endif
