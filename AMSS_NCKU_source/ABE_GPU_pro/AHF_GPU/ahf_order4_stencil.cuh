// CPU fd_grid.h ORDER4 macros scale after the unit-spacing stencil,
// including each intermediate sigma derivative in the mixed derivative.
__device__ double cpu_sigma(const int *st,const double *h,int row,int D,int R,double inv) {
    const double p2=h[st[row*D+R+2]],p1=h[st[row*D+R+1]];
    const double m1=h[st[row*D+R-1]],m2=h[st[row*D+R-2]];
    const double sum=__dadd_rn(__dmul_rn(8.0/12.0,p1-m1),__dmul_rn(1.0/12.0,m2-p2));
    return inv*sum;
}
__device__ void angular(DeviceSurface s,int q,double a[5]) {
    const int R=s.molecule_radius,D=s.molecule_diameter;
    const int *st=s.stencils+q*D*D;
    const double *h=s.ghost_h;
    const double dr=s.points[q].dr,ds=s.points[q].ds;
    // The packed stencil and CPU fd_grid.h both use the same R=2, D=5 footprint.
    const double rp2=h[st[(R+2)*D+R]],rp1=h[st[(R+1)*D+R]];
    const double r0=h[st[R*D+R]];
    const double rm1=h[st[(R-1)*D+R]],rm2=h[st[(R-2)*D+R]];
    const double sp2=h[st[R*D+R+2]],sp1=h[st[R*D+R+1]];
    const double sm1=h[st[R*D+R-1]],sm2=h[st[R*D+R-2]];
    a[0]=dr*__dadd_rn(__dmul_rn(8.0/12.0,rp1-rm1),__dmul_rn(1.0/12.0,rm2-rp2));
    a[1]=cpu_sigma(st,h,R,D,R,ds);
    a[2]=(dr*dr)*__dsub_rn(__dadd_rn(__dmul_rn(-(30.0/12.0),r0),__dmul_rn(16.0/12.0,rm1+rp1)),__dmul_rn(1.0/12.0,rm2+rp2));
    const double mp2=cpu_sigma(st,h,R+2,D,R,ds);
    const double mp1=cpu_sigma(st,h,R+1,D,R,ds);
    const double mm1=cpu_sigma(st,h,R-1,D,R,ds);
    const double mm2=cpu_sigma(st,h,R-2,D,R,ds);
    a[3]=dr*__dadd_rn(__dmul_rn(8.0/12.0,mp1-mm1),__dmul_rn(1.0/12.0,mm2-mp2));
    a[4]=(ds*ds)*__dsub_rn(__dadd_rn(__dmul_rn(-(30.0/12.0),r0),__dmul_rn(16.0/12.0,sm1+sp1)),__dmul_rn(1.0/12.0,sm2+sp2));
}
