#pragma once
#include <vector>
#include <cuda_runtime.h>
namespace AHFinderDirect { class patch_system; }
namespace ahf_diagnostics {
/** Static angle/topology metadata: deliberately contains no origin or radius. */
struct Point { int family=0,ghost=0; double first=0,second=0,sign=0,normalizer=0,dr=0,ds=0,cr=0,cs=0; };
struct GhostOperation { int destination,begin,end; };
struct Patch { int begin,end; double dr,ds; };
struct Arc { int plane,axis,begin,end; double delta; };
struct HostPlan {
 int n=0,ng=0,symmetry=0;
 std::vector<Point> points; std::vector<int> stencils,sources,arc_points;
 std::vector<double> weights,arc_coefficients;
 std::vector<GhostOperation> operations; std::vector<Patch> patches; std::vector<Arc> arcs;
};
/** Build once on a disposable original CPU patch system, fixed ORDER4 radius2. */
HostPlan make_plan(AHFinderDirect::patch_system &ps,int symmetry);
struct DeviceSurface { int molecule_radius=2,molecule_diameter=5; const int *stencils=nullptr; const Point *points=nullptr; double *ghost_h=nullptr; };
struct DevicePlan {
 int n=0,ng=0,symmetry=0,np=0,no=0,na=0;
 Point *points=nullptr; int *stencils=nullptr,*sources=nullptr,*arc_points=nullptr;
 double *weights=nullptr,*arc_coefficients=nullptr;
 GhostOperation *operations=nullptr; Patch *patches=nullptr; Arc *arcs=nullptr;
 double *ghost_h=nullptr,*work=nullptr,*integrals=nullptr;
};
/** Upload metadata and allocate reusable device scratch; no runtime field inputs.
 * DevicePlan owns allocations: do not copy the owning handle. Synchronize outstanding
 * calls before release_plan; allocation/upload/free are outside the hot-path API.
 */
DevicePlan upload_plan(const HostPlan &plan);
void release_plan(DevicePlan &plan) noexcept;
/** Asynchronous pure-device diagnostics. Inputs h[n], origin[3], physical35[n][35].
 * Outputs slots[30]: centroid xyz; quadrupole xx,xy,xz,yy,yz,zz;
 * min/max/mean radius; min/max x,y,z; circumference xy,xz,yz;
 * area,irreducible_mass,areal_radius; Px,Py,Pz,Sx,Sy,Sz.
 * physical35 is physical metric at offsets 0,4,8,12,16,20; K at 28..33;
 * trK at 34. No host transfer/fallback. Use --fmad=false for CUDA and
 * -ffp-contract=off for the original CPU oracle; no fast-math.
 * Requires finite positive h and finite positive-definite physical metric.
 * Calls for a given plan must be serialized on one stream (scratch is shared).
 */
void compute(DevicePlan &plan,const double *h,const double *origin,const double *physical35,double *slots,cudaStream_t stream=0);
/** Sequence-faithful integration from point-major xx,xy,xz,yy,yz,zz scratch.
 * Does not reconstruct moments or mutate the input scratch; remaining arithmetic
 * and P/S computation are identical to clean compute. Original NaNs retained.
 */
void compute_from_scratch(DevicePlan &plan,const double *h,const double *origin,const double *physical35,const double *scratch6,double *slots,cudaStream_t stream=0);
/** Original baseline signature; oldh[n], gauge[n][4] are runtime device inputs.
 * Exact-pole CPU NaNs intentionally retained; no production finite-pole policy.
 */
void signature(DevicePlan &plan,const double *h,const double *physical35,const double *oldh,const double *gauge,double dt,double *result,cudaStream_t stream=0);
}
