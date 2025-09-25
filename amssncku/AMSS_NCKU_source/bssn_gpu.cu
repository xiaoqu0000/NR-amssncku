// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <sys/time.h>
#include <cuda.h>
//#include "cutil.h"
#include <nvrtc.h>
#include <cuda_runtime.h>
using namespace std;

//includes, bssn
#include "gpu_mem.h"
#include "bssn_gpu.h"
#ifdef RESULT_CHECK
#include <fstream>
#endif

void compare_result_gpu(int ftag1,double * datac,int data_num){
	double * data = (double*)malloc(sizeof(double)*data_num);
	cudaMemcpy(data, datac, data_num * sizeof(double), cudaMemcpyDeviceToHost);
	compare_result(ftag1,data,data_num);
	free(data);
}

__global__ void test_const_address(double * testd){
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	if(_t == 0)
		testd[0] = F1o3;
}

__global__ void enforce_ga(double * trA){
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	//int ps; //TOTRY: i,j,k; double value;
	
	while(_t < _3D_SIZE[0])
	{
		  M_ gxx[_t] = M_ dxx[_t] + 1;
		  M_ gyy[_t] = M_ dyy[_t] + 1;
		  M_ gzz[_t] = M_ dzz[_t] + 1;
		// for M_ g;
		  M_ gupzz[_t] =  M_ gxx[_t] * M_ gyy[_t] * M_ gzz[_t] + M_ gxy[_t] * M_ gyz[_t] * M_ gxz[_t] + M_ gxz[_t] * M_ gxy[_t] * M_ gyz[_t] - 
		           M_ gxz[_t] * M_ gyy[_t] * M_ gxz[_t] - M_ gxy[_t] * M_ gxy[_t] * M_ gzz[_t] - M_ gxx[_t] * M_ gyz[_t] * M_ gyz[_t]; 
		
		  M_ gupzz[_t] = 1.0 / pow( M_ gupzz[_t] , F1o3 ) ;
		  
		  M_ gxx[_t] = M_ gxx[_t] * M_ gupzz[_t]; 
		  M_ gxy[_t] = M_ gxy[_t] * M_ gupzz[_t]; 
		  M_ gxz[_t] = M_ gxz[_t] * M_ gupzz[_t]; 
		  M_ gyy[_t] = M_ gyy[_t] * M_ gupzz[_t]; 
		  M_ gyz[_t] = M_ gyz[_t] * M_ gupzz[_t]; 
		  M_ gzz[_t] = M_ gzz[_t] * M_ gupzz[_t]; 
		
		  M_ dxx[_t] = M_ gxx[_t] - 1;
		  M_ dyy[_t] = M_ gyy[_t] - 1;
		  M_ dzz[_t] = M_ gzz[_t] - 1;
		// for A  ;
		
		  M_ gupxx[_t] =   ( M_ gyy[_t] * M_ gzz[_t] - M_ gyz[_t] * M_ gyz[_t] );
		  M_ gupxy[_t] = - ( M_ gxy[_t] * M_ gzz[_t] - M_ gyz[_t] * M_ gxz[_t] );
		  M_ gupxz[_t] =   ( M_ gxy[_t] * M_ gyz[_t] - M_ gyy[_t] * M_ gxz[_t] );
		  M_ gupyy[_t] =   ( M_ gxx[_t] * M_ gzz[_t] - M_ gxz[_t] * M_ gxz[_t] );
		  M_ gupyz[_t] = - ( M_ gxx[_t] * M_ gyz[_t] - M_ gxy[_t] * M_ gxz[_t] );
		  M_ gupzz[_t] =   ( M_ gxx[_t] * M_ gyy[_t] - M_ gxy[_t] * M_ gxy[_t] );
		
		  trA[_t] =         M_ gupxx[_t] *M_ Axx[_t] + M_ gupyy[_t] * M_ Ayy[_t] + M_ gupzz[_t] * M_ Azz[_t] 
		       + 2 * (M_ gupxy[_t] *M_ Axy[_t] + M_ gupxz[_t] *M_ Axz[_t] + M_ gupyz[_t] * M_ Ayz[_t]);
		
		  M_ Axx[_t] = M_ Axx[_t] - F1o3 * M_ gxx[_t] * trA[_t];
		  M_ Axy[_t] = M_ Axy[_t] - F1o3 * M_ gxy[_t] * trA[_t];
		  M_ Axz[_t] = M_ Axz[_t] - F1o3 * M_ gxz[_t] * trA[_t];
		  M_ Ayy[_t] = M_ Ayy[_t] - F1o3 * M_ gyy[_t] * trA[_t];
		  M_ Ayz[_t] = M_ Ayz[_t] - F1o3 * M_ gyz[_t] * trA[_t];
		  M_ Azz[_t] = M_ Azz[_t] - F1o3 * M_ gzz[_t] * trA[_t];
		//-------------------
		_t += STEP_SIZE;
	}
}

inline void sub_enforce_ga(int matrix_size){
	double * trA = M_ chin1;
	enforce_ga<<<GRID_DIM,BLOCK_DIM>>>(trA);
	cudaMemset(trA,0,matrix_size * sizeof(double));
	cudaThreadSynchronize(); 
	
	//cudaMemset(Mh_ gupxx,0,matrix_size * sizeof(double));
	//trA gxx,gyy,gzz gupxx,gupxy,gupxz,gupyy,gupyz,gupzz
	
}
__device__ volatile unsigned int global_count = 0;
__global__ void test_init_matrix(){
	int tid = blockIdx.x*blockDim.x+threadIdx.x;
	int curr = tid;
	while(curr < _3D_SIZE[2])
	{
		metac.fh[curr] = 0;
	  	curr += STEP_SIZE;
	}
	curr = tid;
	while(curr < _3D_SIZE[0])
	{
		metac.betaxx[curr] = 0;
		metac.betaxy[curr] = 0;
		metac.betaxz[curr] = 0;
	  	curr += STEP_SIZE;
	}
}
__global__ void init_matrix(double * mat){
	int tid = blockIdx.x*blockDim.x+threadIdx.x;
	int curr = tid;
	while(curr < _3D_SIZE[0])
	{
		mat[curr] = 0;
	  	curr += STEP_SIZE;
	}
}
__global__ void init_3_matrixs(double * mat1,double* mat2,double *mat3){
	int tid = blockIdx.x*blockDim.x+threadIdx.x;
	int curr = tid;
	while(curr < _3D_SIZE[0])
	{
		mat1[curr] = 0;
		mat2[curr] = 0;
		mat3[curr] = 0;
	  	curr += STEP_SIZE;
	}
}
__global__ void init_matrix_fh(double * mat){
	int tid = blockIdx.x*blockDim.x+threadIdx.x;
	int curr = tid;
	while(curr < _3D_SIZE[2])
	{
		mat[curr] = 0;
	  	curr += STEP_SIZE;
	}
}


__global__ void sub_symmetry_bd_partF(int ord, double * func, double *funcc)
{
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps; //TOTRY: i,j,k; double value;
	
	while(curr < _3D_SIZE[0])
	{
		int k = curr / _2D_SIZE[0];
		ps = curr - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
		int j = ps / ex_c[0];
		int i = ps  - (j * ex_c[0]);  //= ps % ex_c[0];
		
		funcc[i+ ord + (ord +j)* _1D_SIZE[ord] + (k + ord) * _2D_SIZE[ord]] = func[curr];
	
		curr += STEP_SIZE;
	}
	
}

#ifdef Vertex
__global__ void sub_symmetry_bd_partI(int ord, double * func, double * funcc,double S1){
	//for i
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps;
	int m;
	while(curr < (ex_c[1]+ord)*(ex_c[2]+ord) ){
		m =  ord * 2;
		ps = curr * _1D_SIZE[ord];
		for(int i = 0;i < ord; ++i){
			funcc[ps] = funcc [ps + m] * S1;
			ps ++;
			m -= 2;
		}
		curr+= STEP_SIZE;
	}
	__syncthreads();
}
__global__ void sub_symmetry_bd_partJ(int ord,double * func, double * funcc,double S2){
	//for j
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps;
	int m;
	
	while(curr < (ex_c[0]+ord)*(ex_c[2]+ord))
	{
		m = 2 * ord;
		ps = (curr/_1D_SIZE[ord])*_2D_SIZE[ord] + (curr % _1D_SIZE[ord]);
		for(int i = 0;i<ord;++i){
			funcc[ps] = funcc[ps + _1D_SIZE[ord] * m]* S2;
			ps += _1D_SIZE[ord];
			m -= 2;
		}
		curr+=STEP_SIZE;
	}
}
__global__ void sub_symmetry_bd_partK(int ord,double * func, double * funcc,double S3){
	//for k
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int m ;
	while(curr < _2D_SIZE[ord]){
		m = ord * 2;
		for(int i = 0;i < ord;++i){
			funcc[curr] = funcc[curr + _2D_SIZE[ord] * m] * S3;
			curr += _2D_SIZE[ord];
			m -= 2;
		}
		curr+= STEP_SIZE;
	}
	
	//TOTALLY use only 1 cycle
	/*curr = tid + _2D_SIZE[2];
	while(curr < 2 * _2D_SIZE[2]){
		funcc[curr] = funcc [curr + _2D_SIZE[2]*2] * SoA[3];
		curr+= STEP_SIZE;
	}*/

}

#else //ifdef Vertex
#ifdef Cell
__global__ void sub_symmetry_bd_partI(int ord, double * func, double * funcc,double S1){
	//for i
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps;
	int m;
	while(curr < (ex_c[1]+ord)*(ex_c[2]+ord) ){
		m =  ord * 2 - 1;  //changed from "m =  ord * 2"
		ps = curr * _1D_SIZE[ord];
		for(int i = 0;i < ord; ++i){
			funcc[ps] = funcc [ps + m] * S1;
			ps ++;
			m -= 2;
		}
		curr+= STEP_SIZE;
	}
	__syncthreads();
}

__global__ void sub_symmetry_bd_partJ(int ord,double * func, double * funcc,double S2){
	//for j
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps;
	int m;
	
	while(curr < (ex_c[0]+ord)*(ex_c[2]+ord))
	{
		m = 2 * ord - 1;
		ps = (curr/_1D_SIZE[ord])*_2D_SIZE[ord] + (curr % _1D_SIZE[ord]);
		for(int i = 0;i<ord;++i){
			funcc[ps] = funcc[ps + _1D_SIZE[ord] * m]* S2;
			ps += _1D_SIZE[ord];
			m -= 2;
		}
		curr+=STEP_SIZE;
	}
}

__global__ void sub_symmetry_bd_partK(int ord,double * func, double * funcc,double S3){
	//for k
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int m ;
	while(curr < _2D_SIZE[ord]){
		m = ord * 2 - 1;
		for(int i = 0;i < ord;++i){
			funcc[curr] = funcc[curr + _2D_SIZE[ord] * m] * S3;
			curr += _2D_SIZE[ord];
			m -= 2;
		}
		curr+= STEP_SIZE;
	}
}
#endif //ifdef Cell
#endif //ifdef Vertex
inline void sub_symmetry_bd(int ord,double * func, double * funcc,double * SoA){
	sub_symmetry_bd_partF<<<GRID_DIM,BLOCK_DIM>>>(ord,func,funcc);
	cudaThreadSynchronize();
	sub_symmetry_bd_partI<<<GRID_DIM,BLOCK_DIM>>>(ord,func,funcc,SoA[0]);
	cudaThreadSynchronize();
	sub_symmetry_bd_partJ<<<GRID_DIM,BLOCK_DIM>>>(ord,func,funcc,SoA[1]);
	cudaThreadSynchronize();
	sub_symmetry_bd_partK<<<GRID_DIM,BLOCK_DIM>>>(ord,func,funcc,SoA[2]);
	cudaThreadSynchronize();
}


__global__ void sub_fdderivs_part1(double * f,double *fh,double *fxx,double *fxy,double *fxz,double *fyy,double *fyz,double *fzz)
 {	
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps; //TOTRY: i,j,k; double value;
	
	while(curr < _3D_SIZE[0])
	{
		int k = curr / _2D_SIZE[0];
		ps = curr - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
		int j = ps / ex_c[0];
		int i = ps  - (j * ex_c[0]);
		
		if(k == ex_c[2]-1 || i == ex_c[0]-1 || j == ex_c[1]-1){
			curr += STEP_SIZE;	
			continue;
		}
		else
		{
			//xx
			if(i+2 <= ijk_max[0] && i-2 >= ijk_min[0]){
				fxx[curr] = Fdxdx*(-_FH2_(i,(j+2),(k+2))+16*_FH2_((i+1),(j+2),(k+2))-30*_FH2_((i+2),(j+2),(k+2)) 
			            -_FH2_((i+4),(j+2),(k+2))+16*_FH2_((i+3),(j+2),(k+2))		);
							
			}
			else if(i+1 <= ijk_max[0] && i-1 >= ijk_min[0]){
				fxx[curr] = Sdxdx*(_FH2_((i+1),(j+2),(k+2))-2*_FH2_((i+2),(j+2),(k+2)) 
			            +_FH2_(i+3,(j+2),(k+2))              );
			}
			//zz--
			if(k+2 <= ijk_max[2] && k-2 >= ijk_min[2]){
				fzz[curr] = Fdzdz * (-_FH2_((i+2),(j+2),k) + 16 *_FH2_((i+2),(j+2),(k+1))- 30*_FH2_((i+2),(j+2),(k+2)) 
		            -_FH2_((i+2),(j+2),(k+4))+ 16*_FH2_((i+2),(j+2),(k+3))              );
			}
			else if(k+1 <= ijk_max[2] && k-1 >= ijk_min[2]){
				fzz[curr] = Sdzdz*(_FH2_((i+2),(j+2),(k+1))- 2 * _FH2_((i+2),(j+2),(k+2)) 
			            + _FH2_((i+2),(j+2),(k+3))              );
			}

			//yy--
			if(j+2 <= ijk_max[1] && j-2 >= ijk_min[1]){
			    fyy[curr] = Fdydy*(-_FH2_((i+2),j,(k+2))+16*_FH2_((i+2),(j+1),(k+2))-30*_FH2_((i+2),(j+2),(k+2)) 
			            -_FH2_((i+2),(j+4),(k+2))+16*_FH2_((i+2),(j+3),(k+2))              );
			}
			else if(j+1 <= ijk_max[1] && j-1 >= ijk_min[1]){
				fyy[curr] = Sdydy*(_FH2_((i+2),(j+1),(k+2))-2*_FH2_((i+2),(j+2),(k+2)) 
			            +_FH2_((i+2),(j+3),(k+2))              );
			}

			
			
			//xy
		    if(i+2 <= ijk_max[0] && i-2 >= ijk_min[0] && j+2 <= ijk_max[1] && j-2 >= ijk_min[1])	
		   		fxy[curr] = Fdxdy*((_FH2_(i,j,(k+2))-8*_FH2_((i+1),j,(k+2))+8*_FH2_((i+3),j,(k+2))-_FH2_((i+4),j,(k+2))) 
		                       -8 *(_FH2_(i,(j+1),(k+2))-8*_FH2_((i+1),(j+1),(k+2))+8*_FH2_((i+3),(j+1),(k+2))-_FH2_((i+4),(j+1),(k+2)))  
		                       +8 *(_FH2_(i,(j+3),(k+2))-8*_FH2_((i+1),(j+3),(k+2))+8*_FH2_((i+3),(j+3),(k+2))-_FH2_((i+4),(j+3),(k+2)))  
		                       -    (_FH2_(i,(j+4),(k+2))-8*_FH2_((i+1),(j+4),(k+2))+8*_FH2_((i+3),(j+4),(k+2))-_FH2_((i+4),(j+4),(k+2))));

		   	else if(i+1 <= ijk_max[0] && i-1 >= ijk_min[0] && j+1 <= ijk_max[1] && j-1 >= ijk_min[1])
		                
		   		fxy[curr] = Sdxdy*(_FH2_((i+1),(j+1),(k+2))-_FH2_((i+3),(j+1),(k+2))-_FH2_((i+1),(j+3),(k+2))+_FH2_((i+3),(j+3),(k+2)));
			//xz
		    if(i+2 <= ijk_max[0] && i-2 >= ijk_min[0] && k+2 <= ijk_max[2] && k-2 >= ijk_min[2])
		   		fxz[curr] = Fdxdz*(     (_FH2_(i,(j+2),k)-8*_FH2_((i+1),(j+2),k)+8*_FH2_((i+3),(j+2),k)-_FH2_((i+4),(j+2),k)) 
		                       -8 *(_FH2_(i,(j+2),(k+1))-8*_FH2_((i+1),(j+2),(k+1))+8*_FH2_((i+3),(j+2),(k+1))-_FH2_((i+4),(j+2),(k+1))) 
		                       +8 *(_FH2_(i,(j+2),(k+3))-8*_FH2_((i+1),(j+2),(k+3))+8*_FH2_((i+3),(j+2),(k+3))-_FH2_((i+4),(j+2),(k+3)))  
		                       -    (_FH2_(i,(j+2),(k+4))-8*_FH2_((i+1),(j+2),(k+4))+8*_FH2_((i+3),(j+2),(k+4))-_FH2_((i+4),(j+2),(k+4))));
		                       
			else if(i+1 <= ijk_max[0] && i-1 >= ijk_min[0] && k+1 <= ijk_max[2] && k-1 >= ijk_min[2])
					fxz[curr] = Sdxdz*(_FH2_((i+1),(j+2),(k+1))-_FH2_((i+3),(j+2),(k+1))-_FH2_((i+1),(j+2),(k+3))+_FH2_((i+3),(j+2),(k+3)));
			//yz
			if(j+2 <= ijk_max[1] && j-2 >= ijk_min[1] && k+2 <= ijk_max[2] && k-2 >= ijk_min[2])
					fyz[curr] = Fdydz*(     (_FH2_((i+2),j,k)-8*_FH2_((i+2),(j+1),k)+8*_FH2_((i+2),(j+3),k)-_FH2_((i+2),(j+4),k))  
		                       -8 *(_FH2_((i+2),j,(k+1))-8*_FH2_((i+2),(j+1),(k+1))+8*_FH2_((i+2),(j+3),(k+1))-_FH2_((i+2),(j+4),(k+1)))  
		                       +8 *(_FH2_((i+2),j,(k+3))-8*_FH2_((i+2),(j+1),(k+3))+8*_FH2_((i+2),(j+3),(k+3))-_FH2_((i+2),(j+4),(k+3)))  
		                       -    (_FH2_((i+2),j,(k+4))-8*_FH2_((i+2),(j+1),(k+4))+8*_FH2_((i+2),(j+3),(k+4))-_FH2_((i+2),(j+4),(k+4))));
		                       
			else if(j+1 <= ijk_max[1] && j-1 >= ijk_min[1] && k+1 <= ijk_max[2] && k-1 >= ijk_min[2])
					fyz[curr] = Sdydz*(_FH2_((i+2),(j+1),(k+1))-_FH2_((i+2),(j+3),(k+1))-_FH2_((i+2),(j+1),(k+3))+_FH2_((i+2),(j+3),(k+3)));
			
			curr += STEP_SIZE;
		}
	}
 	
 	__syncthreads();
 }

inline void sub_fdderivs(double * f,double *fh,double *fxx,double *fxy,double *fxz,double *fyy,double *fyz,double *fzz,double* SoA)
{
	sub_symmetry_bd(2,f,fh,SoA);
	cudaMemset(fxx,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fxy,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fxz,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fyy,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fyz,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fzz,0,_3D_SIZE[0] * sizeof(double));
	cudaThreadSynchronize(); 
	sub_fdderivs_part1<<<GRID_DIM,BLOCK_DIM>>>(f,fh,fxx,fxy,fxz,fyy,fyz,fzz);
	cudaThreadSynchronize(); 
}

__global__ void sub_fderivs_part1(double * f,double * fh,double *fx,double *fy,double *fz  )
 {	
	int curr = blockIdx.x*blockDim.x+threadIdx.x;
	int ps; //TOTRY: i,j,k; double value;
	
	while(curr < _3D_SIZE[0])
	{
		int k = curr / _2D_SIZE[0];
		ps = curr - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
		int j = ps / ex_c[0];
		int i = ps  - (j * ex_c[0]);
		
		if(k == ex_c[2]-1 || i == ex_c[0]-1 || j == ex_c[1]-1){
			curr += STEP_SIZE;	
			continue;
		}
		
			//X--
			if(i+2 <= ijk_max[0] && i-2 >= ijk_min[0])
				fx[curr] = d12dxyz[0]*(fh[i+(j+2)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]] - 
								8*fh[i+1+(j+2)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	+
								8*fh[i+3+(j+2)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	-
								fh[i+4+(j+2)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	);
							
			else if(i+1 <= ijk_max[0] && i-1 >= ijk_min[0])
				fx[curr] = d2dxyz[0]*(-fh[i+1+(j+2)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	+
								fh[i+3+(j+2)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	);
			//Y--
			if(j+2 <= ijk_max[1] && j-2 >= ijk_min[1])
	      		fy[curr]=d12dxyz[1]*(fh[i+2+j*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]-
	      						8*fh[i+2+(j+1)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	+
	      						8*fh[i+2+(j+3)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]	-
	      						fh[i+2+(j+4)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]);
	      						
	    	else if(j+1 <= ijk_max[1] && j-1 >= ijk_min[1])
	     		fy[curr]=d2dxyz[1]*(-fh[i+2+(j+1)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]		+
	     						fh[i+2+(j+3)*_1D_SIZE[2]+(k+2)*_2D_SIZE[2]]);
	     	//Z--
      
	     	if(k+2 <= ijk_max[2] && k-2 >= ijk_min[2])
	      	   fz[curr]=d12dxyz[2]*(	fh[i+2+(j+2)*_1D_SIZE[2]+k    *_2D_SIZE[2]]		-
	      						 8*		fh[i+2+(j+2)*_1D_SIZE[2]+(k+1)*_2D_SIZE[2]]		+
	      						 8*		fh[i+2+(j+2)*_1D_SIZE[2]+(k+3)*_2D_SIZE[2]]		-
	      								fh[i+2+(j+2)*_1D_SIZE[2]+(k+4)*_2D_SIZE[2]]);
	      					
	    	else if(k+1 <= ijk_max[2] && k-1 >= ijk_min[2])
	      		fz[curr]=d2dxyz[2]*(-fh[i+2+(j+2)*_1D_SIZE[2]+(k+1)*_2D_SIZE[2]]+
	      						fh[i+2+(j+2)*_1D_SIZE[2]+(k+3)*_2D_SIZE[2]]);
      
			curr += STEP_SIZE;
		
	}
 }
 
inline void sub_fderivs(double * f,double * fh,double *fx,double *fy,double *fz,double * SoA)
{
	sub_symmetry_bd(2,f,fh,SoA);
	
	cudaMemset(fx,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fy,0,_3D_SIZE[0] * sizeof(double));
	cudaMemset(fz,0,_3D_SIZE[0] * sizeof(double));
	
	cudaThreadSynchronize(); 
	sub_fderivs_part1<<<GRID_DIM,BLOCK_DIM>>>(f,fh,fx,fy,fz);
	cudaThreadSynchronize();
}

__global__ void computeRicci_part1(double * dst)
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{	
		dst[_t] = M_ gupxx [_t]* M_ fxx [_t]+ M_ gupyy[_t]* M_ fyy[_t]+ M_ gupzz[_t]* M_ fzz[_t]+ 
         ( M_ gupxy[_t]* M_ fxy[_t]+ M_ gupxz[_t]* M_ fxz[_t]+ M_ gupyz[_t]* M_ fyz[_t]) * 2;
	
		_t += STEP_SIZE;
	}
}

 inline void computeRicci(double * src,double* dst,double * SoA, Meta* meta)
{
	sub_fdderivs(src,Mh_ fh,Mh_ fxx,Mh_ fxy,Mh_ fxz,Mh_ fyy,Mh_ fyz,Mh_ fzz,SoA);
	cudaThreadSynchronize();
	computeRicci_part1<<<GRID_DIM,BLOCK_DIM>>>(dst);
	cudaThreadSynchronize();
	
}/*Exception*/

__global__ void sub_kodis_part1(double *f,double *fh,double *f_rhs)
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	int ps; //TOTRY: i,j,k; double value;
	double inc_f_rhs;
	while(_t < _3D_SIZE[0])
	{
		int k = _t / _2D_SIZE[0];
		ps = _t - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
		int j = ps / ex_c[0];
		int i = ps  - (j * ex_c[0]);
		
		if(k == ex_c[2]-1 && i == ex_c[0]-1 && j == ex_c[1]-1){
			_t += STEP_SIZE;
			continue;
		}

		if(i-3 >= ijk_min3[0] && i+3 <= ijk_max[0] && 
		    j-3 >= ijk_min3[1] && j+3 <= ijk_max[1] && 
		    k-3 >= ijk_min3[2] && k+3 <= ijk_max[2])
		{ 
		// x direction
		   inc_f_rhs =  		( (_FH3_(i,(j+3),(k+3))+_FH3_((i+6),(j+3),(k+3))) - 
		                6*(_FH3_((i+1),(j+3),(k+3))+_FH3_((i+5),(j+3),(k+3))) + 
		                15*(_FH3_((i+2),(j+3),(k+3))+_FH3_((i+4),(j+3),(k+3))) - 
		                20* _FH3_((i+3),(j+3),(k+3))          ) /dX;
		                
       					
		// y direction

	   		inc_f_rhs +=		(  (_FH3_((i+3),j,(k+3))+_FH3_((i+3),(j+6),(k+3))) - 
		                6*(_FH3_((i+3),(j+1),(k+3))+_FH3_((i+3),(j+5),(k+3))) + 
		                15*(_FH3_((i+3),(j+2),(k+3))+_FH3_((i+3),(j+4),(k+3))) - 
		                20* _FH3_((i+3),(j+3),(k+3))             )/dY;
		                
		// z direction
					
		   	inc_f_rhs +=		 ( (_FH3_((i+3),(j+3),k)+_FH3_((i+3),(j+3),(k+6))) - 
		            6*(_FH3_((i+3),(j+3),(k+1))+_FH3_((i+3),(j+3),(k+5))) + 
		            15*(_FH3_((i+3),(j+3),(k+2))+_FH3_((i+3),(j+3),(k+4))) - 
		            20* _FH3_((i+3),(j+3),(k+3))             )/dZ;
		   inc_f_rhs *= eps_c;
		   inc_f_rhs /= 64;         
		   f_rhs[_t] += inc_f_rhs;  //be careful the mark is "+=" not "==" !
		}
	   
		_t += STEP_SIZE;
	}
}

inline void sub_kodis(double *f,double *fh,double *f_rhs,double *SoA)
{
	sub_symmetry_bd(3,f,fh,SoA);
	cudaThreadSynchronize();
	sub_kodis_part1<<<GRID_DIM,BLOCK_DIM>>>(f,fh,f_rhs);
	cudaThreadSynchronize();
}
 
__global__ void  sub_lopsided_part1(double *f,double* fh,double *f_rhs,double *Sfx,double *Sfy,double *Sfz)
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	int ps; //TOTRY: i,j,k; double value;
	
	while(_t < _3D_SIZE[0])
	{
		int k = _t / _2D_SIZE[0];
		ps = _t - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
		int j = ps / ex_c[0];
		int i = ps  - (j * ex_c[0]);
		
		if(k < ex_c[2]-1 && i < ex_c[0]-1 && j < ex_c[1]-1){
			// x direction   
		    if(Sfx[_t] >= 0 && i+3 <= ijk_max[0] && i-1 >= ijk_min2[0])
		     f_rhs[_t]=f_rhs[_t]+                                                   
		                  Sfx[_t]*d12dxyz[0]*(-3*_FH3_((i+2),(j+3),(k+3))-10*_FH3_((i+3),(j+3),(k+3))+18*_FH3_((i+4),(j+3),(k+3)) 
		                                    -6*_FH3_((i+5),(j+3),(k+3))+    _FH3_((i+6),(j+3),(k+3)));

		     else if(Sfx[_t] <= 0 && i-3 >= ijk_min2[0] && i+1 <= ijk_max[0])
		     f_rhs[_t]=f_rhs[_t]-                                                   
		                  Sfx[_t]*d12dxyz[0]*(-3*_FH3_((i+4),(j+3),(k+3))-10*_FH3_((i+3),(j+3),(k+3))+18*_FH3_((i+2),(j+3),(k+3)) 
		                                    -6*_FH3_((i+1),(j+3),(k+3))+    _FH3_(i,(j+3),(k+3)));

		     else if(i+2 <= ijk_max[0] && i-2 >= ijk_min2[0])


		     f_rhs[_t]=f_rhs[_t]+                                                           
		                  Sfx[_t]*d12dxyz[0]*(_FH3_((i+1),(j+3),(k+3))-8*_FH3_((i+2),(j+3),(k+3))+8*_FH3_((i+4),(j+3),(k+3))-_FH3_((i+5),(j+3),(k+3)));

		     else if(i+1 <= ijk_max[0] && i-1 >= ijk_min2[0])

		     f_rhs[_t]=f_rhs[_t] + Sfx[_t]*d2dxyz[0]*(-_FH3_((i+2),(j+3),(k+3))+_FH3_((i+4),(j+3),(k+3)));


			// y direction   
		    if(Sfy[_t] >= 0 && j+3 <= ijk_max[1] && j-1 >= ijk_min2[1])

		     f_rhs[_t]=f_rhs[_t]+                                                   
		                  Sfy[_t]*d12dxyz[1]*(-3*_FH3_((i+3),(j+2),(k+3))-10*_FH3_((i+3),(j+3),(k+3))+18*_FH3_((i+3),(j+4),(k+3)) 
		                                    -6*_FH3_((i+3),(j+5),(k+3))+    _FH3_((i+3),(j+6),(k+3)));

		    else if(Sfy[_t] <= 0 && j-3 >= ijk_min2[1] && j+1 <= ijk_max[1])
		     f_rhs[_t]=f_rhs[_t]-                                                   
		                  Sfy[_t]*d12dxyz[1]*(-3*_FH3_((i+3),(j+4),(k+3))-10*_FH3_((i+3),(j+3),(k+3))+18*_FH3_((i+3),(j+2),(k+3)) 
		                                    -6*_FH3_((i+3),(j+1),(k+3))+    _FH3_((i+3),j,(k+3)));

		    else if(j+2 <= ijk_max[1] && j-2 >= ijk_min2[1])

		     f_rhs[_t]=f_rhs[_t]+                                                            
		                  Sfy[_t]*d12dxyz[1]*(_FH3_((i+3),(j+1),(k+3))-8*_FH3_((i+3),(j+2),(k+3))+8*_FH3_((i+3),(j+4),(k+3))-_FH3_((i+3),(j+5),(k+3)));

		    else if(j+1 <= ijk_max[1] && j-1 >= ijk_min2[1])

		     f_rhs[_t]=f_rhs[_t] + Sfy[_t]*d2dxyz[1]*(-_FH3_((i+3),(j+2),(k+3))+_FH3_((i+3),(j+4),(k+3)));
		     

			// z direction   
		    if(Sfz[_t] >= 0 && k+3 <= ijk_max[2] && k-1 >= ijk_min2[2])
			//         v
			// D f = ------[ - 3f    - 10f  + 18f    - 6f     + f     ]
			//  i     12dx       i-v      i      i+v     i+2v    i+3v
		     f_rhs[_t]=f_rhs[_t]+                                                   
		                  Sfz[_t]*d12dxyz[2]*(-3*_FH3_((i+3),(j+3),(k+2))-10*_FH3_((i+3),(j+3),(k+3))+18*_FH3_((i+3),(j+3),(k+4)) 
		                                    -6*_FH3_((i+3),(j+3),(k+5))+    _FH3_((i+3),(j+3),(k+6)));

		    else if(Sfz[_t] <= 0 && k-3 >= ijk_min2[2] && k+1 <= ijk_max[2])
		     f_rhs[_t]=f_rhs[_t]-                                                   
		                  Sfz[_t]*d12dxyz[2]*(-3*_FH3_((i+3),(j+3),(k+4))-10*_FH3_((i+3),(j+3),(k+3))+18*_FH3_((i+3),(j+3),(k+2)) 
		                                    -6*_FH3_((i+3),(j+3),(k+1))+    _FH3_((i+3),(j+3),k));

		     else if(k+2 <= ijk_max[2] && k-2 >= ijk_min2[2])

		     f_rhs[_t]=f_rhs[_t]+                                                            
		                  Sfz[_t]*d12dxyz[2]*(_FH3_((i+3),(j+3),(k+1))-8*_FH3_((i+3),(j+3),(k+2))+8*_FH3_((i+3),(j+3),(k+4))-_FH3_((i+3),(j+3),(k+5)));

		     else if(k+1 <= ijk_max[2] && k-1 >= ijk_min2[2])

		     f_rhs[_t]=f_rhs[_t]+Sfz[_t]*d2dxyz[2]*(-_FH3_((i+3),(j+3),(k+2))+_FH3_((i+3),(j+3),(k+4)));
		}
		//-------------------
		_t += STEP_SIZE;
	}
}


inline void  sub_lopsided(double *f,double*fh,double *f_rhs,double *Sfx,double *Sfy,double *Sfz,double *SoA){
	sub_symmetry_bd(3,f,fh,SoA);
	cudaThreadSynchronize(); 
	sub_lopsided_part1<<<GRID_DIM,BLOCK_DIM>>>(f,fh,f_rhs,Sfx,Sfy,Sfz);
	cudaThreadSynchronize(); 
}

__global__ void compute_rhs_bssn_part1() 
{
	int tid = blockIdx.x*blockDim.x+threadIdx.x;
	int curr = tid;
	while(curr < _3D_SIZE[0])
	{
		metac.alpn1[curr] = metac.Lap[curr] + 1; 
	  	metac.chin1[curr] = metac.chi[curr] + 1;
	  	metac.gxx[curr] = metac.dxx[curr] + 1;
	  	metac.gyy[curr] = metac.dyy[curr] + 1;
	  	metac.gzz[curr] = metac.dzz[curr] + 1;
	  	
	  	curr += STEP_SIZE;
	}
}

__global__ void compute_rhs_bssn_part2() 
{
	//__shared__ int judge = 1;
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{

		M_ div_beta[_t] = M_ betaxx[_t] + M_ betayy[_t] + M_ betazz[_t];
		M_ chi_rhs[_t] = F2o3 *M_ chin1[_t]*( M_ alpn1[_t] * M_ trK[_t] - M_ div_beta[_t] ); //rhs[_t] for M_ chi

		M_ gxx_rhs[_t] = - 2 * M_ alpn1[_t] * M_ Axx[_t]   -  F2o3 * M_ gxx[_t]* M_ div_beta[_t]        + 
		            2 *(  M_ gxx[_t]* M_ betaxx[_t]+   M_ gxy[_t]* M_ betayx[_t]+   M_ gxz[_t]* M_ betazx[_t]);
		M_ gyy_rhs[_t] = - 2 * M_ alpn1[_t] * M_ Ayy[_t]   -  F2o3 * M_ gyy[_t]* M_ div_beta[_t]        + 
		            2 *(  M_ gxy[_t]* M_ betaxy[_t]+   M_ gyy[_t]* M_ betayy[_t]+   M_ gyz[_t]* M_ betazy[_t]);

		M_ gzz_rhs[_t] = - 2 * M_ alpn1[_t] * M_ Azz[_t]   -  F2o3 * M_ gzz[_t]* M_ div_beta[_t]        + 
		            2 *(  M_ gxz[_t]* M_ betaxz[_t]+   M_ gyz[_t]* M_ betayz[_t]+   M_ gzz[_t]* M_ betazz[_t]);

		M_ gxy_rhs[_t] = - 2 * M_ alpn1[_t] * M_ Axy[_t]   +  F1o3 * M_ gxy[_t]   * M_ div_beta[_t]     + 
		                    M_ gxx[_t]* M_ betaxy[_t]                 +   M_ gxz[_t]* M_ betazy[_t]+ 
		                                    M_ gyy[_t]* M_ betayx[_t]+   M_ gyz[_t]* M_ betazx[_t]  
		                                                -   M_ gxy[_t]* M_ betazz[_t];

		M_ gyz_rhs[_t] = - 2 * M_ alpn1[_t] * M_ Ayz[_t]   +  F1o3 * M_ gyz[_t]   * M_ div_beta[_t]     + 
		                    M_ gxy[_t]* M_ betaxz[_t]+   M_ gyy[_t]* M_ betayz[_t]                 + 
		                    M_ gxz[_t]* M_ betaxy[_t]                 +   M_ gzz[_t]* M_ betazy[_t]  
		                                                -   M_ gyz[_t]* M_ betaxx[_t];

		M_ gxz_rhs[_t] = - 2 * M_ alpn1[_t] * M_ Axz[_t]   +  F1o3 * M_ gxz[_t]   * M_ div_beta[_t]     + 
		                    M_ gxx[_t]* M_ betaxz[_t]+   M_ gxy[_t]* M_ betayz[_t]                 + 
		                                    M_ gyz[_t]* M_ betayx[_t]+   M_ gzz[_t]* M_ betazx[_t]  
		                                                -   M_ gxz[_t]* M_ betayy[_t];    //rhs[_t] for gij

		// invert tilted metric
		M_ gupzz[_t]=  M_ gxx[_t]* M_ gyy[_t]* M_ gzz[_t]+ M_ gxy[_t]* M_ gyz[_t]* M_ gxz[_t]+ M_ gxz[_t]* M_ gxy[_t]* M_ gyz[_t]- 
		        M_ gxz[_t]* M_ gyy[_t]* M_ gxz[_t]- M_ gxy[_t]* M_ gxy[_t]* M_ gzz[_t]- M_ gxx[_t]* M_ gyz[_t]* M_ gyz[_t];
		M_ gupxx[_t]=   ( M_ gyy[_t]* M_ gzz[_t]- M_ gyz[_t]* M_ gyz[_t]) / M_ gupzz[_t];
		M_ gupxy[_t]= - ( M_ gxy[_t]* M_ gzz[_t]- M_ gyz[_t]* M_ gxz[_t]) / M_ gupzz[_t];
		M_ gupxz[_t]=   ( M_ gxy[_t]* M_ gyz[_t]- M_ gyy[_t]* M_ gxz[_t]) / M_ gupzz[_t];
		M_ gupyy[_t]=   ( M_ gxx[_t]* M_ gzz[_t]- M_ gxz[_t]* M_ gxz[_t]) / M_ gupzz[_t];
		M_ gupyz[_t]= - ( M_ gxx[_t]* M_ gyz[_t]- M_ gxy[_t]* M_ gxz[_t]) / M_ gupzz[_t];
		M_ gupzz[_t]=   ( M_ gxx[_t]* M_ gyy[_t]- M_ gxy[_t]* M_ gxy[_t]) / M_ gupzz[_t];
		//if(threadIdx.x == 0){
		//	judge = co_c;
		//}
		//__syncthreads();

		if(co_c == 0)
		{
		// M_ Gam^i_Res = M_ Gam^i + M_ gup^ij_,j
		M_ Gmx_Res[_t] = M_ Gamx[_t] - (M_ gupxx[_t]*(M_ gupxx[_t]*M_ gxxx[_t]+M_ gupxy[_t]*M_ gxyx[_t]+M_ gupxz[_t]*M_ gxzx[_t])
		                +M_ gupxy[_t]*(M_ gupxx[_t]*M_ gxyx[_t]+M_ gupxy[_t]*M_ gyyx[_t]+M_ gupxz[_t]*M_ gyzx[_t])
		                +M_ gupxz[_t]*(M_ gupxx[_t]*M_ gxzx[_t]+M_ gupxy[_t]*M_ gyzx[_t]+M_ gupxz[_t]*M_ gzzx[_t])
		                +M_ gupxx[_t]*(M_ gupxy[_t]*M_ gxxy[_t]+M_ gupyy[_t]*M_ gxyy[_t]+M_ gupyz[_t]*M_ gxzy[_t])
		                +M_ gupxy[_t]*(M_ gupxy[_t]*M_ gxyy[_t]+M_ gupyy[_t]*M_ gyyy[_t]+M_ gupyz[_t]*M_ gyzy[_t])
		                +M_ gupxz[_t]*(M_ gupxy[_t]*M_ gxzy[_t]+M_ gupyy[_t]*M_ gyzy[_t]+M_ gupyz[_t]*M_ gzzy[_t])
		                +M_ gupxx[_t]*(M_ gupxz[_t]*M_ gxxz[_t]+M_ gupyz[_t]*M_ gxyz[_t]+M_ gupzz[_t]*M_ gxzz[_t])
		                +M_ gupxy[_t]*(M_ gupxz[_t]*M_ gxyz[_t]+M_ gupyz[_t]*M_ gyyz[_t]+M_ gupzz[_t]*M_ gyzz[_t])
		                +M_ gupxz[_t]*(M_ gupxz[_t]*M_ gxzz[_t]+M_ gupyz[_t]*M_ gyzz[_t]+M_ gupzz[_t]*M_ gzzz[_t]));
		M_ Gmy_Res[_t] = M_ Gamy[_t] - (M_ gupxx[_t]*(M_ gupxy[_t]*M_ gxxx[_t]+M_ gupyy[_t]*M_ gxyx[_t]+M_ gupyz[_t]*M_ gxzx[_t])
		                +M_ gupxy[_t]*(M_ gupxy[_t]*M_ gxyx[_t]+M_ gupyy[_t]*M_ gyyx[_t]+M_ gupyz[_t]*M_ gyzx[_t])
		                +M_ gupxz[_t]*(M_ gupxy[_t]*M_ gxzx[_t]+M_ gupyy[_t]*M_ gyzx[_t]+M_ gupyz[_t]*M_ gzzx[_t])
		                +M_ gupxy[_t]*(M_ gupxy[_t]*M_ gxxy[_t]+M_ gupyy[_t]*M_ gxyy[_t]+M_ gupyz[_t]*M_ gxzy[_t])
		                +M_ gupyy[_t]*(M_ gupxy[_t]*M_ gxyy[_t]+M_ gupyy[_t]*M_ gyyy[_t]+M_ gupyz[_t]*M_ gyzy[_t])
		                +M_ gupyz[_t]*(M_ gupxy[_t]*M_ gxzy[_t]+M_ gupyy[_t]*M_ gyzy[_t]+M_ gupyz[_t]*M_ gzzy[_t])
		                +M_ gupxy[_t]*(M_ gupxz[_t]*M_ gxxz[_t]+M_ gupyz[_t]*M_ gxyz[_t]+M_ gupzz[_t]*M_ gxzz[_t])
		                +M_ gupyy[_t]*(M_ gupxz[_t]*M_ gxyz[_t]+M_ gupyz[_t]*M_ gyyz[_t]+M_ gupzz[_t]*M_ gyzz[_t])
		                +M_ gupyz[_t]*(M_ gupxz[_t]*M_ gxzz[_t]+M_ gupyz[_t]*M_ gyzz[_t]+M_ gupzz[_t]*M_ gzzz[_t]));
		M_ Gmz_Res[_t] = M_ Gamz[_t] - (M_ gupxx[_t]*(M_ gupxz[_t]*M_ gxxx[_t]+M_ gupyz[_t]*M_ gxyx[_t]+M_ gupzz[_t]*M_ gxzx[_t])
		                +M_ gupxy[_t]*(M_ gupxz[_t]*M_ gxyx[_t]+M_ gupyz[_t]*M_ gyyx[_t]+M_ gupzz[_t]*M_ gyzx[_t])
		                +M_ gupxz[_t]*(M_ gupxz[_t]*M_ gxzx[_t]+M_ gupyz[_t]*M_ gyzx[_t]+M_ gupzz[_t]*M_ gzzx[_t])
		                +M_ gupxy[_t]*(M_ gupxz[_t]*M_ gxxy[_t]+M_ gupyz[_t]*M_ gxyy[_t]+M_ gupzz[_t]*M_ gxzy[_t])
		                +M_ gupyy[_t]*(M_ gupxz[_t]*M_ gxyy[_t]+M_ gupyz[_t]*M_ gyyy[_t]+M_ gupzz[_t]*M_ gyzy[_t])
		                +M_ gupyz[_t]*(M_ gupxz[_t]*M_ gxzy[_t]+M_ gupyz[_t]*M_ gyzy[_t]+M_ gupzz[_t]*M_ gzzy[_t])
		                +M_ gupxz[_t]*(M_ gupxz[_t]*M_ gxxz[_t]+M_ gupyz[_t]*M_ gxyz[_t]+M_ gupzz[_t]*M_ gxzz[_t])
		                +M_ gupyz[_t]*(M_ gupxz[_t]*M_ gxyz[_t]+M_ gupyz[_t]*M_ gyyz[_t]+M_ gupzz[_t]*M_ gyzz[_t])
		                +M_ gupzz[_t]*(M_ gupxz[_t]*M_ gxzz[_t]+M_ gupyz[_t]*M_ gyzz[_t]+M_ gupzz[_t]*M_ gzzz[_t]));
		}//if(co == 0)

		// second kind of connection
		M_ Gamxxx[_t]=HALF*( M_ gupxx[_t]*M_ gxxx[_t]+ M_ gupxy[_t]*(2*M_ gxyx[_t]- M_ gxxy[_t]) + M_ gupxz[_t]*(2*M_ gxzx[_t]- M_ gxxz[_t]));
		M_ Gamyxx[_t]=HALF*( M_ gupxy[_t]*M_ gxxx[_t]+ M_ gupyy[_t]*(2*M_ gxyx[_t]- M_ gxxy[_t]) + M_ gupyz[_t]*(2*M_ gxzx[_t]- M_ gxxz[_t]));
		M_ Gamzxx[_t]=HALF*( M_ gupxz[_t]*M_ gxxx[_t]+ M_ gupyz[_t]*(2*M_ gxyx[_t]- M_ gxxy[_t]) + M_ gupzz[_t]*(2*M_ gxzx[_t]- M_ gxxz[_t]));

		M_ Gamxyy[_t]=HALF*( M_ gupxx[_t]*(2*M_ gxyy[_t]- M_ gyyx[_t]) + M_ gupxy[_t]*M_ gyyy[_t]+ M_ gupxz[_t]*(2*M_ gyzy[_t]- M_ gyyz[_t]));
		M_ Gamyyy[_t]=HALF*( M_ gupxy[_t]*(2*M_ gxyy[_t]- M_ gyyx[_t]) + M_ gupyy[_t]*M_ gyyy[_t]+ M_ gupyz[_t]*(2*M_ gyzy[_t]- M_ gyyz[_t]));
		M_ Gamzyy[_t]=HALF*( M_ gupxz[_t]*(2*M_ gxyy[_t]- M_ gyyx[_t]) + M_ gupyz[_t]*M_ gyyy[_t]+ M_ gupzz[_t]*(2*M_ gyzy[_t]- M_ gyyz[_t]));

		M_ Gamxzz[_t]=HALF*( M_ gupxx[_t]*(2*M_ gxzz[_t]- M_ gzzx[_t]) + M_ gupxy[_t]*(2*M_ gyzz[_t]- M_ gzzy[_t]) + M_ gupxz[_t]*M_ gzzz[_t]);
		M_ Gamyzz[_t]=HALF*( M_ gupxy[_t]*(2*M_ gxzz[_t]- M_ gzzx[_t]) + M_ gupyy[_t]*(2*M_ gyzz[_t]- M_ gzzy[_t]) + M_ gupyz[_t]*M_ gzzz[_t]);
		M_ Gamzzz[_t]=HALF*( M_ gupxz[_t]*(2*M_ gxzz[_t]- M_ gzzx[_t]) + M_ gupyz[_t]*(2*M_ gyzz[_t]- M_ gzzy[_t]) + M_ gupzz[_t]*M_ gzzz[_t]);

		M_ Gamxxy[_t]=HALF*( M_ gupxx[_t]*M_ gxxy[_t]+ M_ gupxy[_t]*M_ gyyx[_t]+ M_ gupxz[_t]*( M_ gxzy[_t]+ M_ gyzx[_t]- M_ gxyz[_t]) );
		M_ Gamyxy[_t]=HALF*( M_ gupxy[_t]*M_ gxxy[_t]+ M_ gupyy[_t]*M_ gyyx[_t]+ M_ gupyz[_t]*( M_ gxzy[_t]+ M_ gyzx[_t]- M_ gxyz[_t]) );
		M_ Gamzxy[_t]=HALF*( M_ gupxz[_t]*M_ gxxy[_t]+ M_ gupyz[_t]*M_ gyyx[_t]+ M_ gupzz[_t]*( M_ gxzy[_t]+ M_ gyzx[_t]- M_ gxyz[_t]) );

		M_ Gamxxz[_t]=HALF*( M_ gupxx[_t]*M_ gxxz[_t]+ M_ gupxy[_t]*( M_ gxyz[_t]+ M_ gyzx[_t]- M_ gxzy[_t]) + M_ gupxz[_t]*M_ gzzx[_t]);
		M_ Gamyxz[_t]=HALF*( M_ gupxy[_t]*M_ gxxz[_t]+ M_ gupyy[_t]*( M_ gxyz[_t]+ M_ gyzx[_t]- M_ gxzy[_t]) + M_ gupyz[_t]*M_ gzzx[_t]);
		M_ Gamzxz[_t]=HALF*( M_ gupxz[_t]*M_ gxxz[_t]+ M_ gupyz[_t]*( M_ gxyz[_t]+ M_ gyzx[_t]- M_ gxzy[_t]) + M_ gupzz[_t]*M_ gzzx[_t]);

		M_ Gamxyz[_t]=HALF*( M_ gupxx[_t]*( M_ gxyz[_t]+ M_ gxzy[_t]- M_ gyzx[_t]) + M_ gupxy[_t]*M_ gyyz[_t]+ M_ gupxz[_t]*M_ gzzy[_t]);
		M_ Gamyyz[_t]=HALF*( M_ gupxy[_t]*( M_ gxyz[_t]+ M_ gxzy[_t]- M_ gyzx[_t]) + M_ gupyy[_t]*M_ gyyz[_t]+ M_ gupyz[_t]*M_ gzzy[_t]);
		M_ Gamzyz[_t]=HALF*( M_ gupxz[_t]*( M_ gxyz[_t]+ M_ gxzy[_t]- M_ gyzx[_t]) + M_ gupyz[_t]*M_ gyyz[_t]+ M_ gupzz[_t]*M_ gzzy[_t]);
		// Raise indices of \tilde A_{ij} and store in R_ij

		M_ Rxx[_t]=    M_ gupxx[_t]* M_ gupxx[_t]* M_ Axx[_t]+ M_ gupxy[_t]* M_ gupxy[_t]* M_ Ayy[_t]+ M_ gupxz[_t]* M_ gupxz[_t]* M_ Azz[_t]+ 
		    2*(M_ gupxx[_t]* M_ gupxy[_t]* M_ Axy[_t]+ M_ gupxx[_t]* M_ gupxz[_t]* M_ Axz[_t]+ M_ gupxy[_t]* M_ gupxz[_t]* M_ Ayz[_t]);

		M_ Ryy[_t]=    M_ gupxy[_t]* M_ gupxy[_t]* M_ Axx[_t]+ M_ gupyy[_t]* M_ gupyy[_t]* M_ Ayy[_t]+ M_ gupyz[_t]* M_ gupyz[_t]* M_ Azz[_t]+ 
		    2*(M_ gupxy[_t]* M_ gupyy[_t]* M_ Axy[_t]+ M_ gupxy[_t]* M_ gupyz[_t]* M_ Axz[_t]+ M_ gupyy[_t]* M_ gupyz[_t]* M_ Ayz[_t]);

		M_ Rzz[_t]=    M_ gupxz[_t]* M_ gupxz[_t]* M_ Axx[_t]+ M_ gupyz[_t]* M_ gupyz[_t]* M_ Ayy[_t]+ M_ gupzz[_t]* M_ gupzz[_t]* M_ Azz[_t]+ 
		    2*(M_ gupxz[_t]* M_ gupyz[_t]* M_ Axy[_t]+ M_ gupxz[_t]* M_ gupzz[_t]* M_ Axz[_t]+ M_ gupyz[_t]* M_ gupzz[_t]* M_ Ayz[_t]);

		M_ Rxy[_t]=    M_ gupxx[_t]* M_ gupxy[_t]* M_ Axx[_t]+ M_ gupxy[_t]* M_ gupyy[_t]* M_ Ayy[_t]+ M_ gupxz[_t]* M_ gupyz[_t]* M_ Azz[_t]+ 
		        (M_ gupxx[_t]* M_ gupyy[_t]      + M_ gupxy[_t]* M_ gupxy[_t])* M_ Axy[_t]                      + 
		        (M_ gupxx[_t]* M_ gupyz[_t]      + M_ gupxz[_t]* M_ gupxy[_t])* M_ Axz[_t]                      + 
		        (M_ gupxy[_t]* M_ gupyz[_t]      + M_ gupxz[_t]* M_ gupyy[_t])* M_ Ayz[_t];

		M_ Rxz[_t]=    M_ gupxx[_t]* M_ gupxz[_t]* M_ Axx[_t]+ M_ gupxy[_t]* M_ gupyz[_t]* M_ Ayy[_t]+ M_ gupxz[_t]* M_ gupzz[_t]* M_ Azz[_t]+ 
		        (M_ gupxx[_t]* M_ gupyz[_t]      + M_ gupxy[_t]* M_ gupxz[_t])* M_ Axy[_t]                      + 
		        (M_ gupxx[_t]* M_ gupzz[_t]      + M_ gupxz[_t]* M_ gupxz[_t])* M_ Axz[_t]                      + 
		        (M_ gupxy[_t]* M_ gupzz[_t]      + M_ gupxz[_t]* M_ gupyz[_t])* M_ Ayz[_t];

		M_ Ryz[_t]=    M_ gupxy[_t]* M_ gupxz[_t]* M_ Axx[_t]+ M_ gupyy[_t]* M_ gupyz[_t]* M_ Ayy[_t]+ M_ gupyz[_t]* M_ gupzz[_t]* M_ Azz[_t]+ 
		        (M_ gupxy[_t]* M_ gupyz[_t]      + M_ gupyy[_t]* M_ gupxz[_t])* M_ Axy[_t]                      + 
		        (M_ gupxy[_t]* M_ gupzz[_t]      + M_ gupyz[_t]* M_ gupxz[_t])* M_ Axz[_t]                      + 
		        (M_ gupyy[_t]* M_ gupzz[_t]      + M_ gupyz[_t]* M_ gupyz[_t])* M_ Ayz[_t];

		// Right hand side for M_ Gam^i without shift terms...

		M_ Gamx_rhs[_t] = - 2 * (   M_ Lapx[_t] * M_ Rxx[_t]+   M_ Lapy[_t] * M_ Rxy[_t]+   M_ Lapz[_t] * M_ Rxz[_t]) + 
		    2 * M_ alpn1[_t] * (                                                
		    -F3o2/M_ chin1[_t] * (   M_ chix[_t] * M_ Rxx[_t]+   M_ chiy[_t] * M_ Rxy[_t]+   M_ chiz[_t] * M_ Rxz[_t]) - 
		            M_ gupxx[_t]* (   F2o3 * M_ Kx[_t]  +  8 * PI * M_ Sx[_t]            ) - 
		            M_ gupxy[_t]* (   F2o3 * M_ Ky[_t]  +  8 * PI * M_ Sy[_t]           ) - 
		            M_ gupxz[_t]* (   F2o3 * M_ Kz[_t]  +  8 * PI * M_ Sz[_t]            ) + 
		                    M_ Gamxxx[_t]* M_ Rxx[_t]+ M_ Gamxyy[_t]* M_ Ryy[_t]+ M_ Gamxzz[_t]* M_ Rzz[_t]  + 
		            2 * ( M_ Gamxxy[_t]* M_ Rxy[_t]+ M_ Gamxxz[_t]* M_ Rxz[_t]+ M_ Gamxyz[_t]* M_ Ryz[_t]) );

		M_ Gamy_rhs[_t] = - 2 * (   M_ Lapx[_t] * M_ Rxy[_t]+   M_ Lapy[_t] * M_ Ryy[_t]+   M_ Lapz[_t] * M_ Ryz[_t]) + 
		    2 * M_ alpn1[_t] * (                                                
		    -F3o2/M_ chin1[_t] * (   M_ chix[_t] * M_ Rxy[_t]+  M_ chiy[_t] * M_ Ryy[_t]+    M_ chiz[_t] * M_ Ryz[_t]) - 
		            M_ gupxy[_t]* (   F2o3 * M_ Kx[_t]  +  8 * PI * M_ Sx[_t]            ) - 
		            M_ gupyy[_t]* (   F2o3 * M_ Ky[_t]  +  8 * PI * M_ Sy[_t]            ) - 
		            M_ gupyz[_t]* (   F2o3 * M_ Kz [_t] +  8 * PI * M_ Sz[_t]            ) + 
		                    M_ Gamyxx[_t]* M_ Rxx[_t]+ M_ Gamyyy[_t]* M_ Ryy[_t]+ M_ Gamyzz[_t]* M_ Rzz[_t]  + 
		            2 * ( M_ Gamyxy[_t]* M_ Rxy[_t]+ M_ Gamyxz[_t]* M_ Rxz[_t]+ M_ Gamyyz[_t]* M_ Ryz[_t]) );

		M_ Gamz_rhs[_t] = - 2 * (   M_ Lapx[_t] * M_ Rxz[_t]+   M_ Lapy[_t] * M_ Ryz[_t]+   M_ Lapz[_t] * M_ Rzz[_t]) + 
		    2 * M_ alpn1[_t] * (                                                
		    -F3o2/M_ chin1[_t] * (   M_ chix[_t] * M_ Rxz[_t]+  M_ chiy[_t] * M_ Ryz[_t]+    M_ chiz[_t] * M_ Rzz[_t]) - 
		            M_ gupxz[_t]* (   F2o3 * M_ Kx[_t]  +  8 * PI * M_ Sx[_t]            ) - 
		            M_ gupyz[_t]* (   F2o3 * M_ Ky[_t]  +  8 * PI * M_ Sy[_t]            ) - 
		            M_ gupzz[_t]* (   F2o3 * M_ Kz[_t]  +  8 * PI * M_ Sz[_t]           ) + 
		                    M_ Gamzxx[_t]* M_ Rxx[_t]+ M_ Gamzyy[_t]* M_ Ryy[_t]+ M_ Gamzzz[_t]* M_ Rzz[_t]  + 
		            2 * ( M_ Gamzxy[_t]* M_ Rxy[_t]+ M_ Gamzxz[_t]* M_ Rxz[_t]+ M_ Gamzyz[_t]* M_ Ryz[_t]) );
			
		_t += STEP_SIZE;
	}
}

__global__ void compute_rhs_bssn_part3() 
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
		M_ fxx [_t]= M_ gxxx[_t]+ M_ gxyy[_t]+ M_ gxzz[_t];
		M_ fxy[_t]= M_ gxyx[_t]+ M_ gyyy[_t]+ M_ gyzz[_t];
		M_ fxz[_t]= M_ gxzx[_t]+ M_ gyzy[_t]+ M_ gzzz[_t];

		M_ Gamxa[_t]=       M_ gupxx [_t]* M_ Gamxxx [_t]+ M_ gupyy[_t]* M_ Gamxyy[_t]+ M_ gupzz[_t]* M_ Gamxzz[_t]+ 
		    2*( M_ gupxy[_t]* M_ Gamxxy[_t]+ M_ gupxz[_t]* M_ Gamxxz[_t]+ M_ gupyz[_t]* M_ Gamxyz[_t]);
		M_ Gamya[_t]=       M_ gupxx [_t]* M_ Gamyxx [_t]+ M_ gupyy[_t]* M_ Gamyyy[_t]+ M_ gupzz[_t]* M_ Gamyzz[_t]+ 
		    2*( M_ gupxy[_t]* M_ Gamyxy[_t]+ M_ gupxz[_t]* M_ Gamyxz[_t]+ M_ gupyz[_t]* M_ Gamyyz[_t]);
		M_ Gamza[_t]=       M_ gupxx [_t]* M_ Gamzxx [_t]+ M_ gupyy[_t]* M_ Gamzyy[_t]+ M_ gupzz[_t]* M_ Gamzzz[_t]+ 
		    2*( M_ gupxy[_t]* M_ Gamzxy[_t]+ M_ gupxz[_t]* M_ Gamzxz[_t]+ M_ gupyz[_t]* M_ Gamzyz[_t]);



		M_ Gamx_rhs[_t] =  M_ Gamx_rhs[_t] +  F2o3 *  M_ Gamxa[_t]* M_ div_beta[_t]       - 
		                M_ Gamxa[_t]* M_ betaxx [_t]- M_ Gamya[_t]* M_ betaxy[_t]- M_ Gamza[_t]* M_ betaxz[_t] + 
		        F1o3 * (M_ gupxx [_t]* M_ fxx [_t]   + M_ gupxy[_t]* M_ fxy[_t]   + M_ gupxz[_t]* M_ fxz[_t]   ) + 
		                M_ gupxx [_t]* M_ gxxx [_t]  + M_ gupyy[_t]* M_ gyyx [_t]  + M_ gupzz[_t]* M_ gzzx [_t]   + 
		        2 * (M_ gupxy[_t]* M_ gxyx [_t]  + M_ gupxz[_t]* M_ gxzx [_t]  + M_ gupyz[_t]* M_ gyzx [_t] );

		M_ Gamy_rhs[_t] =               M_ Gamy_rhs[_t] +  F2o3 *  M_ Gamya[_t]* M_ div_beta[_t]       - 
		                M_ Gamxa[_t]* M_ betayx [_t]- M_ Gamya[_t]* M_ betayy[_t]- M_ Gamza[_t]* M_ betayz[_t] + 
		        F1o3 * (M_ gupxy[_t]* M_ fxx [_t]   + M_ gupyy[_t]* M_ fxy[_t]   + M_ gupyz[_t]* M_ fxz[_t]   ) + 
		                M_ gupxx [_t]* M_ gxxy[_t]  + M_ gupyy[_t]* M_ gyyy[_t]  + M_ gupzz[_t]* M_ gzzy[_t]   + 
		        2 * (M_ gupxy[_t]* M_ gxyy[_t]  + M_ gupxz[_t]* M_ gxzy[_t]  + M_ gupyz[_t]* M_ gyzy[_t] );

		M_ Gamz_rhs[_t] =               M_ Gamz_rhs[_t] +  F2o3 *  M_ Gamza[_t]* M_ div_beta[_t]       - 
		                M_ Gamxa[_t]* M_ betazx [_t]- M_ Gamya[_t]* M_ betazy[_t]- M_ Gamza[_t]* M_ betazz[_t] + 
		        F1o3 * (M_ gupxz[_t]* M_ fxx [_t]   + M_ gupyz[_t]* M_ fxy[_t]   + M_ gupzz[_t]* M_ fxz[_t]   ) + 
		                M_ gupxx [_t]* M_ gxxz[_t]  + M_ gupyy[_t]* M_ gyyz[_t]  + M_ gupzz[_t]* M_ gzzz[_t]   + 
		        2 * (M_ gupxy[_t]* M_ gxyz[_t]  + M_ gupxz[_t]* M_ gxzz[_t]  + M_ gupyz[_t]* M_ gyzz[_t] )  ;  //rhs M_ for M_ Gam^i

		//first kind of connection stored in M_ gij,k
		M_ gxxx [_t]= M_ gxx [_t]* M_ Gamxxx [_t]+ M_ gxy[_t]* M_ Gamyxx [_t]+ M_ gxz[_t]* M_ Gamzxx[_t];
		M_ gxyx [_t]= M_ gxx [_t]* M_ Gamxxy[_t]+ M_ gxy[_t]* M_ Gamyxy[_t]+ M_ gxz[_t]* M_ Gamzxy[_t];
		M_ gxzx [_t]= M_ gxx [_t]* M_ Gamxxz[_t]+ M_ gxy[_t]* M_ Gamyxz[_t]+ M_ gxz[_t]* M_ Gamzxz[_t];
		M_ gyyx [_t]= M_ gxx [_t]* M_ Gamxyy[_t]+ M_ gxy[_t]* M_ Gamyyy[_t]+ M_ gxz[_t]* M_ Gamzyy[_t];
		M_ gyzx [_t]= M_ gxx [_t]* M_ Gamxyz[_t]+ M_ gxy[_t]* M_ Gamyyz[_t]+ M_ gxz[_t]* M_ Gamzyz[_t];
		M_ gzzx [_t]= M_ gxx [_t]* M_ Gamxzz[_t]+ M_ gxy[_t]* M_ Gamyzz[_t]+ M_ gxz[_t]* M_ Gamzzz[_t];
		M_ gxxy[_t]= M_ gxy[_t]* M_ Gamxxx [_t]+ M_ gyy[_t]* M_ Gamyxx [_t]+ M_ gyz[_t]* M_ Gamzxx[_t];
		M_ gxyy[_t]= M_ gxy[_t]* M_ Gamxxy[_t]+ M_ gyy[_t]* M_ Gamyxy[_t]+ M_ gyz[_t]* M_ Gamzxy[_t];
		M_ gxzy[_t]= M_ gxy[_t]* M_ Gamxxz[_t]+ M_ gyy[_t]* M_ Gamyxz[_t]+ M_ gyz[_t]* M_ Gamzxz[_t];
		M_ gyyy[_t]= M_ gxy[_t]* M_ Gamxyy[_t]+ M_ gyy[_t]* M_ Gamyyy[_t]+ M_ gyz[_t]* M_ Gamzyy[_t];
		M_ gyzy[_t]= M_ gxy[_t]* M_ Gamxyz[_t]+ M_ gyy[_t]* M_ Gamyyz[_t]+ M_ gyz[_t]* M_ Gamzyz[_t];
		M_ gzzy[_t]= M_ gxy[_t]* M_ Gamxzz[_t]+ M_ gyy[_t]* M_ Gamyzz[_t]+ M_ gyz[_t]* M_ Gamzzz[_t];
		M_ gxxz[_t]= M_ gxz[_t]* M_ Gamxxx [_t]+ M_ gyz[_t]* M_ Gamyxx [_t]+ M_ gzz[_t]* M_ Gamzxx[_t];
		M_ gxyz[_t]= M_ gxz[_t]* M_ Gamxxy[_t]+ M_ gyz[_t]* M_ Gamyxy[_t]+ M_ gzz[_t]* M_ Gamzxy[_t];
		M_ gxzz[_t]= M_ gxz[_t]* M_ Gamxxz[_t]+ M_ gyz[_t]* M_ Gamyxz[_t]+ M_ gzz[_t]* M_ Gamzxz[_t];
		M_ gyyz[_t]= M_ gxz[_t]* M_ Gamxyy[_t]+ M_ gyz[_t]* M_ Gamyyy[_t]+ M_ gzz[_t]* M_ Gamzyy[_t];
		M_ gyzz[_t]= M_ gxz[_t]* M_ Gamxyz[_t]+ M_ gyz[_t]* M_ Gamyyz[_t]+ M_ gzz[_t]* M_ Gamzyz[_t];
		M_ gzzz[_t]= M_ gxz[_t]* M_ Gamxzz[_t]+ M_ gyz[_t]* M_ Gamyzz[_t]+ M_ gzz[_t]* M_ Gamzzz[_t];
		
	  	
	  	_t += STEP_SIZE;
	}
}

__global__ void compute_rhs_bssn_part4() 
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
		M_ Rxx [_t]=     - HALF *M_ Rxx [_t]                                  + 
		        M_ gxx [_t]* M_ Gamxx[_t] +M_ gxy[_t]* M_ Gamyx [_t]  +   M_ gxz[_t]* M_ Gamzx [_t]+ 
		        M_ Gamxa[_t]*M_ gxxx [_t]+  M_ Gamya[_t]*M_ gxyx [_t]+  M_ Gamza[_t]*M_ gxzx [_t] + 
		M_ gupxx [_t]*(                                                  
		    2*(M_ Gamxxx [_t]*M_ gxxx [_t]+ M_ Gamyxx [_t]*M_ gxyx [_t]+ M_ Gamzxx [_t]*M_ gxzx[_t]) + 
		        M_ Gamxxx [_t]*M_ gxxx [_t]+ M_ Gamyxx [_t]*M_ gxxy[_t]+ M_ Gamzxx [_t]*M_ gxxz[_t])+ 
		M_ gupxy[_t]*(                                                  
		    2*(M_ Gamxxx [_t]*M_ gxyx [_t]+ M_ Gamyxx [_t]*M_ gyyx [_t]+ M_ Gamzxx [_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gxxx [_t]+ M_ Gamyxy[_t]*M_ gxyx [_t]+ M_ Gamzxy[_t]*M_ gxzx[_t]) + 
		        M_ Gamxxy[_t]*M_ gxxx [_t]+ M_ Gamyxy[_t]*M_ gxxy[_t]+ M_ Gamzxy[_t]*M_ gxxz[_t] + 
		        M_ Gamxxx [_t]*M_ gxyx [_t]+ M_ Gamyxx [_t]*M_ gxyy[_t]+ M_ Gamzxx [_t]*M_ gxyz[_t])+ 
		M_ gupxz[_t]*(                                                  
		    2*(M_ Gamxxx [_t]*M_ gxzx [_t]+ M_ Gamyxx [_t]*M_ gyzx [_t]+ M_ Gamzxx [_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gxxx [_t]+ M_ Gamyxz[_t]*M_ gxyx [_t]+ M_ Gamzxz[_t]*M_ gxzx[_t]) + 
		        M_ Gamxxz[_t]*M_ gxxx [_t]+ M_ Gamyxz[_t]*M_ gxxy[_t]+ M_ Gamzxz[_t]*M_ gxxz[_t] + 
		        M_ Gamxxx [_t]*M_ gxzx [_t]+ M_ Gamyxx [_t]*M_ gxzy[_t]+ M_ Gamzxx [_t]*M_ gxzz[_t])+ 
		M_ gupyy[_t]*(                                                  
		    2*(M_ Gamxxy[_t]*M_ gxyx [_t]+ M_ Gamyxy[_t]*M_ gyyx [_t]+ M_ Gamzxy[_t]*M_ gyzx[_t]) + 
		        M_ Gamxxy[_t]*M_ gxyx [_t]+ M_ Gamyxy[_t]*M_ gxyy[_t]+ M_ Gamzxy[_t]*M_ gxyz[_t])+ 
		M_ gupyz[_t]*(                                                  
		    2*(M_ Gamxxy[_t]*M_ gxzx [_t]+ M_ Gamyxy[_t]*M_ gyzx [_t]+ M_ Gamzxy[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gxyx [_t]+ M_ Gamyxz[_t]*M_ gyyx [_t]+ M_ Gamzxz[_t]*M_ gyzx[_t]) + 
		        M_ Gamxxz[_t]*M_ gxyx [_t]+ M_ Gamyxz[_t]*M_ gxyy[_t]+ M_ Gamzxz[_t]*M_ gxyz[_t] + 
		        M_ Gamxxy[_t]*M_ gxzx [_t]+ M_ Gamyxy[_t]*M_ gxzy[_t]+ M_ Gamzxy[_t]*M_ gxzz[_t])+ 
		M_ gupzz[_t]*(                                                  
		    2*(M_ Gamxxz[_t]*M_ gxzx [_t]+ M_ Gamyxz[_t]*M_ gyzx [_t]+ M_ Gamzxz[_t]*M_ gzzx[_t]) + 
		        M_ Gamxxz[_t]*M_ gxzx [_t]+ M_ Gamyxz[_t]*M_ gxzy[_t]+ M_ Gamzxz[_t]*M_ gxzz[_t]);

		M_ Ryy[_t]=     - HALF *M_ Ryy[_t]                                  + 
		            M_ gxy[_t]* M_ Gamxy[_t]+ M_ gyy[_t]* M_ Gamyy[_t] + M_ gyz[_t]* M_ Gamzy[_t]  + 
		            M_ Gamxa[_t]*M_ gxyy[_t]+  M_ Gamya[_t]*M_ gyyy[_t]+  M_ Gamza[_t]*M_ gyzy[_t] + 
		M_ gupxx [_t]*(                                                  
		    2*(M_ Gamxxy[_t]*M_ gxxy[_t]+ M_ Gamyxy[_t]*M_ gxyy[_t]+ M_ Gamzxy[_t]*M_ gxzy[_t]) + 
		        M_ Gamxxy[_t]*M_ gxyx [_t]+ M_ Gamyxy[_t]*M_ gxyy[_t]+ M_ Gamzxy[_t]*M_ gxyz[_t])+ 
		M_ gupxy[_t]*(                                                  
		    2*(M_ Gamxxy[_t]*M_ gxyy[_t]+ M_ Gamyxy[_t]*M_ gyyy[_t]+ M_ Gamzxy[_t]*M_ gyzy[_t] + 
		        M_ Gamxyy[_t]*M_ gxxy[_t]+ M_ Gamyyy[_t]*M_ gxyy[_t]+ M_ Gamzyy[_t]*M_ gxzy[_t]) + 
		        M_ Gamxyy[_t]*M_ gxyx [_t]+ M_ Gamyyy[_t]*M_ gxyy[_t]+ M_ Gamzyy[_t]*M_ gxyz[_t] + 
		        M_ Gamxxy[_t]*M_ gyyx [_t]+ M_ Gamyxy[_t]*M_ gyyy[_t]+ M_ Gamzxy[_t]*M_ gyyz[_t])+ 
		M_ gupxz[_t]*(                                                  
		    2*(M_ Gamxxy[_t]*M_ gxzy[_t]+ M_ Gamyxy[_t]*M_ gyzy[_t]+ M_ Gamzxy[_t]*M_ gzzy[_t] + 
		        M_ Gamxyz[_t]*M_ gxxy[_t]+ M_ Gamyyz[_t]*M_ gxyy[_t]+ M_ Gamzyz[_t]*M_ gxzy[_t]) + 
		        M_ Gamxyz[_t]*M_ gxyx [_t]+ M_ Gamyyz[_t]*M_ gxyy[_t]+ M_ Gamzyz[_t]*M_ gxyz[_t] + 
		        M_ Gamxxy[_t]*M_ gyzx [_t]+ M_ Gamyxy[_t]*M_ gyzy[_t]+ M_ Gamzxy[_t]*M_ gyzz[_t])+ 
		M_ gupyy[_t]*(                                                  
		    2*(M_ Gamxyy[_t]*M_ gxyy[_t]+ M_ Gamyyy[_t]*M_ gyyy[_t]+ M_ Gamzyy[_t]*M_ gyzy[_t]) + 
		        M_ Gamxyy[_t]*M_ gyyx [_t]+ M_ Gamyyy[_t]*M_ gyyy[_t]+ M_ Gamzyy[_t]*M_ gyyz[_t])+ 
		M_ gupyz[_t]*(                                                  
		    2*(M_ Gamxyy[_t]*M_ gxzy[_t]+ M_ Gamyyy[_t]*M_ gyzy[_t]+ M_ Gamzyy[_t]*M_ gzzy[_t] + 
		        M_ Gamxyz[_t]*M_ gxyy[_t]+ M_ Gamyyz[_t]*M_ gyyy[_t]+ M_ Gamzyz[_t]*M_ gyzy[_t]) + 
		        M_ Gamxyz[_t]*M_ gyyx [_t]+ M_ Gamyyz[_t]*M_ gyyy[_t]+ M_ Gamzyz[_t]*M_ gyyz[_t] + 
		        M_ Gamxyy[_t]*M_ gyzx [_t]+ M_ Gamyyy[_t]*M_ gyzy[_t]+ M_ Gamzyy[_t]*M_ gyzz[_t])+ 
		M_ gupzz[_t]*(                                                  
		    2*(M_ Gamxyz[_t]*M_ gxzy[_t]+ M_ Gamyyz[_t]*M_ gyzy[_t]+ M_ Gamzyz[_t]*M_ gzzy[_t]) + 
		        M_ Gamxyz[_t]*M_ gyzx [_t]+ M_ Gamyyz[_t]*M_ gyzy[_t]+ M_ Gamzyz[_t]*M_ gyzz[_t]);

		M_ Rzz[_t]=     - HALF *M_ Rzz[_t]                                  + 
		            M_ gxz[_t]* M_ Gamxz[_t] +M_ gyz[_t]* M_ Gamyz[_t] +   M_ gzz[_t]* M_ Gamzz[_t] + 
		            M_ Gamxa[_t]*M_ gxzz[_t]+  M_ Gamya[_t]*M_ gyzz[_t]+  M_ Gamza[_t]*M_ gzzz[_t] + 
		M_ gupxx [_t]*(                                                  
		    2*(M_ Gamxxz[_t]*M_ gxxz[_t]+ M_ Gamyxz[_t]*M_ gxyz[_t]+ M_ Gamzxz[_t]*M_ gxzz[_t]) + 
		        M_ Gamxxz[_t]*M_ gxzx [_t]+ M_ Gamyxz[_t]*M_ gxzy[_t]+ M_ Gamzxz[_t]*M_ gxzz[_t])+ 
		M_ gupxy[_t]*(                                                  
		    2*(M_ Gamxxz[_t]*M_ gxyz[_t]+ M_ Gamyxz[_t]*M_ gyyz[_t]+ M_ Gamzxz[_t]*M_ gyzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxxz[_t]+ M_ Gamyyz[_t]*M_ gxyz[_t]+ M_ Gamzyz[_t]*M_ gxzz[_t]) + 
		        M_ Gamxyz[_t]*M_ gxzx [_t]+ M_ Gamyyz[_t]*M_ gxzy[_t]+ M_ Gamzyz[_t]*M_ gxzz[_t] + 
		        M_ Gamxxz[_t]*M_ gyzx [_t]+ M_ Gamyxz[_t]*M_ gyzy[_t]+ M_ Gamzxz[_t]*M_ gyzz[_t])+ 
		M_ gupxz[_t]*(                                                  
		    2*(M_ Gamxxz[_t]*M_ gxzz[_t]+ M_ Gamyxz[_t]*M_ gyzz[_t]+ M_ Gamzxz[_t]*M_ gzzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxxz[_t]+ M_ Gamyzz[_t]*M_ gxyz[_t]+ M_ Gamzzz[_t]*M_ gxzz[_t]) + 
		        M_ Gamxzz[_t]*M_ gxzx [_t]+ M_ Gamyzz[_t]*M_ gxzy[_t]+ M_ Gamzzz[_t]*M_ gxzz[_t] + 
		        M_ Gamxxz[_t]*M_ gzzx [_t]+ M_ Gamyxz[_t]*M_ gzzy[_t]+ M_ Gamzxz[_t]*M_ gzzz[_t])+ 
		M_ gupyy[_t]*(                                                  
		    2*(M_ Gamxyz[_t]*M_ gxyz[_t]+ M_ Gamyyz[_t]*M_ gyyz[_t]+ M_ Gamzyz[_t]*M_ gyzz[_t]) + 
		        M_ Gamxyz[_t]*M_ gyzx [_t]+ M_ Gamyyz[_t]*M_ gyzy[_t]+ M_ Gamzyz[_t]*M_ gyzz[_t])+ 
		M_ gupyz[_t]*(                                                  
		    2*(M_ Gamxyz[_t]*M_ gxzz[_t]+ M_ Gamyyz[_t]*M_ gyzz[_t]+ M_ Gamzyz[_t]*M_ gzzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxyz[_t]+ M_ Gamyzz[_t]*M_ gyyz[_t]+ M_ Gamzzz[_t]*M_ gyzz[_t]) + 
		        M_ Gamxzz[_t]*M_ gyzx [_t]+ M_ Gamyzz[_t]*M_ gyzy[_t]+ M_ Gamzzz[_t]*M_ gyzz[_t] + 
		        M_ Gamxyz[_t]*M_ gzzx [_t]+ M_ Gamyyz[_t]*M_ gzzy[_t]+ M_ Gamzyz[_t]*M_ gzzz[_t])+ 
		M_ gupzz[_t]*(                                                  
		    2*(M_ Gamxzz[_t]*M_ gxzz[_t]+ M_ Gamyzz[_t]*M_ gyzz[_t]+ M_ Gamzzz[_t]*M_ gzzz[_t]) + 
		        M_ Gamxzz[_t]*M_ gzzx [_t]+ M_ Gamyzz[_t]*M_ gzzy[_t]+ M_ Gamzzz[_t]*M_ gzzz[_t]);

		M_ Rxy[_t]= HALF*(     -M_ Rxy[_t]                                  + 
		            M_ gxx [_t]* M_ Gamxy[_t]+   M_ gxy[_t]* M_ Gamyy[_t]+M_ gxz[_t]* M_ Gamzy[_t] + 
		            M_ gxy[_t]* M_ Gamxx [_t]+   M_ gyy[_t]* M_ Gamyx [_t]+M_ gyz[_t]* M_ Gamzx [_t] + 
		            M_ Gamxa[_t]*M_ gxyx [_t]+  M_ Gamya[_t]*M_ gyyx [_t]+  M_ Gamza[_t]*M_ gyzx [_t] + 
		            M_ Gamxa[_t]*M_ gxxy[_t]+  M_ Gamya[_t]*M_ gxyy[_t]+  M_ Gamza[_t]*M_ gxzy[_t])+ 
		M_ gupxx [_t]*(                                                  
		        M_ Gamxxx [_t]*M_ gxxy[_t]+ M_ Gamyxx [_t]*M_ gxyy[_t]+ M_ Gamzxx [_t]*M_ gxzy[_t] + 
		        M_ Gamxxy[_t]*M_ gxxx [_t]+ M_ Gamyxy[_t]*M_ gxyx [_t]+ M_ Gamzxy[_t]*M_ gxzx [_t] + 
		        M_ Gamxxx [_t]*M_ gxyx [_t]+ M_ Gamyxx [_t]*M_ gxyy[_t]+ M_ Gamzxx [_t]*M_ gxyz[_t])+ 
		M_ gupxy[_t]*(                                                  
		        M_ Gamxxx [_t]*M_ gxyy[_t]+ M_ Gamyxx [_t]*M_ gyyy[_t]+ M_ Gamzxx [_t]*M_ gyzy[_t] + 
		        M_ Gamxxy[_t]*M_ gxyx [_t]+ M_ Gamyxy[_t]*M_ gyyx [_t]+ M_ Gamzxy[_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gxyx [_t]+ M_ Gamyxy[_t]*M_ gxyy[_t]+ M_ Gamzxy[_t]*M_ gxyz[_t] + 
		        M_ Gamxxy[_t]*M_ gxxy[_t]+ M_ Gamyxy[_t]*M_ gxyy[_t]+ M_ Gamzxy[_t]*M_ gxzy[_t] + 
		        M_ Gamxyy[_t]*M_ gxxx [_t]+ M_ Gamyyy[_t]*M_ gxyx [_t]+ M_ Gamzyy[_t]*M_ gxzx [_t] + 
		        M_ Gamxxx [_t]*M_ gyyx [_t]+ M_ Gamyxx [_t]*M_ gyyy[_t]+ M_ Gamzxx [_t]*M_ gyyz[_t])+ 
		M_ gupxz[_t]*(                                                  
		        M_ Gamxxx [_t]*M_ gxzy[_t]+ M_ Gamyxx [_t]*M_ gyzy[_t]+ M_ Gamzxx [_t]*M_ gzzy[_t] + 
		        M_ Gamxxy[_t]*M_ gxzx [_t]+ M_ Gamyxy[_t]*M_ gyzx [_t]+ M_ Gamzxy[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gxyx [_t]+ M_ Gamyxz[_t]*M_ gxyy[_t]+ M_ Gamzxz[_t]*M_ gxyz[_t] + 
		        M_ Gamxxz[_t]*M_ gxxy[_t]+ M_ Gamyxz[_t]*M_ gxyy[_t]+ M_ Gamzxz[_t]*M_ gxzy[_t] + 
		        M_ Gamxyz[_t]*M_ gxxx [_t]+ M_ Gamyyz[_t]*M_ gxyx [_t]+ M_ Gamzyz[_t]*M_ gxzx [_t] + 
		        M_ Gamxxx [_t]*M_ gyzx [_t]+ M_ Gamyxx [_t]*M_ gyzy[_t]+ M_ Gamzxx [_t]*M_ gyzz[_t])+ 
		M_ gupyy[_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxyy[_t]+ M_ Gamyxy[_t]*M_ gyyy[_t]+ M_ Gamzxy[_t]*M_ gyzy[_t] + 
		        M_ Gamxyy[_t]*M_ gxyx [_t]+ M_ Gamyyy[_t]*M_ gyyx [_t]+ M_ Gamzyy[_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gyyx [_t]+ M_ Gamyxy[_t]*M_ gyyy[_t]+ M_ Gamzxy[_t]*M_ gyyz[_t])+ 
		M_ gupyz[_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxzy[_t]+ M_ Gamyxy[_t]*M_ gyzy[_t]+ M_ Gamzxy[_t]*M_ gzzy[_t] + 
		        M_ Gamxyy[_t]*M_ gxzx [_t]+ M_ Gamyyy[_t]*M_ gyzx [_t]+ M_ Gamzyy[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gyyx [_t]+ M_ Gamyxz[_t]*M_ gyyy[_t]+ M_ Gamzxz[_t]*M_ gyyz[_t] + 
		        M_ Gamxxz[_t]*M_ gxyy[_t]+ M_ Gamyxz[_t]*M_ gyyy[_t]+ M_ Gamzxz[_t]*M_ gyzy[_t] + 
		        M_ Gamxyz[_t]*M_ gxyx [_t]+ M_ Gamyyz[_t]*M_ gyyx [_t]+ M_ Gamzyz[_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gyzx [_t]+ M_ Gamyxy[_t]*M_ gyzy[_t]+ M_ Gamzxy[_t]*M_ gyzz[_t])+ 
		M_ gupzz[_t]*(                                                  
		        M_ Gamxxz[_t]*M_ gxzy[_t]+ M_ Gamyxz[_t]*M_ gyzy[_t]+ M_ Gamzxz[_t]*M_ gzzy[_t] + 
		        M_ Gamxyz[_t]*M_ gxzx [_t]+ M_ Gamyyz[_t]*M_ gyzx [_t]+ M_ Gamzyz[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gyzx [_t]+ M_ Gamyxz[_t]*M_ gyzy[_t]+ M_ Gamzxz[_t]*M_ gyzz[_t]);

		M_ Rxz[_t]= HALF*(     -M_ Rxz[_t]                                  + 
		            M_ gxx [_t]* M_ Gamxz[_t]+ M_ gxy[_t]* M_ Gamyz[_t]+M_ gxz[_t]* M_ Gamzz[_t]   + 
		            M_ gxz[_t]* M_ Gamxx [_t]+ M_ gyz[_t]* M_ Gamyx [_t]+M_ gzz[_t]* M_ Gamzx [_t]   + 
		            M_ Gamxa[_t]*M_ gxzx [_t]+  M_ Gamya[_t]*M_ gyzx [_t]+  M_ Gamza[_t]*M_ gzzx [_t] + 
		            M_ Gamxa[_t]*M_ gxxz[_t]+  M_ Gamya[_t]*M_ gxyz[_t]+  M_ Gamza[_t]*M_ gxzz[_t])+ 
		M_ gupxx [_t]*(                                                  
		        M_ Gamxxx [_t]*M_ gxxz[_t]+ M_ Gamyxx [_t]*M_ gxyz[_t]+ M_ Gamzxx [_t]*M_ gxzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxxx [_t]+ M_ Gamyxz[_t]*M_ gxyx [_t]+ M_ Gamzxz[_t]*M_ gxzx [_t] + 
		        M_ Gamxxx [_t]*M_ gxzx [_t]+ M_ Gamyxx [_t]*M_ gxzy[_t]+ M_ Gamzxx [_t]*M_ gxzz[_t])+ 
		M_ gupxy[_t]*(                                                  
		        M_ Gamxxx [_t]*M_ gxyz[_t]+ M_ Gamyxx [_t]*M_ gyyz[_t]+ M_ Gamzxx [_t]*M_ gyzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxyx [_t]+ M_ Gamyxz[_t]*M_ gyyx [_t]+ M_ Gamzxz[_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gxzx [_t]+ M_ Gamyxy[_t]*M_ gxzy[_t]+ M_ Gamzxy[_t]*M_ gxzz[_t] + 
		        M_ Gamxxy[_t]*M_ gxxz[_t]+ M_ Gamyxy[_t]*M_ gxyz[_t]+ M_ Gamzxy[_t]*M_ gxzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxxx [_t]+ M_ Gamyyz[_t]*M_ gxyx [_t]+ M_ Gamzyz[_t]*M_ gxzx [_t] + 
		        M_ Gamxxx [_t]*M_ gyzx [_t]+ M_ Gamyxx [_t]*M_ gyzy[_t]+ M_ Gamzxx [_t]*M_ gyzz[_t])+ 
		M_ gupxz[_t]*(                                                  
		        M_ Gamxxx [_t]*M_ gxzz[_t]+ M_ Gamyxx [_t]*M_ gyzz[_t]+ M_ Gamzxx [_t]*M_ gzzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxzx [_t]+ M_ Gamyxz[_t]*M_ gyzx [_t]+ M_ Gamzxz[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gxzx [_t]+ M_ Gamyxz[_t]*M_ gxzy[_t]+ M_ Gamzxz[_t]*M_ gxzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxxz[_t]+ M_ Gamyxz[_t]*M_ gxyz[_t]+ M_ Gamzxz[_t]*M_ gxzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxxx [_t]+ M_ Gamyzz[_t]*M_ gxyx [_t]+ M_ Gamzzz[_t]*M_ gxzx [_t] + 
		        M_ Gamxxx [_t]*M_ gzzx [_t]+ M_ Gamyxx [_t]*M_ gzzy[_t]+ M_ Gamzxx [_t]*M_ gzzz[_t])+ 
		M_ gupyy[_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxyz[_t]+ M_ Gamyxy[_t]*M_ gyyz[_t]+ M_ Gamzxy[_t]*M_ gyzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxyx [_t]+ M_ Gamyyz[_t]*M_ gyyx [_t]+ M_ Gamzyz[_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gyzx [_t]+ M_ Gamyxy[_t]*M_ gyzy[_t]+ M_ Gamzxy[_t]*M_ gyzz[_t])+ 
		M_ gupyz[_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxzz[_t]+ M_ Gamyxy[_t]*M_ gyzz[_t]+ M_ Gamzxy[_t]*M_ gzzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxzx [_t]+ M_ Gamyyz[_t]*M_ gyzx [_t]+ M_ Gamzyz[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gyzx [_t]+ M_ Gamyxz[_t]*M_ gyzy[_t]+ M_ Gamzxz[_t]*M_ gyzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxyz[_t]+ M_ Gamyxz[_t]*M_ gyyz[_t]+ M_ Gamzxz[_t]*M_ gyzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxyx [_t]+ M_ Gamyzz[_t]*M_ gyyx [_t]+ M_ Gamzzz[_t]*M_ gyzx [_t] + 
		        M_ Gamxxy[_t]*M_ gzzx [_t]+ M_ Gamyxy[_t]*M_ gzzy[_t]+ M_ Gamzxy[_t]*M_ gzzz[_t])+ 
		M_ gupzz[_t]*(                                                  
		        M_ Gamxxz[_t]*M_ gxzz[_t]+ M_ Gamyxz[_t]*M_ gyzz[_t]+ M_ Gamzxz[_t]*M_ gzzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxzx [_t]+ M_ Gamyzz[_t]*M_ gyzx [_t]+ M_ Gamzzz[_t]*M_ gzzx [_t] + 
		        M_ Gamxxz[_t]*M_ gzzx [_t]+ M_ Gamyxz[_t]*M_ gzzy[_t]+ M_ Gamzxz[_t]*M_ gzzz[_t]);

		M_ Ryz[_t]= HALF*(     -M_ Ryz[_t]                                  + 
		            M_ gxy[_t]* M_ Gamxz[_t]+M_ gyy[_t]* M_ Gamyz[_t]+M_ gyz[_t]* M_ Gamzz[_t]    + 
		            M_ gxz[_t]* M_ Gamxy[_t]+M_ gyz[_t]* M_ Gamyy[_t]+M_ gzz[_t]* M_ Gamzy[_t]    + 
		            M_ Gamxa[_t]*M_ gxzy[_t]+  M_ Gamya[_t]*M_ gyzy[_t]+  M_ Gamza[_t]*M_ gzzy[_t] + 
		            M_ Gamxa[_t]*M_ gxyz[_t]+  M_ Gamya[_t]*M_ gyyz[_t]+  M_ Gamza[_t]*M_ gyzz[_t])+ 
		M_ gupxx [_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxxz[_t]+ M_ Gamyxy[_t]*M_ gxyz[_t]+ M_ Gamzxy[_t]*M_ gxzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxxy[_t]+ M_ Gamyxz[_t]*M_ gxyy[_t]+ M_ Gamzxz[_t]*M_ gxzy[_t] + 
		        M_ Gamxxy[_t]*M_ gxzx [_t]+ M_ Gamyxy[_t]*M_ gxzy[_t]+ M_ Gamzxy[_t]*M_ gxzz[_t])+ 
		M_ gupxy[_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxyz[_t]+ M_ Gamyxy[_t]*M_ gyyz[_t]+ M_ Gamzxy[_t]*M_ gyzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxyy[_t]+ M_ Gamyxz[_t]*M_ gyyy[_t]+ M_ Gamzxz[_t]*M_ gyzy[_t] + 
		        M_ Gamxyy[_t]*M_ gxzx [_t]+ M_ Gamyyy[_t]*M_ gxzy[_t]+ M_ Gamzyy[_t]*M_ gxzz[_t] + 
		        M_ Gamxyy[_t]*M_ gxxz[_t]+ M_ Gamyyy[_t]*M_ gxyz[_t]+ M_ Gamzyy[_t]*M_ gxzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxxy[_t]+ M_ Gamyyz[_t]*M_ gxyy[_t]+ M_ Gamzyz[_t]*M_ gxzy[_t] + 
		        M_ Gamxxy[_t]*M_ gyzx [_t]+ M_ Gamyxy[_t]*M_ gyzy[_t]+ M_ Gamzxy[_t]*M_ gyzz[_t])+ 
		M_ gupxz[_t]*(                                                  
		        M_ Gamxxy[_t]*M_ gxzz[_t]+ M_ Gamyxy[_t]*M_ gyzz[_t]+ M_ Gamzxy[_t]*M_ gzzz[_t] + 
		        M_ Gamxxz[_t]*M_ gxzy[_t]+ M_ Gamyxz[_t]*M_ gyzy[_t]+ M_ Gamzxz[_t]*M_ gzzy[_t] + 
		        M_ Gamxyz[_t]*M_ gxzx [_t]+ M_ Gamyyz[_t]*M_ gxzy[_t]+ M_ Gamzyz[_t]*M_ gxzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxxz[_t]+ M_ Gamyyz[_t]*M_ gxyz[_t]+ M_ Gamzyz[_t]*M_ gxzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxxy[_t]+ M_ Gamyzz[_t]*M_ gxyy[_t]+ M_ Gamzzz[_t]*M_ gxzy[_t] + 
		        M_ Gamxxy[_t]*M_ gzzx [_t]+ M_ Gamyxy[_t]*M_ gzzy[_t]+ M_ Gamzxy[_t]*M_ gzzz[_t])+ 
		M_ gupyy[_t]*(                                                  
		        M_ Gamxyy[_t]*M_ gxyz[_t]+ M_ Gamyyy[_t]*M_ gyyz[_t]+ M_ Gamzyy[_t]*M_ gyzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxyy[_t]+ M_ Gamyyz[_t]*M_ gyyy[_t]+ M_ Gamzyz[_t]*M_ gyzy[_t] + 
		        M_ Gamxyy[_t]*M_ gyzx [_t]+ M_ Gamyyy[_t]*M_ gyzy[_t]+ M_ Gamzyy[_t]*M_ gyzz[_t])+ 
		M_ gupyz[_t]*(                                                  
		        M_ Gamxyy[_t]*M_ gxzz[_t]+ M_ Gamyyy[_t]*M_ gyzz[_t]+ M_ Gamzyy[_t]*M_ gzzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxzy[_t]+ M_ Gamyyz[_t]*M_ gyzy[_t]+ M_ Gamzyz[_t]*M_ gzzy[_t] + 
		        M_ Gamxyz[_t]*M_ gyzx [_t]+ M_ Gamyyz[_t]*M_ gyzy[_t]+ M_ Gamzyz[_t]*M_ gyzz[_t] + 
		        M_ Gamxyz[_t]*M_ gxyz[_t]+ M_ Gamyyz[_t]*M_ gyyz[_t]+ M_ Gamzyz[_t]*M_ gyzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxyy[_t]+ M_ Gamyzz[_t]*M_ gyyy[_t]+ M_ Gamzzz[_t]*M_ gyzy[_t] + 
		        M_ Gamxyy[_t]*M_ gzzx [_t]+ M_ Gamyyy[_t]*M_ gzzy[_t]+ M_ Gamzyy[_t]*M_ gzzz[_t])+ 
		M_ gupzz[_t]*(                                                  
		        M_ Gamxyz[_t]*M_ gxzz[_t]+ M_ Gamyyz[_t]*M_ gyzz[_t]+ M_ Gamzyz[_t]*M_ gzzz[_t] + 
		        M_ Gamxzz[_t]*M_ gxzy[_t]+ M_ Gamyzz[_t]*M_ gyzy[_t]+ M_ Gamzzz[_t]*M_ gzzy[_t] + 
		        M_ Gamxyz[_t]*M_ gzzx [_t]+ M_ Gamyyz[_t]*M_ gzzy[_t]+ M_ Gamzyz[_t]*M_ gzzz[_t]);

		_t += STEP_SIZE;
	}
}
__global__ void compute_rhs_bssn_part5() 
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
		M_ fxx [_t]=M_ fxx [_t]- M_ Gamxxx [_t]* M_ chix [_t]- M_ Gamyxx [_t]* M_ chiy[_t]- M_ Gamzxx [_t]* M_ chiz[_t];
		M_ fxy[_t]=M_ fxy[_t]- M_ Gamxxy[_t]* M_ chix [_t]- M_ Gamyxy[_t]* M_ chiy[_t]- M_ Gamzxy[_t]* M_ chiz[_t];
		M_ fxz[_t]=M_ fxz[_t]- M_ Gamxxz[_t]* M_ chix [_t]- M_ Gamyxz[_t]* M_ chiy[_t]- M_ Gamzxz[_t]* M_ chiz[_t];
		M_ fyy[_t]=M_ fyy[_t]- M_ Gamxyy[_t]* M_ chix [_t]- M_ Gamyyy[_t]* M_ chiy[_t]- M_ Gamzyy[_t]* M_ chiz[_t];
		M_ fyz[_t]=M_ fyz[_t]- M_ Gamxyz[_t]* M_ chix [_t]- M_ Gamyyz[_t]* M_ chiy[_t]- M_ Gamzyz[_t]* M_ chiz[_t];
		M_ fzz[_t]=M_ fzz[_t]- M_ Gamxzz[_t]* M_ chix [_t]- M_ Gamyzz[_t]* M_ chiy[_t]- M_ Gamzzz[_t]* M_ chiz[_t];
		// M_ Store D^l D_l M_ chi - 3/(2*M_ chi) D^l M_ chi D_l M_ chi inM_ f[_t]

		M_ f[_t] =       M_ gupxx [_t]* (M_ fxx [_t]- F3o2/M_ chin1[_t] * M_ chix [_t]* M_ chix [_t]) + 
		        M_ gupyy[_t]* (M_ fyy[_t]- F3o2/M_ chin1[_t] * M_ chiy[_t]* M_ chiy[_t]) + 
		        M_ gupzz[_t]* (M_ fzz[_t]- F3o2/M_ chin1[_t] * M_ chiz[_t]* M_ chiz[_t]) + 
		    2 *M_ gupxy[_t]* (M_ fxy[_t]- F3o2/M_ chin1[_t] * M_ chix [_t]* M_ chiy[_t]) + 
		    2 *M_ gupxz[_t]* (M_ fxz[_t]- F3o2/M_ chin1[_t] * M_ chix [_t]* M_ chiz[_t]) + 
		    2 *M_ gupyz[_t]* (M_ fyz[_t]- F3o2/M_ chin1[_t] * M_ chiy[_t]* M_ chiz[_t]);
		// M_ Add M_ chi part toM_ Ricci tensor:

		M_ Rxx [_t]=M_ Rxx [_t]+ (M_ fxx [_t]- M_ chix[_t]*M_ chix[_t]/M_ chin1[_t]/2 +M_ gxx [_t]*M_ f[_t])/M_ chin1[_t]/2;
		M_ Ryy[_t]=M_ Ryy[_t]+ (M_ fyy[_t]- M_ chiy[_t]*M_ chiy[_t]/M_ chin1[_t]/2 +M_ gyy[_t]*M_ f[_t])/M_ chin1[_t]/2;
		M_ Rzz[_t]=M_ Rzz[_t]+ (M_ fzz[_t]- M_ chiz[_t]*M_ chiz[_t]/M_ chin1[_t]/2 +M_ gzz[_t]*M_ f[_t])/M_ chin1[_t]/2;
		M_ Rxy[_t]=M_ Rxy[_t]+ (M_ fxy[_t]- M_ chix[_t]*M_ chiy[_t]/M_ chin1[_t]/2 +M_ gxy[_t]*M_ f[_t])/M_ chin1[_t]/2;
		M_ Rxz[_t]=M_ Rxz[_t]+ (M_ fxz[_t]- M_ chix[_t]*M_ chiz[_t]/M_ chin1[_t]/2 +M_ gxz[_t]*M_ f[_t])/M_ chin1[_t]/2;
		M_ Ryz[_t]=M_ Ryz[_t]+ (M_ fyz[_t]- M_ chiy[_t]*M_ chiz[_t]/M_ chin1[_t]/2 +M_ gyz[_t]*M_ f[_t])/M_ chin1[_t]/2;
		
			
		_t += STEP_SIZE;
	}
}

__global__ void compute_rhs_bssn_part6() 
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
		M_ gxxx [_t]= (M_ gupxx [_t]* M_ chix [_t]+M_ gupxy[_t]* M_ chiy[_t]+M_ gupxz[_t]* M_ chiz[_t])/M_ chin1[_t];
		M_ gxxy[_t]= (M_ gupxy[_t]* M_ chix [_t]+M_ gupyy[_t]* M_ chiy[_t]+M_ gupyz[_t]* M_ chiz[_t])/M_ chin1[_t];
		M_ gxxz[_t]= (M_ gupxz[_t]* M_ chix [_t]+M_ gupyz[_t]* M_ chiy[_t]+M_ gupzz[_t]* M_ chiz[_t])/M_ chin1[_t];
		// nowM_ get physical second kind of connection
		M_ Gamxxx [_t]= M_ Gamxxx [_t]- ( (M_ chix [_t]+ M_ chix[_t])/M_ chin1[_t] -M_ gxx [_t]*M_ gxxx [_t])*HALF;
		M_ Gamyxx [_t]= M_ Gamyxx [_t]- (                     -M_ gxx [_t]*M_ gxxy[_t])*HALF;
		M_ Gamzxx [_t]= M_ Gamzxx [_t]- (                     -M_ gxx [_t]*M_ gxxz[_t])*HALF;
		M_ Gamxyy[_t]= M_ Gamxyy[_t]- (                     -M_ gyy[_t]*M_ gxxx [_t])*HALF;
		M_ Gamyyy[_t]= M_ Gamyyy[_t]- ( (M_ chiy[_t]+ M_ chiy[_t])/M_ chin1[_t] -M_ gyy[_t]*M_ gxxy[_t])*HALF;
		M_ Gamzyy[_t]= M_ Gamzyy[_t]- (                     -M_ gyy[_t]*M_ gxxz[_t])*HALF;
		M_ Gamxzz[_t]= M_ Gamxzz[_t]- (                     -M_ gzz[_t]*M_ gxxx [_t])*HALF;
		M_ Gamyzz[_t]= M_ Gamyzz[_t]- (                     -M_ gzz[_t]*M_ gxxy[_t])*HALF;
		M_ Gamzzz[_t]= M_ Gamzzz[_t]- ( (M_ chiz[_t]+ M_ chiz[_t])/M_ chin1[_t] -M_ gzz[_t]*M_ gxxz[_t])*HALF;
		M_ Gamxxy[_t]= M_ Gamxxy[_t]- (  M_ chiy[_t]       /M_ chin1[_t] -M_ gxy[_t]*M_ gxxx [_t])*HALF;
		M_ Gamyxy[_t]= M_ Gamyxy[_t]- (         M_ chix [_t]/M_ chin1[_t] -M_ gxy[_t]*M_ gxxy[_t])*HALF;
		M_ Gamzxy[_t]= M_ Gamzxy[_t]- (                     -M_ gxy[_t]*M_ gxxz[_t])*HALF;
		M_ Gamxxz[_t]= M_ Gamxxz[_t]- (  M_ chiz[_t]       /M_ chin1[_t] -M_ gxz[_t]*M_ gxxx [_t])*HALF;
		M_ Gamyxz[_t]= M_ Gamyxz[_t]- (                     -M_ gxz[_t]*M_ gxxy[_t])*HALF;
		M_ Gamzxz[_t]= M_ Gamzxz[_t]- (         M_ chix [_t]/M_ chin1[_t] -M_ gxz[_t]*M_ gxxz[_t])*HALF;
		M_ Gamxyz[_t]= M_ Gamxyz[_t]- (                     -M_ gyz[_t]*M_ gxxx [_t])*HALF;
		M_ Gamyyz[_t]= M_ Gamyyz[_t]- (  M_ chiz[_t]       /M_ chin1[_t] -M_ gyz[_t]*M_ gxxy[_t])*HALF;
		M_ Gamzyz[_t]= M_ Gamzyz[_t]- (         M_ chiy[_t]/M_ chin1[_t] -M_ gyz[_t]*M_ gxxz[_t])*HALF;

		M_ fxx [_t]=M_ fxx [_t]- M_ Gamxxx[_t]*M_ Lapx [_t]- M_ Gamyxx[_t]*M_ Lapy[_t]- M_ Gamzxx[_t]*M_ Lapz[_t];
		M_ fyy[_t]=M_ fyy[_t]- M_ Gamxyy[_t]*M_ Lapx [_t]- M_ Gamyyy[_t]*M_ Lapy[_t]- M_ Gamzyy[_t]*M_ Lapz[_t];
		M_ fzz[_t]=M_ fzz[_t]- M_ Gamxzz[_t]*M_ Lapx [_t]- M_ Gamyzz[_t]*M_ Lapy[_t]- M_ Gamzzz[_t]*M_ Lapz[_t];
		M_ fxy[_t]=M_ fxy[_t]- M_ Gamxxy[_t]*M_ Lapx [_t]- M_ Gamyxy[_t]*M_ Lapy[_t]- M_ Gamzxy[_t]*M_ Lapz[_t];
		M_ fxz[_t]=M_ fxz[_t]- M_ Gamxxz[_t]*M_ Lapx [_t]- M_ Gamyxz[_t]*M_ Lapy[_t]- M_ Gamzxz[_t]*M_ Lapz[_t];
		M_ fyz[_t]=M_ fyz[_t]- M_ Gamxyz[_t]*M_ Lapx [_t]- M_ Gamyyz[_t]*M_ Lapy[_t]- M_ Gamzyz[_t]*M_ Lapz[_t];

		// store D^i D_i Lap in M_ trK_rhs[_t] upto M_ chi
		M_ trK_rhs[_t] =   M_ gupxx [_t]*M_ fxx [_t]+M_ gupyy[_t]*M_ fyy[_t]+M_ gupzz[_t]*M_ fzz[_t]+ 
		2* (M_ gupxy[_t]*M_ fxy[_t]+M_ gupxz[_t]*M_ fxz[_t]+M_ gupyz[_t]*M_ fyz[_t]);
		// M_ Add lapse and M_ S_ij parts toM_ Ricci tensor:
		
		//follow bam code
		M_ S[_t] =  M_ chin1[_t] * ( M_ gupxx[_t] * M_ Sxx[_t] + M_ gupyy[_t] * M_ Syy[_t] + M_ gupzz[_t] * M_ Szz[_t] + 

     2 * ( M_ gupxy[_t] * M_ Sxy[_t] + M_ gupxz[_t] * M_ Sxz[_t] + M_ gupyz[_t] * M_ Syz[_t] ) );
     

M_  f[_t] = F2o3 * M_ trK[_t] * M_ trK[_t] -(

       M_ gupxx[_t] * ( 

       M_ gupxx[_t] * M_ Axx[_t] * M_ Axx[_t] + M_ gupyy[_t] * M_ Axy[_t] * M_ Axy[_t] + M_ gupzz[_t] * M_ Axz[_t] * M_ Axz[_t] + 

       2 * (M_ gupxy[_t] * M_ Axx[_t] * M_ Axy[_t] + M_ gupxz[_t] * M_ Axx[_t] * M_ Axz[_t] + M_ gupyz[_t] * M_ Axy[_t] * M_ Axz[_t]) ) + 

       M_ gupyy[_t] * ( 

       M_ gupxx[_t] * M_ Axy[_t] * M_ Axy[_t] + M_ gupyy[_t] * M_ Ayy[_t] * M_ Ayy[_t] + M_ gupzz[_t] * M_ Ayz[_t] * M_ Ayz[_t] + 

       2 * (M_ gupxy[_t] * M_ Axy[_t] * M_ Ayy[_t] + M_ gupxz[_t] * M_ Axy[_t] * M_ Ayz[_t] + M_ gupyz[_t] * M_ Ayy[_t] * M_ Ayz[_t]) ) + 

       M_ gupzz[_t] * ( 

       M_ gupxx[_t] * M_ Axz[_t] * M_ Axz[_t] + M_ gupyy[_t] * M_ Ayz[_t] * M_ Ayz[_t] + M_ gupzz[_t] * M_ Azz[_t] * M_ Azz[_t] + 

       2 * (M_ gupxy[_t] * M_ Axz[_t] * M_ Ayz[_t] + M_ gupxz[_t] * M_ Axz[_t] * M_ Azz[_t] + M_ gupyz[_t] * M_ Ayz[_t] * M_ Azz[_t]) ) + 

       2 * ( 

       M_ gupxy[_t] * ( 

       M_ gupxx[_t] * M_ Axx[_t] * M_ Axy[_t] + M_ gupyy[_t] * M_ Axy[_t] * M_ Ayy[_t] + M_ gupzz[_t] * M_ Axz[_t] * M_ Ayz[_t] + 

       M_ gupxy[_t] * (M_ Axx[_t] * M_ Ayy[_t] + M_ Axy[_t] * M_ Axy[_t]) + 

       M_ gupxz[_t] * (M_ Axx[_t] * M_ Ayz[_t] + M_ Axz[_t] * M_ Axy[_t]) + 

       M_ gupyz[_t] * (M_ Axy[_t] * M_ Ayz[_t] + M_ Axz[_t] * M_ Ayy[_t]) ) + 

       M_ gupxz[_t] * ( 

       M_ gupxx[_t] * M_ Axx[_t] * M_ Axz[_t] + M_ gupyy[_t] * M_ Axy[_t] * M_ Ayz[_t] + M_ gupzz[_t] * M_ Axz[_t] * M_ Azz[_t] + 

       M_ gupxy[_t] * (M_ Axx[_t] * M_ Ayz[_t] + M_ Axy[_t] * M_ Axz[_t]) + 

       M_ gupxz[_t] * (M_ Axx[_t] * M_ Azz[_t] + M_ Axz[_t] * M_ Axz[_t]) + 

       M_ gupyz[_t] * (M_ Axy[_t] * M_ Azz[_t] + M_ Axz[_t] * M_ Ayz[_t]) ) + 

       M_ gupyz[_t] * ( 

       M_ gupxx[_t] * M_ Axy[_t] * M_ Axz[_t] + M_ gupyy[_t] * M_ Ayy[_t] * M_ Ayz[_t] + M_ gupzz[_t] * M_ Ayz[_t] * M_ Azz[_t] + 

       M_ gupxy[_t] * (M_ Axy[_t] * M_ Ayz[_t] + M_ Ayy[_t] * M_ Axz[_t]) + 

       M_ gupxz[_t] * (M_ Axy[_t] * M_ Azz[_t] + M_ Ayz[_t] * M_ Axz[_t]) + 

       M_ gupyz[_t] * (M_ Ayy[_t] * M_ Azz[_t] + M_ Ayz[_t] * M_ Ayz[_t]) ) )) -16 * PI * M_ rho[_t] + 8 * PI * M_ S[_t];
       

  M_ f[_t] = - F1o3 *(  M_ gupxx[_t] * M_ fxx[_t] + M_ gupyy[_t] * M_ fyy[_t] + M_ gupzz[_t] * M_ fzz[_t] + 

        2* ( M_ gupxy[_t] * M_ fxy[_t] + M_ gupxz[_t] * M_ fxz[_t] + M_ gupyz[_t] * M_ fyz[_t] ) + M_ alpn1[_t] / M_ chin1[_t] * M_ f[_t]);

  

  M_ fxx[_t] = M_ alpn1[_t] * (M_ Rxx[_t] - 8 * PI * M_ Sxx[_t]) - M_ fxx[_t];

  M_ fxy[_t] = M_ alpn1[_t] * (M_ Rxy[_t] - 8 * PI * M_ Sxy[_t]) - M_ fxy[_t];

  M_ fxz[_t] = M_ alpn1[_t] * (M_ Rxz[_t] - 8 * PI * M_ Sxz[_t]) - M_ fxz[_t];

  M_ fyy[_t] = M_ alpn1[_t] * (M_ Ryy[_t] - 8 * PI * M_ Syy[_t]) - M_ fyy[_t];

  M_ fyz[_t] = M_ alpn1[_t] * (M_ Ryz[_t] - 8 * PI * M_ Syz[_t]) - M_ fyz[_t];

  M_ fzz[_t] = M_ alpn1[_t] * (M_ Rzz[_t] - 8 * PI * M_ Szz[_t]) - M_ fzz[_t];
		/*
		M_ fxx [_t]= M_ alpn1[_t]* (M_ Rxx [_t]- 8 * PI * M_ Sxx[_t]) -M_ fxx[_t];
		M_ fxy[_t]= M_ alpn1[_t]* (M_ Rxy[_t]- 8 * PI * M_ Sxy[_t]) -M_ fxy[_t];
		M_ fxz[_t]= M_ alpn1[_t]* (M_ Rxz[_t]- 8 * PI * M_ Sxz[_t]) -M_ fxz[_t];
		M_ fyy[_t]= M_ alpn1[_t]* (M_ Ryy[_t]- 8 * PI * M_ Syy[_t]) -M_ fyy[_t];
		M_ fyz[_t]= M_ alpn1[_t]* (M_ Ryz[_t]- 8 * PI * M_ Syz[_t]) -M_ fyz[_t];
		M_ fzz[_t]= M_ alpn1[_t]* (M_ Rzz[_t]- 8 * PI * M_ Szz[_t]) -M_ fzz[_t];

		// Compute trace-free part (note: M_ chi^-1 and M_ chi cancel//):

		M_ f[_t] = F1o3 *( M_ gupxx [_t]*M_ fxx [_t]+M_ gupyy[_t]*M_ fyy[_t]+M_ gupzz[_t]*M_ fzz[_t]+ 
		2* (M_ gupxy[_t]*M_ fxy[_t]+M_ gupxz[_t]*M_ fxz[_t]+M_ gupyz[_t]*M_ fyz[_t]) );
		*/
		M_ Axx_rhs[_t] =M_ fxx [_t]-M_ gxx [_t]*M_ f[_t];
		M_ Ayy_rhs[_t] =M_ fyy[_t]-M_ gyy[_t]*M_ f[_t];
		M_ Azz_rhs[_t] =M_ fzz[_t]-M_ gzz[_t]*M_ f[_t];
		M_ Axy_rhs[_t] =M_ fxy[_t]-M_ gxy[_t]*M_ f[_t];
		M_ Axz_rhs[_t] =M_ fxz[_t]-M_ gxz[_t]*M_ f[_t];
		M_ Ayz_rhs[_t] =M_ fyz[_t]-M_ gyz[_t]*M_ f[_t];

		// Now: store M_ A_il M_ A^l_j intoM_ fij:

		M_ fxx [_t]=      M_ gupxx [_t]* M_ Axx [_t]* M_ Axx [_t]+M_ gupyy[_t]* M_ Axy[_t]* M_ Axy[_t]+M_ gupzz[_t]* M_ Axz[_t]* M_ Axz[_t]+ 
		2 * (M_ gupxy[_t]* M_ Axx [_t]* M_ Axy[_t]+M_ gupxz[_t]* M_ Axx [_t]* M_ Axz[_t]+M_ gupyz[_t]* M_ Axy[_t]* M_ Axz[_t]);

		M_ fyy[_t]=      M_ gupxx [_t]* M_ Axy[_t]* M_ Axy[_t]+M_ gupyy[_t]* M_ Ayy[_t]* M_ Ayy[_t]+M_ gupzz[_t]* M_ Ayz[_t]* M_ Ayz[_t]+ 
		2 * (M_ gupxy[_t]* M_ Axy[_t]* M_ Ayy[_t]+M_ gupxz[_t]* M_ Axy[_t]* M_ Ayz[_t]+M_ gupyz[_t]* M_ Ayy[_t]* M_ Ayz[_t]);

		M_ fzz[_t]=      M_ gupxx [_t]* M_ Axz[_t]* M_ Axz[_t]+M_ gupyy[_t]* M_ Ayz[_t]* M_ Ayz[_t]+M_ gupzz[_t]* M_ Azz[_t]* M_ Azz[_t]+ 
		2 * (M_ gupxy[_t]* M_ Axz[_t]* M_ Ayz[_t]+M_ gupxz[_t]* M_ Axz[_t]* M_ Azz[_t]+M_ gupyz[_t]* M_ Ayz[_t]* M_ Azz[_t]);

		M_ fxy[_t]=      M_ gupxx [_t]* M_ Axx [_t]* M_ Axy[_t]+M_ gupyy[_t]* M_ Axy[_t]* M_ Ayy[_t]+M_ gupzz[_t]* M_ Axz[_t]* M_ Ayz[_t]+ 
		        M_ gupxy[_t]*(M_ Axx [_t]* M_ Ayy[_t]+ M_ Axy[_t]* M_ Axy[_t])                            + 
		        M_ gupxz[_t]*(M_ Axx [_t]* M_ Ayz[_t]+ M_ Axz[_t]* M_ Axy[_t])                            + 
		        M_ gupyz[_t]*(M_ Axy[_t]* M_ Ayz[_t]+ M_ Axz[_t]* M_ Ayy[_t]);
		M_ fxz[_t]=      M_ gupxx [_t]* M_ Axx [_t]* M_ Axz[_t]+M_ gupyy[_t]* M_ Axy[_t]* M_ Ayz[_t]+M_ gupzz[_t]* M_ Axz[_t]* M_ Azz[_t]+ 
		        M_ gupxy[_t]*(M_ Axx [_t]* M_ Ayz[_t]+ M_ Axy[_t]* M_ Axz[_t])                            + 
		        M_ gupxz[_t]*(M_ Axx [_t]* M_ Azz[_t]+ M_ Axz[_t]* M_ Axz[_t])                            + 
		        M_ gupyz[_t]*(M_ Axy[_t]* M_ Azz[_t]+ M_ Axz[_t]* M_ Ayz[_t]);
		M_ fyz[_t]=      M_ gupxx [_t]* M_ Axy[_t]* M_ Axz[_t]+M_ gupyy[_t]* M_ Ayy[_t]* M_ Ayz[_t]+M_ gupzz[_t]* M_ Ayz[_t]* M_ Azz[_t]+ 
		        M_ gupxy[_t]*(M_ Axy[_t]* M_ Ayz[_t]+ M_ Ayy[_t]* M_ Axz[_t])                            + 
		        M_ gupxz[_t]*(M_ Axy[_t]* M_ Azz[_t]+ M_ Ayz[_t]* M_ Axz[_t])                            + 
		        M_ gupyz[_t]*(M_ Ayy[_t]* M_ Azz[_t]+ M_ Ayz[_t]* M_ Ayz[_t]);

		M_ f[_t] = M_ chin1[_t];
		// store D^i D_i Lap in M_ trK_rhs[_t]
		M_ trK_rhs[_t] =M_ f[_t]*M_ trK_rhs[_t];
		    
		M_ Axx_rhs[_t] =          M_ f[_t] * M_ Axx_rhs[_t]+ M_ alpn1[_t]* (M_ trK[_t]* M_ Axx [_t]- 2 *M_ fxx[_t])  + 
		    2 * (  M_ Axx [_t]* M_ betaxx [_t]+   M_ Axy[_t]* M_ betayx [_t]+   M_ Axz[_t]* M_ betazx [_t])- 
		        F2o3 * M_ Axx [_t]* M_ div_beta[_t];

		M_ Ayy_rhs[_t] =          M_ f[_t] * M_ Ayy_rhs[_t]+ M_ alpn1[_t]* (M_ trK[_t]* M_ Ayy[_t]- 2 *M_ fyy[_t])  + 
		    2 * (  M_ Axy[_t]* M_ betaxy[_t]+   M_ Ayy[_t]* M_ betayy[_t]+   M_ Ayz[_t]* M_ betazy[_t])- 
		        F2o3 * M_ Ayy[_t]* M_ div_beta[_t];

		M_ Azz_rhs[_t] =          M_ f[_t] * M_ Azz_rhs[_t]+ M_ alpn1[_t]* (M_ trK[_t]* M_ Azz[_t]- 2 *M_ fzz[_t])  + 
		    2 * (  M_ Axz[_t]* M_ betaxz[_t]+   M_ Ayz[_t]* M_ betayz[_t]+   M_ Azz[_t]* M_ betazz[_t])- 
		        F2o3 * M_ Azz[_t]* M_ div_beta[_t];

		M_ Axy_rhs[_t] =          M_ f[_t] * M_ Axy_rhs[_t]+ M_ alpn1[_t]*( M_ trK[_t]* M_ Axy[_t] - 2 *M_ fxy[_t])+ 
		            M_ Axx [_t]* M_ betaxy[_t]                 +   M_ Axz[_t]* M_ betazy[_t] + 
		                                M_ Ayy[_t]* M_ betayx [_t]+   M_ Ayz[_t]* M_ betazx [_t] + 
		        F1o3 * M_ Axy[_t]* M_ div_beta[_t]               -   M_ Axy[_t]* M_ betazz[_t];

		M_ Ayz_rhs[_t] =          M_ f[_t] * M_ Ayz_rhs[_t]+ M_ alpn1[_t]*( M_ trK[_t]* M_ Ayz[_t] - 2 *M_ fyz[_t])+ 
		            M_ Axy[_t]* M_ betaxz[_t]+   M_ Ayy[_t]* M_ betayz[_t]                  + 
		            M_ Axz[_t]* M_ betaxy[_t]                 +   M_ Azz[_t]* M_ betazy[_t] + 
		        F1o3 * M_ Ayz[_t]* M_ div_beta[_t]               -   M_ Ayz[_t]* M_ betaxx[_t];

		M_ Axz_rhs[_t] =          M_ f[_t] * M_ Axz_rhs[_t]+ M_ alpn1[_t]*( M_ trK[_t]* M_ Axz[_t] - 2 *M_ fxz[_t])+ 
		            M_ Axx [_t]* M_ betaxz[_t]+   M_ Axy[_t]* M_ betayz[_t]                  + 
		                                M_ Ayz[_t]* M_ betayx [_t]+   M_ Azz[_t]* M_ betazx [_t] + 
		        F1o3 * M_ Axz[_t]* M_ div_beta[_t]               -   M_ Axz[_t]* M_ betayy[_t]  ;   //rhsM_ for M_ Aij

		// Compute trace of M_ S_ij

		M_ S[_t] = M_ f[_t] * (M_ gupxx [_t]* M_ Sxx [_t]+M_ gupyy[_t]* M_ Syy[_t]+M_ gupzz[_t]* M_ Szz[_t]+ 
		2 * (M_ gupxy[_t]* M_ Sxy[_t]+M_ gupxz[_t]* M_ Sxz[_t]+M_ gupyz[_t]* M_ Syz[_t]) );

		M_ trK_rhs[_t] = - M_ trK_rhs[_t] + M_ alpn1[_t]*( F1o3 * M_ trK[_t]* M_ trK[_t]        + 
		        M_ gupxx [_t]*M_ fxx [_t]+M_ gupyy[_t]*M_ fyy[_t]+M_ gupzz[_t]*M_ fzz[_t]  + 
		2 * (M_ gupxy[_t]*M_ fxy[_t]+M_ gupxz[_t]*M_ fxz[_t]+M_ gupyz[_t]*M_ fyz[_t]) + 
		4 * PI * ( M_ rho[_t] + M_ S[_t] )) ;                               //rhsM_ for M_ trK[_t]

		////////M_ gauge variable part

		M_ Lap_rhs[_t] = -2*M_ alpn1[_t] * M_ trK[_t];
		
#if (GAUGE == 0)
		M_ betax_rhs[_t] =0.75*M_ dtSfx[_t];
		M_ betay_rhs[_t] =0.75*M_ dtSfy[_t];
		M_ betaz_rhs[_t] =0.75*M_ dtSfz[_t];

		M_ dtSfx_rhs[_t] = M_ Gamx_rhs[_t] -2*M_ dtSfx[_t];
		M_ dtSfy_rhs[_t] = M_ Gamy_rhs[_t] -2*M_ dtSfy[_t];
		M_ dtSfz_rhs[_t] = M_ Gamz_rhs[_t] -2*M_ dtSfz[_t];
		
#elif (GAUGE == 1)
		M_ betax_rhs[_t] =M_ Gamx[_t] - 2 * M_ betax[_t] ;

  		M_ betay_rhs[_t] =M_ Gamy[_t] - 2 * M_ betay[_t] ;

  		M_ betaz_rhs[_t] =M_ Gamz[_t] - 2 * M_ betaz[_t] ;

		M_ dtSfx_rhs[_t] = 0;
		M_ dtSfy_rhs[_t] = 0;
		M_ dtSfz_rhs[_t] = 0;
		
#elif (GAUGE == 2 || GAUGE == 3)

  M_ betax_rhs[_t] = 0.75* M_ dtSfx[_t];

  M_ betay_rhs[_t] = 0.75* M_ dtSfy[_t];

  M_ betaz_rhs[_t] = 0.75* M_ dtSfz[_t];
  	
#elif (GAUGE == 6)
		if(BHN==2)
		{
			int k = _t / _2D_SIZE[0];
			int ps = _t - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
			int j = ps / ex_c[0];
			int i = ps  - (j * ex_c[0]);
			
			r1 = ( pow2((Porg[0]-X[i]))+ pow2((Porg[1]-Y[j]))+ pow2((Porg[2]-Z[k])) ) / 

          		( pow2((Porg[0]-Porg[3]))+ pow2((Porg[1]-Porg[4])) + pow2((Porg[2]-Porg[5])) );
          

		    r2 = ( pow2((Porg[3]-X[i])) + pow2((Porg[4]-Y[j])) + pow2((Porg[5]-Z[k])) )/ 

		          ( pow2((Porg[0]-Porg[3])) + pow2((Porg[1]-Porg[4])) + pow2((Porg[2]-Porg[5])) );
          

     		reta[i+ j*_1D_SIZE[0]+ k*_2D_SIZE[0] ] = A + C1/(1 + 12 * r1) + C2/(1 + 12 *r2);
		}//BHN == 2
		
		M_ betax_rhs[_t] = 0.75*M_ dtSfx[_t];

		M_ betay_rhs[_t] = 0.75*M_ dtSfy[_t];

		M_ betaz_rhs[_t] = 0.75*M_ dtSfz[_t];



		M_ dtSfx_rhs[_t] = M_ Gamx_rhs[_t] - M_ reta[_t] * M_ dtSfx[_t];

		M_ dtSfy_rhs[_t] = M_ Gamy_rhs[_t] - M_ reta[_t] * M_ dtSfy[_t];

		M_ dtSfz_rhs[_t] = M_ Gamz_rhs[_t] - M_ reta[_t] * M_ dtSfz[_t];
		
#elif (GAUGE == 7)
		if(BHN==2){
			int k = _t / _2D_SIZE[0];
			int ps = _t - (_2D_SIZE[0] * k); //TOTRY: = curr % _2D_SIZE[0];
			int j = ps / ex_c[0];
			int i = ps  - (j * ex_c[0]);
			
			r1 = ( pow2((Porg[0]-X[i])) + pow2((Porg[1]-Y[j])) + pow2((Porg[2]-Z[k])) )/ 

          		( pow2((Porg[0]-Porg[3])) + pow2((Porg[1]-Porg[4])) + pow2((Porg[2]-Porg[5])) );
          

		    r2 = ( pow2((Porg[3]-X[i])) + pow2((Porg[4]-Y[j])) + pow2((Porg[5]-Z[k])) )/ 

		          ( pow2((Porg[0]-Porg[3])) + pow2((Porg[1]-Porg[4])) + pow2((Porg[2]-Porg[5])) );
		          

		    M_ reta[_t][i+ j*_1D_SIZE[0]+ k*_2D_SIZE[0] ] = A + C1* exp(-12 *r1) + C2*exp(- 12*r2);
		}//BHN ==2
		
		M_ betax_rhs[_t] = 0.75*M_ dtSfx[_t];

		M_ betay_rhs[_t] = 0.75*M_ dtSfy[_t];

		M_ betaz_rhs[_t] = 0.75*M_ dtSfz[_t];



		M_ dtSfx_rhs[_t] = M_ Gamx_rhs[_t] - M_ reta[_t]*M_ dtSfx[_t];

		M_ dtSfy_rhs[_t] = M_ Gamy_rhs[_t] - M_ reta[_t]*M_ dtSfy[_t];

		M_ dtSfz_rhs[_t] = M_ Gamz_rhs[_t] - M_ reta[_t]*M_ dtSfz[_t];

#endif  //if (GAUGE == ?)

		_t += STEP_SIZE;
	}
}	

__global__ void compute_rhs_bssn_part6_gauge() 
{	
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
#if (GAUGE == 2)
	M_ reta[_t] =  M_ gupxx[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfx_rhs[_t] +  M_ gupyy[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfy_rhs[_t] +  M_ gupzz[_t] *  M_ dtSfz_rhs[_t] *  M_ dtSfz_rhs[_t] + 

       2 * ( M_ gupxy[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfy_rhs[_t] +  M_ gupxz[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfz_rhs[_t] +  M_ gupyz[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfz_rhs[_t]);
       

   M_ reta[_t] = 1.13 / 2 * sqrt( M_ reta[_t]/M_ chin1[_t])/ pow2( ( 1-sqrt(M_ chin1[_t]) ) );
  

   M_ dtSfx_rhs[_t] =  M_ Gamx_rhs[_t] -  M_ reta[_t]* M_ dtSfx[_t];

   M_ dtSfy_rhs[_t] =  M_ Gamy_rhs[_t] -  M_ reta[_t]* M_ dtSfy[_t];

   M_ dtSfz_rhs[_t] =  M_ Gamz_rhs[_t] -  M_ reta[_t]* M_ dtSfz[_t];

#elif (GAUGE == 3)
	M_ reta[_t] =  M_ gupxx[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfx_rhs[_t] +  M_ gupyy[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfy_rhs[_t]
					 +  M_ gupzz[_t] *  M_ dtSfz_rhs[_t] *  M_ dtSfz_rhs[_t] + 

       				2 * ( M_ gupxy[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfy_rhs[_t] +  
       				M_ gupxz[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfz_rhs[_t] +  
       				M_ gupyz[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfz_rhs[_t]);
       

	M_ reta[_t] = 1.13/2 * sqrt( M_ reta[_t]/ M_ chin1[_t])/ pow2((1-M_ chin1[_t]));

	M_ dtSfx_rhs[_t] =  M_ Gamx_rhs[_t] -  M_ reta[_t]* M_ dtSfx[_t];

	M_ dtSfy_rhs[_t] =  M_ Gamy_rhs[_t] -  M_ reta[_t]* M_ dtSfy[_t];

	M_ dtSfz_rhs[_t] =  M_ Gamz_rhs[_t] -  M_ reta[_t]* M_ dtSfz[_t];
	
#elif (GAUGE == 4)
	M_ reta[_t] =  M_ gupxx[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfx_rhs[_t] +  M_ gupyy[_t] *  M_ dtSfy_rhs[_t] *
				  M_ dtSfy_rhs[_t] +  M_ gupzz[_t] *  M_ dtSfz_rhs[_t] *  M_ dtSfz_rhs[_t] + 

       				2 * ( M_ gupxy[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfy_rhs[_t] +  M_ gupxz[_t] *  
       				M_ dtSfx_rhs[_t] *  M_ dtSfz_rhs[_t] +  M_ gupyz[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfz_rhs[_t]);
       

	M_ reta[_t] = 1.13 / 2 * sqrt( M_ reta[_t]/M_ chin1[_t])/ pow( (1-sqrt(M_ chin1[_t])));


	M_ betax_rhs[_t] = 0.75* M_ Gamx[_t] -  M_ reta[_t]*M_ betax[_t];

	M_ betay_rhs[_t] = 0.75* M_ Gamy[_t] -  M_ reta[_t]*M_ betay[_t];

	M_ betaz_rhs[_t] = 0.75* M_ Gamz[_t] -  M_ reta[_t]*M_ betaz[_t];
	
#elif (GAUGE == 5)
	M_ reta[_t] =  M_ gupxx[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfx_rhs[_t] +  M_ gupyy[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfy_rhs[_t] +  M_ gupzz[_t] *  M_ dtSfz_rhs[_t] *  M_ dtSfz_rhs[_t] + 

       2 * ( M_ gupxy[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfy_rhs[_t] +  M_ gupxz[_t] *  M_ dtSfx_rhs[_t] *  M_ dtSfz_rhs[_t] +  M_ gupyz[_t] *  M_ dtSfy_rhs[_t] *  M_ dtSfz_rhs[_t]);
       

	M_ reta[_t] = 1.13 / 2 * sqrt( M_ reta[_t]/M_ chin1)/ pow( (1-M_ chin1[_t]) );

	M_ betax_rhs[_t] = 0.75* M_ Gamx[_t] -  M_ reta[_t]*M_ betax[_t];

	M_ betay_rhs[_t] = 0.75* M_ Gamy[_t] -  M_ reta[_t]*M_ betay[_t];

	M_ betaz_rhs[_t] = 0.75* M_ Gamz[_t] -  M_ reta[_t]*M_ betaz[_t];



	M_ dtSfx_rhs[_t] = 0;

	M_ dtSfy_rhs[_t] = 0;

	M_ dtSfz_rhs[_t] = 0;
#endif
		_t += STEP_SIZE;
	}
}
__global__ void compute_rhs_bssn_part7() 
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
		M_ ham_Res[_t] =   M_ gupxx [_t]* M_ Rxx [_t]+ M_ gupyy[_t]* M_ Ryy[_t]+ M_ gupzz[_t]* M_ Rzz[_t]+ 
						2* ( M_ gupxy[_t]* M_ Rxy[_t]+ M_ gupxz[_t]* M_ Rxz[_t]+ M_ gupyz[_t]* M_ Ryz[_t]);

		M_ ham_Res[_t] = M_ chin1[_t]*M_ ham_Res[_t] + F2o3 * M_ trK[_t] * M_ trK[_t] -(
					M_ gupxx [_t]* ( 
					M_ gupxx [_t]* M_ Axx [_t]* M_ Axx [_t]+ M_ gupyy[_t]* M_ Axy[_t]* M_ Axy[_t]+ M_ gupzz[_t]* M_ Axz[_t]* M_ Axz[_t]+ 
					2 * (M_ gupxy[_t]* M_ Axx [_t]* M_ Axy[_t]+ M_ gupxz[_t]* M_ Axx [_t]* M_ Axz[_t]+ M_ gupyz[_t]* M_ Axy[_t]* M_ Axz[_t]) ) + 
					M_ gupyy[_t]* ( 
					M_ gupxx [_t]* M_ Axy[_t]* M_ Axy[_t]+ M_ gupyy[_t]* M_ Ayy[_t]* M_ Ayy[_t]+ M_ gupzz[_t]* M_ Ayz[_t]* M_ Ayz[_t]+ 
					2 * (M_ gupxy[_t]* M_ Axy[_t]* M_ Ayy[_t]+ M_ gupxz[_t]* M_ Axy[_t]* M_ Ayz[_t]+ M_ gupyz[_t]* M_ Ayy[_t]* M_ Ayz[_t]) ) + 
					M_ gupzz[_t]* ( 
					M_ gupxx [_t]* M_ Axz[_t]* M_ Axz[_t]+ M_ gupyy[_t]* M_ Ayz[_t]* M_ Ayz[_t]+ M_ gupzz[_t]* M_ Azz[_t]* M_ Azz[_t]+ 
					2 * (M_ gupxy[_t]* M_ Axz[_t]* M_ Ayz[_t]+ M_ gupxz[_t]* M_ Axz[_t]* M_ Azz[_t]+ M_ gupyz[_t]* M_ Ayz[_t]* M_ Azz[_t]) ) + 
					2 * ( 
					M_ gupxy[_t]* ( 
					M_ gupxx [_t]* M_ Axx [_t]* M_ Axy[_t]+ M_ gupyy[_t]* M_ Axy[_t]* M_ Ayy[_t]+ M_ gupzz[_t]* M_ Axz[_t]* M_ Ayz[_t]+ 
					M_ gupxy[_t]* (M_ Axx [_t]* M_ Ayy[_t]+ M_ Axy[_t]* M_ Axy[_t]) + 
					M_ gupxz[_t]* (M_ Axx [_t]* M_ Ayz[_t]+ M_ Axz[_t]* M_ Axy[_t]) + 
					M_ gupyz[_t]* (M_ Axy[_t]* M_ Ayz[_t]+ M_ Axz[_t]* M_ Ayy[_t]) ) + 
					M_ gupxz[_t]* ( 
					M_ gupxx [_t]* M_ Axx [_t]* M_ Axz[_t]+ M_ gupyy[_t]* M_ Axy[_t]* M_ Ayz[_t]+ M_ gupzz[_t]* M_ Axz[_t]* M_ Azz[_t]+ 
					M_ gupxy[_t]* (M_ Axx [_t]* M_ Ayz[_t]+ M_ Axy[_t]* M_ Axz[_t]) + 
					M_ gupxz[_t]* (M_ Axx [_t]* M_ Azz[_t]+ M_ Axz[_t]* M_ Axz[_t]) + 
					M_ gupyz[_t]* (M_ Axy[_t]* M_ Azz[_t]+ M_ Axz[_t]* M_ Ayz[_t]) ) + 
					M_ gupyz[_t]* ( 
					M_ gupxx [_t]* M_ Axy[_t]* M_ Axz[_t]+ M_ gupyy[_t]* M_ Ayy[_t]* M_ Ayz[_t]+ M_ gupzz[_t]* M_ Ayz[_t]* M_ Azz[_t]+ 
					M_ gupxy[_t]* (M_ Axy[_t]* M_ Ayz[_t]+ M_ Ayy[_t]* M_ Axz[_t]) + 
					M_ gupxz[_t]* (M_ Axy[_t]* M_ Azz[_t]+ M_ Ayz[_t]* M_ Axz[_t]) + 
					M_ gupyz[_t]* (M_ Ayy[_t]* M_ Azz[_t]+ M_ Ayz[_t]* M_ Ayz[_t]) ) ))- 16 * PI * M_ rho[_t];
		
		_t += STEP_SIZE;
	}
}
__global__ void compute_rhs_bssn_part8() 
{
	int _t = blockIdx.x*blockDim.x+threadIdx.x;
	while(_t < _3D_SIZE[0])
	{
		M_ gxxx [_t]= M_ gxxx [_t]- (  M_ Gamxxx [_t]* M_ Axx [_t]+ M_ Gamyxx [_t]* M_ Axy[_t]+ M_ Gamzxx [_t]* M_ Axz[_t]
		                + M_ Gamxxx [_t]* M_ Axx [_t]+ M_ Gamyxx [_t]* M_ Axy[_t]+ M_ Gamzxx [_t]* M_ Axz[_t]) - M_ chix[_t]*M_ Axx[_t]/M_ chin1[_t];
		                
		M_ gxyx [_t]= M_ gxyx [_t]- (  M_ Gamxxy[_t]* M_ Axx [_t]+ M_ Gamyxy[_t]* M_ Axy[_t]+ M_ Gamzxy[_t]* M_ Axz[_t]
		                + M_ Gamxxx [_t]* M_ Axy[_t]+ M_ Gamyxx [_t]* M_ Ayy[_t]+ M_ Gamzxx [_t]* M_ Ayz[_t]) - M_ chix[_t]*M_ Axy[_t]/M_ chin1[_t];
		                
		M_ gxzx [_t]= M_ gxzx [_t]- (  M_ Gamxxz[_t]* M_ Axx [_t]+ M_ Gamyxz[_t]* M_ Axy[_t]+ M_ Gamzxz[_t]* M_ Axz[_t]
		                + M_ Gamxxx [_t]* M_ Axz[_t]+ M_ Gamyxx [_t]* M_ Ayz[_t]+ M_ Gamzxx [_t]* M_ Azz[_t]) - M_ chix[_t]*M_ Axz[_t]/M_ chin1[_t];
		                
		M_ gyyx [_t]= M_ gyyx [_t]- (  M_ Gamxxy[_t]* M_ Axy[_t]+ M_ Gamyxy[_t]* M_ Ayy[_t]+ M_ Gamzxy[_t]* M_ Ayz[_t]
		                + M_ Gamxxy[_t]* M_ Axy[_t]+ M_ Gamyxy[_t]* M_ Ayy[_t]+ M_ Gamzxy[_t]* M_ Ayz[_t]) - M_ chix[_t]*M_ Ayy[_t]/M_ chin1[_t];
		                
		M_ gyzx [_t]= M_ gyzx [_t]- (  M_ Gamxxz[_t]* M_ Axy[_t]+ M_ Gamyxz[_t]* M_ Ayy[_t]+ M_ Gamzxz[_t]* M_ Ayz[_t]
		                + M_ Gamxxy[_t]* M_ Axz[_t]+ M_ Gamyxy[_t]* M_ Ayz[_t]+ M_ Gamzxy[_t]* M_ Azz[_t]) - M_ chix[_t]*M_ Ayz[_t]/M_ chin1[_t];
		                
		M_ gzzx [_t]= M_ gzzx [_t]- (  M_ Gamxxz[_t]* M_ Axz[_t]+ M_ Gamyxz[_t]* M_ Ayz[_t]+ M_ Gamzxz[_t]* M_ Azz[_t]
		                + M_ Gamxxz[_t]* M_ Axz[_t]+ M_ Gamyxz[_t]* M_ Ayz[_t]+ M_ Gamzxz[_t]* M_ Azz[_t]) - M_ chix[_t]*M_ Azz[_t]/M_ chin1[_t];
		                
		M_ gxxy[_t]= M_ gxxy[_t]- (  M_ Gamxxy[_t]* M_ Axx [_t]+ M_ Gamyxy[_t]* M_ Axy[_t]+ M_ Gamzxy[_t]* M_ Axz[_t]
		                + M_ Gamxxy[_t]* M_ Axx [_t]+ M_ Gamyxy[_t]* M_ Axy[_t]+ M_ Gamzxy[_t]* M_ Axz[_t]) - M_ chiy[_t]*M_ Axx[_t]/M_ chin1[_t];
		                
		M_ gxyy[_t]= M_ gxyy[_t]- (  M_ Gamxyy[_t]* M_ Axx [_t]+ M_ Gamyyy[_t]* M_ Axy[_t]+ M_ Gamzyy[_t]* M_ Axz[_t]
		                + M_ Gamxxy[_t]* M_ Axy[_t]+ M_ Gamyxy[_t]* M_ Ayy[_t]+ M_ Gamzxy[_t]* M_ Ayz[_t]) - M_ chiy[_t]*M_ Axy[_t]/M_ chin1[_t];
		                
		M_ gxzy[_t]= M_ gxzy[_t]- (  M_ Gamxyz[_t]* M_ Axx [_t]+ M_ Gamyyz[_t]* M_ Axy[_t]+ M_ Gamzyz[_t]* M_ Axz[_t]
		                + M_ Gamxxy[_t]* M_ Axz[_t]+ M_ Gamyxy[_t]* M_ Ayz[_t]+ M_ Gamzxy[_t]* M_ Azz[_t]) - M_ chiy[_t]*M_ Axz[_t]/M_ chin1[_t];
		                
		M_ gyyy[_t]= M_ gyyy[_t]- (  M_ Gamxyy[_t]* M_ Axy[_t]+ M_ Gamyyy[_t]* M_ Ayy[_t]+ M_ Gamzyy[_t]* M_ Ayz[_t]
		                + M_ Gamxyy[_t]* M_ Axy[_t]+ M_ Gamyyy[_t]* M_ Ayy[_t]+ M_ Gamzyy[_t]* M_ Ayz[_t]) - M_ chiy[_t]*M_ Ayy[_t]/M_ chin1[_t];
		                
		M_ gyzy[_t]= M_ gyzy[_t]- (  M_ Gamxyz[_t]* M_ Axy[_t]+ M_ Gamyyz[_t]* M_ Ayy[_t]+ M_ Gamzyz[_t]* M_ Ayz[_t]
		                + M_ Gamxyy[_t]* M_ Axz[_t]+ M_ Gamyyy[_t]* M_ Ayz[_t]+ M_ Gamzyy[_t]* M_ Azz[_t]) - M_ chiy[_t]*M_ Ayz[_t]/M_ chin1[_t];
		                
		M_ gzzy[_t]= M_ gzzy[_t]- (  M_ Gamxyz[_t]* M_ Axz[_t]+ M_ Gamyyz[_t]* M_ Ayz[_t]+ M_ Gamzyz[_t]* M_ Azz[_t]
		                + M_ Gamxyz[_t]* M_ Axz[_t]+ M_ Gamyyz[_t]* M_ Ayz[_t]+ M_ Gamzyz[_t]* M_ Azz[_t]) - M_ chiy[_t]*M_ Azz[_t]/M_ chin1[_t];
		                
		M_ gxxz[_t]= M_ gxxz[_t]- (  M_ Gamxxz[_t]* M_ Axx [_t]+ M_ Gamyxz[_t]* M_ Axy[_t]+ M_ Gamzxz[_t]* M_ Axz[_t]
		                + M_ Gamxxz[_t]* M_ Axx [_t]+ M_ Gamyxz[_t]* M_ Axy[_t]+ M_ Gamzxz[_t]* M_ Axz[_t]) - M_ chiz[_t]*M_ Axx[_t]/M_ chin1[_t];
		                
		M_ gxyz[_t]= M_ gxyz[_t]- (  M_ Gamxyz[_t]* M_ Axx [_t]+ M_ Gamyyz[_t]* M_ Axy[_t]+ M_ Gamzyz[_t]* M_ Axz[_t]
		                + M_ Gamxxz[_t]* M_ Axy[_t]+ M_ Gamyxz[_t]* M_ Ayy[_t]+ M_ Gamzxz[_t]* M_ Ayz[_t]) - M_ chiz[_t]*M_ Axy[_t]/M_ chin1[_t];
		                
		M_ gxzz[_t]= M_ gxzz[_t]- (  M_ Gamxzz[_t]* M_ Axx [_t]+ M_ Gamyzz[_t]* M_ Axy[_t]+ M_ Gamzzz[_t]* M_ Axz[_t]
		                + M_ Gamxxz[_t]* M_ Axz[_t]+ M_ Gamyxz[_t]* M_ Ayz[_t]+ M_ Gamzxz[_t]* M_ Azz[_t]) - M_ chiz[_t]*M_ Axz[_t]/M_ chin1[_t];
		                
		M_ gyyz[_t]= M_ gyyz[_t]- (  M_ Gamxyz[_t]* M_ Axy[_t]+ M_ Gamyyz[_t]* M_ Ayy[_t]+ M_ Gamzyz[_t]* M_ Ayz[_t]
		                + M_ Gamxyz[_t]* M_ Axy[_t]+ M_ Gamyyz[_t]* M_ Ayy[_t]+ M_ Gamzyz[_t]* M_ Ayz[_t]) - M_ chiz[_t]*M_ Ayy[_t]/M_ chin1[_t];
		                
		M_ gyzz[_t]= M_ gyzz[_t]- (  M_ Gamxzz[_t]* M_ Axy[_t]+ M_ Gamyzz[_t]* M_ Ayy[_t]+ M_ Gamzzz[_t]* M_ Ayz[_t]
		                + M_ Gamxyz[_t]* M_ Axz[_t]+ M_ Gamyyz[_t]* M_ Ayz[_t]+ M_ Gamzyz[_t]* M_ Azz[_t]) - M_ chiz[_t]*M_ Ayz[_t]/M_ chin1[_t];
		                
		M_ gzzz[_t]= M_ gzzz[_t]- (  M_ Gamxzz[_t]* M_ Axz[_t]+ M_ Gamyzz[_t]* M_ Ayz[_t]+ M_ Gamzzz[_t]* M_ Azz[_t]
		                + M_ Gamxzz[_t]* M_ Axz[_t]+ M_ Gamyzz[_t]* M_ Ayz[_t]+ M_ Gamzzz[_t]* M_ Azz[_t]) - M_ chiz[_t]*M_ Azz[_t]/M_ chin1[_t];
		                
		M_ movx_Res[_t] = M_ gupxx[_t]*M_ gxxx [_t]+ M_ gupyy[_t]*M_ gxyy[_t]+ M_ gupzz[_t]*M_ gxzz[_t]
		        +M_ gupxy[_t]*M_ gxyx [_t]+ M_ gupxz[_t]*M_ gxzx [_t]+ M_ gupyz[_t]*M_ gxzy[_t]
		        +M_ gupxy[_t]*M_ gxxy[_t]+ M_ gupxz[_t]*M_ gxxz[_t]+ M_ gupyz[_t]*M_ gxyz[_t];
		M_ movy_Res[_t] = M_ gupxx[_t]*M_ gxyx [_t]+ M_ gupyy[_t]*M_ gyyy[_t]+ M_ gupzz[_t]*M_ gyzz[_t]
		        +M_ gupxy[_t]*M_ gyyx [_t]+ M_ gupxz[_t]*M_ gyzx [_t]+ M_ gupyz[_t]*M_ gyzy[_t]
		        +M_ gupxy[_t]*M_ gxyy[_t]+ M_ gupxz[_t]*M_ gxyz[_t]+ M_ gupyz[_t]*M_ gyyz[_t];
		        
		M_ movz_Res[_t] = M_ gupxx[_t]*M_ gxzx [_t]+ M_ gupyy[_t]*M_ gyzy[_t]+ M_ gupzz[_t]*M_ gzzz[_t]
		        +M_ gupxy[_t]*M_ gyzx [_t]+ M_ gupxz[_t]*M_ gzzx [_t]+ M_ gupyz[_t]*M_ gzzy[_t]
		        +M_ gupxy[_t]*M_ gxzy[_t]+ M_ gupxz[_t]*M_ gxzz[_t]+ M_ gupyz[_t]*M_ gyzz[_t];

		M_ movx_Res[_t] = M_ movx_Res[_t] - F2o3*M_ Kx [_t]- 8*PI*M_ Sx[_t];
		M_ movy_Res[_t] = M_ movy_Res[_t] - F2o3*M_ Ky[_t]- 8*PI*M_ Sy[_t];
		M_ movz_Res[_t] = M_ movz_Res[_t] - F2o3*M_ Kz[_t]- 8*PI*M_ Sz[_t];
		
		_t += STEP_SIZE;
	}
}



__global__ void device_test(double * result, double * Xt){	
	/*result[0] = MAXSIZE;
	result[1] = STEP;
	result[2] = ex_c[0];
	result[3] = ex_c[1];
	result[4] = ex_c[2];
	result[5] = Xt[0];
	result[6] = Xt[1];
	result[7] = metac.X[0];
	result[8] = metac.X[1];	*/
	
	result[0] = metac.gzz[0];
	result[1] = metac.gzz[1];
	result[2] = metac.gzz[2];
	result[3] = metac.gyy[0];
	result[4] = metac.gyy[1];
	result[5] = metac.gyy[2];
	result[6] = _3D_SIZE[0];
	result[7] = STEP_SIZE;
	result[8] = blockDim.x * gridDim.x;
}

void destroy_meta(Meta *meta)
{
	/*
	if(Mh_ X) CUDA_SAFE_CALL(cudaFree(Mh_ X));
	if(Mh_ Y) CUDA_SAFE_CALL(cudaFree(Mh_ Y));
	if(Mh_ Z) CUDA_SAFE_CALL(cudaFree(Mh_ Z));
	if(Mh_ chi) CUDA_SAFE_CALL(cudaFree(Mh_ chi));
	if(Mh_ dxx) CUDA_SAFE_CALL(cudaFree(Mh_ dxx));
	if(Mh_ dyy) CUDA_SAFE_CALL(cudaFree(Mh_ dyy));
	if(Mh_ dzz) CUDA_SAFE_CALL(cudaFree(Mh_ dzz));
	if(Mh_ trK) CUDA_SAFE_CALL(cudaFree(Mh_ trK));
	if(Mh_ gxy) CUDA_SAFE_CALL(cudaFree(Mh_ gxy));
	if(Mh_ gxz) CUDA_SAFE_CALL(cudaFree(Mh_ gxz));
	if(Mh_ gyz) CUDA_SAFE_CALL(cudaFree(Mh_ gyz));
	if(Mh_ Axx) CUDA_SAFE_CALL(cudaFree(Mh_ Axx));
	if(Mh_ Axy) CUDA_SAFE_CALL(cudaFree(Mh_ Axy));
	if(Mh_ Axz) CUDA_SAFE_CALL(cudaFree(Mh_ Axz));
	if(Mh_ Ayz) CUDA_SAFE_CALL(cudaFree(Mh_ Ayz));
	if(Mh_ Ayy) CUDA_SAFE_CALL(cudaFree(Mh_ Ayy));
	if(Mh_ Azz) CUDA_SAFE_CALL(cudaFree(Mh_ Azz));
	if(Mh_ Gamx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamx));
	if(Mh_ Gamy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamy));
	if(Mh_ Gamz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamz));
	if(Mh_ Lap) CUDA_SAFE_CALL(cudaFree(Mh_ Lap));
	if(Mh_ betax) CUDA_SAFE_CALL(cudaFree(Mh_ betax));
	if(Mh_ betay) CUDA_SAFE_CALL(cudaFree(Mh_ betay));
	if(Mh_ betaz) CUDA_SAFE_CALL(cudaFree(Mh_ betaz));
	if(Mh_ dtSfx) CUDA_SAFE_CALL(cudaFree(Mh_ dtSfx));
	if(Mh_ dtSfy) CUDA_SAFE_CALL(cudaFree(Mh_ dtSfy));
	if(Mh_ dtSfz) CUDA_SAFE_CALL(cudaFree(Mh_ dtSfz));
	if(Mh_ chi_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ chi_rhs));
	if(Mh_ trK_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ trK_rhs));
	if(Mh_ gxy_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ gxy_rhs));
	if(Mh_ gxz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ gxz_rhs));
	if(Mh_ gyz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ gyz_rhs));
	if(Mh_ Axx_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Axx_rhs));
	if(Mh_ Axy_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Axy_rhs));
	if(Mh_ Axz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Axz_rhs));
	if(Mh_ Ayz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Ayz_rhs));
	if(Mh_ Ayy_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Ayy_rhs));
	if(Mh_ Azz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Azz_rhs));
	if(Mh_ Gamx_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Gamx_rhs));
	if(Mh_ Gamy_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Gamy_rhs));
	if(Mh_ Gamz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Gamz_rhs));
	if(Mh_ Lap_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ Lap_rhs));
	if(Mh_ betax_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ betax_rhs));
	if(Mh_ betay_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ betay_rhs));  
	if(Mh_ betaz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ betaz_rhs));
	if(Mh_ dtSfx_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ dtSfx_rhs));
	if(Mh_ dtSfy_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ dtSfy_rhs));
	if(Mh_ dtSfz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ dtSfz_rhs));
	if(Mh_ rho) CUDA_SAFE_CALL(cudaFree(Mh_ rho));
	if(Mh_ Sx) CUDA_SAFE_CALL(cudaFree(Mh_ Sx));
	if(Mh_ Sy) CUDA_SAFE_CALL(cudaFree(Mh_ Sy));
	if(Mh_ Sz) CUDA_SAFE_CALL(cudaFree(Mh_ Sz));
	if(Mh_ Sxx) CUDA_SAFE_CALL(cudaFree(Mh_ Sxx));
	if(Mh_ Sxy) CUDA_SAFE_CALL(cudaFree(Mh_ Sxy));
	if(Mh_ Sxz) CUDA_SAFE_CALL(cudaFree(Mh_ Sxz));
	if(Mh_ Syz) CUDA_SAFE_CALL(cudaFree(Mh_ Syz));
	if(Mh_ Syy) CUDA_SAFE_CALL(cudaFree(Mh_ Syy));
	if(Mh_ Szz) CUDA_SAFE_CALL(cudaFree(Mh_ Szz));
	if(Mh_ Gamxxx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxxx));
	if(Mh_ Gamxxy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxxy));
	if(Mh_ Gamxxz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxxz));
	if(Mh_ Gamxyy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxyy));
	if(Mh_ Gamxyz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxyz));
	if(Mh_ Gamxzz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxzz));
	if(Mh_ Gamyxx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyxx));
	if(Mh_ Gamyxy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyxy));
	if(Mh_ Gamyxz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyxz));
	if(Mh_ Gamyyy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyyy));
	if(Mh_ Gamyyz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyyz));
	if(Mh_ Gamyzz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyzz));
	if(Mh_ Gamzxx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzxx));
	if(Mh_ Gamzxy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzxy));
	if(Mh_ Gamzxz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzxz));
	if(Mh_ Gamzyz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzyz));
	if(Mh_ Gamzyy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzyy));
	if(Mh_ Gamzzz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzzz));
	if(Mh_ Rxx) CUDA_SAFE_CALL(cudaFree(Mh_ Rxx));
	if(Mh_ Rxy) CUDA_SAFE_CALL(cudaFree(Mh_ Rxy));
	if(Mh_ Rxz) CUDA_SAFE_CALL(cudaFree(Mh_ Rxz));
	if(Mh_ Ryy) CUDA_SAFE_CALL(cudaFree(Mh_ Ryy));
	if(Mh_ Ryz) CUDA_SAFE_CALL(cudaFree(Mh_ Ryz));
	if(Mh_ Rzz) CUDA_SAFE_CALL(cudaFree(Mh_ Rzz));
	if(Mh_ ham_Res) CUDA_SAFE_CALL(cudaFree(Mh_ ham_Res));
	if(Mh_ movx_Res) CUDA_SAFE_CALL(cudaFree(Mh_ movx_Res));
	if(Mh_ movy_Res) CUDA_SAFE_CALL(cudaFree(Mh_ movy_Res));
	if(Mh_ movz_Res) CUDA_SAFE_CALL(cudaFree(Mh_ movz_Res));
	if(Mh_ Gmx_Res) CUDA_SAFE_CALL(cudaFree(Mh_ Gmx_Res));
	if(Mh_ Gmy_Res) CUDA_SAFE_CALL(cudaFree(Mh_ Gmy_Res));
	if(Mh_ Gmz_Res) CUDA_SAFE_CALL(cudaFree(Mh_ Gmz_Res));
	if(Mh_ gxx) CUDA_SAFE_CALL(cudaFree(Mh_ gxx));
	if(Mh_ gyy) CUDA_SAFE_CALL(cudaFree(Mh_ gyy));
	if(Mh_ gzz) CUDA_SAFE_CALL(cudaFree(Mh_ gzz));
	if(Mh_ chix) CUDA_SAFE_CALL(cudaFree(Mh_ chix));
	if(Mh_ chiy) CUDA_SAFE_CALL(cudaFree(Mh_ chiy));
	if(Mh_ chiz) CUDA_SAFE_CALL(cudaFree(Mh_ chiz));
	if(Mh_ gxxx) CUDA_SAFE_CALL(cudaFree(Mh_ gxxx));
	if(Mh_ gxyx) CUDA_SAFE_CALL(cudaFree(Mh_ gxyx));
	if(Mh_ gxzx) CUDA_SAFE_CALL(cudaFree(Mh_ gxzx));
	if(Mh_ gyyx) CUDA_SAFE_CALL(cudaFree(Mh_ gyyx));
	if(Mh_ gyzx) CUDA_SAFE_CALL(cudaFree(Mh_ gyzx));
	if(Mh_ gzzx) CUDA_SAFE_CALL(cudaFree(Mh_ gzzx));
	if(Mh_ gxxy) CUDA_SAFE_CALL(cudaFree(Mh_ gxxy));
	if(Mh_ gxyy) CUDA_SAFE_CALL(cudaFree(Mh_ gxyy));
	if(Mh_ gxzy) CUDA_SAFE_CALL(cudaFree(Mh_ gxzy));
	if(Mh_ gyyy) CUDA_SAFE_CALL(cudaFree(Mh_ gyyy));
	if(Mh_ gyzy) CUDA_SAFE_CALL(cudaFree(Mh_ gyzy));
	if(Mh_ gzzy) CUDA_SAFE_CALL(cudaFree(Mh_ gzzy));
	if(Mh_ gxxz) CUDA_SAFE_CALL(cudaFree(Mh_ gxxz));
	if(Mh_ gxyz) CUDA_SAFE_CALL(cudaFree(Mh_ gxyz));
	if(Mh_ gxzz) CUDA_SAFE_CALL(cudaFree(Mh_ gxzz));
	if(Mh_ gyyz) CUDA_SAFE_CALL(cudaFree(Mh_ gyyz));
	if(Mh_ gyzz) CUDA_SAFE_CALL(cudaFree(Mh_ gyzz));
	if(Mh_ gzzz) CUDA_SAFE_CALL(cudaFree(Mh_ gzzz));
	if(Mh_ Lapx) CUDA_SAFE_CALL(cudaFree(Mh_ Lapx));
	if(Mh_ Lapy) CUDA_SAFE_CALL(cudaFree(Mh_ Lapy));
	if(Mh_ Lapz) CUDA_SAFE_CALL(cudaFree(Mh_ Lapz));
	if(Mh_ betaxx) CUDA_SAFE_CALL(cudaFree(Mh_ betaxx));
	if(Mh_ betaxy) CUDA_SAFE_CALL(cudaFree(Mh_ betaxy));
	if(Mh_ betaxz) CUDA_SAFE_CALL(cudaFree(Mh_ betaxz));
	if(Mh_ betayy) CUDA_SAFE_CALL(cudaFree(Mh_ betayy));
	if(Mh_ betayz) CUDA_SAFE_CALL(cudaFree(Mh_ betayz));
	if(Mh_ betazz) CUDA_SAFE_CALL(cudaFree(Mh_ betazz));
	if(Mh_ betayx) CUDA_SAFE_CALL(cudaFree(Mh_ betayx));
	if(Mh_ betazy) CUDA_SAFE_CALL(cudaFree(Mh_ betazy));
	if(Mh_ betazx) CUDA_SAFE_CALL(cudaFree(Mh_ betazx));
	if(Mh_ Kx) CUDA_SAFE_CALL(cudaFree(Mh_ Kx));
	if(Mh_ Ky) CUDA_SAFE_CALL(cudaFree(Mh_ Ky));
	if(Mh_ Kz) CUDA_SAFE_CALL(cudaFree(Mh_ Kz));
	if(Mh_ Gamxx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxx));
	if(Mh_ Gamxy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxy));
	if(Mh_ Gamxz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxz));
	if(Mh_ Gamyy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyy));
	if(Mh_ Gamyz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyz));
	if(Mh_ Gamzz) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzz));
	if(Mh_ Gamyx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamyx));
	if(Mh_ Gamzy) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzy));
	if(Mh_ Gamzx) CUDA_SAFE_CALL(cudaFree(Mh_ Gamzx));
	if(Mh_ div_beta) CUDA_SAFE_CALL(cudaFree(Mh_ div_beta));
	if(Mh_ S) CUDA_SAFE_CALL(cudaFree(Mh_ S));
	if(Mh_ f) CUDA_SAFE_CALL(cudaFree(Mh_ f));
	if(Mh_ fxx) CUDA_SAFE_CALL(cudaFree(Mh_ fxx));
	if(Mh_ fxy) CUDA_SAFE_CALL(cudaFree(Mh_ fxy));
	if(Mh_ fxz) CUDA_SAFE_CALL(cudaFree(Mh_ fxz));
	if(Mh_ fyy) CUDA_SAFE_CALL(cudaFree(Mh_ fyy));
	if(Mh_ fyz) CUDA_SAFE_CALL(cudaFree(Mh_ fyz));
	if(Mh_ fzz) CUDA_SAFE_CALL(cudaFree(Mh_ fzz));
	if(Mh_ gupxx) CUDA_SAFE_CALL(cudaFree(Mh_ gupxx));
	if(Mh_ gupxy) CUDA_SAFE_CALL(cudaFree(Mh_ gupxy));
	if(Mh_ gupxz) CUDA_SAFE_CALL(cudaFree(Mh_ gupxz));
	if(Mh_ gupyy) CUDA_SAFE_CALL(cudaFree(Mh_ gupyy));
	if(Mh_ gupyz) CUDA_SAFE_CALL(cudaFree(Mh_ gupyz));
	if(Mh_ gupzz) CUDA_SAFE_CALL(cudaFree(Mh_ gupzz));
	if(Mh_ Gamxa) CUDA_SAFE_CALL(cudaFree(Mh_ Gamxa));
	if(Mh_ Gamya) CUDA_SAFE_CALL(cudaFree(Mh_ Gamya));
	if(Mh_ Gamza) CUDA_SAFE_CALL(cudaFree(Mh_ Gamza));
	if(Mh_ alpn1) CUDA_SAFE_CALL(cudaFree(Mh_ alpn1));
	if(Mh_ chin1) CUDA_SAFE_CALL(cudaFree(Mh_ chin1));
	if(Mh_ fh) CUDA_SAFE_CALL(cudaFree(Mh_ fh));
	if(Mh_ fh2) CUDA_SAFE_CALL(cudaFree(Mh_ fh2));
	if(Mh_ gxx_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ gxx_rhs));
	if(Mh_ gyy_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ gyy_rhs));
	if(Mh_ gzz_rhs) CUDA_SAFE_CALL(cudaFree(Mh_ gzz_rhs));
	*/
	
	if(Mh_ X) cudaFree(Mh_ X);
	if(Mh_ Y) cudaFree(Mh_ Y);
	if(Mh_ Z) cudaFree(Mh_ Z);
	if(Mh_ chi) cudaFree(Mh_ chi);
	if(Mh_ dxx) cudaFree(Mh_ dxx);
	if(Mh_ dyy) cudaFree(Mh_ dyy);
	if(Mh_ dzz) cudaFree(Mh_ dzz);
	if(Mh_ trK) cudaFree(Mh_ trK);
	if(Mh_ gxy) cudaFree(Mh_ gxy);
	if(Mh_ gxz) cudaFree(Mh_ gxz);
	if(Mh_ gyz) cudaFree(Mh_ gyz);
	if(Mh_ Axx) cudaFree(Mh_ Axx);
	if(Mh_ Axy) cudaFree(Mh_ Axy);
	if(Mh_ Axz) cudaFree(Mh_ Axz);
	if(Mh_ Ayz) cudaFree(Mh_ Ayz);
	if(Mh_ Ayy) cudaFree(Mh_ Ayy);
	if(Mh_ Azz) cudaFree(Mh_ Azz);
	if(Mh_ Gamx) cudaFree(Mh_ Gamx);
	if(Mh_ Gamy) cudaFree(Mh_ Gamy);
	if(Mh_ Gamz) cudaFree(Mh_ Gamz);
	if(Mh_ Lap) cudaFree(Mh_ Lap);
	if(Mh_ betax) cudaFree(Mh_ betax);
	if(Mh_ betay) cudaFree(Mh_ betay);
	if(Mh_ betaz) cudaFree(Mh_ betaz);
	if(Mh_ dtSfx) cudaFree(Mh_ dtSfx);
	if(Mh_ dtSfy) cudaFree(Mh_ dtSfy);
	if(Mh_ dtSfz) cudaFree(Mh_ dtSfz);
	if(Mh_ chi_rhs) cudaFree(Mh_ chi_rhs);
	if(Mh_ trK_rhs) cudaFree(Mh_ trK_rhs);
	if(Mh_ gxy_rhs) cudaFree(Mh_ gxy_rhs);
	if(Mh_ gxz_rhs) cudaFree(Mh_ gxz_rhs);
	if(Mh_ gyz_rhs) cudaFree(Mh_ gyz_rhs);
	if(Mh_ Axx_rhs) cudaFree(Mh_ Axx_rhs);
	if(Mh_ Axy_rhs) cudaFree(Mh_ Axy_rhs);
	if(Mh_ Axz_rhs) cudaFree(Mh_ Axz_rhs);
	if(Mh_ Ayz_rhs) cudaFree(Mh_ Ayz_rhs);
	if(Mh_ Ayy_rhs) cudaFree(Mh_ Ayy_rhs);
	if(Mh_ Azz_rhs) cudaFree(Mh_ Azz_rhs);
	if(Mh_ Gamx_rhs) cudaFree(Mh_ Gamx_rhs);
	if(Mh_ Gamy_rhs) cudaFree(Mh_ Gamy_rhs);
	if(Mh_ Gamz_rhs) cudaFree(Mh_ Gamz_rhs);
	if(Mh_ Lap_rhs) cudaFree(Mh_ Lap_rhs);
	if(Mh_ betax_rhs) cudaFree(Mh_ betax_rhs);
	if(Mh_ betay_rhs) cudaFree(Mh_ betay_rhs);  
	if(Mh_ betaz_rhs) cudaFree(Mh_ betaz_rhs);
	if(Mh_ dtSfx_rhs) cudaFree(Mh_ dtSfx_rhs);
	if(Mh_ dtSfy_rhs) cudaFree(Mh_ dtSfy_rhs);
	if(Mh_ dtSfz_rhs) cudaFree(Mh_ dtSfz_rhs);
	if(Mh_ rho) cudaFree(Mh_ rho);
	if(Mh_ Sx) cudaFree(Mh_ Sx);
	if(Mh_ Sy) cudaFree(Mh_ Sy);
	if(Mh_ Sz) cudaFree(Mh_ Sz);
	if(Mh_ Sxx) cudaFree(Mh_ Sxx);
	if(Mh_ Sxy) cudaFree(Mh_ Sxy);
	if(Mh_ Sxz) cudaFree(Mh_ Sxz);
	if(Mh_ Syz) cudaFree(Mh_ Syz);
	if(Mh_ Syy) cudaFree(Mh_ Syy);
	if(Mh_ Szz) cudaFree(Mh_ Szz);
	if(Mh_ Gamxxx) cudaFree(Mh_ Gamxxx);
	if(Mh_ Gamxxy) cudaFree(Mh_ Gamxxy);
	if(Mh_ Gamxxz) cudaFree(Mh_ Gamxxz);
	if(Mh_ Gamxyy) cudaFree(Mh_ Gamxyy);
	if(Mh_ Gamxyz) cudaFree(Mh_ Gamxyz);
	if(Mh_ Gamxzz) cudaFree(Mh_ Gamxzz);
	if(Mh_ Gamyxx) cudaFree(Mh_ Gamyxx);
	if(Mh_ Gamyxy) cudaFree(Mh_ Gamyxy);
	if(Mh_ Gamyxz) cudaFree(Mh_ Gamyxz);
	if(Mh_ Gamyyy) cudaFree(Mh_ Gamyyy);
	if(Mh_ Gamyyz) cudaFree(Mh_ Gamyyz);
	if(Mh_ Gamyzz) cudaFree(Mh_ Gamyzz);
	if(Mh_ Gamzxx) cudaFree(Mh_ Gamzxx);
	if(Mh_ Gamzxy) cudaFree(Mh_ Gamzxy);
	if(Mh_ Gamzxz) cudaFree(Mh_ Gamzxz);
	if(Mh_ Gamzyz) cudaFree(Mh_ Gamzyz);
	if(Mh_ Gamzyy) cudaFree(Mh_ Gamzyy);
	if(Mh_ Gamzzz) cudaFree(Mh_ Gamzzz);
	if(Mh_ Rxx) cudaFree(Mh_ Rxx);
	if(Mh_ Rxy) cudaFree(Mh_ Rxy);
	if(Mh_ Rxz) cudaFree(Mh_ Rxz);
	if(Mh_ Ryy) cudaFree(Mh_ Ryy);
	if(Mh_ Ryz) cudaFree(Mh_ Ryz);
	if(Mh_ Rzz) cudaFree(Mh_ Rzz);
	if(Mh_ ham_Res) cudaFree(Mh_ ham_Res);
	if(Mh_ movx_Res) cudaFree(Mh_ movx_Res);
	if(Mh_ movy_Res) cudaFree(Mh_ movy_Res);
	if(Mh_ movz_Res) cudaFree(Mh_ movz_Res);
	if(Mh_ Gmx_Res) cudaFree(Mh_ Gmx_Res);
	if(Mh_ Gmy_Res) cudaFree(Mh_ Gmy_Res);
	if(Mh_ Gmz_Res) cudaFree(Mh_ Gmz_Res);
	if(Mh_ gxx) cudaFree(Mh_ gxx);
	if(Mh_ gyy) cudaFree(Mh_ gyy);
	if(Mh_ gzz) cudaFree(Mh_ gzz);
	if(Mh_ chix) cudaFree(Mh_ chix);
	if(Mh_ chiy) cudaFree(Mh_ chiy);
	if(Mh_ chiz) cudaFree(Mh_ chiz);
	if(Mh_ gxxx) cudaFree(Mh_ gxxx);
	if(Mh_ gxyx) cudaFree(Mh_ gxyx);
	if(Mh_ gxzx) cudaFree(Mh_ gxzx);
	if(Mh_ gyyx) cudaFree(Mh_ gyyx);
	if(Mh_ gyzx) cudaFree(Mh_ gyzx);
	if(Mh_ gzzx) cudaFree(Mh_ gzzx);
	if(Mh_ gxxy) cudaFree(Mh_ gxxy);
	if(Mh_ gxyy) cudaFree(Mh_ gxyy);
	if(Mh_ gxzy) cudaFree(Mh_ gxzy);
	if(Mh_ gyyy) cudaFree(Mh_ gyyy);
	if(Mh_ gyzy) cudaFree(Mh_ gyzy);
	if(Mh_ gzzy) cudaFree(Mh_ gzzy);
	if(Mh_ gxxz) cudaFree(Mh_ gxxz);
	if(Mh_ gxyz) cudaFree(Mh_ gxyz);
	if(Mh_ gxzz) cudaFree(Mh_ gxzz);
	if(Mh_ gyyz) cudaFree(Mh_ gyyz);
	if(Mh_ gyzz) cudaFree(Mh_ gyzz);
	if(Mh_ gzzz) cudaFree(Mh_ gzzz);
	if(Mh_ Lapx) cudaFree(Mh_ Lapx);
	if(Mh_ Lapy) cudaFree(Mh_ Lapy);
	if(Mh_ Lapz) cudaFree(Mh_ Lapz);
	if(Mh_ betaxx) cudaFree(Mh_ betaxx);
	if(Mh_ betaxy) cudaFree(Mh_ betaxy);
	if(Mh_ betaxz) cudaFree(Mh_ betaxz);
	if(Mh_ betayy) cudaFree(Mh_ betayy);
	if(Mh_ betayz) cudaFree(Mh_ betayz);
	if(Mh_ betazz) cudaFree(Mh_ betazz);
	if(Mh_ betayx) cudaFree(Mh_ betayx);
	if(Mh_ betazy) cudaFree(Mh_ betazy);
	if(Mh_ betazx) cudaFree(Mh_ betazx);
	if(Mh_ Kx) cudaFree(Mh_ Kx);
	if(Mh_ Ky) cudaFree(Mh_ Ky);
	if(Mh_ Kz) cudaFree(Mh_ Kz);
	if(Mh_ Gamxx) cudaFree(Mh_ Gamxx);
	if(Mh_ Gamxy) cudaFree(Mh_ Gamxy);
	if(Mh_ Gamxz) cudaFree(Mh_ Gamxz);
	if(Mh_ Gamyy) cudaFree(Mh_ Gamyy);
	if(Mh_ Gamyz) cudaFree(Mh_ Gamyz);
	if(Mh_ Gamzz) cudaFree(Mh_ Gamzz);
	if(Mh_ Gamyx) cudaFree(Mh_ Gamyx);
	if(Mh_ Gamzy) cudaFree(Mh_ Gamzy);
	if(Mh_ Gamzx) cudaFree(Mh_ Gamzx);
	if(Mh_ div_beta) cudaFree(Mh_ div_beta);
	if(Mh_ S) cudaFree(Mh_ S);
	if(Mh_ f) cudaFree(Mh_ f);
	if(Mh_ fxx) cudaFree(Mh_ fxx);
	if(Mh_ fxy) cudaFree(Mh_ fxy);
	if(Mh_ fxz) cudaFree(Mh_ fxz);
	if(Mh_ fyy) cudaFree(Mh_ fyy);
	if(Mh_ fyz) cudaFree(Mh_ fyz);
	if(Mh_ fzz) cudaFree(Mh_ fzz);
	if(Mh_ gupxx) cudaFree(Mh_ gupxx);
	if(Mh_ gupxy) cudaFree(Mh_ gupxy);
	if(Mh_ gupxz) cudaFree(Mh_ gupxz);
	if(Mh_ gupyy) cudaFree(Mh_ gupyy);
	if(Mh_ gupyz) cudaFree(Mh_ gupyz);
	if(Mh_ gupzz) cudaFree(Mh_ gupzz);
	if(Mh_ Gamxa) cudaFree(Mh_ Gamxa);
	if(Mh_ Gamya) cudaFree(Mh_ Gamya);
	if(Mh_ Gamza) cudaFree(Mh_ Gamza);
	if(Mh_ alpn1) cudaFree(Mh_ alpn1);
	if(Mh_ chin1) cudaFree(Mh_ chin1);
	if(Mh_ fh) cudaFree(Mh_ fh);
	if(Mh_ fh2) cudaFree(Mh_ fh2);
	if(Mh_ gxx_rhs) cudaFree(Mh_ gxx_rhs);
	if(Mh_ gyy_rhs) cudaFree(Mh_ gyy_rhs);
	if(Mh_ gzz_rhs) cudaFree(Mh_ gzz_rhs);
	
#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)
	// if(Mh_ reta) CUDA_SAFE_CALL(cudaFree(Mh_ reta));
	if(Mh_ reta) cudaFree(Mh_ reta);

#endif
	
	//if(Mh_ other_int) cudaFree(Mh_ other_int);
	//if(Mh_ other_double) cudaFree(Mh_ other_double);
	//cout<<"Address of meta:"<<&meta<<endl;
	
	/*if(meta){
		cout<<"Free gpu meta."<<endl;
		free(meta);
	}*/
}

/*void fetch_data(Meta *meta, int matrix_size)
{
	

}*/

int gpu_rhs(int calledby, int mpi_rank, int *ex, double &T,double *X, double *Y, double *Z,                                     
               double *chi, double *  trK ,                                             
               double *dxx , double *  gxy    ,double *gxz ,double * dyy,double *gyz,double *dzz,     
               double *Axx ,   double *Axy ,  double * Axz ,  double * Ayy ,  double * Ayz , double * Azz,     
               double *Gamx ,  double *Gamy ,  double *Gamz ,                                  
               double *Lap ,  double *betax ,  double *betay ,  double *betaz ,                       
               double *dtSfx,  double *dtSfy ,  double *dtSfz ,                                  
               double *chi_rhs, double *  trK_rhs,                                             
               double *gxx_rhs,  double * gxy_rhs,  double * gxz_rhs,  double * gyy_rhs,   double *gyz_rhs,  double * gzz_rhs, 
               double *Axx_rhs, double *  Axy_rhs,  double * Axz_rhs,  double * Ayy_rhs,   double *Ayz_rhs,   double *Azz_rhs, 
           	   double *Gamx_rhs, double * Gamy_rhs, double * Gamz_rhs,                                 
               double *Lap_rhs,  double *betax_rhs, double * betay_rhs, double * betaz_rhs,                    
               double *dtSfx_rhs, double * dtSfy_rhs, double * dtSfz_rhs,                              
               double *rho,double *Sx,double *Sy,double *Sz,double *Sxx,
               double *Sxy,double *Sxz,double *Syy,double *Syz,double *Szz,                           
               double *Gamxxx,double *Gamxxy,double *Gamxxz,double *Gamxyy,double *Gamxyz,double *Gamxzz,                      
               double *Gamyxx,double *Gamyxy,double *Gamyxz,double *Gamyyy,double *Gamyyz,double *Gamyzz,                      
               double *Gamzxx,double *Gamzxy,double *Gamzxz,double *Gamzyy,double *Gamzyz,double *Gamzzz,                      
               double *Rxx,double *Rxy,double *Rxz,double *Ryy,double *Ryz,double *Rzz,                                        
               double *ham_Res, double *movx_Res, double *movy_Res,double * movz_Res, 
               double * Gmx_Res, double *Gmy_Res,double * Gmz_Res ,
               int & Symmetry,int &Lev, double &eps, int &co)
{
	//#1------------init gpu meta data---------------------
	//cout<<"init GPU meta data\n";

#ifdef DEVICE_ID  
  	// which device to use
  	cudaSetDevice(DEVICE_ID);
#endif

#ifdef DEVICE_ID_BY_PID
	pid_t pid = getpid();
	cudaSetDevice(pid % 2);
	cout<<"My pid= "<<pid<<endl;
#endif

#ifdef DEVICE_ID_BY_MPI_RANK
	cudaSetDevice(mpi_rank % 2);
#endif

#ifdef TIMING  
  struct timeval tvStart, tvEnd;
  struct timeval tv1, tv2;
  gettimeofday(&tvStart, NULL );
  gettimeofday(&tv1, NULL );
#endif

	//int dim = 3;
	int matrix_size = ex[0] * ex[1] * ex[2];
	Meta met;
	Meta * meta = &met;
	
	/*
	//#1--------------------init_gpu_meta(meta,matrix_size)---------------------------

	//1.1 inout
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ X), ex[0] * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Y), ex[1] * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Z), ex[2] * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ chi), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ trK), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Axx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Axy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Axz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ayz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ayy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Azz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Lap), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betax), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betay), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betaz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dtSfx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dtSfy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dtSfz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ chi_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ trK_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxx_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxy_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyy_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gzz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Axx_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Axy_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Axz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ayz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ayy_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Azz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamx_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamy_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Lap_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betax_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betay_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betaz_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dtSfx_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dtSfy_rhs), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ dtSfz_rhs), matrix_size * sizeof(double)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ rho), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Sx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Sy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Sz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Sxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Sxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Sxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Syz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Syy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Szz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Rxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Rxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Rxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ryy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ryz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Rzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ ham_Res), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ movx_Res), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ movy_Res), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ movz_Res), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gmx_Res), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gmy_Res), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gmz_Res), matrix_size * sizeof(double)));

	//1.2 local Data
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ chix), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ chiy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ chiz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxyx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxzx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyyx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyzx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gzzx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxzy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyzy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gzzy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gxzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gyzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gzzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Lapx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Lapy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Lapz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betaxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betaxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betaxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betayy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betayz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betazz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betayx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betazy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ betazx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Kx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Ky), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Kz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamyx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamzx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ div_beta), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ S), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ f), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gupxx), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gupxy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gupxz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gupyy), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gupyz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ gupzz), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamxa), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamya), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ Gamza), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ alpn1), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ chin1), matrix_size * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fh), (ex[0]+2)*(ex[1]+2)*(ex[2]+2) * sizeof(double)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&(Mh_ fh2), (ex[0]+3)*(ex[1]+3)*(ex[2]+3) * sizeof(double)));
	*/
	
	//#1--------------------init_gpu_meta(meta,matrix_size)---------------------------

	//1.1 inout
	cudaMalloc((void**)&(Mh_ X), ex[0] * sizeof(double));
	cudaMalloc((void**)&(Mh_ Y), ex[1] * sizeof(double));
	cudaMalloc((void**)&(Mh_ Z), ex[2] * sizeof(double));
	cudaMalloc((void**)&(Mh_ chi), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ trK), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Axx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Axy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Axz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ayz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ayy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Azz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Lap), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betax), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betay), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betaz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dtSfx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dtSfy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dtSfz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ chi_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ trK_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxx_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxy_rhs), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ gyy_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxz_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyz_rhs), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ gzz_rhs), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ Axx_rhs), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ Axy_rhs), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ Axz_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ayz_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ayy_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Azz_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamx_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamy_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamz_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Lap_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betax_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betay_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betaz_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dtSfx_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dtSfy_rhs), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ dtSfz_rhs), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ rho), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Sx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Sy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Sz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Sxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Sxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Sxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Syz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Syy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Szz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxzz), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ Gamyxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Rxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Rxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Rxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ryy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ryz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Rzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ ham_Res), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ movx_Res), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ movy_Res), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ movz_Res), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gmx_Res), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gmy_Res), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gmz_Res), matrix_size * sizeof(double));


	//1.2 local Data
	cudaMalloc((void**)&(Mh_ gxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ chix), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ chiy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ chiz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxyx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxzx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyyx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyzx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gzzx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxzy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyzy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gzzy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gxzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gyzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gzzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Lapx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Lapy), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ Lapz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betaxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betaxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betaxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betayy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betayz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betazz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betayx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betazy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ betazx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Kx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Ky), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Kz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamyx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamzx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ div_beta), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ S), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ f), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ fxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ fxy), matrix_size * sizeof(double));

	cudaMalloc((void**)&(Mh_ fxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ fyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ fyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ fzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gupxx), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gupxy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gupxz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gupyy), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gupyz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ gupzz), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamxa), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamya), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ Gamza), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ alpn1), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ chin1), matrix_size * sizeof(double));
	cudaMalloc((void**)&(Mh_ fh), (ex[0]+2)*(ex[1]+2)*(ex[2]+2) * sizeof(double));
	cudaMalloc((void**)&(Mh_ fh2), (ex[0]+3)*(ex[1]+3)*(ex[2]+3) * sizeof(double));
	
	#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5 || GAUGE == 6 || GAUGE == 7)

  		cudaMalloc((void**)&(Mh_ reta), matrix_size * sizeof(double));

	#endif
	  
//2 ----------------Copy Data to Device------------------
	cudaMemcpy(Mh_ X, X, ex[0] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Y, Y, ex[1] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Z, Z, ex[2] * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ chi, chi, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ dxx, dxx, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ dyy, dyy, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ dzz, dzz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ trK, trK, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ gxy, gxy, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ gxz, gxz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ gyz, gyz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Axx, Axx, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Axy, Axy, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Axz, Axz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Ayz, Ayz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Ayy, Ayy, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Azz, Azz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Gamx, Gamx, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Gamy, Gamy, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Gamz, Gamz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ betax, betax, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ betay, betay, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ betaz, betaz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ Lap, Lap, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ dtSfx, dtSfx, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ dtSfy, dtSfy, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemcpy(Mh_ dtSfz, dtSfz, matrix_size * sizeof(double), cudaMemcpyHostToDevice);
	cudaMemset(Mh_ rho,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Sxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Sxy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Sxz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Syz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Syy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Szz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Sx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Sy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Sz,0,matrix_size * sizeof(double));

	//init local var
	cudaMemset(Mh_ gxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ chix,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ chiy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ chiz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxyx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxzx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyyx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyzx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gzzx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxxy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxyy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxzy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyyy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyzy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gzzy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxxz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxyz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gxzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyyz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gyzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gzzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Lapx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Lapy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Lapz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betaxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betaxy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betaxz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betayy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betayz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betazz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betayx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betazy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ betazx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Kx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Ky,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Kz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamxy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamxz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamyy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamyz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamyx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamzy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamzx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ div_beta,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ S,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ f,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ fxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ fxy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ fxz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ fyy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ fyz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ fzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gupxx,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gupxy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gupxz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gupyy,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gupyz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ gupzz,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamxa,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamya,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ Gamza,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ alpn1,0,matrix_size * sizeof(double));
	cudaMemset(Mh_ chin1,0,matrix_size * sizeof(double));
	
	
	double sss[3] = {1,1,1};
	double aas[3] = {-1,-1,1};
	double asa[3] = {-1,1,-1};
	double saa[3] = {1,-1,-1};
	double ass[3] = {-1,1,1};
	double sas[3] = {1,-1,1};
	double ssa[3] = {1,1,-1};

//3 --------------------Init constant memory---------------------
#if (GAUGE == 6 || GAUGE == 7)
	int BHN_h;

	double Porg_h[9];

	double Mass_h[3];
//call getpbh(BHN,Porg,Mass)
#ifdef fortran1

	getpbh

#endif	

#ifdef fortran2

	GETPBH

#endif

#ifdef fortran3

	getpbh_

#endif
	(BHN_h,double Porg_h,double Mass_h);
	
	cudaMemcpyToSymbol(BHN,&BHN_h, sizeof(int));
	cudaMemcpyToSymbol(Porg,&Porg_h,9 * sizeof(double));
	cudaMemcpyToSymbol(Mass,&Mass_h,3 * sizeof(double));
	
	double tmp_con = Mass[0] + Mass[1]; //t = M
	cudaMemcpyToSymbol(M, &tmp_con, sizeof(double));
	tmp_con = 2 / tmp_con;  //t = A
	cudaMemcpyToSymbol(A, &tmp_con, sizeof(double));
	
	double tmp_con2 = 1/Mass[0] - tmp_con;
	cudaMemcpyToSymbol(C1, &tmp_con2, sizeof(double));
	double tmp_con2 = 1/Mass[1] - tmp_con;
	cudaMemcpyToSymbol(C2, &tmp_con2, sizeof(double));
		

#endif//if (GAUGE == 6 || GAUGE == 7)
//3.1-----for compute_rhs_bssn---------
	//cout<<"Size of Meta:"<<sizeof(Meta)<<endl;
	cudaMemcpyToSymbol(metac,meta, sizeof(Meta));
	cudaMemcpyToSymbol(ex_c,ex, 3*sizeof(int));
	cudaMemcpyToSymbol(T_c,&T, sizeof(double));
	cudaMemcpyToSymbol(Symmetry_c,&Symmetry, sizeof(int));
	cudaMemcpyToSymbol(Lev_c,&Lev, sizeof(int));
	cudaMemcpyToSymbol(co_c,&co, sizeof(int));
	cudaMemcpyToSymbol(eps_c,&eps, sizeof(double));
	
	double F1o3h  = 1.0;	F1o3h /= 3.0;
	double F2o3h  = 2.0;	F2o3h /= 3.0;
	double F1o6h  = 1.0;	F1o6h /= 6.0;
	double PIh = M_PI;
	int step = GRID_DIM * BLOCK_DIM;
	double dXh = X[1] - X[0];
	double dYh = Y[1] - Y[0];
	double dZh = Z[1] - Z[0];
	
	cudaMemcpyToSymbol(F1o3,&F1o3h, sizeof(double));
	cudaMemcpyToSymbol(F2o3,&F2o3h, sizeof(double));
	cudaMemcpyToSymbol(F1o6,&F1o6h, sizeof(double));
	cudaMemcpyToSymbol(PI,&PIh, sizeof(double));
	cudaMemcpyToSymbol(STEP_SIZE,&step, sizeof(int));
	cudaMemcpyToSymbol(dX,&dXh, sizeof(double));
	cudaMemcpyToSymbol(dY,&dYh, sizeof(double));
	cudaMemcpyToSymbol(dZ,&dZh, sizeof(double));
	
	int _1d_size[4];
	int _2d_size[4];
	int _3d_size[4];
	for(int i = 0;i<4;++i){
		_1d_size[i] = ex[0] + i;
		_2d_size[i] = _1d_size[i] * (ex[1]+i);
		_3d_size[i] = _2d_size[i] * (ex[2]+i);
		//cout<<_1d_size[i]<<' '<<_2d_size[i]<<' '<<_3d_size[i]<<endl;
	}
	cudaMemcpyToSymbol(_1D_SIZE,_1d_size, 4*sizeof(int));
	cudaMemcpyToSymbol(_2D_SIZE,_2d_size, 4*sizeof(int));
	cudaMemcpyToSymbol(_3D_SIZE,_3d_size, 4*sizeof(int));

	
//3.2--------for fderivs------------
	int ijkmax_h[3] = {ex[0]-1,ex[1]-1,ex[2]-1};
	int ijkmin_h[3] = {0,0,0};
	int ijkmin2_h[3] = {0,0,0};
	int ijkmin3_h[3] = {0,0,0};
	
	double abs[3] = {X[0],Y[0],Z[0]};
	for(int i = 0;i<3;++i){
		if(abs[i] < 0) abs[i] = -abs[i]; 
	}
  	if(Symmetry > 1 && abs[0] < dXh) {ijkmin_h[0] = -2; ijkmin2_h[0] = -3;}
  	if(Symmetry > 1 && abs[1] < dYh) {ijkmin_h[1] = -2; ijkmin2_h[1] = -3;}
  	if(Symmetry > 0 && abs[2] < dZh) {ijkmin_h[2] = -2; ijkmin2_h[2] = -3;}
  	
  	if(Symmetry > 2 && abs[0] < dXh) {ijkmin3_h[0] = -3;}
  	if(Symmetry > 2 && abs[1] < dYh) {ijkmin3_h[1] = -3;}
  	if(Symmetry > 0 && abs[2] < dZh) {ijkmin3_h[2] = -3;}
  	
  	cudaMemcpyToSymbol(ijk_max,ijkmax_h,3*sizeof(int));
  	cudaMemcpyToSymbol(ijk_min,ijkmin_h,3*sizeof(int));
  	cudaMemcpyToSymbol(ijk_min2,ijkmin2_h,3*sizeof(int));
  	cudaMemcpyToSymbol(ijk_min3,ijkmin3_h,3*sizeof(int));
	
	double d12dxyz_h[3] = {1.0,1.0,1.0};
	double d2dxyz_h[3] = {1.0,1.0,1.0};
	d12dxyz_h[0] /= 12; d12dxyz_h[1] /= 12; d12dxyz_h[2] /= 12;
	d12dxyz_h[0] /= dXh; d12dxyz_h[1] /= dYh; d12dxyz_h[2] /= dZh;
	d2dxyz_h[0] /= 2; d2dxyz_h[1] /= 2; d2dxyz_h[2] /= 2;
	d2dxyz_h[0] /= dXh; d2dxyz_h[1] /= dYh; d2dxyz_h[2] /= dZh;
	
	cudaMemcpyToSymbol(d12dxyz,d12dxyz_h,3*sizeof(double));
	cudaMemcpyToSymbol(d2dxyz,d2dxyz_h,3*sizeof(double));
	
//3.3--------for fdderivs------------
	double Sdxdxh =  1.0 /( dXh * dXh ); 
	double Sdydyh =  1.0 /( dYh * dYh );
	double Sdzdzh =  1.0 /( dZh * dZh );
	double Fdxdxh = 1.0 / 12.0 /( dXh * dXh );
	double Fdydyh = 1.0 / 12.0 /( dYh * dYh );
	double Fdzdzh = 1.0 / 12.0 /( dZh * dZh );
	double Sdxdyh = 1.0/4.0 /( dXh * dYh );
	double Sdxdzh = 1.0/4.0 /( dXh * dZh );
	double Sdydzh = 1.0/4.0 /( dYh * dZh );
	double Fdxdyh = 1.0/144.0 /( dXh * dYh );
	double Fdxdzh = 1.0/144.0 /( dXh * dZh );
	double Fdydzh = 1.0/144.0 /( dYh * dZh );
	cudaMemcpyToSymbol(Sdxdx,&Sdxdxh,sizeof(double));
	cudaMemcpyToSymbol(Sdydy,&Sdydyh,sizeof(double));
	cudaMemcpyToSymbol(Sdzdz,&Sdzdzh,sizeof(double));
	cudaMemcpyToSymbol(Sdxdy,&Sdxdyh,sizeof(double));
	cudaMemcpyToSymbol(Sdxdz,&Sdxdzh,sizeof(double));
	cudaMemcpyToSymbol(Sdydz,&Sdydzh,sizeof(double));
	cudaMemcpyToSymbol(Fdxdx,&Fdxdxh,sizeof(double));
	cudaMemcpyToSymbol(Fdydy,&Fdydyh,sizeof(double));
	cudaMemcpyToSymbol(Fdzdz,&Fdzdzh,sizeof(double));
	cudaMemcpyToSymbol(Fdxdy,&Fdxdyh,sizeof(double));
	cudaMemcpyToSymbol(Fdxdz,&Fdxdzh,sizeof(double));
	cudaMemcpyToSymbol(Fdydz,&Fdydzh,sizeof(double));

//3.4---------for lopsided---------------------------


#ifdef TIMING1
	cudaThreadSynchronize();
	gettimeofday(&tv2, NULL);
   	cout<<"TIME USED"<<TimeBetween(tv1, tv2)<<endl; 
#endif	
	//cout<<"GPU meta data ready.\n";
	
	cudaThreadSynchronize();

//--------------test constant memory address & value--------------
/*	double rank = mpi_rank;
	cudaMemcpyToSymbol(F1o3,&rank, sizeof(double));
	double ctest1 = -1;
	double * ctest = &ctest1;
	double * ctest_d;
	cudaMalloc((void**)&ctest_d,sizeof(double));
	test_const_address<<<1,1>>>(ctest_d);
	cudaMemcpy(ctest, ctest_d, sizeof(double), cudaMemcpyDeviceToHost);
	cout<<"My rank is: "<<rank<<" "<<"const value is: "<<ctest[0]<<" const Address is:"<<&F1o3<<endl;
	cudaFree(ctest_d);
*/
//-------------get device info-------------------------------------
/*	int deviceCount; cudaGetDeviceCount(&deviceCount);
	cout<<"myrank is: "<<mpi_rank<<" deviceCount is:"<<deviceCount<<endl;
*/
//#4-----------------------calculate------------------------------
	//4.0------enforce_ga---------
	//sub_enforce_ga(matrix_size);
	//4.1-----compute rhs---------
	compute_rhs_bssn_part1<<<GRID_DIM,BLOCK_DIM>>>();
	cudaThreadSynchronize();

	sub_fderivs(Mh_ betax,Mh_ fh,Mh_ betaxx,Mh_ betaxy,Mh_ betaxz,ass);
	sub_fderivs(Mh_ betay,Mh_ fh,Mh_ betayx,Mh_ betayy,Mh_ betayz,sas);
	sub_fderivs(Mh_ betaz,Mh_ fh,Mh_ betazx,Mh_ betazy,Mh_ betazz,ssa);
	sub_fderivs(Mh_ chi,Mh_ fh,Mh_ chix,Mh_ chiy,Mh_ chiz, sss);
	sub_fderivs(Mh_ Lap,Mh_ fh,Mh_ Lapx,Mh_ Lapy,Mh_ Lapz, sss);
	sub_fderivs(Mh_ trK,Mh_ fh,Mh_ Kx,Mh_ Ky,Mh_ Kz, sss);
	sub_fderivs(Mh_ dxx,Mh_ fh,Mh_ gxxx,Mh_ gxxy,Mh_ gxxz, sss);
	sub_fderivs(Mh_ dyy,Mh_ fh,Mh_ gyyx,Mh_ gyyy,Mh_ gyyz, sss);
	sub_fderivs(Mh_ dzz,Mh_ fh,Mh_ gzzx,Mh_ gzzy,Mh_ gzzz, sss);
	sub_fderivs(Mh_ gxy,Mh_ fh,Mh_ gxyx,Mh_ gxyy,Mh_ gxyz, aas);
	sub_fderivs(Mh_ gxz,Mh_ fh,Mh_ gxzx,Mh_ gxzy,Mh_ gxzz, asa);
	sub_fderivs(Mh_ gyz,Mh_ fh,Mh_ gyzx,Mh_ gyzy,Mh_ gyzz, saa);
  	
  	compute_rhs_bssn_part2<<<GRID_DIM,BLOCK_DIM>>>();
	cudaThreadSynchronize();
	
	sub_fdderivs(Mh_ betax,Mh_ fh,Mh_ gxxx,Mh_ gxyx,Mh_ gxzx,Mh_ gyyx,Mh_ gyzx,Mh_ gzzx,ass);
	sub_fdderivs(Mh_ betay,Mh_ fh,Mh_ gxxy,Mh_ gxyy,Mh_ gxzy,Mh_ gyyy,Mh_ gyzy,Mh_ gzzy,sas);
	sub_fdderivs(Mh_ betaz,Mh_ fh,Mh_ gxxz,Mh_ gxyz,Mh_ gxzz,Mh_ gyyz,Mh_ gyzz,Mh_ gzzz,ssa);        
	sub_fderivs( Mh_ Gamx, Mh_ fh,Mh_ Gamxx, Mh_ Gamxy, Mh_ Gamxz,ass);
	sub_fderivs( Mh_ Gamy, Mh_ fh,Mh_ Gamyx, Mh_ Gamyy, Mh_ Gamyz,sas);
	sub_fderivs( Mh_ Gamz, Mh_ fh,Mh_ Gamzx, Mh_ Gamzy, Mh_ Gamzz,ssa);
	
	compute_rhs_bssn_part3<<<GRID_DIM,BLOCK_DIM>>>();
	cudaThreadSynchronize();
	
	computeRicci(Mh_ dxx,Mh_ Rxx,sss, meta);
	computeRicci(Mh_ dyy,Mh_ Ryy,sss, meta);
	computeRicci(Mh_ dzz,Mh_ Rzz,sss, meta);
	computeRicci(Mh_ gxy,Mh_ Rxy,aas, meta);
	computeRicci(Mh_ gxz,Mh_ Rxz,asa, meta);
	computeRicci(Mh_ gyz,Mh_ Ryz,saa, meta);
	
	cudaThreadSynchronize();
	
	compute_rhs_bssn_part4<<<GRID_DIM,BLOCK_DIM>>>();
	cudaThreadSynchronize();
	
	sub_fdderivs(Mh_ chi,Mh_ fh,Mh_ fxx,Mh_ fxy,Mh_ fxz,Mh_ fyy,Mh_ fyz,Mh_ fzz,sss);
	
	compute_rhs_bssn_part5<<<GRID_DIM,BLOCK_DIM>>>();
	cudaThreadSynchronize();
	
	sub_fdderivs(Mh_ Lap,Mh_ fh,Mh_ fxx,Mh_ fxy,Mh_ fxz,Mh_ fyy,Mh_ fyz,Mh_ fzz,sss);
	
	compute_rhs_bssn_part6<<<GRID_DIM,BLOCK_DIM>>>();
	cudaThreadSynchronize();
	
#if (GAUGE == 2 || GAUGE == 3 || GAUGE == 4 || GAUGE == 5)
	sub_fderivs(Mh_ chi,Mh_ fh, Mh_ dtSfx_rhs, Mh_ dtSfy_rhs, Mh_ dtSfz_rhs,sss);
	compute_rhs_bssn_part6_gauge<<<GRID_DIM,BLOCK_DIM>>>();
#endif	

	sub_lopsided(Mh_ gxx,Mh_ fh2,Mh_ gxx_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ gxy,Mh_ fh2,Mh_ gxy_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,aas);
	sub_lopsided(Mh_ gxz,Mh_ fh2,Mh_ gxz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,asa);
	sub_lopsided(Mh_ gyy,Mh_ fh2,Mh_ gyy_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ gyz,Mh_ fh2,Mh_ gyz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,saa);
	sub_lopsided(Mh_ gzz,Mh_ fh2,Mh_ gzz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ Axx,Mh_ fh2,Mh_ Axx_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ Axy,Mh_ fh2,Mh_ Axy_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,aas);
	sub_lopsided(Mh_ Axz,Mh_ fh2,Mh_ Axz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,asa);
	sub_lopsided(Mh_ Ayy,Mh_ fh2,Mh_ Ayy_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ Ayz,Mh_ fh2,Mh_ Ayz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,saa);
	sub_lopsided(Mh_ Azz,Mh_ fh2,Mh_ Azz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ chi,Mh_ fh2,Mh_ chi_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ trK,Mh_ fh2,Mh_ trK_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	sub_lopsided(Mh_ Gamx,Mh_ fh2,Mh_ Gamx_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,ass);
	sub_lopsided(Mh_ Gamy,Mh_ fh2,Mh_ Gamy_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sas);
	sub_lopsided(Mh_ Gamz,Mh_ fh2,Mh_ Gamz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,ssa);
	sub_lopsided(Mh_ Lap,Mh_ fh2,Mh_ Lap_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sss);
	
#if (GAUGE == 0 || GAUGE == 1 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)

  	sub_lopsided(Mh_ betax,Mh_ fh2,Mh_ betax_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,ass);
	sub_lopsided(Mh_ betay,Mh_ fh2,Mh_ betay_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sas);
	sub_lopsided(Mh_ betaz,Mh_ fh2,Mh_ betaz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,ssa);

#endif
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)	
	sub_lopsided(Mh_ dtSfx,Mh_ fh2,Mh_ dtSfx_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,ass);
	sub_lopsided(Mh_ dtSfy,Mh_ fh2,Mh_ dtSfy_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,sas);
	sub_lopsided(Mh_ dtSfz,Mh_ fh2,Mh_ dtSfz_rhs,Mh_ betax,Mh_ betay,Mh_ betaz,ssa);
#endif	
	if(eps > 0){
		sub_kodis(Mh_ chi,Mh_ fh2, Mh_ chi_rhs,sss);
		sub_kodis(Mh_ trK,Mh_ fh2, Mh_ trK_rhs,sss);
		sub_kodis(Mh_ dxx,Mh_ fh2, Mh_ gxx_rhs,sss);
		sub_kodis(Mh_ gxy,Mh_ fh2, Mh_ gxy_rhs,aas);
		sub_kodis(Mh_ gxz,Mh_ fh2, Mh_ gxz_rhs,asa);
		sub_kodis(Mh_ dyy,Mh_ fh2, Mh_ gyy_rhs,sss);
		sub_kodis(Mh_ gyz,Mh_ fh2, Mh_ gyz_rhs,saa);
		sub_kodis(Mh_ dzz,Mh_ fh2, Mh_ gzz_rhs,sss);
		sub_kodis(Mh_ Axx,Mh_ fh2, Mh_ Axx_rhs,sss);
		sub_kodis(Mh_ Axy,Mh_ fh2, Mh_ Axy_rhs,aas);
		sub_kodis(Mh_ Axz,Mh_ fh2, Mh_ Axz_rhs,asa);
		sub_kodis(Mh_ Ayy,Mh_ fh2, Mh_ Ayy_rhs,sss);
		sub_kodis(Mh_ Ayz,Mh_ fh2, Mh_ Ayz_rhs,saa);
		sub_kodis(Mh_ Azz,Mh_ fh2, Mh_ Azz_rhs,sss);
		sub_kodis(Mh_ Gamx,Mh_ fh2, Mh_ Gamx_rhs,ass);
		sub_kodis(Mh_ Gamy,Mh_ fh2, Mh_ Gamy_rhs,sas);
		sub_kodis(Mh_ Gamz,Mh_ fh2, Mh_ Gamz_rhs,ssa);
		
		sub_kodis(Mh_ Lap,Mh_ fh2, Mh_ Lap_rhs,sss);
		sub_kodis(Mh_ betax,Mh_ fh2, Mh_ betax_rhs,ass);
		sub_kodis(Mh_ betay,Mh_ fh2, Mh_ betay_rhs,sas);
		sub_kodis(Mh_ betaz,Mh_ fh2, Mh_ betaz_rhs,ssa);
		
#if (GAUGE == 0 || GAUGE == 2 || GAUGE == 3 || GAUGE == 6 || GAUGE == 7)		
		sub_kodis(Mh_ dtSfx,Mh_ fh2, Mh_ dtSfx_rhs,ass);
		sub_kodis(Mh_ dtSfy,Mh_ fh2, Mh_ dtSfy_rhs,sas);
		sub_kodis(Mh_ dtSfz,Mh_ fh2, Mh_ dtSfz_rhs,ssa);
#endif

	}
	
	if(co == 0){
		compute_rhs_bssn_part7<<<GRID_DIM,BLOCK_DIM>>>();
		cudaThreadSynchronize();

		sub_fderivs(Mh_ Axx,Mh_ fh,Mh_ gxxx,Mh_ gxxy,Mh_ gxxz,sss);
		sub_fderivs(Mh_ Axy,Mh_ fh,Mh_ gxyx,Mh_ gxyy,Mh_ gxyz,aas);
		sub_fderivs(Mh_ Axz,Mh_ fh,Mh_ gxzx,Mh_ gxzy,Mh_ gxzz,asa);
		sub_fderivs(Mh_ Ayy,Mh_ fh,Mh_ gyyx,Mh_ gyyy,Mh_ gyyz,sss);
		sub_fderivs(Mh_ Ayz,Mh_ fh,Mh_ gyzx,Mh_ gyzy,Mh_ gyzz,saa);
		sub_fderivs(Mh_ Azz,Mh_ fh,Mh_ gzzx,Mh_ gzzy,Mh_ gzzz,sss);
		compute_rhs_bssn_part8<<<GRID_DIM,BLOCK_DIM>>>();
		cudaThreadSynchronize();
	}

#if (ABV == 1)
	cout<<"TODO: bssn_gpu.cu::2373 (ABV == 1)"<<endl;
#endif
//5---------------------------get result----------------------------
	/*cudaMemcpy(chi, Mh_ chi, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(dxx, Mh_ dxx, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(dyy, Mh_ dyy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(dzz, Mh_ dzz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(Lap, Mh_ Lap, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(betax, Mh_ betax, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(betay, Mh_ betay, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	cudaMemcpy(betaz, Mh_ betaz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);*/
	if(calledby == CALLED_BY_STEP)
	{	
		cudaMemcpy(chi_rhs, Mh_ chi_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(trK_rhs, Mh_ trK_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(gxx_rhs, Mh_ gxx_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(gxy_rhs, Mh_ gxy_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(gxz_rhs, Mh_ gxz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(gyy_rhs, Mh_ gyy_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(gyz_rhs, Mh_ gyz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(gzz_rhs, Mh_ gzz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Axx_rhs, Mh_ Axx_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Axy_rhs, Mh_ Axy_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Axz_rhs, Mh_ Axz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ayy_rhs, Mh_ Ayy_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ayz_rhs, Mh_ Ayz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Azz_rhs, Mh_ Azz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamx_rhs, Mh_ Gamx_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamy_rhs, Mh_ Gamy_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamz_rhs, Mh_ Gamz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Lap_rhs, Mh_ Lap_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(betax_rhs, Mh_ betax_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(betay_rhs, Mh_ betay_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(betaz_rhs, Mh_ betaz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(dtSfx_rhs, Mh_ dtSfx_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(dtSfy_rhs, Mh_ dtSfy_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(dtSfz_rhs, Mh_ dtSfz_rhs, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	}
	else if(calledby == CALLED_BY_CONSTRAINT)
	{
		cudaMemcpy(Gamxxx, Mh_ Gamxxx, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamxxy, Mh_ Gamxxy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamxxz, Mh_ Gamxxz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamxyy, Mh_ Gamxyy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamxyz, Mh_ Gamxyz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamxzz, Mh_ Gamxzz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamyxx, Mh_ Gamyxx, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamyxy, Mh_ Gamyxy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamyxz, Mh_ Gamyxz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamyyy, Mh_ Gamyyy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamyyz, Mh_ Gamyyz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamyzz, Mh_ Gamyzz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamzxx, Mh_ Gamzxx, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamzxy, Mh_ Gamzxy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamzxz, Mh_ Gamzxz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamzyy, Mh_ Gamzyy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamzyz, Mh_ Gamzyz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gamzzz, Mh_ Gamzzz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Rxx, Mh_ Rxx, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Rxy, Mh_ Rxy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Rxz, Mh_ Rxz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ryy, Mh_ Ryy, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Ryz, Mh_ Ryz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Rzz, Mh_ Rzz, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(ham_Res, Mh_ ham_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(movx_Res, Mh_ movx_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(movy_Res, Mh_ movy_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(movz_Res, Mh_ movz_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gmx_Res, Mh_ Gmx_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gmy_Res, Mh_ Gmy_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
		cudaMemcpy(Gmz_Res, Mh_ Gmz_Res, matrix_size * sizeof(double), cudaMemcpyDeviceToHost);
	}

//-----------------------------------------------------
//-------------------FOR GPU TEST----------------------
//-----------------------------------------------------
#ifdef TIMING
	cudaThreadSynchronize();
	gettimeofday(&tv2, NULL);
   	cout<<"MPI rank is: "<<mpi_rank<<" GPU TIME is"<<TimeBetween(tv1, tv2)<<" (s)."<<endl; 
#endif


	destroy_meta(meta);

	
	return 0;//TODO return
}
