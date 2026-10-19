#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
#include <cassert>
using namespace std;

#include "Block.h"
#include "misc.h"

#ifdef USE_GPU
#include "gpu_manager.h"
#include "gpu_info.h"
#endif
#include "fmisc.h"

#ifdef USE_GPU
constexpr int pad_x_left = GPU_PAD_X_LEFT;
#endif

Block::Block(int DIM, int *shapei, double *bboxi, int ranki, int ingfsi, int fngfsi, int levi, const int cgpui) : rank(ranki), ingfs(ingfsi), fngfs(fngfsi), lev(levi), cgpu(cgpui)
{
	for (int i = 0; i < dim; i++)
		X[i] = 0;

	if (DIM != dim)
	{
		cout << "dimension is not consistent in Block construction" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	bool flag = false;
	for (int i = 0; i < dim; i ++) {
		shape[i] = shapei[i];
		if (shape[i] <= 0)
			flag = true;
		bbox[i] = bboxi[i];
		bbox[dim + i] = bboxi[dim + i];
	}
#ifdef USE_GPU
	int raw_shape_pad_x = pad_x_left + shape[0] + ghost_width /* + GPU_DIM_X - 1 */;
	shape_pad[0] = (raw_shape_pad_x + GPU_ALIGNMENT_ASYNC_DOUBLE - 1) / GPU_ALIGNMENT_ASYNC_DOUBLE * GPU_ALIGNMENT_ASYNC_DOUBLE;
	shape_pad[1] = shape[1] + 2 * ghost_width + GPU_DIM_Y - 1;
	shape_pad[2] = shape[2] + 2 * ghost_width + GPU_DIM_Z - 1;
	// cout << "shape_pad: (" << shape_pad[0] << "," << shape_pad[1] << "," << shape_pad[2] << ")" << endl;
#endif

	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (flag)
	{
		cout << "myrank: " << myrank << ", on rank: " << rank << endl;
		cout << "error shape in Block construction: (" << shape[0] << "," << shape[1] << "," << shape[2] << ")" << endl;
		cout << "box boundary: (" << bbox[0] << ":" << bbox[3] << "," << bbox[1] << ":" << bbox[4] << "," << bbox[2] << ":" << bbox[5] << ")" << endl;
		cout << "belong to level " << lev << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	if (myrank == rank)
	{
		for (int i = 0; i < dim; i++)
		{
				X[i] = new double[shape[i]];
#ifdef Vertex
#ifdef Cell
#error "Both Cell and Vertex are defined"
#endif
				double h = (bbox[dim + i] - bbox[i]) / (shape[i] - 1);
				for (int j = 0; j < shape[i]; j++)
					X[i][j] = bbox[i] + j * h;
#else
#ifdef Cell
				double h = (bbox[dim + i] - bbox[i]) / shape[i];
				for (int j = 0; j < shape[i]; j++)
					X[i][j] = bbox[i] + (j + 0.5) * h;
#else
#error "Not define Vertex nor Cell"
#endif
#endif
#ifdef USE_GPU
            d_X[i] = GPUManager::getInstance().allocate_device_memory<double>(shape[i], 0, false);
            GPUManager::getInstance().sync_to_gpu(X[i], d_X[i], shape[i]);
#endif
		}

		int nn = shape[0] * shape[1] * shape[2];
		fgfs = new double *[fngfs];
#ifdef USE_GPU
		int nn_pad = shape_pad[0] * shape_pad[1] * shape_pad[2] + GPU_DIM_X - 1;

		// initialize GPU shadow pointers and valid flags
		d_fgfs = new double *[fngfs];
		d_fgfs_base = new double *[fngfs];
		cpu_valid = new bool[fngfs];
		gpu_valid = new bool[fngfs];

		int offset = (ghost_width * shape_pad[1] + ghost_width) * shape_pad[0] + pad_x_left;
		d_fgfs_transfer_buffer = GPUManager::getInstance().allocate_device_memory<double>(nn, 0, false);
#endif

		for (int i = 0; i < fngfs; i++) {
			fgfs[i] = (double *)malloc(sizeof(double) * nn);
			if (!(fgfs[i]))
			{
				cout << "on node#" << rank << ", out of memory when constructing Block." << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			// memset(fgfs[i], 0, sizeof(double) * nn); // ?
#ifdef USE_GPU
            d_fgfs_base[i] = GPUManager::getInstance().allocate_device_memory<double>(nn_pad, 0, false);
			d_fgfs[i] = d_fgfs_base[i] + offset; // point to the start of the actual data region, skipping ghost zones
			// assert(reinterpret_cast<std::uintptr_t>(d_fgfs[i]) % 32 == 0);
			// assert(reinterpret_cast<std::uintptr_t>(d_fgfs_base[i]) % 128 == 0);
            cpu_valid[i] = true;
            gpu_valid[i] = false;
#endif
		}

		igfs = new int *[ingfs];
		for (int i = 0; i < ingfs; i++)
		{
			igfs[i] = (int *)malloc(sizeof(int) * nn);
			if (!(igfs[i]))
			{
				cout << "on node#" << rank << ", out of memory when constructing Block." << endl;
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
			// memset(igfs[i], 0, sizeof(int) * nn);
		}
	}

#ifdef USE_GPU
	stream = GPUManager::getInstance().get_stream();
	GPUManager::getInstance().synchronize_memory();
#endif
}
Block::~Block()
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == rank) {
		int nn = shape[0] * shape[1] * shape[2];
#ifdef USE_GPU
		int nn_pad = shape_pad[0] * shape_pad[1] * shape_pad[2];
#endif
		for (int i = 0; i < dim; i ++) {
#ifdef USE_GPU
			GPUManager::getInstance().free_device_memory(d_X[i], shape[i], 0, false);
#endif
			delete[] X[i];
		}
		for (int i = 0; i < ingfs; i++)
			free(igfs[i]);
		delete[] igfs;
		for (int i = 0; i < fngfs; i ++) {
			free(fgfs[i]);
#ifdef USE_GPU
			GPUManager::getInstance().free_device_memory(d_fgfs_base[i], nn_pad, 0, false);
#endif
		}
#ifdef USE_GPU
		GPUManager::getInstance().free_device_memory(d_fgfs_transfer_buffer, nn, 0, false);
#endif
		delete[] fgfs;
#ifdef USE_GPU
		delete[] d_fgfs;
		delete[] d_fgfs_base;
		delete[] cpu_valid;
		delete[] gpu_valid;
#endif
		X[0] = X[1] = X[2] = 0;
		igfs = 0;
		fgfs = 0;
#ifdef USE_GPU
		d_fgfs = 0;
		d_fgfs_base = 0;
		cpu_valid = 0;
		gpu_valid = 0;
#endif
	}
#ifdef USE_GPU
	GPUManager::getInstance().synchronize_memory();
#endif
}
void Block::checkBlock()
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0)
	{
		cout << "belong to level " << lev << endl;
		cout << "shape: [";
		for (int i = 0; i < dim; i++)
		{
			cout << shape[i];
			if (i < dim - 1)
				cout << ",";
			else
				cout << "]";
		}
		cout << " resolution: [";
		for (int i = 0; i < dim; i++)
		{
			cout << getdX(i);
			if (i < dim - 1)
				cout << ",";
			else
				cout << "]" << endl;
		}
		cout << "locate on node " << rank << ", at (includes ghost zone):" << endl;
		cout << "(";
		for (int i = 0; i < dim; i++)
		{
			cout << bbox[i] << ":" << bbox[dim + i];
			if (i < dim - 1)
				cout << ",";
			else
				cout << ")" << endl;
		}
		cout << "has " << ingfs << " int type grids functions," << fngfs << " double type grids functions" << endl;
	}
}
double Block::getdX(int dir)
{
	if (dir < 0 || dir >= dim)
	{
		cout << "Block::getdX: error input dir = " << dir << ", this Block has direction (0," << dim - 1 << ")" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	double h;
#ifdef Vertex
#ifdef Cell
#error "Both Cell and Vertex are defined"
#endif
	if (shape[dir] == 1)
	{
		cout << "Block::getdX: for direction " << dir << ", this Block has only one point. Can not determine dX for vertex center grid." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	h = (bbox[dim + dir] - bbox[dir]) / (shape[dir] - 1);
#else
#ifdef Cell
	h = (bbox[dim + dir] - bbox[dir]) / shape[dir];
#else
#error "Not define Vertex nor Cell"
#endif
#endif
	return h;
}
void Block::swapList(MyList<var> *VarList1, MyList<var> *VarList2, int myrank)
{
	if (rank == myrank)
	{
		MyList<var> *varl1 = VarList1, *varl2 = VarList2;
		while (varl1 && varl2)
		{
			misc::swap<double *>(fgfs[varl1->data->sgfn], fgfs[varl2->data->sgfn]);
#ifdef USE_GPU
			misc::swap<double *>(d_fgfs[varl1->data->sgfn], d_fgfs[varl2->data->sgfn]);
			misc::swap<double *>(d_fgfs_base[varl1->data->sgfn], d_fgfs_base[varl2->data->sgfn]);
			misc::swap<bool>(cpu_valid[varl1->data->sgfn], cpu_valid[varl2->data->sgfn]);
			misc::swap<bool>(gpu_valid[varl1->data->sgfn], gpu_valid[varl2->data->sgfn]);
#endif
			varl1 = varl1->next;
			varl2 = varl2->next;
		}
		if (varl1 || varl2)
		{
			cout << "error in Block::swaplist, var lists does not match." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
	}
}


#ifdef USE_GPU
void Block::require_on_gpu(int var_index) {
    if (!gpu_valid[var_index] && cpu_valid[var_index]) {
        int nn = shape[0] * shape[1] * shape[2];
        GPUManager::getInstance().sync_to_gpu(fgfs[var_index], d_fgfs_transfer_buffer, nn);
		GPUManager::getInstance().unpack_to_pad(d_fgfs_transfer_buffer, d_fgfs[var_index], shape, shape_pad);
        gpu_valid[var_index] = true; 
    }
}

void Block::require_on_cpu(int var_index) {
    if (!cpu_valid[var_index] && gpu_valid[var_index]) {
        int nn = shape[0] * shape[1] * shape[2];
        GPUManager::getInstance().pack_to_pad(d_fgfs[var_index], d_fgfs_transfer_buffer, shape, shape_pad);
        GPUManager::getInstance().sync_to_cpu(fgfs[var_index], d_fgfs_transfer_buffer, nn);
        cpu_valid[var_index] = true; 
    }
}

void Block::mark_gpu_modified(int var_index) {
    gpu_valid[var_index] = true;
    cpu_valid[var_index] = false; 
}

void Block::mark_cpu_modified(int var_index) {
    cpu_valid[var_index] = true;
    gpu_valid[var_index] = false; 
}

void Block::move_to_gpu(MyList<var> *VarList) {
	MyList<var> *iter = VarList;
	while (iter) {
		int var_index = iter->data->sgfn;
		mark_cpu_modified(var_index);
		require_on_gpu(var_index);
		iter = iter->next;
	}
}

void Block::move_to_cpu(MyList<var> *VarList) {
	MyList<var> *iter = VarList;
	while (iter) {
		int var_index = iter->data->sgfn;
		mark_gpu_modified(var_index);
		require_on_cpu(var_index);
		iter = iter->next;
	}
}
#endif
