#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "Block.h"
#include "misc.h"

#include "gpu_runtime_coordinator.h"

namespace {

int current_rank_id() {
	int myrank = 0;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	return myrank;
}

int block_point_count(const Block& block) {
	return block.shape[0] * block.shape[1] * block.shape[2];
}

void abort_block_error(const char* message) {
	cout << message << endl;
	MPI_Abort(MPI_COMM_WORLD, 1);
}

void initialize_coordinate_axis(Block& block, int axis) {
	block.X[axis] = new double[block.shape[axis]];
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
	if (block.shape[axis] == 1)
		abort_block_error("Block coordinate initialization: Vertex grid needs at least two points");
	const double spacing = (block.bbox[dim + axis] - block.bbox[axis]) / (block.shape[axis] - 1);
	for (int point_index = 0; point_index < block.shape[axis]; ++point_index)
		block.X[axis][point_index] = block.bbox[axis] + point_index * spacing;
#else
#ifdef Cell
	const double spacing = (block.bbox[dim + axis] - block.bbox[axis]) / block.shape[axis];
	for (int point_index = 0; point_index < block.shape[axis]; ++point_index)
		block.X[axis][point_index] = block.bbox[axis] + (point_index + 0.5) * spacing;
#else
#error Not define Vertex nor Cell
#endif
#endif
	block.d_X[axis] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(block.shape[axis]);
	GpuRuntimeCoordinator::instance().copy_to_device(block.X[axis], block.d_X[axis], block.shape[axis]);
}

void initialize_field_storage(Block& block, int point_count) {
	block.fgfs = new double *[block.fngfs];
	block.d_fgfs = new double *[block.fngfs];
	block.cpu_valid = new bool[block.fngfs];
	block.gpu_valid = new bool[block.fngfs];

	for (int field_index = 0; field_index < block.fngfs; ++field_index) {
		block.fgfs[field_index] = (double *)malloc(sizeof(double) * point_count);
		if (!(block.fgfs[field_index])) {
			cout << "on node#" << block.rank << ", out of memory when constructing Block." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		memset(block.fgfs[field_index], 0, sizeof(double) * point_count);
		block.d_fgfs[field_index] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(point_count);
		block.cpu_valid[field_index] = true;
		block.gpu_valid[field_index] = false;
	}
}

void initialize_index_storage(Block& block, int point_count) {
	block.igfs = new int *[block.ingfs];
	for (int field_index = 0; field_index < block.ingfs; ++field_index) {
		block.igfs[field_index] = (int *)malloc(sizeof(int) * point_count);
		if (!(block.igfs[field_index])) {
			cout << "on node#" << block.rank << ", out of memory when constructing Block." << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		memset(block.igfs[field_index], 0, sizeof(int) * point_count);
	}
}

void release_coordinate_storage(Block& block) {
	for (int axis = 0; axis < dim; ++axis) {
		GpuRuntimeCoordinator::instance().release_device_buffer(block.d_X[axis], block.shape[axis]);
		delete[] block.X[axis];
	}
	block.X[0] = block.X[1] = block.X[2] = 0;
}

void release_index_storage(Block& block) {
	for (int field_index = 0; field_index < block.ingfs; ++field_index)
		free(block.igfs[field_index]);
	delete[] block.igfs;
	block.igfs = 0;
}

void release_field_storage(Block& block, int point_count) {
	for (int field_index = 0; field_index < block.fngfs; ++field_index) {
		free(block.fgfs[field_index]);
		GpuRuntimeCoordinator::instance().release_device_buffer(block.d_fgfs[field_index], point_count);
	}
	delete[] block.fgfs;
	delete[] block.d_fgfs;
	delete[] block.cpu_valid;
	delete[] block.gpu_valid;
	block.fgfs = 0;
	block.d_fgfs = 0;
}

void mark_storage_state(bool& cpu_valid, bool& gpu_valid, bool mark_cpu) {
	if (mark_cpu) {
		cpu_valid = true;
		gpu_valid = false;
	} else {
		gpu_valid = true;
		cpu_valid = false;
	}
}

void transfer_var_list(Block& block, MyList<var> *var_list, bool to_gpu) {
	MyList<var> *iter = var_list;
	while (iter) {
		const int var_index = iter->data->sgfn;
		if (to_gpu) {
			block.mark_cpu_modified(var_index);
			block.require_on_gpu(var_index);
		} else {
			block.mark_gpu_modified(var_index);
			block.require_on_cpu(var_index);
		}
		iter = iter->next;
	}
}

} 

Block::Block(int DIM, int *shapei, double *bboxi, int ranki, int ingfsi, int fngfsi, int levi, const int cgpui) : rank(ranki), ingfs(ingfsi), fngfs(fngfsi), lev(levi), cgpu(cgpui)
{
	for (int i = 0; i < dim; i++)
		X[i] = 0;

	if (DIM != dim)
		abort_block_error("dimension is not consistent in Block construction");

	bool flag = false;
	for (int i = 0; i < dim; i++)
	{
		shape[i] = shapei[i];
		if (shape[i] <= 0)
			flag = true;
		bbox[i] = bboxi[i];
		bbox[dim + i] = bboxi[dim + i];
	}

	const int myrank = current_rank_id();
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
		for (int axis = 0; axis < dim; ++axis)
			initialize_coordinate_axis(*this, axis);

		const int point_count = block_point_count(*this);
		initialize_field_storage(*this, point_count);
		initialize_index_storage(*this, point_count);
	}

	stream = GpuRuntimeCoordinator::instance().acquire_stream();
	GpuRuntimeCoordinator::instance().wait_for_memory_ops();
}
Block::~Block()
{
	const int myrank = current_rank_id();
	if (myrank == rank) {
		const int point_count = block_point_count(*this);
		release_coordinate_storage(*this);
		release_index_storage(*this);
		release_field_storage(*this, point_count);
	}
	GpuRuntimeCoordinator::instance().wait_for_memory_ops();
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
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
	if (shape[dir] == 1)
	{
		cout << "Block::getdX: for direction " << dir << ", this Block has only one point. Can not determine dX for vertex center grid." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return (bbox[dim + dir] - bbox[dir]) / (shape[dir] - 1);
#else
#ifdef Cell
	return (bbox[dim + dir] - bbox[dir]) / shape[dir];
#else
#error Not define Vertex nor Cell
#endif
#endif
}
void Block::swapList(MyList<var> *VarList1, MyList<var> *VarList2, int myrank)
{
	if (rank == myrank)
	{
		MyList<var> *varl1 = VarList1, *varl2 = VarList2;
		while (varl1 && varl2)
		{
			misc::swap<double *>(fgfs[varl1->data->sgfn], fgfs[varl2->data->sgfn]);
			misc::swap<double *>(d_fgfs[varl1->data->sgfn], d_fgfs[varl2->data->sgfn]);
			misc::swap<bool>(cpu_valid[varl1->data->sgfn], cpu_valid[varl2->data->sgfn]);
			misc::swap<bool>(gpu_valid[varl1->data->sgfn], gpu_valid[varl2->data->sgfn]);
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


void Block::require_on_gpu(int var_index) {
    if (!gpu_valid[var_index] && cpu_valid[var_index]) {
		GpuRuntimeCoordinator::instance().copy_to_device(fgfs[var_index], d_fgfs[var_index], block_point_count(*this), stream);
		gpu_valid[var_index] = true;
	}
}

void Block::require_on_cpu(int var_index) {
	if (!cpu_valid[var_index] && gpu_valid[var_index]) {
		GpuRuntimeCoordinator::instance().copy_to_host(fgfs[var_index], d_fgfs[var_index], block_point_count(*this), stream);
		cpu_valid[var_index] = true;
	}
}

void Block::mark_gpu_modified(int var_index) {
	mark_storage_state(cpu_valid[var_index], gpu_valid[var_index], false);
}

void Block::mark_cpu_modified(int var_index) {
	mark_storage_state(cpu_valid[var_index], gpu_valid[var_index], true);
}

void Block::move_to_gpu(MyList<var> *VarList) {
	transfer_var_list(*this, VarList, true);
}

void Block::move_to_cpu(MyList<var> *VarList) {
	transfer_var_list(*this, VarList, false);
}
