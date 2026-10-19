
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <new>
using namespace std;

#include "misc.h"
#include "MPatch.h"
#include "Parallel.h"
#include "fmisc.h"
#include "gpu_runtime_coordinator.h"

namespace {

double patch_spacing(const double* bbox, const int* shape, int axis, const char* owner) {
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
	if (shape[axis] == 1)
	{
		cout << owner << ": for direction " << axis << ", this Patch has only one point. Can not determine dX for vertex center grid." << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	return (bbox[dim + axis] - bbox[axis]) / (shape[axis] - 1);
#else
#ifdef Cell
	return (bbox[dim + axis] - bbox[axis]) / shape[axis];
#else
#error Not define Vertex nor Cell
#endif
#endif
}

void set_block_interp_bounds(const Block* block, const Patch* patch, const double* DH, double* llb, double* uub) {
	for (int i = 0; i < dim; i++)
	{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
		llb[i] = (feq(block->bbox[i], patch->bbox[i], DH[i] / 2)) ? block->bbox[i] + patch->lli[i] * DH[i] : block->bbox[i] + (ghost_width - 0.5) * DH[i];
		uub[i] = (feq(block->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2)) ? block->bbox[dim + i] - patch->uui[i] * DH[i] : block->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
		llb[i] = (feq(block->bbox[i], patch->bbox[i], DH[i] / 2)) ? block->bbox[i] + patch->lli[i] * DH[i] : block->bbox[i] + ghost_width * DH[i];
		uub[i] = (feq(block->bbox[dim + i], patch->bbox[dim + i], DH[i] / 2)) ? block->bbox[dim + i] - patch->uui[i] * DH[i] : block->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
	}
}

}

Patch::Patch(int DIM, int *shapei, double *bboxi, int levi, bool buflog, int Symmetry)
		: lev(levi),
			d_interp_local_shellf(0),
			d_interp_local_weight(0),
			d_interp_global_weight(0),
			d_interp_err_flag(0),
			h_interp_shellf(0),
			h_interp_weight(0),
			interp_shellf_capacity(0),
			interp_weight_capacity(0),
			interp_host_shellf_capacity(0),
			interp_host_weight_capacity(0)
{

	int hbuffer_width = buffer_width;
	if (lev == 0)
		hbuffer_width = CS_width; // specific for shell-box coulping

	if (DIM != dim)
	{
		cout << "dimension is not consistent in Patch construction" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	for (int i = 0; i < dim; i++)
	{
		shape[i] = shapei[i];
		bbox[i] = bboxi[i];
		bbox[dim + i] = bboxi[dim + i];
		lli[i] = uui[i] = 0;
		if (buflog)
		{
			double DH;
			DH = patch_spacing(bbox, shape, i, "Patch construction");
			uui[i] = hbuffer_width;
			bbox[dim + i] = bbox[dim + i] + uui[i] * DH;
			shape[i] = shape[i] + uui[i];
		}
	}

	if (buflog)
	{
		if (DIM != 3)
		{
			cout << "Symmetry in Patch construction only support 3 yet but dim = " << DIM << endl;
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		double tmpb, DH;
		if (Symmetry > 0)
		{
			DH = patch_spacing(bbox, shape, 2, "Patch construction");
			tmpb = Mymax(0, bbox[2] - hbuffer_width * DH);
			lli[2] = int((bbox[2] - tmpb) / DH + 0.4);
			bbox[2] = bbox[2] - lli[2] * DH;
			shape[2] = shape[2] + lli[2];
			if (lli[2] < hbuffer_width)
			{
				if (feq(bbox[2], 0, DH / 2))
					lli[2] = 0;
				else
				{
					cout << "Code mistake for lli[2] = " << lli[2] << ", bbox[2] = " << bbox[2] << endl;
					MPI_Abort(MPI_COMM_WORLD, 1);
				}
			}
			if (Symmetry > 1)
			{
				DH = patch_spacing(bbox, shape, 0, "Patch construction");
				tmpb = Mymax(0, bbox[0] - hbuffer_width * DH);
				lli[0] = int((bbox[0] - tmpb) / DH + 0.4);
				bbox[0] = bbox[0] - lli[0] * DH;
				shape[0] = shape[0] + lli[0];
				if (lli[0] < hbuffer_width)
				{
					if (feq(bbox[0], 0, DH / 2))
						lli[0] = 0;
					else
					{
						cout << "Code mistake for lli[0] = " << lli[0] << ", bbox[0] = " << bbox[0] << endl;
						MPI_Abort(MPI_COMM_WORLD, 1);
					}
				}
				DH = patch_spacing(bbox, shape, 1, "Patch construction");
				tmpb = Mymax(0, bbox[1] - hbuffer_width * DH);
				lli[1] = int((bbox[1] - tmpb) / DH + 0.4);
				bbox[1] = bbox[1] - lli[1] * DH;
				shape[1] = shape[1] + lli[1];
				if (lli[1] < hbuffer_width)
				{
					if (feq(bbox[1], 0, DH / 2))
						lli[1] = 0;
					else
					{
						cout << "Code mistake for lli[1] = " << lli[1] << ", bbox[1] = " << bbox[1] << endl;
						MPI_Abort(MPI_COMM_WORLD, 1);
					}
				}
			}
			else
			{
				for (int i = 0; i < 2; i++)
				{
					DH = patch_spacing(bbox, shape, i, "Patch construction");
					lli[i] = hbuffer_width;
					bbox[i] = bbox[i] - lli[i] * DH;
					shape[i] = shape[i] + lli[i];
				}
			}
		}
		else
		{
			for (int i = 0; i < dim; i++)
			{
				DH = patch_spacing(bbox, shape, i, "Patch construction");
				lli[i] = hbuffer_width;
				bbox[i] = bbox[i] - lli[i] * DH;
				shape[i] = shape[i] + lli[i];
			}
		}
	}

	blb = ble = 0;
}
Patch::~Patch()
{
	GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_local_shellf, interp_shellf_capacity);
	GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_local_weight, interp_weight_capacity);
	GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_global_weight, interp_weight_capacity);
	GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_err_flag, 1);
	delete[] h_interp_shellf;
	delete[] h_interp_weight;
}

void Patch::ensure_interp_device_cache(int point_count, int variable_count, cudaStream_t stream)
{
	const int required_shellf = point_count * variable_count;
	if (required_shellf > interp_shellf_capacity)
	{
		GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_local_shellf, interp_shellf_capacity, stream);
		d_interp_local_shellf = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(required_shellf, stream, false);
		interp_shellf_capacity = required_shellf;
	}

	if (point_count > interp_weight_capacity)
	{
		GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_local_weight, interp_weight_capacity, stream);
		GpuRuntimeCoordinator::instance().release_device_buffer(d_interp_global_weight, interp_weight_capacity, stream);
		d_interp_local_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(point_count, stream, false);
		d_interp_global_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(point_count, stream, false);
		interp_weight_capacity = point_count;
	}

	if (!d_interp_err_flag)
	{
		d_interp_err_flag = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(1, stream, false);
	}
}

void Patch::ensure_interp_host_cache(int point_count, int variable_count)
{
	const int required_shellf = point_count * variable_count;
	if (required_shellf > interp_host_shellf_capacity)
	{
		delete[] h_interp_shellf;
		h_interp_shellf = new double[required_shellf];
		interp_host_shellf_capacity = required_shellf;
	}

	if (point_count > interp_host_weight_capacity)
	{
		delete[] h_interp_weight;
		h_interp_weight = new int[point_count];
		interp_host_weight_capacity = point_count;
	}
}
// buflog 1: with buffer points; 0 without
void Patch::checkPatch(bool buflog)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0)
	{
		if (buflog)
		{
			cout << " belong to level " << lev << endl;
			cout << " shape: [";
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
			cout << " range:" << "(";
			for (int i = 0; i < dim; i++)
			{
				cout << bbox[i] << ":" << bbox[dim + i];
				if (i < dim - 1)
					cout << ",";
				else
					cout << ")" << endl;
			}
		}
		else
		{
			cout << " belong to level " << lev << endl;
			cout << " shape: [";
			for (int i = 0; i < dim; i++)
			{
				cout << shape[i] - lli[i] - uui[i];
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
			cout << " range:" << "(";
			for (int i = 0; i < dim; i++)
			{
				cout << bbox[i] + lli[i] * getdX(i) << ":" << bbox[dim + i] - uui[i] * getdX(i);
				if (i < dim - 1)
					cout << ",";
				else
					cout << ")" << endl;
			}
		}
	}
}
void Patch::checkPatch(bool buflog, const int out_rank)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == out_rank)
	{
		cout << " out_rank = " << out_rank << endl;
		if (buflog)
		{
			cout << " belong to level " << lev << endl;
			cout << " shape: [";
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
			cout << " range:" << "(";
			for (int i = 0; i < dim; i++)
			{
				cout << bbox[i] << ":" << bbox[dim + i];
				if (i < dim - 1)
					cout << ",";
				else
					cout << ")" << endl;
			}
		}
		else
		{
			cout << " belong to level " << lev << endl;
			cout << " shape: [";
			for (int i = 0; i < dim; i++)
			{
				cout << shape[i] - lli[i] - uui[i];
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
			cout << " range:" << "(";
			for (int i = 0; i < dim; i++)
			{
				cout << bbox[i] + lli[i] * getdX(i) << ":" << bbox[dim + i] - uui[i] * getdX(i);
				if (i < dim - 1)
					cout << ",";
				else
					cout << ")" << endl;
			}
		}
	}
}

void Patch::Interp_Points(MyList<var> *VarList,
                          int NN, double **XX,
                          double *Shellf, int Symmetry)
{
    // NOTE: we do not Synchnize variables here, make sure of that before calling this routine
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    int ordn = 2 * ghost_width;
    MyList<var> *varl;
    int num_var = 0;
    varl = VarList;
    while (varl)
    {
        num_var++;
        varl = varl->next;
    }

    double *shellf;
    shellf = new double[NN * num_var];
    memset(shellf, 0, sizeof(double) * NN * num_var);

    int *weight;
    weight = new int[NN];
    memset(weight, 0, sizeof(int) * NN);

    double *DH, *llb, *uub;
    DH = new double[dim];

    for (int i = 0; i < dim; i++)
    {
        DH[i] = getdX(i);
    }
    llb = new double[dim];
    uub = new double[dim];

    // --------------------------------------------------------------------------
    // 保留 CPU 端的安全检查（只执行一次，开销极低）
    for (int j = 0; j < NN; j++) 
    {
        for (int i = 0; i < dim; i++)
        {
            if (isnan(XX[i][j])) {
                puts("Error: Patch::Interp_Points encounters NaN in input point coordinates.");
            }
            if (myrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
            {
                cout << "Patch::Interp_Points: point (";
                for (int k = 0; k < dim; k++)
                {
                    cout << XX[k][j];
                    if (k < dim - 1)
                        cout << ",";
                    else
                        cout << ") is out of current Patch." << endl;
                }
                MPI_Abort(MPI_COMM_WORLD, 1);
            }
        }
    }

    // ================== GPU 零拷贝重构部分 开始 ==================

    // 1. 分配 GPU 内存并将所有插值点坐标一次性拷贝至 GPU
	double *d_XX[3] = {nullptr, nullptr, nullptr};
	for (int i = 0; i < dim; i++) {
		d_XX[i] = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(NN);
		GpuRuntimeCoordinator::instance().copy_to_device(XX[i], d_XX[i], NN);
	}

	double *d_shellf = GpuRuntimeCoordinator::instance().allocate_device_buffer<double>(NN * num_var);

	int *d_weight = GpuRuntimeCoordinator::instance().allocate_device_buffer<int>(NN);

	GpuRuntimeCoordinator::instance().wait_for_memory_ops();

    // 2. Block-Centric 遍历：以 Block 为单位启动 Kernel
    MyList<Block> *Bp = blb;
    while (Bp)
    {
        Block *BP = Bp->data;

	        if (myrank == BP->rank)
	        {
	            // 计算当前 Block 包含幽灵区(ghost zones)的 Bounding Box
	            set_block_interp_bounds(BP, this, DH, llb, uub);

            // 提取多维参数（展平传递给 Kernel 以避免隐式 Host 指针问题）
            int shape_0 = BP->shape[0];
            int shape_1 = (dim > 1) ? BP->shape[1] : 1;
            int shape_2 = (dim > 2) ? BP->shape[2] : 1;

            double llb_0 = llb[0], llb_1 = (dim > 1) ? llb[1] : 0.0, llb_2 = (dim > 2) ? llb[2] : 0.0;
            double uub_0 = uub[0], uub_1 = (dim > 1) ? uub[1] : 0.0, uub_2 = (dim > 2) ? uub[2] : 0.0;

            varl = VarList;
            int k = 0;
            while (varl) 
            {
                gpu_global_interp_launch(
					BP->stream,
                    NN, dim, 
                    d_XX[0], d_XX[1], d_XX[2],
                    shape_0, shape_1, shape_2,
                    BP->d_X[0], BP->d_X[1], BP->d_X[2],
                    BP->d_fgfs[varl->data->sgfn],
                    llb_0, llb_1, llb_2, 
                    uub_0, uub_1, uub_2, 
                    ordn, varl->data->SoA[0], varl->data->SoA[1], varl->data->SoA[2], 
                    Symmetry, k, num_var, d_shellf, d_weight
                );
                varl = varl->next;
                k++;
            }
        }
        if (Bp == ble) break;
        Bp = Bp->next;
    }

	GpuRuntimeCoordinator::instance().wait_for_device();
	GpuRuntimeCoordinator::instance().copy_to_host(shellf, d_shellf, NN * num_var);
	GpuRuntimeCoordinator::instance().copy_to_host(weight, d_weight, NN);

    // 4. 清理 GPU 临时内存
	for (int i = 0; i < dim; i++) {
		GpuRuntimeCoordinator::instance().release_device_buffer(d_XX[i], NN);
	}
	GpuRuntimeCoordinator::instance().release_device_buffer(d_shellf, NN * num_var);
	GpuRuntimeCoordinator::instance().release_device_buffer(d_weight, NN);
    // ================== GPU 零拷贝重构部分 结束 ==================

    MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    int *Weight;
    Weight = new int[NN];
    MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

    for (int i = 0; i < NN; i++)
    {
        if (Weight[i] > 1)
        {
            if (myrank == 0)
                cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
            for (int j = 0; j < num_var; j++)
                Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
        }
        else if (Weight[i] == 0 && myrank == 0)
        {
            cout << "ERROR: Patch::Interp_Points fails to find point (";
            for (int j = 0; j < dim; j++)
            {
                cout << XX[j][i];
                if (j < dim - 1)
                    cout << ",";
                else
                    cout << ")";
            }
            cout << " on Patch (";
            for (int j = 0; j < dim; j++)
            {
                cout << bbox[j] << "+" << lli[j] * getdX(j);
                if (j < dim - 1)
                    cout << ",";
                else
                    cout << ")--";
            }
            cout << "(";
            for (int j = 0; j < dim; j++)
            {
                cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
                if (j < dim - 1)
                    cout << ",";
                else
                    cout << ")" << endl;
            }
            cout << "splited domains:" << endl;
            {
                MyList<Block> *Bp = blb;
                while (Bp)
                {
                    Block *BP = Bp->data;

	                    set_block_interp_bounds(BP, this, DH, llb, uub);
                    cout << "(";
                    for (int j = 0; j < dim; j++)
                    {
                        cout << llb[j] << ":" << uub[j];
                        if (j < dim - 1)
                            cout << ",";
                        else
                            cout << ")" << endl;
                    }
                    if (Bp == ble)
                        break;
                    Bp = Bp->next;
                }
            }
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }

    delete[] shellf;
    delete[] weight;
    delete[] Weight;
    delete[] DH;
    delete[] llb;
    delete[] uub;
	GpuRuntimeCoordinator::instance().wait_for_memory_ops();
}

// void Patch::Interp_Points(MyList<var> *VarList,
// 													int NN, double **XX,
// 													double *Shellf, int Symmetry)
// {
// 	// NOTE: we do not Synchnize variables here, make sure of that before calling this routine
// 	int myrank;
// 	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

// 	int ordn = 2 * ghost_width;
// 	MyList<var> *varl;
// 	int num_var = 0;
// 	varl = VarList;
// 	while (varl)
// 	{
// 		num_var++;
// 		varl = varl->next;
// 	}

// 	double *shellf;
// 	shellf = new double[NN * num_var];
// 	memset(shellf, 0, sizeof(double) * NN * num_var);

// 	// we use weight to monitor code, later some day we can move it for optimization
// 	int *weight;
// 	weight = new int[NN];
// 	memset(weight, 0, sizeof(int) * NN);

// 	double *DH, *llb, *uub;
// 	DH = new double[dim];

// 	for (int i = 0; i < dim; i++)
// 	{
// 		DH[i] = getdX(i);
// 	}
// 	llb = new double[dim];
// 	uub = new double[dim];

// 	GpuBatchInterp gpu_manager;
// 	for (int j = 0; j < NN; j++) // run along points
//     {
//         double pox[dim];
//         for (int i = 0; i < dim; i++)
//         {
//             pox[i] = XX[i][j];
// 			if (isnan(pox[i])) {
// 				puts("Error: Patch::Interp_Points encounters NaN in input point coordinates.");
// 			}
//             // [Original Logic] 边界检查，保持不变
//             if (myrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
//             {
//                 cout << "Patch::Interp_Points: point (";
//                 for (int k = 0; k < dim; k++)
//                 {
//                     cout << XX[k][j];
//                     if (k < dim - 1)
//                         cout << ",";
//                     else
//                         cout << ") is out of current Patch." << endl;
//                 }
//                 MPI_Abort(MPI_COMM_WORLD, 1);
//             }
//         }

//         MyList<Block> *Bp = blb;
//         bool notfind = true;
//         while (notfind && Bp) // run along Blocks
//         {
//             Block *BP = Bp->data;

//             bool flag = true;
//             // [Original Logic] Bounding Box 复杂判断，保持不变
//             for (int i = 0; i < dim; i++)
//             {
//                 llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
//                 uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];

//                 if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
//                 {
//                     flag = false;
//                     break;
//                 }
//             }

//             if (flag)
//             {
//                 notfind = false;
//                 if (myrank == BP->rank) // [Check] 只处理属于本 rank 的 Block
//                 {
//                     // ---> interpolation [Refactored for GPU]
//                     varl = VarList;
//                     int k = 0;
//                     while (varl) // run along variables
//                     {
//                         // 计算 Shellf 数组的目标索引
//                         int target_idx = j * num_var + k;

//                         // [GPU] 添加任务到队列，而不是立即计算
//                         // 注意：这里我们传入了 Block 的几何信息、场数据指针、点坐标和目标索引
//                         gpu_manager.add_task(
//                             BP->shape, 
//                             BP->X, 
//                             varl->data->SoA,           // SoA
//                             BP->fgfs[varl->data->sgfn], // 当前变量的 Field Data
//                             pox, 
//                             ordn, 
//                             Symmetry, 
//                             target_idx                  // 计算结果应该填回 Shellf 的哪个位置
//                         );

//                         varl = varl->next;
//                         k++;
//                     }
//                     // 标记权重 (CPU 侧逻辑保持不变，因为我们确信这个点被找到了)
//                     weight[j] = 1; 
//                 }
//             }
//             if (Bp == ble)
//                 break;
//             Bp = Bp->next;
//         }
//     }

//     // --------------------------------------------------------------------------
//     // [New Code] 循环结束，触发 GPU 批量计算并写回结果
//     // 这一步会执行：拷贝任务 -> Kernel 计算 -> 拷贝结果 -> 填入 Shellf
//     gpu_manager.execute(shellf);

// 	MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
// 	int *Weight;
// 	Weight = new int[NN];
// 	MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

// 	//  misc::tillherecheck("print me");

// 	for (int i = 0; i < NN; i++)
// 	{
// 		if (Weight[i] > 1)
// 		{
// 			if (myrank == 0)
// 				cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
// 			for (int j = 0; j < num_var; j++)
// 				Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
// 		}
// 		else if (Weight[i] == 0 && myrank == 0)
// 		{
// 			cout << "ERROR: Patch::Interp_Points fails to find point (";
// 			for (int j = 0; j < dim; j++)
// 			{
// 				cout << XX[j][i];
// 				if (j < dim - 1)
// 					cout << ",";
// 				else
// 					cout << ")";
// 			}
// 			cout << " on Patch (";
// 			for (int j = 0; j < dim; j++)
// 			{
// 				cout << bbox[j] << "+" << lli[j] * getdX(j);
// 				if (j < dim - 1)
// 					cout << ",";
// 				else
// 					cout << ")--";
// 			}
// 			cout << "(";
// 			for (int j = 0; j < dim; j++)
// 			{
// 				cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
// 				if (j < dim - 1)
// 					cout << ",";
// 				else
// 					cout << ")" << endl;
// 			}
// 			cout << "splited domains:" << endl;
// 			{
// 				MyList<Block> *Bp = blb;
// 				while (Bp)
// 				{
// 					Block *BP = Bp->data;

// 					for (int i = 0; i < dim; i++)
// 					{
// 						llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
// 						uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
// 					}
// 					cout << "(";
// 					for (int j = 0; j < dim; j++)
// 					{
// 						cout << llb[j] << ":" << uub[j];
// 						if (j < dim - 1)
// 							cout << ",";
// 						else
// 							cout << ")" << endl;
// 					}
// 					if (Bp == ble)
// 						break;
// 					Bp = Bp->next;
// 				}
// 			}
// 			MPI_Abort(MPI_COMM_WORLD, 1);
// 		}
// 	}

// 	delete[] shellf;
// 	delete[] weight;
// 	delete[] Weight;
// 	delete[] DH;
// 	delete[] llb;
// 	delete[] uub;
// }
void Patch::Interp_Points(MyList<var> *VarList,
													int NN, double **XX,
													double *Shellf, int Symmetry, MPI_Comm Comm_here)
{
	// NOTE: we do not Synchnize variables here, make sure of that before calling this routine
	int myrank, lmyrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_rank(Comm_here, &lmyrank);

	int ordn = 2 * ghost_width;
	MyList<var> *varl;
	int num_var = 0;
	varl = VarList;
	while (varl)
	{
		num_var++;
		varl = varl->next;
	}

	double *shellf;
	shellf = new double[NN * num_var];
	memset(shellf, 0, sizeof(double) * NN * num_var);

	// we use weight to monitor code, later some day we can move it for optimization
	int *weight;
	weight = new int[NN];
	memset(weight, 0, sizeof(int) * NN);

	double *DH, *llb, *uub;
	DH = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		DH[i] = getdX(i);
	}
	llb = new double[dim];
	uub = new double[dim];

	for (int j = 0; j < NN; j++) // run along points
	{
		double pox[dim];
		for (int i = 0; i < dim; i++)
		{
			pox[i] = XX[i][j];
			if (lmyrank == 0 && (XX[i][j] < bbox[i] + lli[i] * DH[i] || XX[i][j] > bbox[dim + i] - uui[i] * DH[i]))
			{
				cout << "Patch::Interp_Points: point (";
				for (int k = 0; k < dim; k++)
				{
					cout << XX[k][j];
					if (k < dim - 1)
						cout << ",";
					else
						cout << ") is out of current Patch." << endl;
				}
				MPI_Abort(MPI_COMM_WORLD, 1);
			}
		}

		MyList<Block> *Bp = blb;
		bool notfind = true;
		while (notfind && Bp) // run along Blocks
		{
			Block *BP = Bp->data;

			bool flag = true;
				set_block_interp_bounds(BP, this, DH, llb, uub);
				for (int i = 0; i < dim; i++)
				{
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
					if (XX[i][j] - llb[i] < -DH[i] / 2 || XX[i][j] - uub[i] > DH[i] / 2)
					{
						flag = false;
					break;
				}
			}

			if (flag)
			{
				notfind = false;
				if (myrank == BP->rank)
				{
					//---> interpolation
					varl = VarList;
					int k = 0;
					while (varl) // run along variables
					{
						//              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
						//	  		                                    pox,ordn,varl->data->SoA,Symmetry);
						f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[j * num_var + k],
														pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
						varl = varl->next;
						k++;
					}
					weight[j] = 1;
				}
			}
			if (Bp == ble)
				break;
			Bp = Bp->next;
		}
	}

	MPI_Allreduce(shellf, Shellf, NN * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);
	int *Weight;
	Weight = new int[NN];
	MPI_Allreduce(weight, Weight, NN, MPI_INT, MPI_SUM, Comm_here);

	//  misc::tillherecheck("print me");
	//  if(lmyrank == 0) cout<<"myrank = "<<myrank<<"print me"<<endl;

	for (int i = 0; i < NN; i++)
	{
		if (Weight[i] > 1)
		{
			if (lmyrank == 0)
				cout << "WARNING: Patch::Interp_Points meets multiple weight" << endl;
			for (int j = 0; j < num_var; j++)
				Shellf[j + i * num_var] = Shellf[j + i * num_var] / Weight[i];
		}
#if 0 // for not involved levels, this may fail     
		 else if(Weight[i] == 0 && lmyrank == 0)
		 {
			 cout<<"ERROR: Patch::Interp_Points fails to find point (";
			 for(int j=0;j<dim;j++)
			 {
		cout<<XX[j][i];
		if(j<dim-1) cout<<",";
		else        cout<<")";
			 }
			 cout<<" on Patch (";
			 for(int j=0;j<dim;j++)
			 {
		cout<<bbox[j]<<"+"<<lli[j]*getdX(j);
		if(j<dim-1) cout<<",";
		else        cout<<")--";
			 }
			 cout<<"(";
			 for(int j=0;j<dim;j++)
			 {
		cout<<bbox[dim+j]<<"-"<<uui[j]*getdX(j);
		if(j<dim-1) cout<<",";
		else        cout<<")"<<endl;
			 }
#if 0
			 checkBlock();
#else
	cout<<"splited domains:"<<endl;
	{
		 MyList<Block> *Bp=blb;
		 while(Bp)
		 {
	Block *BP=Bp->data;

	for(int i=0;i<dim;i++)
	{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif    
					llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +(ghost_width-0.5)*DH[i];
					uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-(ghost_width-0.5)*DH[i];
#else
#ifdef Cell
					llb[i] = (feq(BP->bbox[i]    ,bbox[i]    ,DH[i]/2)) ? BP->bbox[i]+lli[i]*DH[i]     : BP->bbox[i]    +ghost_width*DH[i];
					uub[i] = (feq(BP->bbox[dim+i],bbox[dim+i],DH[i]/2)) ? BP->bbox[dim+i]-uui[i]*DH[i] : BP->bbox[dim+i]-ghost_width*DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif 
	}       
			 cout<<"(";
			 for(int j=0;j<dim;j++)
			 {
		cout<<llb[j]<<":"<<uub[j];
		if(j<dim-1) cout<<",";
		else        cout<<")"<<endl;
			 }
	if(Bp == ble) break;
	Bp=Bp->next;
		 }
	}
#endif       
			 MPI_Abort(MPI_COMM_WORLD,1);
		 }
#endif
	}

	delete[] shellf;
	delete[] weight;
	delete[] Weight;
	delete[] DH;
	delete[] llb;
	delete[] uub;
}
void Patch::checkBlock()
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	if (myrank == 0)
	{
		MyList<Block> *BP = blb;
		while (BP)
		{
			BP->data->checkBlock();
			if (BP == ble)
				break;
			BP = BP->next;
		}
	}
}
double Patch::getdX(int dir)
{
	if (dir < 0 || dir >= dim)
	{
		cout << "Patch::getdX: error input dir = " << dir << ", this Patch has direction (0," << dim - 1 << ")" << endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	double h;
	h = patch_spacing(bbox, shape, dir, "Patch::getdX");
	return h;
}
bool Patch::Interp_ONE_Point(MyList<var> *VarList, double *XX,
														 double *Shellf, int Symmetry)
{
	// NOTE: we do not Synchnize variables here, make sure of that before calling this routine
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int ordn = 2 * ghost_width;
	MyList<var> *varl;
	int num_var = 0;
	varl = VarList;
	while (varl)
	{
		num_var++;
		varl = varl->next;
	}

	double *shellf;
	shellf = new double[num_var];
	memset(shellf, 0, sizeof(double) * num_var);

	double *DH, *llb, *uub;
	DH = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		DH[i] = getdX(i);
	}
	llb = new double[dim];
	uub = new double[dim];

	double pox[dim];
	for (int i = 0; i < dim; i++)
	{
		pox[i] = XX[i];
		// has excluded the buffer points
		if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
		{
			delete[] shellf;
			delete[] DH;
			delete[] llb;
			delete[] uub;
			return false; // out of current patch,
										// remember to delete the allocated arrays before return!!!
		}
	}

	MyList<Block> *Bp = blb;
	bool notfind = true;
	while (notfind && Bp) // run along Blocks
	{
		Block *BP = Bp->data;

		bool flag = true;
			set_block_interp_bounds(BP, this, DH, llb, uub);
			for (int i = 0; i < dim; i++)
			{
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
				if (XX[i] - llb[i] < -DH[i] / 2 || XX[i] - uub[i] > DH[i] / 2)
				{
					flag = false;
				break;
			}
		}

		if (flag)
		{
			notfind = false;
			if (myrank == BP->rank)
			{
				//---> interpolation
				varl = VarList;
				int k = 0;
				while (varl) // run along variables
				{
					//              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
					//	  		                                    pox,ordn,varl->data->SoA,Symmetry);
					f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[k],
													pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
					varl = varl->next;
					k++;
				}
			}
		}
		if (Bp == ble)
			break;
		Bp = Bp->next;
	}

	if (notfind && myrank == 0)
	{
		cout << "ERROR: Patch::Interp_Points fails to find point (";
		for (int j = 0; j < dim; j++)
		{
			cout << XX[j];
			if (j < dim - 1)
				cout << ",";
			else
				cout << ")";
		}
		cout << " on Patch (";
		for (int j = 0; j < dim; j++)
		{
			cout << bbox[j] << "+" << lli[j] * getdX(j);
			if (j < dim - 1)
				cout << ",";
			else
				cout << ")--";
		}
		cout << "(";
		for (int j = 0; j < dim; j++)
		{
			cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
			if (j < dim - 1)
				cout << ",";
			else
				cout << ")" << endl;
		}
		cout << "splited domains:" << endl;
		{
			MyList<Block> *Bp = blb;
			while (Bp)
			{
				Block *BP = Bp->data;

					set_block_interp_bounds(BP, this, DH, llb, uub);
				cout << "(";
				for (int j = 0; j < dim; j++)
				{
					cout << llb[j] << ":" << uub[j];
					if (j < dim - 1)
						cout << ",";
					else
						cout << ")" << endl;
				}
				if (Bp == ble)
					break;
				Bp = Bp->next;
			}
		}
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Allreduce(shellf, Shellf, num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	delete[] shellf;
	delete[] DH;
	delete[] llb;
	delete[] uub;

	return true;
}
bool Patch::Interp_ONE_Point(MyList<var> *VarList, double *XX,
														 double *Shellf, int Symmetry, MPI_Comm Comm_here)
{
	// NOTE: we do not Synchnize variables here, make sure of that before calling this routine
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	int ordn = 2 * ghost_width;
	MyList<var> *varl;
	int num_var = 0;
	varl = VarList;
	while (varl)
	{
		num_var++;
		varl = varl->next;
	}

	double *shellf;
	shellf = new double[num_var];
	memset(shellf, 0, sizeof(double) * num_var);

	double *DH, *llb, *uub;
	DH = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		DH[i] = getdX(i);
	}
	llb = new double[dim];
	uub = new double[dim];

	double pox[dim];
	for (int i = 0; i < dim; i++)
	{
		pox[i] = XX[i];
		// has excluded the buffer points
		if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
		{
			delete[] shellf;
			delete[] DH;
			delete[] llb;
			delete[] uub;
			return false; // out of current patch,
										// remember to delete the allocated arrays before return!!!
		}
	}

	MyList<Block> *Bp = blb;
	bool notfind = true;
	while (notfind && Bp) // run along Blocks
	{
		Block *BP = Bp->data;

		bool flag = true;
		for (int i = 0; i < dim; i++)
		{
// NOTE: our dividing structure is (exclude ghost)
// -1 0
//       1  2
// so (0,1) does not belong to any part for vertex structure
// here we put (0,0.5) to left part and (0.5,1) to right part
// BUT for cell structure the bbox is (-1.5,0.5) and (0.5,2.5), there is no missing region at all
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
			llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
			uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
			llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
			uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
			if (XX[i] - llb[i] < -DH[i] / 2 || XX[i] - uub[i] > DH[i] / 2)
			{
				flag = false;
				break;
			}
		}

		if (flag)
		{
			notfind = false;
			if (myrank == BP->rank)
			{
// test old code
#if 0
#define floorint(a) ((a) < 0 ? int(a) - 1 : int(a))
//---> interpolation
								int ixl,iyl,izl,ixu,iyu,izu;
				double Delx,Dely,Delz;

		ixl = 1+floorint((pox[0]-BP->X[0][0])/DH[0]);
		 	iyl = 1+floorint((pox[1]-BP->X[1][0])/DH[1]);
		 	izl = 1+floorint((pox[2]-BP->X[2][0])/DH[2]);

		int nn=ordn/2;

		ixl = ixl-nn;
		iyl = iyl-nn;
		izl = izl-nn;
		 
		int tmi;
		tmi = (Symmetry==2)?-1:0;
		if(ixl<tmi) ixl=tmi;
		 	if(iyl<tmi) iyl=tmi;
		tmi = (Symmetry>0)?-1:0;
		 	if(izl<tmi) izl=tmi;
			
		 	if(ixl+ordn>BP->shape[0]) ixl=BP->shape[0]-ordn;
		 	if(iyl+ordn>BP->shape[1]) iyl=BP->shape[1]-ordn;
		 	if(izl+ordn>BP->shape[2]) izl=BP->shape[2]-ordn;
// support cell center
		if(ixl>=0) Delx = ( pox[0] - BP->X[0][ixl] )/ DH[0];
		else       Delx = ( pox[0] + BP->X[0][0] )/ DH[0];
								if(iyl>=0) Dely = ( pox[1] - BP->X[1][iyl] )/ DH[1];
		else       Dely = ( pox[1] + BP->X[1][0] )/ DH[1];
								if(izl>=0) Delz = ( pox[2] - BP->X[2][izl] )/ DH[2];
		else       Delz = ( pox[2] + BP->X[2][0] )/ DH[2];
//change to fortran index
								ixl++;
		 	iyl++;
		 	izl++;
		 	ixu = ixl + ordn - 1;
		 	iyu = iyl + ordn - 1;
		 	izu = izl + ordn - 1;
				varl=VarList;
		int j=0;
				while(varl)
		{
								 f_interp_2(BP->shape,BP->fgfs[varl->data->sgfn],shellf[j],ixl,ixu,iyl,iyu,izl,izu,Delx,Dely,Delz,
																		 ordn,varl->data->SoA,Symmetry);
		 varl=varl->next;
		 j++;
		} //varl
#else
				//---> interpolation
				varl = VarList;
				int k = 0;
				while (varl) // run along variables
				{
					//              shellf[j*num_var+k] = Parallel::global_interp(dim,BP->shape,BP->X,BP->fgfs[varl->data->sgfn],
					//	  		                                    pox,ordn,varl->data->SoA,Symmetry);
					f_global_interp(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], shellf[k],
													pox[0], pox[1], pox[2], ordn, varl->data->SoA, Symmetry);
					varl = varl->next;
					k++;
				}
#endif
			}
		}
		if (Bp == ble)
			break;
		Bp = Bp->next;
	}

	if (notfind && myrank == 0)
	{
		cout << "ERROR: Patch::Interp_Points fails to find point (";
		for (int j = 0; j < dim; j++)
		{
			cout << XX[j];
			if (j < dim - 1)
				cout << ",";
			else
				cout << ")";
		}
		cout << " on Patch (";
		for (int j = 0; j < dim; j++)
		{
			cout << bbox[j] << "+" << lli[j] * getdX(j);
			if (j < dim - 1)
				cout << ",";
			else
				cout << ")--";
		}
		cout << "(";
		for (int j = 0; j < dim; j++)
		{
			cout << bbox[dim + j] << "-" << uui[j] * getdX(j);
			if (j < dim - 1)
				cout << ",";
			else
				cout << ")" << endl;
		}
#if 0
			 checkBlock();
#else
		cout << "splited domains:" << endl;
		{
			MyList<Block> *Bp = blb;
			while (Bp)
			{
				Block *BP = Bp->data;

				for (int i = 0; i < dim; i++)
				{
#ifdef Vertex
#ifdef Cell
#error Both Cell and Vertex are defined
#endif
					llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + (ghost_width - 0.5) * DH[i];
					uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - (ghost_width - 0.5) * DH[i];
#else
#ifdef Cell
					llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? BP->bbox[i] + lli[i] * DH[i] : BP->bbox[i] + ghost_width * DH[i];
					uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? BP->bbox[dim + i] - uui[i] * DH[i] : BP->bbox[dim + i] - ghost_width * DH[i];
#else
#error Not define Vertex nor Cell
#endif
#endif
				}
				cout << "(";
				for (int j = 0; j < dim; j++)
				{
					cout << llb[j] << ":" << uub[j];
					if (j < dim - 1)
						cout << ",";
					else
						cout << ")" << endl;
				}
				if (Bp == ble)
					break;
				Bp = Bp->next;
			}
		}
#endif
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	MPI_Allreduce(shellf, Shellf, num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

	delete[] shellf;
	delete[] DH;
	delete[] llb;
	delete[] uub;

	return true;
}
// find maximum of abstract value, XX store position for maximum, Shellf store maximum themselvs
void Patch::Find_Maximum(MyList<var> *VarList, double *XX,
												 double *Shellf)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	MyList<var> *varl;
	int num_var = 0;
	varl = VarList;
	while (varl)
	{
		num_var++;
		varl = varl->next;
	}

	double *shellf, *xx;
	shellf = new double[num_var];
	xx = new double[dim * num_var];
	memset(shellf, 0, sizeof(double) * num_var);
	memset(xx, 0, sizeof(double) * dim * num_var);

	double *DH;
	int *llb, *uub;
	DH = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		DH[i] = getdX(i);
	}

	llb = new int[dim];
	uub = new int[dim];

	MyList<Block> *Bp = blb;
	while (Bp) // run along Blocks
	{
		Block *BP = Bp->data;

		if (myrank == BP->rank)
		{

			for (int i = 0; i < dim; i++)
			{
				llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
				uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
			}

			varl = VarList;
			int k = 0;
			double tmp, tmpx[dim];
			while (varl) // run along variables
			{
				f_find_maximum(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], tmp, tmpx, llb, uub);
				if (tmp > shellf[k])
				{
					shellf[k] = tmp;
					for (int i = 0; i < dim; i++)
						xx[dim * k + i] = tmpx[i];
				}
				varl = varl->next;
				k++;
			}
		}

		if (Bp == ble)
			break;
		Bp = Bp->next;
	}

	struct mloc
	{
		double val;
		int rank;
	};

	mloc *IN, *OUT;
	IN = new mloc[num_var];
	OUT = new mloc[num_var];
	for (int i = 0; i < num_var; i++)
	{
		IN[i].val = shellf[i];
		IN[i].rank = myrank;
	}

	MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

	for (int i = 0; i < num_var; i++)
	{
		Shellf[i] = OUT[i].val;
		if (myrank != OUT[i].rank)
			for (int k = 0; k < 3; k++)
				xx[3 * i + k] = 0;
	}

	MPI_Allreduce(xx, XX, dim * num_var, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	delete[] IN;
	delete[] OUT;
	delete[] shellf;
	delete[] xx;
	delete[] DH;
	delete[] llb;
	delete[] uub;
}
void Patch::Find_Maximum(MyList<var> *VarList, double *XX,
												 double *Shellf, MPI_Comm Comm_here)
{
	int myrank;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

	MyList<var> *varl;
	int num_var = 0;
	varl = VarList;
	while (varl)
	{
		num_var++;
		varl = varl->next;
	}

	double *shellf, *xx;
	shellf = new double[num_var];
	xx = new double[dim * num_var];
	memset(shellf, 0, sizeof(double) * num_var);
	memset(xx, 0, sizeof(double) * dim * num_var);

	double *DH;
	int *llb, *uub;
	DH = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		DH[i] = getdX(i);
	}

	llb = new int[dim];
	uub = new int[dim];

	MyList<Block> *Bp = blb;
	while (Bp) // run along Blocks
	{
		Block *BP = Bp->data;

		if (myrank == BP->rank)
		{

			for (int i = 0; i < dim; i++)
			{
				llb[i] = (feq(BP->bbox[i], bbox[i], DH[i] / 2)) ? lli[i] : ghost_width;
				uub[i] = (feq(BP->bbox[dim + i], bbox[dim + i], DH[i] / 2)) ? uui[i] : ghost_width;
			}

			varl = VarList;
			int k = 0;
			double tmp, tmpx[dim];
			while (varl) // run along variables
			{
				f_find_maximum(BP->shape, BP->X[0], BP->X[1], BP->X[2], BP->fgfs[varl->data->sgfn], tmp, tmpx, llb, uub);
				if (tmp > shellf[k])
				{
					shellf[k] = tmp;
					for (int i = 0; i < dim; i++)
						xx[dim * k + i] = tmpx[i];
				}
				varl = varl->next;
				k++;
			}
		}

		if (Bp == ble)
			break;
		Bp = Bp->next;
	}

	struct mloc
	{
		double val;
		int rank;
	};

	mloc *IN, *OUT;
	IN = new mloc[num_var];
	OUT = new mloc[num_var];
	for (int i = 0; i < num_var; i++)
	{
		IN[i].val = shellf[i];
		IN[i].rank = myrank;
	}

	MPI_Allreduce(IN, OUT, num_var, MPI_DOUBLE_INT, MPI_MAXLOC, Comm_here);

	for (int i = 0; i < num_var; i++)
	{
		Shellf[i] = OUT[i].val;
		if (myrank != OUT[i].rank)
			for (int k = 0; k < 3; k++)
				xx[3 * i + k] = 0;
	}

	MPI_Allreduce(xx, XX, dim * num_var, MPI_DOUBLE, MPI_SUM, Comm_here);

	delete[] IN;
	delete[] OUT;
	delete[] shellf;
	delete[] xx;
	delete[] DH;
	delete[] llb;
	delete[] uub;
}
// if the given point locates in the present Patch return true
// otherwise return false
bool Patch::Find_Point(double *XX)
{
	double *DH;
	DH = new double[dim];

	for (int i = 0; i < dim; i++)
	{
		DH[i] = getdX(i);
	}

	for (int i = 0; i < dim; i++)
	{
		// has excluded the buffer points
		if (XX[i] < bbox[i] + lli[i] * DH[i] - DH[i] / 100 || XX[i] > bbox[dim + i] - uui[i] * DH[i] + DH[i] / 100)
		{
			delete[] DH;
			return false; // out of current patch,
										// remember to delete the allocated arrays before return!!!
		}
	}

	delete[] DH;

	return true;
}
