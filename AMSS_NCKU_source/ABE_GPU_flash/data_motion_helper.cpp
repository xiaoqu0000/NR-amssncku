#include "data_motion_helper.h"

#include "MyList.h"
#include "MPatch.h"
#include "Block.h"
#include "var.h"
#include "bssn_gpu_class.h"

namespace {

template<typename TransferFn>
void for_each_local_block(MyList<Patch>* patch_list, int mpi_rank, MyList<var>* variable_list, TransferFn transfer_fn) {
    for (MyList<Patch>* patch_cursor = patch_list; patch_cursor; patch_cursor = patch_cursor->next) {
        for (MyList<Block>* block_cursor = patch_cursor->data->blb; block_cursor; block_cursor = block_cursor->next) {
            Block* block = block_cursor->data;
            if (mpi_rank == block->rank) {
                transfer_fn(block, variable_list);
            }
        }
    }
}

void wait_for_local_block_streams(MyList<Patch>* patch_list, int mpi_rank) {
    for (MyList<Patch>* patch_cursor = patch_list; patch_cursor; patch_cursor = patch_cursor->next) {
        for (MyList<Block>* block_cursor = patch_cursor->data->blb; block_cursor; block_cursor = block_cursor->next) {
            Block* block = block_cursor->data;
            if (mpi_rank == block->rank) {
                CUDA_CHECK(cudaStreamSynchronize(block->stream));
            }
            if (block_cursor == patch_cursor->data->ble) {
                break;
            }
        }
    }
}

}  // namespace

namespace DataMotionHelper {

void pull_all_data_to_host(MyList<Patch>* patch_list, int mpi_rank, MyList<var>* variable_list) {
    for_each_local_block(patch_list, mpi_rank, variable_list, [](Block* block, MyList<var>* variables) {
        block->move_to_cpu(variables);
    });
    wait_for_local_block_streams(patch_list, mpi_rank);
}

void push_all_data_to_device(cgh* grid_hierarchy, int mpi_rank, MyList<var>* variable_list) {
    for (int level_index = 0; level_index < grid_hierarchy->levels; ++level_index) {
        DataMotionHelper::push_all_data_to_device(grid_hierarchy->PatL[level_index], mpi_rank, variable_list);
    }
}

void push_all_data_to_device(MyList<Patch>* patch_list, int mpi_rank, MyList<var>* variable_list) {
    for_each_local_block(patch_list, mpi_rank, variable_list, [](Block* block, MyList<var>* variables) {
        block->move_to_gpu(variables);
    });
    wait_for_local_block_streams(patch_list, mpi_rank);
}

void pull_all_data_to_host(cgh* grid_hierarchy, int mpi_rank, MyList<var>* variable_list) {
    for (int level_index = 0; level_index < grid_hierarchy->levels; ++level_index) {
        DataMotionHelper::pull_all_data_to_host(grid_hierarchy->PatL[level_index], mpi_rank, variable_list);
    }
}

}  // namespace DataMotionHelper
