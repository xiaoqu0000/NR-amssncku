#ifndef DATA_MOTION_HELPER_H
#define DATA_MOTION_HELPER_H

#include "MyList.h"
#include "MPatch.h"
#include "Block.h"
#include "var.h"
#include "bssn_gpu_class.h"

namespace DataMotionHelper {
    void push_all_data_to_device(MyList<Patch>* patch_list, int mpi_rank, MyList<var>* variable_list);
    void push_all_data_to_device(cgh* grid_hierarchy, int mpi_rank, MyList<var>* variable_list);
    void pull_all_data_to_host(MyList<Patch>* patch_list, int mpi_rank, MyList<var>* variable_list);
    void pull_all_data_to_host(cgh* grid_hierarchy, int mpi_rank, MyList<var>* variable_list);
}

#endif
