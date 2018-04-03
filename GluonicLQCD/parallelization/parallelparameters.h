#ifndef PARALLELPARAMETERS_H
#define PARALLELPARAMETERS_H

#include <mpi.h>

namespace Parallel {
class ParallelParameters {
public:
    // MPI Groups
    static MPI_Group ACTIVE_GROUP,WORLD_GROUP;
    static MPI_Comm ACTIVE_COMM;

    // Variable for storing if processor is active(will always seek maximum possible of cores divisible to 2^n
    static bool active;

    ParallelParameters();
    ~ParallelParameters();
};
}

#endif // PARALLELPARAMETERS_H
