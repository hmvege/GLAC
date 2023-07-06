#include "parallelparameters.h"

// MPI Groups
MPI_Group Parallel::ParallelParameters::ACTIVE_GROUP = MPI_GROUP_NULL;
MPI_Group Parallel::ParallelParameters::WORLD_GROUP = MPI_GROUP_NULL;
MPI_Comm Parallel::ParallelParameters::ACTIVE_COMM = MPI_COMM_NULL;

// Variable for storing if processor is active(will always seek maximum possible of cores divisible to 2^n
bool Parallel::ParallelParameters::active = true;

Parallel::ParallelParameters::ParallelParameters()
{
}

Parallel::ParallelParameters::~ParallelParameters()
{
}
