#include "parallelparameters.h"

// MPI Groups
MPI_Group Parallel::ParallelParameters::ACTIVE_GROUP;
MPI_Group Parallel::ParallelParameters::WORLD_GROUP;
MPI_Comm Parallel::ParallelParameters::ACTIVE_COMM;

// Variable for storing if processor is active(will always seek maximum possible of cores divisible to 2^n
bool Parallel::ParallelParameters::active = true;

Parallel::ParallelParameters::ParallelParameters()
{
}

Parallel::ParallelParameters::~ParallelParameters()
{
}
