/*!
 * \class ParallelParameters
 *
 * \brief ParallelParameters holds two groups and one communicator.
 *
 * MPI Groups used for communications. Needed in case some processors remain inactive(i.e. we are unable to have zero in remainder hven dividing by 2^n).
 *
 * \author Mathias M. Vege
 * \version 1.0
 * \date 2017-2019
 * \copyright MIT Licence
 */
#ifndef PARALLELPARAMETERS_H
#define PARALLELPARAMETERS_H

#include <mpi.h>

/*!
 * \namespace Parallel
 *
 * \brief Parallel contains all of the relevant methods for communicating between the lattices.
 */

namespace Parallel
{
class ParallelParameters
{
public:
    //! The WORLD_GROUP is the communicator used all processors allocated to the program.
    static MPI_Group WORLD_GROUP;
    /*!
     * The ACTIVE_GROUP is the communicator used for the processors involved directly in the calculation.
     *
     * This is due to the fact that when allocating processors to the program, it may exceed what is the optimal number of processors for the specified sub-lattice geoemtry.
     */
    static MPI_Group ACTIVE_GROUP;
    //! The communicator ACTIVE_COMM is the communicator for the ACTIVE_GROUP.
    static MPI_Comm ACTIVE_COMM;

    //! Variable for storing if processor is active(will always seek maximum possible of cores divisible to 2^n).
    static bool active;

    ParallelParameters();
    ~ParallelParameters();
};
}

#endif // PARALLELPARAMETERS_H
