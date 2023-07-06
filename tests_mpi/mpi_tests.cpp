#include <lattice/test_lattice_operations.h>
#include <parallelization/communicator.h>

#include <iomanip>
#include <iostream>

int main(int numberOfArguments, char* cmdLineArguments[])
{
  Parallel::Communicator::init(numberOfArguments, cmdLineArguments);

  std::cout << "Starting MPI tests\n";

  bool passed = true;

  if (Parallel::ParallelParameters::active)
  {
    passed |= testLatticeOperations();
  }

  return passed && MPI_Finalize();
}