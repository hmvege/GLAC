#include <mpi.h>

#include <catch2/catch_all.hpp>

TEST_CASE("MPI Test", "[mpi]")
{
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (rank == 0)
  {
    // Test something on the master process
    std::cout << "Ok from 0\n";
    REQUIRE(1 + 1 == 2);
  }
  else
  {
    // Test something on the worker processes
    std::cout << "Ok from others\n";
    REQUIRE(2 + 2 == 4);
  }
}