#include <mpi.h>

#include <catch2/catch_all.hpp>

TEST_CASE("MPI Test", "[mpi]")
{
  SECTION("Sub-sec1")
  {
    SECTION("Subsub")
    {
      int rank;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);

      double b = 3.14;
      CAPTURE(b);

      INFO("Running test for rank " + rank);

      if (rank == 0)
      {
        // Test something on the master process
        std::cout << "Ok from 0\n";
        REQUIRE(1 + 1 == 2);
      }
      else if (rank == 2)
      {
        std::cout << "Not Ok from 2\n";
        const int x = 1;
        REQUIRE_FALSE(1 + x == 3);
      }
      else
      {
        // Test something on the worker processes
        std::cout << "Ok from others\n";
        REQUIRE(2 + 2 == 4);
      }
    }
  }

  SECTION("Sub-sec2")
  {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    double b = 3.14;
    CAPTURE(b);
    // REQUIRE(1 + 1 == 4);

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
}

// TEST_CASE("MPI counter test", "[mpi]")
// {
//   int rank;
//   MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//   double b = 3.14 * 100;
//   CAPTURE(b);
//   REQUIRE(1 + 1 == 40);
// }