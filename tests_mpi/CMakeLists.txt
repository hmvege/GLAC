# Add the test files
set(MPI_TEST_SOURCES
    # lattice/test_lattice.cpp
    lattice/test_actions.cpp
    # lattice/action_test.cpp
    # lattice/observable_test.cpp
    # io/io_test.cpp
    mpi_tests.cpp
)

# Add the test executable(These tests can use the Catch2-provided main)
add_executable(mpi_tests ${MPI_TEST_SOURCES})
target_compile_features(mpi_tests PRIVATE cxx_std_17)

# Link the necessary libraries
target_link_libraries(mpi_tests PRIVATE math actions MPI::MPI_CXX Catch2::Catch2)
target_include_directories(mpi_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Registers CMake tests
include(CTest)

# Allow changing number of ranks without editing CMake
set(MPI_TEST_NPROCS 4 CACHE STRING "Number of ranks for mpi_tests")

add_test(NAME mpi_tests_np${MPI_TEST_NPROCS}
  COMMAND
    ${MPIEXEC_EXECUTABLE}
    ${MPIEXEC_NUMPROC_FLAG} ${MPI_TEST_NPROCS}
    ${MPIEXEC_PREFLAGS}
    $<TARGET_FILE:mpi_tests>
    ${MPIEXEC_POSTFLAGS}
)

# # Avoid running multiple mpirun jobs in parallel (optional)
# set_tests_properties(mpi_tests_np${MPI_TEST_NPROCS} PROPERTIES RUN_SERIAL TRUE)
