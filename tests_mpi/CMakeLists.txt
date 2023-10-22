find_package(Catch2 3 REQUIRED)

# Add the test files
set(MPI_TEST_SOURCES
    # mpi_tests.cpp
    # lattice/test_lattice.cpp
    lattice/test_actions.cpp
    # lattice/action_test.cpp
    # lattice/observable_test.cpp
    # io/io_test.cpp
)

# Add the test executable(These tests can use the Catch2-provided main)
add_executable(mpi_tests ${MPI_TEST_SOURCES} mpi_tests.cpp)

# Link the necessary libraries
target_link_libraries(mpi_tests PRIVATE math actions ${MPI_LIBRARIES} Catch2::Catch2)
target_include_directories(mpi_tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
