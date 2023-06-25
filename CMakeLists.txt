cmake_minimum_required(VERSION 3.10)

# Set the project name, description and programming language
project(GLAC VERSION 1.0
DESCRIPTION "Lattice Quantum Chromo Dynamics program for generating pure gauge field configurations and applying gradient flow on them"
LANGUAGES CXX)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Add compile options
# add_compile_options(-O3 -march=native -Werror -Wall -Wextra)

# Set the compilers to mpicc and mpicxx
set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_C_COMPILER mpicc)

# Set the output directory for the binary files
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Include directories
include_directories(${CMAKE_SOURCE_DIR}/lib)

# Add subdirectories
add_subdirectory(${CMAKE_SOURCE_DIR}/src/math)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/actions)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/config)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/flow)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/io)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/observables)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/parallelization)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/tests_old)
add_subdirectory(tests)
add_subdirectory(tests_mpi)

# Add the executable
add_executable(GLAC ${CMAKE_SOURCE_DIR}/src/main.cpp ${CMAKE_SOURCE_DIR}/src/system.cpp)

# Get the MPI compile and link flags
execute_process(COMMAND mpicxx --showme:compile OUTPUT_VARIABLE MPI_COMPILE_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND mpicxx --showme:link OUTPUT_VARIABLE MPI_LINK_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)

# Set the MPI compile and link flags
set_target_properties(GLAC PROPERTIES
COMPILE_FLAGS "${MPI_COMPILE_FLAGS}"
LINK_FLAGS "${MPI_LINK_FLAGS}")

# Link the libraries to the executable
# target_link_libraries(GLAC math actions config flow io observables parallelization tests)
target_link_libraries(GLAC config flow io math observables parallelization actions tests_old)

# Add a custom target to run the tests
add_custom_target(run_tests COMMAND tests)
add_custom_target(run_mpi_tests COMMAND mpi_tests)