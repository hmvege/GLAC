cmake_minimum_required(VERSION 3.20)  # 3.20+ is plenty
project(GLAC VERSION 1.0
  DESCRIPTION "LQCD program for generating pure gauge configs and gradient flow"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# MPI discovery (no compiler overrides)
find_package(MPI 3.0 REQUIRED COMPONENTS CXX)

# Optional: your include dir
include_directories(${CMAKE_SOURCE_DIR}/lib)

# Catch2 (as you had it)
option(CATCH_BUILD_TESTING "" OFF)
option(CATCH_BUILD_EXAMPLES "" OFF)
option(CATCH_BUILD_EXTRA_TESTS "" OFF)
option(CATCH_INSTALL_DOCS "" OFF)
option(CATCH_INSTALL_EXTRAS "" OFF)
add_subdirectory(external/Catch2)

# Your subdirs
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

# Executable
add_executable(GLAC
  ${CMAKE_SOURCE_DIR}/src/main.cpp
  ${CMAKE_SOURCE_DIR}/src/system.cpp)

target_compile_features(GLAC PRIVATE cxx_std_17)
# if(CMAKE_BUILD_TYPE STREQUAL "Release")
#   target_compile_options(GLAC PRIVATE -O3 -march=native)
# endif()
# target_compile_options(GLAC PRIVATE -Werror -Wall -Wextra)

# Link MPI the right way
target_link_libraries(GLAC PRIVATE MPI::MPI_CXX
  config flow io math observables parallelization actions tests_old)

# If you want the runtime dir standardized
set_target_properties(GLAC PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# CTest integration for MPI tests (example)
include(CTest)
# add_test(NAME glac_np4 COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4
#   ${MPIEXEC_PREFLAGS} $<TARGET_FILE:GLAC> ${MPIEXEC_POSTFLAGS})
