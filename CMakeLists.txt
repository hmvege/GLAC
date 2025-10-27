cmake_minimum_required(VERSION 3.20)  # 3.20+ is plenty
project(GLAC VERSION 1.0
  DESCRIPTION "LQCD program for generating pure gauge configs and gradient flow"
  LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# MPI discovery
find_package(MPI 3.0 REQUIRED COMPONENTS CXX)

# JSON discovery
add_subdirectory(external/nlohmann_json EXCLUDE_FROM_ALL)
# find_package(nlohm^ann_json 3.12.0 REQUIRED)

# GLAC
add_subdirectory(${CMAKE_SOURCE_DIR}/src)

# Testing
# Catch2 pre-amble
option(CATCH_BUILD_TESTING "" OFF)
option(CATCH_BUILD_EXAMPLES "" OFF)
option(CATCH_BUILD_EXTRA_TESTS "" OFF)
option(CATCH_INSTALL_DOCS "" OFF)
option(CATCH_INSTALL_EXTRAS "" OFF)
add_subdirectory(external/Catch2)
add_subdirectory(tests)
add_subdirectory(tests_mpi)