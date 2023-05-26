# List the source files
set(SOURCES
    communicator.cpp
    index.cpp
    neighbourlist.cpp
    neighbours.cpp
    parallelparameters.cpp
)

# Add the executable
add_library(parallelization ${SOURCES})

message(STATUS "CMAKE_SOURCE_DIR: ${CMAKE_SOURCE_DIR}")

target_include_directories(parallelization PUBLIC ${CMAKE_SOURCE_DIR}/src/)