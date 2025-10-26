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

target_link_libraries(parallelization PUBLIC MPI::MPI_CXX)
target_include_directories(parallelization PUBLIC ${CMAKE_SOURCE_DIR}/src/)