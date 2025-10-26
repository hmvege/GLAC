set(SOURCES
    fieldio.cpp
    observablesio.cpp
)

add_library(io ${SOURCES})
target_link_libraries(io PUBLIC MPI::MPI_CXX)
target_include_directories(io PUBLIC ${CMAKE_SOURCE_DIR}/src/)