set(SOURCES
    flow.cpp
)

add_library(flow ${SOURCES})

target_link_libraries(flow PUBLIC MPI::MPI_CXX)
target_include_directories(flow PUBLIC ${CMAKE_SOURCE_DIR}/src/)