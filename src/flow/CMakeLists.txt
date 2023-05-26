set(SOURCES
    flow.cpp
)

add_library(flow ${SOURCES})

target_include_directories(flow PUBLIC ${CMAKE_SOURCE_DIR}/src/)