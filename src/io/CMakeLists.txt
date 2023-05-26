set(SOURCES
    fieldio.cpp
    observablesio.cpp
)

add_library(io ${SOURCES})

target_include_directories(io PUBLIC ${CMAKE_SOURCE_DIR}/src/)