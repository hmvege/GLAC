set(SOURCES
    configloader.cpp
    parameters.cpp
    sysprint.cpp
)

add_library(config ${SOURCES})

target_include_directories(config PUBLIC ${CMAKE_SOURCE_DIR}/src/)

target_link_libraries(config parallelization)