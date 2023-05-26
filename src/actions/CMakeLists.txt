set(SOURCES
    action.cpp
    wilsonexplicitder.cpp
    wilsongaugeaction.cpp
)

add_library(actions ${SOURCES})

target_include_directories(actions PUBLIC ${CMAKE_SOURCE_DIR}/src/)

target_link_libraries(actions math parallelization config)