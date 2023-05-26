set(SOURCES
    correlator.cpp
    energydensity.cpp
    latticeactionchargedensity.cpp
    mastersampler.cpp
    mastersamplertopcxyz.cpp
    plaquette.cpp
    supersampler.cpp
    topologicalcharge.cpp
    tools/observablestorer.cpp
)

add_library(observables ${SOURCES})

target_include_directories(observables PUBLIC ${CMAKE_SOURCE_DIR}/src/)