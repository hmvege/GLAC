set(SOURCES
    performancetests.cpp
    testsuite.cpp
    unit_tests/actiontests.cpp
    unit_tests/complexoperations.cpp
    unit_tests/iotests.cpp
    unit_tests/latticeoperations.cpp
    unit_tests/observabletests.cpp
    unit_tests/randommatrixtests.cpp
    unit_tests/su2operations.cpp
    unit_tests/su3exptest.cpp
    unit_tests/su3operations.cpp
    unit_tests/testcore.cpp
)

add_library(tests ${SOURCES})

target_include_directories(tests PUBLIC ${CMAKE_SOURCE_DIR}/src/)

target_link_libraries(tests flow io math observables parallelization actions)