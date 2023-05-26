# List the source files
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

# # Add the executable
# add_executable(tests ${SOURCES})
add_library(tests ${SOURCES})

# # Link the libraries to the executable
# target_link_libraries(tests actions config flow io math observables parallelization)