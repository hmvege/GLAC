find_package(Catch2 3 REQUIRED)

# Add the test files
set(TEST_SOURCES
    math/su3_test.cpp
    math/su2_test.cpp
    math/complex_test.cpp
    math/random_matrix_test.cpp
    math/exp_test.cpp
    math/taylor_exp_test.cpp
    math/functions_test.cpp
)

set(TEST_UTILS
    utils/helper_math.cpp
)

# Add the test executable(These tests can use the Catch2-provided main)
add_executable(tests ${TEST_SOURCES} ${TEST_UTILS})

# Link the necessary libraries
target_link_libraries(tests PRIVATE math config Catch2::Catch2WithMain)
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/utils)