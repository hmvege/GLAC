find_package(Catch2 3 REQUIRED)

# Add the test files
set(TEST_SOURCES
    math/su3_test.cpp
    # config/config_test.cpp
)

set(TEST_UTILS
    utils/helper_comparers.h
    utils/helper_generators.h
    utils/helper_math.h
)

# Add the test executable(These tests can use the Catch2-provided main)
add_executable(tests ${TEST_SOURCES})

# Link the necessary libraries
target_link_libraries(tests PRIVATE math Catch2::Catch2WithMain)
