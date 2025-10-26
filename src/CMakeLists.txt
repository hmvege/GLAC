# Adding GLAC libraries
add_subdirectory(math)
add_subdirectory(actions)
add_subdirectory(config)
add_subdirectory(flow)
add_subdirectory(io)
add_subdirectory(observables)
add_subdirectory(parallelization)
add_subdirectory(tests_old)

# Executable
add_executable(GLAC
  main.cpp
  system.cpp)

target_compile_features(GLAC PRIVATE cxx_std_17)
# if(CMAKE_BUILD_TYPE STREQUAL "Release")
#   target_compile_options(GLAC PRIVATE -O3 -march=native)
# endif()
# target_compile_options(GLAC PRIVATE -Werror -Wall -Wextra)

# Link MPI the right way
target_link_libraries(GLAC PRIVATE MPI::MPI_CXX
  config flow io math observables parallelization actions tests_old)

set_target_properties(GLAC PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
