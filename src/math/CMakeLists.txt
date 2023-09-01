set(SOURCES
    complex.cpp
    exponentiation/expluscher.cpp
    exponentiation/su3exp.cpp
    exponentiation/taylorexp.cpp
    matrices/su2.cpp
    matrices/su3.cpp
    matrices/su3matrixgenerator.cpp
)

add_library(math ${SOURCES})

target_include_directories(math PUBLIC ${CMAKE_SOURCE_DIR}/src/)