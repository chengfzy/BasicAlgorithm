cmake_minimum_required(VERSION 3.12)
project(BasicAlgorithm)
set(CMAKE_CXX_STANDARD 17)


# projects build setting
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# sub-project
add_subdirectory(HammingDistance)
add_subdirectory(CubicSpline)
add_subdirectory(PolyRoots)
