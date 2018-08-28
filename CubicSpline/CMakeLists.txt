cmake_minimum_required(VERSION 3.12)
project(CubicSpline)
set(CMAKE_CXX_STANDARD 17)


# libraries
find_package(Eigen3 REQUIRED)                           # Eigen
include_directories(${EIGEN3_INCLUDE_DIR})
set(THIRD_LIBS ${THIRD_LIBS} gflags glog)


# add example
add_subdirectory(tkSpline)
add_subdirectory(MySpline)