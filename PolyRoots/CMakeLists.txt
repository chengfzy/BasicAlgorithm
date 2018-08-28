project(PolyRoots)
set(CMAKE_CXX_STANDARD 17)


# libraries
find_package(Eigen3 REQUIRED)                           # Eigen
include_directories(${EIGEN3_INCLUDE_DIR})
set(THIRD_LIBS ${THIRD_LIBS} gflags glog)


# header and source files
include_directories(./include)
aux_source_directory(./src FILE_SRCS)

# build the library
add_executable(${PROJECT_NAME} ${FILE_SRCS})
target_link_libraries(${PROJECT_NAME} ${THIRD_LIBS})