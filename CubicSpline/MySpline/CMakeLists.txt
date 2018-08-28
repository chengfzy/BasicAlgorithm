project(MySpline)

# header and source files
include_directories(./include)
aux_source_directory(./src FILE_SRCS)

# build the library
add_executable(${PROJECT_NAME} ${FILE_SRCS})
target_link_libraries(${PROJECT_NAME} ${THIRD_LIBS})