project(tkSpline)

# header and source files
aux_source_directory(./ FILE_SRCS)

# build the library
add_executable(${PROJECT_NAME} ${FILE_SRCS})
target_link_libraries(${PROJECT_NAME} ${THIRD_LIBS})