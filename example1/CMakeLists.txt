cmake_minimum_required(VERSION 3.25.0)
project(GridEngine VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(
    SOURCE_FILE_LIST
    example1.h
    main.c
    example1.c
    compress_matrix.c
)

set(SUNDIALS_DIR ${PROJECT_SOURCE_DIR}/../../sundials)
set(SUITE_SPARSE_DIR ${PROJECT_SOURCE_DIR}/../../SuiteSparse)


add_link_options(
    -static 
    -static-libgcc 
    -static-libstdc++  # Needed to avoid user having ligbcc runtime
)

add_executable(example1 ${SOURCE_FILE_LIST})

target_include_directories(
    example1 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${SUNDIALS_DIR}/include
    ${SUNDIALS_DIR}/build/include
    ${SUITE_SPARSE_DIR}/SuiteSparse_config
    ${SUITE_SPARSE_DIR}/AMD/include
    ${SUITE_SPARSE_DIR}/COLAMD/include
    ${SUITE_SPARSE_DIR}/BTF/include
    ${SUITE_SPARSE_DIR}/KLU/include
)

target_link_libraries(example1 ${SUNDIALS_DIR}/build/bin/Release/sundials_core.lib)
target_link_libraries(example1 ${SUNDIALS_DIR}/build/bin/Release/sundials_idas.lib)
target_link_libraries(example1 ${SUNDIALS_DIR}/build/bin/Release/sundials_sunlinsolklu.lib)