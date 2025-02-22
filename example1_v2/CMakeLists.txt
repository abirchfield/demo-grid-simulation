cmake_minimum_required(VERSION 3.25.0)
project(Examples VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(
    SOURCE_FILE_LIST
    main.cpp
    example1_v2.cpp
)

#set(GITHUB_DIR ${PROJECT_SOURCE_DIR}/../..)
set(GITHUB_DIR "C:/Users/abirchfield/Documents/GitHub")
set(SUNDIALS_DIR ${GITHUB_DIR}/sundials)
set(SUITE_SPARSE_DIR ${GITHUB_DIR}/SuiteSparse)

add_link_options(
    -static 
    -static-libgcc 
    -static-libstdc++  # Needed to avoid user having ligbcc runtime
)

add_executable(example1_v2 ${SOURCE_FILE_LIST})

target_include_directories(
    example1_v2 PRIVATE
    ${CMAKE_SOURCE_DIR}
    ${SUNDIALS_DIR}/include
    ${SUNDIALS_DIR}/build/include
    ${SUITE_SPARSE_DIR}/SuiteSparse_config
    ${SUITE_SPARSE_DIR}/AMD/include
    ${SUITE_SPARSE_DIR}/COLAMD/include
    ${SUITE_SPARSE_DIR}/BTF/include
    ${SUITE_SPARSE_DIR}/KLU/include
)

target_link_libraries(example1_v2 ${SUNDIALS_DIR}/build/bin/Release/sundials_core.lib)
target_link_libraries(example1_v2 ${SUNDIALS_DIR}/build/bin/Release/sundials_idas.lib)
target_link_libraries(example1_v2 ${SUNDIALS_DIR}/build/bin/Release/sundials_sunlinsolklu.lib)