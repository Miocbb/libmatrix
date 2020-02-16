cmake_minimum_required(VERSION 3.8)

# compiler setting: use intel C/C++ compiler.
find_program(ICC icc)
if (ICC)
    message(STATUS "icc found: ${ICC}")
    set(CMAKE_C_COMPILER "${ICC}" CACHE FILEPATH "Intel C compiler" FORCE)
else()
    message(WARNING "icc not found.")
endif()
find_program(ICPC icpc)
if (ICPC)
    message(STATUS "icpc found: ${ICPC}")
    set(CMAKE_CXX_COMPILER "${ICPC}" CACHE FILEPATH "Intel C++ compiler" FORCE)
else()
    message(WARNING "icpc not found.")
endif()

# C/C++ standard setting
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(project_name matrix)
project(${project_name} C CXX)

# blas library
find_library(
    BLAS_LIBRARY
    NAMES blas
    PATHS "/usr/lib")
if (BLAS_LIBRARY)
    message(STATUS "libblas: ${BLAS_LIBRARY}")
else()
    message(FATAL_ERROR "libblas not found")
endif()

# lapack library
find_library(
    LAPACK_LIBRARY
    NAMES lapack
    PATHS "/usr/lib")
if (BLAS_LIBRARY)
    message(STATUS "liblapack: ${LAPACK_LIBRARY}")
else()
    message(FATAL_ERROR "liblapack not found")
endif()

file(GLOB Headers
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h")
file(GLOB Sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.c")

#if(IS_DIRECTORY $ENV{HOME}/local/googletest)
#    message(STATUS "found googletest source folder.")
#    add_subdirectory($ENV{HOME}/local/googletest $ENV{HOME}/local/googletest/build)
#else()
#    message(FATAL_ERROR "googletest source folder not found.
#        Run `git clone https://github.com/google/googletest.git` under path `$ENV{HOME}/local/`")
#endif()

# enable openmp or not.
if (OPENMP)
    message(STATUS "trying to use openmp.")
    add_compile_definitions(USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "openmp for CXX is found and enabled.")
    else()
        message(WARNING "Fail to find openmp for CXX. Will compile WITHOUT openmp for CXX.")
    endif()
else()
    message(STATUS "openmp for CXX is not used. If you want to enable openmp, use `cmake -DOPENMP=on`.")
endif()

add_library(
    ${project_name} SHARED
    ${Sources}
    ${Headers})

target_link_libraries(
    ${project_name}
    PUBLIC
    ${BLAS_LIBRARY}
    ${LAPACK_LIBRARY})

if (OPENMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(
            ${project_name}
            PUBLIC
            OpenMP::OpenMP_CXX)
    endif()
endif()

target_include_directories(
    ${project_name}
    PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

# check if gtest_main target is existed.
if (NOT TARGET gtest_main)
    enable_testing()

    # copied from google-test readme.md. See
    # https://github.com/google/googletest/blob/master/googletest/README.md.
    # {{{
    # Download and unpack googletest at configure time
    configure_file(cmake/cmake_gtest.in ${CMAKE_CURRENT_BINARY_DIR}/download/googletest/CMakeLists.txt)
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/download/googletest )
    if(result)
      message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
      RESULT_VARIABLE result
      WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/download/googletest )
    if(result)
      message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory(${CMAKE_SOURCE_DIR}/extern/googletest)
endif()

add_subdirectory(test)
