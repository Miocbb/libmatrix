cmake_minimum_required(VERSION 3.8)

# need blas library
find_library(
    BLAS_LIBRARY
    NAMES blas
    PATHS "/usr/lib")
if (BLAS_LIBRARY)
    message(STATUS "libblas: ${BLAS_LIBRARY}")
else()
    message(FATAL_ERROR "libblas not found")
endif()

# need lapack library
find_library(
    LAPACK_LIBRARY
    NAMES lapack
    PATHS "/usr/lib")
if (BLAS_LIBRARY)
    message(STATUS "liblapack: ${LAPACK_LIBRARY}")
else()
    message(FATAL_ERROR "liblapack not found")
endif()

# ==> build matrix <==
file(GLOB SOURCES "*.cpp")
add_library(
    ${PROJECT_MATRIX} SHARED
    ${SOURCES})

target_link_libraries(
    ${PROJECT_MATRIX}
    PUBLIC
    ${BLAS_LIBRARY}
    ${LAPACK_LIBRARY})

if (OPENMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(
            ${PROJECT_MATRIX}
            PUBLIC
            OpenMP::OpenMP_CXX)
    endif()
endif()

set(MATRIX_PUBLIC_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/include")
target_include_directories(
    ${PROJECT_MATRIX}
    PUBLIC
    ${MATRIX_PUBLIC_INCLUDE_DIR})

# ==> install matrix <==
install(TARGETS ${PROJECT_MATRIX} DESTINATION lib)
install(DIRECTORY "${MATRIX_PUBLIC_INCLUDE_DIR}/matrix" DESTINATION include)
