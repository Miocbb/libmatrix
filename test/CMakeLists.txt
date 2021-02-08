cmake_minimum_required(VERSION 3.8)

set(PROJECT_MATRIX_TEST matrix_test)
project(${PROJECT_MATRIX_TEST} C CXX)

find_package (Eigen3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 library found.")
else()
    message(FETAL_ERROR "Can not found library Eigen3, which will be used for testing.\
        Test code will be skipped until Eigen3 is installed.")
    return()
endif()

file(GLOB PROJECT_SOURCES
    "*.cpp"
    "*.cc")

add_executable(
    ${PROJECT_MATRIX_TEST}
    ${PROJECT_SOURCES})

target_link_libraries(
    ${PROJECT_MATRIX_TEST}
    PUBLIC
    ${PROJECT_MATRIX}
    gtest_main
    Eigen3::Eigen)

add_test(
    NAME ${PROJECT_MATRIX_TEST}
    COMMAND ${PROJECT_MATRIX_TEST})
