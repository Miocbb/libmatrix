cmake_minimum_required(VERSION 3.8)

set(project_name matrix_test)
project(${project_name} C CXX)

find_package (Eigen3 REQUIRED NO_MODULE)
if(Eigen3_FOUND)
    message(STATUS "Eigen3 library found.")
else()
    message(WARNING "Can not found library Eigen3, which will be used for testing.\
        Test code will be skipped until Eigen3 is installed.")
    return()
endif()

file(GLOB Sources
    "*.cpp"
    "*.cc")
file(GLOB Headers
    "*.hpp"
    "*.h")

add_executable(
    ${project_name}
    ${Sources}
    ${Headers})

target_link_libraries(
    ${project_name} PUBLIC
    matrix
    gtest_main
    Eigen3::Eigen)

add_test(
    NAME ${project_name}
    COMMAND ${project_name})
