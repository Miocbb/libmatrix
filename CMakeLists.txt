cmake_minimum_required(VERSION 3.8)

# set project matrix and its name
set(PROJECT_MATRIX matrix)
project(${PROJECT_MATRIX} C CXX)

# C/C++ standard setting
set(CMAKE_C_STANDARD 99)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# ==> build matrix <==
add_subdirectory("${CMAKE_SOURCE_DIR}/src/matrix")

# ==> build tests <==
option(BUILD_TEST "build test or not." OFF)
if (BUILD_TEST)
    enable_testing()
    # Download googletest library
    set(GOOGLE_TEST_SOURCE_DIR "${CMAKE_SOURCE_DIR}/external/googletest")
    if(NOT EXISTS "${GOOGLE_TEST_SOURCE_DIR}")
        message(STATUS "Download external library: googletest.")
        configure_file(cmake/cmake_gtest.in download/googletest/CMakeLists.txt)
        execute_process(COMMAND "${CMAKE_COMMAND}" -G "${CMAKE_GENERATOR}" .
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/googletest"
        )
        execute_process(COMMAND "${CMAKE_COMMAND}" --build .
            WORKING_DIRECTORY "${CMAKE_BINARY_DIR}/download/googletest"
        )
    endif()

    # ==> build google test <==
    # Prevent overriding the parent project's compiler/linker
    # settings on Windows
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    # Add googletest directly to our build. This defines
    # the gtest and gtest_main targets.
    add_subdirectory("${GOOGLE_TEST_SOURCE_DIR}")

    # ==> build test <==
    add_subdirectory(test)
else()
    message(STATUS "Note: Tests are disabled. Use \"-DBUILD_TEST=On\" to build tests.")
endif()
