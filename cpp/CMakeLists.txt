cmake_minimum_required(VERSION 3.16)
project(WheelLegCpp LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 3.3 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    include(FetchContent)
    message(STATUS "Eigen3 not found; downloading Eigen 3.4.0")
    FetchContent_Declare(
        eigen
        URL https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
        DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )
    FetchContent_GetProperties(eigen)
    if(NOT eigen_POPULATED)
        FetchContent_Populate(eigen)
    endif()
    add_library(Eigen3::Eigen INTERFACE IMPORTED)
    set_target_properties(Eigen3::Eigen PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${eigen_SOURCE_DIR}"
    )
endif()

add_library(wheelleg_core
    math_utils.cpp
    model.cpp
    dYdt.cpp
    lqr.cpp
)
target_include_directories(wheelleg_core PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")
target_link_libraries(wheelleg_core PUBLIC Eigen3::Eigen)
target_compile_options(wheelleg_core PRIVATE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall;-Wextra;-Wpedantic>
)
target_compile_definitions(wheelleg_core PUBLIC
    WHEELLEG_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/.."
    WHEELLEG_CPP_DIR="${CMAKE_CURRENT_SOURCE_DIR}"
)

add_executable(wheelleg_open_loop main.cpp)
target_link_libraries(wheelleg_open_loop PRIVATE wheelleg_core)

add_executable(wheelleg_lqr lqr_control.cpp)
target_link_libraries(wheelleg_lqr PRIVATE wheelleg_core)

include(CTest)
if(BUILD_TESTING)
    add_executable(wheelleg_regression regression_test.cpp)
    target_link_libraries(wheelleg_regression PRIVATE wheelleg_core)
    add_test(NAME open_loop COMMAND wheelleg_open_loop)
    add_test(NAME lqr_closed_loop COMMAND wheelleg_lqr)
    add_test(NAME numerical_regression COMMAND wheelleg_regression)
endif()
