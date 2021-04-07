project(spline)

cmake_minimum_required(VERSION 3.1)

find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

include_directories(include)

set(SOURCES_BSPLINE src/bspline.cpp)
set(SOURCES_SPLINE src/spline.cpp)

add_executable(spline ${SOURCES_SPLINE})
add_executable(bspline ${SOURCES_BSPLINE})
