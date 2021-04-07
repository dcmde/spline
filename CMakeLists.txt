project(spline)

cmake_minimum_required(VERSION 3.1)

include_directories(include)

set(SOURCES_BSPLINE src/bspline.cpp)
set(SOURCES_SPLINE src/spline.cpp)

add_executable(bspline ${SOURCES_BSPLINE})
add_executable(spline ${SOURCES_SPLINE})
