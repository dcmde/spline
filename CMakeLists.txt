project(spline)

cmake_minimum_required(VERSION 3.1)

include_directories(include)

set(SOURCES src/main.cpp)

add_executable(${PROJECT_NAME} ${SOURCES})
