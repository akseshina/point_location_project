cmake_minimum_required(VERSION 3.7)
project(point_location_project)

set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-O0")
SET(CMAKE_C_FLAGS "-O0")

set(SOURCE_FILES main.cpp DCEL.h read_file.h point.h double_arithmetics.h search_structure.h geom_primitives.h)
add_executable(point_location_project ${SOURCE_FILES})