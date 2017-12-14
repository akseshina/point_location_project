cmake_minimum_required(VERSION 3.7)
project(point_location_project)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES main.cpp DCEL.h read_file.h)
add_executable(point_location_project ${SOURCE_FILES})