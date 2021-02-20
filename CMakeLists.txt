cmake_minimum_required(VERSION 3.17)
project(LAB2)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)
include_directories(.idea)

add_executable(LAB2
        .idea/.gitignore
        .idea/workspace.xml
        FORSYTHE.h
        main.cpp
        MATRIX.CPP
        MATRIX.H)
