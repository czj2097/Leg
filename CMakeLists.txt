cmake_minimum_required(VERSION 2.8)

project(LegTest)
include_directories(.)
add_compile_options(-std=c++11)
add_executable(LegTest main.cpp LegKinematics.h LegKinematics.cpp TimeOptimalMotion.cpp)
