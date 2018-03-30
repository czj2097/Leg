cmake_minimum_required(VERSION 2.8)

project(LegTest)
include_directories(.)
add_executable(LegTest main.cpp LegKinematics.h LegKinematics.cpp)
