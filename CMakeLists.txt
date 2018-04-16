cmake_minimum_required(VERSION 2.8)

project(LegTest)
include_directories(.)
set(CMAKE_CXX_FLAGS "-std=c++11")
#add_compile_options(-std=c++11)
set(LIBRARY_OUTPUT_PATH lib)
add_library(time_optimal TimeOptimalMotion.cpp)

add_executable(LegTest main.cpp)
target_link_libraries(LegTest time_optimal)

set(CMAKE_INSTALL_PREFIX "/usr/TimeOptimal")
install(TARGETS time_optimal ARCHIVE DESTINATION lib)
install(FILES TimeOptimalMotion.h DESTINATION include)

add_custom_target(uninstall rm -r /usr/TimeOptimal)
