cmake_minimum_required(VERSION 3.18)
project(LagrangianGVS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)                       # specify C++ standard
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -march=native -g -fno-omit-frame-pointer")
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 3.3 REQUIRED)   # find Eigen

add_executable(Run main.cpp)
target_link_libraries(Run Eigen3::Eigen)