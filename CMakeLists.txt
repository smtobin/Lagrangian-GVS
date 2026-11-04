cmake_minimum_required(VERSION 3.18)
project(LagrangianGVS LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)                       # specify C++ standard
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -march=native -g -fno-omit-frame-pointer")
set(CMAKE_BUILD_TYPE Release)

find_package(Eigen3 3.3 REQUIRED)   # find Eigen

# find VTK
find_package(VTK COMPONENTS
    CommonColor
    CommonCore
    FiltersSources
    IOImage
    IOGeometry
    InteractionStyle
    RenderingContextOpenGL2
    RenderingCore
    RenderingFreeType
    RenderingGL2PSOpenGL2
    RenderingOpenGL2
)

if (NOT VTK_FOUND)
  message(FATAL_ERROR "Unable to find the VTK build folder.")
endif()

add_executable(Run main.cpp)
target_link_libraries(Run Eigen3::Eigen ${VTK_LIBRARIES})

vtk_module_autoinit(
  TARGETS Run
  MODULES ${VTK_LIBRARIES}
)