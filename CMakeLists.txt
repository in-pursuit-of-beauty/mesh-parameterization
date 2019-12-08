cmake_minimum_required(VERSION 3.1)
project(add-texcoords)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

option(LIBIGL_USE_STATIC_LIBRARY  "Use libigl as static library"  OFF)
option(LIBIGL_WITH_OPENGL         "Use OpenGL"                    OFF)
option(LIBIGL_WITH_OPENGL_GLFW    "Use GLFW"                      OFF)

find_package(LIBIGL REQUIRED QUIET)

# Add the project files
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/")
file(GLOB SRCFILES src/*.cpp)

add_executable(${PROJECT_NAME} "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp" ${SRCFILES})
target_link_libraries(${PROJECT_NAME} igl::core)
