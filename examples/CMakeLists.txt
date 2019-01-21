cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

# disable GLFW docs, examples and tests
# see http://www.glfw.org/docs/latest/build_guide.html
set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)

add_subdirectory(submodules/glfw)
include_directories(..)
include_directories(submodules)
include_directories(submodules/glfw/include)

project (dj_fft)

add_executable(img_to_fft2d img_to_fft2d.cpp)

add_executable(gpu_acceleration gpu_acceleration.cpp submodules/glad/glad.c)
target_link_libraries(gpu_acceleration glfw)
