cmake_minimum_required (VERSION 3.1)
set (CMAKE_CXX_STANDARD 11)

include_directories(..)

project (dj_fft)
add_executable(img_to_fft2d img_to_fft2d.cpp)
