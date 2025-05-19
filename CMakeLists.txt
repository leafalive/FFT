cmake_minimum_required(VERSION 3.28)
project(FFT)

set(CMAKE_CXX_STANDARD 14)
set(ArrayFire_DIR "D:/Program Files/ArrayFire/v3/cmake")

find_package(Eigen3 CONFIG REQUIRED)
find_package(ArrayFire REQUIRED)

# The Unified backend lets you choose the backend at runtime.
# To use the Unified backend, do the following:
add_executable(FFT main.cpp
        fft.cpp
        fft.cpp
        fft.h
        fft.h)
target_link_libraries(FFT ArrayFire::af Eigen3::Eigen)

