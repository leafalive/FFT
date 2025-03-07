cmake_minimum_required(VERSION 3.28)
project(FFT)

set(CMAKE_CXX_STANDARD 14)

find_package(Eigen3 CONFIG REQUIRED)
add_executable(FFT main.cpp
        fft.cpp
        fft.cpp
        fft.h
        fft.h)

target_link_libraries(FFT PRIVATE Eigen3::Eigen)

