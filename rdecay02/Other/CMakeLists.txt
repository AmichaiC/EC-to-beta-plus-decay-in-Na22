cmake_minimum_required(VERSION 3.10)

project(plot_gamma_histogram)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find ROOT package
find_package(ROOT REQUIRED)

include(${ROOT_USE_FILE})

add_executable(plot_gamma_histogram plot_gamma_histogram.cpp)

# Link the ROOT libraries
target_link_libraries(plot_gamma_histogram ${ROOT_LIBRARIES})
