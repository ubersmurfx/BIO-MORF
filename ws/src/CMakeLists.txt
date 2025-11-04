cmake_minimum_required(VERSION 3.10)
project(FishSimulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Qt5 REQUIRED COMPONENTS Core Widgets)

# QCustomPlot
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/qcustomplot)

add_executable(fish_simulation main.cpp qcustomplot.cpp)

target_link_libraries(fish_simulation Qt5::Core Qt5::Widgets)