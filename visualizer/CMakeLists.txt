cmake_minimum_required(VERSION 3.16)
project(visualizer LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(visualizer src/main.cpp)

target_link_libraries(visualizer PRIVATE
    diskpack
    Boost::program_options
)
