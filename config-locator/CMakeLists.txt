cmake_minimum_required(VERSION 3.16)
project(config_locator VERSION 1.0 LANGUAGES CXX)

find_package(Boost REQUIRED COMPONENTS program_options)

add_executable(config_locator
    src/main.cpp
)

target_link_libraries(config_locator PRIVATE
    diskpack
    Boost::program_options
)