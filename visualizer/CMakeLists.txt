cmake_minimum_required(VERSION 3.16)

project(visualization_tool LANGUAGES CXX)

# Dependencies
find_package(Boost REQUIRED COMPONENTS program_options)

# Executable
add_executable(visualizer
    src/main.cpp
)

# Link libraries
target_link_libraries(visualizer
    PRIVATE
        diskpack
        Boost::program_options
)

# C++ standard
target_compile_features(visualizer PRIVATE cxx_std_17)