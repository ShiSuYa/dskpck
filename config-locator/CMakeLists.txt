cmake_minimum_required(VERSION 3.16)

project(config_locator LANGUAGES CXX)

# Dependencies
find_package(Boost REQUIRED COMPONENTS program_options)

# Executable
add_executable(config-locator
    src/main.cpp
)

# Link libraries
target_link_libraries(config-locator
    PRIVATE
        diskpack
        Boost::program_options
)

# C++ standard
target_compile_features(config-locator PRIVATE cxx_std_17)