cmake_minimum_required(VERSION 3.16)

project(diskpack_project LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Output directory
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Optimization flags
if(CMAKE_BUILD_TYPE STREQUAL "Release")
    add_compile_options(-O3)
else()
    add_compile_options(-O0 -g)
endif()

# Submodules
add_subdirectory(diskpack)
add_subdirectory(visualizer)
add_subdirectory(config-locator)