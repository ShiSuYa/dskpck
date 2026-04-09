cmake_minimum_required(VERSION 3.16)

project(diskpack LANGUAGES CXX)

find_package(Boost REQUIRED)

# Library
add_library(diskpack
    src/geometry.cpp
    src/corona.cpp
    src/generator.cpp
    src/codec.cpp
    src/search.cpp
)

# Include directories
target_include_directories(diskpack
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${Boost_INCLUDE_DIRS}
)

target_link_libraries(diskpack PUBLIC
        Boost::headers
)

# Compile options
target_compile_features(diskpack PUBLIC cxx_std_17)