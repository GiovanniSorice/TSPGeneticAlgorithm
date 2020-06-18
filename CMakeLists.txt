cmake_minimum_required(VERSION 3.16)

project(TSPGeneticAlgorithm VERSION 1.0.0 LANGUAGES CXX)

# find package needed by the project
find_package(Threads REQUIRED)

# if not specified force Relase build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -O0 -fstack-check -fstack-protector-strong -fno-omit-frame-pointer -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3") # To view vectorization -ftree-vectorize -fopt-info-vec

set(CMAKE_CXX_STANDARD 17)

add_library(TSPGeneticAlgorithm
            geneticAlgorithm.h
            TSPGeneticAlgorithm.cpp
            TSPGeneticAlgorithm.h
            undirectedGraph.cpp
            undirectedGraph.h
            graph.h)

target_link_libraries(TSPGeneticAlgorithm
        Threads::Threads
        )

set_target_properties(TSPGeneticAlgorithm PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)

# create main target
add_executable(main main.cpp)
set_target_properties(main PROPERTIES
        CXX_STANDARD 17
        LINKER_LANGUAGE CXX)
target_link_libraries(main TSPGeneticAlgorithm)
