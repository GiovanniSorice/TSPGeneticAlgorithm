cmake_minimum_required(VERSION 3.11)

project(TSPGeneticAlgorithm VERSION 1.0.0 LANGUAGES CXX)

# find package needed by the project
find_package(Threads REQUIRED)

# if not specified force Relase build type
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif ()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -O0 -fsanitize=thread -fstack-check -fstack-protector-strong -fno-omit-frame-pointer -g -DDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3  -ftree-vectorize -pthread") # To view vectorization -ftree-vectorize -fopt-info-vec use-missed -note -all

set(CMAKE_CXX_STANDARD 17)

add_library(TSPGeneticAlgorithm
        src/geneticAlgorithm/geneticAlgorithm.h
        src/geneticAlgorithm/TSPGeneticAlgorithm.cpp
        src/geneticAlgorithm/TSPGeneticAlgorithm.h
        src/graph/undirectedGraph.cpp
        src/graph/undirectedGraph.h
        src/graph/graph.h
        src/geneticAlgorithm/TSPGeneticAlgorithmOMP.cpp
        src/geneticAlgorithm/TSPGeneticAlgorithmOMP.h
        src/geneticAlgorithm/TSPGeneticAlgorithmST.cpp
        src/geneticAlgorithm/TSPGeneticAlgorithmST.h src/geneticAlgorithm/TSPGeneticAlgorithmFF.cpp src/geneticAlgorithm/TSPGeneticAlgorithmFF.h)

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
