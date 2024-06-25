#[===[
/* CCSMSM solver, that solves the Cardinality-Constrained Submodular Monotone
   Subset Maximization problem.
   Copyright (C) 2024  Henning Woydt

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or any
   later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
==============================================================================*/
]] <- this is not the closing
]===]

cmake_minimum_required(VERSION 3.24)
project(CCSMSM)

set(CMAKE_CXX_STANDARD 17)

# Google Test
include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/03597a01ee50ed33e9dfd640b249b4be3799d395.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)
enable_testing()

# Compiler Flags
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++17 -O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-std=c++17 -g3 -O0 -Wall -Wextra -pedantic")

file(GLOB blossom5_files "src/3rd_party/blossom5/*.cpp" "src/3rd_party/blossom5/*.h" "src/3rd_party/blossom5/GEOM/*.cpp" "src/3rd_party/blossom5/GEOM/*.h" "src/3rd_party/blossom5/MinCost/*.cpp" "src/3rd_party/blossom5/MinCost/*.h")
list(FILTER blossom5_files EXCLUDE REGEX "example.cpp")

# Main executable
add_executable(CCSMSM
        main.cpp
        ${blossom5_files}
        src/structures/StructureInterface.h
        src/structures/Graph.h
        src/structures/Matrix.h
        src/structures/Clustering.h
        src/structures/FacilityLocation.h
        src/structures/WeightedCoverage.h
        src/structures/BipartiteInfluence.h
        src/algorithms/BFSolver.h
        src/algorithms/Solver.h
        src/algorithms/CandidateManager.h
        src/algorithms/Heap.h
        src/algorithms/Cache.h
        src/algorithms/PairwiseUB.h
        src/algorithms/StatisticProfiler.h
        src/util/utility.h
        src/util/utility.cpp
        src/util/AlgorithmConfiguration.h
        src/util/CommandLineParser.h)

# Test executable
add_executable(CCSMSM_Test
        ${blossom5_files}
        src/structures/StructureInterface.h
        src/structures/Graph.h
        src/structures/Matrix.h
        src/structures/Clustering.h
        src/structures/FacilityLocation.h
        src/structures/WeightedCoverage.h
        src/structures/BipartiteInfluence.h
        src/algorithms/BFSolver.h
        src/algorithms/Solver.h
        src/algorithms/CandidateManager.h
        src/algorithms/Heap.h
        src/algorithms/Cache.h
        src/algorithms/PairwiseUB.h
        src/algorithms/StatisticProfiler.h
        src/util/utility.h
        src/util/utility.cpp
        src/util/AlgorithmConfiguration.h
        src/util/CommandLineParser.h
        tests/BruteForceTests.cpp
        tests/PartialDominatingSetTests.cpp
        tests/GroupClosenessCentralityTests.cpp
        tests/EuclidianClusteringTests.cpp
        tests/FacilityLocationTests.cpp
        tests/ManhattanClusteringTests.cpp
        tests/PartialVertexSetTests.cpp
        tests/WeightedCoverageTests.cpp
        tests/BipartiteInfluenceTests.cpp
        tests/MiscTest.cpp)

# IPO Support
include(CheckIPOSupported)
check_ipo_supported(RESULT supported OUTPUT error)
if (supported)
    message(STATUS "IPO / LTO enabled")
    set_property(TARGET CCSMSM PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
    set_property(TARGET CCSMSM_Test PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
else ()
    message(STATUS "IPO / LTO not supported: <${error}>")
endif ()

# Link to Google Test
target_link_libraries(CCSMSM_Test GTest::gtest_main)
include(GoogleTest)
gtest_discover_tests(CCSMSM_Test)
