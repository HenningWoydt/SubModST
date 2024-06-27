#[===[
/* SubModST solver, that solves the Cardinality-Constrained Submodular Monotone
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
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

project(SubModST CXX)
set(PROJECT_VENDOR "Henning Woydt")
set(PROJECT_CONTACT "henning.woydt@informatik.uni-heidelberg.de")
set(PROJECT_URL "https://github.com/HenningWoydt/SubModST")
set(PROJECT_DESCRIPTION "SubModST: Solving Cardinality-Constrained Submodular Monotone Subset Maximization problem")
set(PROJECT_VERSION "1.0.0")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


# Compiler Flags
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -finline-functions -fwhole-program -ftree-vectorize -flto -funroll-loops -falign-loops -march=native -Wall -Wextra -pedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -Wall -Wextra -pedantic")

# Find all source and header files
file(GLOB_RECURSE SubModST_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
file(GLOB_RECURSE SubModST_HEADERS CONFIGURE_DEPENDS "src/*.h")

# Find all test files
file(GLOB_RECURSE SubModST_Test_SOURCES CONFIGURE_DEPENDS "tests/*.cpp")
file(GLOB_RECURSE SubModST_Test_HEADERS CONFIGURE_DEPENDS "tests/*.h")

# Find BlossomV files
file(GLOB_RECURSE blossom5_files "extern/blossom5/*.cpp" "extern/blossom5/*.h")
list(FILTER blossom5_files EXCLUDE REGEX "example.cpp")


# Main executable
add_executable(submodst
        main.cpp
        ${blossom5_files}
        ${SubModST_HEADERS}
        ${SubModST_SOURCES})

add_subdirectory(external_tools/googletest)
include(gmock)
enable_testing()

# Test executable
add_executable(submodst_gtest
        ${blossom5_files}
        ${SubModST_HEADERS}
        ${SubModST_SOURCES}
        ${SubModST_Test_HEADERS}
        ${SubModST_Test_SOURCES})

# Link to Google Test
target_link_libraries(submodst_gtest gtest gtest_main)
