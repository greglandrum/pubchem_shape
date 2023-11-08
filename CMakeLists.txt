cmake_minimum_required(VERSION 3.10)

project(pubchem_shape)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


find_package(RDKit REQUIRED)
add_library(pubchem_align3d SHARED shape_functions1.cpp shape_functions2.cpp shape_neighbor.cpp)

add_executable(rdkit_test rdkit_align3d_test.cpp )
target_link_libraries(rdkit_test pubchem_align3d RDKit::FileParsers)