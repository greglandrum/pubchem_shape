cmake_minimum_required(VERSION 3.10)

project(pubchem_shape)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)


Include(FetchContent)

FetchContent_Declare(
  Catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG        v3.4.0 # or a later release
)

FetchContent_MakeAvailable(Catch2)

find_package(RDKit REQUIRED)
add_library(pubchem_align3d SHARED shape_functions1.cpp shape_functions2.cpp shape_neighbor.cpp)


add_library(rdkit_shape SHARED rdkit_shape.cpp)
target_link_libraries(rdkit_shape pubchem_align3d RDKit::SmilesParse RDKit::SubstructMatch)
target_include_directories(rdkit_shape PUBLIC ${RDKIT_INCLUDE_DIR} ${Boost_INCLUDE_DIR})

add_subdirectory(Wrap)

add_executable(rdkit_test rdkit_align3d_test.cpp )
target_link_libraries(rdkit_test pubchem_align3d RDKit::FileParsers)

add_executable(sdf_align sdf_align.cpp )
target_link_libraries(sdf_align rdkit_shape pubchem_align3d RDKit::SubstructMatch RDKit::FileParsers)

add_executable(tests test.cpp)
target_link_libraries(tests PRIVATE Catch2::Catch2WithMain rdkit_shape RDKit::FileParsers )


