cmake_minimum_required(VERSION 3.20)
project(DipoleCell)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-march=native -mtune=native)
message(STATUS "Using -march=native -mtune=native")

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
        message(STATUS "Using STL debug mode on Clang")
        add_compile_definitions(_LIBCPP_DEBUG=1)
    elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        message(STATUS "Using STL debug mode on GCC")
        add_compile_definitions(_GLIBCXX_DEBUG)
    else()
        message(AUTHOR_WARNING "Building Debug but not using STL debug mode since compiler not known to support it")
    endif()

    # elseif (CMAKE_BUILD_TYPE MATCHES "Release")
endif()

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(TBB REQUIRED)
find_package(Boost REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)
find_package(OpenMP REQUIRED)

set(common_source_files_all_targets
    CommonTypes.cpp CommonTypes.h
    RandomHelper.cpp RandomHelper.h pcg-cpp-0.98/pcg_random.hpp pcg-cpp-0.98/pcg_uint128.hpp pcg-cpp-0.98/pcg_extras.hpp
    SimulationBox.cpp SimulationBox.h
    ModelCellParameter.h CellBase.cpp CellBase.h Cell.cpp Cell.h
    ModelTissueParameter.h Tissue.cpp Tissue.h Initialization.cpp Initialization.h
    InitialCondition.cpp InitialCondition.h)
set(common_linked_libraries_all_targets TBB::tbb OpenMP::OpenMP_CXX Eigen3::Eigen Boost::boost)

add_executable(DipoleCell main.cpp PostProcessingPybind11.cpp PostProcessingPybind11.h ${common_source_files_all_targets})
target_link_libraries(DipoleCell PUBLIC ${common_linked_libraries_all_targets} pybind11::embed)

# add_custom_command(TARGET DipoleCell PRE_BUILD
# COMMAND ${CMAKE_COMMAND} -E copy
# ${CMAKE_SOURCE_DIR}/PostProcessing.py
# ${CMAKE_BINARY_DIR})
configure_file(${CMAKE_SOURCE_DIR}/PostProcessing.py ${CMAKE_BINARY_DIR}/PostProcessing.py COPYONLY)

# pybind11_add_module(DipoleCellPybind11 TissuePybind11.cpp ${common_source_files_all_targets})
# target_link_libraries(DipoleCellPybind11 PUBLIC ${common_linked_libraries_all_targets})
