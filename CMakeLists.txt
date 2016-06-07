cmake_minimum_required(VERSION 2.8)

project(pb_solvers)

set( CMAKE_VERBOSE_MAKEFILE false )
set(PBSOLV_ROOT ${PROJECT_SOURCE_DIR})
set(PBSOLV_BUILD ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PBSOLV_BUILD}/bin)
set(EXECUTABLE_INSTALL_PATH bin)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -funroll-loops")
elseif(COMPILER_SUPPORTS_CXX0X)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x -funroll-loops")
else()
  message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()


enable_testing()

add_subdirectory(pb_shared/src)
add_subdirectory(pb_shared/unittest)

add_subdirectory(pbam/src)
add_subdirectory(pbam/unittest)

add_subdirectory(pbsam/src)
add_subdirectory(pbsam/unittest)

file(MAKE_DIRECTORY build/bin)
