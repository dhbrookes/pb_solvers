cmake_minimum_required(VERSION 3.12)
project(pb_solvers)
set(CMAKE_VERBOSE_MAKEFILE false)

################################################
################################################

# ---[ Using cmake scripts and modules
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
include(ExternalProject)
include(cmake/Utils.cmake)
include(cmake/Dependencies.cmake)
include(cmake/Summary.cmake)

set(PBSOLV_ROOT ${PROJECT_SOURCE_DIR})
set(PBSOLV_BUILD ${PROJECT_BINARY_DIR})
set(EXECUTABLE_OUTPUT_PATH ${PBSOLV_BUILD}/bin)
set(EXECUTABLE_INSTALL_PATH bin)

option(ENABLE_PBAM "Main flag for enabling PBAM" OFF)
option(ENABLE_PBSAM "Main flag for enabling PBSAM" OFF)
option(ENABLE_PB_TESTING "Testing in pb_solvers" OFF)
option(ENABLE_PBAM_APBS "pbam APBS" OFF)
option(ENABLE_PBSAM_APBS "pbsam APBS" OFF)


################################################
##### Need a compiler for C++11
################################################

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

################################################
### Compiler flags for optimization
################################################
if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3") ## Optimize
endif()

################################################
###### OMP paralleization
################################################

option(ENABLE_OPENMP "Enable OpenMP parallelism" OFF)

if (ENABLE_OPENMP)
  find_package(OpenMP)
  if (OPENMP_FOUND)
    add_definitions(-D__OMP)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  endif()
endif()

################################################
################################################
### For APBS Sphinx
################################################
option(ENABLE_PBAM_SPHINX "pbam SPHINX" OFF)
if(ENABLE_PBAM_SPHINX)

  #if(APPLE)
  #  list(APPEND CMAKE_MODULE_LINKER_FLAGS "-Wl,-undefined,dynamic_lookup")
  #endif()
  #if(WIN32)
  #  list(APPEND CMAKE_MODULE_LINKER_FLAGS "/LIBPATH:C:\Python34\libs")
  #endif()

  #list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/../../../cmake")
  find_package(Python3 3.4 REQUIRED COMPONENTS Interpreter Development)
  #find_package(Cython REQUIRED)
 
  include_directories(${Python3_INCLUDE_DIRS})
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pbam_py.c"
                     COMMAND "${PYTHON_EXECUTABLE}" -m cython
                     -o "${CMAKE_CURRENT_BINARY_DIR}/pbam_py.c"
                     "${CMAKE_CURRENT_SOURCE_DIR}/pbam_sph.pyx"
                     MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/pbam_sph.pyx")
endif()


################################################
###### Testing
################################################
if(ENABLE_PB_TESTING)
  find_package(GTest REQUIRED)
  find_package(Threads REQUIRED)


endif()


################################################
###### Sub directories for all of source
################################################
if(ENABLE_PB_TESTING)
  add_subdirectory(pb_shared/unittest)
endif()

if(ENABLE_PBAM)
  add_subdirectory(pbam/src)
  if(ENABLE_PB_TESTING)
    add_subdirectory(pbam/unittest)
  endif()
endif()

if(ENABLE_PBSAM)
  add_subdirectory(pbsam/src)
  if(ENABLE_PB_TESTING)
    add_subdirectory(pbsam/unittest)
  endif()
endif()

################################################
##### Make a directory for the executable
################################################
#file(MAKE_DIRECTORY build/bin)

pbsam_print_configuration_summary()
