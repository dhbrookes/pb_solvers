cmake_minimum_required(VERSION 3.0)

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

add_executable(pbsam
	pbsam/BesselCalc.cpp
	pbsam/BesselCalc.h
	pbsam/Constants.cpp
	pbsam/Constants.h
	pbsam/main.cpp
	pbsam/MyMatrix.h
	pbsam/readutil.h
	pbsam/setup.cpp
	pbsam/setup.h
	pbsam/SHCalc.cpp
	pbsam/SHCalc.h
	pbsam/System.cpp
	pbsam/System.h
	pbsam/util.h
	)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

install(TARGETS pbsam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
