find_package(Threads REQUIRED)

# Enable ExternalProject CMake module
include(ExternalProject)

# Download and install GoogleTest
ExternalProject_Add(
    gtest
    URL https://googletest.googlecode.com/files/gtest-1.7.0.zip
    PREFIX ${CMAKE_CURRENT_BINARY_DIR}/gtest
    # Disable install step
    INSTALL_COMMAND ""
)

# Create a libgtest target to be used as a dependency by test programs
add_library(libgtest IMPORTED STATIC GLOBAL)
add_dependencies(libgtest gtest)

# Set gtest properties
ExternalProject_Get_Property(gtest source_dir binary_dir)
set_target_properties(libgtest PROPERTIES
    "IMPORTED_LOCATION" "${binary_dir}/libgtest.a"
    "IMPORTED_LINK_INTERFACE_LIBRARIES" "${CMAKE_THREAD_LIBS_INIT}"
)

include_directories("${source_dir}/include")

file(GLOB SRCS *.cpp)

include_directories("../src")

ADD_EXECUTABLE(pbamtest 
  ../src/ASolver.cpp
  ../src/BD.cpp
  ../src/BesselCalc.cpp
  ../src/Constants.cpp
  ../src/Electrostatics.cpp
  ../src/EnergyForce.cpp
  ../src/ReExpCalc.cpp
  ../src/SHCalc.cpp
  ../src/System.cpp
  ../src/setup.cpp
  main.cpp
  )
  
TARGET_LINK_LIBRARIES(pbamtest
    libgtest
)

add_test(NAME pbamtest
         COMMAND pbamtest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
  
install(TARGETS pbamtest RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
