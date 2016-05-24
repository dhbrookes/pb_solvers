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
#    "INTERFACE_INCLUDE_DIRECTORIES" "${source_dir}/include"
)
# I couldn't make it work with INTERFACE_INCLUDE_DIRECTORIES
include_directories("${source_dir}/include")

file(GLOB SRCS *.cpp)

include_directories("../pbam")

ADD_EXECUTABLE(pbamtest 
  ../pbam/ASolver.cpp
  ../pbam/ASolver.h
  ../pbam/BD.cpp
  ../pbam/BD.h
  ../pbam/BesselCalc.cpp
  ../pbam/BesselCalc.h
  ../pbam/CMakeLists.txt
  ../pbam/Constants.cpp
  ../pbam/Constants.h
  ../pbam/Electrostatics.cpp
  ../pbam/Electrostatics.h
  ../pbam/EnergyForce.cpp
  ../pbam/EnergyForce.h
  ../pbam/MyMatrix.h
  ../pbam/ReExpCalc.cpp
  ../pbam/ReExpCalc.h
  ../pbam/SHCalc.cpp
  ../pbam/SHCalc.h
  ../pbam/System.cpp
  ../pbam/System.h
  ../pbam/readutil.h
  ../pbam/setup.cpp
  ../pbam/setup.h
  ../pbam/util.h
  ASolverTest.h
  ASolverUnitTest.h
  BDUnitTest.h
  BesselCalcTest.h
  BesselCalcUnitTest.h
  ConstantsTest.h
  ConstantsUnitTest.h
  ElectrostaticsUnitTest.h
  EnergyForceTest.h
  EnergyForceUnitTest.h
  MyMatrixUnitTest.h
  ReExpCalcTest.h
  ReExpCalcUnitTest.h
  SHCalcTest.h
  SHCalcUnitTest.h
  SystemUnitTest.h
  TerminateUnitTest.h
  ThreeBodyUnitTest.h
  main.cpp
  readutilUnitTest.h
  setupUnitTest.h
  tester.h
  utilUnitTest.h
  )
  
TARGET_LINK_LIBRARIES(pbamtest
    libgtest
)

add_test(NAME pbamtest
         COMMAND pbamtest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
  
install(TARGETS pbamtest RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
