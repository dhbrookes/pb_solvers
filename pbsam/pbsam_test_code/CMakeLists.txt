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

include_directories("../pbsam")

ADD_EXECUTABLE(pbsamtest 
  ../pbsam/BesselCalc.cpp
  ../pbsam/BesselCalc.h
  ../pbsam/Constants.cpp
  ../pbsam/Constants.h
  ../pbsam/MyMatrix.h
  ../pbsam/readutil.h
  ../pbsam/ReExpCalc.cpp
  ../pbsam/ReExpCalc.h
  ../pbsam/setup.cpp
  ../pbsam/setup.h
  ../pbsam/SHCalc.cpp
  ../pbsam/SHCalc.h
  ../pbsam/Solver.cpp
  ../pbsam/Solver.h
  ../pbsam/System.cpp
  ../pbsam/System.h
  ../pbsam/util.h
  BesselCalcTest.h
  BesselCalcUnitTest.h
  CMakeLists.txt
  ConstantsTest.h
  ConstantsUnitTest.h
  MyMatrixUnitTest.h
  SHCalcTest.h
  SHCalcUnitTest.h
  SystemUnitTest.h
  main.cpp
  readutilUnitTest.h
  setupUnitTest.h
  tester.h
  utilUnitTest.h
  )

TARGET_LINK_LIBRARIES(pbsamtest
    libgtest
)

add_test(NAME pbamtest
         COMMAND pbamtest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

install(TARGETS pbsamtest RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
