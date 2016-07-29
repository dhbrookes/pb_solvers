# Set gtest properties
ExternalProject_Get_Property(gtest source_dir binary_dir)
include_directories("${source_dir}/include")

file(GLOB SRCS *.cpp)

include_directories("../src")
include_directories("../../pb_shared/src")
include_directories("../../pb_shared/unittest")

ADD_EXECUTABLE(pbsamtest 
  ../../pb_shared/src/BesselCalc.cpp
  ../../pb_shared/src/Constants.cpp
  ../src/Electrostatics.cpp
  ../src/Gradsolvmat.cpp
  ../../pb_shared/src/ReExpCalc.cpp
  ../../pb_shared/src/SHCalc.cpp
  ../../pb_shared/src/setup.cpp
  ../src/PBSAM.cpp
  ../src/Solver.cpp
  ../src/Solvmat.cpp
  ../src/System.cpp
  ../src/TMatrix.cpp
  main.cpp
  )

TARGET_LINK_LIBRARIES(pbsamtest
    libgtest
)

add_test(NAME pbamtest
         COMMAND pbamtest)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbsamtest RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
