# Set gtest properties
ExternalProject_Get_Property(gtest source_dir binary_dir)
include_directories("${source_dir}/include")

file(GLOB SRCS *.cpp)

include_directories("../src")
include_directories("../../pb_shared/src")
include_directories("../../pb_shared/unittest")
include_directories("../../pb_wrap/src")

ADD_EXECUTABLE(pbsamtest
  ../../pb_shared/src/BaseElectro.cpp
  ../../pb_shared/src/BaseSys.cpp
  ../src/BDSAM.cpp
  ../../pb_shared/src/BesselCalc.cpp
  ../../pb_shared/src/Constants.cpp
  ../src/ElectrostaticsSAM.cpp
  ../src/Gradsolvmat.cpp
  ../src/PBSAM.cpp
  ../src/PhysCalcSAM.cpp
  ../../pb_shared/src/ReExpCalc.cpp
  ../../pb_shared/src/SHCalc.cpp
  ../../pb_shared/src/setup.cpp
  ../src/Solver.cpp
  ../src/Solvmat.cpp
  ../src/SystemSAM.cpp
  ../src/TMatrix.cpp
  main.cpp
  )

add_test(NAME pbamtest
         COMMAND pbamtest)

TARGET_LINK_LIBRARIES(pbsamtest
                      libgtest
                      ${PBSAM_LINKER_LIBS})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbsamtest RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
