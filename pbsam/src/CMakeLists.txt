
include_directories("../../pb_shared/src")
include_directories("./")

add_executable(pbsam
  ../../pb_shared/src/BesselCalc.cpp
  ../../pb_shared/src/Constants.cpp
  Electrostatics.cpp
  Gradsolvmat.cpp  
  main.cpp
  PBSAM.cpp 
  ../../pb_shared/src/ReExpCalc.cpp
  ../../pb_shared/src/setup.cpp
  ../../pb_shared/src/SHCalc.cpp
  Solver.cpp
  Solvmat.cpp
  System.cpp
  TMatrix.cpp
  )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbsam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
