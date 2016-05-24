add_executable(pbam
  ASolver.cpp
  ASolver.h
  BD.cpp
  BD.h
  BesselCalc.cpp
  BesselCalc.h
  Constants.cpp
  Constants.h
  Electrostatics.cpp
  Electrostatics.h
  EnergyForce.cpp
  EnergyForce.h
  main.cpp
  MyMatrix.h
  readutil.h
  ReExpCalc.cpp
  ReExpCalc.h
  setup.cpp
  setup.h
  SHCalc.cpp
  SHCalc.h
  System.cpp
  System.h
  util.h
  )

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

install(TARGETS pbam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
