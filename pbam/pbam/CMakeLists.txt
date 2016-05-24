project(pbam)

option(ENABLE_PBAM_APBS "pbam APBS" OFF)
if(ENABLE_PBAM_APBS)
   add_definitions(-DPBAM_APBS)
   include_directories(
         ../../../src/ 
         ../../../tools/windows/configuration\ headers/
         ../../../externals/fetk/maloc/src/base
         ../../../externals/fetk/maloc/src/vsys
         ../../../externals/fetk/maloc/src/vsh
         ../../../externals/fetk/maloc/src/psh
         )
endif()


# C++ code for standalone
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


####################
# APBS components?
####################

#
#  build pbam library
#
add_library(PBAMLib 
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

#
#  build standalone C example
#
add_executable(PBAMWrap
  runPBAMWrap.c
)


target_link_libraries(PBAMWrap PBAMLib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)
install(TARGETS pbam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
