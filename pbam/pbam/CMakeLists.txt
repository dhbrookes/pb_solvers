project(pbam)

option(ENABLE_PBAM_APBS "pbam APBS" OFF)
if(ENABLE_PBAM_APBS)
   add_definitions(-DPBAM_APBS)
   include_directories(
         ../../../src/ 
         ../../../src/generic
         ../../../tools/windows/configuration\ headers/
         )
endif()


# C++ code for standalone
add_executable(pbam
  ASolver.cpp
  BD.cpp
  BesselCalc.cpp
  Constants.cpp
  Electrostatics.cpp
  EnergyForce.cpp
  main.cpp
  PBAM.cpp
  ReExpCalc.cpp
  setup.cpp
  SHCalc.cpp
  System.cpp
  )


####################
# APBS components?
####################

#
#  build pbam library
#
add_library(PBAMLib 
  ASolver.cpp
  BD.cpp
  BesselCalc.cpp
  Constants.cpp
  Electrostatics.cpp
  EnergyForce.cpp
  PBAM.cpp
  PBAMWrap.cpp
  ReExpCalc.cpp
  setup.cpp
  SHCalc.cpp
  System.cpp
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
