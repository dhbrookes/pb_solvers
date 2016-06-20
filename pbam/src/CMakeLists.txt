project(pbam)

option(ENABLE_PBAM_APBS "pbam APBS" OFF)
if(ENABLE_PBAM_APBS)
   add_definitions(-DPBAM_APBS)
   include_directories(
         ../../../../src/ 
         ../../../../tools/windows/configuration\ headers/
         ../../../../externals/fetk/maloc/src/base
         ../../../../externals/fetk/maloc/src/vsys
         ../../../../externals/fetk/maloc/src/vsh
         ../../../../externals/fetk/maloc/src/psh
         )
endif()

include_directories("../../pb_shared/src")
include_directories("./")

# C++ code for standalone
add_executable(pbam
  ASolver.cpp
  BD.cpp
  ../../pb_shared/src/BesselCalc.cpp
  ../../pb_shared/src/Constants.cpp
  Electrostatics.cpp
  EnergyForce.cpp
  main.cpp
  PBAM.cpp
  ../../pb_shared/src/ReExpCalc.cpp
  ../../pb_shared/src/setup.cpp
  ../../pb_shared/src/SHCalc.cpp
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
  ../../pb_shared/src/BesselCalc.cpp
  ../../pb_shared/src/Constants.cpp
  Electrostatics.cpp
  EnergyForce.cpp
  PBAM.cpp
  PBAMWrap.cpp
  ../../pb_shared/src/ReExpCalc.cpp
  ../../pb_shared/src/setup.cpp
  ../../pb_shared/src/SHCalc.cpp
  System.cpp
  )

#
#  build standalone C example
#
add_executable(PBAMWrap
  runPBAMWrap.c
)


target_link_libraries(PBAMWrap PBAMLib)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
