project(pbam)

include_directories("../../pb_shared/src")
include_directories("../../pb_wrap/src")
include_directories("./")

# C++ code for standalone
add_executable(pbam
  ASolver.cpp
  ../../pb_shared/src/BaseSys.cpp
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


################################################
###### APBS components
################################################
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

  #  build pbam library
  add_library(PBAMLib 
              ASolver.cpp
              ../../pb_shared/src/BaseSys.cpp
              BD.cpp
              ../../pb_shared/src/BesselCalc.cpp
              ../../pb_shared/src/Constants.cpp
              Electrostatics.cpp
              EnergyForce.cpp
              PBAM.cpp
              ../../pb_wrap/src/PBAMWrap.cpp
              ../../pb_shared/src/ReExpCalc.cpp
              ../../pb_shared/src/setup.cpp
              ../../pb_shared/src/SHCalc.cpp
              System.cpp
              )
  
  #  build standalone C example
  add_executable(PBAMWrap
                 ../../pb_wrap/src/runPBAMWrap.c
                 )

  target_link_libraries(PBAMWrap PBAMLib)
endif()


################################################
###### Sphinx
################################################
if(ENABLE_PBAM_SPHINX)
  include_directories("${PYTHON_INCLUDE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pbam_py.c"
                     COMMAND "${PYTHON_EXECUTABLE}" -m cython
                     -o "${CMAKE_CURRENT_BINARY_DIR}/pbam_py.c"
                     "${CMAKE_CURRENT_SOURCE_DIR}/pbam_sph.pyx"
                     MAIN_DEPENDENCY "${CMAKE_CURRENT_SOURCE_DIR}/pbam_sph.pyx")
  
  
  add_library(pbam_sph MODULE
              ../../pb_wrap/src/pbam_py.c
              ASolver.cpp
              BD.cpp
              ../../pb_shared/src/BesselCalc.cpp
              ../../pb_shared/src/Constants.cpp
              Electrostatics.cpp
              EnergyForce.cpp
              PBAM.cpp
              ../../pb_wrap/src/PBAMWrap.cpp
              ReExpCalc.cpp
              ../../pb_shared/src/setup.cpp
              ../../pb_shared/src/SHCalc.cpp
              System.cpp)
  
  set_target_properties(pbam_sph PROPERTIES PREFIX "")
endif()

################################################
###### Install
################################################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
