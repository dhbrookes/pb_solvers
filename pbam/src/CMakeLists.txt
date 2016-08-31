project(pbam)

include_directories("../../pb_shared/src")
include_directories("../../pb_wrap/src")
include_directories("./")

# C++ code for standalone
add_executable(pbam
  ASolver.cpp
  ../../pb_shared/src/BaseSys.cpp
  BDAM.cpp
  ../../pb_shared/src/BesselCalc.cpp
  ../../pb_shared/src/Constants.cpp
  ElectrostaticsAM.cpp
  PhysCalcAM.cpp
  main.cpp
  PBAM.cpp
  ../../pb_shared/src/ReExpCalc.cpp
  ../../pb_shared/src/setup.cpp
  ../../pb_shared/src/SHCalc.cpp
  SystemAM.cpp
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
              BDAM.cpp
              ../../pb_shared/src/BesselCalc.cpp
              ../../pb_shared/src/Constants.cpp
              ElectrostaticsAM.cpp
              PhysCalcAM.cpp
              PBAM.cpp
              ../../pb_wrap/src/PBAMWrap.cpp
              ../../pb_shared/src/ReExpCalc.cpp
              ../../pb_shared/src/setup.cpp
              ../../pb_shared/src/SHCalc.cpp
              SystemAM.cpp
              )
  
endif()

################################################
###### Sphinx
################################################
if(ENABLE_PBAM_SPHINX)
  include_directories("${PYTHON_INCLUDE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}")
  
  add_library(pbam_sph MODULE
              ../../pb_wrap/src/pbam_py.c
              ASolver.cpp
              BDAM.cpp
              ../../pb_shared/src/BesselCalc.cpp
              ../../pb_shared/src/Constants.cpp
              ElectrostaticsAM.cpp
              PhysCalcAM.cpp
              PBAM.cpp
              ../../pb_wrap/src/PBAMWrap.cpp
              ReExpCalc.cpp
              ../../pb_shared/src/setup.cpp
              ../../pb_shared/src/SHCalc.cpp
              SystemAM.cpp)
  
  set_target_properties(pbam_sph PROPERTIES PREFIX "")
endif()

################################################
###### Install
################################################

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
