set(pb_shared_dir ${pb_solvers_SOURCE_DIR}/pb_shared/src/)
set(pb_wrap_dir ${pb_solvers_SOURCE_DIR}/pb_wrap/src)

include_directories(
  ${pb_shared_dir}
  ${pb_wrap_dir}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if(ENABLE_PBAM_APBS)
  add_definitions(-DPBAM_APBS)
  include_directories(
    ${apbs_SOURCE_DIR}/src
    ${Capbs_BINARY_DIR}/src
    ${FETK_INSTALL_PREFIX}/include
  )
endif()

#  build pbam library
add_library(PBAMLib 
  ASolver.cpp
  ${pb_shared_dir}/BaseBD.cpp
  ${pb_shared_dir}/BaseElectro.cpp
  ${pb_shared_dir}/BaseSys.cpp
  BDAM.cpp
  ${pb_shared_dir}/BesselCalc.cpp
  ${pb_shared_dir}/Constants.cpp
  ElectrostaticsAM.cpp
  PhysCalcAM.cpp
  PBAM.cpp
  ${pb_shared_dir}/ReExpCalc.cpp
  ${pb_shared_dir}/setup.cpp
  ${pb_shared_dir}/SHCalc.cpp
  SystemAM.cpp
)
if(WIN32)
  target_sources(PBAMLib PUBLIC ${pb_shared_dir}/drand48.cpp)
endif()

# C++ code for standalone
add_executable(pbam main.cpp)
target_link_libraries(pbam PBAMLib)

################################################
###### Sphinx
################################################
if(ENABLE_PBAM_SPHINX)
  include_directories(${PYTHON_INCLUDE_DIR})
  
  add_library(pbam_sph MODULE
    ${pb_wrap_dir}/pbam_py.c
    ${pb_wrap_dir}/PBAMWrap.cpp
  )
  target_link_libraries(pbam_sph PBAMLib)
  
  set_target_properties(pbam_sph PROPERTIES PREFIX "")
endif()

################################################
###### Install
################################################

#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbam RUNTIME)
install(TARGETS PBAMLib LIBRARY)
if(ENABLE_PBAM_SPHINX)
  #install(TARGETS pbam_sph LIBRARY)
endif()
