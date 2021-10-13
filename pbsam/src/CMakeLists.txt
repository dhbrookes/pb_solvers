set(pb_shared_dir ${pb_solvers_SOURCE_DIR}/pb_shared/src/)
set(pb_wrap_dir ${pb_solvers_SOURCE_DIR}/pb_wrap/src)

include_directories(
  ${pb_shared_dir}
  ${pb_wrap_dir}
  ${CMAKE_CURRENT_SOURCE_DIR}
)

if(ENABLE_PBSAM_APBS)
  add_definitions(-DPBSAM_APBS)
  include_directories(
    ${APBS_SOURCE_DIR}/src
    ${APBS_BINARY_DIR}/src
    ${fetk_SOURCE_DIR}/maloc/src/vsys
    ${fetk_BINARY_DIR}/maloc/src
  )
  list(APPEND PBSAMLib_src
    ${APBS_SOURCE_DIR}/src/generic/vatom.c
#    ${fetk_SOURCE_DIR}/maloc/src/vsys/vmem.c
  )
  list(APPEND PBSAM_LINKER_LIBS maloc)
  # NOTE: at the moment this will fail because we don't necessarily have vmem.c (if we're using a built version of FETK) 
  #if(APPLE)
  #  LIST(APPEND PBSAMLib_src
  #    ${apbs_SOURCE_DIR}/src/generic/vatom.c
  #    ../../../../externals/fetk/maloc/src/vsys/vmem.c
  #  )
  #endif()
  list(APPEND PBSAMLib_src
    ${pb_wrap_dir}/PBSAMWrap.cpp
  )
endif()

#  build pbsam library
list(APPEND PBSAMLib_src
  ${pb_shared_dir}/BaseBD.cpp
  ${pb_shared_dir}/BaseElectro.cpp
  ${pb_shared_dir}/BaseSys.cpp
  BDSAM.cpp
  ${pb_shared_dir}/BesselCalc.cpp
  ${pb_shared_dir}/Constants.cpp
  ElectrostaticsSAM.cpp
  Gradsolvmat.cpp
  PBSAM.cpp
  PhysCalcSAM.cpp
  ${pb_shared_dir}/ReExpCalc.cpp
  ${pb_shared_dir}/setup.cpp
  ${pb_shared_dir}/SHCalc.cpp
  Solver.cpp
  Solvmat.cpp
  SystemSAM.cpp
  TMatrix.cpp
)
if(WIN32)
  list(APPEND PBSAMLib_src ${pb_shared_dir}/drand48.cpp)
endif()

add_library(PBSAMLib ${PBSAMLib_src})
target_link_libraries(PBSAMLib ${PBSAM_LINKER_LIBS})

# build pbsam executable
add_executable(pbsam main.cpp)
target_link_libraries(pbsam PBSAMLib)


################################################
#set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbsam RUNTIME)
install(TARGETS PBSAMLib LIBRARY)
