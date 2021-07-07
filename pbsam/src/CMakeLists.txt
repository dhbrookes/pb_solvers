
include_directories("../../pb_shared/src")
include_directories("../../pb_wrap/src")

include_directories("./")

add_executable(pbsam
               ../../pb_shared/src/BaseBD.cpp
               ../../pb_shared/src/BaseElectro.cpp
               ../../pb_shared/src/BaseSys.cpp
               BDSAM.cpp
               ../../pb_shared/src/BesselCalc.cpp
               ../../pb_shared/src/Constants.cpp
               ElectrostaticsSAM.cpp
               Gradsolvmat.cpp
               main.cpp
               PBSAM.cpp
               PhysCalcSAM.cpp
               ../../pb_shared/src/ReExpCalc.cpp
               ../../pb_shared/src/setup.cpp
               ../../pb_shared/src/SHCalc.cpp
               Solver.cpp
               Solvmat.cpp
               SystemSAM.cpp
               TMatrix.cpp
               )

if(WIN32)
  target_sources(pbsam PUBLIC ../../pb_shared/src/drand48.cpp)
endif()

################################################
################################################
##### For APBS build of PBSAM
################################################
option(ENABLE_PBSAM_APBS "pbsam APBS" OFF)
if(ENABLE_PBSAM_APBS)
  add_definitions(-DPBSAM_APBS)
  include_directories(
    ../../../../src/
    ../../../../build/src
    ../../../../externals/fetk/maloc/src/base
    ../../../../externals/fetk/maloc/src/vsys
    ../../../../externals/fetk/maloc/src/vsh
    ../../../../externals/fetk/maloc/src/psh
    ##maybe?
    ../../../../src/generic
    )

  LIST(APPEND PBSAMLib_src
    ../../pb_shared/src/BaseBD.cpp
    ../../pb_shared/src/BaseElectro.cpp
    ../../pb_shared/src/BaseSys.cpp
    BDSAM.cpp
    ../../pb_shared/src/BesselCalc.cpp
    ../../pb_shared/src/Constants.cpp
    ElectrostaticsSAM.cpp
    Gradsolvmat.cpp
    PBSAM.cpp
    PhysCalcSAM.cpp
    ../../pb_wrap/src/PBSAMWrap.cpp
    ../../pb_shared/src/ReExpCalc.cpp
    ../../pb_shared/src/setup.cpp
    ../../pb_shared/src/SHCalc.cpp
    Solver.cpp
    Solvmat.cpp
    SystemSAM.cpp
    TMatrix.cpp
    )

  if(APPLE AND BUILD_SHARED_LIBS)
    LIST(APPEND PBSAMLib_src
      ../../../../src/generic/vatom.c
      ../../../../externals/fetk/maloc/src/vsys/vmem.c
      )
  endif()
  
  #  build pbsam library
  add_library(PBSAMLib ${PBSAMLib_src})

    TARGET_LINK_LIBRARIES(PBSAMLib ${PBSAM_LINKER_LIBS})

  if(WIN32)
    target_sources(PBSAMLib PUBLIC ../../pb_shared/src/drand48.cpp)
  endif()

  MESSAGE( STATUS "linkers: " ${APBS_LIBS})
  list(APPEND APBS_LIBS ${PBSAM_LINKER_LIBS})
  MESSAGE( STATUS "linkers: " ${APBS_LIBS})

endif()

################################################
TARGET_LINK_LIBRARIES( pbsam ${PBSAM_LINKER_LIBS})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
install(TARGETS pbsam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
