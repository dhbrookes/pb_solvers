# Required for static linking and exported to PBSAM.cmake
set(PBSAM_LINKER_LIBS "")

# ---[ BLAS

# if BLA_VENDOR is set
#  find required
# else (BLA_VENDOR not set)
#   find ACML
#   if found
#     ...
#   endif
#   find MKL
#   if found
#     ...
#   endif
#   find Apple
#   if found
#     ...
#   endif
#   find NAS (veclib)
#   if found
#     ...
#   endif
#   find generic required
# endif

if(NOT BLA_VENDOR)
  set(BLA_VENDOR "OpenBLAS")
  message(STATUS "BLAS vendor was not set; looking for OpenBLAS")
endif()

find_package(BLAS REQUIRED)
message(STATUS "Found requested BLAS: <${BLA_VENDOR}>")

if(BLA_VENDOR MATCHES "^ACML[A-Za-z0-9_]*")
  add_definitions(-D__ACML)
elseif(BLA_VENDOR MATCHES "^Intel[A-Za-z0-9_]*")
  add_definitions(-D__MKL)
elseif(BLA_VENDOR STREQUAL "Apple")
  add_definitions(-D__MACOS)
elseif(BLA_VENDOR STREQUAL "NAS")
  add_definitions(-D__XCODE)
elseif(BLA_VENDOR STREQUAL "OpenBLAS")

endif()

add_definitions(-D__LAU)
#include_directories(${???})
list(APPEND PBSAM_LINKER_LIBS ${BLAS_LIBRARIES})


#if(NOT APPLE)
#  set(BLAS "Open" CACHE STRING "Selected BLAS library")
#  set_property(CACHE BLAS PROPERTY STRINGS "Atlas;Open;MKL")

#  if(BLAS STREQUAL "Open" OR BLAS STREQUAL "open")
#    find_package(OpenBLAS)
#    if( OpenBLAS_FOUND)
#      include_directories(${OpenBLAS_INCLUDE_DIR})
#      list(APPEND PBSAM_LINKER_LIBS ${OpenBLAS_LIB})
#      add_definitions(-D__MKL)
#      add_definitions(-D__LAU)
#    endif()
#  endif()
#elseif(APPLE)
#  find_package(vecLib)
#  if(vecLib_FOUND)
#    include_directories(${vecLib_INCLUDE_DIR})
#    MESSAGE( STATUS "Dependencies: " ${vecLib_INCLUDE_DIR} )
#    list(APPEND PBSAM_LINKER_LIBS ${vecLib_LINKER_LIBS})
#    add_definitions(-D__MACOS)
#    add_definitions(-D__LAU)
#  endif()
# MESSAGE( STATUS "sys root " ${CMAKE_OSX_SYSROOT} )
#endif()

MESSAGE( STATUS "linkers: " ${PBSAM_LINKER_LIBS} )
