if(ENABLE_PBSAM)
  # Required for static linking and exported to PBSAM.cmake
  set(PBSAM_LINKER_LIBS "")

  # ---[ BLAS

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
    message(FATAL_ERROR "Cannot currently use Mac-native BLAS; recommend OpenBLAS")
    add_definitions(-D__MACOS)
    find_file(cblas_h cblas.h)
  elseif(BLA_VENDOR STREQUAL "NAS")
    add_definitions(-D__XCODE)
  elseif(BLA_VENDOR STREQUAL "OpenBLAS")

  endif()

  add_definitions(-D__LAU)
  #include_directories(${???})
  list(APPEND PBSAM_LINKER_LIBS ${BLAS_LIBRARIES})

  MESSAGE( STATUS "linkers: " ${PBSAM_LINKER_LIBS} )
endif(ENABLE_PBSAM)
