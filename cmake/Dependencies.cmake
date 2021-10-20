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
    add_definitions(-D__MACOS)
    message(STATUS "BLAS_LINKER_FLAGS: ${BLAS_LINKER_FLAGS}")
    find_file(cblas_h cblas.h)
    message(STATUS "cblas.h: ${cblas_h}")
    get_target_property(include_dirs_blas BLAS::BLAS INCLUDE_DIRECTORIES)
    get_target_property(ifc_include_dirs_blas BLAS::BLAS INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(imported_loc BLAS::BLAS IMPORTED_LOCATION)
    get_target_property(imported_libname BLAS::BLAS IMPORTED_LIBNAME)
    message(STATUS "include dirs: ${include_dirs_blas}")
    message(STATUS "ifc include dirs: ${ifc_include_dirs_blas}")
    message(STATUS "imported_loc: ${imported_loc}")
    message(STATUS "imported_libname: ${imported_libname}")
  elseif(BLA_VENDOR STREQUAL "NAS")
    add_definitions(-D__XCODE)
  elseif(BLA_VENDOR STREQUAL "OpenBLAS")

  endif()

  add_definitions(-D__LAU)
  #include_directories(${???})
  list(APPEND PBSAM_LINKER_LIBS ${BLAS_LIBRARIES})

  MESSAGE( STATUS "linkers: " ${PBSAM_LINKER_LIBS} )
endif(ENABLE_PBSAM)
