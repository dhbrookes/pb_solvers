# Find the vecLib libraries as part of Accelerate.framework 
#  or as standalone framework
#
# The following are set after configuration is done:
#  VECLIB_FOUND
#  vecLib_INCLUDE_DIR
#  vecLib_LINKER_LIBS


if(NOT APPLE)
  return()
endif()

set(__vlib_incl_suff "Frameworks/vecLib.framework/Versions/Current/Headers")
set(__accel_incl_suff "Frameworks/Accelerate.framework/Versions/Current/")
find_path(vecLib_INCLUDE_DIR clapack.h
          PATHS /System/Library/${__veclib_include_suffix}
                /System/Library/${__accel_incl_suff}${__vlib_incl_suff}
         )


if (${vecLib_INCLUDE_DIR} MATCHES "vecLib_INCLUDE_DIR-NOTFOUND")
  set(__vlib_incl_suff "Frameworks/vecLib.framework/Versions/Current/")
  find_path(vecLib_INCLUDE_DIR libBLAS.dylib
            PATHS /System/Library/${__vlib_incl_suff}
                  /System/Library/${__accel_incl_suff}${__vlib_incl_suff}
           )

  MESSAGE( STATUS "incl 1: " ${vecLib_INCLUDE_DIR} )
endif()  

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(vecLib DEFAULT_MSG vecLib_INCLUDE_DIR)

if(VECLIB_FOUND)
  if(vecLib_INCLUDE_DIR MATCHES "^/System/Library/Frameworks/vecLib.framework.*")
    set(vecLib_LINKER_LIBS -lcblas "-framework vecLib")
    message(STATUS "Found standalone vecLib.framework")
  else()
    set(vecLib_LINKER_LIBS -lcblas "-framework Accelerate")
    message(STATUS "Found vecLib as part of Accelerate.framework")
  endif()

  mark_as_advanced(vecLib_INCLUDE_DIR)
endif()
