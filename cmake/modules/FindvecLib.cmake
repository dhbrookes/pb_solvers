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

include(FindPackageHandleStandardArgs)

set(__vlib_incl_suff "Frameworks/vecLib.framework/Versions/Current/")
set(__accel_incl_suff "/System/Library/Frameworks/Accelerate.framework/")
find_path(vecLib_INCLUDE_DIR clapack.h
          PATHS /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
                /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX.sdk/System/Library/Frameworks/Accelerate.framework/Versions/Current/Frameworks/vecLib.framework/Headers/
                #/System/Library/${__vlib_incl_suff}
                #${__accel_incl_suff}${__vlib_incl_suff}
                #${CMAKE_OSX_SYSROOT}/${__accel_incl_suff}${__vlib_incl_suff}/Headers/
                NO_DEFAULT_PATH
         )

if (${vecLib_INCLUDE_DIR} MATCHES "vecLib_INCLUDE_DIR-NOTFOUND")
  find_path(vecLib_INCLUDE_DIR libBLAS.dylib
            PATHS /System/Library/${__accel_incl_suff}
           )
endif()

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
