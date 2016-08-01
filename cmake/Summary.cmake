################################################################################################
# pbsam status report function.
# Automatically align right column and selects text based on condition.
# Usage:
#   pbsam_status(<text>)
#   pbsam_status(<heading> <value1> [<value2> ...])
#   pbsam_status(<heading> <condition> THEN <text for TRUE> ELSE <text for FALSE> )
function(pbsam_status text)
  set(status_cond)
  set(status_then)
  set(status_else)

  set(status_current_name "cond")
  foreach(arg ${ARGN})
    if(arg STREQUAL "THEN")
      set(status_current_name "then")
    elseif(arg STREQUAL "ELSE")
      set(status_current_name "else")
    else()
      list(APPEND status_${status_current_name} ${arg})
    endif()
  endforeach()

  if(DEFINED status_cond)
    set(status_placeholder_length 23)
    string(RANDOM LENGTH ${status_placeholder_length} ALPHABET " " status_placeholder)
    string(LENGTH "${text}" status_text_length)
    if(status_text_length LESS status_placeholder_length)
      string(SUBSTRING "${text}${status_placeholder}" 0 ${status_placeholder_length} status_text)
    elseif(DEFINED status_then OR DEFINED status_else)
      message(STATUS "${text}")
      set(status_text "${status_placeholder}")
    else()
      set(status_text "${text}")
    endif()

    if(DEFINED status_then OR DEFINED status_else)
      if(${status_cond})
        string(REPLACE ";" " " status_then "${status_then}")
        string(REGEX REPLACE "^[ \t]+" "" status_then "${status_then}")
        message(STATUS "${status_text} ${status_then}")
      else()
        string(REPLACE ";" " " status_else "${status_else}")
        string(REGEX REPLACE "^[ \t]+" "" status_else "${status_else}")
        message(STATUS "${status_text} ${status_else}")
      endif()
    else()
      string(REPLACE ";" " " status_cond "${status_cond}")
      string(REGEX REPLACE "^[ \t]+" "" status_cond "${status_cond}")
      message(STATUS "${status_text} ${status_cond}")
    endif()
  else()
    message(STATUS "${text}")
  endif()
endfunction()


################################################################################################
# Function for fetching pbsam version from git and headers
# Usage:
#   pbsam_extract_pbsam_version()
function(pbsam_extract_pbsam_version)
  set(pbsam_GIT_VERSION "unknown")
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
                    ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE
                    WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
                    OUTPUT_VARIABLE pbsam_GIT_VERSION
                    RESULT_VARIABLE __git_result)
    if(NOT ${__git_result} EQUAL 0)
      set(pbsam_GIT_VERSION "unknown")
    endif()
  endif()

  set(pbsam_GIT_VERSION ${pbsam_GIT_VERSION} PARENT_SCOPE)
  set(pbsam_VERSION "<TODO> (pbsam doesn't declare its version in headers)" PARENT_SCOPE)

  # pbsam_parse_header(${pbsam_INCLUDE_DIR}/pbsam/version.hpp pbsam_VERSION_LINES pbsam_MAJOR pbsam_MINOR pbsam_PATCH)
  # set(pbsam_VERSION "${pbsam_MAJOR}.${pbsam_MINOR}.${pbsam_PATCH}" PARENT_SCOPE)

  # or for #define pbsam_VERSION "x.x.x"
  # pbsam_parse_header_single_define(pbsam ${pbsam_INCLUDE_DIR}/pbsam/version.hpp pbsam_VERSION)
  # set(pbsam_VERSION ${pbsam_VERSION_STRING} PARENT_SCOPE)

endfunction()


################################################################################################
# Prints accumulated pbsam configuration summary
# Usage:
#   pbsam_print_configuration_summary()

function(pbsam_print_configuration_summary)
  pbsam_extract_pbsam_version()
  set(pbsam_VERSION ${pbsam_VERSION} PARENT_SCOPE)

  pbsam_merge_flag_lists(__flags_rel CMAKE_CXX_FLAGS_RELEASE CMAKE_CXX_FLAGS)
  pbsam_merge_flag_lists(__flags_deb CMAKE_CXX_FLAGS_DEBUG   CMAKE_CXX_FLAGS)

  pbsam_status("")
  pbsam_status("******************* pbsam Configuration Summary *******************")
  pbsam_status("General:")
  pbsam_status("  Version           :   ${pbsam_TARGET_VERSION}")
  pbsam_status("  Git               :   ${pbsam_GIT_VERSION}")
  pbsam_status("  System            :   ${CMAKE_SYSTEM_NAME}")
  pbsam_status("  C++ compiler      :   ${CMAKE_CXX_COMPILER}")
  pbsam_status("  Release CXX flags :   ${__flags_rel}")
  pbsam_status("  Debug CXX flags   :   ${__flags_deb}")
  pbsam_status("  Build type        :   ${CMAKE_BUILD_TYPE}")
  pbsam_status("")
  pbsam_status("  BUILD_SHARED_LIBS :   ${BUILD_SHARED_LIBS}")
  pbsam_status("  CPU_ONLY          :   ${CPU_ONLY}")
  pbsam_status("Dependencies:")
  pbsam_status("  BLAS              : " APPLE THEN "Yes (vecLib)" ELSE "Yes (${BLAS})")
# if(BUILD_docs)
#   pbsam_status("Documentaion:")
#   pbsam_status("  Doxygen           :" DOXYGEN_FOUND THEN "${DOXYGEN_EXECUTABLE} (${DOXYGEN_VERSION})" ELSE "No")
#   pbsam_status("  config_file       :   ${DOXYGEN_config_file}")

#   pbsam_status("")
# endif()
  pbsam_status("Install:")
  pbsam_status("  Install path      :   ${CMAKE_INSTALL_PREFIX}")
  pbsam_status("")
endfunction()
