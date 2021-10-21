include_directories("${pb_solvers_SOURCE_DIR}/pbsam/src")
include_directories("${pb_solvers_SOURCE_DIR}/pb_shared/src")
include_directories("${pb_solvers_SOURCE_DIR}/pb_shared/unittest")

add_definitions("-DPBSAM_TEST_FILES_PATH=${pb_solvers_SOURCE_DIR}/pbsam/pbsam_test_files/gtest/")
message(STATUS "%%%%%%% Adding definition: -DPBSAM_TEST_FILES_PATH=${pb_solvers_SOURCE_DIR}/pbsam/pbsam_test_files/gtest/")

add_executable(pbsamtest 
  main.cpp
)
  
target_link_libraries(pbsamtest 
  PBSAMLib 
  gtest_main
  ${PBSAM_LINKER_LIBS}
)

gtest_discover_tests(pbsamtest)

install(TARGETS pbsamtest)
