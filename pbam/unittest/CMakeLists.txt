include_directories("${pb_solvers_SOURCE_DIR}/pbam/src")
include_directories("${pb_solvers_SOURCE_DIR}/pb_shared/src")
include_directories("${pb_solvers_SOURCE_DIR}/pb_shared/unittest")

add_definitions("-DPBAM_TEST_FILES_PATH=${pb_solvers_SOURCE_DIR}/pbam/pbam_test_files/gtest/")

add_executable(pbamtest 
  main.cpp
)
  
target_link_libraries(pbamtest 
  PBAMLib 
  gtest_main
)

gtest_discover_tests(pbamtest)

install(TARGETS pbamtest)
