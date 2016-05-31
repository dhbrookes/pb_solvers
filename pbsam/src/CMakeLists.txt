
add_executable(pbsam
	BesselCalc.cpp
	Constants.cpp
	main.cpp
	ReExpCalc.cpp
	setup.cpp
	SHCalc.cpp
  Solver.cpp
	System.cpp
	)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

install(TARGETS pbsam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/bin)
