
add_executable(pbsam
	BesselCalc.cpp
	BesselCalc.h
	Constants.cpp
	Constants.h
	main.cpp
	MyMatrix.h
	readutil.h
	ReExpCalc.cpp
	ReExpCalc.h
	setup.cpp
	setup.h
	SHCalc.cpp
	SHCalc.h
  Solver.cpp
  Solver.h
	System.cpp
	System.h
	util.h
	)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/../bin)

install(TARGETS pbsam RUNTIME DESTINATION ${CMAKE_BINARY_DIR}/../bin)
