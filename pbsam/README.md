PB-SAM code: Poisson-Boltzmann Semi-Analytical Method
============

Welcome to the home for the [PB-SAM code](https://github.com/dhbrookes/pb_solvers/tree/master/pbsam)!

## Building

To build this branch, follow the instructions below:

~~~
cd [pb_solvers_home]
mkdir build
cd build
cmake ..
make pbsam pbsamtest
~~~

This will make executables for pbsam and pbsamtest.  
Both executables can be found in the build directory subdirectory
`bin/`.

## Running

To run, simply execute the following on the command line:

### For PBSAM

~~~
./build/bin/pbsam run.inp
~~~


### For Google Tests

~~~
./build/bin/pbsamtest
~~~
