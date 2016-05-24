PB-AM code: Poisson-Boltzmann Semi-Analytical Method
============

Welcome to the home for the [PB-AM code](https://github.com/davas301/pb_solvers/tree/pbsam_dev)!

## Building


To build this branch, follow the instructions below:

~~~
mkdir build
cd build
cmake ..
make
~~~

This will make executables for pbam and pbamtest. Right now, the pbam 
Both executables can be found in the build directory subdirectories
`pbam/` and `pbam_test_code/`

## Running

To run, simply execute the following on the command line:

### For PBAM


~~~
./build/pbam/pbam
~~~


### For Google Tests

~~~
cd build # if not already in build
./build/pbsamtest/pbsamtest
~~~

Note that some tests are looking for paths relative to build for files,
so they will probably fail if you are not there!
