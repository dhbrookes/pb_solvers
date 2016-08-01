PB-AM code: Poisson-Boltzmann Analytical Method
============

Welcome to the home for the [PB-AM code](https://github.com/davas301/pb_solvers)!

## Building

To build this method, follow the instructions below:

~~~
cd [pb_solvers_home]
mkdir build
cd build
cmake ..
make pbam pbamtest
~~~

This will make executables for pbam and pbamtest. Right now, the pbam 
Both executables can be found in the build subdirectory `bin/`

## Running

To run, simply execute the following on the command line:

### For PBAM

~~~
./build/bin/pbam runfile.inp
~~~


### For Google Tests

~~~
cd build # if not already in build
./build/bin/pbsamtest
~~~

Note that some tests are looking for paths relative to build for files,
so they will probably fail if you are not there!
