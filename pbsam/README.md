PB-SAM code: Poisson-Boltzmann Semi-Analytical Method
============

Welcome to the home for the [PB-SAM code](http://www.poissonboltzmann.org)!

## Building


To build this branch, follow the instructions below:

~~~
mkdir build
cd build
cmake ..
make
~~~

This will make executables for pbsam and pbsamtest. Right now, the pbsam 
executable gets put in a bin directory but for some reason the pbsamtest
does not? Both executables can be found in the build directory subdirectories
`pbsam/` and `pbsamtest/`

## Running

To run, simply execute the following on the command line:

### For PBSAM


~~~
./build/pbsam/pbsam
~~~


### For Google Tests

~~~
./build/pbsamtest/pbsamtest
~~~
