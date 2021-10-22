PB-AM code: Poisson-Boltzmann Analytical Method
============

Welcome to the home for the [PB-AM code](https://github.com/dhbrookes/pb_solvers)!

## Building

To build this method with testing enabled, follow the instructions below:

~~~
cd [pb_solvers_home]
mkdir build
cd build
cmake -DENABLE_PBAM=ON -DENABLE_PBSAM=OFF -DENABLE_PB_TESTING=ON ..
make
~~~

This will make executables for pbam and pbamtest. 
Both executables can be found in the build subdirectory `bin/`

If the user is building with macOS, then `vecLib` will be used to provide BLAS libraries.
If the user is building on another system, the default `cmake` command assumes that OpenBLAS is installed on the user's system, and uses OpenBLAS to provide the BLAS libraries.  

If you would like to instead  specify  either the MKL (Intel Math Kernel Library) or Atlas BLAS libraries when building on a non-macOS system, you can specify the `-DBLAS` flag when running the `cmake` command:

`cmake -DBLAS="Atlas" ..` for Atlas BLAS libraries

`cmake -DBLAS="MKL" ..` for MKL BLAS libraries

For information on installing the OpenBLAS libraries, see http://www.openblas.net.

## Running

To run, simply execute the following on the command line:

### For PBAM

~~~
./build/bin/pbam runfile.inp
~~~


### For Google Tests

~~~
cd [pb_solvers_home]/build
ctest
~~~
