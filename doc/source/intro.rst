

Getting Started
===============

Welcome to the documentation for the PB-[S]AM software suite! The theory underlying this code can is presented in [LoHe06]_, [YaHe10]_, [YaHe13]_

Installation
------------

To install the software suite, first clone the github repository::

$ git clone https://github.com/davas301/pb_solvers.git

Next, navigate to the :code:`pb_solvers` directory and follow the instructions below::

$ mkdir build
$ cd build
$ cmake ..
$ make pbam pbsam

This will make executables for pbam and pbsam. Both executables can then be found in :code:`./build/bin/`

Running
-------

To run the programs, simply execute the following on the command line::

For PB-AM
^^^^^^^^^
:code:`$ ./build/bin/pbam <input_file>`

For PB-SAM
^^^^^^^^^^
:code:`$ ./build/bin/pbsam <input_file>`

where :code:`<input_file>` is a PB-[S]AM input file described in detail in the proceeding sections of this documentation



