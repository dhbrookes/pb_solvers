
PB-SAM Program details
=======================

Some details about the code: it is written in C++11 and it
has some linear algebra utilities and OMP parallelism implemented.

Below is a flow diagram of the program from call to finish.

.. image:: ../images/pbsam_flow.png
   :width: 85%

Additionally, some of the pre-run components are quite time-consuming
to compute, so if they are not included as inputs, they will be printed to 
the working directory. They can then be used for future runs, forgoing
the need to re-calculate them.

