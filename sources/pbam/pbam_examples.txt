
PB-AM Examples
==============

This is a few examples of input files
and their expected outputs for each type of PB-AM run.
The following runs are described in the previous section,
and examples of their inputs and outputs follow:

- :ref:`physlab`
- :ref:`electlab`
- :ref:`dynamlab`

.. _physlab:

Physical Calculation run
------------------------


Input Files
^^^^^^^^^^^

name:  ``run.energyforce.inp``

.. code-block:: bash

  runtype energyforce
  runname energyforce.2sp.jmol.out
  
  units jmol
  salt 0.01
  temp 353
  idiel 4 
  sdiel 78
  
  attypes 1
  type 1 2
  pqr 1 single_charge.pqr
  xyz 1 positions_2.xyz

  
The files for PQR and XYZ are:

name:  ``single_charge.pqr`` 

.. code-block:: bash

  ATOM      1  N   NTR     0       1.000   0.000   0.000 -1.0000 3.7300
  ATOM      1  N   NTR     0       0.000   1.000   0.000 -1.0000 6.3200


name:  ``positions_2.xyz`` 

.. code-block:: bash

  -10.0  23.4  -8.7
    0.0   0.0  -2.5


To run:

.. code-block:: bash

  $$ ../../bin/pbam run.energyforce.inp


Output Files
^^^^^^^^^^^^

And the resulting file:

name: ``energyforce.2sp.jmol.out`` 

.. code-block:: bash

  My units are Joules/Mol
  MOLECULE #1
          POSITION: [-10, 23.4, -8.7]
          ENERGY: 1328.86
          FORCE: 1.19858e+08, [-36.6012 1.19858e+08 -1.65408e-05]
          TORQUE: 1.78426e+06, [1.28425 1.78426e+06 1.9978e-06]
  MOLECULE #2
          POSITION: [0, 0, -2.5]
          ENERGY: 1328.86
          FORCE: 1.19858e+08, [36.6012 -1.19858e+08 1.65408e-05]
          TORQUE: 1.78354e+06, [1.28372 1.78354e+06 1.99699e-06]


.. _electlab:

Electrostatics run
------------------


Input Files
^^^^^^^^^^^

name:  ``run.electrostatic.inp``

.. code-block:: bash

  runtype electrostatics 140
  runname electrostatic
  
  units kT
  salt 0.01
  temp 298
  idiel 4 
  sdiel 78
  
  dx out.dx
  
  3dmap electro_map.out
  
  gridct 2
  grid2D 1 out.x.0.dat x 0
  grid2D 2 out.x.-1.dat x -1
  
  attypes 2
  type 1 2
  pqr 1 single_charge.pqr
  xyz 1 positions_2.xyz
  
  type 2 2
  pqr 2 pos_charge.pqr
  xyz 2 positions_pos.xyz



The files for PQR and XYZ files are:

name:  ``single_charge.pqr``

.. code-block:: bash

  ATOM      1  N   NTR     0       0.000   1.000   0.000  4.0000 0.3200
  ATOM      1  N   NTR     0       0.000   0.000  -1.000  4.0000 0.3200
  ATOM      1  X   CEN     0       0.000   0.000   0.000  0.0000 2.0000




name:  ``positions_2.xyz``

.. code-block:: bash

    0.0   0.0  -5.0
    0.0   0.0   5.0



name:  ``pos_charge.pqr``

.. code-block:: bash

  ATOM      1  N   NTR     0       0.000   1.000   0.000 -4.0000 0.3200
  ATOM      1  N   NTR     0       0.000   0.000  -1.000 -4.0000 0.3200
  ATOM      1  X   CEN     0       0.000   0.000   0.000  0.0000 2.0000




name:  ``positions_pos.xyz``

.. code-block:: bash

    0.0   5.0   0.0
    0.0  -5.0   0.0



To run:

.. code-block:: bash

  $$ ../../bin/pbam run.electrostatic.inp



Output Files
^^^^^^^^^^^^

And the resulting files:

name: ``out.dx``

.. code-block:: bash

  # Data from PBAM Electrostat run
  # My runname is out.dx and units kT/e
  object 1 class gridpositions counts 140 140 140
  origin -4 -9 -9
  delta 0.0571429 0.0e+00 0.0e+00
  delta 0.0e00 0.128571 0.0e+00
  delta 0.0e00 0.0e+00 0.128571
  object 2 class gridconnections counts 140 140 140
  object 3 class array type double rank 0 items 2744000 data follows
  2.7203115e-01  3.0271755e-01  3.3459723e-01  
  3.6769040e-01  4.0201595e-01  4.3759129e-01 
  .....
  -1.3185519e-01  -1.5849252e-01  -1.8359631e-01
  -2.0722087e-01  -2.2942006e-01  -2.5024714e-01
  -2.6975467e-01  -2.8799442e-01
  attribute "dep" string "positions"
  object "regular positions regular connections" class field
  component "positions" value 1
  component "connections" value 2
  component "data" value 3



name: ``electro_map.out``

.. code-block:: bash

  # Data from PBAM Electrostat run
  # My runname is electro_map.out and units kT/e
  grid 10 10 10
  origin -4 -9 -9
  delta 0.8 1.8 1.8
    0.00825   0.00006  -2.90002 -5.899956 
    0.00822   0.00071  -2.90002 -5.902602 



name: ``out.x.0.dat``

.. code-block:: bash

  # Data from PBAM Electrostat run
  # My runname is out.x.0.dat
  units kT
  grid 140 140 
  axis x 0 
  origin -9 -9
  delta 0.128571 0.128571
  maxmin 39.23 -39.23
     0.3605004     0.4030045     0.4474874     0.4940082     0.5426260     0.5933995


.. _dynamlab:

Dynamics run
------------


Input Files
^^^^^^^^^^^

name:  ``run.dynamics.inp``

.. code-block:: bash

  runtype dynamics 2
  runname dyn_cont_barn
  
  salt 0.01
  temp 298
  idiel 4 
  sdiel 78
    
  termct 1
  termcombine or
  term 1 contact 2.5 1 2
  
  attypes 2
  type 1 2 move 0.015 0.000045
  pqr 1 1BRS_chainA.pqr
  xyz 1 1 pos_1_1.xyz
  xyz 1 2 pos_1_2.xyz
  
  type 2 2 move 0.015 0.000045
  pqr 2 1BRS_chainD.pqr
  xyz 2 1 pos_2_1.xyz
  xyz 2 2 pos_2_2.xyz



The files for PQR (first 5 lines) and XYZ files for the first trajectories are:

name:  ``1BRS_chainA.pqr``

.. code-block:: bash

  ATOM   1700  N    ALA B   1      20.757 52.394 30.692     0.1414  1.8240
  ATOM   1702  CA   ALA B   1      20.602 52.680 29.268     0.0962  1.9080
  ATOM   1703  C    ALA B   1      19.286 52.138 28.675     0.6163  1.9080
  ATOM   1704  O    ALA B   1      18.578 51.351 29.318    -0.5722  1.6612
  ATOM   1705  CB   ALA B   1      21.739 52.033 28.476    -0.0597  1.9080
  



name:  ``pos_1_1.xyz``

.. code-block:: bash

  61.25 61.25 61.25
  -26.25 61.25 -26.25



name:  ``1BRS_chainD.pqr``

.. code-block:: bash

  ATOM      1  N    LYS D   1      48.330 40.393  9.798     0.0966  1.8240
  ATOM      2  CA   LYS D   1      47.401 39.287  9.370    -0.0015  1.9080
  ATOM      3  C    LYS D   1      47.507 38.911  7.890     0.7214  1.9080
  ATOM      4  O    LYS D   1      47.126 39.582  6.905    -0.6013  1.6612
  ATOM      5  CB   LYS D   1      45.995 39.632  9.817     0.0212  1.9080
  



name:  ``pos_2_1.xyz``

.. code-block:: bash

  -26.25 61.25 61.25
  61.25 -26.25 61.25



To run:

.. code-block:: bash

  $$ ../../bin/pbam run.dynamics.inp


Output Files
^^^^^^^^^^^^

And the resulting files:

name: ``dyn_cont_barn_[traj#].xyz`` VMD readable XYZ file 
that shows the trajectory of molecules in the system. The 
time that is snapshot was printed from is given on the 
same line as the word Atom. The atoms of your input file are 
currently labeled N, and the coarse-grain center is labeled "X" 
in the first column of the XYZ file.

.. code-block:: bash

  3135
  Atoms. Timestep (ps): 0
  N   -7.241   -0.530   18.703
  N   -6.015   -0.503   17.910
  N   -5.784    0.840   17.188
  N   -6.682    1.690   17.128
  N   -6.066   -1.580   16.827
  N   -7.519   -1.481   18.863
  N   -7.084   -0.079   19.584



name: ``dyn_cont_barn_[traj\#].dat`` Statistics from simulation 
printed out at the same time as each XYZ snapshot. The energy 
is not computed and should be ignored.

.. code-block:: bash

  My units are Internal. Time (ps) 500.4
  MOLECULE #1
      POSITION: [0, 0, 0]
      ENERGY: 0
      FORCE: 3.39124e-06, [1.69863e-06 2.07547e-06 6.5356e-07]
      TORQUE: 2.55224e-05, [-2.11728e-05 1.00774e-05 3.08631e-05]
  MOLECULE #2
      POSITION: [87.211, 43.861, 21.691]
      ENERGY: 0
      FORCE: 3.65373e-06, [-1.87502e-06 -2.21744e-06 -7.27314e-07]
      TORQUE: 1.91656e-05, [8.14396e-06 -1.22678e-05 1.56284e-05]



name: ``dyn_nam_barn.stat`` Details about how each simulation has 
terminated and the time at which this occurred.

.. code-block:: bash

  Molecule type 1 has fulfilled condition: r >= 500.00;    at time (ps) 1.32367e+06
  Molecule type 1 has fulfilled condition: r >= 500.00;    at time (ps) 1.15712e+06
  System has fulfilled condition: Type 0 and Type 1 are within  2.50;  at time (ps) 1.90603e+06
  Molecule type 1 has fulfilled condition: r >= 500.00;    at time (ps) 2.18533e+06
  System has fulfilled condition: Type 0 and Type 1 are within  2.50;  at time (ps) 1.59066e+06




