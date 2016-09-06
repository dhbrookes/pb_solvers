
PB-SAM Examples
==================

This is a few examples of input files
and their expected outputs for each type of PB-SAM run. Like PB-AM
runs can be one of three options:

- :ref:`pbsam_enfo`
- :ref:`pbsam_elec`
- :ref:`pbsam_dyn`

Because PB-SAM is a bit more theoretically involved than PB-AM, it is 
time saving to save some files with program intermediates that 
can be later used as program inputs. They are described in the following section.


PB-SAM intermediates
---------------------

When running PB-SAM, a few intermediate quantities are generated
and saved to files. These include

- :ref:`pqrlab`
- :ref:`imatlab`
- :ref:`explab`

.. _pqrlab:

Coarse-Grain PQR files
^^^^^^^^^^^^^^^^^^^^^^^

The program uses the vertex file output from MSMS
to coarse-grain the system. If the input file does not
contain any CG centers (indicated by ``CEN`` keyword in 
the PQR file, then the MC CG process is run, and the output
is stored to a file called ``[pqr input]_cg.pqr``, where 
``[pqr input]`` is the input file name with the last four characters
removed (presumably '.pqr').

.. code-block:: bash

  pqr 1 barnase.pqr
  surf 1 barnase.vert

**Later use:** If you wish to run the system again, you can 
change the filename of the PQR file from the original to
the new ``[pqr]_cg.pqr``, and the CG process will be skipped.

.. code-block:: bash

  pqr 1 barnase_cg.pqr # _cg file has replaced pqr
  #surf 1 barnase.vert # commented out, no longer needed

.. _imatlab:

Imat: Surface integral files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the system has been coarse-grained, surface integral
matrices are generated for each CG sphere. These are binary
files, and will be named as follows: ``[pqr input]sph[#].bin``,
where like the CG output file, the pqr input has removed the last
four characters, and the ``#`` indicating the number of the sphere
within the molecule (zero-based).

**Later use:** If you wish to run the system again, you
can add the flag ``imat`` into the input file, with the prefix
``[pqr input]sph`` for the molecule. The program will append the sphere
numbers after the ``sph``. Input file invocation is given below.
Please note that if system conditions (dielectric
constants, temperature, salt concentration) are changed, the ``imat`` 
files do not need to be regenerated.

.. code-block:: bash

  pqr 1 barnase.pqr
  imat 1 barnasesph 

.. _explab:

Exp: Expansion files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After the system has been coarse-grained and surface integral
generated, the program will then perform a self-polarization for
each molecule type designated in the input. Once the self-polarization
has completed, the multipole expansions **H** and **F** representing
the effective charge distribution on the molecular surface are printed
to files, called: ``[pqr input].H.[#].exp`` and ``[pqr input].F.[#].exp``.
where like the CG output file, the pqr input has removed the last
four characters, and the ``#`` indicating the number of the sphere
within the molecule (zero-based).

**Later use:** If you wish to run the system again, you
can add the flag ``exp`` into the input file, with the prefix
``[pqr input]`` for the molecule. The program will append letter H or F
and the sphere numbers. Please note that if system conditions (dielectric
constants, temperature, salt concentration) are changed, the ``exp`` 
files should be regenerated.

.. code-block:: bash

  pqr 1 barnase.pqr
  exp 1 barnase

.. _pbsam_enfo:

Physical Calculation run
------------------------


Input Files
^^^^^^^^^^^

name:  ``run.energyforce.inp``

.. code-block:: bash

  runtype energyforce

  
The files for PQR and XYZ are:

name:  ``single_charge.pqr`` 

.. code-block:: bash

  ATOM      1  N   NTR     0       1.000   0.000   0.000 -1.0000 3.7300


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

.. _pbsam_elec:

Electrostatics run
------------------


Input Files
^^^^^^^^^^^

name:  ``run.electrostatic.inp``

.. code-block:: bash

  runtype electrostatics 150
  runname barnase.out
  
  units kT
  salt 0.01
  temp 298
  idiel 2 
  sdiel 78

  dx barnase.dx
  
  3dmap barnase_map.out
  
  gridct 2
  grid2D 1 barnase.y.0.dat y 0
  grid2D 2 barnase.z.0.dat z 0
  
  attypes 1
  type 1 1
  pqr 1 barnase.pqr  
  #pqr 1 barnase_cg.pqr #For when CG process is done
  #imat 1 barnasesph    #For when imat calculation is done
  #exp 1 barnase        #For when self-polarization is done
  xyz 1 zero.xyz


The files for PQR and XYZ files are:

name:  ``barnase.pqr``

.. code-block:: bash

  ATOM   1700  N    ALA B   1      20.757 52.394 30.692     0.1414  1.8240
  ATOM   1702  CA   ALA B   1      20.602 52.680 29.268     0.0962  1.9080
  ATOM   1703  C    ALA B   1      19.286 52.138 28.675     0.6163  1.9080
  ATOM   1704  O    ALA B   1      18.578 51.351 29.318    -0.5722  1.6612
  ATOM   1705  CB   ALA B   1      21.739 52.033 28.476    -0.0597  1.9080

name:  ``zero.xyz``

.. code-block:: bash

    0.0   0.0   0.0

name: ``barnase.vert``

.. code-block:: bash

  ATOM      1  O    CEN     1      15.773   43.159   13.061  0.00  12.5813
  ATOM      2  O    CEN     2      29.519   49.345   20.442  0.00  10.3065
  ATOM      3  O    CEN     3      13.647   40.795    6.493  0.00  11.5453
  ATOM      4  O    CEN     4      10.614   49.582   13.541  0.00  10.4199


To run:

.. code-block:: bash

  $$ ../../bin/pbsam run.electrostatic.inp

.. _pbsam_dyn:

Output Files
^^^^^^^^^^^^

And the resulting files:

name: ``barnase.dx``

.. code-block:: bash

  # Data from PBAM Electrostat run
  # My runname is out.dx and units kT/e

name: ``barnase_map.out``

.. code-block:: bash

  # Data from PBAM Electrostat run


name: ``barnase.y.0.dat``

.. code-block:: bash

  # Data from PBAM Electrostat run


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
  term 1 contact contact.dat
  
  attypes 2
  type 1 2 move 0.015 0.000045
  pqr 1 1BRS_chainA.pqr
  surf 1 1BRS_chainA.vert
  xyz 1 1 pos_1_1.xyz
  xyz 1 2 pos_1_2.xyz
  
  type 2 2 move 0.015 0.000045
  pqr 2 1BRS_chainD.pqr
  surf 2 1BRS_chainD.vert
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




