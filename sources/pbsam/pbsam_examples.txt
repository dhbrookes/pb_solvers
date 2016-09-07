
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
  runname gly_0.05M.out
    
  units kT
  salt 0.05
  temp 298.15
  idiel 4
  sdiel 80
  
  attypes 1
  type 1 2
  tolsp 2.0
  surf 1 gly.vert #comment when system is CGed
  #imat 1 glyph   #uncomment when system has been run once
  #exp 1 gly      #uncomment when system has been run once
  pqr 1 gly.pqr
  xyz 1 zero.xyz

  
The files for PQR, Surf, and XYZ are:

name:  ``gly.pqr`` ( the first few lines)

.. code-block:: bash

  ATOM       1    C   A  0        0.0000   0.0000   0.0000  -0.1550 1.8700   
  ATOM       2    C   A  0        1.5000   0.0000   0.0000   0.6500 1.8700   
  ATOM       3    O   A  0        2.1870   0.9730   0.0000  -0.5330 1.7600   
  ATOM       4    O   A  0        1.9820  -1.2390   0.0110  -0.4280 1.5200   
  ATOM       5    C   A  0        3.3800  -1.4580   0.0810   0.1240 1.8700   
  ATOM       6    C   A  0        3.8360  -1.6880   1.5170   0.1220 1.8700   
  ATOM       7    C   A  0        5.1860  -2.3750   1.5470   0.1430 1.8700   
  ATOM       8    O   A  0        5.4810  -2.6680   2.9010  -0.4260 1.5200   
  ATOM       9    C   A  0        6.6410  -3.2580   3.1620   0.6450 1.8700   
  ATOM      10    C   A  0        6.8260  -3.4910   4.6340  -0.1550 1.8700 

name:  ``gly.vert`` ( the first few lines)

.. code-block:: bash

  # MSMS solvent excluded surface vertices for gly.xyzr
  #vertex #sphere density probe_r
      642      29  3.00  1.50
     -0.041    -1.794    -0.525    -0.022    -0.960    -0.280       0       1  2 
     -0.248    -1.747    -0.847     0.116    -0.991    -0.065       0      16  2 
      0.947    -2.243    -0.470    -0.681    -0.661    -0.317       0       4  2 
      0.101    -0.896    -1.738     0.433    -0.217    -0.875       0      16  2 
      0.416    -0.678    -1.692     0.223    -0.362    -0.905       0       1  2 
      1.084    -0.678    -1.692    -0.223    -0.362    -0.905       0       2  2 
      1.084     0.906    -1.582    -0.223     0.485    -0.846       0       2  2 


name:  ``zero.xyz`` 

.. code-block:: bash

    0.0   0.0   0.0
   12.0  12.0  12.0


To run:

.. code-block:: bash

  $$ ../../bin/pbsam run.energyforce.inp


Output Files
^^^^^^^^^^^^

And the resulting file:

name: ``gly_0.05M.out``

.. code-block:: bash

  My units are kT. Time: 0
  Molecule #1
      POSITION: [0, 0, 0]
      ENERGY: 6.17661e-05
      FORCE: 0.00072349, [-0.000537635 -0.000423847 -0.000233967]
      TORQUE: 2.03503e-06, [-4.31343e-05 -0.000822915 0.00078854]
  Molecule #2
      POSITION: [12, 12, 12] 
      ENERGY: 6.21059e-05
      FORCE: 0.000737173, [0.000535151 0.000445966 0.000241146]
    TORQUE: 8.2822e-06, [0.00196746 0.00132961 -0.00398844]

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
  surf 1 barnase.vert
  #pqr 1 barnase_cg.pqr #For when CG process is done
  #imat 1 barnasesph    #For when imat calculation is done
  #exp 1 barnase        #For when self-polarization is done
  xyz 1 zero.xyz


The files for PQR, Surf, and XYZ files are:

name:  ``barnase.pqr`` ( the first few lines)

.. code-block:: bash

  ATOM   1700  N    ALA B   1      20.757 52.394 30.692     0.1414  1.8240
  ATOM   1702  CA   ALA B   1      20.602 52.680 29.268     0.0962  1.9080
  ATOM   1703  C    ALA B   1      19.286 52.138 28.675     0.6163  1.9080
  ATOM   1704  O    ALA B   1      18.578 51.351 29.318    -0.5722  1.6612
  ATOM   1705  CB   ALA B   1      21.739 52.033 28.476    -0.0597  1.9080


name: ``barnase.vert`` ( the first few lines)

.. code-block:: bash

  # MSMS solvent excluded surface vertices for barnase.xyzr
  #vertex #sphere density probe_r
     5720     878  1.00  1.50
      5.097    50.485    18.262     0.322    -0.456     0.830       0     123  2 
      5.549    50.063    18.030     0.021    -0.174     0.984       0     121  2 
      6.503    50.902    19.073    -0.615    -0.734     0.289       0     133  2 
      5.437    49.956    18.007    -0.035    -0.228     0.973       0     121  2 
      4.986    50.378    18.239     0.266    -0.509     0.818       0     123  2 
      5.273    49.355    17.993     0.074     0.173     0.982       0     122  2 
      3.731    49.067    15.692    -0.890    -0.007    -0.457       0     122  2 

name:  ``zero.xyz``

.. code-block:: bash

    0.0   0.0   0.0

To run:

.. code-block:: bash

  $$ ../../bin/pbsam run.electrostatic.inp


Output Files
^^^^^^^^^^^^

And the resulting files:

name: ``barnase.dx``

.. code-block:: bash

  # Data from PBSAM Electrostat run
  # My runname is barnase.dx and units kT/e
  object 1 class gridpositions counts 100 100 100
  origin -25.025 -24.4258 -30.4642
  delta 0.538326 0.0e+00 0.0e+00
  delta 0.0e00 0.493468 0.0e+00
  delta 0.0e00 0.0e+00 0.563884
  object 2 class gridconnections counts 100 100 100
  object 3 class array type double rank 0 items 1000000 data follows
   0.003659521  0.003697636  0.003732662  0.003764229  0.003791946 
   0.003815395  0.003834137  0.003847709  0.003855628  0.003857388 
   0.003852465  0.003840319  0.003820396  0.003792134  0.003754962 
   0.003708309  0.003651608  0.003584305  0.003505857  0.003415750 
   0.003313498  0.003198656  0.003070826  0.002929665  0.002774897 

name: ``barnase_map.out``

.. code-block:: bash

  # Data from PBSAM Electrostat run
  # My runname is barnase_map.out and units kT
  grid 100 100 100
  origin -25.025 -24.4258 -30.4642
  delta 0.538326 0.493468 0.563884
     4.8667332   -1.1809119   -8.3553659    0.2419499 
     5.1270905   -4.6114465   -7.2974789    0.2407265 
     5.5570112   -3.1729867   -7.2710317    0.2437944 
     5.7783599   -1.6880478   -7.2445845    0.2337809 
     5.9996470   -4.5953014   -6.0809087    0.2450571 

name: ``barnase.y.0.dat``

.. code-block:: bash

  # Data from PBSAM Electrostat run
  # My runname is barnase.y.0.dat
  units kT/e
  grid 100 100
  axis y -0.245894
  origin -25.025 -30.4642
  delta 0.538326 0.563884
  maxmin 0.350066 -0.311879
     0.0074041     0.0076604     0.0079257     0.0081997     0.0084821 


.. _pbsam_dyn:

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
  term 1 contact contact.barn_bars

  attypes 2
  type 1 2 move 0.015 0.000045
  pqr 1 barnase.pqr
  #pqr 1 barnase_cg.pqr # for after CG
  surf 1 barnase.vert
  xyz 1 1 pos_1_1.xyz
  xyz 1 2 pos_1_2.xyz
  
  type 2 2 move 0.015 0.000045
  pqr 2 barstar.pqr
  #pqr 2 barstar_cg.pqr  # for after CG
  surf 2 barstar.vert
  xyz 2 1 pos_2_1.xyz
  xyz 2 2 pos_2_2.xyz
    

The files for PQR (first 5 lines) and XYZ files for the first trajectories are:

name:  ``barnase.pqr`` ( the first few lines)

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



name:  ``barstar.pqr`` ( the first few lines)

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


name: ``contact.barn_bars``

.. code-block:: bash

  1   872  2  1208    3.0
  1  1565  2   538    3.2 
  1   894  2   541    2.8 
  1   862  2   566    2.9 
  1   425  2   671    3.0 
  1  1242  2   474    2.7 
  1  1249  2   631    2.5 
  1  1248  2   683    3.1


To run:

.. code-block:: bash

  $$ ../../bin/pbsam run.dynamics.inp


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




