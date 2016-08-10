
PB-SAM Inputs
=============

Here are the inputs that are used in a PB-SAM run, and
their meanings. There is much overlap between PB-AM inputs and PB-SAM, but 
the full PB-SAM list is provided here. The general format of the input is:

``keyword value``

The keyword that is most important is the `runtype` keyword.
This specifies which type of output you are interested in generating.

The current options are: `energyforce`, `electrostatics` 
and `dynamics`. The keywords for each specific run follow 
the next section.

General input file parameters
-----------------------------

+-------------+--------------------+--------------------------------------------------------+
| *Keyword*   |  *Parameters*      |  *Description*                                         |
|             |                    |                                                        |
+=============+====================+========================================================+
| runname     | `name`             | `name` is desired internal name of this run.           |
+-------------+--------------------+--------------------------------------------------------+
| attypes     | `numtypes`         | Set the number of different atom                       |
|             |                    |                                                        |
|             |                    | types to `numtypes`                                    |
+-------------+--------------------+--------------------------------------------------------+
| pqr         | `idx`  `fpath`     | The molecule index for this xyz file.                  |
|             |                    |                                                        |
|             |                    | Provide input PQR file at `fpath`                      |
+-------------+--------------------+--------------------------------------------------------+
| xyz         | `idx`  `fpath`     | The molecule index for this xyz file.                  |
|             |                    |                                                        |
|             |                    | Provide input XYZ file at `fpath`                      |
+-------------+--------------------+--------------------------------------------------------+
| surf        | `idx`  `fpath`     | The molecule index for this surface/vertex file.       |
|             |                    |                                                        |
|             |                    | Provide input solvent vertex file at `fpath`. See      |
|             |                    |                                                        |
|             |                    | details of the vertex file below.                      |
+-------------+--------------------+--------------------------------------------------------+
| tolsp       | `tolerance`        | The desired tolerance for coarse-graining, indicates   |
|             |                    |                                                        |
|             |                    | how far beyond the surface the CG spheres may extend   |
|             |                    |                                                        |
|             |                    | Larger tolsp yields a coarser molecular representation |
+-------------+--------------------+--------------------------------------------------------+
|  transrot   | `idx`  `fpath`     | The molecule type index for this translation/rotation  | 
|             |                    |                                                        |
|             |                    | file that can be input instead of a xyz file.          |
|             |                    |                                                        |
|             |                    | Provide input file at `fpath`                          |
+-------------+--------------------+--------------------------------------------------------+
|  randorient |                    | If you want your molecules to be randomly              |
|             |                    |                                                        |
|             |                    | rotated, use this flag                                 |
+-------------+--------------------+--------------------------------------------------------+
|  units      | `units`            | Set the units of output to `units`.                    |
|             |                    |                                                        |
|             |                    | The current options are: `jmol` (Joules/mole),         |
|             |                    |                                                        |
|             |                    | `kT` (kT/e) and `kcalmol` (kCal/mole).                 |
+-------------+--------------------+--------------------------------------------------------+
|  salt       | `con`              | Set salt concentration in the system to `con`          |
+-------------+--------------------+--------------------------------------------------------+
|  temp       | `T`                | Set system temperature to `T`                          |
+-------------+--------------------+--------------------------------------------------------+
|  idiel      | `ival`             | Set the interior dielectric constant to `ival`         |
+-------------+--------------------+--------------------------------------------------------+
|  sdiel      | `sval`             | Set the solvent dielectric constant to `sval`          |
+-------------+--------------------+--------------------------------------------------------+
|  pbc        | `boxlength`        | Set size of periodic box to `boxlength`                |
+-------------+--------------------+--------------------------------------------------------+
|  random     | `seed`             | Seed the random number generator with `seed`           |
+-------------+--------------------+--------------------------------------------------------+
|  type       |   `idx` `ct`       | Set attributes of an atom type, where `idx` is the     | 
|             |                    |                                                        |
|             |   `movetype` `dtr` | integer id of this type, which can be 1 to `numtypes`  |
|             |                    |                                                        |
|             |   `drot`           | (above). `ct` is the number of atoms of this type      |
|             |                    |                                                        |
|             |                    | in the system and `movetype` describes the way this    |
|             |                    |                                                        |
|             |                    | type is allowed to move in a dynamics run (`move`,     |
|             |                    |                                                        |
|             |                    | `rot`, or `stat`). If `movetype` is `move` , then a    |
|             |                    |                                                        |
|             |                    | translational diffusion coefficient `dtr` and a        |
|             |                    |                                                        |
|             |                    | rotational diffusional coefficient `drot` are          |
|             |                    |                                                        |
|             |                    | required. If `movetype` is `rot` then just `drot`      |
|             |                    |                                                        |
|             |                    | is required.                                           |
+-------------+--------------------+--------------------------------------------------------+

Physical calculations
---------------------

Will calculate the interaction energy, forces and torques
for all molecules in the system. 

+-------------+--------------------+--------------------------------------------------------+
| *Keyword*   |  *Parameters*      |  *Description*                                         |
|             |                    |                                                        |
+=============+====================+========================================================+
| runtype     | `energyforce`      | This will invoke the calculation of energies, forces   |
|             |                    |                                                        |
|             | `outfilename`      | and torques for the input system.                      |
+-------------+--------------------+--------------------------------------------------------+

Electrostatic Potential Visualization
-------------------------------------

Will compute the electrostatic potential (ESP)
and print out 2D and 3D geometries for visualization.
Post-processing visualization can be performed in VMD
and with the python scripts provided in the repository.

+-------------+--------------------+--------------------------------------------------------+
| *Keyword*   |  *Parameters*      |  *Description*                                         |
|             |                    |                                                        |
+=============+====================+========================================================+
| runtype     | `electrostatics`   | Will write the results of electrostatics calculations  |
|             |                    |                                                        |
|             | `gridpts`          | for every 3D grid point to `fname`, in the same output |
|             |                    |                                                        |
|             |                    | format as an APBS dx file.                             |
+-------------+--------------------+--------------------------------------------------------+
| dx          | `energyforce`      | Will write the results of electrostatics calculations  |
|             |                    |                                                        |
|             | `outfilename`      | for every 3D grid point to `fname`, in the same output |
|             |                    |                                                        |
|             |                    | format as an APBS dx file.                             |
+-------------+--------------------+--------------------------------------------------------+
| 3dmap       | `energyforce`      | Will write the results of electrostatics calculations  |
|             |                    |                                                        |
|             | `outfilename`      | for points on the surface of molecules of the system   |
+-------------+--------------------+--------------------------------------------------------+
| gridct      | `ct`               | `ct` is the number of 2D grids to output.              |
+-------------+--------------------+--------------------------------------------------------+
| grid2d      | `fname` `axis`     | Set attributes of a grid output where `idx` is the     |
|             |                    |                                                        |
|             | `val`              | integer id of this grid, which can be 1 to `ct` (above)|
|             |                    |                                                        |
|             |                    | Write output of calculations for a cross section       |
|             |                    |                                                        |
|             |                    | along `axis` (*x*, *y*, or *z*) at `value`             |
+-------------+--------------------+--------------------------------------------------------+


Brownian Dynamics Simulations
-----------------------------

Will perform a brownian dynamics simulation
for the given system. The user may select and 
combine a variety of termination conditions 
depending on the system and desired results.

+-------------+--------------------+--------------------------------------------------------+
| *Keyword*   |  *Parameters*      |  *Description*                                         |
|             |                    |                                                        |
+=============+====================+========================================================+
| runtype     | `dynamics`         | Will perform a brownian dynamics run. A directory where|
|             |                    |                                                        |
|             | `outname` `ntraj`  | trajectory information will be stored in and the number|
|             |                    |                                                        |
|             |                    | of trajectories is required.                           |
+-------------+--------------------+--------------------------------------------------------+
|  termct     | `ct`               | `ct` is the number of termination conditions.          |
+-------------+--------------------+--------------------------------------------------------+
|  termcombine| `andor`            | How termination conditions will be combined. `andor`   |
|             |                    |                                                        |
|             |                    | should be *and* or *or*. Default is *or*.              |
+-------------+--------------------+--------------------------------------------------------+
|  term       | `idx` `type` `val` | Set attributes of a termination condition where `idx`  |
|             |                    |                                                        |
|             | `mols`             | is the integer id of this condition, which can be 1 to |
|             |                    |                                                        |
|             |                    | `ct` (above). `type` can be *time*,  *x<=*, *y<=*,     |
|             |                    |                                                        |
|             |                    | *z<=*, or *r<=* (or the *>=* equivalents), `val`       |
|             |                    |                                                        |
|             |                    | is the value where the simulation terminates. `mols`   |
|             |                    |                                                        |
|             |                    | is a whitespace-delimited list of molecule indices that|
|             |                    |                                                        |
|             |                    | this condition applies to (*time* requires 0, and all  |
|             |                    |                                                        |
|             |                    | else require 1).                                       |
+-------------+--------------------+--------------------------------------------------------+
|  term `idx` | `confile`          | Set attributes of contact termination condition, where |
|             |                    |                                                        |
|  contact    |                    | `idx` is the integer id of this condition, `confile`   |
|             |                    |                                                        |
|             |                    | is a path to a file containing the contact information.|
|             |                    |                                                        |
|             |                    | See below for more info.                               |
+-------------+--------------------+--------------------------------------------------------+
|  xyz        | `idx` `trajidx`    | `idx` is the molecule index for this xyz file.         |
|             |                    |                                                        |
|             | `fpath`            | Provide input XYZ file at `fpath`. For the             |
|             |                    |                                                        |
|             |                    | dynamics run, a starting configuration is              |
|             |                    |                                                        |
|             |                    | needed for each trajectory for all the molecule        |
|             |                    |                                                        |
|             |                    | types, so there should be `ntraj` xyz lines for        |
|             |                    |                                                        |
|             |                    | each molecule, the trajectory number denoted by        |
|             |                    |                                                        |
|             |                    | `trajidx`.                                             |
+-------------+--------------------+--------------------------------------------------------+



Other input files
-----------------


PQR File
^^^^^^^^
All the options above require a *PQR* file name. A PQR file 
can be generated from a PDB file using the PDB2PQR program, 
available as a web server or for download at: 

| http://nbcr-222.ucsd.edu/pdb2pqr/
| http://www.poissonboltzmann.org/docs/pdb2pqr-installation/ 

|

It may also be formatted manually. The general format of a PQR 
file is as follows, and is whitespace-delimited: 

``recName  serial  atName  resName  chainID  resNum  X  Y  Z  charge rad``

===============  ==========================================================
Parameter        Description
===============  ==========================================================
``recName``      A string that should either be ATOM or HETATM.
---------------  ----------------------------------------------------------
``serial``       An integer that provides the atom index 
---------------  ----------------------------------------------------------
``atName``       A string that provides the atom name.
---------------  ----------------------------------------------------------
``resName``      A string that provides the residue name. 
---------------  ----------------------------------------------------------
``chainID``      An optional string that provides the chain ID of the atom.
---------------  ----------------------------------------------------------
``resNumber``    An integer that provides the residue index.
---------------  ----------------------------------------------------------
``X Y Z``        Three floats that provide the atomic coordinates.
---------------  ----------------------------------------------------------
``charge``       A float that provides the atomic charge (in electrons). 
---------------  ----------------------------------------------------------
``Rad``          A float that provides the atomic radius (in A).
===============  ==========================================================



XYZ File
^^^^^^^^

The *XYZ* file simply specifies the desired molecule 
centers for a given molecule type. 

| ``mol1X  mol1Y  mol1Z``
| ``mol2X  mol2Y  mol2Z``
| ``mol3X  mol3Y  mol3Z``

|

Translation/Rotation File
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Translation/Rotation* Instead of a XYZ file, one can input a file 
specifying the translations and rotations that should be applied to
each molecule of a particular type. For these files, we follow 
the PDB standard for rotation matrices and translation vectors,
which is as follows: 

| ``mol1 rot_1_11 rot_1_12 rot_1_13 trans_1_1``
| ``mol1 rot_1_21 rot_1_22 rot_1_23 trans_1_2``
| ``mol1 rot_1_31 rot_1_32 rot_1_33 trans_1_3``
| ``mol2 rot_2_11 rot_2_12 rot_2_13 trans_2_1``
| ``mol2 rot_2_21 rot_2_22 rot_2_23 trans_2_2``
| ``mol2 rot_2_31 rot_2_32 rot_2_33 trans_2_3``

|

where ``mol1`` and ``mol2`` are indices of the molecule of 
the type this file applies to, ``rot_i_jk`` is the ``j,k`` index
of the rotation matrix for molecule ``i`` and ``trans_i_j`` 
is the ``j`` th element in the translation vector for molecule ``i``.


Contact File
^^^^^^^^^^^^

*Contact* files describe contacts between two molecular types. 
Generally this information is used to determine if a dynamics
simulation should be terminated (e.g. terminate a simulation after two 
proteins have docked). The contact file contains lines with the format: 

``moltype1  at1 moltype2 at2 dist``

where ``moltype1`` and ``moltype2`` are indices of the 
molecular types, ``at1`` is the index of an atom from the first
molecular type, ``at2`` is the index of an atom from the second 
molecular type and ``dist`` is the maximum distance between
the two atoms that defines the contact.  Note that because of the
coarse-graining, the program will identify the CG sphere closest
to the contact atom, and use the surface-to-surface distance of those 
CG spheres to compare against the reported ``dist``.


Vertex/Surface File
^^^^^^^^^^^^^^^^^^^^^^

As part of the coarse-graining process a definition of the molecular
surface is necessary. For this we have historically used the program
MSMS_ by M. Sanner, or on the online web server_

.. _MSMS: http://mgltools.scripps.edu/packages/MSMS

.. _server: http://mgl.scripps.edu/people/sanner/html/msms_server.html

If using the command line tool, after downloading it for the correct platform, 
it can be run as follows on the command line. It requires an xyzr file as input, which
is the xyz coordinates of each atom of the system followed by the vDW radius. This
information can all be found in the PQR file.

``./msms.system -if [filename].xyzr -of [outfile]``

This will produce a \*.face file and a \*.vert file, of which the \*.vert is needed. 
The vertex file is given as follows: 

.. code-block:: bash

    1669      95  3.00  1.50
   2.965    12.871    -1.084    -0.751    -0.636    -0.175       0      81  2
   3.241    11.952    -0.817    -0.936    -0.024    -0.353       0      69  2
   3.026    11.791    -0.439    -0.792     0.084    -0.604       0      79  2
   4.481    14.391    -3.026    -0.879    -0.246    -0.409       0      73  2
   5.413    15.674    -0.948    -0.337     0.499     0.798       0      73  2
   4.478    15.093    -0.297     0.286     0.886     0.365       0      81  2
   4.930    15.004    -0.240    -0.015     0.945     0.326       0      71  2
   4.072    13.663     0.763    -0.465     0.242     0.852       0      71  2

Where the first line is the number of vertex points, followed by information 
on the density of the surface, and the lines that follow indicate the cartesian 
locations of each vertex point, followed by the unit norm of the surface. 
This vertex file is used to coarse-grain the molecule.

