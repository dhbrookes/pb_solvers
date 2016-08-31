
PB-AM Inputs
============

Here are the inputs that are used in a PB-AM run, and
their meanings. The general format of the input is:

``keyword <value>``

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
| runname     | `<name>`           | `name` is desired internal name of this run.           |
+-------------+--------------------+--------------------------------------------------------+
| attypes     | `<numtypes>`       | Set the number of different atom                       |
|             |                    |                                                        |
|             |                    | types to `numtypes`                                    |
+-------------+--------------------+--------------------------------------------------------+
| pqr         | `<idx>`  `<fpath>` | The molecule index for this xyz file.                  |
|             |                    |                                                        |
|             |                    | Provide input PQR file at `fpath`                      |
+-------------+--------------------+--------------------------------------------------------+
| xyz         | `<idx>`  `<fpath>` | The molecule index for this xyz file.                  |
|             |                    |                                                        |
|             |                    | Provide input XYZ file at `fpath`                      |
+-------------+--------------------+--------------------------------------------------------+
|  transrot   | `<idx>`  `<fpath>` | The molecule type index for this translation/rotation  | 
|             |                    |                                                        |
|             |                    | file that can be input instead of a xyz file.          |
|             |                    |                                                        |
|             |                    | Provide input file at `fpath`                          |
+-------------+--------------------+--------------------------------------------------------+
|  randorient |                    | If you want your molecules to be randomly              |
|             |                    |                                                        |
|             |                    | rotated, use this flag                                 |
+-------------+--------------------+--------------------------------------------------------+
|  units      | `<units>`          | Set the units of output to `units`.                    |
|             |                    |                                                        |
|             |                    | The current options are: `jmol` (Joules/mole),         |
|             |                    |                                                        |
|             |                    | `kT` (kT/e) and `kcalmol` (kCal/mole).                 |
+-------------+--------------------+--------------------------------------------------------+
|  salt       | `<con>`            | Set salt concentration in the system to `con`          |
+-------------+--------------------+--------------------------------------------------------+
|  temp       | `<T>`              | Set system temperature to `T`                          |
+-------------+--------------------+--------------------------------------------------------+
|  idiel      | `<ival>`           | Set the interior dielectric constant to `ival`         |
+-------------+--------------------+--------------------------------------------------------+
|  sdiel      | `<sval>`           | Set the solvent dielectric constant to `sval`          |
+-------------+--------------------+--------------------------------------------------------+
|  pbc        | `<boxlength>`      | Set size of periodic box to `boxlength`                |
+-------------+--------------------+--------------------------------------------------------+
|  random     | `<seed>`           | Seed the random number generator with `seed`           |
+-------------+--------------------+--------------------------------------------------------+
|  type       | `<idx>` `<ct>`     | Set attributes of an atom type, where `idx` is the     | 
|             |                    |                                                        |
|             | `<mvtype>` `<dtr>` | integer id of this type, which can be 1 to `numtypes`  |
|             |                    |                                                        |
|             | `<drot>`           | (above). `ct` is the number of atoms of this type      |
|             |                    |                                                        |
|             |                    | in the system and `mvtype` describes the way this      |
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
| runtype     | `<energyforce>`    | This will invoke the calculation of energies, forces   |
|             |                    |                                                        |
|             | `<outfilename>`    | and torques for the input system.                      |
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
| runtype     | `<electrostatics>` | Will write the results of electrostatics calculations  |
|             |                    |                                                        |
|             | `<gridpts>`        | for every 3D grid point to `fname`, in the same output |
|             |                    |                                                        |
|             |                    | format as an APBS dx file.                             |
+-------------+--------------------+--------------------------------------------------------+
| dx          | `<energyforce>`    | Will write the results of electrostatics calculations  |
|             |                    |                                                        |
|             | `<outfilename>`    | for every 3D grid point to `fname`, in the same output |
|             |                    |                                                        |
|             |                    | format as an APBS dx file.                             |
+-------------+--------------------+--------------------------------------------------------+
| 3dmap       | `<energyforce>`    | Will write the results of electrostatics calculations  |
|             |                    |                                                        |
|             | `<outfilename>`    | for points on the surface of molecules of the system   |
+-------------+--------------------+--------------------------------------------------------+
| gridct      | `<ct>`             | `ct` is the number of 2D grids to output.              |
+-------------+--------------------+--------------------------------------------------------+
| grid2d      | `<fname>` `<axis>` | Set attributes of a grid output where `idx` is the     |
|             |                    |                                                        |
|             | `<val>`            | integer id of this grid, which can be 1 to `ct` (above)|
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
| runtype     | `<dynamics>`       | Will perform a brownian dynamics run. A directory where|
|             |                    |                                                        |
|             | `<outnm>` `<ntraj>`| trajectory information will be stored in and the number|
|             |                    |                                                        |
|             |                    | of trajectories is required.                           |
+-------------+--------------------+--------------------------------------------------------+
|  termct     | `<ct>`             | `ct` is the number of termination conditions.          |
+-------------+--------------------+--------------------------------------------------------+
|  termcombine| `<andor>`          | How termination conditions will be combined. `andor`   |
|             |                    |                                                        |
|             |                    | should be *and* or *or*. Default is *or*.              |
+-------------+--------------------+--------------------------------------------------------+
|  term       | `<idx>` `<type>`   | Set attributes of a termination condition where `idx`  |
|             |                    |                                                        |
|             | `<val>` `<mols>`   | is the integer id of this condition, which can be 1 to |
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
|  term `idx` | `<confile>` `<pad>`| Set attributes of contact termination condition, where |
|             |                    |                                                        |
|  contact    |                    | `idx` is the integer id of this condition, `confile`   |
|             |                    |                                                        |
|             |                    | is a path to a file containing the contact information,|
|             |                    |                                                        |
|             |                    | and `pad` specifies a correction for the case when the |
|             |                    |                                                        |
|             |                    | contact distance cannot be reached due to the spherical|
|             |                    |                                                        |
|             |                    | assumption of the model. See below for more info.      |
+-------------+--------------------+--------------------------------------------------------+
|  xyz        | `<idx>` `<trajidx>`| `idx` is the molecule index for this xyz file.         |
|             |                    |                                                        |
|             | `<fpath>`          | Provide input XYZ file at `fpath`. For the             |
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

.. code-block:: bash

    mol1X  mol1Y  mol1Z
    mol2X  mol2Y  mol2Z
    mol3X  mol3Y  mol3Z

|

Translation/Rotation File
^^^^^^^^^^^^^^^^^^^^^^^^^^

*Translation/Rotation* Instead of a XYZ file, one can input a file 
specifying the translations and rotations that should be applied to
each molecule of a particular type. For these files, we follow 
the PDB standard for rotation matrices and translation vectors,
which is as follows: 

.. code-block:: bash

    mol1 rot_1_11 rot_1_12 rot_1_13 trans_1_1  
    mol1 rot_1_21 rot_1_22 rot_1_23 trans_1_2  
    mol1 rot_1_31 rot_1_32 rot_1_33 trans_1_3  
    mol2 rot_2_11 rot_2_12 rot_2_13 trans_2_1  
    mol2 rot_2_21 rot_2_22 rot_2_23 trans_2_2  
    mol2 rot_2_31 rot_2_32 rot_2_33 trans_2_3  

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
the two atoms that defines the contact.  Note that sometimes 
these distances cannot be reached due to the assumption in this model that
the molecule is spherical. To correct for this case, one must 
specify a "pad"  distance that is defined as the maximum distance between
the radial projections of the atoms onto the surface of their 
respective spheres that defines a contact.
