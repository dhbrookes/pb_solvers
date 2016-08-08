# Scripts

Here is a list of files included in the scripts directory for use with PB-[S]AM code.

## Creating PB-[S]AM inputs

In the `setup` directory, there are the following scripts to
help generate inputs for PB-[S]AM. These scripts are rough,
but they do help.

### make_starting_xyz.py

For BD trajectories, you will probably want to run multiple times to
generate good statistics. We have added some options for generating 
multiple input xyz files for a BD NAM and a MFPT simulation.

There are two types of starting configurations to be generated, 
but both are called in a similar way:

For sample_grid, the program will create many different random configurations
in a cubic grid for an MFPT simulation. It is invoked as follows:

`python make_starting_xyz.py sample_grid [ntraj] [avg_rag] [xyz_name] [ntypes] [typ1_ct] ... [typn_ct]`

Where [ntraj] is the number of xyz files you want generated,
[xyz_name] is the prefix of the xyz files to be written.
They will be terminated by `_typeno_trajno_.xyz`.
[ntypes] is the number of different molecule types,
followed by [typ1_ct] [typ2_ct] ... [typn_ct] the number
of molecules of each type. [avg_rad] is the average radius of
the system. 

For random_rad, which is used in a NAM simulation, a random position
is generated for one molecule at a distance [rad] from the origin.

`python make_starting_xyz.py random_rad [ntraj] [rad] [xyz_name]`

Where [ntraj] is the number of xyz files to be printed. [rad]
is the distance from the origin that this point will be, and 
[xyz_name] is the prefix of the xyz file to be generated.

### pdb_rot_trans.py

Given a list of rotational and translational matrices (taken from REMARK 350
of a PDB file, rotate and translate a PDB or a PQR.

`python pdb_rot_trans.py [PDB or PQR fname] [Rot and trans fname]`

Where the rotational, translational file contains 3 rows for each 
repeat unit, each line containing four columns, three rotational 
columns followed by one translational.

The program outputs a pdb or pqr with the name [PDB or PQR fname].rottrans
The input PDB or PQR must be called *.pqr or *.pdb


## Plotting scripts

In the directory `plotting` there are scripts to generate 2 and 3D plots
in python. 

### plot_3D_surf.py

This file will print a 3D representation of the potential at the surface of
a molecule, given a PB-[S]AM generated map file.

`python plot_3D_surf.py [3D map name] [figure_name]`

It take as inputs the map file name, followed by the
prefix of the desired figure name. It will print out
two images of the molecule. If you want to look at the 3D plot and generate
other images, change line number 93 from `plt.close()` to `plt.show()`

### plot_2D_potential.py

This file will print a 2D representation of the potential of the system
a molecule, given a PB-[S]AM generated 2D potential file.

`python plot_2D_potential.py [2D pot name] [figure_name]`

It take as inputs the map file name, followed by the
prefix of the desired figure name. It will print out
a 2D plot of the potential at a given cross-section.

### plot_dx_to_2D_potential.py

Take a dx file and generate a 2D plot from it.

`python plot_dx_to_2D_potential.py [dx file] [2D outfile] [x/y/z CX] [CX loc] [center of mass] [CG molecule radius]`

The [dx file] and the [2D outfile] are pretty self-explanatory. 
[x/y/z CX] is the cartesian dimension in which you want the cross section, and
[CX loc] is the value in Angstroms for this cross-section.
[Center of mass] is the location of the center of mass in Angstroms, separated by 
spaces. Finally, the [CG molecule radius] is the radius of the molecule coarse-grained
from the PB-AM run.

## Analysis

After the systems are run, other analysis beyond plots may be performed. For now they 
are as follows:


