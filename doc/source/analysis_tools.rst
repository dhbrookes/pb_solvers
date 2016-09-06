
Analysis Tools
==============

Viewing PQR in VMD
------------------

First, open VMD and load the PQR file by selecting File > New Molecule on the toolbar and entering the path to the PQR file in the Filename section. Then select Graphics > Representations on the toolbar.  In the Selected Atoms selection of the pop-up, type::

 not name X

This tells VMD to remove the coarse grained spheres from the current representation.  Change the Coloring Method to Charge and the Drawing Method
to VDW. Now, to add the coarse grained spheres back into the image, select the Create Rep button and in the Selected Atoms section type::

 name X

Change the Drawing Method of the new representation to VDW and the Material to Transparent. You should now see a collection of charges inside
transparent coarse grained spheres. The final view of the Graphical Representations screen and the ouput image after this process are shown below, for both PB-AM and PB-SAM.

.. image:: images/vmd_graph_rep_pqr.png
   :width: 30%

.. image:: images/vmd_4sp.png
   :width: 66%

PB-SAM DX:

.. image:: images/vmd_barn_bars_rep.png
   :width: 35%

.. image:: images/vmd_barn_bars.png
   :width: 56%


Viewing Electrostatics in VMD
-----------------------------

To view the electrostatic results, first follow the steps above to load the PQR file. Then load the .dx file by selecting File > New Molecule on the toolbar, using the Load Files For toggle to select the previously loaded PQR and then entering the path to the .dx file in the Filename section. We will now use the .dx file to draw isosurfaces representing the surface in which the system has a selected charge. First, open the Graphical Representation screen again and select the Create Rep button. Now change the Drawing Method to Isosurface. Now in the new Draw toggle, select Solid Surface. Move the Value bar to change the charge that the isosurface represents. Change the Coloring Method to ColorID and then select the color by entering a ColorID number
in the box that appears next to the Coloring Method toggle. You may also toggle through the options by expanding this box. You may add an arbritary number
of isosurfaces by again pressing Create Rep and choosing a new Value and ColorID. The final view of the Graphical Representations screen and the ouput image after this process are shown below, for both PB-AM and PB-SAM

PB-AM DX: 

.. image:: images/vmd_graph_rep_dx.png
   :width: 30%

.. image:: images/vmd_4sp_dx.png
   :width: 66%

PB-SAM DX:

.. image:: images/vmd_barn_bars_rep_dx.png
   :width: 35%

.. image:: images/vmd_barn_bars_dx.png
   :width: 56%

For PB-SAM, in addition to using dx files for viewing ESP isosurfaces, they can also be used to 
visualize the potential on the molecule surfaces. An excellent tutorial for this 
procedure is given here_, starting at step 15.

.. _here: http://www.poissonboltzmann.org/examples/comp_tut/#VMD

.. image:: images/vmd_barn_bars_rep_surf.png
   :width: 25%

.. image:: images/vmd_barn_bars_surf1.png
   :width: 36%

.. image:: images/vmd_barn_bars_surf2.png
   :width: 36%


2D ESP Plots
------------

If one chooses the 'electrostatic' run type in PB-[S]AM, then a file for each specified 
cross section will be output with the following format::

	# Data from PBAM Electrostat run
	# My runname is barnase.x.0.dat
	units kT
	grid 200 200
	axis x 0
	origin -51.2204	-51.2204
	delta 0.512204	0.512204
	maxmin	1.35396	0
		0.2352107	0.2360552	0.2368904	0.2377159

In our :code:`scripts` directory we have provided a python script for plotting these potential
files. Simply call::

	$ python plot_2D_potential.py	<file_name>	<out_name>

where :code:`<file_name>` is an output file with format detailed above and :code:`<out_name>` 
is the desired prefix of the output file. This script creates a JPG file of the cross sectioned 
potential. An example of this output is shown below.

PB-AM 2D plot:

.. image:: images/pot_x_0.jpg
   :width: 70%
   :align: center

PB-SAM 2D plot:

.. image:: images/pbsam_barstar.png
   :width: 70%
   :align: center

3D ESP plots
------------

Additionally, if one chooses the 'electrostatic' run type in PB-[S]AM, then a .dx file will 
be created that maps the electrostatic potential to a specified number of grid points. This 
.dx file has the following format::

	# Data from PBAM Electrostat run
	# My runname is electro.map.out and units kT
	grid 10 10 10
	origin -4 -9 -9
	delta 0.8 1.8 1.8
		0.00000	0.00000	-2.90000	-5.899581

In our :code:`scripts` directory we have provided a python script for plotting these potential
files. Simply call::

	$ python plot_3d_surf.py <dx_file> <out_name> 

which will create a number of 3D potential JPG images from various viewing angles.  An example
output image is shown below. 

PB-AM 3D plot:

.. image:: images/4sp_surf.jpg
   :width: 70%
   :align: center

PB-SAM 3D plot:

.. image:: images/barstar_180.png
   :width: 70%
   :align: center

