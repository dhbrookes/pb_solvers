
PB-AM
=====

PB-AM is an analytical solution to the linearized Poisson-Boltzmann 
equation for multiple spherical objects of arbitrary charge distribution 
in an ionic solution. The solution can be reduced to a simple system 
of equations as follows:

\\[ A = \\Gamma \\cdot (\\Delta \\cdot T \\cdot A + E) \\]

Where \\(A^{(i)}\\) represents the effective multipole expansion 
of the charge distributions of molecule \\(i\\). \\(E^{(i)}\\) is 
the free charge distribution of molecule \\(i\\).  \\(\\Gamma\\) is 
a dielectric boundary-crossing operator, \\(\\Delta\\) is a cavity 
polarization operator, \\(T\\) an operator that transforms the 
multipole expansion to a local coordinate frame.  More details on 
the method are available in [LoHe06]_. Once \\(A^{(i)}\\) has been 
solved, through an iterative SCF method, physical properties of the 
system can be computed, as detailed in the next section.

Physical calculations
---------------------

Interaction energies
^^^^^^^^^^^^^^^^^^^^^

From the above formulation, computation of the interaction energy 
(\\(\\Omega^{(i)}\\)) for molecule i, is given as follows:

\\[\\Omega^{(i)}=\\frac{1}{\\epsilon_s} \\sum_{j \\ne i}^N \\left \\langle  T \\cdot A^{(j)} ,  A^{(i)} \\right \\rangle \\]

Where \\(\\langle . . . \\rangle\\) denotes an inner product.

Forces and Torques
^^^^^^^^^^^^^^^^^^

When energy is computed, forces follow as:

\\[ \\textbf{F}^{(i)} = \\nabla_i \\Omega^{(i)}=\\frac{1}{\\epsilon_s} [ \\langle \\nabla_i \\,T \\cdot A^{(j)} ,  A^{(i)} \\rangle +  \\langle T \\cdot A^{(j)} ,   \\nabla_i \\, A^{(i)} \\rangle ]\\]

The method to calculate the torque \\(\\boldsymbol{\\tau}^{(i)}\\) on 
molecule is outside the scope of this manual, but is discussed extensively in [LoHe06]_
