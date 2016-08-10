
PB-SAM
===========

PB-SAM is a semi-analytical solution to the linearized Poisson-Boltzmann 
equation for multiple molecules of arbitrary charge distribution 
in an ionic solution. The solution is an extension of the analytical method,
leveraging Fast-Multipole methods as well as boundary elements. Each molecule is
coarse-grained as a system of overlapping spheres, whose surface charges are represented
by the multipole expansions \\(H^{(i)}\\) and \\(F^{(i)}\\). To solve for the potential,
the following interactions are considered:

 - Intra-molecular interactions between overlapping spheres are treated numerically
 - Intra-molecular interactions between non-overlapping spheres are treated analytically
 - Inter-molecular interactions between spheres on different molecules

With these interactions, the multipole expansions are solved with an iterative 
SCF method, briefly given as

\\[ H^{(i,k)} = I_{E}^{(i,k)} \\cdot \\left ( H^{(i,k)} + F^{(i,k)} + T \\cdot H^{(j,l)} \\right ) \\]
\\[ F^{(i,k)} = I_{E}^{(i,k)} \\cdot \\left ( H^{(i,k)} + F^{(i,k)} + T \\cdot F^{(j,l)} \\right ) \\]

For details on the method, please see [YaHe10]_ and [YaHe13]_.

Physical calculations
---------------------

Interaction energies
^^^^^^^^^^^^^^^^^^^^^

From the above formulation, computation of the interaction energy 
(\\(\\Omega^{(i)}\\)) for molecule i, is given as a sum of all the interactions
of spheres \\(k\\) within it with all external spheres (in a simplified form) as follows:

\\[\\Omega^{(i)}=\\frac{1}{\\epsilon_s} \\sum_{k \\, in\\, i} \\sum_{j \\ne i}^N \\sum_{l\\, in \\, j} \\left \\langle  T \\cdot H^{(j,l)} ,  H^{(i,k)} \\right \\rangle \\]

Where \\(\\langle . . . \\rangle\\) denotes an inner product.

Forces and Torques
^^^^^^^^^^^^^^^^^^

When energy is computed, forces follow as:

\\[ \\textbf{F}^{(i)} = \\nabla_i \\Omega^{(i)}=\\frac{1}{\\epsilon_s} [ \\langle \\nabla_i \\,T \\cdot H^{(j,l)} ,  H^{(i,k)} \\rangle +  \\langle T \\cdot H^{(j,l)} ,   \\nabla_i \\, H^{(i,k)} \\rangle ]\\]

The method to calculate the torque \\(\\boldsymbol{\\tau}^{(i)}\\) on 
molecule is outside the scope of this manual, but is discussed extensively in [YaHe13]_
