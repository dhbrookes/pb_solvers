
PBAM
====

PB-AM is an analytical solution to the linearized Poisson-Boltzmann equation for multiple spherical objects of arbitrary charge distribution in an ionic solution. The solution can be reduced to a simple system of equations as follows:

\\[ A = \\Gamma \\cdot (\\Delta \\cdot T \\cdot A + E) \\]

Where \\(A^{(i)}\\) represents the effective multipole expansion of the charge distributions of molecule \\(i\\). \\(E^{(i)}\\) is the free charge distribution of molecule \\(i\\). \\(\Gamma\\) is a dielectric boundary-crossing operator, \\(\Delta\\) is a cavity polarization operator, \\(T\\) an operator that transforms the multipole expansion to a local coordinate frame.  More details on the method are available in [LoHe06]_. Once \\(A^{(i)}\\) has been solved, through an iterative SCF method, physical properties of the system can be computed, as detailed in the next section.