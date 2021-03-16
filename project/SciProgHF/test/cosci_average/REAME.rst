COSCI states of the V atom
===========================

i) In this test we demostrate obtaining the converged molecular spinors (MS) of the V high-spin open-shell atom in averaged
configuration.
The complicated shell structure of V(Z=23, ^4 F_3/2), [Ar]3d(3)4s(2), requires the starting set of MSs
from the V(5+) closed-shell cation and, in the following SCF step of the neutral V atom,
applying the non-dynamical overlap selection. 
Within converged molecular spinors, the averaged 3d(3) open-shell in V is 
under the 4s(2) closed shell, and both are the same fermion symmetry (E1g).

ii) Another purpose of this test is checking the averaged COSCI energy (-948.0620823505), 
which must be identical (within available decimal places) with the total SCF energy.
This is because the diagonalization performed in the COSCI step does not change
the value of the trace \Sum_i <D_i | H | D_i>, 
where the i index runs over all determinants D contained in the open-shell expansion and H is the Hamiltonian. 
If these two values differ this signifies either a bug or an approximation that is made in the integral transformation. 

