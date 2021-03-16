Molecular mean-field Hamiltonian schemes
========================================

This test calculates the spin-orbit splitting in the ground state of the FO molecule 
using various molecular-mean-field models to the full Dirac-Coulomb-Gaunt Hamiltonian.

This test shows the use of 4 of the approximate Hamiltonian models 
to the full 4c-Dirac-Coulomb-Gaunt Hamiltonian (which is at present not available in Dirac at the post-HF level)
which are described in the paper [Sikkema2009]_ .
 
Following the nomenclature given in Table I in the above paper the Hamiltonian models tested here are:

::

 ^4DCG*          ==> input name: 4dcg_x.inp
 ^4DCG**         ==> input name: 4dcg_xx.inp
 ^2DCG^M(U_mol)  ==> input name: 2dcg_mmf.inp
 ^2DCG^A(U_nuc)  ==> input name: 2dcg_x2c.inp


