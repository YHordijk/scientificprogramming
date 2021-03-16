Energy and dipole moment
========================

This test calculates the energy and the dipole moment
of the water molecule at different levels of theory.
Note that we consistently leave out the SS integrals.

The small run is ONLY for testing pruposes : it is NOT a good approximation
to leave out SL integrals without any compensation and also the virtual space 
is too much truncated !

We have also added the so-called spin-orbit correlated calculations, where spin-free molecular
orbitals are read in and the full relativistic Hamiltonian is used at the correlated
leverl (MP2, CC).

Values of the correlated dipole moment with the Dirac-Coulomb Hamiltonian:

MP2:     2.02142292 D
SO-MP2:  2.02142759 D
