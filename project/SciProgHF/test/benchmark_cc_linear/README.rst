Coupled Cluster Calculations
============================

This test calculates the energy of the doubly ionized nitrogen molecule  using Coupled-Cluster (CC) theory.
Note that we consistently leave out the SS integrals.

The nitrogen molecule, N2(2+), is closed shell system. All electrons are taken for the correlation.

The test is meant to check the automatic detection of symmetry AND the implementation of linear symmetry.
If the homunuclear diatomic molecule is not in the position, 
from which the program can find out the the symmetry, and we want
to impose the linear symmetry, we have to add the "Dinfh" label into molecule input file.

We assume that the employed correlation consistent basis sets (ccpVDZ, ccpVQZ, ccpVTZ) are available
in your standard path. These sets can be
obtained from the DALTON distribution and may later be included in the DIRAC distribution as well.


