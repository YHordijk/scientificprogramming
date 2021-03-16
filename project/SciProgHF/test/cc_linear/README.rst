Coupled Cluster Calculations
============================

Setting the linear symmetry
---------------------------
This test calculates the energy of the doubly ionized nitrogen molecule using Coupled-Cluster (CC) theory.
Note that we consistently leave out the SS integrals.

The nitrogen molecule, N2, is closed shell system. 
Four lowest electrons -  sigma(1s) and sigma(1s)* - are frozen and remaining 10 valence electrons
is included in the correlation.

This test is meant to check the automatic detection of symmetry AND the implementation of the linear symmetry.

If the homunuclear diatomic molecule is not in the position, 
from which the program can find out the the symmetry, and we want
to impose the linear symmetry, we have to add the "Dinfh" label into molecule input file.

We assume that the employed correlation consistent basis sets (ccpVDZ, ccpVQZ, ccpVTZ) are available
in your standard path. These sets can be
obtained from the DALTON distribution and may later be included in the DIRAC distribution as well.

X2C-Coupled Cluster (re)start
-----------------------------
The second test is dedicated for (re)starting X2C-Coupled Cluster calculations after the (X2c) SCF and MOLTRA steps.

Note that in the X2c-CC start you have to tell program not to perform picture change transformation of four-component property
operators sitting in the AOPROPER file.

