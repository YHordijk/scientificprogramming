Symmetry recognition in mol files
=================================

The test is meant to check the automatic detection of molecular symmetry in mol input files and calculates the nonrelativistic energy.

We assume that the basis sets are available in your standard path.
These sets can be obtained from the Dalton distribution and may later be included in the DIRAC distribution as well.

Note that a coarse DFT grid for LDA is used and the SCF calculation is stopped after 4 iterations.

Examples
--------
The following table contains few examples of
the symmetry detection from the xyz coordinates in mol files.


==============    =================
Molecule          Detected symmetry
==============    =================
CH4               T(d)  ->  D2
H2O               C(2v) ->  C2v 
SF6               O(h)  ->  D2h
trans CHF=CHF     C(2h) ->  C2h
==============    =================

Few more examples can be added here to demonstrate (and to test) the wider variety of symmetry detections.
