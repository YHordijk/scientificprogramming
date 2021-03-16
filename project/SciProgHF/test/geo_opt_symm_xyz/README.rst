Geometry optimization tests
===========================

Collection of geometry optimizations tests of symmetric molecules 
with xyz molecule inputs.

Testsa are aimed at checking the symmetry detection and proper geometry displacements
of simple symmetric molecules.

Resulting converged geometries must be of the given symmetry.

CH4
---

Has the T(d) point group symmetry.
Sets the D2 computational symmetry.
For getting  the converged geometry, one has to apply the .NORTSD keyword.


NH3
---

Has the C(3v) point group symmetry.
Sets the Cs computational symmetry.

Allene
------

Has the D(2d) point group symmetry.
Molecule can not detect symmetry and during optimization 
it gets wrong coordinates resulting in SCF nonconvergence.


SF5Cl
-----

Has the C(4v) point group symmetry.
Sets the C2v computational symmetry.


SF6
-----

Has the O(h) point group symmetry.
Sets the D2h computational symmetry.

