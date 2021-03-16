Geometry optimization in xyz coordinates
========================================

Acetone
-------
The testing example - the acetone molecule,
with the input geometry in the xyz format,
gives identical (within set up limits of accuracy) geometries and
final energies with and without imposed symmetry for 
the xyz-input reading. These two energies are identical with 
the resulting energy from the geometry optimization with the mol-input file
with the symmetry (see the neighbouring test geom_opt).

Total energies of acetone after geom.optimizations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
-189.536032850805 a.u. (xyz; with symmetry C(2v) -> C2v)
-189.536032850926 a.u. (xyz; no symmetry)
-189.536032850804 a.u. (mol; C2v symmetry; test geom_opt)

Water
-----
The water molecule has the starting geometry of a bent molecule.
If it would be the linear starting structure, the geometry optimization
would fail.

Again, both symmetry imposed and symmetry not imposed results give identical
converged geometries (in limits of given accuracy) and the total energies.
