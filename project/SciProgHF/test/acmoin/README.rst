ACMOIN test
===========

This test demonstrates the export and import of coefficients
between different symmetries.

The test considers a highly fictitious system, namely ammonia
in D_3h symmetry at c = 20.0 a.u.

The C3 axis is along the z-axis and we would like to analyze
the structure of the molecule by projection analysis using
atomic orbitals will well-defined mj-values with respect to
this axis. This is achieved by calculating the H and N atoms in
linear symmetry, exporting the coefficients to C1 symmetry (ACMOUT)
and them importing them to Cs symmetry (ACOMIN).

To check the import of coefficients for nitrogen, the SCF calculation
is repeated in Cs, and should then converge immediately to the right
energy and orbitals.

This test is run both for spinfree and full relativistic two-component (X2c) Hamiltonians.

