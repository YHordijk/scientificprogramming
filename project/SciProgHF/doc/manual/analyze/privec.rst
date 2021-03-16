:orphan:
 

star(PRIVEC)

Print vectors read from DFCOEF.


keyword(VECPRI)

For each fermion irrep, give an :ref:`orbital_strings` of orbitals to print.

*Default:* Print the occupied electronic solutions.


keyword(PRICMP)

Separate control for printing of large and small components.

*Default:* Print large components.

::

    .PRICMP
     1 0

Print only small components:

::

    .PRICMP
     0 1


keyword(AOLAB)

Print vectors in AO basis.

*Default:* Print vectors in SO basis.


keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

