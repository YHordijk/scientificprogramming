:orphan:
 

star(EXPECTATION VALUES)

This section gives directives for the calculation of expectation values.

keyword(OPERATOR)

Specification of general one-electron operators.

See the :ref:`one_electron_operators` section for more information and explicit
examples.

**Advanced options**

keyword(ORBANA)

Analyze contributions to an expectation value from the individual
orbitals.

keyword(PROJEC)

Decompose an expectation value in contribution from fragment orbitals.
This action presupposes that a projection analysis has been performed,
see :ref:`*PROJECTION`
section for
details.

**Programmers options**

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0


