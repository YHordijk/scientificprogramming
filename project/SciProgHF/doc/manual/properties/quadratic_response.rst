:orphan:
 

star(QUADRATIC RESPONSE)

This section gives directives for the calculation of quadratic response
functions :cite:`Saue2002a`.


General control statements
==========================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

Definition of the quadratic response function
=============================================

keyword(DIPLEN)

Specification of dipole operators for A, B, and C
(see :ref:`one_electron_operators` for details).

keyword(A OPERATOR)

Specification of the A operator
(see :ref:`one_electron_operators` for details).

keyword(B OPERATOR)

Specification of the B operator
(see :ref:`one_electron_operators` for details).

keyword(C OPERATOR)

Specification of the C operator
(see :ref:`one_electron_operators` for details).


keyword(B FREQ)

Specify frequencies of operator B.

*Example:* 3 different frequencies.

::

    .B FREQ
     3
     0.001
     0.002
     0.01

*Default:* Static case.

::

    .B FREQ
     1
     0.0

keyword(C FREQ)

Specify frequencies of operator C (see :ref:`QUADRATIC_RESPONSE_.B FREQ`).

keyword(ALLCMB)

Evaluate all nonzero quadratic response functions and thereby
disregarding analysis of overall permutational symmetry.

*Default:* Evaluate only unique, nonzero, response functions.

Excited state properties
========================

**This page describes unreleased functionality. The keywords may not be
available in your version of DIRAC.**

First order properties of excited states can be computed from the
quadratic response function.

keyword(EXCPRP)

Give the number of "left" and "right" states in each boson symmetry.

*Example*:

::

    .EXCPRP
    5 5 5 5
    0 0 0 0

Compute the excited state expectation values \|\\langle i\|\\hat{A}\|i\\rangle|,
where i goes from 1 to 5 in each symmetry (four symmetries in this case). The
zeros can be substituted with positive integers to generate transition state moments \|\\langle i\|\\hat{A}\|j\\rangle|.

Control variational parameters
==============================

keyword(SKIPEE)

Exclude all rotations between occupied positive-energy and virtual
positive-energy orbitals.

keyword(SKIPEP)

Exclude all rotations between occupied positive-energy and virtual
negative-energy orbitals.

Control reduced equations
=========================

keyword(MAXITR)

Maximum number of iterations.

*Default:*

::

    .MAXITR
     30

keyword(MAXRED)

Maximum dimension of matrix in reduced system.

*Default:*

::

    .MAXRED
     100

keyword(THRESH)

Threshold for convergence of reduced system.

*Default:*

::

    .THRESH
     1.0D-5

Control integral contributions
==============================

The user is encouraged to experiment with these options since they may
have an important effect on run time.

keyword(INTFLG)

Specify what two-electron integrals to include
(default: :ref:`HAMILTONIAN_.INTFLG` under :ref:`**HAMILTONIAN`).

keyword(CNVINT)

Set threshold for convergence before adding SL and SS integrals to
SCF-iterations.

*2 (real) Arguments:*

::

    .CNVINT
     CNVXQR(1) CNVXQR(2)

*Default:* Very large numbers.

keyword(ITRINT)

Set the number of iterations before adding SL and SS integrals to
SCF-iterations.

*Default:*

::

    .ITRINT
     1 1

Control trial vectors
=====================

keyword(XQRNRM)

Normalize trial vectors. Using normalized trial vectors will reduce
efficiency of screening.

*Default:* Use un-normalized vectors.

Advanced/debug flags
====================

keyword(NOPREC)

No preconditioning of initial trial vectors.

*Default:* Preconditioning of trial vectors.

keyword(RESFAC)

New trial vector will be generated only for variational parameter
classes whose residual has a norm that is larger than a fraction
1/RESFAC of the maximum norm.

*Default:*

::

    .RESFAC
     1000.0

