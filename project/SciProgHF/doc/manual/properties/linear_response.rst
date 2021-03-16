:orphan:
 

star(LINEAR RESPONSE)

Linear Response module written by T. Saue and H. J. Aa. Jensen :cite:`Saue2003`.


General control statements
==========================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

Definition of the linear response function
==========================================

keyword(A OPERATOR)

Specification of the A operator
(see :ref:`one_electron_operators` for details).

keyword(B OPERATOR)

Specification of the B operator
(see :ref:`one_electron_operators` for details).

keyword(OPERATORS)

Specification of both the A and B operators
(see :ref:`one_electron_operators` for details).

keyword(EPOLE)

Specification of electric Cartesian multipole operators of order L .
Specify order.

*Example:* Electric dipole operators

::
      .EPOLE
      1

keyword(MPOLE)

Specification of magnetic Cartesian multipole operators of order L .
Specify order.
*Example:* Magnetic dipole operators

::
      .MPOLE
      1


keyword(TRIAB)

Enforce triangularity of response function.

.. math::

  \langle\langle A; B \rangle\rangle_\omega = \langle\langle B; A \rangle\rangle_\omega

Only one function is calculated.

*Default:* Deactivated.

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

keyword(IMAGIN)

Employ imaginary frequencies.

keyword(ALLCMB)

Form all possible combinations
even if imaginary.

*Default:* Deactivated.

keyword(UNCOUP)

Uncoupled calculation.

Control variational parameters
==============================

keyword(OCCUP)

For each fermion ircop give an :ref:`orbital_strings`
of inactive orbitals
to include in the linear response calculation.

keyword(VIRTUA)

For each fermion ircop give an :ref:`orbital_strings`
of virtual orbitals
to include in the linear response calculation.

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
     200

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

keyword(REAXVC)

Read solution vectors from file XVCFIL

::

    .REAXVC
    XVCFIL

*Default:* No restart on solution vectors. The file has to have six
characters. Make sure there is no blank character in front of the file
name.

For a restart on solution vectors it is useful to set

::

    .REAXVC
    XVCFIL
    .ITRINT
     0 0

otherwise LS-integrals (and SS-integrals) are switched on later and one
may first iterate away and then back to a possibly converged response
vector.

Often you have a converged SCF wave function along with a response
vector. In this case make sure that

::

    **DIRAC
    #.WAVE FUNCTION

is commented out. Make then also sure that you use the DFCOEF file which
has been obtained in the *same* calculation as the response vector file.
Otherwise you may observe more response solver iterations than
necessary.

keyword(XLRNRM)

Normalize trial vectors. Using normalized trial vectors will reduce
efficiency of screening. CLARIFY!

*Default:* Use un-normalized vectors.

Analysis
========

keyword(ANALYZ)

The linear response function is obtained by contracting the property gradient :math:`\boldsymbol{E}_{A}^{[1]}` associated with property :math:`A` with the solution vector
:math:`\boldsymbol{X}_{B}` associated with property :math:`B`:

.. math::

   \langle\langle\hat{A};\hat{B}\rangle\rangle_{\omega_{b}}=\boldsymbol{E}_{A}^{[1]\dagger}\boldsymbol{X}_{B}

The indices of the two vectors run over orbital rotation indices, that is, one occupied and one virtual index. This analysis shows the most important contributions from orbital pairs to
the dot product as well as the most important occupied and virtual orbitals. More information about these orbitals can then be gleaned from Mulliken population analysis.

keyword(ANATHR)

This keyword adjusts the number of terms shown in the above analysis. The threshold is in terms of percentage value to the total value of the linear response function.

*Default:* :math:`2.0`.

Advanced/debug flags
====================

keyword(E2CHEK)

Generate a complete set of trial vector which implicitly allows the
explicit construction of the electronic Hessian. Only to be used for
small systems !

keyword(ONLYSF)

Only call FMOLI in sigmavector routine: only generate one-index
transformed Fock matrix :cite:`Saue2003`.

keyword(ONLYSG)

Only call FMOLI in sigmavector routine: 2-electron Fock matrices using
one-index transformed densities :cite:`Saue2003`.

keyword(STERNHEIM)

Set diagonal elements of orbital part of Hessian equal to

.. math::

  -2 m c^2

for rotations between occupied positive-energy and
virtual negative-energy orbitals.

*Default:* Deactivated.

keyword(STERNC)

(Sternheim complement) allows to separate basis set incompleteness from
the replacement of an inner sum over negative-energy orbitals only by
the full sum. In order to benefit from this functionality (only for
specialists !), you should run with print level 2 under properties.

Then you can do a sequence of calculations: 1) .SKIPEP 2) .STERNH 3)
.STERNC The diamagnetic contribution of 1) is the non-relativistic
expectation value, whereas 2) is the Sternheim approximation, that is
replacing orbital energy differences with

.. math::

  -2 m c^2

With no
basis set incompleteness the sum of the diamagnetic contribution 2) and
the paramagnetic contribution 3) should equal the diamagnetic
contribution of 1).

*Default:* Deactivated.

keyword(COMPRESSION)

Reduce number of orbital variation parameters by checking corresponding
elements of gradient vector against a threshold. This may reduce memory.

*Default:* No compression.

::

    .COMPRESSION
     0.0

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

