:orphan:

.. _DIRRCI:

==========================
DIRRCI -- Direct CI module
==========================

*Relativistic RASCI module* written by Lucas Visscher

This section allows for Kramers unrestricted relativistic Configuration Interaction
calculations, Ref. :cite:`MOLFDIR`.

Note that when the namelist :ref:`goscip`
is also present in
the input, **GOSCIP**  will be called
instead of **DIRRCI** !

The input should be given in namelist form.


&CIROOT -- Select the state to converge on
==========================================

IREPNA
------

Abelian symmetry group of the desired state(s).

*Default: First abelian symmetry in the list*

NROOTS
------

Number of states to optimize on.

*Default:*

::

    NROOTS=1


**Advanced options**

ISTART
------

Start vector method.

COSCI start vectors:

::

    ISTART=1

Determinant with lowest eigenvalue:

::

    ISTART=2

First (reference) determinant:

::

    ISTART=3

If the COSCI start vectors are not available, then the default will be
2.

NSEL(NROOTS)
------------

Rank number of states to be optimized.

*Default:*

::

    NSEL=1,2,3,...

SELECT
------

Select wave functions on basis of largest overlap with the start wave
function.

*Default:*

::

    SELECT=F


&DIRECT -- Convergence control
==============================

CONVERE
-------

Convergence threshold for energy.

*Default:*

::

    CONVERR=1.0D-9

MAXITER
-------

Maximum number of direct CI iterations.

*Default:*

::

    MAXITER=10

**Advanced options**

CONVERR
-------

Convergence threshold for residual vector.

*Default:*

::

    CONVERR=1.0D-10

RESTART
-------

Restart on CI-vectors present in MRCFINV.

*Default:*

::

    RESTART=F

CPUMAX
------

Maximum amount of CPU-seconds to be used.

*Default:*

::

    CPUMAX=604800


&LEADDET -- Analyze the CI wave function
========================================

**Advanced options**

GETDET
------

Get the list of dominant determinants.

*Default:*

::

    GETDET=T

COMIN
-----

Print contributions of determinants only if the square of the
coefficients is larger than COMIN.

*Default:*

::

    COMIN=0.1


&OPTIM -- Fine tuning of the algorithm
======================================

**Programmers options**

IGENEX
------

Write coupling coefficients to file (default):

::

    IGENEX=2

Calculate coupling coefficients when needed:

::

    IGENEX=1


&RASORB -- Specify the type of CI and the active space
======================================================

NELEC
-----

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELEC=0

NRAS1
-----

Number of spinors in the RAS1 space for each abelian irrep.

*Default:*

::

    NRAS1=NSYMRP*0

NRAS2
-----

Number of spinors in the RAS2 space for each abelian irrep.

*Default:*

::

    NRAS2=NSYMRP*0

MAXH1
-----

Maximum number of holes in RAS1 spinors.

*Default:*

::

    MAXH1=0

MAXE3
-----

Maximum number of electrons in RAS3 spinors.

*Default:*

::

    MAXE3=0


&CIFOPR -- Specify different options in CI Property Module
==========================================================

PROPER
------

Activate the property module to evaluate expectation values of :ref:`one_electron_operators` defined inside

**PRPTRA** under **MOLTRA**, over the DIRRCI wavefunction :cite:`Nayak2006`, :cite:`Nayak2007`, :cite:`Nayak2009`.

*Default:*

::

     PROPER=F

NEOPER
------

Define the number of property operators one needs to calculate expectation values.


*Default:*

::

     NEOPER=0

NAMEE
-----

Mention the names of :ref:`one_electron_operators` as defined inside **PRPTRA**.

For example, Electric dipole moment and magnetic hyperfine structure constants

can be defined as follows. Of coures, the given names are users own choice.

::

     NAMEE='Z-DIP','X1-HYP','Y1-HYP','Z1-HYP'
