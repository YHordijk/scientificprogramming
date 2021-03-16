:orphan:
 

.. _goscip:

======================
GOSCIP -- COSCI module
======================

General Open Shell CI Program written by Olivier Visser,  Ref. :cite:`Visser1992`.

This module allows for small full configuration interaction calculations. 
It is invoked from within :ref:`DIRRCI` if the namelist GOSCIP is present in the input, or if :ref:`*RESOLVE` is specified when doing open shell HF calculations.

The input is given in namelist form.

Historical note: This program was originally written for the :cite:`MOLFDIR` suite and was later also included in DIRAC.


&GOSCIP
=======

Specify the CI space.

NELEC
-----

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELEC=0

**NOTE: in the developer's version this line reads as**

NELACT
------

Number of electrons (excluding frozen core electrons).

*Default:*

::

    NELACT=0

**Programmers options**

IPRNT
-----

Print level.

*Default:*

::

    IPRNT=0

.. _popana:

&POPANA
=======

Analyze the CI wave function.

**Advanced options**

THRESH
------

Print only determinants with coefficients higher than THRESH.

*Default:*

::

    THRESH=1.0D-3

DEGEN
-----

Threshold for when several states are considered to be degenerate.

*Default:*

::

    DEGEN=1.0D-10

SELPOP
------

Select only states with relative energies lower than SELPOP.

*Default:*

::

    SELPOP=1.0D2

