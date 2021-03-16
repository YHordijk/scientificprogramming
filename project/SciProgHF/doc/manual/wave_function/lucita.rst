:orphan:
 

star(LUCITA)

Spin-free relativistic GASCI module written by Jeppe Olsen, adaptation to DIRAC by Timo Fleig, parallelization by Stefan
Knecht and Hans Joergen Aagaard Jensen.

LUCITA is a string-based Hamiltonian-direct configuration interaction (CI) program, based on the LUCIA code :cite:`Olsen1990`.
It is capable of doing efficient CI computations at arbitrary excitation level, e.g. FCI, SDCI, RASCI, and MRCI using
general active spaces. The code is interfaced :cite:`Fleig2005` to molecular integrals obtained in the spin-free Dirac
formalism and uses non-relativistic point group symmetry. It is implemented as a full parallel version :cite:`Knecht2008`.

A central feature of the program is the Generalized Active Space (GAS) concept, in which the underlying total orbital
space is subdivided into a basically arbitrary number (save for an upper limit) of subspaces with arbitrary occupation
constraints. This is the most general approach to orbital space subdivisions. The program uses DIRAC orbitals from either
a closed- or an open-shell calculation in a spin-free relativistic formalism or the non-relativistic Lévy-Leblond
formalism (see :ref:`**HAMILTONIAN`).

The technical limitations are roughly set by several 100 million determinants in the CI expansion on PCs and common computing clusters and several billions of determinants on supercomputers with ample memory.
The program also computes 1- and, if requested 2-particle densities from optimized CI wave functions. The density matrices may be printed along with natural orbital occupations and the corresponding eigenvectors (NOs).

**General information**
=======================

Please add the lines

::

    .SCHEME
     4

to the :ref:`**MOLTRA` input section for the integral transformation.

**Mandatory keywords**
======================

keyword(INIWFC)

specifies the initial SCF wave function. 

*closed-shell:*

::

    .INIWFC
     DHFSCF

*open-shell:*

::

    .INIWFC
     OSHSCF

keyword(CITYPE)

Type of CI calculation. Typical multi-reference CI calculations should be defined by RASCI or GASCI.

valid choices are: FCI, SDCI, SDTQ, RASCI, GASCI. If you choose RASCI or GASCI please take a look at the *RASCI/GASCI*
keyword section below.

*example:*

::

    .CITYPE
     FCI

keyword(MULTIP)

specifies the target state spin multiplicity (2S+1) of the desired eigenstates.

*example for singlet state(s) (S=0):*

::

    .MULTIP
     1

**Optional keywords**
=====================

keyword(DENSI)

Level of computed density matrices. *default:* one-elctron density matrices and natural orbital occupation numbers.
To compute one- and two-particle density matrices:

::

    .DENSI
     2

keyword(NROOTS)

Number of states(roots) to optimize. *default*:

::

    .NROOTS
     1

keyword(SYMMET)

Spatial symmetry of the desired state(s). This is the boson symmetry label referring to an irrep ordering as defined by
the group generators. *default*:

::

    .SYMMET
     1


keyword(MAXITR)

maximum number of CI iterations. *default*:

::

    .MAXITR
     100

keyword(INACTI)

Inactive orbitals per boson symmetry, separated by commas. This keyword is not allowed in connection with the citype
choice :ref:`LUCITA_.CITYPE` GASCI. *default*: all orbitals active. 
example for :math:`D_{2h}`:

::

    .INACTI
     1,1,0,0,1,0,0,2

keyword(RSTRCI)

Restart CI from vector(s) on file LUCVECT resp. LUCVECT.0 (parallel calculation; until including DIRAC v13.1):

::

    .RSTRCI
     1

*default*: no restart.

::

    .RSTRCI
     0

Convergence threshold for energy (double-precision value):

::

    .CONVER
     x.xd-0x

*default*: 1.0d-08

::

    .CONVER
     1.0d-08

star(LUCITA GASCI)

**Specific input -- GASCI**
===========================

keyword(NACTEL)

Number of active electrons. *default*: none. example:

::

    .NACTEL
     10

keyword(GASSHE)

Number and specification of GAS orbitals. Line with the number of GA spaces used (1-7), 
followed by one line per GAS with number of orbitals per boson symmetry, separated by commas. *default*: none.
example for 2 GAS spaces and active orbital distribution when running in the :math:`D_{2h}` symmetry:

::

    .GASSHE
     2
     2,0,4,4,8,2,1,1
     8,2,6,6,19,5,3,3

keyword(GASSPC)

Number and specification of sequential CI calculations. Line with the number of CI calculations with given GA spaces
(currently only 1 is allowed), followed by one line per GAS with 2 numbers each: The first gives the minimal number of
accumulated electrons after this GAS, the second the corresponding maximum number, separated by blanks (defining
occupation constraints of each GAS). *default*: None. example for an excitation pattern in connection with the 2-GAS
input from :ref:`LUCITA_GASCI_.GASSHE` above and 10 active electrons:

::

    .GASSPC
     1
     8  10
     10 10

star(LUCITA RASCI)

**Specific input -- RASCI**
===========================

keyword(NACTEL)

Number of active electrons. *default*: none. example:

::

    .NACTEL
     10

keyword(RAS1)

RAS1 specification and maximum number of holes. Line with orbitals per boson symmetry, separated by commas, 
followed by a line with the maximum number of holes in RAS1. *default*: none. 
Example for an active orbital distribution when running in the D:math:`_{2h}` symmetry and max 2 holes:

::

    .RAS1
     2,0,4,4,2,2,1,1
     2


keyword(RAS2)

RAS2 specification. Line with orbitals per boson symmetry, separated by commas. 
Example for an active orbital distribution 
when running in the D:math:`_{2h}` symmetry:


::

    .RAS2
     1,0,1,1,1,1,0,1

keyword(RAS3)

RAS3 specification and maximum number of electrons. Line with orbitals per boson symmetry, separated by commas, followed by a line with the maximum number of electrons in RAS3.
Example for an active orbital distribution when running in the D:math:`_{2h}` symmetry and max 2 electrons in RAS3:

::

    .RAS3
     2,0,4,4,2,2,1,1
     2

