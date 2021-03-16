:orphan:

.. _DIRAC_computer_exercises:

########################
DIRAC computer exercises
########################

The following exercises are to be done in groups of 2-3 students each covering different groups/atoms of the periodic table of elements.
At the end of the day, each group should give a short presentation on their findings to the other groups.

In order to run a DIRAC calculation you need an input file of arbitrary name, e.g. :file:`dirac.inp` which will be setup in the
following and which contains the computational directives. In addition, a second file is needed providing the molecular structure.
DIRAC can handle two formats but for convenience we here focus on the standard **XYZ** molecular input file, 
named for example :file:`molecule.xyz`.
The run-command inside your cluster run script looks then as follows::

  $ pam --noarch --inp=dirac.inp --mol=molecule.xyz

where the DIRAC output is re-directed to the file :file:`dirac_molecule.out`. 
Make sure that you rename the output file with a proper name
if you want to take a look at what you did back home.

More information on running DIRAC calculations can be found in the :ref:`FirstSteps` section of this documentation.

If not stated otherwise we will always use a Gaussian nuclear model in DIRAC.

*****************
Exercises - Day 3
*****************

In the following exercises we will do our (first) calculations in DIRAC looking at the effect on the total energy and the
valence-shell spin-orbit splitting of a heavy-element atom while picking a different nuclear
model, varying the actual speed of light as well as going to the non-relativistic limit.

::

   the contender: Pb atom


A general :file:`.xyz` file (see :ref:`molecule_using_xyz`) for an atom reads (replace here XX with Pb)::

  1
  XX atom
  XX 0.0 0.0 0.0


Atomic basis set calculations
=============================

In order to get started we calculate the ground state configuration and lowest excited states 
within the p :math:`^2` manifold ( :math:`^3` P, :math:`^1` D, :math:`^1` S) of the Pb 6p block atom, 
using DIRAC and a double-:math:`\zeta` basis set.

The corresponding DIRAC input file (for downloading :download:`Pb.dirac.inp <Pb.dirac.inp>`) is:

.. literalinclude:: Pb.dirac.inp

.. note:: We explicitly give the electronic occupation (.CLOSED SHELL / .OPEN SHELL) which will active the average-of-configuration SCF. DIRAC has an automatic occupation mode but one should always carefully check the output for correctness.

.. note:: The convergence threshold (.EVCCNV) is reduced in order to speed-up the exercises. Note the keyword ".DOSSSS" which activates the (SS|SS) integral classes to be explicitly calculated.

.. note:: ".RESOLVE" activates the open-shell state analysis module and the line "&GOSCIP IPRNT=5 &END" raises the print flag within this module.

Nuclear models in DIRAC
-----------------------

DIRAC offers two different nuclear models, run the p-block atom for the point-nucleus model as well as with a finite nucleus (Gaussian model).

* point nucleus::

   **INTEGRALS
   .NUCMOD
    1

* finite nucleus (Gaussian, default in DIRAC except for ``.NONREL`` and ``.ECP``)::

   **INTEGRALS
   .NUCMOD
    2

**Question**:

1. Do you observe the same trends (total energy, spin-orbit splitting) as with the numerical atomic GRASP calculations on Day 1?

Speed of light and the non-relativistic limit
---------------------------------------------

* vary speed of light *c* between *Z* < *c* < :math:`\infty`::

   **GENERAL
   .CVALUE
    137.0d0


* run a true non-relativistic calculation::

   **HAMILTONIAN
   .NONREL

**Question**:

1. How does the latter calculation compare with a calculation using a large *c* value, *c* :math:`>> 137`?

2. Does the p-shell become degenerate in the non-relativistic limit?

*****************
Exercises - Day 4
*****************


The Hamiltonian families
========================

DIRAC offers a variety of Hamiltonians, ranging from non-relativistic, scalar-relativistic to "fully" relativistic. In
the following we will examine the dependence of the ground state configuration, spin-orbit splitting and lowest excited
state manifold of Pb on the choice of Hamiltonian.

* non-relativistic Hamiltonian(s):

::

   **HAMILTONIAN
   .NONREL

::

   **HAMILTONIAN
   .LEVY-LEBLOND

* scalar-relativistic Hamiltonians (just a selection, look at the DIRAC manual page for even more choices):

::

   **HAMILTONIAN
   .SPINFREE     ! Ken Dyall's spinfree Hamiltonian

::

   **HAMILTONIAN
   .DKH2         ! spinfree Douglas-Kroll-Hess 2nd order
   .SPINFREE

::

   **HAMILTONIAN
   .X2C          ! spinfree exact two-component
   .SPINFREE

::

   **HAMILTONIAN
   .BSS          ! spinfree exact two-component
   109
   .SPINFREE

::

   **HAMILTONIAN
   .BSS          ! spinfree Douglas-Kroll-Hess 2nd order
   102
   .SPINFREE

* relativistic (including spin-orbit coupling) Hamiltonians:

::

   **HAMILTONIAN
           ! Dirac-Coulomb Hamiltonian in which (SS|SS) integrals are neglected
           !and replaced by an interatomic SS correction (calculated as a
           ! classical repulsion term of (tabulated) small component atomic charges)

::

   **HAMILTONIAN
   .LVNEW  ! Dirac-Coulomb Hamiltonian in which (SS|SS) integrals are neglected
           ! and replaced by an interatomic SS correction (with atomic small component
           ! charge calculated via a Mulliken analysis instead of a table look-up

::

   **HAMILTONIAN
   .DOSSSS ! the "full" Dirac-Coulomb Hamiltonian


::

   **HAMILTONIAN
   .DKH2  ! Douglas-Kroll-Hess second order;
          ! implies inclusion of 2-electron spin-same-orbit corrections (AMFI)


::

   **HAMILTONIAN
   .X2C   ! "exact" two-component Hamiltonian;
          ! implies the inclusion of 2-electron spin-same-orbit corrections (AMFI).

::

   **HAMILTONIAN
   .X2C
   .GAUNT ! include 2-electron spin-same and spin-other-orbit corrections.

::

   **HAMILTONIAN
   .X2C
   .NOAMFI ! neglect 2-electron spin-orbit corrections.

::

   **HAMILTONIAN
   .DKH2
   .NOAMFI ! neglect 2-electron spin-orbit corrections.

.. note:: We don't include the zeroth-order approximation (ZORA) to the Dirac-Coulomb Hamiltonian as in its present implementation in DIRAC it only works
 for closed-shell systems. 
 In case you are interested to run ZORA calculations have a look at the :ref:`HAMILTONIAN_.ZORA` keyword in this documentation.

* *special* cases - the molecular-mean-field Hamiltonians:

.. note:: At present DIRAC does not have the functionalities to transform the "Gaunt integrals" to the molecular basis needed for post-Hartree-Fock methods.  The following Hamiltonian combinations allow to approximate the Gaunt interaction in a "molecular-mean-field" fashion.  We can likewise use this approach to approximate the Dirac-Coulomb Hamiltonian ["mean-field" for (SS|LL) + (SS|SS) integrals]. In fact, there are even more molecular-mean-field Hamiltonian available (X2Cmmf) but current implementation restrictions allow those only for the Coupled Cluster module (see Day 4).

::

   **HAMILTONIAN  ! molecular mean-field approximation for the Dirac-Coulomb Hamiltonian: 4DC**
   .DOSSSS
   **MOLTRA
   .INTFL2
   1 1 1 0
   .INTFL4
   1 0 0 0

::

   **HAMILTONIAN ! molecular mean-field approximation for the Dirac-Coulomb-Gaunt Hamiltonian: 4DCG**
   .DOSSSS
   .GAUNT
   **MOLTRA
   .INTFL2
   1 1 1 1
   .INTFL4
   1 0 0 0

**Questions**:

1. How do the two-component (X2C, DKH2) and four-component results compare?

2. Does the inclusion of spin-other-orbit 2-electron interaction (Gaunt contributions) alter the spin-orbit splitting/excitation energies?

3. Is the inclusion of 2-electron spin-orbit corrections (here via "AMFI") important in 2-component spin-orbit calculations (X2C, DKH2)?

4. How do the non-relativistic and various scalar-relativistic ("SPINFREE") results compare?

5. There are two exact-two-component spinfree Hamiltonians listed. Do you obtain slightly different results and if, why?

6. Is the molecular-mean-field approach a reasonable approximation to the "full" 4-component Hamiltonian?

**Further reading and literature**

a. :cite:`Dyall1994`

b. :cite:`Visscher2000`

c. :cite:`Ilias2007`

d. :cite:`Sikkema2009`

*****************
Exercises - Day 5
*****************


The valence ground- and excited state manifold of Se :math:`_2`
===============================================================

After the intense study of atoms with GRASP and DIRAC, we are ready to run our first molecular calculations. In this
exercise we look at the zero-field splitting of the ground state and valence excitation manifold of Se :math:`_2` which is determined by (in approximate spin-orbit
free notation) :math:`(\pi^*)^2` manifold: :math:`^3\Sigma^-_g`,  :math:`^1\Delta_g`,  :math:`^1\Sigma^+_g`.

**Task**:

1. Draw a molecular orbital diagram for the valence part of the Se :math:`_2` molecule. Hint: the valence configuration of Se is p :math:`^4`.

The xyz-molecular input file looks as follows (:download:`Se2.xyz  <Se2.xyz>`):

.. literalinclude:: Se2.xyz

The SCF input reads as (:download:`Se2.dirac.inp  <Se2.dirac.inp>`)

.. literalinclude:: Se2.dirac.inp

Repeat the calculation with the "X2C", "X2C+Gaunt" and "Gaunt molecular-mean-field" Hamiltonian.

**Questions**

1. How does the zero-field splitting change with the inclusion of the Gaunt interaction?

2. Does the X2C Hamiltonian provide accurate results (compared to the fully relativistic calculation)?

A selection of properties in DIRAC
==================================

In our final part of the computer exercises we will look at some property calculations in DIRAC.

NMR shieldings
--------------

The first part concerns the calculation of NMR shieldings at the DFT level using the recently implemented 
simple magnetic balance scheme and the use of London
orbitals in comparison with RKB and URKB calculations.

An illustrative tutorial including a brief introduction to the theory and a guideline for the exercise 
can be found on the DIRAC tutorial web page, :ref:`NMRShieldSimpleMagBal`.

Density at the nucleus and the effective density
------------------------------------------------

The second part illustrates the calculation of the density at the nucleus in comparison with the effective density
(average density over the finite nucleus). We will look at the molecule HI (:download:`HI.xyz <HI.xyz>`)

.. literalinclude:: HI.xyz


The dirac input for the density at the nucleus then reads as (:download:`HI.dc.dirac.inp <HI.dc.dirac.inp>`)

.. literalinclude:: HI.dc.dirac.inp

Let's repeat this exercise with the X2C Hamiltonian, including AMFI corrections (:download:`HI.x2c.dirac.inp <HI.x2c.dirac.inp>`)

.. literalinclude:: HI.x2c.dirac.inp

As the last exercise in this course we will look at picture change errors when transforming from a 4-component to a
2-component formalism. This illustrates the importance of transforming the property operator when you change the
"picture" you are working in! In DIRAC we have the option to **NOT** transform the property operators when doing an X2C
calculation. The keyword is ".NOPCTR" under "**PROPERTIES**".
The input file is (for download, :download:`HI.x2c.nopctr.dirac.inp <HI.x2c.nopctr.dirac.inp>`):

.. literalinclude:: HI.x2c.nopctr.dirac.inp

**Questions**:

1. Compare the effective density and density at the iodine center? What do you observe? Do they agree?

2. Picture-change errors: are they large for both, iodine and hydrogen?

3. Picture-change errors: is there a difference in order of magnitude of the effects for the dipole moment compared with the densities, and if so why?

