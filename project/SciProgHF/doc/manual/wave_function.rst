:orphan:
 

starstar(WAVE FUNCTION)


**Get the wave function**
=========================

This section allows the specification of which wave function module(s)
to activate. By default no modules are activated. To activate any of
these modules you must also specify :ref:`DIRAC_.WAVE FUNCTION`
under :ref:`**DIRAC`,
otherwise this input is not read.

Note that the order below specifies the order in which the different
modules are called if you ask for more than one.


keyword(SCF)

Activates the Hartree-Fock/Kohn-Sham module.

Specification of the SCF module can be given in the :ref:`*SCF` subsection.

If :ref:`HAMILTONIAN_.DFT` has been specified
under :ref:`**HAMILTONIAN`, then a
Kohn-Sham calculation will be performed, otherwise a Hartree-Fock
calculation will be performed.


keyword(RESOLVE)

Resolve open-shell states: do a small CI calculation to get the
individual energies of the states present in an
average-of-configurations open-shell Hartree-Fock calculation (see
:ref:`*RESOLVE`).

keyword(COSCI)

Activates advanced COSCI method, see :ref:`*COSCI`.


keyword(MP2)

Activates :ref:`*MP2CAL`.


keyword(MVO)

Calculate modified virtual orbitals (see :ref:`*MVOCAL`).

keyword(MP2 NO)

Activates the :ref:`*MP2 NO` module to calculate MP2 natural orbitals.

keyword(RELCCSD)

Activates the  :ref:`**RELCC`
(and the :ref:`**MOLTRA`
module to get 4-index
transformed integrals).

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
:ref:`**MOLTRA` section.


keyword(RELADC)

Activates the :ref:`RELADC` and
calculates the single and double ionization spectra by the (A)lgebraic (D)iagrammatic
(C)onstruction ADC. Also activates the :ref:`**MOLTRA`
module to get 4-index
transformed integrals.


keyword(POLPRP)

Activates the :ref:`**POLPRP` module for
calculation of the excitation spectrum by the strict or extended second order
(A)lgebraic (D)iagrammatic (C)onstruction ADC. Also activates the :ref:`**MOLTRA`
module to get 4-index transformed integrals.


keyword(DIRRCI)

Activates the MOLFDIR CI module (and also the
:ref:`**MOLTRA`  module to get 4-index transformed integrals).

Specification of input for the MOLFDIR CI module is given in
the :ref:`DIRRCI`  and :ref:`GOSCIP` sections.

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
:ref:`**MOLTRA`  section.


keyword(LUCITA)

Activates the :ref:`*LUCITA` (and
the :ref:`**MOLTRA` module to get 4-index
transformed integrals).

By default, molecular orbitals with orbital energy between -10 and +20
hartree (a.u.) are included, this can be modified in the
:ref:`**MOLTRA` section.


keyword(EXACC)

Activates the  :ref:`**EXACC` module, the new coupled cluster implementation
based on the ExaTensor library.


**Pre-SCF orbital manipulations**
=================================

keyword(REORDER MO)

Interchange initial molecular orbitals prior to the SCF-calculation. The
start orbitals from DFCOEF are read and reordered.

For each fermion irrep give the new order of orbitals.

*Example:*

::

    .REORDER MO'S
    1..8,10,9


keyword(ORBROT)

Jacobi rotations between pairs of orbitals.

On the line following the keyword, give first the rotation angle, then
on the following line(s) for each fermion irrep, give an :ref:`orbital_strings`
of orbitals to
rotate.


**Post-SCF orbital manipulations**
==================================

keyword(POST SCF REORDER MO)

Interchange converged molecular orbitals. The orbitals from DFCOEF are
read and reordered just before exiting the SCF subroutine.

For each fermion irrep give the new order of orbitals.

*Example:*

::

    .POST DHF REORDER MO'S
    1..8,10,9


keyword(PHCOEF)

Phase adjustment of coefficients DFCOEF: make the largest element of a
given orbital real and positive.


keyword(KRCI)

Activates the :ref:`*KRCI` module
for the calculation of ground and excited states at the relativistic CI level.

keyword(KRMCSCF)

Activates the :ref:`*KRMCSCF`
module for the optimization of ground and excited states (in other than the ground state symmetry) 
at the relativistic MCSCF level.

keyword(LAPLCE)

Activates the :ref:`*LAPLCE`
module to compute weights for Laplace transformation of orbital
energy denominators with the algorithm of Helmich-Paris.
No subsequent calculations, only output of the Laplace points and weights.
