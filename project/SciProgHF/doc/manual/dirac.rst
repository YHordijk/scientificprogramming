:orphan:
 

starstar(DIRAC)

By default DIRAC only activates the generation of the basis set (e.g.  kinetic
balance conditions for small component functions) and the one-electron modules.
Two-electron integrals over the atomic basis functions will be calculated when
needed by the job modules given below.

We recommend that you use the `Dalton <http://www.daltonprogram.org>`_ program
package if you want to e.g. save the two-electron integrals to disk for another
purpose, this is not possible with DIRAC.


keyword(WAVE FUNCTION)

Activates the wave function module(s).  This activates the reading of the
:ref:`**WAVE FUNCTION` section, where the desired wave function type(s) must be
specified.  If you read in converged MO coefficients from the file DFCOEF and
you want to skip the SCF step and proceed directly to properties or the
post-SCF step, then a trick to do this is to comment out (or remove) this
keyword.


keyword(ANALYZE)

Activates the Hartree--Fock or Kohn--Sham analysis module.  This activates the
reading of the :ref:`**ANALYZE` section.


keyword(PROPERTIES)

Activates the property module (which will call the integral module for property
integrals).  This activates the reading of the :ref:`**PROPERTIES` section.


keyword(OPTIMIZE)

Activates the geometry optimization.  This activates the reading of the
:ref:`*OPTIMIZE` subsection.


keyword(4INDEX)

Explicitly activates the transformation of integrals to molecular orbital
basis.  This activates the reading of the :ref:`**MOLTRA` section.

These transformed integrals are currently only used by the :ref:`**RELCC`,
:ref:`RELADC`, :ref:`**POLPRP`, :ref:`DIRRCI`, and :ref:`*LUCITA` modules, and if one of these
three modules are requested under :ref:`**WAVE FUNCTION`, then this flag is
automatically activated unless .NO4INDEX is specified in this input module.

By default, molecular orbitals with orbital energy between -10 and +20 hartree
are included, this can be modified in the :ref:`**MOLTRA` section.


keyword(NO4INDEX)

Do not automatically activate integral transformation to molecular
orbital basis if any of
:ref:`**RELCC`,
:ref:`RELADC`,
:ref:`**POLPRP`,
:ref:`DIRRCI`, and
:ref:`*LUCITA`
modules
are requested under :ref:`**WAVE FUNCTION`.

This keyword is utilized when repeating correlated CC or CI calculations
(with different parameters for instance) based on saved files after the
integral transformation.


keyword(TITLE)

Title line (max. 50 characters).
Example::

  .TITLE
   my first DIRAC calculation


keyword(INPTEST)

Input test - no job modules are called, only verification of DIRAC input files.
It is often useful to start a new set of calculations with an input test in
order to check that input file processing is correct before submitting your
(long-term run) job.


keyword(ONLY INTEGRALS)

Stop after the calculation of the one-electron integrals for the Hamiltonian
and the one-electron integrals specified under :ref:`**INTEGRALS`.  The integrals
are written to disk.


keyword(NOSFCR)

.. warning:: documentation missing


keyword(MINIMI)

.. warning:: documentation missing


keyword(XMLOUT)

Create the output xml-file containing selected data (in development).
