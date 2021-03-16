:orphan:
 

starstar(ANALYZE)

Wave function analyzing tool.

This section is aimed to analyze the final Hartree-Fock wave function by
using one or more analysis modules. By default none of them is
activated.


keyword(PRIVEC)

Print vectors. Activates the :ref:`*PRIVEC` subsection.


keyword(MULPOP)

Perform Mulliken population analysis :cite:`Mulliken1955`.
Continues to :ref:`*MULPOP` subsection.


keyword(PROJECTION)

Perform projection analysis :cite:`Faegri2001` Activates the :ref:`*PROJECTION`
subsection.


keyword(DENSITY)

Write density to a formatted file in Gaussian cube format. Activates the
:ref:`*DENSITY` subsection.


keyword(LOCALIZATION)

Localize orbitals using the Pipek-Mezey criterion :cite:`Dubillard2006` .

Continues to the activated :ref:`*LOCALIZATION`
subsection. Note that in the present implementation this only works in 
:math:`C_1` symmetry.

