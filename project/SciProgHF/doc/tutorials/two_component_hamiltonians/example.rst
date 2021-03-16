:orphan:

Spin-orbit states from the COSCI method
=======================================

This tutorial demonstrates the importance of the effective mean-field spin-orbit screening
on spin-orbit states of open-shell systems. Several two-component Hamiltonians are employed.

Spin-orbit states of the F atom
-------------------------------

In the DIRAC test we calculate the energy difference between spin-orbit splitted states of the :math:`^{2}P` state of Fluorine,
using the COSCI wavefunction and with several different Hamiltonians. 
All input files for download (together with output files) are in the corresponding test directory of DIRAC, :file:`test/cosci_energy`.

The following table shows the energy difference betweem :math:`X ^{2}P_{3/2}` and :math:`A ^{2}P_{1/2}` states:

+-------------------+-------------------+
| Hamiltonian       |   Splitting/cm-1  |
+===================+===================+
| DC                |     434.511758    |
+-------------------+-------------------+
| BSS+MFSSO         |     438.792872    |
+-------------------+-------------------+
| BSS_RKB+MFSSO(*)  |     438.793184    |
+-------------------+-------------------+
| DKH2+MFSSO        |     438.792782    |
+-------------------+-------------------+
| BSSsfBSO1+MFSSO   |     438.868634    |
+-------------------+-------------------+
| DKH2sfBSO1+MFSSO  |     438.868738    |
+-------------------+-------------------+
| BSSsfESO1+MFSSO   |     438.866098    |
+-------------------+-------------------+
| DKH2sfESO1+MFSSO  |     438.866201    |
+-------------------+-------------------+
| BSS               |     583.459766    |
+-------------------+-------------------+
| BSS_RKB(**)       |     583.459995    |
+-------------------+-------------------+
| DKH2              |     583.459700    |
+-------------------+-------------------+
| BSSsfESO1         |     583.533060    |
+-------------------+-------------------+
| DKH2sfESO1        |     583.533187    |
+-------------------+-------------------+
| BSSsfBSO1         |     583.535908    |
+-------------------+-------------------+
| DKH2sfBSO1        |     583.536036    |
+-------------------+-------------------+
| DC2BSS_RKB(DF)    |     585.906861    |
+-------------------+-------------------+

 (*) Known as X2C.
 (**) Known as X2C-NOAMFI.

Calculated values can be devided into two categories: those with the mean-field spin-orbit term (MFSSO)
and those without.
Results matching the four-component Dirac-Coulomb (DC) Hamiltonian are those
containing the MFSSO screening term.

For more information, see Refs. [Ilias2001]_, [Ilias2007]_.

Spin-orbit states of the :math:`Rn^{77+}` cation
-------------------------------------------------

Let us proceed with the isoelectronic, but heavier system:  the Fluorine-like (9 electrons), highly charged :math:`Rn^{77+}` cation (Z=86).
All input files for download (together with output files) are in the corresponding test directory of DIRAC, :file:`test/cosci_energy`.
Calculated energy differences between the ground, :math:`X ^{2}P_{3/2}`, and the first excited state, :math:`A ^{2}P_{1/2}`,
are in the following table:

+-------------------+-------------------+
| Hamiltonian       |   Splitting/eV    |
+===================+===================+
| DC                |   3700.081        |
+-------------------+-------------------+
| BSS+MFSSO         |   3796.844        |
+-------------------+-------------------+
| DKH2+MFSSO        |   3777.837        |
+-------------------+-------------------+
| DC2BSS_RKB(DF)    |   3810.190        |
+-------------------+-------------------+
| BSS               |   3808.859        |
+-------------------+-------------------+
| BSS_RKB (*)       |   3810.273        |
+-------------------+-------------------+
| DKH2              |   3790.044        |
+-------------------+-------------------+
| DKH2sfBSO1+MFSSO  |   4047.324        |
+-------------------+-------------------+
| DKH2sfBSO1        |   4056.349        |
+-------------------+-------------------+

 (*) Known as X2C-NOAMFI.


Excercises
----------

1. Why is the MFSSO term more important for the ligher element (F) than for the  heavy :math:`Rn^{77+}` ? 

2. The one-electron spin-orbit term, SO1, is sufficient for representing spin-orbital effects in the Flourine atom,
   but not of the `Rn^{77+}` cation. Why ?

3. For the Flourine atom, increase the speed of light (:ref:`GENERAL_.CVALUE`) in four-component calculations to emulate non-relativistic description.
   What is the effect on the spin-orbit splitting ? 
   What artificial value of the speed of light generates the DC-SCF energy identical with nonrelativistic SCF energy up to 5 decimal places ?

4. To "increase" relativistic effects in Flourine, decrease the speed of light in four-component calculations.
   How does it affect the spin-orbit splitting ?

5. Change the symmetry from D2h to automatic symmetry detection in the F mol file and
   add molecular spinors analysis to the input file (:ref:`**ANALYZE`).
   Identify molecular spinors (orbitals) of Flourine according to the extra quantum number in linear symmetry.

