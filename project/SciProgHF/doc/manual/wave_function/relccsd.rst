:orphan:
 

.. _relccsd:

============================================
Kramers unrestricted Coupled-Cluster methods
============================================

See :cite:`Visscher1996` for the initial CCSD(T) implementation. The Fock-space CC 
implementation is described in :cite:`Visscher2001`.

starstar(RELCC)

Specification of reference determinant, type of calculation, and general
settings.

keyword(FOCKSP)

Activate the Fock space module. This option should be used for
multireference calculations. See further :ref:`*CCFSPC`.
Because the first sector of Fock space gives the same result as a 
regular CCSD calculation, the latter calculation is switched off. If you
do wish to perform also a regular CC calculation (e.g. to get the CCSD(T) energy)
you need to activate this explicitly via :ref:`RELCC_.ENERGY` (see below).

keyword(ENERGY)

Activate the energy calculation. This is the preferred option for
calculations on closed shell or simple open shell systems and need not
be specified explictly in such cases. For Fock space calculations the
keyword switches on a separate single reference calculation done prior
to the FS calculation.

*Default:*

::

    Perform energy calculation.

keyword(GRADIENT)

Calculate the effective 1-particle density matrix. This option can be used to calculate
molecular properties. For closed-shell MP2 see :cite:`vanStralen2005`. For closed-shell CC 
see :cite:`Shee2016`.

This calculation generates the file CCDENS, which contains the MP2 or CC ground-state density 
matrix in AO basis. In this release, CCDENS is used by the property module to calculate 
ground-state expectation values. 

If saved, CCDENS can be used in a property calculation (see :ref:`PROPERTIES_.RDCCDENS`) 
without the need to invoke this module.


*Default:*

::

    No gradient calculation.

keyword(EOMCC)

Activate the equation of motion (EOMCC) module. This option can be used to obtain
excitation energies (EOM-EE), singly ionized (EOM-IP) and electron attached (EOM-EA)
states starting from a closed shell reference and at the CCSD level, by solving for
the right-hand eigenvalues and eigenvectors of the appropriate similarity-transformed 
Hamiltonian.  

Unlike in Fock-space calculations, where the full model Hamiltonian for each sector
is diagonalized, in EOMCC the user must specify the desired number of roots per symmetry,
(see :ref:`*EOMCC` for further details) and these are determined via matrix-free 
diagonalization procedures such as the Davidson method (see :ref:`*CCDIAG` for the 
options available). 

*Default:*

::

    No EOM-CC calculation.

keyword(EOMPROP)

.. warning:: Development version only. 

Activate the calculation of excited-state first-order properties for the EOMCC module for
CCSD wavefunctions, by solving for the left-hand eigenvectors of the corresponding 
similarity-transformed Hamiltonian. It should be specified along with  :ref:`RELCC_.EOMCC` 
as it requires the right-hand eigenvectors and, as such, the same type of calculation 
(EOM-EE/IP/EA) will be performed. See :ref:`*EOMPROP` for further details.

*Default:*

::

    No EOM-CC excited-state expectation value calculation.

keyword(NELEC)

Number of correlated electrons. This variable determines the reference determinant to be used in the exponential expansion of the wave function. 
Since the default values correspond to the information passed on by the MOLTRA code on basis of the Hartree-Fock occupations 
and chosen range of active orbitals in MOLTRA, **for CLOSED SHELL systems there is usually no need to specify this variable manually**. 
If you do chose to specify this manually, you should make sure to count electrons carefully as the numbers relate to the number of correlated
electrons rather than the total number. This number may therefore change if you change tresholds in the integral transformation.

For OPEN SHELL systems it is usually better to employ the multireference Fock space approach, except for simple cases such as 
a high-spin open shell state. In this case the single reference CCSD(T) ansatz works rater well. For such cases it is easier to use 
the :ref:`RELCC_.NELEC_OPEN` keyword that allows to just specify the distribution of open shell electrons over the irreps, but for backwards
compatibility we retain the older NELEC option as well.

*Arguments:*

::

    Integer (NELEC(I),I=1,NFSYM*2).

*Default:*

::

     Number of correlated electrons in closed shells in these irreps (written by **MOLTRA).

keyword(NELEC_OPEN)

Distribution of correlated open shell electrons over the irreps.  This determines the reference determinant for open shell single reference calculations.

This input should always be given for average-of-configuration SCF calculations that are followed by a CC calculation. The total number of electrons that is
given should correspond to the number of open shell electrons.

*Arguments:*

::

    Integer (NELEC_OPEN(I),I=1,NFSYM*2).

*Default:*

::

     Zero. (note that this default leads to wrong results because a closed shell ion will be calculated if no input is given).

keyword(NEL_F1)

Number of electrons in the *gerade* irreps of the Abelian symmetry group. This works like the NELEC keyword, but uses the 
supersymmetry possible for linear groups in which irreps are ordered as 1/2, -1/2, 3/2, -3/2, 5/2, ... (with the number
being the m_j value). 

keyword(NEL_F2)

Number of electrons in the *ungerade* irreps of the Abelian symmetry group. 

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(COUNTMEM)

Stop CC module after counting the total memory demand. Needs only MRCONEE. 

keyword(TIMING)

Print detailed timing information.

*Default:*

::

    Only limited timing information is printed.

keyword(DEBUG)

Print debug information.

*Default:*

::

    Debug information is not printed.

keyword(RESTART)

Reuses information from prior calculations to resume the calculation from the last sucessfully recorded checkpoint. In order for a restart to be possible, the appropriate files from both the 4-index transformation - e.g. MRCONEE, MDCINT (and MDCINX* for parallel runs) and, if the transformation of property operators was requested, MDPROP - and from the previous coupled cluster run (ft.* for the sorting of integrals, MCCRES* for the amplitudes) have to be present at the scratch directory. Given the size and number of these files, this means in practice one must request the the scratch directory is kept at the end of the runs (see the pam script options for details)

*Default:*

::

    No restart is performed.

keyword(NOSORT)

Forces the code to skip the sorting of integrals coming from the 4-index transformation into six integral classes. Using this option without the integral sorting stap having been properly carried out may produce incorrect results. 

*Default:*

::

    Sorting of integrals into classes is performed.


star(CCENER)

Covers options related to energy.

keyword(NOMP2)

Deactivate MP2 calculation.

keyword(NOSD)

Deactivate CCSD calculation.

keyword(NOSDT)

Deactivate the calculation of perturbative triples. This is potentially
useful when running into memory problems for very big calculations and will also
save some CPU time.

keyword(MAXIT)

Set maximum number of iterations allowed to solve the CC equations.

keyword(MAXDIM)

Set maximum number of amplitude vectors used in the DIIS extrapolation.

keyword(NTOL)

Specify requested convergence (10^-NTOL) in the amplitudes.

keyword(NOSING)

Eliminate T1 amplitudes in the calculation (only interesting for test
purposes, this gives no computational speed-up).

keyword(NODOUB)

Eliminate T2 amplitudes in the calculation (only interesting for test
purposes, this gives no computational speed-up).
Deactivate contribution from doubles; corresponds to a CCS calculation.

star(CCFOPR)

Calculate first-order properties (expectation values) for the MP2 and CCSD wave
function.

keyword(MP2G)
Use MP2 wave function. 

keyword(CCSDG)
Use CCSD wave function.

keyword(NATORB)

Calculate natural orbitals (currently only for MP2 density matrix)

keyword(RELAXED)

Use orbital-relaxed density matrix (currently only for MP2).
The current default for MP2 and CC wave functions is to use the unrelaxed 
density matrix. This is computationally less expensive, but also less accurate.



star(EOMCC)

This menu controls the parameters for the definition of the EOM model used and the number of roots per symmetry. 

Currently the implementation supports the excitation energy (EE), single electron attachment (EA) and detachment (IP) models for CCSD wavefunctions only. 

Note that there is no default, if no options are selected no calculations will be performed. 

keyword(EE)

Selects excitation energy calculations. 

Expects two integers, the first specifying the state symmetry number and the second the number of states for this symmetry.
The keyword should be repeated multiple times if different symmetries are desired. The  implementation does not allow for mixing
the different EOM models in the same run.

NB: finding the state symmetry number for EE may require some experimentation, the totally symmetry irrep is always number 1,
but the order of the other boson irreps depends on the double group that is used and how the molecule is oriented.

*Example:* requesting two excitation energies for symmetry 1, four for symmetry 3, and one for symmetry 8::

      .EE
      1 2
      .EE
      3 4
      .EE
      8 1

*Default:*

::

     No excitation energies requested


keyword(EA)

Selects single electron attachment calculations. 

Expects two integers, the first specifying the state symmetry number and the second the number of states for this symmetry.
The keyword should be repeated multiple times if different symmetries are desired. The  implementation does not allow for mixing
the different EOM models in the same run.

NB: the state symmetry number for EA depends on the ordering of irreps for the spinors, this information is printed at
the beginning of the RELCCSD output.


*Example:* requesting two electron attachment energies for symmetry 1, four for symmetry 3, and one for symmetry 8::

      .EA
      1 2
      .EA
      3 4
      .EA
      8 1

*Default:*

::

     No excitation energies requested

keyword(IP)

Selects electron detachment calculations. 

Expects two integers, the first specifying the state symmetry number and the second the number of states for this symmetry.
The keyword should be repeated multiple times if different symmetries are desired. The  implementation does not allow for mixing
the different EOM models in the same run.

NB: the state symmetry number for IP depends on the ordering of irreps for the spinors, this information is printed at
the beginning of the RELCCSD output.


*Example:* requesting two single electron detachment energies for symmetry 1, four for symmetry 3, and one for symmetry 8::

      .IP
      1 2
      .IP
      3 4
      .IP
      8 1

*Default:*

::

     No excitation energies requested

star(EOMPROP)

keyword(EXCPRP)

.. warning:: Development version only. 

number of excited states to calculate expectation values for, on a given symmetry.

Arguments: two integers, the first specifying the state symmetry and the second the number of states for this symmetry

star(CCDIAG)

This menu controls the parameters to the iterative Davidson diagonalization procedure.

Note that the implementation currently does not support the saving of eigenvectors, so no restarts are possible.

keyword(CONVERG)

Sets the convergence threshold on the norm of the residual vectors

*Default:*

::

     .CONVERG
     1.0E-8

keyword(MAXSIZE)

Sets the maximum size of the subspace

*Default:*

::

     .MAXSIZE
     128

keyword(MAXITER)

Sets the maximum number of iterations

*Default:*

::

     .MAXITR
     80

keyword(TRV_I)
`
Creates trial vectors based on unit vectors sorted by lowest energy of the diagonal of the similarity-transoformed Hamitonian (pivoting).


*Default:*

::

    Enabled

keyword(TRV_FULLMATRIX)

Creates unit trial vectors from a unit matrix. This roughly means that the first trial vectors will be those for high-lying states. 

If combined with a number of EE/EA/IP roots set to -1 and a maximum number of iterations set to 1, the full similarity transformed Hamiltonian for a given symmetry will be diagonalized. 
This yields the full spectrum, but for anything other than CC singles will require significant amounts of memory, so in practice such a combination is useful for debug purposes only.

*Default:*

::

    Disabled 

keyword(TRV_CCS)

Creates trial vectors from the eigenvectors of the singles-singles EOM-EE/IP/EA blocks.

keyword(TRV_R2L) 
keyword(TRV_NOR2L)

Control the use the right-hand side (R) vectors as guess for the left-hand side (L). If disabled, (.TRV_NOR2L) trial vector generation will be controlled by the options above.

*Default:*

::

    .TRV_R2L

keyword(OVERLAP) 
keyword(NOOVERLAP)

Controls the use of root following via overlap in sorting the left/right-hand side eigenvalues and eigenvectors.

*Default:*

::

    .NOOVERLAP


keyword(DEBUG)

Enables debug printout. Massive output.

*Default:*

::

    Does not print out debug information



star(CCPROJ)

This menu controls the use of projection operators in the CC calculations.

These projectors are used to remove terms from the ground state CC equations (frozen core), the EOM-CC similarity transformed Hamiltonian (CVS and restricted excitation window), or a combination of both.

keyword(FCORE)

Defines which active occupied spinors to consider as part of frozen core within the ground-state CC calculation.

Input is an threshold energy value in atomic units. T amplitudes contaning spinors with energies below such value are set to zero.

Note this value does not have to be the same as for the definition of the core spinors for core-valence separation. Assuring the consistency of these choices is up to the user.

 
*Example:* requests that all spinors with energies below -10.0 a.u. are frozen::

      .FCORE
      -10.0

*Default:*

::

     No frozen core spinors


keyword(CVS_CORE)

Defines which active occupied spinors to consider as core spinors for EOM-CC calculations employing the core-valence separation (CVS) framework. To be used with EOM-IP or EOM-EE, in order to target high-lying excited/ionized states, namely for core spectra.

Input is an threshold energy value in atomic units. 

Right/Left-hand side trial vector and similarity transformed Hamiltonian matrix elements representing excited determinants containing at least one occupied spinor with energies below such value are set to zero.

Note this value does not have to be the same as for the definition of the frozen core. Assuring the consistency of these choices is up to the user. 

In addition to this, users should make sure to enable root following (see :ref:`CCDIAG_.OVERLAP`) if they want to target states with dominant singly ionized or singly excited character.


*Example:* requests that all spinors with energies below -10.0 a.u. are considered to be core::

      .CVS_CORE
      -10.0

*Default:*

::

     No core-valence separation is enforced


keyword(NODOCC)

Requests that, when using core-valence separation,  Right/Left-hand side trial vector and similarity transformed Hamiltonian matrix elements representing excited determinants containing both core occupied spinors to be projected out. 


*Example:* requests that double core excited determinants are projected out::

      .NODOCC

*Default:*

::

     Double core excited determinants are retained in the calculations


keyword(REW_OCC)
keyword(REW_VIRT)

Defines a set or active occupied and/or virtual spinors as part of restricted window for EOM-CC calculations. Can be used with any EOM variant. 

If an occupied/virtual restricted window is defined, the Right/Left-hand side trial vector and similarity transformed Hamiltonian matrix elements representing excited determinants containing at least one occupied/virtual  spinor within the window(s) are retained in the calculation, the rest being projected out.

Input for each keyword is a pair of values defining thethreshold energy value in atomic units, one per line. Note that the two keywords can be used independently.

Right/Left-hand side trial vector and similarity transformed Hamiltonian matrix elements representing excited determinants containing at least one occupied spinor with energies below such value are set to zero.

The option :ref:`CCPROJ_.REW_OCC` is incompatible with :ref:`CCPROJ_.CVS_CORE`. It can however be used with :ref:`CCPROJ_.FCORE`, but assuring the consistency of these choices is up to the user.

In addition to this, users should make sure to enable root following (see under :ref:`*EOMCC`) if they want to target states with dominant singly ionized or singly excited character.


*Example:* requests two restricted windows, one for occupied spinors with energies between -200.0 a.u. and -150.0 a.u. and another for virtual spinors with energies between 0.001 a.u. and 35.0 a.u.::

      .REW_OCC      
       -200.0   
       -150.0
      .REW_VIRT      
       0.001 
       35.0

*Default:*

::

     No restricted window is defined, all occupied and/or virtual spinors are active










star(CCFSPC)

Perform a Fock space MRCC calculation. Fock space allows variable
particle number. Sectors :math:`(m,n)` in Fock space corresponds to
:math:`N+n-m` - electron states obtained by the generation of
:math:`m` holes (electron removal) and :math:`n` particles
(electron attachment)  with respect to
a closed-shell reference determinant, the :math:`(0,0)` sector. 
Within each specified sector (and the lower ones), an effective Hamiltonian
is built and diagonalized to give CC energies for a set of
states.

keyword(DOIH)

Use the Intermediate Hamiltonian formalism in which an auxiliary space
is used to prevent the "intruder state" problem. Default: IH formalism
not used.

keyword(DOEA)

Calculate electron affinities (add one electron to the reference state,
allowing occupation of the active virtual orbitals, corresponding to the
:math:`(0,1)` sector).

keyword(DOIE)

Calculate ionization energies (remove one electron from the reference
state, allowing depletion of the active occupied orbitals, corresponding to the
:math:`(1,0)` sector).

keyword(DOEA2)

Calculate second electron affinities (add two electrons to the reference
state, allowing occupation of the active virtual orbitals, corresponding to the
:math:`(0,2)` sector).

keyword(DOIE2)

Calculate second ionization energies (remove two electrons from the
reference state, allowing depletion of the active occupied orbitals, corresponding to the :math:`(2,0)` sector).

keyword(DOEXC)

Calculate excitation energies (allow excitation from the set of active
occupied orbitals to the set of active virtual orbitals, corresponding to the
:math:`(1,1)` sector).

keyword(NACTH)

Specification of the set of active hole orbitals (from which
ionization/excitation takes place)

keyword(NACTP)

Specification of the set of active particle orbitals (to which electron
attachment/excitation takes place)

keyword(MAXIT)

Maximum number of iterations allowed to solve the FSCC equations

keyword(MAXDIM)

Set maximum number of amplitude vectors used in the DIIS extrapolation.

keyword(NTOL)

Specify requested convergence (10^-NTOL) in the amplitudes.

keyword(GESTAT)

Specify the state number in the last active sector to pick the energy 
from (remember to account for degeneracies) for a state-specific FSCC 
geometry optimization based on a numerical gradient. 

star(CCIH)

Options for intermediate hamiltonian in FSCC.

keyword(EHMIN)

Minimum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
secundary (Q) space and do not contribute to the model space.

low limit of orbital energies of active occupied
orbitals, which constitute the secondary Pi space. Could be
used in (1,0), (2,0) and (1,1) sectors. Arguments: real.

keyword(EHMAX)

Maximum orbital energy of occupied orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
primary (Pm) space.

This is upper limit of one-electronic energies of
active occupied orbitals, which constitute the secondary Pi space.
Could be used in (1,0), (2,0) and (1,1) sectors. Arguments: real.

keyword(EPMIN)

Minimum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies lower than this energy are taken in the
primary (Pm) space.

This is the low limit of orbital energies of active virtual
orbitals, which constitute the secondary Pi space.
Could be used in (0,1), (0,2) and (1,1) sectors. Arguments: real.

keyword(EPMAX)

Maximum orbital energy of virtual orbitals forming the auxiliary (Pi)
space. Orbitals with energies higher than this energy are taken in the
secundary (Q) space and do not contribute to the model space.

This is the upper limit of one-electronic energies of active virtual
orbitals, which constitute the secondary Pi space. Could be
used in (0,1), (0,2) and (1,1) sectors. Arguments: real.
 
Other Intermediate Hamiltonian (IH) input parameters
====================================================
 
For *experts* only.

Following keywords belong to the CCIH namelist section.
 
keyword(IHSCHEME)

Choose particular IH scheme. Arguments: Integer IHSCHEME =  1,  or 2.
 
The IHSCHEME=1 corresponds to the extrapolated IH
(XIH) approach, described in the paper :cite:`Eliav2005`.

Main idea: proper modification of the energetic denominators,
containing problematic Pi  -> Q transition.
The original denominator 1/(E_Pi  -  E_Q) , used during CC
iterations, is substituted by the following expression (1)

.. math::

 \frac{(1-[\frac{AIH*SHIFT}{(E_{Pi} - E_{Q} + SHIFT)}]^{NIH})}{\frac{(1-AIH*SHIFT}{(E_{Pi}  -  E_{Q} + SHIFT))}},

here AIH, SHIFT,NIH are parameters, specially chosen for
overcoming of the intruder states problem. These parameters
could be used in the procedure of the extrapolation of
the "exact" effective Hamiltonian solutions from
corresponding IH CC energies and wave functions.
 
The IHSCHEME=2 corresponds to the simplified IH-2 approach, described in
the paper :cite:`Landau2004`.
 
Here the problematic denominators :math:`1/(E_{Pi}  -  E_{Q})` are
substituted simply by the factor 0.
 
Default: IHSCHEME = 2
 
Next key options are used only in case of XIH (IHSCHEME = 1).
 
keyword(SHIFTH11)

Energy shift for the one-electronic excitations in (1,0) sector. Arguments: real.

keyword(SHIFTH12)

Energy shift for the two-electronic excitations in (1,0) sector. Arguments: real.

keyword(SHIFTH2)

Energy shift for the two-electronic excitations in (2,0) sector. Arguments: real.

keyword(SHIFTP11)

Energy shift for the one-electronic excitations in (0,1) sector. Arguments: real.

keyword(SHIFTP12)

Energy shift for the two-electronic excitations in (0,1) sector. Arguments: real.

keyword(SHIFTP2)

Energy shift for the two-electronic excitations in (0,2) sector. Arguments: real.
Usually we choose the approximate difference between the highest orbital energy belonging to Pi and the lowest
orbital energy belonging to the Pm space. 
Works only with the old style of RELCC input.

keyword(SHIFTHP)

Energy shift for the two-electronic excitations in (1,1) sector. Arguments: real

keyword(AIH)

Compensation factor, used in expression (1). Arguments: real positive, not greater then 1.0.
 
keyword(NIH)

Compensation power, used in expression (1). Arguments: integer.
 
In the case of the limit: AIH=1.0 and NIH -> "infinity" ( NIH>100, in practice)
we have so called "full compensation" method, corresponding to the
extrapolation of the effective Hamiltonian from the intermediate one.

star(CCSORT)

Specialist options related to the sorting of two-electron integrals
and the calculation of the reference Fock matrix.

keyword(NORECMP)

Do not recompute the Fock matrix, but assume a diagonal matrix with the orbital energies 
taken from the SCF program on the dioagonal. This is usually not recommended as the latter
correspond to a restricted open shell expression and RELCCSD uses an unrestricted formalism.
For closed shell systems the two expressions are identical and this option merely suppresses
a build-in check on the accuracy of transformed integrals.

keyword(USEOE)

Ignore recomputed Fock matrix and use orbital energies supplied by the SCF program. This option is
sometimes useful for degenerate open shell cases in which case the perturbation theory for the unrestricted formalism
is not invariant for rotations among degenerate orbitals. It should only change the outcome
of the [T], (T) and -T energy corrections.

star(CCRESTART)

Control parameters for the restart option. The default behavior of the restart option is to verify whether in the RELCCSD the successive 
checkpoints have been passed, and restart the calculation at the first one which is not flagged as "Completed, restartable". 

For example, one would have 

::

 Status of the calculations
 Integral sort # 1 :                   Completed, restartable
 Integral sort # 2 :                   Completed, restartable
 Fock matrix build :                   Completed, restartable
 MP2 energy calculation :              Completed, restartable
 CCSD energy calculation :             Completed, restartable
 CCSD(T) energy calculation :          Completed, restartable

if a calculations has successfully has passed through all checkpoints.

For the single-reference calculations the restart is straightforward and in general no additional keywords are necessary. Users must nevertheless be careful that the restart is performed for exactly the same calculation, as there are no interla consistency checks in place in the case of a restart: this means that, for example, if either the geometry or the number of corrlated electrons or virtual spinors has been changed after the initial calculations, if the RESTART option is on the code will proceed with the old checkpoint data (sorted integrals, etc) and any results will be erroneous. 

Fock-space calculations can also be restarted, but for these additional care must be taken:
(1) One must redo the integral sorting steps if the model (P) or correlation (Q) spaces change dimension. 
A change in the values that define the partition between main (Pm) and intermediate (Pi) model spaces, on the other hand, 
does not require the integral sorting to be redone.
(2) One has to specifically indentify which which sectors of Fock space are to be (re)calculated and which should be skipped. Some examples of such a procedure can be found in the test set. 

Important: for EOMCC calculations it is currently only possible to restart from ground-state calculations. In this case, the ".REDOCSD" keyword *must* be used, otherwise the code will produce incorrect results.


keyword(UNCONVERGED)

Forces results from unconverged iterative procedures to be considered as converged. 

*Default:*

::

    False

keyword(FORCER)

Forces restart even if setup is potentially different (number of electrons, active/inactive spinors  etc).

*Default:*

::
      
      False
    
keyword(REDOCCSD)

Forces the iterative procedure to solve the CCSD (or CCD or CCS) equations to be performed, even if in a prior run it has been marked as completed.

*Default:*

::

      False 
     

keyword(REDOSECT)

In the case of fock-space calculations, indicate which sectors are to be relcalculated during the restart.

This keyword expects two lines as input; on the first line, an integer specifying how many sectors are recalculated
and in the second line a list of the sectors in question in the RELCC notation. 

::

    .REDOSECT
     2
     00 01

will redo sectors 0h0p and 0h1p 

*Default:*

::

    All sectors are recalculated 

keyword(SKIPSECT)

In the case of fock-space calculations, indicate which sectors are to be skipped during the restart

This keyword expects two lines as input; on the first line, an integer specifying how many sectors are recalculated
and in the second line a list of the sectors in question in the RELCC notation.

::

    .SKIPSECT
     2
     00 01

will bypass the calculation of sectors 0h0p and 0h1p.


*Default:*

::

    No sectors are skipped


keyword(REDOSORT)

Forces the sorting into the six integral classes to be performed again, even if prior sorting was completed. 

*Default:*

::

    False. integrals are not resorted if the prior sorting was completed.


