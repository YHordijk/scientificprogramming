:orphan:
 

starstar(MOLTRA)

Integral transformation module.

This module transforms integrals from the scalar AO basis to integrals in the
molecular spinor basis. The algorithm is not yet described in detail but it may
be helpful to read Ref. :cite:`Visscher2002` to find information about the double
quaternion formalism that is used. This module is usually invoked by one of the
correlation modules but can also be run as a stand-alone module by specifying
:ref:`DIRAC_.4INDEX` in the :ref:`**DIRAC` section of input. The input serves to
change the default active space and to specify special wishes (e.g modification
of :ref:`MOLTRA_.INTFLG`).

Three separate transformations are considered:

-  Transformation of a core Fock matrix, giving effective one-electron
   matrices.
-  Transformation of two-electron Coulomb integrals.
-  Transformation of integrals over general (one-electron) operators.

The output is written to the files MDCINT, MRCONEE, and MDPROP, respectively,
and can be read in by the :ref:`**RELCC` and :ref:`DIRRCI` and programs to
perform correlated calculations and by the :ref:`RELADC` and :ref:`**POLPRP` modules
to perform propagator calculations. Some common input needs to be specified both in this
section and in the :ref:`DIRRCI` sections.

keyword(ACTIVE)

Specify the active set of spinors.

*Arguments:* One line with an :ref:`orbital_strings` for each fermion irrep.

*Default:*

::

    energy -10.0 20.0 1.0

**Advanced options**

keyword(THROUT)

Threshold value for writing transformed 2-electron integrals to file. A
higher value will decrease the size of the MDCINT files but will
decrease the accuracy of the calculation.

*Default:*

::

    .THROUT
     1.0D-14

keyword(PRPTRA)

Perform transformation of property integrals. The type(s) of integrals
that are to be transformed can be specified with the .OPERATOR keyword.

*Default:* Do not transform property integrals.

keyword(SCHEME)

Selects the integral transformation algorithm ("scheme"). Several
schemes have been implemented, but users are to consider only "scheme 4"
and "scheme 6" for production calculations. Other schemes refer to
outdated to special-purpose algorithms, and should be used with care.
The most suitable scheme to use depends to a significant extent on the
type of calculation in question (size of the active space etc), and on
the computational resources used (compilers, amount of memory/scratch
disk space, whether or not disks are shared to local to compute nodes
etc).

*Scheme 6 (Default):*

::

    .SCHEME
     6

Scheme 6 is designed to decrease I/O operations at the expense of more
communication between processes, something which is more often than not
advantageous in today's architectures. It should be noted, however, that
the significant amount of communication between processes happens for
the so-called "1HT" step ("1HT" standing for first half-transformation,
since it is where the first two indexes of the electron repulsion
integrals are transformed to MO basis and these intermediate,
half-transformed entities are stored on disk), with no communication on
the "2HT" step (which will yield the final integrals in MO basis). Due
to the large amount of half-transformed entities, this scheme may demand
large scratch spaces (e.g. a Dirac-Coulomb calculation on water with the
aug-cc-pV5Z basis, where the full spinor space is active and LL and SL
integrals are transformed can reach nearly 1Tb of disk use) but with the
advantage that the full size of the problem is exactly divided between
processes - so the larger the number of MPI processes, the less disk
space is used per process.

 *Scheme 4:*

::

    .SCHEME
     4

Scheme 4 is designed to reduce the communication between different
processors at the expense of increasing I/O operations. This increase in
I/O operations stems from the integrals being replicated for all MPI
processes, so that doubling the number of MPI processes doubles the
total disk usage but not the usage per MPI process.

keyword(HTSORT)

In connection to the scheme 6 the
user can also choose whether or perform an intermediate sorting step
prior to the second half-transformation step, via the keyword :ref:`MOLTRA_.HTSORT`.
In this step the 1HT entities on disk are rearranged from
their "natural" order (arising from the two-electron integral algorithm)
to an order which will minimize the amount of I/O operations taking
place at the 2HT step.

 *HT sort disabled (Default):*

::

    .SCHEME
     6
    .HTSORT
     0

Sorting of HT intermediates is
disabled by default. Without such sorting, scheme 6's performance can be
negatively affected on computer systems where reading data is relatively
expensive (one example is for IBM p4/p5/p6 systems and the XLF compiler)
even for relatively modest calculations Dirac-Coulomb calculations
(where LL and SL integrals are considered). However, in two-component
calculations (where only LL integrals are transformed and therefore
there are much less 1HT entities on disk) of similar or larger size, in
terms of active space, good performance can still be achieved in the
same architecture.

 *HT sort (strategy #1) enabled:*

::

    .SCHEME
     6
    .HTSORT
     1

This sort strategy minimizes the
amount of read/write operations by taking up as much 1HT data in memory
and sorting it before writing it to disk again. This provides efficiency
for systems (such as the IBM systems mentioned above) where reading is
expensive, at the expense of larger memory usage (rule of thumb: the
extra memory required \*per MPI process\* will be about 1-2% of the
\*total disk\* space taken up by the 1HT quantities - so if that is 1Tb,
the extra memory used will be of of about 1-2Gb; see the calculation
output for more details).

keyword(INTFLG)

Specify which classes of integrals should contribute to the transformed
two-electron integrals and the effective (core) Hamiltonian. Default is
to follow the definitions given in \*\*HAMILTONIAN. If you want
different integral classes in the 2- and 4-index transformation, then
use the :ref:`MOLTRA_.INTFL2` and :ref:`MOLTRA_.INTFL4` keywords below
(default: :ref:`HAMILTONIAN_.INTFLG` from :ref:`**HAMILTONIAN`).

keyword(INTFL2)

Specify what two-electron integrals to include in the 2-index
transformation ,i.e., in the two-electron part of the Fock core matrix
(default: :ref:`HAMILTONIAN_.INTFLG` from :ref:`**HAMILTONIAN`).

keyword(INTFL4)

Specify what two-electron integrals to include in the 4-index transformation
(default: :ref:`HAMILTONIAN_.INTFLG` from :ref:`**HAMILTONIAN`).

keyword(CORE)

Specify the frozen core spinors,

For each fermion irrep, give an :ref:`orbital_strings`
of active orbitals.

*Default:* All occupied orbitals that are not active.

keyword(CORE2)

Specify the frozen core spinors for open (average of configuration) shell,

In the first line give the number of electrons in the open shell.
For each fermion irrep, give an :ref:`orbital_strings`
of active orbitals. ATTENTION: A continous set of orbitals is assumed.

*Default:* No frozen open shell. 

keyword(POSITRONS)

Allow negative energy ("positronic") spinors in active set. If this
keyword is not specified then only electronic spinors will be included
and any positronic spinors specified with :ref:`MOLTRA_.ACTIVE` are
ignored.

keyword(NO2IND)

Skip the 2-index transformation of the effective Fock matrix.

keyword(NO4IND)

Skip the 4-index transformation.

keyword(SCREEN)

Screening threshold in 4-index transformation (a negative value disables
screening).

*Default:*

::

    .SCREEN
     1.0D-14

keyword(ASCII)

Write integrals to the ASCII file MO\_integrals.txt to facilitate
interfaces with correlation codes. Warning: only meant for initial
interfacing as this format is very inefficient.

*Default:* Do not write the ASCII file.

keyword(RCORBS)

Recanonize orbitals before transforming.

*Default:* Do not recanonize orbitals before transforming.

**Programmers options**

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(2INDEX)

Ranges of active orbitals in 2-index transformation module specified for
index 1 and 2 and fermion irrep 1 and 2.

*Default:* Ranges set by :ref:`MOLTRA_.ACTIVE`.

keyword(4INDEX)

Ranges of active orbitals in 4-index transformation module specified for
index 1 to 4 and fermion irrep 1 and 2.

*Default:* Ranges set by :ref:`MOLTRA_.ACTIVE`.

keyword(MDCINT)

Write 4-index transformed integrals to MDCINT file (default).

keyword(NOMDCI)

Skip writing 4-index transformed integrals to MDCINT file. Instead, keep 4IND* files
and write file 4INDINFO, which contains info about how to read the 4IND* files.

keyword(SCATTER)

Scatter 4-index transformed integrals to all nodes (default).

keyword(NOSCAT)

Do not scatter 4-index transformed integrals to all nodes

keyword(MOFILE)

Name of file with MO coefficients to be used for integral transformation.
Valid options: "DFCOEF", "KRMCSCF", "KRMCOLD", and "UNKNOWN".
The "UNKNOWN" option means use the first one found of "KRMCSCF", "KRMCOLD", "DFCOEF" in this order.

*Default:*

::

    .MOFILE
    UNKNOWN


keyword(PAR4BAS)

Specify number of batches for each node in 4-index transformation. Only
working for parallel distributions.

*Default:*

::

    .PAR4BAS
     -1

star(PRPTRA)

Property integrals transformation
---------------------------------

This section defines the property operators that should to be transformed. Such
integrals are only used in an experimental section of :ref:`**RELCC`, and this
subsection can thus usually be omitted.

**Programmers options**

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(OPERATOR)

Specification of general one-electron operators.

See the :ref:`one_electron_operators` section for more information and explicit
examples.

