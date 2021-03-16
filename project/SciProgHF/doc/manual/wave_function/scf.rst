:orphan:
 

star(SCF)

This section gives directives for Hartree-Fock and Kohn-Sham calculations.
Kohn-Sham calculations are activated by invoking the keyword
:ref:`HAMILTONIAN_.DFT` under :ref:`**HAMILTONIAN`.

Open-shell calculations correspond to either
average-of-configurations (:ref:`SCF_.AOC`) (see :ref:`aoc` for theory) or fractional occupation (:ref:`SCF_.FOCC`).
The former is the default for
Hartree-Fock calculations, whereas the latter is default for
Kohn-Sham calculations. Note that average-of-configurations
Kohn-Sham calculations are not well defined.

Occupation
==========

keyword(CLOSED SHELL)

For each fermion irrep give the number of closed shell electrons.

The specification of the closed shell electrons is simple. For
symmetry groups *without* inversion symmetry, there is only one
fermion irrep, and you need only to specify the number of
electrons.

For symmetry groups *with* inversion symmetry, you need to specify
the distribution of the electrons in the two fermion irreps :cite:`Saue2000`.

keyword(OPEN SHELL)

Specification of open shell(s).

For each open shell give the number of electrons and the number of
active spinors.

*Short example:*

::

    .OPEN SHELL
     1
     5/0,6

1 open shell with 5 electrons in 6 spinors (= 3 Kramers pairs) in
irrep 2 (the *ungerade* one). Thus, the fractional occupation is 5/6.

The open shell module in DIRAC is based on
average-of-configurations :cite:`Thyssen1998` .
The simplest case is one electron in two spinors (= one Kramers
pair). For this special case the average-of-configuration
calculation gives the same result as the usual restricted
open-shell Hartree-Fock. For all other cases the calculation gives
the average energy of many states.

Note that the order of closed and open shells are assumed to be as
in the following scheme:

**Virtual orbitals**
Not occupied

...
**Open shell 2**
Fractionally occupied

**Open shell 1**
Fractionally occupied

**Closed shell**

Doubly occupied;
that is, the lowest molecular orbitals are doubly occupied, the
next ones are occupied with the electrons of open shell no. 1,
etc.

Other orderings can be achieved by using :ref:`WAVE_FUNCTION_.REORDER MO`
and :ref:`SCF_.OVLSEL`.

To get the energies of the individual states present in the
average-of-configurations, specify :ref:`WAVE_FUNCTION_.RESOLVE` (see also
:ref:`*RESOLVE`).

To get the energies of (some) of the individual states present in
the average of configurations, you can use the
:ref:`GOSCIP`, the :ref:`DIRRCI`, or the :ref:`*LUCITA`.

keyword(BOSONS)

Occupation of boson irreps in spin-free calculation. For example,
for the D2h symmetry eight numbers in subsequent line, for the C2v
symmetry there four occupation numbers in line.

keyword(MJSELE)

In the case of linear supersymmetry give occupation for each :math:`M_J` value.

keyword(KPSELE)

In the case of **atomic supersymmetry** give occupation for each :math:`\kappa` - value.
This option also works in linear supersymmetry when a single atomic center is combined with a ghost atom.

The format of :ref:`SCF_.KPSELE` is illustrated for the case of Uranium (:math:`[Rn]5f^36d^17s^2`). We first provide the
usual specification of closed and open shells

::

  .CLOSED SHELL
  44 44
  .OPEN SHELL
   2
   3/0,14
   1/10,0

The closed shells are those of Radon as well as the outer :math:`7s` shell of Uranium. However, their presence will
lead to convergence problems because by default orbitals are ordered according to energy, but also with closed shells
before open ones. This means that the outer :math:`7s` shell of Uranium will end up amongst the open-shell orbitals, whereas
some :math:`6d` orbitals end up being defined as closed shell, thus creating havoc. This is completely avoided by in addition
giving the occupation in terms of :math:`\kappa` - values as shown below

::
   
   .KPSELEC
   7                              # Number of the Kappa-splitted orbital
    -1   1  -2   2  -3   3  -4    # Values of Kappa: s   p-  p+  d-  d+  f-  f+
    14  10  20  12  18   6   8    # Number of the electrons in the closed orbitals
     0   0   0   0   0   6   8    # Number of the orbitals of open shell 1 (5f^3)
     0   0   0   4   6   0   0    # Number of the orbitals of open shell 2 (6d^1)

   
keyword(AOC)

Average-of-configuration calculation (default for open-shell
Hartree-Fock).

keyword(FOCC)

Fractional occupation (default for open-shell
Kohn-Sham) CLARIFY

:ref:`SCF_.FOCC` calculations are less memory-intensive than
:ref:`SCF_.AOC` calculations. In the latter case one additional
AO-Fock matrix is generated for each open shell.

:ref:`SCF_.FOCC` calculations are therefore an interesting option
for generating start orbitals for MCSCF as well as initial
convergence in open-shell Hartree-Fock.

keyword(AUTOCC)

Program is allowed to change occupation during SCF cycles. This is
deactivated by default. However, the program will still try to do
an automatic initial occupation if neither
:ref:`SCF_.CLOSED SHELL` nor :ref:`SCF_.OPEN SHELL` is given.

Trial function
==============

An SCF-calculation (HF or DFT) may be initiated in four different
ways:


-  Using **MO coefficients** from a previous calculation.
-  Using coefficients obtained by diagonalization of the
   one-electron Fock matrix: the **bare nucleus approach**.
-  The **corrected bare nucleus approach**. There are two flavors:

   - SCRPOT: sum of atomic LDA potentials, generated by :cite:`GRASP` (default)
   - BNCORR: bare nucleus potential corrected with screening factors based on Slaters rules

-  Using the **two-electron Fock matrix** from a previous
   calculation; this may be thought of as starting from a converged
   SCF potential.
-  Using an **atomic start** based on densities from atomic SCF runs for the individual centers, see e.g. :cite:`vanLenthe2006` .
-  Using an **extended Hückel start** based on *atomic fragments*

The default is to start from MO coefficients if the file DFCOEF is
present. Otherwise the corrected bare nucleus approach (SCRPOT) is followed.
In all cases linear dependencies are removed in the zeroth
iteration.

keyword(ATOMST)

Start first SCF iteration with a molecular density matrix constructed from
atomic densities.  The keyword ``ATOMST`` is followed by input for each atomic
type. The details, orbital strings (see :ref:`orbital_strings` for the syntax)
and occupation, usually correspond to those of the atomic runs, but the user
may modify this at will.  The syntax is explained in the parenthesis "" for
each atomic type but we highly recommend to carefully check the tutorial
example :ref:`atomic_start_guess`. Please note that the order of atoms
corresponds to the order they appear in the molecule file.

::

  .ATOMST
  "SCF coefficients file name (6 characters)" "integer specifying # of occupation patterns, here: 2"
  orbital occupation string #1 for atomic type 1
  occupation (real*8 value in the range of 0.0d0 - 1.0d0)
  orbital occupation string #2 for atomic type 1
  occupation (real*8 value in the range of 0.0d0 - 1.0d0)
  "SCF coefficients file name (6 characters)" "integer specifying # of occupation patterns, here: 1"
  orbital occupation string #1 for atomic type 2
  occupation (real*8 value in the range of 0.0d0 - 1.0d0)
  ...

keyword(AD HOC)

Start first SCF iteration with orbitals generated from an extended Hückel calculation using pre-calculated orbitals 
for each constituent atomic type of the molecule. 

::

   .AD HOC
   "SCF coefficients file name (6 characters)" 
   orbital occupation string for atomic type 1
   "SCF coefficients file name (6 characters)" 
   orbital occupation string for atomic type 2
   ...

The order of atomic types follows that of the input. Presently, this functionality only works without symmetry.

keyword(HUCPAR)
Modify the the Wolfsberg-Helmholtz constant :math:`K` of the extended Hückel calculation.
*Default:* 1.75

keyword(SCRPOT)
Default start procedure: use sum of atomic potentials generated using LDA on Hartree-Fock
densities (numerical 4C, generated by :cite:`GRASP`).

keyword(BNCORR)
Old version of the start potential.
Two-electron repulsion is estimated via nuclear-attraction type integrals:

.. math::

 \langle X_{A} \vert \sum_{C} \frac{-Z_{C} \cdot \sum_j a_j e^{(-\alpha^C_j r_{C}^{2})}}{r_{C}} \vert X_{B} \rangle, \ \ \ X = L,S

The coefficients *a*\ :sub:`*j*`\  and the exponents :math:`\alpha^{C_j}`
in this expression are chosen according to Slater's rules to obtain
an approximate atomic electronic density for the initial guess. For
example, with one heavy element and without this correction (that
is, with the bare nucleus Hamiltonian) all electrons will end up on
that heavy element in the initial guess!

keyword(NOBNCR)

Switch off all bare nucleus corrections (SCRPOT or BNCORR).

keyword(BNC_FORCE)

Force employing the bare nucleus correction (BNC). This keyword is worth when the calculated system
is highly positively charged what makes (from defined charge value) switching off the default BNC.
The BNC can help to achieve better convergence also for non-neutral systems.

keyword(FOMOUT)

Print Fock MO matrices for diagonalization (according to symmetries) into own formatted files. 
Programmer's option suitable for testing. Only in for the linear symmetry.

keyword(TRIVEC)

Start SCF-iterations from the vector file.

keyword(TRIFCK)

Start SCF-iterations from the two-electron Fock matrix from
previous calculation (stored on file DFFCK2).

keyword(MOSTART)

This keyword collects most start guess possibilities. It is followed by a second line specifying start guess. The available options are:

- Bare nucleus start::

    .MOSTART
    BARNUC

- Bare nucleus correction, using sum of atomic potentials generated using LDA on Hartree-Fock densities (numerical 4C, generated by :cite:`GRASP`) ::

    .MOSTART
    SCRPOT

- Bare nucleus potential corrected with screening factors based on Slater's rules ::

    .MOSTART
    BNCORR

- Start SCF-iterations from the vector file DFCOEF ::

    .MOSTART
    TRIVEC

- Start SCF-iterations from the two-electron Fock matrix from previous calculation (stored on file DFFCK2) ::

    .MOSTART
    TRIFCK


Convergence criteria
====================

Three different criteria for convergence may be chosen:


-  The norm of the DIIS error vector
   :math:`\mathbf{e} = [\mathbf{F}, \mathbf{D}]` (in MO basis). This
   corresponds to the norm of the electronic gradient and is the
   recommended convergence criterion. When you are only interested in
   the energy :ref:`SCF_.EVCCNV` = 1.0D-5 is usually sufficient.
   For properties and correlated methods you should converge to
   :ref:`SCF_.EVCCNV` = 1.0D-9. Large negative energy eigenvalues
   lead to a loss of precision that might lead to convergence
   problems. Remember also that a too loose screening threshold (too
   many integrals neglected) will hinder convergence. You should
   modify
   :ref:`TWOINT_.SCREEN` under
   :ref:`*TWOINT` if you
   modify :ref:`SCF_.EVCCNV` or one of the other two convergence
   criteria.
-  The difference in total energy between two consecutive
   iterations.
-  The largest absolute difference in the total Fock matrix between
   two consecutive iterations.

The change in total energy is approximately the square of the
largest element in the error vector or the largest change in the
Fock matrix. The default is convergence on electronic gradient with
1.0D-6 as threshold. Alternatively, the iterations will stop at the
maximum number of iterations.

Sometimes it may happen that the specified convergence criterion is
too tight for the given basis set and/or other input parameters. In
this case one needs to decide whether one should proceed with
post-HF steps (like correlation calculations) or not. The program
decides this by looking at a secondary convergence criterion that
gives the *allowed* convergence. This value is by default the same
as first or *desired* convergence criterion but can be made lower
to make sure that a calculation does not abort when the convergence
is slightly above threshold.

For more detailed help see SCF help on convergence troubleshooting  and related pages.

keyword(MAXITR)

Maximum number of SCF iterations.

*Default:*

::

    .MAXITR
     50

When restarting SCF itrations from previous molecular orbitals file
(DFCMO or formatted DFPCMO), we recommend to decrease maximum
number of iterations together with readjusting desired and allowed
convergence thresholds. By properly set desired and allowed
thresholds one can have exact number of iterations specified by
*.MAXITR*.

keyword(EVCCNV)

Converge on error vector (electronic gradient).

*2 (real) Arguments:*

::

    .EVCCNV
     desired threshold allowed threshold

keyword(ERGCNV)

Threshold for convergence on total energy.

*2 (real) Arguments:*

::

    .ERGCNV
     desired threshold allowed threshold

keyword(FCKCNV)

Converge on largest absolute change in Fock matrix.

*2 (real) Arguments:*

::

    .FCKCNV
     desired threshold allowed threshold

Note that the *allowed* threshold may be omitted. It is then made
equal to the *desired* threshold.

Convergence acceleration
========================

It is imperative to keep the number of SCF iterations at a minimum.
This may be achieved by convergence acceleration schemes:


-  **Damping:** The simplest scheme is damping of the Fock matrix
   that may remove oscillations. In :math:`n + 1` iteration the Fock matrix
   to be diagonalized is:
   :math:`\mathbf{F}\' = (1-c) \mathbf{F}_{n+1} + c \mathbf{F}_n`,
   where :math:`c` is the damping factor.

-  **DIIS:** Direct inversion of iterative subspaces, Refs. :cite:`Pulay1980` , :cite:`Pulay1982` , :cite:`Hamilton1986`,
   may be thought of as generalized damping involving Fock matrices from many iterations.
   Damping factors are obtained by solving a simple matrix equation
   involving the B-matrix constructed from error vectors (approximate
   gradients). Linear dependent columns in the B-matrix is removed.

In DIRAC DIIS takes precedence over damping.

keyword(DIISTH)

Change the default convergence threshold for initiation of DIIS,
based on largest element of error vector.

*Default:*

::

    .DIISTH
     a very large number

keyword(MXDIIS)

Maximum dimension of B-matrix in the DIIS module.

*Default:*

::

    .MXDIIS
     10

keyword(DIISMO)

Activate DIIS in orthogonal basis (MO) with the error vector as
described above.

keyword(DIISAO)

Activate
`DALTON <http://www.kjemi.uio.no/software/dalton/dalton.html>`_-like
DIIS using AO-basis. The error vector is

.. math::

 {\mathbf{e}}={\mathbf{f}}-{\mathbf{f}}^{\dagger} 

where the term :math:`\mathbf{f}` is given by

.. math::

   \mathbf{f}=\mathbf{C}^{\dagger}\cdot\mathbf{S}_{AO}\cdot \left[ \mathbf{D}^{C}_{AO}\cdot\mathbf{F}^{D}_{AO}+\sum_{O\in\mathcal{O}}f_{O}\cdot\mathbf{D}^{O}_{{AO}}\cdot ( \mathbf{F}^{D}_{{AO}}+(a_{O}-1)\mathbf{Q}^{V,O}_{{AO}} ) \right] \mathbf{C}

keyword(NODIIS)

Do not perform DIIS. The default is to activate
:ref:`SCF_.DIISMO` for closed-shell calculations, and to
activate :ref:`SCF_.DIISAO` for average-of-configurations
calculations.

keyword(DAMPFC)

Change the default damping factor.

*Default:*

::

    .DAMPFC
     0.25

keyword(NODAMP)

Do not perform damping of the Fock matrix. Damping is activated by
default, but DIIS takes precedence. In case all columns in the
B-matrix is removed by linear dependency, damping is activated.

Level shifts
============

keyword(LSHIFT)

Activate level shift (for virtual orbitals). Followed by a real
argument (level shift).

keyword(OLEVEL)

Activate level shift (for open-shell orbitals). Followed by a real
argument (level shift), one line for each open shell {Please give
example}.

keyword(OPENFACTOR)

Change the default factor on an open-shell diagonal contribution to the Fock matrix (see :ref:`aoc` for theory). A factor of one corresponds to a Koopmans interpretation of the orbital energies. However, experience shows convergence can be improved by tuning this factor. DIRAC therefore presently employs a default factor of 1/2.

*Default:*

::

    .OPENFAC
     0.5

2nd-order optimization
======================

keyword(2NDOPT)

The default SCF of DIRAC uses only gradient information. By adding this keyword 2nd-order optimization, 
using both gradient and Hessian information, is activated in case the regular SCF does not converge.
This scheme is computationally more expensive and so far only available for closed-shell Hartree-Fock.


State selection
===============

Convergence can be improved by selection of vectors based on
overlap with vectors from a previous iteration. This method may
also be used for convergence to some excited state.

If dynamic overlap selection is used, the vector set from the
previous iteration is used as the criterion. For the first
iteration either restart vectors or vectors generated by the bare
nucleus approach (not*recommended) are used.*

If :ref:`SCF_.NODYNSEL` is given, either the restart vectors
or the bare nucleus vectors are used, i.e. the overlap selection
vectors are *not* updated in each iteration. Please note, that
overlap selection based on vectors from the bare nucleus approach
is not recommended.

Overlap selection is very useful together with
:ref:`WAVE_FUNCTION_.REORDER MO`.
This will reshuffle the vectors within the restart coefficients.

Example: First one might do a open shell calculation on Boron, this
would give the *P*\ :sub:`1 / 2`\  state. But if we restart on the
*P*\ :sub:`1 / 2`\  coefficients, interchange the
*p*\ :sub:`1 / 2`\  with the *p*\ :sub:`3 / 2`\  orbitals, and
request overlap selection, we can converge to the
*P*\ :sub:`3 / 2`\  state.

There also exists a keyword for reordering the converged SCF
orbitals. This is useful for reordering the orbitals for the
4-index transformation and subsequent correlation calculations
(CCSD, CI etc.) (see :ref:`WAVE_FUNCTION_.POST SCF REORDER MO`).

keyword(OVLSEL)

Activate dynamic overlap selection. The default is no overlap
selection.

keyword(NODYNSEL)

No dynamic update of overlap selection vectors. The default is
dynamic update.

 .. note:: Overlap selection is nowadays marketed hard as MOM (Maximum Orbital Method, see :cite:`Gilbert_JPCA2008`), but this method has been included in DIRAC for at least two decades and goes back to the pioneering work of `Paul Bagus <http://cascam.unt.edu/people/psbagus.htm>`_ It was used in :cite:`Bagus_JCP1971`, but not reported explicitly. However, it is for instance documented in the `1975 manual of the ALCHEMY program <http://k-sek01.t-komazawa.ac.jp/msekiya/alchemy/scfm.pdf>`_ (On pdf page 15 you find a description of keyword MOORDR using a "maximum overlap criterion").


Iteration speedup
=================

The total run time may be reduced significantly by reducing the
number of integrals to be processed in each iteration:


-  **Screening on integrals:** Thresholds may be set to eliminate
   integrals below the threshold value, see :cite:`Saue1997`.
   . The threshold for LL
   integrals is set in the basis file, but this threshold may be
   adjusted for SL and SS integrals by threshold factors set in the
   :ref:`\*\*INTEGRALS` section.

-  **Screening on density:** In direct mode further reductions are
   obtained by screening on the density matrix as well, see Ref. :cite:`Saue1997`.
   This becomes even more
   effective if one employs *differential densities*, that is
   :math:`\Delta \mathbf{D} = \mathbf{D}_{n+1} - \alpha \cdot \mathbf{D}_n`.
   The default value for :math:`alpha` is
   :math:`\alpha=\frac{ \mathbf{D}_{n+1} \cdot \mathbf{D}_n }{ \mathbf{D}_{n}  \cdot \mathbf{D}_n }`
   which corresponds to a Gram-Schmidt orthogonalization. As SCF
   converges, :math:`alpha` goes towards 1, but  :math:`alpha` can also explicitly be set equal
   to 1 with :ref:`SCF_.FIXDIF`.


-  **Neglect of integrals:** The number of integrals to be
   processed may be reduced even further by adding SL and SS integrals
   only at an advanced stage in the SCF-iterations, as determined
   either by the number of iterations or by energy convergence. The
   latter takes precedence over the former.

keyword(NODSCF)

Do not perform SCF-iterations with differential density matrix.

*Default:* Use differential density matrix in direct SCF.

keyword(FIXDIF)

Set :math:`alpha` equal to 1.

keyword(CNVINT)

Set thresholds for convergence before adding SL and SS/GT integrals to
SCF-iterations. 

*2 (real) Arguments:*

::

    .CNVINT
     CNVINT(1) CNVINT(2)

keyword(ITRINT)

Set the number of iterations before adding SL and SS/GT integrals to
SCF-iterations.

*Default:*

::

    .ITRINT
     1 1

keyword(INTFLG)

Specify what two-electron integrals to include (see
:ref:`HAMILTONIAN_.INTFLG` under :ref:`**HAMILTONIAN`).

*Default:* :ref:`HAMILTONIAN_.INTFLG` from :ref:`**HAMILTONIAN`.

Output control
==============

keyword(PRINT)

General print level for the SCF method. For instance, value of 2
prints eignevalues during each iteration.

*Default:*

::

    .PRINT
     0

keyword(EIGPRI)

Controls the print-out of positive energy and negative energy
eigenvalues (1 = on; 0 = off).

*Default:* Only the positive energy eigenvalues are printed.

::

    .EIGPRI
     1 0

Eliminating/freezing orbitals
=============================

In studies of electronic structure it may be of interest to
eliminate or freeze certain orbitals. This option is furthermore
useful for convergence, in particular to excited electronic states.
A simple case is the thallium atom. The ground state
\ :sup:`2`\ *P*\ :sub:`1 / 2`\  has the electronic configuration
[Xe]|4f^{14}5d^{10}6s^{2}6p^1\_{1/2}|. The first excited state
\ :sup:`2`\ *P*\ :sub:`3 / 2`\  with the electronic configuration
[Xe]4f^{14}5d^{10}6s^{2}6p^1\_{3/2}| can easily be accessed by
first calculating the ground state, then eliminating the
6p^1\_{1/2}| from the ensuing calculation. In the final
calculation the excited state is relaxed using overlap selection.
The use of frozen orbitals is demonstrated in test/33.frozen: When
the geometry of the water molecule is optimized with the oxygen 1s
and 2s orbitals frozen, a bond angle of 96.242 degrees is found,
contrary to the 90 degrees one might have expected when s-p
hybridization is thus blocked :cite:`Kutzelnigg1984`.

The elimination of orbitals is achieved by projecting the selected
orbitals out of the transformation matrix to orthonormal basis. The
selected orbitals can be expressed either in the full molecular
basis or in the basis set of the chosen fragment. In the latter
case, the same set of fragment orbitals can in the case of atomic
fragments be used at different molecular geometries. One may even
perform a geometry optimization, but only using the numerical
gradient. When freezing orbitals the selected orbitals are first
eliminated from the transformation to orthonormal basis, but then
reintroduced in the backtransformation step. They will appear in
the output with zero orbital eigenvalues. Note that when freezing
orbitals the orbitals to be eliminated must be specified as well.
The frozen orbitals must be a subset of the eliminated orbitals.

Fragments are defined with respect to the list of
symmetry-independent atoms appearing in the DIRAC basis file:
Consider the water molecule in the full *C*\ :sub:`2*v*`\
symmetry. Then there are two fragments: the oxygen atom and the
H\ :sub:`2`\  moiety. However, with no symmetry there will be three
fragments: the oxygen atom and the two hydrogen atoms. At the
moment there are no orthonormalization of fragments on different
fragments and so in practice one should only use orbitals from one
fragment.

keyword(PROJECT)

Eliminate orbitals by projecting them out of the transformation
matrix to orthonormal basis.

*Arguments:* Number of fragments (NPRJREF).

Then for each fragment read name PRJFIL of the coefficient file
followed by the number of symmetry-independent nuclei in this
fragment followed by an
:ref:`orbital_strings` of
selected orbitals for each fermion irrep.

::

    DO J = 1,NPRJREF
      READ(LUCMD,'(A6)') PRJFIL(J)
      READ(LUCMD,'(I6)') NPRJNUC(J)
      READ(LUCMD,'(A72)') (VCPROJ(I,J),I=1,NFSYM)
    ENDDO

keyword(PRJTHR)

Smallest norm accepted when eliminating orbitals.

*Default:*

::

    .PRJTHR
     1.0D-10

keyword(OWNBAS)

Eliminated/frozen orbitals are given in the fragment basis. Note
that the list of fragments is assumed to follow the list of
symmetry independent nuclei in the DIRAC basis file.

keyword(FROZEN)

Freeze orbitals. This keyword must only be used in conjunction with the keyword
:ref:`SCF_.PROJECT`. The latter keyword eliminates selected orbitals from the
variational space. From the list of eliminated orbitals, the user can select those
that should be added to the coefficient array after diagonalization of the Fock matrix
in orthonormal basis. To do so, the user must, for each fermion irrep, give an :ref:`orbital_strings`
giving the position of the eliminated orbital in the coefficient array. If the value zero
is given, it means that this eliminated orbital is definitely out and not kept frozen.

