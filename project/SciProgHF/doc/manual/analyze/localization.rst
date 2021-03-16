:orphan:
 

star(LOCALIZATION)

Performs molecular orbitals localization using the Pipek-Mezey
criterion. The implementation uses an exponential parametrization and a
Trus region minimization method instead of the Edminston-Ruedenberg
algorithm. It works with the *C*\ :sub:`1`\ symmetry and for the occupied orbitals.


keyword(PRJLOC)

Localization from projection analysis instead of Mulliken analysis. [Default: False]::

  .PRJLOC


keyword(OWNBAS)

Calculate fragments in their own basis. For now it is the only option which can be
used with projection analysis option (.PRJLOC) [Default: False]::

  .OWNBAS

keyword(VECREF)

First give number of fragments to project onto. Then for each fragment give
filename of MO coefficients and string specifying which orbitals will be
used in the projection (:ref:`orbital_strings`)::

  .VECREF
   3
   DFO1AT
   1..6
   DFH1AT
   1
   DFH2AT
   1


keyword(HESLOC)

Define model of the Hessian [Default: FULL].

Use full Hessian::

  .HESLOC
   FULL

Use diagonal approximation of the Hessian. If using this approximation,
there is danger of converging to the stationary point instead of minima::

  .HESLOC
   DIAG

First the diagonal approximation of the Hessian will be used 
and after reaching convergence criterion the localization will 
switch to the construction of the full Hessian::

  .HESLOC
   COMB

keyword(CHECK)

This keyword can be used only in combination with .HESLOC = COMB.
The full Hessian will be calculated only at the last iteration.
This way one can check if the procedure converged to the minimum.


keyword(THFULL)

Convergence threshold for the localization process. It applies on the functional
value of the localization criteria when full Hessian is calculated. It is used in 
.HESLOC = FULL scheme and in the second part of the .HESLOC = COMB scheme.
When not specified gradient criterion or criterion of the number of negative
eigenvalues equals zero will be used::

  .THFULL
   1.0D-13


keyword(THGRAD)

Convergence threshold for the localization process. It applies on the gradient
of the localization functional. It is used in .HESLOC = FULL or DIAG scheme and in the
second part of the .HESLOC = COMB scheme. When not specified the functional value 
criterion or criterion of the number of negative eigenvalues equals zero will be used::

  .THGRAD
   1.0D-17


keyword(THDIAG)

Convergence threshold for the localization process. It applies on the functional 
value of the localization criteria when diagonal Hessian is calculated.
It can be used with the .HESLOC = DIAG keyword.
If .HESLOC = COMB it is used to switch from the first to the second stage of the
convergence process::

  .THDIAG
   1.0D-05


keyword(SELECT)

Select subset of MOs for localization (use :ref:`orbital_strings` syntax). 
Currently only Pipek-Mezey localization is implemented, therefore localize
only occupied orbitals. Localization of virtual orbitals is not recommended.

::

  .SELECT
   1..5


keyword(MAXITR)

Maximum number of iterations. Default::

  .MAXITR
   100


keyword(PRINT)

Set print level. Default::

  .PRINT
   1

