:orphan:
 

star(OPTIMIZE)

Geometry optimization directives.

This section controls the geometry optimization. The geometry optimization
algorithm is based on the one of `Dalton <http://www.daltonprogram.org>`_. The
directives are therefore similar except for these few changes:

-  No second-order algorithms are available since the molecular Hessian
   is not implemented.
-  A few new keywords have been introduced (see below).

If :ref:`**PROPERTIES` or :ref:`**ANALYZE` is specified in the job input section
together with :ref:`*OPTIMIZE`, then the property or analysis module is called
in each optimization iteration and at the converged geometry.

Depending on convergence criteria for energies and gradients the values of the
thresholds may be automatically adjusted.


keyword(NO SKIP)

Don't skip SSLL or SSSS gradient contributions when small.
This keyword forces the LS and SS two-electron gradient to always be evaluated
in all geometry iterations (depending on the integral flag).

The SL and SS two-electron integral contributions to the gradient are normally
skipped if their contribution is estimated to be small. An "empirical" estimate
of the norms of the LS and SS two-electron gradients based on the norm of the
LL two-electron gradient.  This trick is by default activated in a geometry
optimization, since when the current geometry is far away from the equilibrium,
e.g. the norm of the gradient is, say, 1.0, then there is no need to calculate
the LS and/or SS two-electron gradient because they have a norm of, say, 0.001.


keyword(NUMGRA)

Force the use of numerical gradient in geometry optimization.


keyword(TWOGRD)

Include LL, SL, and SS integral contributions to the gradient
(1 = on; 0 = off). The default is to turn all on::

  .TWOGRD
   1 1 1


keyword(BAKER)

Baker's convergence criteria :cite:`Baker1993` will be used.
The minimum is then said to be
found when the largest element of the gradient vector (in Cartesian or redundant
internal coordinates) falls below 3.0D-4  and either the energy change from the last
iteration is less than 1.0D-6 or the largest element of the predicted step vector is
less 3.0D-4.


keyword(GRADIENT)

Convergence threshold for the gradient. The default is 1.0D-5.


keyword(ENERGY)

Convergence threshold for the energy. The default is 1.0D-6.


keyword(STEP THRESHOLD)

Convergence threshold for the step. The default is 1.0D-5.


keyword(PRINT)

Set print level for this module. Read one more line containing print level. Default
value is 0, any value higher than 12 gives debugging level output.


keyword(MAX IT)

Read the maximum number of geometry iterations. The default value is 25.


keyword(TRUSTR)

Set initial trust radius for calculation. This will also be the maximum step length for
the first iteration. The trust radius is updated after each iteration depending on the
ratio between predicted and actual energy change. The default trust radius is 0.5 a.u.


keyword(TR FAC)

::

 READ(LUCMD,*) TRSTIN, TRSTDE

Read two factors that will be used for increasing and decreasing the trust radius
respectively. Default values are 1.2 and 0.7.


keyword(TR LIM)

::

 READ(LUCMD,*) RTENBD, RTENGD, RTRJMN, RTRJMX

Read four limits for the ratio between the actual and predicted energies. This ratio
indicates how good the step is?that is, how accurately the quadratic model describes
the true energy surface. If the ratio is below RTRJMN or above RTRJMX, the step is
rejected. With a ratio between RTRJMN and RTENBD, the step is considered bad an
the trust radius decreased to less than the step length. Ratios between RTENBD and
RTENGD are considered satisfactory, the trust radius is set equal to the norm of the
step. Finally ratios above RTENGD (but below RTRJMX) indicate a good step, and the
trust radius is given a value larger than the step length. The amount the trust radius
is increased or decreased can be adjusted with .TR FAC. The default values of RTENBD,
RTENGD, RTRJMN and RTRJMX are 0.4, 0.8, -0.1 and 3.0 respectively.


keyword(MAX RE)

::

 READ(LUCMD,*) MAXREJ

Read maximum number of rejected steps in each iterations, default is 3.


keyword(NOTRUS)

Turns off the trust radius, so that a full Newton step is taken in each iteration.
This should be used with caution, as global convergence is no longer guaranteed. If
long steps are desired, it is safer to adjust the initial trust radius and the limits for
the actual/predicted energy ratio.


keyword(CONDIT)

::

 READ (LUCMD,*) ICONDI

Set the number of convergence criteria that should be fulfilled before convergence
occurs. There are three different convergence thresholds, one for the energy, one for
the gradient norm and one for the step norm. The possible values for this variable
is therefore between 1 and 3. Default is 2. The three convergence thresholds can be
adjusted with the keywords :ref:`OPTIMIZE_.ENERGY`, :ref:`OPTIMIZE_.GRADIENT` and 
:ref:`OPTIMIZE_.STEP THRESHOLD`.


keyword(NOBREA)

Disables breaking of symmetry. The geometry will be optimized within the given
symmetry, even if a non-zero molecular Hessian index is found. The default is to let
the symmetry be broken until a minimum is found with a molecular Hessian index of
zero. This option only has effect when second-order methods are used.


keyword(SP BAS)

::

 READ(LUCMD,*) SPBSTX

Read a string containing the name of a basis set. When the geometry has converged,
a single-point energy will be calculated using this basis set.


keyword(PREOPT)

::

 READ (LUCMD,*) NUMPRE
 DO I = 1, NUMPRE
 READ (LUCMD,*) PREBTX(I)
 END DO

First we read the number of basis sets that should be used for preoptimization, then
we read those basis set names as strings. These sets will be used for optimization in
the order they appear in the input. One should therefore place the smaller basis at
the top. After the preoptimization, optimization is performed with the basis specified
in the molecule input file.


keyword(VISUAL)

Specifies that the molecule should be visualized, writing a VRML file of the molecular
geometry. No optimization will be performed when this keyword is given. See also
related keywords .VR-BON, .VR-COR, .VR-EIG and .VR-VIB.


keyword(VRML) 

Specifies that the molecule should be visualized. VRML files describing both the
initial and final geometry will be written (as initial.wrl and final.wrl). The file
final.wrl is updated in each iteration, so that it always reflects the latest geometry.
See also related keywords .VR-BON, .VR-COR, .VR-EIG and .VR-VIB.


keyword(SYMTHR)

::

 READ(LUCMD,*) THRSYM

Determines the gradient threshold (in a.u.) for breaking of the symmetry. That is,
if the index of the molecular molecular Hessian is non-zero when the gradient norm
drops below this value, the symmetry is broken to avoid unnecessary iterations within
the wrong symmetry. This option only applies to second-order methods and when the
keyword .NOBREA is not present. The default value of this threshold is 0.005.


keyword(TRSTRG)

Specifies that the level-shifted trust region method should be used to control
the step. This is the default, so the keyword is actually redundant at the moment.
Alternative step control methods are :ref:`OPTIMIZE_.RF` and :ref:`OPTIMIZE_.GDIIS`.


keyword(VR-BON)

Only has effect together with .VRML or .VISUAL. Specifies that the VRML files
should include bonds between nearby atoms. The bonds are drawn as grey cylinders,
making it easier to see the structure of the molecule. If .VR-BON is omitted, only the
spheres representing the different atoms will be drawn.


keyword(VR-EIG)

Only has effect together with .VRML or .VISUAL. Specifies that the eigenvectors of
the molecule (that is the eigenvectors of the molecular Hessian, which differs from the
normal modes as they are not mass-scaled) should be visualized. These are written
to the files eigv_###.wrl.


keyword(INITHE)

Specifies that the initial molecular Hessian should be calculated (analytical 
molecular Hessian), thus yielding a first step that is identical to that of second-order 
methods. This provides an excellent starting point for first-order methods, but should only
be used when the molecular Hessian can be calculated within a reasonable amount of
time. It has only effect for first-order methods and overrides the keywords .INITEV
and .INIRED. It has no effect when .HESFIL has been specified.


keyword(INITEV)

::

 READ(LUCMD,*) EVLINI 

The default initial molecular Hessian for first-order minimizations
is the identity matrix when Cartesian coordinates are used, and a diagonal
matrix when redundant internal coordinates are used. If .INITEV is used, all the 
diagonal elements (and therefore the eigenvalues) are set equal to the value EVLINI. This
option only has effect when first-order methods are used and .INITHE and .HESFIL
are non-present.


keyword(HESFIL)

Specifies that the initial molecular Hessian should be read from the file DALTON.HES(?).
This applies to first-order methods, and the Hessian in the file must have the correct
dimensions. This option overrides other options for the initial Hessian.
Each time a Hessian is calculated or updated, it?s written to this file (in Cartesian
coordinates). If an optimization is interrupted, it can be restarted with the last geometry
and the molecular Hessian in DALTON.HES, minimizing the loss of information. 
Another useful possibility is to transfer the molecular Hessian from a calculation on the
same molecule with another (smaller) basis and/or a cheaper wave function. Finally,
one can go in and edit the file directly to set up a specific force field.


keyword(REJINI)

Specifies that the molecular Hessian should be reinitialized after every rejected
step, as a rejected step indicates that the molecular Hessian models the true potential
surface poorly. Only applies to first-order methods.


keyword(STEEPD)

Specifies that the first-order steepest descent method should be used. No update
is done on the molecular Hessian, so the optimization will be guided by the gradient
alone. The ?pure? steepest descent method is obtained when the molecular Hessian
is set equal to the identity matrix. Each step will then be the negative of the gradient
vector, and the convergence towards the minimum will be extremely slow. However,
this option can be combined with other initial molecular Hessians in Cartesian or
redundant internal coordinates, giving a method where the main feature is the lack
of molecular Hessian updates (static molecular Hessian).


keyword(RANKON)

Specifies that a first-order method with the rank one update formula should be
used for optimization. This updating is also referred to as symmetric rank one (SR1)
or Murtagh-Sargent (MS).


keyword(PSB)

Specifies that a first-order method with the Powell-Symmetric-Broyden (PSB) update
formula should be used for optimization.


keyword(DFP)

Specifies that a first-order method with the Davidon-Fletcher-Powell (DFP) update
formula should be used for optimization. May be used for both minimizations and
transition state optimizations.


keyword(BFGS)

Specifies the use of a first-order method with the Broyden-Fletcher-Goldfarb-Shanno
(BFGS) update formula for optimization. This is the preferred first-order method for
minimizations, as this update is able to maintain a positive definite Hessian. Note that
this also makes it unsuitable for transitions state optimization (where one negative
eigenvalue is sought).


keyword(NEWTON)

Specifies that a second-order Newton method should be used for optimization
that is, the analytical molecular Hessian will be calculated at every geometry. By
default the level-shifted trust region method will be used, but it is possible to override
this by using one of the two keywords :ref:`OPTIMIZE_.RF` or :ref:`OPTIMIZE_.GDIIS`.
Not implemented in DIRAC.


keyword(QUADSD)

Use the 2nd order quadratic steepest descent method 
(missing in Dalton manual).


keyword(SCHLEG)

Specifies that a first-order method with Schlegel's updating scheme [dalton ref 152] should
be used. This makes use of all previous displacements and gradients, not just the last,
to update the molecular Hessian.


keyword(HELLMA)

Use gradients and Hessians calculated using the Hellmann-Feynman approximation.
Currently not working properly.


keyword(M-BFGS)

A list of old geometries and gradients are kept. At each new point, displacements
and gradient difference for the last few steps are calculated, and all of these are then
used to sequentially update the molecular Hessian, the most weight being given to
the last displacement and gradient difference. Each update is done using the BFGS
formula, and it?s thus only suitable for minimizations. Only applies to first-order
methods.


keyword(CARTES)

Indicates that Cartesian coordinates should be used in the optimization. This is
the default for second-order methods.


keyword(REDINT)

Specifies that redundant internal coordinates should be used in the optimization.
This is the default for first-order methods.


keyword(INIRED)

Use a simple model Hessian [dalton ref 18] diagonal in redundant internal coordinates as the
initial Hessian. All diagonal elements are determined based on an extremely simplified
molecular mechanics model, yet this model provides Hessians that are good starting
points for most systems, thus avoiding any calculation of the exact Hessian. This is
the default for first-order methods.


keyword(1STORD)

Use default first-order method. This means that the BFGS update will be used,
and that the optimization is carried out in redundant internal coordinates. Same
effect as the combination of the two keywords :ref:`OPTIMIZE_.BFGS` and :ref:`OPTIMIZE_.REDINT`.
Since the BFGS method ensures a positive definite Hessian,
the :ref:`OPTIMIZE_.BOFILL` optimization method is used
by default in case of searches for transition states.


keyword(2NDORD)

Use default second-order method. Molecular Hessians will be calculated at every
geometry. The level-shifted Newton method and Cartesian coordinates are used.
Identical to specifying the keywords :ref:`OPTIMIZE_.NEWTON` and :ref:`OPTIMIZE_.CARTES`.
Not implemented in DIRAC.


keyword(GRDINI)

Specifies that the molecular Hessian should be reinitialized every time the norm of
the gradient is larger than norm of the gradient two iterations earlier. This keyword
should only be used when it?s difficult to obtain a good approximation to the molecular
Hessian during optimization. Only applies to first-order methods.

keyword(DISPLA)

::

 READ (LUCMD,*) DISPLA

Read one more line containing the norm of the displacement vector to be used during
numerical evaluation of the molecular gradient, as is needed when doing geometry
optimizations with CI or MP2 wave functions. Default is 0.001 a.u.


keyword(CONSTR)

::

 READ (LUCMD, *) NCON
 DO I = 1, NCON
 READ(LUCMD,*) ICON
 ICNSTR(ICON) = 1
 END DO

Request a constrained geometry optimization. Only works when using redundant
internal coordinates. The number of primitive coordinates that should be frozen has
to be specified (NCON), then a list follows with the individual coordinate numbers.
The coordinate numbers can be found by first running Dalton(DIRAC?) with the .FINDRE
keyword. Any number of bonds, angles and dihedral angles may be frozen. NOTE:
Symmetry takes precedence over constraints, if you e.g. want to freeze just one of
several symmetric bonds, symmetry must be lowered or switched off.


keyword(MODHES)

Determine a new model molecular Hessian (see .INIMOD) at every geometry without
doing any updating. The model is thus used in much the same manner as an
exact molecular Hessian, though it is obviously only a relatively crude approximation
to the analytical molecular Hessian.


keyword(REMOVE)

::

 READ (LUCMD, *) NREM
 DO I = 1, NREM
 READ(LUCMD,*) IREM
 ICNSTR(IREM) = 2
 END DO

Only has effect when using redundant internal coordinates. Specifies internal 
coordinates that should be removed. The input is identical to the one for .CONSTRAINT,
that is one has to specify the number of coordinates that should be removed, then
the number of each of those internal coordinates. The coordinate numbers can first
be determined by running with .FINDRE set.
Removing certain coordinates can sometimes be useful in speeding up constrained
geometry optimization, as certain coordinates sometimes "struggle" against the 
constraints. See also .NODIHE.


keyword(INIMOD)

Use a simple model Hessian [dalton ref 18] diagonal in redundant internal coordinates as the
initial Hessian. All diagonal elements are determined based on an extremely simplified
molecular mechanics model, yet this model provides Hessians that are good starting
points for most systems, thus avoiding any calculation of the exact Hessian. This is
the default for first-order methods.


keyword(FINDRE)

Determines the redundant internal coordinate system then quits without doing an
actual calculation. Useful for setting up constrained geometry optimizations, where
the numbers of individual primitive internal coordinates are needed.


keyword(CMBMOD)

Uses a combination of the BFGS update and the model molecular Hessian 
(diagonal in redundant internal coordinates). The two have equal weight in the first
iteration of the geometry optimization, then for each subsequent iteration the weight
of the model Hessian is halved. Only suitable for minimizations.


keyword(RF)

Use the rational function method [dalton ref.12] instead of level-shifted Newton which is the
default. The RF method is often slightly faster than the level-shifted Newton, but
also slightly less robust.
For saddle point optimizations there's a special partitioned rational function method
(used automatically when both .RF and .SADDLE are set). However, this method is
both slower and less stable than the default trust-region level-shifted image method
(which is the default).


keyword(GDIIS)

Use the Geometrical DIIS[dalton ref 151] algorithm to control the step. Works in much the
same way as DIIS for wave functions. However, the rational function and level-shifted
Newton methods are generally more robust and more efficient. Can only be used for
minimizations.


keyword(DELINT)

Use delocalized internal coordinates. These are built up as non-redundant linear
combinations of the redundant internal coordinates. Performance is more or less the
same as for the redundant internals, but the transformation of displacements (step)
is slightly less stable.

keyword(NODIHE)

No dihedral angles will be used as coordinates (just bonds and angles).


keyword(VR-COR)

Draws x-, y- and z-axis in the VRML scenes with geometries. Somewhat useful if
one is struggling to build a reasonable geometry by adjusting coordinates manually.


keyword(VR-VIB)

Similar to .VR-EIG, but more useful as it draws the actual normal mode vectors
(the mass-weighted eigenvectors). These are written to the files norm_###.wrl.
Keyword only has effect when a vibrational analysis has been requested.


keyword(VR-SYM)

Draws in all symmetry elements of the molecule as vectors (rotational axes) and
semi-transparent planes (mirror planes).


keyword(M-PSB)

This identical to .M-BFGS, except the PSB formula is used for the updating. Only
applies to first-order methods, but it can be used for both minimizations and saddle
point optimizations.


keyword(LINE S)

Turns on line searching, using a quartic polynomial. By default this is turned off,
as there seems to be no gain in efficiency. Can only be used for minimizations.


keyword(SADDLE)

Indicates that a saddle point optimization should be performed rather than a
minimization. The default method is to calculate the molecular Hessian analytically
at the initial geometry, then update it using Bofill's update. The optimization is
performed in redundant internal coordinates and using the trust-region level-shifted
image method to control the step. That is by default all the keywords .INITHE,
.FILL and .REDINT are already set, but this can of course be overridden by specifying
other keywords. If locating the desired transition state is difficult, and provided
analytical molecular Hessians are available, it may sometimes be necessary to use the
.NEWTON keyword so that molecular Hessians are calculated at every geometry.


keyword(MODE)

::

 READ(LUCMD,*) NSPMOD

Only has effect when doing saddle point optimizations. Determines which molecular
Hessian eigenmode should be maximized (inverted in the image method). By default
this is the mode corresponding to the lowest eigenvalue, i.e. mode 1. If an optimization
does not end up at the correct transition state, it may be worthwhile following other
modes (only the lower ones are usually interesting).


keyword(BOFILL)

Bofill's update [dalton ref 20] is the default updating scheme for transition state optimizations.
It's a linear combination of the symmetric rank one and the PSB updating
schemes, automatically giving more weight to PSB whenever the rank one potentially
is numerically unstable.


keyword(LINDH)

Use original Roland Lindh r_ref with .MODHES
(missing in Dalton manual).


keyword(GRD IN)

Specifies that the molecular Hessian should be reinitialized every time the norm of
the gradient is larger than norm of the gradient two iterations earlier. This keyword
should only be used when it?s difficult to obtain a good approximation to the molecular
Hessian during optimization. Only applies to first-order methods.


keyword(GRD SCREEN)

::

 READ (LUCMD,*) SCRGRD

Read the screening threshold for gradient calculations
(missing in Dalton manual).


keyword(NOAUX)

Only has effect when using redundant internal coordinates. The default for minimizations
is to add auxiliary bonds between atoms that are up to two and half times
further apart then regular (chemical) bonds. This increases the redundancy of the
coordinate system, but usually speeds up the geometry optimization slightly. .NOAUX
turns this off. For saddle point optimizations and constrained geometry optimization
this is off by default (cannot be switched on).


keyword(BFGSR1)

Use a linear combination of the BFGS and the symmetric rank one updating
schemes in the same fashion as Bofill's update. Only suitable for minimizations.


keyword(IPRGRD)

Print level for DIRAC molecular gradient evaluation (missing in Dalton manual).

