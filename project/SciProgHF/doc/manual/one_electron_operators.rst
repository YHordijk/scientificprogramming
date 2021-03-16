:orphan:
 

.. _one_electron_operators:

One-electron operators
======================

Syntax for the specification of one-electron operators
------------------------------------------------------

A general one-electron property operator in 4-component calculations is generated from
linear combinations of the basic form:

.. math::

  \hat{P} = f M_{4 \times 4} \hat{\Omega}

with the scalar factor :math:`f` and the scalar operator :math:`\hat{\Omega}`, and where

.. math::

  M_{4 \times 4}

is one of the following :math:`4 \times 4` matrices:

.. math::

  I_{4 \times 4}, \gamma_5, \beta \gamma_5,

  i\alpha_x, i\alpha_y, i\alpha_z

  \Sigma_x, \Sigma_y,  \Sigma_z

  \beta \Sigma_x, \beta \Sigma_y, \beta \Sigma_z

  i \beta \alpha_x, i \beta \alpha_y, i \beta \alpha_z

One thing to notice is that an imaginary :math:`i` is added to the time-antisymmetric Dirac :math:`\boldsymbol{\alpha}` - matrices and their derivatives to make them time symmetric and hence fit into the
quaternion symmetry scheme of DIRAC (see :cite:`Saue1999` and :cite:`Salek2005` for more information).

Operator types
--------------

There are 21 basic operator types used in DIRAC, listed in this Table:

===========  ============================================================================= ===============
**Keyword**  **Operator form**                                                             **Nr. factors**
===========  ============================================================================= ===============
DIAGONAL     :math:`f I_{4 \times 4} \Omega`                                               1
XALPHA       :math:`f \alpha_x \Omega`                                                     1
YALPHA       :math:`f \alpha_y \Omega`                                                     1
ZALPHA       :math:`f \alpha_z \Omega`                                                     1
XAVECTOR     :math:`f_1 \alpha_y \Omega_z - f_2 \alpha_z \Omega_y`                         2
YAVECTOR     :math:`f_1 \alpha_z \Omega_x - f_2 \alpha_x \Omega_z`                         2
ZAVECTOR     :math:`f_1 \alpha_x \Omega_y - f_2 \alpha_y \Omega_x`                         2
ALPHADOT     :math:`f_1 \alpha_x \Omega_x + f_2 \alpha_y \Omega_y + f_3 \alpha_z \Omega_z` 3
GAMMA5       :math:`f \gamma_5 \Omega`                                                     1
XSIGMA       :math:`f \Sigma_x \Omega`                                                     1
YSIGMA       :math:`f \Sigma_y \Omega`                                                     1
ZSIGMA       :math:`f \Sigma_z \Omega`                                                     1
XBETASIG     :math:`f \beta \Sigma_x \Omega`                                               1
YBETASIG     :math:`f \beta \Sigma_y \Omega`                                               1
ZBETASIG     :math:`f \beta \Sigma_z \Omega`                                               1
XiBETAAL     :math:`f i \beta \alpha_x \Omega`                                             1
YiBETAAL     :math:`f i \beta \alpha_y \Omega`                                             1
ZiBETAAL     :math:`f i \beta \alpha_z \Omega`                                             1
BETA         :math:`f \beta \Omega`                                                        1
SIGMADOT     :math:`f_1 \Sigma_x \Omega_x + f_2 \Sigma_y \Omega_y + f_3 \Sigma_z \Omega_z` 1
iBETAGAMMA5  :math:`f i \beta \gamma_5 \Omega`                                             1
===========  ============================================================================= ===============


Operator specification
----------------------

Operators are specified by the keyword :ref:`HAMILTONIAN_.OPERATOR` with the following
arguments::

  .OPERATOR
   'operator name'
   operator type keyword
   operator labels for each component
   FACTORS
   factors for each component
   CMULT
   COMFACTOR
   common factor for all components

Note that the arguments following the keyword :ref:`HAMILTONIAN_.OPERATOR` must start with
a blank. The arguments are optional, except for the operator label.
Component factors as well as the common factor are all one if not specified.


List of one-electron operators
------------------------------

+-------------+----------------------------------+----------------+----------------+----------------------------------+
| **Operator**| **Description**                  | **Symmetry**   | **Components** | **Operators**                    |
| **label**   |                                  |                |                |                                  |
+=============+==================================+================+================+==================================+
| MOLFIELD    | Nuclear attraction integrals     | Symmetric      | MOLFIELD       | :math:`\Omega_1 = \sum_K V_{iK}` |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| OVERLAP     | Overlap integrals                | Symmetric      | OVERLAP        | :math:`\Omega_1 = 1`             |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| BETAMAT     | Overlap integrals, only SS-block | Symmetric      | BETAMAT        | :math:`\Omega_1 = 1`             |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| DIPLEN      | Dipole length integrals          | Symmetric      | XDIPLEN        | :math:`\Omega_1 = x`             |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | YDIPLEN        | :math:`\Omega_2 = y`             |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | ZDIPLEN        | :math:`\Omega_3 = z`             |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| DIPVEL      | Dipole velocity integrals        | Anti-symmetric | XDIPVEL        |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | YDIPVEL        |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | ZDIPVEL        |                                  |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| QUADRUP     | Quadrupole moments integrals     | Symmetric      | XXQUADRU       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | XYQUADRU       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | XZQUADRU       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | YYQUADRU       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | YZQUADRU       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | ZZQUADRU       |                                  |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| SPNORB      | Spatial spin-orbit integrals     | Anti-symmetric | X1SPNORB       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | Y1SPNORB       |                                  |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | Z1SPNORB       |                                  |
+-------------+----------------------------------+----------------+----------------+----------------------------------+
| SECMOM      | Second moments integrals         | Symmetric      | XXSECMOM       | :math:`\Omega_1 = xx`            |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | XYSECMOM       | :math:`\Omega_2 = xy`            |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | XZSECMOM       | :math:`\Omega_3 = xz`            |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | YYSECMOM       | :math:`\Omega_4 = yy`            |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | YZSECMOM       | :math:`\Omega_5 = yz`            |
|             |                                  |                +----------------+----------------------------------+
|             |                                  |                | ZZSECMOM       | :math:`\Omega_6 = zz`            |
+-------------+----------------------------------+----------------+----------------+----------------------------------+


=========== ====================================================================================================================
**Keyword** **Description**
=========== ====================================================================================================================
THETA       Traceless theta quadrupole integrals
CARMOM      Cartesian moments integrals, symmetric integrals, (l + 1)(l + 2)/2 components ( l = i + j + k)
SPHMOM      Spherical moments integrals (real combinations), symmetric integrals, (2l + 1) components ( m = +0, -1, +1, ..., +l)
SOLVENT     Electronic solvent integrals
FERMI C     One-electron Fermi contact integrals
PSO         Paramagnetic spin-orbit integrals
SPIN-DI     Spin-dipole integrals
DSO         Diamagnetic spin-orbit integrals
SDFC        Spin-dipole + Fermi contact integrals
HDO         Half-derivative overlap integrals
S1MAG       Second order contribution from overlap matrix to magnetic properties
ANGLON      Angular momentum around the nuclei
ANGMOM      Electronic angular momentum around the molecular center of mass
LONMOM      London orbital contribution to angular momentum
MAGMOM      One-electron contributions to magnetic moment
KINENER     Electronic kinetic energy
DSUSNOL     Diamagnetic susceptibility without London contribution
DSUSLH      Angular London orbital contribution to diamagnetic susceptibility
DIASUS      Angular London orbital contribution to diamagnetic susceptibility
NUCSNLO     Nuclear shielding integrals without London orbital contribution
NUCSLO      London orbital contribution to nuclear shielding tensor integrals
NUCSHI      Nuclear shielding tensor integrals
NEFIELD     Electric field at the individual nuclei
ELFGRDC     Electric field gradient at the individual nuclei, cartesian
ELFGRDS     Electric field gradient at the individual nuclei, spherical
S1MAGL      Bra-differentiation of overlap matrix with respect to magnetic field
S1MAGR      Ket-differentiation of overlap matrix with respect to magnetic field
HDOBR       Ket-differentiation of HDO-integrals with respect to magnetic field
NUCPOT      Potential energy of the interaction of electrons with individual nuclei, divided by the nuclear charge
HBDO        Half B-differentiated overlap matrix
SQHDO       Half-derivative overlap integrals not to be antisymmetrized
DSUSCGO     Diamagnetic susceptibility with common gauge origin
NSTCGO      Nuclear shielding integrals with common gauge origin
EXPIKR      Cosine and sine integrals
MASSVEL     Mass velocity integrals
DARWIN      Darwin type integrals
CM1         First order magnetic field derivatives of electric field
CM2         Second order magnetic field derivatives of electric field
SQHDOR      Half-derivative overlap integrals not to be anti-symmetrized
SQOVLAP     Second order derivatives overlap integrals
=========== ====================================================================================================================


Examples of using various operators
-----------------------------------

We give here several concrete examples on how to construct operators for
various properties.

Kinetic part of the Dirac Hamiltonian
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The kinetic part of the Dirac Hamiltonian may be specified by::

  .OPERATOR
   'Kin energy'
   ALPHADOT
   XDIPVEL
   YDIPVEL
   ZDIPVEL
   COMFACTOR
   -68.51799475

where -68.51799475 is :math:`-c/2`. 

The speed of light :math:`c` is an important parameter in relativistic
theory, but its explicit value in atomic units not necessarily remembered.
A simpler way to specify the kinetic energy operator is therefore::


  .OPERATOR
   'Kin energy'
   ALPHADOT
   XDIPVEL
   YDIPVEL
   ZDIPVEL
   CMULT
   COMFACTOR
   -0.5

where the keyword *CMULT* assures multiplication of the common factor -0.5 by :math:`c`.
This option has the further advantage that *CMULT* follows any user-specified modification 
of the speed of light, as provided by :ref:`GENERAL_.CVALUE`.

XAVECTOR
~~~~~~~~

Another example::

  .OPERATOR
   'B_x'
   XAVECTOR
   ZDIPLEN
   YDIPLEN
   CMULT
   COMFACTOR
   -0.5


The program will assume all operators to be Hermitian and will therefore insert
an imaginary phase *i* if necessary (applies to antisymmetric scalar
operators).

If no other arguments are given, the program assumes the operator to be a
diagonal operator and expects the operator name to be the component label, for
instance::

  .OPERATOR
   OVERLAP


Dipole moment as finite field perturbation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Another example is the finite perturbation calculation with the :math:`\hat{z}`
dipole length operator added to the Hamiltonian (don't forget to decrease the
symmetry of your system):

.. math::

    \hat{H} = \hat{H}_0 + 0.01 \cdot \hat{z}

::

  .OPERATOR
   ZDIPLEN
   COMFACTOR
   0.01


Fermi-contact integrals
~~~~~~~~~~~~~~~~~~~~~~~

Here is an example where the Fermi-contact (FC) integrals for a certain nucleus
are added to the Hamiltonian in a finite-field calculation.  Let's assume you
are looking at a PbX dimer (order in the .mol file: 1.  Pb, 2. X) and you want
to add to the Dirac-Coulomb :ref:`**HAMILTONIAN` the FC integrals for the Pb
nucleus as a perturbation with a given field-strength (FACTORS).

**Important note:** The raw density values obtained after the fit of
your finite-field energies need to be scaled by
:math:`\frac{1}{\frac{4 \pi g_{e}}{3}} = \frac{1}{8.3872954891254192}`,
a factor that originates from the definition of the operator for
calculating the density at the nucleus::

  **HAMILTONIAN
  .OPERATOR
   'Density at nucleus'
   DIAGONAL
   'FC Pb 01'
   FACTORS
   -0.000000001

Here is next example of how-to calculate the electron density at the
nucleus as an expectation value :math:`\langle 0 \vert \delta(r-R) \vert 0
\rangle` for a Dirac-Coulomb HF wave function including a decomposition of the
molecular orbital contributions to the density::

  **DIRAC
  .WAVE FUNCTION
  .PROPERTIES
  **HAMILTONIAN
  **WAVE FUNCTION
  .SCF
  **PROPERTIES
  .RHONUC
  *EXPECTATION
  .ORBANA
  *END OF


Cartesian moment expectation value
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the following example I am calculating a cartesian moment expectation value
:math:`\langle 0 \vert x^1 y^2 z^3 \vert 0 \rangle` for a Levy-Leblond HF wave
function::

  **DIRAC
  .WAVE FUNCTION
  .PROPERTIES
  **HAMILTONIAN
  .LEVY-LEBLOND
  **WAVE FUNCTION
  .SCF
  **PROPERTIES
  *EXPECTATION
  .OPERATOR
   CM010203
  *END OF
