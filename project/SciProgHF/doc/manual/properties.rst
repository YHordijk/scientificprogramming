:orphan:
 

starstar(PROPERTIES)

This section allows for the evaluation of a large number of molecular
properties. Available properties include:

-  Expectation values (e.g. dipole moment and electric field gradients).

-  Linear response properties (e.g. polarizability and NMR parameters).

-  Quadratic response properties (e.g. hyperpolarizabilities).

For convenience some common properties can be specified directly in this
section, which means that the user in principle does not need to know
how they are calculated. Note, however, that response functions are by
default static, but frequencies can be added in the relevant
subsections.

Properties which are not predefined must be specified in detail in the
relevant input section (see :ref:`one_electron_operators`).

By default no properties are calculated.

General control statements
==========================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(ABUNDANCIES)

For properties that make reference to isotopes, give threshold level (in
% abundance) for isotopes to print.

*Default:*

::

    .ABUNDANCIES
     1.0

keyword(RKBIMP)

Import coefficients calculated with restricted kinetic balance (RKB) in
a calculation using unrestricted kinetic balance (UKB). This option is a
simple way to generated restricted magnetic balance for the calculation
of NMR shieldings. This option works in the general SO case, but not in
the spinfree case since spinfree calculations are not possible with UKB.

keyword(NOPCTR)

In two-component infinite-order relativistic calculations (with :ref:`HAMILTONIAN_.X2C`) take only LL block
of four-component property operators to avoid the picture change transformation.
Experimental option, use with care.

keyword(RDCCDENS)

Activates the reading of the file CCDENS obtained from a previous CC calculation with either the 
:ref:`WAVE_FUNCTION_.RELCCSD` or with :ref:`WAVE_FUNCTION_.EXACC` modules.  It is not necessary to 
use this keyword in runs in which the correlated and  property modules are both activated.

CCDENS is not saved by pam, so unless the scratch directory from the previous calculation is
kept (see --keep_scratch in pam), you should retrieve it after a correlated calculation, e.g. ::

   pam --get=CCDENS ...

and then copy it back for the property calculation e.g. ::

   pam --put=CCDENS ...



Predefined electric properties
==============================

keyword(DIPOLE)

Evaluate the electronic electric dipole moment

Expectation values: :math:`\langle\hat{\mu}_\alpha\rangle=-e\langle r_\alpha\rangle`

keyword(QUADRUPOLE)

Evaluate the electronic traceless electric quadrupole moment

Expectation values: :math:`\langle\Theta_{\alpha\beta}\rangle=-e\frac {3}{2}\langle r_\alpha r_\beta-\frac{1}{3}\delta_{\alpha\beta}r^2\rangle`)

keyword(EFG)

Evaluate electric field gradients at nuclear positions, see :cite:`Visscher_JCP1998` .

Electronic contribution to center :math:`K` (expectation values) :

.. math::

    \phi^{[2]el}_{\alpha\beta}(\mathbf{R}_K)=\frac{-e}{4\pi\varepsilon_0}\left<\frac{3r_{K;\alpha}r_{K;\beta}-\delta_{\alpha\beta}r_K^2}{r_K^5}\right>

Nuclear contributions to center :math:`K` :

.. math::

   \phi^{[2]nuc}_{\alpha\beta}(\mathbf{R}_K)=\sum_{A\ne K}\frac{Z_Ae}{4\pi\varepsilon_0}\left[\frac{3r_{KA;\alpha}r_{KA;\beta}-\delta_{\alpha\beta}r_{KA}^2}{r_{KA}^5}\right]

Results are also reported with respect to a principal axis system for each center.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.


keyword(NQCC)

Evaluate nuclear quadrupole coupling constants (NQCC) (expectation values).
The NQCC is formally defined as

.. math::

   \frac{e^2qQ}{h}

where :math:`Q` is the electric quadrupole moment of the nucleus and :math:`q=\phi^{[2]}_{zz}/e` (in the principal axis system) is the field gradient. 
The NQCC may be extracted from experiment, whereas electronic structure calculations may provide the field gradient :math:`q`.
The two quantities are related as

.. math::

   \mbox{NQCC [in MHz] } = 234.9647\ \times\ Q\mbox{ [in b] }\ \times\ q\mbox{ [in atomic units }E_h/ea_0^2\mbox{ ]}

The calculations proceed similar to :ref:`PROPERTIES_.EFG`. The total electric field gradients
for each center are transformed to a principal axis system for which

.. math::

   |\phi^{[2]}_{zz}|\ge|\phi^{[2]}_{yy}|\ge|\phi^{[2]}_{xx}|

DIRAC reports the more general expressions

.. math::

   \mbox{NQCC}_{\alpha\alpha}\mbox{ [in MHz] } = 234.9647\ \times\ Q\mbox{ [in b] }\ \times\ \phi^{[2]}_{\alpha\alpha}/e \mbox{ [in atomic units }E_h/ea_0^2\mbox{ ]}
   
The asymmetry factor is defined as

.. math::

   \eta = \frac{\phi^{[2]}_{xx}-\phi^{[2]}_{yy}}{\phi^{[2]}_{zz}}


Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(POLARIZABILITY)

Evaluate the electronic dipole polarizability tensor, see :cite:`Saue2003` (HF) and :cite:`Salek2005` (DFT).

Linear response function: :math:`\quad\alpha_{\alpha\beta}(-\omega;\omega)=\langle\langle\hat{\mu}_{\alpha};\hat{\mu}_{\beta}\rangle\rangle_{\omega}`

keyword(FIRST ORDER HYPERPOLARIZABILITY)

Evaluate static electronic dipole first-order hyperpolarizability tensor, see :cite:`Norman_JCP2004` (HF) and :cite:`Henriksson:2008` (DFT).

Quadratic response function: :math:`\quad\beta_{\alpha\beta\gamma}(-\omega_\sigma;\omega_1,\omega_2)=\langle\langle\hat{\mu}_{\alpha};\hat{\mu}_{\beta},\hat{\mu}_{\gamma}\rangle\rangle_{\omega_1,\omega_2}`


Results are also given for the static electronic
dipole polarizability.

keyword(TWO-PHOTON)

Evaluate two-photon absorption cross sections :cite:`Henriksson:2005`, obtained as a first-order residue of the first-order hyperpolarizability. Give
the number of desired states in each boson symmetry. Cannot be specified
in combination with other quadratic response calculations.

*Example:* Point group with four boson irreps, (e.g. :math:`C_{2v}`)

::

    .TWO-PHOTON
     5 5 5 0

Predefined magnetic properties
==============================

keyword(NMR)

Evaluate nuclear magnetic shieldings and indirect spin-spin couplings (linear response functions), see :cite:`Visscher_jcc1999` and :cite:`Ilias2009`.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

See below for advice on calculation of diamagnetic terms.

keyword(SHIELDING)

Evaluate nuclear magnetic shieldings (linear response), see :cite:`Visscher_jcc1999` and :cite:`Ilias2009` .

Elements of the shielding tensor for center :math:`K` are given by

.. math::

   \sigma_{K;\mu\nu}=\frac{\partial^2}{\partial m_{K;\mu}\partial B_{0;\nu}}\langle\langle\hat{h}^{hfs}_{K};\hat{h}^Z\rangle\rangle_0

where appears the relativistic hyperfine operator

.. math::

   \hat{h}^{hfs}_{K}=-\sum_i\mathbf{m}_K\cdot\hat{\mathbf{B}}^{el}_{K}(i);\quad \mathbf{B}^{el}_{K}(i)=-\frac{1}{4\pi\varepsilon_0 c^2}\frac{\mathbf{r}_{iK}\times ec\boldsymbol{\boldsymbol{\alpha}}}{r_{iK}^3},

expressed in terms of the nuclear magnetic dipole :math:`\mathbf{m}_K` and the operator :math:`\hat{\mathbf{B}}^{el}_{K}` giving the magnetic field due to the electrons at the nuclear position, and the relativistic Zeeman operator

.. math::

   \hat{h}^Z=-\hat{\mathbf{m}}_e^{[1]}\cdot\mathbf{B}_0;
   \quad\hat{\mathbf{m}}^{[1]}_e=-\sum_i\frac{e}{2}(\mathbf{r}_{iG}\times c\boldsymbol{\alpha}(i)),

expressed in terms of the operator :math:`\hat{\mathbf{m}}_e^{[1]}` associated with the magnetic dipole moment of the electrons and the external magnetic field :math:`\mathbf{B}_0`. Note reference to the gauge origin :math:`G`.

Note that :ref:`PROPERTIES_.PRINT` 2 gives the full tensor and longer output. The :ref:`PROPERTIES_.PRINT` 4 gives the raw values in symmetry coordinates as well.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(MAGNET)

Evaluate the (static) magnetizablity tensor :cite:`Ilias2013`

Lineare response function

.. math::

      \xi_{K;\alpha\beta}= - \frac{\partial^2}{\partial B_{\alpha}\partial B_{\beta}}\langle\langle\hat{h}^Z;\hat{h}^Z\rangle\rangle_0

where appears the relativistic Zeeman operator

.. math::

   \hat{h}^Z=-\hat{\mathbf{m}}_e^{[1]}\cdot\mathbf{B}_0;
   \quad\hat{\mathbf{m}}^{[1]}_e=-\sum_i\frac{e}{2}(\mathbf{r}_{iG}\times c\boldsymbol{\alpha}(i)),

expressed in terms of the operator :math:`\hat{\mathbf{m}}_e^{[1]}` associated with the magnetic dipole moment of the electrons and the external magnetic field :math:`\mathbf{B}_0`. Note reference to the gauge origin :math:`G`.

keyword(ROTG)

Evaluate rotational g-tensors: linear response and nuclear contributions, see :cite:`Aucar_JCP2014`.

Elements of the rotational g-tensor are given by

.. math::

   g_{\mu\nu} = g^{nuc}_{\mu\nu} + g^{elec}_{\mu\nu}

with

.. math::

   g^{elec}_{\mu\nu}= - \frac{2 m_p}{e} \frac{\partial^2}{\partial L_{\mu}\partial B_{0;\nu}}\langle\langle\hat{h}^{BO};\hat{h}^Z\rangle\rangle_0

where appears the first order correction to the Born-Oppenheimer (BO) approximation

.. math::

   \hat{h}^{BO}=-\boldsymbol{\boldsymbol{\omega}}\cdot\hat{\mathbf{J}}_e;
   \quad\boldsymbol{\boldsymbol{\omega}}=\mathbf{L}\otimes\mathbf{I}^{-1},

expressed in terms of the angular velocity :math:`\boldsymbol{\boldsymbol{\omega}}` associated with the total angular momentum of the electrons :math:`\hat{\mathbf{J}}_e=\hat{\mathbf{L}}_e+\hat{\mathbf{S}}_e`, and the relativistic Zeeman operator

.. math::

   \hat{h}^Z=-\hat{\mathbf{m}}_e^{[1]}\cdot\mathbf{B}_0;
   \quad\hat{\mathbf{m}}^{[1]}_e=-\sum_i\frac{e}{2}(\mathbf{r}_{iG}\times c\boldsymbol{\alpha}(i)),

expressed in terms of the operator :math:`\hat{\mathbf{m}}_e^{[1]}` associated with the magnetic dipole moment of the electrons and the external magnetic field :math:`\mathbf{B}_0`. Note reference to the gauge origin :math:`G`.

Results are dimensionless.

The total g-tensor, as well as its linear response and nuclear contributions are always given separately.

Using :ref:`PROPERTIES_.PRINT` 1, the paramagnetic (e-e) and diamagnetic (e-p) parts of the linear response contributions are given separately, together with results for the :math:`\mathbf{L}` and :math:`\mathbf{S}` parts of the linear response.

keyword(SPIN-SPIN COUPLING)

Evaluate indirect spin-spin couplings :cite:`Visscher_jcc1999` . The indirect spin-spin tensor :math:`J_{KL}` associated with nuclei :math:`K` and :math:`L` may be expressed as

.. math::

   2\pi J_{KL}=\gamma_{K}\gamma_{L}K_{KL}

where appears gyromagnetic ratios :math:`\gamma_K`. The elements of the reduced tensor :math:`K_{KL}` are expressed in terms of linear response functions as

.. math::

   K_{KL:\mu\nu} = \frac{\partial^{2}}{\partial m_{K;\mu}\partial m_{L;\nu}}\langle \langle \hat{h}_{K,\mathrm{rel}}^{\rm hfs}; \hat{h}_{L,\mathrm{rel}}^{\rm hfs}\rangle\rangle _{0}

where appears the relativistic hyperfine operator

.. math::

   \hat{h}^{hfs}_{K}=-\sum_i\mathbf{m}_K\cdot\hat{\mathbf{B}}^{el}_{K}(i);\quad \mathbf{B}^{el}_{K}(i)=-\frac{1}{4\pi\varepsilon_0 c^2}\frac{\mathbf{r}_{iK}\times ec\boldsymbol{\boldsymbol{\alpha}}}{r_{iK}^3},

expressed in terms of the nuclear magnetic dipole :math:`\mathbf{m}_K` and the operator :math:`\hat{\mathbf{B}}^{el}_{K}` giving the magnetic field due to the electrons at the nuclear position.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

The default is to calculate the diamagnetic term via occupied positive energy to virtual negative energy orbital rotations
(also called electron-positron rotations), see :cite:`Aucar1999` for the theory.
The quality of this is very basis set dependent.
It is generally more accurate to use the non-relativistic expectation value expression
for the diamagnetic term, activated with keyword .DSO in this section.
You must also add :ref:`LINEAR_RESPONSE_.SKIPEP` under :ref:`*LINEAR RESPONSE` to 
exclude the diamagnetic term from the linear response calculation.

keyword(DSO)

Evaluate the diamagnetic contribution to indirect spin-spin couplings as
an expectation value of the non-relativistic DSO operator.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(NSTDIAMAGNETIC)

Evaluate the diamagnetic contribution to nuclear magnetic shielding tensor as
an expectation value of the non-relativistic operator.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

Mixed electric and magnetic properties
======================================

keyword(OPTROT)

Calculate optical rotation. The most common experimental setup uses
light with a frequency corresponding to the sodium D-line (589.29 nm).
The optical rotation is reported as the number of 
degrees of rotation of the plane of polarization per mole of sample
for a sample cell of length 1~dm, and at a temperature of 
:math:`25^\circ{\rm C}`

.. math::

   \left[\alpha\right]_D^{25} = -288\cdot 10^{-30}\frac{\pi^2\mathcal{N}a_0^4\omega}{3M}\sum_{\alpha}G^{\prime}_{\alpha\alpha };\quad G^{\prime}_{\alpha\beta}(-\omega;\omega)= -\mbox{Im}\langle \langle \hat{\mu}_{\alpha};\hat{m}_{\beta}\rangle \rangle _\omega

where :math:`M` is the molecular mass in g :math:`\mbox{mol}^{-1}` and :math:`\mathcal{N}` is the number density. 


keyword(VERDET)

Evaluate Verdet constants :cite:`Ekstrom2005` for a dynamic electric field corresponding to Ruby
laser wavelength of 694 nm and a static magnetic field along the
propagation direction of the light beam (in this case, the default
frequencies of the quadratic response function thus become
ω\ :sub:`*B*`\  = 0.0656 and ω\ :sub:`*C*`\  = 0.0).

The Verdet constant is given in terms of quadratic response functions

.. math::

   V(\omega)=\omega C\epsilon_{\alpha\beta\gamma}\mbox{Im}\langle\langle\hat{\mu}_{\alpha};\hat{\mu}_{\beta},\hat{m}_{\gamma}\rangle\rangle_{\omega,0}

where :math:`C=eN/(24c_0\epsilon_0m_e` and :math:`N` is the number density of the gas.  A Verdet calculation cannot be specified in combination with other quadratic
response calculations.

The frequencies can be changed using :ref:`QUADRATIC_RESPONSE_.B FREQ` in :ref:`*QUADRATIC RESPONSE`.


Other predefined properties
===========================

keyword(MOLGRD)

Evaluate the molecular gradient, i.e.

.. math::

  \frac{\partial E}{\partial \mathbf{X}_A}

where :math:`\mathbf{X}_{A}` are the coordinates of the nuclei. This
is an expectation value of one- and two-electron operators. Normally the
molecular gradient evaluation is not invoked explicitly with this
keyword but rather implicitly in the geometry optimization module.

keyword(PVC)

Calculate matrix elements over the nuclear spin-independent
parity-violating operator, e.g. calculate energy differences between
enantiomers, see :cite:`Laerdahl1999` and :cite:`Bast2011` .

The parity violating energy is calculated as an expectation value

.. math::

   E_\text{PV} = \sum_A \langle H_\text{PV}^A \rangle;\quad H_\text{PV}^A= \frac{G_\text{F}}{2 \sqrt{2}} Q_\text{w}^A \sum_i \gamma_5 (i) \rho^A (\mathbf{r}_i)

where appears appears the normalized nuclear charge densities :math:`\rho^A` and the the Fermi coupling constant :math:`G_F=2.222 55\times 10^{-14}E_ha_0^3`.  The weak nuclear charge

.. math::

   Q_\text{w}^A = Z^AC_V^\text{p}+N^AC_V^\text{n}=Z^A (1 - 4 \sin^2 \theta_\text{W}) - N^A,

given in terms of the number of protons and neutrons --- :math:`Z^A` and :math:`N^A` --  in nucleus :math:`A` and the Weinberg angle :math:`\theta_\text{W}` which describes the rotation of :math:`B^0` and :math:`W^0` bosons by spontaneous symmetry breaking to form photons and :math:`Z^0` bosons (DIRAC uses :math:`\sin^2 \theta_\text{W} = 0.2319`). 

keyword(PVCNMR)

Calculate parity-violating contribution to the NMR shielding tensor, see :cite:`Bast:2006`.
Elements of the parity-violating contribution to the shielding tensor for center :math:`K` are
given by

.. math::

   \sigma^{PV}_{K;\mu\nu}=\frac{\partial^2}{\partial m_{K;\mu}\partial B_{0;\nu}}\langle\langle h_\text{PV2}^K;\hat{h}^Z\rangle\rangle_0

where appears the nuclear spin-dependent parity-violating operator

.. math::

   h_\text{PV2}^K = -\frac{G_\text{F}(1-4\sin^2\theta_\text{W})}{\sqrt{2}}\sum_{i}\frac{1}{\gamma_K}\boldsymbol{\alpha}\cdot\mathbf{M}_K\rho_K(\mathbf{r}_i)

where appears the Fermi coupling constant :math:`G_F=2.222 55\times 10^{-14}E_ha_0^3`, the gyromagnetic ratio :math:`\gamma_K`, the nuclear magnetic dipole moment :math:`\mathbf{M}_K=\gamma_K\mathbf{I}_K` and the normalized nuclear charge density :math:`\rho_K`.

The relativistic Zeeman operator

.. math::

   \hat{h}^Z=-\hat{\mathbf{m}}_e^{[1]}\cdot\mathbf{B}_0;
   \quad\hat{\mathbf{m}}^{[1]}_e=-\sum_i\frac{e}{2}(\mathbf{r}_{iG}\times c\boldsymbol{\alpha}(i)),

is expressed in terms of the operator :math:`\hat{\mathbf{m}}_e^{[1]}` associated with the magnetic dipole moment of the electrons and the external magnetic field :math:`\mathbf{B}_0`. Note reference to the gauge origin :math:`G`.

keyword(RHONUC)

Calculate electronic density at the nuclear positions, also known as the contact density :cite:`Knecht2011` and :cite:`Almoukhalalati:2016b`.

It is formally an expectation value

.. math::

   \rho_e^K = -e\langle\delta^3(\mathbf{r}-\mathbf{R}_K)\rangle,

an important observation in view of picture change effects, see  :cite:`Knecht2011` .

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.

keyword(EFFDEN)

Calculate effective electronic density associated with nuclei, see :cite:`Knecht2011`.
This quantity appears in expressions for the Mössbauer isomer shift. Starting from
the electrostatic electron-nucleus interaction

.. math::

   E^{el}\left(R\right)=\int  \rho_e (\mathbf{r})\phi_n(\mathbf{r};R){\rm d}^3 \mathbf{r} ,

we consider the change in the electrostatic energy upon a change of nuclear radius. If we
ignore any change in the electronic density :math:`\rho_e`, we may express this as

.. math::

   \Delta E_{\gamma} = \left.\frac{\partial E^{el}}{\partial R}\right|_{R=R_0}\Delta R = \left[\int \rho_e (\mathbf{r})\frac{\partial\phi_n(\mathbf{r})}{\partial R}{\rm d}^3 \mathbf{r}\right]_{R=R_0} \Delta R = \bar\rho_e\int \left[\frac{\partial\phi_n(\mathbf{r})}{\partial R}{\rm d}^3 \mathbf{r}\right]_{R=R_0} \Delta R

where the effective density :math:`\bar\rho_e` is introduced in the last step. It is often approximated by the contact density, see :ref:`PROPERTIES_.RHONUC` , but this is discouraged since it may introduce errors on the order of 10%.


Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.



keyword(SPIN-ROTATION)

Evaluate nuclear spin-rotation constants: linear response, expectation value and nuclear contributions, see :cite:`Aucar_JCP2012` and :cite:`AucarChap2019`.

Elements of the nuclear spin-rotation tensor for center :math:`K` in a molecule in equilibrium are given by

.. math::

   M_{K;\mu\nu} = M^{nuc}_{K;\mu\nu} + M^{elec}_{K;\mu\nu}

with

.. math::

   M^{elec}_{K;\mu\nu} = - \hslash^2 \frac{\partial^2}{\partial I_{K;\mu}\partial L_{\nu}}\langle\langle\hat{h}^{hfs}_{K};\hat{h}^{BO}\rangle\rangle_0

where appears the relativistic hyperfine operator

.. math::

   \hat{h}^{hfs}_{K}=-\sum_i\mathbf{m}_K\cdot\hat{\mathbf{B}}^{el}_{K}(i);\quad \mathbf{B}^{el}_{K}(i)=-\frac{1}{4\pi\varepsilon_0 c^2}\frac{\mathbf{r}_{iK}\times ec\boldsymbol{\boldsymbol{\alpha}}}{r_{iK}^3},

expressed in terms of the nuclear magnetic dipole :math:`\mathbf{m}_K=\hslash \gamma_K \mathbf{I}_K` and the operator :math:`\hat{\mathbf{B}}^{el}_{K}` giving the magnetic field due to the electrons at the nuclear position, and the first order correction to the Born-Oppenheimer (BO) approximation

.. math::

   \hat{h}^{BO}=-\boldsymbol{\boldsymbol{\omega}}\cdot\hat{\mathbf{J}}_e;
   \quad\boldsymbol{\boldsymbol{\omega}}=\mathbf{L}\otimes\mathbf{I}^{-1},

expressed in terms of the angular velocity :math:`\boldsymbol{\boldsymbol{\omega}}` associated with the total angular momentum of the electrons :math:`\hat{\mathbf{J}}_e=\hat{\mathbf{L}}_e+\hat{\mathbf{S}}_e`. Note that the origin of the orbital angular momentum is the molecular center of mass.


.. note::

   The current implementation gives by default (:ref:`PROPERTIES_.PRINT` values up to 3) results only for molecules in equilibrium.

   
Results are given in kHz, but for particular :ref:`PROPERTIES_.PRINT` values they could also be given in ppm (to compare results with :ref:`PROPERTIES_.SHIELDING`).

The total spin-rotation tensors, as well as their electronic (linear response) and nuclear contributions are always given separately.

Using :ref:`PROPERTIES_.PRINT` 0 or 1, results are only given in kHz, whereas employing :ref:`PROPERTIES_.PRINT` 2 or 3 they are also shown in ppm.

In addition, when :ref:`PROPERTIES_.PRINT` 1 or 3 are used, the paramagnetic-like (e-e) and diamagnetic-like (e-p) parts of the linear response contributions are given separately, together with results for the :math:`\mathbf{L}` and :math:`\mathbf{S}` parts of the linear response.

Finally, employing :ref:`PROPERTIES_.PRINT` 4 expectation value and nuclear contributions are given separately, with the inclusion of Thomas precesion effects, in order to properly include contributions to nuclear spin-rotations out of the equilibrium geometry of the molecular system.

For more details, the user is welcome to see the Tutorial Section.

Atomic centers may be restricted with :ref:`INTEGRALS_.SELECT` under :ref:`**INTEGRALS`.
