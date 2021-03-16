:orphan:

star(EXCITATION ENERGIES)

Calculate excitation energies using time dependent Hartree-Fock or DFT.
The excitation energies are found as the lowest generalized eigenvalues
of the electronic Hessian. DIRAC supports TDDFT kernels from all ground
state functionals included in the code. Currently the iterative
eigenvalue solver may fail to converge more than about twenty roots per
symmetry.

Define excitations and transition moments
=========================================

keyword(EXCITA)

::

    .EXCITA
    SYM N

Number of excitation energies N calculated in boson symmetry no. SYM.
This keyword can be repeated if you want excitation energies in more
than one boson symmetry.

keyword(OPERATOR)

Specification of a transition moment operator (see
:ref:`one_electron_operators` for details). This keyword can be given multiple
times to add more operators.

keyword(EPOLE)

Specification of electric Cartesian multipole operators of order :math:`n`

.. math::

    \hat{Q}_{j_{1}\ldots j_{n}}^{\left[n\right]}=-er_{1}r_{2}\ldots r_{j_{n}}   

for the calculation of transition moments (note that they contribute to one order less in the wave vector). Specify order.

*Example:* Electric dipole operators::
  
      .EPOLE
      1

keyword(MPOLE)

Specification of magnetic Cartesian multipole operators of order :math:`n`

.. math::

   \hat{m}_{j_{1}\ldots j_{n-1};j_{n}}^{\left[n\right]}=\frac{n}{n+1}r_{j_{1}}r_{j_{2}}\ldots r_{j_{n-1}}(\boldsymbol{r}\times\hat{\mathbf{j}})_{j_{n}};\quad\hat{\mathbf{j}}=-ec\boldsymbol{\alpha}


for the calculation of transition moments (note that they contribute to the same order in the wave vector). Specify order.

*Example:* Magnetic dipole operators::

      .MPOLE
      1


keyword(ANALYZE)

Analyze solution vectors and show the most important excitations at the
orbital level.

keyword(INTENS)

Invoke calculation of oscillator strengths to order k in the wave vector.
Default order is zero, which corresponds to the widely used electric-dipole approximation. By default results are given both in the length and the velocity representation, but a selection can be made using keywords :ref:`EXCITATION_ENERGIES_.NOLENR` and :ref:`EXCITATION_ENERGIES_.NOVELR`. Only even orders contribute. For further details, see :cite:`List_JCP2020`

*Example:* ::

     .INTENS
     2


keyword(BED)

Invoke calculation of oscillator strengths using the full operator coupling the molecule to an electromagnetic plane wave. The oscillator strength is then given by

.. math::

   f_{n\leftarrow 0}=\frac{2\omega}{\hbar e^{2}}\left|\langle n|T\left(\omega\right)|0\rangle\right|^2

where appears the effective interaction operator

.. math::

   T\left(\omega\right)=\frac{ec}{\omega}\left(\boldsymbol{\alpha}\cdot\boldsymbol{\epsilon}\right)e^{+i\left(\mathbf{k}\cdot\mathbf{r}\right)}.

In the above expression :math:`\mathbf{k}` refers to the wave vector of length :math:`k=\omega/c` and direction :math:`\mathbf{m}`, whereas the polarization of the electric component is specificed by :math:`\boldsymbol{\epsilon}`.

Since experiment is typically carried out in an isotropic medium, rotational average is performed.
Rather than rotate the molecule we shall rotate the experimental apparatus.
To rotate we use the unit vectors of the spherical coordinates

.. math::

   \begin{array}{lcl} \mathbf{e}_{r} & = & \mathbf{e}_{x}\sin\theta\cos\phi+\mathbf{e}_{y}\sin\theta\sin\phi+\mathbf{e}_{z}\cos\theta\\
   \mathbf{e}_{\theta} & = & \mathbf{e}_{x}\cos\theta\cos\phi+\mathbf{e}_{y}\cos\theta\sin\phi-\mathbf{e}_{z}\sin\theta\\ \mathbf{e}_{\phi} & = & -\mathbf{e}_{x}\sin\phi+\mathbf{e}_{y}\cos\phi \end{array}   

We now choose to align the wave vector with the :math:`\mathbf{e}_r` unit vector, that is

.. math::

   \mathbf{k} = k\mathbf{e}_r

This means that the polarization vector :math:`\boldsymbol{\epsilon}` is in the plane spanned by the
unit vectors :math:`\mathbf{e}_{\theta}` and :math:`\mathbf{e}_{\phi}`. We therefore set

.. math::

   \boldsymbol{\epsilon} = \cos\chi \mathbf{e}_{\theta}+\sin\chi\mathbf{e}_{\phi}

We see the solid angle :math:`\left(\theta,\phi\right)` gives all possible
directions of the wave vector :math:`\mathbf{k}`, wheras the angle
:math:`\chi` provides all possible orientations of the polarization vector
:math:`\boldsymbol{\epsilon}` in the plane perpendicular to :math:`\mathbf{k}`. 

The general expression for the rotational average will be      

.. math::

   \left\langle f\left(\boldsymbol{r}\right)\right\rangle _{\theta,\phi,\chi}=\frac{1}{8\pi^{2}}\int_{0}^{2\pi}\int_{0}^{2\pi}\int_{0}^{\pi}f\left(\boldsymbol{r}\right)\sin\theta d\theta d\phi d\chi.


In our case we have   

.. math::

   \left\langle f_{n\leftarrow 0}\right\rangle _{\theta,\phi,\chi} = \frac{2\omega}{\hbar e^{2}}\left\langle \epsilon_{\alpha}\epsilon_{\beta}\langle n|\frac{ec}{\omega}\alpha_{\alpha}e^{+i\left(\boldsymbol{k}\cdot\boldsymbol{r}\right)}|0\rangle\langle n|\frac{ec}{\omega}\alpha_{\beta}e^{+i\left(\boldsymbol{k}\cdot\boldsymbol{r}\right)}|0\rangle^{\ast}\right\rangle _{\theta,\phi,\chi},

which simplifies to   

.. math::

   \left\langle f_{n\leftarrow 0}\right\rangle _{\theta,\phi,\chi} = \frac{2\omega}{\hbar e^{2}}\left\langle \left\langle\epsilon_{\alpha}\epsilon_{\beta}\right\rangle _{\chi}\langle n|\frac{ec}{\omega}\alpha_{\alpha}e^{+i\left(\boldsymbol{k}\cdot\boldsymbol{r}\right)}|0\rangle\langle n|\frac{ec}{\omega}\alpha_{\beta}e^{+i\left(\boldsymbol{k}\cdot\boldsymbol{r}\right)}|0\rangle^{\ast}\right\rangle _{\theta,\phi},

since only the polarization vectors depend on the angle :math:`\chi`. The average over the angle :math:`\chi` can be expressed compactly in terms of the wave unit vector :math:`\mathbf{m}` (:math:`\mathbf{k}=k\mathbf{m}`)

.. math::

  \left\langle\epsilon_{\alpha}\epsilon_{\beta}\right\rangle _{\chi}=\frac{1}{2}\left(\delta_{\alpha\beta}-m_{\alpha}m_{\beta}\right),

whereas the average over angles :math:`\theta` and :math:`\phi` is handled by `Lebedev quadrature <https://en.wikipedia.org/wiki/Lebedev_quadrature>`_ .

A full account is given in :cite:`List_JCP2020`

keyword(ORIENT)

Specify fixed experimental configuration (no rotational average). The orientation of the wave and polarization vector is given by specification of the angles :math:`\theta`, :math:`\phi` and :math:`\chi`, see the :ref:`EXCITATION_ENERGIES_.BED` keyword for more details. For instance, to specify that the wave vector is along the :math:`z` -axis  and the polarization vector along the :math:`x` - axis, we set

::

    .ORIENT
    0.0 0.0 0.0

keyword(NROTAV)

As described under the .BED keyword,  `Lebedev quadrature <https://en.wikipedia.org/wiki/Lebedev_quadrature>`_ is employed for rotational average. This quadrature over the solid angles can integrate a spherical harmonic to high accuracy with a maximum angular momentum :math:`L_\mbox{max}`. The default value of :math:`L_\mbox{max}` is presently 5, but can be reset with this keyword.

keyword(BEDCON)

Specification of contributions of the full light-matter interaction of order :math:`n` in the
wave vector

.. math::

   \hat{T}_{\mathrm{full}}^{\left[n\right]}(\omega)=\frac{k^{n}}{n!}\frac{d^{n}}{dk^{n}}\left[\frac{e}{\omega}
   \left(c\boldsymbol{\alpha}\cdot\boldsymbol{\epsilon}\right)e^{+i\left(\boldsymbol{k}\cdot\boldsymbol{r}\right)}\right]_{k=0}
   =\frac{e}{\omega}\frac{i^{n}}{n!}\left(c\boldsymbol{\alpha}\cdot\boldsymbol{\epsilon}\right)\left(\boldsymbol{k}\cdot\boldsymbol{r}\right)^{n}

keyword(NOLENR)

Deactivate length representation.

keyword(NOVELR)

Deactivate velocity representation.


Control variational parameters
==============================

keyword(OCCUP)

For each fermion ircop give an :ref:`orbital_strings` of inactive orbitals from
which excitations are allowed. By default excitations from all occupied
orbitals are included in the generalized eigenvalue problem.

Example: ::

    .OCCUP
    1..3
    7,8

This would include excitations from gerade orbitals 1,2,3, and ungerade
orbitals 7 and 8.

keyword(VIRTUA)

For each fermion ircop give an :ref:`orbital_strings`
of virtual orbitals
to which excitations are allowed. By default excitations to all virtal
orbitals are included in the generalized eigenvalue problem.

keyword(SKIPEE)

Exclude all rotations between occupied positive-energy and virtual
positive-energy orbitals.

keyword(SKIPEP)

Exclude all rotations between occupied positive-energy and virtual
negative-energy orbitals.

Control reduced equations
=========================

keyword(MAXITR)

Maximum number of iterations.

*Default:* ::

    .MAXITR
     30

keyword(MAXRED)

Maximum dimension of matrix in reduced system.

*Default:* ::

    .MAXRED
     200

keyword(THRESH)

Threshold for convergence of reduced system.

*Default:* ::

    .THRESH
     1.0D-5

Control integral contributions
==============================

The user is encouraged to experiment with these options since they may
have an important effect on run time.

keyword(INTFLG)

Specify what two-electron integrals to include
(default: :ref:`HAMILTONIAN_.INTFLG` under :ref:`**HAMILTONIAN`).

keyword(CNVINT)

Set threshold for convergence before adding SL and SS integrals to
SCF-iterations.

*2 (real) Arguments:* ::

    .CNVINT
     CNVXQR(1) CNVXQR(2)

*Default:* Very large numbers.

keyword(ITRINT)

Set the number of iterations before adding SL and SS integrals to
SCF-iterations.

*Default:* ::

    .ITRINT
     1 1

Advanced/debug flags
====================

keyword(E2CHEK)

Generate a complete set of trial vector which implicitly allows the
explicit construction of the electronic Hessian. Only to be used for
small systems !

keyword(ONLYSF)

Only call FMOLI in sigmavector routine: only generate one-index
transformed Fock matrix  :cite:`Saue2003`.

keyword(ONLYSG)

Only call FMOLI in sigmavector routine: 2-electron Fock matrices using
one-index transformed densities :cite:`Saue2003`.

keyword(GNOISE)

To test the robustness of property gradients to numerical noise artificial noise is added to
the MO-coefficients. More precisely, the user activates noise and provides a "noise level".
Then, for each element of the coefficient array, a pseudo-random number in the interval (-1,+1]
is selected, multiplied with the "noise level" and added to the element.

*Example:* ::

    .GNOISE
     1.0D-10
