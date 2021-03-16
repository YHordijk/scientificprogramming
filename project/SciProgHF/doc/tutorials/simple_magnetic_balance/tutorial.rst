:orphan:

.. _NMRShieldSimpleMagBal:

===========================================================
Calculation of NMR shieldings using simple magnetic balance
===========================================================

In this tutorial we look at 4-component relativistic calculations of NMR
shielding constants. For best results we recommend to use the scheme of simple
magnetic balance (sMB) in conjunction with London orbitals :cite:`Olejniczak2012`.


Basic theory of magnetic balance
================================

We encourage the reader to consult the original paper :cite:`Olejniczak2012` for
full details; here we just provide some key points.

In the absence of magnetic fields the coupling between the large and small
components of the Dirac equation is

.. math::

  c \psi^S = \frac{1}{2m} \left[1 + \frac{E-V}{2mc^2} \right]^{-1} (\vec \sigma \cdot \vec p) \psi^L

The exact coupling can not be implemented *per se* at the basis set level due to
the energy dependence of the coupling.  One therefore employs the
non-relativistic limit of this coupling

.. math::

  \lim_{c \to \infty} c \psi^S = \frac{1}{2m} (\vec \sigma \cdot \vec p) \psi^L,

valid for positive-energy solutions (and extended nuclear charge
distributions). Construction of the small component basis in this manner
assures the correct representation of the kinetic energy operator in the
non-relativistic limit, and the coupling is therefore referred to as *kinetic
balance* :cite:`Stanton1984`.  In a 2-spinor basis kinetic balance can be
implemented in such a manner that there is a 1:1 ratio between the sizes of the
large and small component basis; this is referred to as *restricted kinetic
balance* (RKB).

.. math::

  \chi_{\mu}^{2c;S} = (\vec \sigma \cdot \vec p) \chi_{\mu}^{2c;L}

In a scalar basis the straightforward implementation of RKB is not possible and
one typically resorts to *unrestricted kinetic balance* (UKB)

.. math::

  \{ \chi_{\mu}^{1c;S} \} = \{ \vec p \chi_{\mu}^{1c;L} \},

basically generating the small component basis functions as derivative of the
large component ones.

Magnetic fields are introduced through minimal substitution

.. math::

  \vec p \to \vec \pi = \vec p + e \vec A

which manifestly modifies the coupling between the large and the small
components. The non-relativistic limit now reads

.. math::

  \lim_{c \to \infty} c \psi^S = \frac{1}{2m} (\vec \sigma \cdot \vec \pi) \psi^L

and is referred to as *magnetic balance* :cite:`Aucar1999`.  Magnetic balance needs
to be taken into account when constructing the first-order correction to the
wave function upon inclusion of an external magnetic field and is not assured
in a RKB basis :cite:`Pecul2004`.  In the cited paper :cite:`Olejniczak2012` we show that
magnetic balance can be obtained in a Gaussian basis combining UKB with London
orbitals. However, UKB generally leads to a small component basis larger than
the large component basis set from which it was generated and may increase
computational cost as well as introduce linear dependencies in the basis set as
compared to RKB. We therefore propose calculations of NMR shielding parameters
by a simple two-step procedure described in the next section.


A sample calculation
====================

We will consider a 4-component relativistic density functional calculation of
the NMR shielding constant of hydrogen fluoride using the KT2 functional
:cite:`Keal2003`.  The molecular input file ``hf.mol`` reads:

.. literalinclude:: hf.mol


Step 1: SCF calculation
-----------------------

The first step is the generation of the Kohn--Sham wave function in a RKB
basis. Although DIRAC employs scalar Gaussian basis functions restricted
kinetic balance is by default imposed when transforming to orthonormal basis
:cite:`Visscher2000`.

The menu file ``scf.inp`` reads:

.. literalinclude:: scf.inp

and we make sure to recover the coefficient file ``DFCOEF`` when running the calculation:

.. literalinclude:: scf.run


Step 2: Response calculation
----------------------------

In the second step we start from the RKB coefficient obtained in the SCF calculation, but DIRAC will add the UKB complement
and
calculate the NMR shielding constant using unrestricted kinetic balance
combined with London orbitals. The menu file ``response.inp`` accordingly reads:

.. literalinclude:: response.inp

and the ``pam`` command is:

.. literalinclude:: response.run


Comparison with RKB calculations
--------------------------------

With the above computational procedure we obtain the following NMR shielding constants::

             isotropic shielding

         ----------------------------

  atom          total          dia         para

  ---------------------------------------------------

  H            29.9502       21.8737        8.0765

  F           415.8316      480.2015      -64.3699

  ---------------------------------------------------

By not adding the UKB complement (deleting .RKBIMP and .URKBAL), that is, using RKB, gives::

              isotropic shielding

          ----------------------------

  atom          total          dia         para

  --------------------------------------------------------

  H            28.8028       20.7264        8.0765

  F           323.2672      387.6459      -64.3787

  --------------------------------------------------------

If we in addition switch off the use of London-orbitals (deleting .LONDON) we obtain::

              isotropic shielding

          ----------------------------

  atom          total          dia         para

  ------------------------------------------------------

  H            28.2920       20.2563        8.0357

  F           322.9065      387.1375      -64.2310

  ------------------------------------------------------

The experimental value for the absolute isotropic shielding available for H in HF is 28.72 ppm :cite:`Chesnut1994`.
