:orphan:

====================================
Excitation energies at the SCF-level
====================================

In this tutorial we will look at the theory behind the calculation of excitation energies at the SCF-level, that is, 
either by time-dependent Hartree-Fock (TD-HF) or TD-DFT. The implementations in DIRAC are presented in :cite:`Bast2009`.

Excitation energies :math:`\omega_m` are found by solving the generalized eigenvalue problem

.. math::

   \left(E_0^{[2]}-\hbar\lambda S^{[2]}\right)X_m=0

where the matrices :math:`E_0^{[2]}` and :math:`S^{[2]}` are the electronic Hessian and generalized metric, respectively.
The eigenvector have a paired structure

.. math::

   \left(\lambda_{+;m}=+\left|\omega_m\right|, X_{+;m} = \right)

.. math::

  \left(\begin{matrix} a & b \\
    c & d \end{matrix}\right)

.. math::

  a & = (x + y)^2 \\
    & = x^2 + 2xy + y^2

