:orphan:

==============================
Molecular rotational g-tensors
==============================

Introduction
------------

In this tutorial we introduce the calculation of molecular rotational g-tensors as given in the DIRAC code, based on the theoretical
developments by I. Agustín Aucar *et al.* For details, you are welcome to consult :cite:`Aucar_JCP2014`.

The molecular rotational g-tensor elements are given by

.. math::

   g_{\alpha\beta} = g_{\alpha\beta}^{nuc} + g_{\alpha\beta}^{elec}

They have two terms: the first of them is independent of the electronic variables, whereas the second one is given by a linear response function.

In tensorial notation, the rotational g-tensor contributions of a molecule are given by:

.. math::

   {\bf g}^{nuc} = m_p \; \left\{ \sum_{M} Z_M \left[ \left( {\bf R}_{M,GO} \cdot {\bf R}_{M,CM} \right) {\bf 1} - {\bf R}_{M,GO} {\bf R}_{M,CM} \right] \right\} \otimes {\bf I}^{-1}

.. math::

   {\bf g}^{elec} = {\bf g}^{LR} = m_p \;  \langle \langle  \left ({\bf r}-{\bf R}_{GO} \times c \, {\bf \alpha}\right) \; ; \; {\bf J}_e\rangle\rangle \; \otimes  \; {\bf I}^{-1}

where :math:`{\bf R}_{M,GO}` is the position of nucleus :math:`M` with respect to the gauge origin; :math:`{\bf R}_{M,CM}` is the position of nucleus :math:`M` with respect to the molecular center of mass; :math:`{\bf R}_{GO}` is the position of the gauge origin; :math:`{\bf I}` is the inertia tensor of the molecule and :math:`{\bf J}_e = \left({\bf r} - {\bf R}_{CM} \right) \times {\bf p}+{\bf S}_e` is the electronic total angular momentum.


Application to the HF molecule
------------------------------

As an example, we show a calculation of the g-tensor of the Hydrogen fluoride molecule. The input file `rotg.inp` is given by

.. literalinclude:: rotg.inp

whereas the molecular input file `HF_cv3z.mol` is

.. literalinclude:: HF_cv3z.mol

The calculation is run using::

   pam --inp=rotg --mol=HF_cv3z


As a result, g-factors are obtained at the coupled Hartree-Fock level. The code also works at the DFT level.

It is also interesting to note that as :ref:`HAMILTONIAN_.URKBAL` is requested in the present calculation, the results will be very close to those obtained using the RKB prescription, because the diamagnetic-like (e-p) contributions to :math:`g_{\alpha\beta}^{LR}` are almost zero (see Eq. 41 of :cite:`Aucar_JCP2014`).


Reading the output file
-----------------------

As the :ref:`PROPERTIES_.PRINT` flag in the input file is set to 1, the results are fully detailed.

The g-factor Hydrogen fluoride molecule will look like:

.. literalinclude:: rotg-HF

One should recall that in the case of linear molecules, as the present one, only one tensor element is printed out (the g-factor) because for these molecules the g-tensor has only two equal and non-zero elements.

As it can be seen, the total g-tensor of the Hydrogen fluoride molecule is given, and then separated in its two terms (nuc and LR). The latter contribution, the linear response function, is further separated in their (e-e) and (e-p) linear response parts, as well as in their :math:`\mathbf{L}` and :math:`\mathbf{S}` parts. It is seen how the (e-p) contribution to the linear response term of the g-factor is almost zero.
