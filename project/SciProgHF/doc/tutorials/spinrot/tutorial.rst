:orphan:

===============================
Nuclear spin-rotation constants
===============================

Introduction
------------

In this tutorial we introduce the calculation of nuclear spin-rotation tensors as given in the DIRAC code, based on the theoretical
developments by I. Agustín Aucar *et al.* For details, you are welcome to consult :cite:`Aucar_JCP2012`.

The nuclear spin-rotation (SR) tensor elements of a nucleus :math:`N` are given by

.. math::

   M_{N,\alpha\beta} = M_{N,\alpha\beta}^{NU} + M_{N,\alpha\beta}^{EV} + M_{N,\alpha\beta}^{LR}

They have three terms: the first of them is independent of the electronic variables, whereas the second and third are given by an expectation value (EV) and a linear response (LR) function, respectively.

In tensorial notation, the nuclear spin-rotation contributions of a nucleus :math:`N` are given, in SI units, by:

.. math::

   {\bf M}_N^{NU} =& \; \frac{1}{4\pi\epsilon_0c^2} \, \frac{e^2 \hslash^2}{2 m_p h} \, g_N  \sum_{M \neq N} \frac{Z_M}{|{\bf R}_{MN}|^3} \\
                   & \hspace{2cm} \left(  \left\{ \left[ {\bf R}_{M,CM} - \left( 1 - \frac{Z_N m_p}{m_N g_N} \right) {\bf R}_{N,CM} \right] \cdot {\bf R}_{MN} \right\} {\bf 1}  \right. \\
                   & \hspace{4cm} \left. - \left[ {\bf R}_{M,CM} - \left( 1 - \frac{Z_N m_p}{m_N g_N} \right) {\bf R}_{N,CM} \right] {\bf R}_{MN} \right)   \otimes {\bf I}^{-1}  \\
                   & \\
                   & \\
   {\bf M}_N^{EV} =& \; \frac{1}{4\pi\epsilon_0c^2} \, \frac{e^2 \hslash^2}{2 m_p h} \, g_N \left( 1 - \frac{Z_N m_p}{m_N g_N} \right) \\
                   & \hspace{2cm} \left[ \left(  \langle 0 | \frac{{\bf r}-{\bf r}_N}{|{\bf r} -{\bf r}_N|^3 } | 0 \rangle \cdot {\bf R}_{N,CM} \right) {\bf 1} - \langle 0 | \frac{{\bf r}-{\bf r}_N}{|{\bf r} -{\bf r}_N|^3 } | 0 \rangle  \; \; {\bf R}_{N,CM}\right] \otimes {\bf I}^{-1} \\
                   & \\
                   & \\
   {\bf M}_N^{LR} =& \; \frac{1}{4\pi\epsilon_0c^2} \, \frac{e^2 \hslash^2}{2 m_p h} \, g_N  \;  \langle\langle \; \left (\frac{{\bf r}-{\bf r}_N}{|{\bf r} -{\bf r}_N|^3 }\times c \, {\bf \alpha}\right) \; ; \; {\bf J}_e \; \rangle\rangle \; \otimes  \; {\bf I}^{-1}

where :math:`{\bf R}_{MN}` is the position of nucleus :math:`M` with respect to the position of nucleus :math:`N`; :math:`{\bf R}_{N,CM}` is the position of nucleus :math:`N` with respect to the molecular center of mass; :math:`{\bf I}` is the inertia tensor of the molecule and :math:`{\bf J}_e = \left({\bf r} - {\bf R}_{CM} \right) \times {\bf p}+{\bf S}_e` is the electronic total angular momentum.

For molecules at their equilibrium geometry, it is found that the sum of the first two terms, :math:`{\bf M}_N^{NU} (eq) + {\bf M}_N^{EV} (eq)`, is equal to a new tensor which is independent of the electronic dynamics (see Eq. 60 of :cite:`Aucar_JCP2012`),

.. math::

   {\bf M}_N^{nuc} =& \; {\bf M}_N^{NU} (eq) \; + \; {\bf M}_N^{EV} (eq) \\
                           =& \; \frac{1}{4\pi\epsilon_0c^2} \, \frac{e^2 \hslash^2}{2 m_p h} \, g_N \sum_{M \neq N} Z_M
                              \left[  \left( {\bf R}_{M,CM}  \cdot \frac{{\bf R}_{MN}}{|{\bf R}_{MN}|^3} \right) {\bf 1}
                                    - {\bf R}_{M,CM} \frac{{\bf R}_{MN}}{|{\bf R}_{MN}|^3} \right]   \otimes {\bf I}^{-1}

Therefore, the electronic dependence of the nuclear spin-rotation tensor of a nucleus :math:`N` in a molecule in equilibrium is completely given by its linear response term, :math:`{\bf M}_N^{LR}` (see Eq. 59 of :cite:`Aucar_JCP2012`).


.. note::
   
   The current implementation gives by default (:ref:`PROPERTIES_.PRINT` values up to 3) results only for molecules in equilibrium.


Application to the HF molecule
------------------------------

As an example, we show a calculation of the SR constant of the fluorine nucleus at the Hydrogen fluoride molecule. The input file `spinrot.inp` is given by

.. literalinclude:: spinrot.inp

whereas the molecular input file `HF_cv3z.mol` is

.. literalinclude:: HF_cv3z.mol

The calculation is run using::

   pam --inp=spinrot --mol=HF_cv3z


As a result, SR constants are obtained at the coupled Hartree-Fock level. The code also works at the DFT level.

It is also interesting to note that as :ref:`HAMILTONIAN_.URKBAL` is requested in the present calculation, the results will be very close to those obtained using the RKB prescription, because the diamagnetic-like (e-p) contributions to :math:`M_{N,\alpha\beta}^{LR}` are almost zero (see Eq. 60 of :cite:`Aucar_JCP2012`).


Reading the output file
-----------------------

As the :ref:`PROPERTIES_.PRINT` flag in the input file is set to 4, the results are fully detailed and given in both kHz and ppm units.

The SR constant of the fluorine nucleus, in kHz, will look like:

.. literalinclude:: SR-HF-khz

One should recall that in the case of linear molecules, as the present one, only one tensor element is printed out (the spin-rotation constant) because for these molecules the SR tensor has only two equal and non-zero elements.

As it can be seen, the total SR constant of the fluorine nucleus is given, and then separated in its two terms (nuc and LR). The latter contribution, the linear response function, is further separated in their (e-e) and (e-p) parts, as well as in their :math:`\mathbf{L}` and :math:`\mathbf{S}` parts. It is seen how the (e-p) contribution to the linear response term of the SR constants is almost zero.

In addition, the results are fully detailed, showing the three terms mentioned above (NU, EV and LR).

Finally, if one is interested in the relationship between the SR constants and :ref:`PROPERTIES_.SHIELDING` (see :cite:`Aucar_JPCL2016` and :cite:`AucarChap2019`), which is given, in SI units, by

.. math::

   {\bf \sigma}_N = \; \frac{1}{4\pi\epsilon_0c^2} \, \frac{e^2}{2} \;  \langle\langle \; \left (\frac{{\bf r}-{\bf r}_N}{|{\bf r} -{\bf r}_N|^3 }\times c \, {\bf \alpha}\right) \; ; \; \left( {\bf r}-{\bf r}_{GO} \times c \, {\bf \alpha}\right) \; \rangle\rangle


then we print also the SR constants in ppm, by multiplying each SR tensor element :math:`M_{N,\alpha\beta}` (given in kHz) by :math:`\frac{m_p I_{\beta \gamma}}{g_N}\,\frac{2\, \pi \, 10^9}{m_e \, \hslash}`, where :math:`m_p` and :math:`m_e` are the proton and electron masses, respectively, :math:`\hslash` is the reduced Planck's constant, and :math:`g_N` is the g-factor of nucleus `N`.


In this way, one obtains the following output:

.. literalinclude:: SR-HF-ppm
