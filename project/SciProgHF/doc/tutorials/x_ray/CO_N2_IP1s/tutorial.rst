:orphan:
 

Core ionization in the CO and N2 molecules
==========================================

Introduction
============

We would like to study K-edge energies of CO and N2 molecules.

Carbon monoxide
===============

A first approximation to the O1s and C1s  ionization energies is provided by Koopmans' theorem.
Starting from the molecular input file `CO.mol` 

.. literalinclude:: CO.mol

and the menu file `CO.inp` 

.. literalinclude:: CO.inp

we first run a Hartree-Fock calculation of the ground electronic state of carbon monoxide::

   pam --inp=CO --mol=CO --get "DFCOEF=cf.CO"

We find that :math:`-\varepsilon (O1s)=20.6838\,E_h=562.8` eV and :math:`-\varepsilon (C1s)=11.3672 \,E_h=309.3` eV. This is significantly off the values of **542.1 eV** and **295.9 eV**, respectively, reported by Kai Siegbahn and co-workers :cite:`Siegbahn1969` .

Koopmans' theorem work quite well for valence excitations because the two major sources, that is, i) lack of electron correlation and ii) lack of orbital relaxation, typically are of about the same magnitude, but of opposite sign and hence tend to cancel out. This is not the case for core ionization, where orbital relaxation is very much more important
than electron correlation. To introduce orbital relaxation of the core-ionized state we carry out a :math:`\Delta SCF`
calculation, as pioneered by Paul Bagus and Henry F. Schaefer :cite:`Bagus_JCP1971`.

In order to calculate the O1s ionized system we set up the menu file `CO_O1s.inp` 

.. literalinclude:: CO_O1s.inp

and use the command::

  pam --inp=CO_O1s --mol=CO --put "cf.CO=DFCOEF"

It is seen that we start from the coefficients of the neutral CO molecule. We then converge to the core-ionized system using reordering and overlap selection: In the input we have specified twelve electrons in six inactive (closed) orbitals, followed
by a single electron in an active (open) orbital. We want the O1s to be the active orbital, but this is not achieved automatically since DIRAC will normally order orbitals according to their energy.
We therefore start be reordering the orbitals such that the O1s orbital from the previous calculation on the neutral system comes out on top of the occupied orbitals. However, this is not enough to converge to the desired state since after the first diagonalization DIRAC will again by default order orbitals according to their energy. This is why we use *overlap selection*, that is, we ask DIRAC to
rather order orbitals according to their overlap with some reference orbitals. By default (dynamic overlap selection) this will be the orbitals from the previous iteration. However, in this case we activate *non-dynamic* overlap selection, which means that we order orbitals according to their overlap with the starting orbitals.

The analogous menu file for C1s ionization is `CO_C1s.inp` 

.. literalinclude:: CO_C1s.inp

which puts the C1s on top of the occupied orbitals, and then use the command::

  pam --inp=CO_C1s --mol=CO --put "cf.CO=DFCOEF"

The HF energies of the neutral, the O1s-ionized and the C1s-ionized systems are -112.857912 :math:`E_h`, -92.939954 :math:`E_h` and -101.933919 :math:`E_h`, respectively, from which we calculate IP(O1s) = 19.9180 :math:`E_h` = 542.0 eV
and IP(O1s) = 10.9240 :math:`E_h` = 297.3 eV, which agrees remarkably well with the experimental values, in particular
for the oxygen K-edge. 

 .. note:: Overlap selection is nowadays marketed hard as MOM (Maximum Orbital Method, see :cite:`Gilbert_JPCA2008`), but this method has been included in DIRAC for at least two decades and goes back to the pioneering work of `Paul Bagus <http://cascam.unt.edu/people/psbagus.htm>`_ It was used in :cite:`Bagus_JCP1971`, but not reported explicitly. However, it is for instance documented in the `1975 manual of the ALCHEMY program <http://k-sek01.t-komazawa.ac.jp/msekiya/alchemy/scfm.pdf>`_ (On pdf page 15 you find a description of keyword MOORDR using a "maximum overlap criterion").

Nitrogen molecule
=================

We now switch to the isoelectronic nitrogen molecule. Starting from the molecular input file `N2.mol` 

.. literalinclude:: N2.mol

and the menu file `N2.inp` 

.. literalinclude:: N2.inp

we run a Hartree-Fock calculation of the ground electronic state of the nitrogen molecule::

   pam --inp=N2 --mol=N2 --get "DFCOEF=cf.N2 DFACMO=ac.N2"

(note that we are also getting the :math:`C_1` coefficients; we shall come back to that). We find that :math:`-\varepsilon (1s\sigma_g)=15.6942\,E_h=427.1` eV and :math:`-\varepsilon (1s\sigma_u)=15.6907 \,E_h=427.0` eV. This is significantly off the value of **409.9 eV** reported by Kai Siegbahn and co-workers :cite:`Siegbahn1969` ; the :math:`1s\sigma_g` and :math:`1s\sigma_u` are indistinguishable in the X-ray PES study (see also :cite:`Wight_JEPRP1972`).

We therefore proceed to :math:`\Delta SCF` calculations: for the :math:`1s\sigma_g` ionization we use the menu file `N2_Kg.inp`

.. literalinclude:: N2_Kg.inp

and the command::

  pam --inp=N2_Kg --mol=N2 --put "cf.N2=DFCOEF"

and equivalent menu file `N2_Ku.inp` for the :math:`1s\sigma_u` ionization:
  
.. literalinclude:: N2_Ku.inp

The HF energies of the neutral, the :math:`1s\sigma_g` -ionized and the :math:`1s\sigma_u` -ionized systems are -109.051015 :math:`E_h`,
-93.626378 :math:`E_h` and -93.630989 :math:`E_h`, respectively, from which we calculate IP(:math:`1s\sigma_g`) = 15.4246 :math:`E_h` = 419.7 eV
and IP(:math:`1s\sigma_u`) = 15.4200 :math:`E_h` = 419.6 eV. However, *this is still far from the experimental value*
of 409.9 eV !

Bagus and Schaefer found a similar discrepancy in their study of the :math:`O_2^+` molecule
:cite:`Bagus_JCP1972` and found that results improved dramatically by *localization* of the core hole,
so let us try this. Localization breaks molecular symmetry, so we go down to :math:`C_1` using the molecule input file `N2_C1.mol` 

.. literalinclude:: N2_C1.mol

This is the reason why we saved the :math:`C_1` coefficients in the initial calculation !
We next set up the menu file `N2loc.inp`

.. literalinclude:: N2loc.inp

for Pipek-Mezey localization. Note that we select only the a little trick: We select only the :math:`1s\sigma_g` and :math:`1s\sigma_u` orbitals for localization. We now localize using the command::

  pam --inp=N2loc --mol=N2_C1 --put "ac.N2=DFCOEF" --get "DFCOEF=ac.N2loc"

We calculate the system with a localized core hole using the menu file  `N2_K.inp`

.. literalinclude:: N2_K.inp

and the command::

  pam --inp=N2_K --mol=N2_C1 --put "ac.N2loc=DFCOEF"

In fact, we can do slightly better because we can go up to :math:`C_{2v}` symmetry using that DIRAC allows
the import of coefficients in :math:`C_1` symmetry from the unformatted file DFACMO to current symmetry.
This assumes that the current symmetry is lower than the symmetry used for obtaining the original coefficients,
which is the case here. We now calculate the system with a localized core hole using the menu file  `N2_KC2v.inp`

.. literalinclude:: N2_KC2v.inp

the molecule input file `N2_C2v.mol`

.. literalinclude:: N2_C2v.mol

and the command::

  pam --inp=N2_KC2v --mol=N2_C2v --put "ac.N2loc=DFACMO"

The calculated energy of the core ionized system is -93.971190  :math:`E_h` and we now obtain an ionization energy
of 15.0978 :math:`E_h` = 410.3 eV, in much better agreement with experiment. The price we had to pay was to break
symmetry. It is possible to keep symmetry, but now at the price of going to more complex wave functions, as shown
by Agren, Bagus and Roos :cite:`Agren_CPL1981`.


