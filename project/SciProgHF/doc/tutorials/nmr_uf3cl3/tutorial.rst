:orphan:


===================================
Case study: NMR shielding of UF3Cl3
===================================

In this tutorial we shall go through the calculation of NMR shieldings of the
UF3Cl3 molecule using the scheme of simple magnetic balance (SMB) in
conjunction with London orbitals :cite:`Olejniczak2012`.


Preparing for the atomic start
==============================

To ensure a smooth convergence of the SCF calculation we employ an atomic start
whereby an initial molecular density matrix is generated as the direct sum of
atomic ones. We therefore have to perform a separate SCF calculation of each
atomic type. At the moment the approach is somewhat cumbersome, but
will be automatized in future versions of DIRAC. Each atom is calculated in
maximal symmetry (linear symmery) and the coefficients then exported to
C1 symmetry using the keyword :ref:`GENERAL_.ACMOUT`.


Uranium atom
------------

We calculate the U+3 cation in the [Rn] 5f3 configuration.
The input file ``UIII.inp`` reads

.. literalinclude:: UIII.inp

and the corresponding molecular input file ``U.mol`` is

.. literalinclude:: U.mol

The minimal pam command is::

  pam --inp=UIII.inp --mol=U.mol --get=DFACMO

We save the C1 coefficients and make a symbolic link for the atomic start::

  mv DFACMO ac.UIII
  ln -s ac.UIII DFUXXX


Fluorine atom
-------------

We calculate the fluorine atom in the [He] 2s2 2p5 ground
state configuration using fractional occupation.
The input file ``F.inp`` reads

.. literalinclude:: F.inp

and the corresponding molecular input file ``F.mol`` is

.. literalinclude:: F.mol

The minimal pam command is::

  pam --inp=F.inp --mol=F.mol --get=DFACMO

We save the C1 coefficients and make a symbolic link for the atomic start::

  mv DFACMO ac.F
  ln -s ac.F DFFXXX


Chlorine atom
-------------

We calculate the chlorine atom in the [Ne] 3s2 3p5 ground
state configuration using fractional occupation.
The input file ``Cl.inp`` reads

.. literalinclude:: Cl.inp

and the corresponding molecular input file ``Cl.mol`` is

.. literalinclude:: Cl.mol

The minimal pam command is::

  pam --inp=Cl.inp --mol=Cl.mol --get=DFACMO

We save the C1 coefficients and make a symbolic link for the atomic start::

  mv DFACMO ac.Cl
  ln -s ac.Cl DFCLXX


Initial SCF calculation
=======================


Atomic Hartree--Fock start
--------------------------

In this case the atomic start at the PBE0 level did not converge. It is then
recommend to proceed via a Hartree--Fock calculation where the HOMO-LUMO gap is
larger. The HF menu file ``SCFatom.inp`` reads

.. literalinclude:: SCFatom.inp

For the atomic start it can be noted that orbitals 1 through 43 of the uranium atom are fully occupied whereas the 5f orbitals (44 through 50) are given an fractional occupation of 3/14.

The corresponding molecular file ``UF3Cl3.mol`` reads

.. literalinclude:: UF3Cl3.mol

The minimal pam command is::

  pam --mol=UF3Cl3.mol --inp=SCFatom.inp --outcmo --copy="DFUXXX DFCLXX DFFXXX"

(I also added --mw=120 as well as --mpi=24)


Initial PBE0 run with restricted kinetic balance
------------------------------------------------

Having good start vectors I next ran the PBE0 SCF calculation using the input ``SCF.inp``:

.. literalinclude:: SCF.inp

I am running with a level shift of 0.2 hartree as well as static overlap selection to ensure convergence. I also use the keyword
:ref:`GRID_.ATSIZE` to determine the partitioning of the molecular volume into atomic ones. The minimal pam command reads::

  pam --inp=SCF.inp --mol=UF3Cl3.mol --incmo --outcmo

This calculations converges smoothly in 20 iterations to::

                                   TOTAL ENERGY
                                   ------------

   Electronic energy                        :    -31794.704487527993

   Other contributions to the total energy

   Nuclear repulsion energy                 :      2038.659981883082
   SS Coulombic correction                  :         0.007588768574

   Sum of all contributions to the energy

   Total energy                             :    -29756.036916876335


Property calculation using London orbitals and simple magnetic balance
======================================================================

Starting from the RKB coefficients we proceed to the actual calculation of NMR
parameters using simple magnetic balance and London orbitals. The input file
``NMR.inp`` reads

.. literalinclude:: NMR.inp

This calculations gives::

  @            isotropic shielding
  @        ----------------------------
  @atom          total          dia         para            skew            span            anis            asym
  @----------------------------------------------------------------------------
  @U         -5418.3074      288.0612    -5706.3685       -0.1585    13571.5391     9641.0019        1.2231
  @Cl1  1    -1023.9194       67.4482    -1091.3676        0.8381     2875.6960     2759.3185        0.1265
  @Cl1  2    -1023.9194       67.4482    -1091.3676        0.8381     2875.6960     2759.3185        0.1265
  @Cl3        -885.2205       85.8218     -971.0423        0.7133     2705.2946     2511.3830        0.2316
  @F1   1     -690.9071      103.1746     -794.0817        0.8742     1380.6609     1337.2555        0.0974
  @F1   2     -690.9071      103.1746     -794.0817        0.8742     1380.6609     1337.2555        0.0974
  @F2         -787.3430       84.2708     -871.6137        0.7175     1704.0367     1583.6803        0.2280
  @----------------------------------------------------------------------------
