Spin-orbit splitting from spin-free orbitals
============================================

The test demonstrates the spin-orbit (SO) splitting, which is generated from spin-free (SF) molecular orbitals.
This is sometimes called the perturbative SO approach.

SCF 
---

At the SCF level, we watch the 2p12, 2p32 spinors of the Rn atom. We generate them from spinfree orbitals
in one iteration using the X2C Hamiltonian, containing  both scalar and spinorbital effects.

One reads SF orbitals from the file and checks orbital energies:

::

  pam --noarch --put="DFPCMO.Rn.x2c_sf.scf.v2z.D2h=DFPCMO" --mol=Rn.v2z.D2h.mol  --inp=Rn.x2c.scf_2fs_1iter.inp

The 2p12,2p32 orbital energies data are in the following table:

           SF            SO from SF           SO
   -570.15736962626  -641.16004038807   -641.78593855634
   -570.15736962625  -541.27583741512   -540.66722888765
   -570.15736962625  -541.27583741511   -540.66722888763

COSCI
-----

In the following case we focus on the Tl tom fine-structure splitting, caused by one electron in the 6p valence shell.
For that we employ the COSCI method, which constructs small determinantal expansion in the space of (active) open-shell orbitals,
giving energy of the the ground 2P12 state (2-fold degenerate) and  of the first excited 2P32 state (4-fold degenerate).

We use the SF orbitals for the open-shell resolve (COSCI) step:

::
  pam --noarch --put="DFPCMO.Tl.x2c_sf.scf.v2z.D2h=DFPCMO" --mol=Tl.v2z.D2h.mol  --inp=Tl.x2c.resolve_2fs.inp

We get the value of 6411.5 cm-1.

However, with the full relativistic Hamiltonian we get the 2P12-2P32 splitting of 7752.4 cm-1.
In the test, this is shown when one employs the full relativistic (SF+SO) molecular orbitals file for the quick open-shell CI:

::
  pam --noarch --mol=Tl.v2z.D2h.mol  --inp=Tl.x2c.resolve_2fs.inp --put="DFPCMO.Tl.x2c.scf.v2z.D2h=DFPCMO"

Note that all calculations are performed in the highest possible symmetry, D2h, for both SF and SO.
