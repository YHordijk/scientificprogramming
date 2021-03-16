DIRAC open-shell tests
========================

This test is dedicated to demonstrate Dirac HF-SCF open-shell capabilities.

We get converged ground and first excited states of 
several open-shell atoms and molecules.

We use .KPSELE for atoms.

For the sake of quick performance we use the STO-2G uncontracted basis sets and 
the 2-component (X2C+AMFI) Hamiltonian.

B 
--

^2P_1/2 (ground state)......2p_1/2(1)2p_3/2(0)
^2P_3/2 (excited state)......2p_1/2(0)2p_3/2(1)

C 
---

LS(low spin).....2p_1/2(2)2p_3/2(0)
HS(high spin)......2p_1/2(1)2p_3/2(1)

N
---

LS(low spin).....2p_1/2(2)2p_3/2(1)
HS(high spin)......2p_1/2(1)2p_3/2(2)

O
--

LS(low spin).....2p_1/2(2)2p_3/2(2)
HS(high spin)......2p_1/2(1)2p_3/2(3)

F
--

^2P_3/2 (ground state)......2p_1/2(2)2p_3/2(3)
^2P_1/2 (excited state)......2p_1/2(1)2p_3/2(4)

CH (^2Pi)
---------

We demonstrate the ^2Pi_1/2-^2Pi_3/2 spin-orbital splitting of the CH molecule. The excited state is
obtained either through the MOs reordering, or with the help of the mJ-occupation scheme.


CO(+)
-----

Starting from the DFCOEF of the neutral CO molecule ground state we are interested in the excited 4sigma-hole state of CO(+).
The calculation with the default SCF setting (OPENFAC=1) converges as desired to the excited 4sigma-hole state. 
However, the calculation with OPENFAC=0, converges to the 5sigma hole state.


