W_S operator in DIRAC
======================

The operator is added as finite-field perturbation at the KU-RelCC level.

This operator contains eitehr PVC of FC representation of the nuclear density.

Test system: MgF, STO-2G uncontracted basis, all (13) correlated electrons


field  symm        SCF                        CCSD(T)
------------------------------------------------------------
0.0    C2v    -287.76439849222800      -288.008945929251695
0.0    C2     -287.76439849222265      -288.008945929251695
              -287.76439849222265      -288.008945929251695
0.0    Cs     -287.76439849222726      -288.008928821452059
0.0    C1     -287.76439849220532      -288.008944099824760 <--- !


1e-4   C2v    -287.76439849222800      -288.008945385656489 (Ws-PVC)
                                       -288.008945385656205 (Ws-FC)

1e-4   C2     -287.76439849222265      -288.008945385656602 (Ws-PVC)
                                       -288.008945385656318 (Ws-FC)

1e-4   Cs     -287.76439849222726      -288.008928821452059 <--- !
1e-4   C1     -287.76439849220532      -288.008943556222960 <--- !



use orbital energies (USEOE) :
0.0    C2v                             -288.00948672977131
0.0    C2                              -288.009486729771766
0.0    Cs                              -288.009486729581226 <---- !
0.0    C1                              -288.009486729771425

1e-4   C2v                           
1e-4   C2                             
1e-4   Cs                             
1e-4   C1                             

