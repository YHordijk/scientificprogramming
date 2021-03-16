W_A operators in DIRAC
======================

The W_A operator of given component (x/y/z-alpha) is added via finite-field perturbation at the KU-RelCC level.

The Wa operator preferably contains PVC representation of the nuclear density.

The test systems are MgF and MgOH molecules in STO-2G uncontracted basis and with all (13) correlated electrons

The "+" field means perturbed energy is higher in absolute value than the unperturbed energy.
The "-" means lower perturbed energy in absolute value than the unperturbed energy.


MgF molecule
------------

For this linear molecule, the operator z-component is zero based on the theory. 
The x and y operator components are identical (x/-y and -x/y) and give nonzero contributions 
for linear systems.
The "+" and "-" components means different (higher and lower wrt absolute value) energies.


oper/ff  symm          CCSD energy           CCSD(T) energy       comment
------------------------------------------------------------------------------------
 n/a     C2v_x_y  -287.998691032182478   -288.008945929244419
z/1e-6   C2v_x_y  -287.998691032182478   -288.008945929244419
z/-1e-6  C2v_x_y  -287.998691032182478   -288.008945929244419 (no impact of z oper)

 n/a     Cs_x     -287.998691025269181   -288.008928821444727
 n/a     Cs_y     -287.998691025269977   -288.008928821445011 (no operator)

z/1e-6   Cs_x     -287.998691025269295   -288.008928821444840 
z/1e-6   Cs_y     -287.998691025270148   -288.008928821445181 
z/-1e-6  Cs_y     -287.998691025270148   -288.008928821445181 (no impact of z oper)

x/1e-6   Cs_x     -287.998691025269181   -288.008928821444727 
y/1e-6   Cs_y     -287.998691025269977   -288.008928821445011 (no impact,symm.cancelation)

y/1e-6   Cs_x     -287.998691019658622   -288.008928816000036 (impact of +y oper on Cs_x,"-" field)
x/-1e-6  Cs_y     -287.998691019659418   -288.008928816000321 (impact of -x oper on Cs_y,"-" field)

x/1e-6   Cs_y     -287.998691030880821   -288.008928826889985 (impact of +x oper on Cs_y,"+" field)
y/-1e-6  Cs_x     -287.998691030879968   -288.008928826889587 (impact of -y oper on Cs_x,"+" field)


MgOH molecule
-------------

As the next example we have the bent MgOH molecule. Clearly, z operator component has no influence 
on resulting CC energy.

The x and y operator components give identical (x/-y and -x/y) nonzero contributions, 
where "+" and "-" components give different (higher and lower) energies.

In the C1 symmetry, however, one does not get "+" and "-"distinguished energies.

oper/ff  symm      CCSD energy           CCSD(T) energy       comment
--------------------------------------------------------------------------------
 n/a     Cs_x   -265.010898953114747   -265.018875137512737  (no operator)
z/+1e-6  Cs_x   -265.010898953114861   -265.018875137512850  (no impact of +z comp. on Cs_x)
z/-1e-6  Cs_x   -265.010898953114861   -265.018875137512850  (no impact of -z comp. on Cs_x)
z/+1e-6  Cs_y   -265.010898953114520   -265.018875137512794  (no impact of +z comp. on Cs_y)

y/+1e-6  Cs_y   -265.010898953114406   -265.018875137512680  (no impact of y comp on Cs_y,symm.cancelation)
x/+1e-6  Cs_x   -265.010898953114747   -265.018875137512737  (no impact of x comp on Cs_x,symm.cancelation)

y/+1e-6  Cs_x   -265.010898948209785   -265.018875132702249  (impact of +y oper on Cs_x,"-" field)
x/-1e-6  Cs_y   -265.010898948209388   -265.018875132702135  (impact of -x oper on Cs_y,"-" field)

x/+1e-6  Cs_y   -265.010898958019652   -265.018875142323452  (impact of +x oper on Cs_y,"+" field)
y/-1e-6  Cs_x   -265.010898958020050   -265.018875142323623  (impact of -y oper on Cs_x,"+" field)

n/a      C1     -265.010898951241131   -265.018885665647588  (no operator)
z/1e-6   C1     -265.010898951241245   -265.018885665647701  (no impact of +z comp. on C1)
y/1e-6   C1     -265.010898947207977   -265.019683350360083  (impact of +y oper on C1,"-" field,CCSD)
x/-1e-6  C1     -265.010898947208034   -265.019683327168991  (impact of -x oper on C1,"-" field,CCSD)

y/-1e-6  C1     -265.010898947208034   -265.019683343151314  (impact of -y oper on C1,"-" field,CCSD)
x/+1e-6  C1     -265.010898947208034   -265.019683336227217  (impact of +x oper on C1,"-" field,CCSD)

