#! /bin/bash
#
for r1 in 1.00 1.05 1.10 1.15 1.20 1.25 1.30 1.35 1.40
do
cat <<EOF > CO2.mol
INTGRL
Water (standard test example, now with symmetry recognition)  
cc-pVDZ basis (the contracted nrel basis)             
C   2    2  X  Y   A
        8.    2
O1     .0000000000         .0000000000        $r1
O2     .0000000000         .0000000000       -$r1
LARGE BASIS cc-pVDZ
        6.    1
C      .0000000000        0.0000000000        0.000
LARGE BASIS cc-pVDZ
FINISH
EOF
pam -noarch -nice -incmo -outcmo -run $r1 CO2.mol CO2.inp
done
exit 0
