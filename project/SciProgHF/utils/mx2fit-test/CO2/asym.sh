#! /bin/bash
#
req=1.16762
for disp in 0.02 0.04 0.06 0.08 0.10
do
r1=$(echo "scale=5; $req+$disp" | bc -l )
r2=$(echo "scale=5; $req-$disp" | bc -l )
cat <<EOF > CO2.mol
INTGRL
Water (standard test example, now with symmetry recognition)  
cc-pVDZ basis (the contracted nrel basis)             
C   2    2  X  Y   A
        8.    2
O1     .0000000000         .0000000000        $r1
O2     .0000000000         .0000000000       -$r2
LARGE BASIS cc-pVDZ
        6.    1
C      .0000000000        0.0000000000        0.000
LARGE BASIS cc-pVDZ
FINISH
EOF
pam -noarch -nice -incmo -outcmo -run ${r1}_${r2} CO2.mol CO2.inp
done
exit 0
