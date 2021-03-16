#! /bin/bash
#
req=1.16762
cnv=$(echo "scale=10; a(1.0)/90.0" | bc -l )
for ang in 180.0 178.0 176.0 174.0 172.0
do
ca=$(echo "scale=10;  c($ang*$cnv)" | bc -l )
sa=$(echo "scale=10;  s($ang*$cnv)" | bc -l )
z1=$(echo "scale=5; $req*$sa" | bc -l )
x1=$(echo "scale=5; $req*$ca" | bc -l )
cat <<EOF > CO2.mol
INTGRL
Water (standard test example, now with symmetry recognition)  
cc-pVDZ basis (the contracted nrel basis)             
C   2    1  Y      A
        8.    2
O1     $x1         .0000000000        $z1
O2     $x1         .0000000000       -$z1
LARGE BASIS cc-pVDZ
        6.    1
C      .0000000000        0.0000000000        0.000
LARGE BASIS cc-pVDZ
FINISH
EOF
pam -noarch -nice -incmo -outcmo -run ang_${ang} CO2.mol CO2.inp
done
exit 0
