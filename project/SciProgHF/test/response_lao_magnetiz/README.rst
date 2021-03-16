Linear response test of magnetizabilities
=========================================

We calculate magnetizabilities (in a.u.) for the LiH molecule with 
London atomic orbitals (LAO) and without - in common gauge origin scheme (CGO).

For comparison of Levy-Leblond data we utilized the Dalton2011 progam suite (daltonprogram.org).

Small testing (uncontracted) basis is employed with point nucleus.

                      Isotropic         Anisotropy
-----------------------------------------------------------------
LL-LAO(GO=0)          -1.9809            -0.4603
LL-LAO(GO=1)          -1.9809            -0.4603
LL-LAO(GO=CM)         -1.9809            -0.4603

DalLAO(GO=CM)         -1.9809            -0.4603
DalLAO(GO=1)          -1.9809            -0.4603
DalLAO(GO=0)           

DalCGO(GO=CM)         -2.4712             0.2752
DalCGO(GO=1)          -68.6812           -68.1754 (1st),99.5966 (2nd)

LL-CGO(GO=1)          -68.6812           -68.1754 (1st),99.5966 (2nd)
LL-CGO(GO=CM)         -2.4712             0.2752

RKB-LAO(GO=CM)        -2.4839             0.2943
RKB-LAO(GO=1)

OpenRSP - LinRSP RKB-LAO(GO=CM) total magnetizabilities comparison:
London Magnetizability = Ebb w: -0.000  0.000    <-- OpenRSP
    -2.582033421    -0.000000000     0.000000000
    -0.000000000    -2.582033421     0.000000000
     0.000000000     0.000000000    -2.287732038
 LinRSP:
    -2.582032043279     -0.000000000000      0.000000000000
    -0.000000000000     -2.582032043279      0.000000000000
     0.000000000000      0.000000000000     -2.287732037645


RKB-CGO(GO=CM)        -2.4713             0.2754
RKB-CGO(GO=0)         -4.2371             2.9241
RKB-CGO(GO=1)

UKB-LAO(GO=1)         -1.7451            -0.8138
UKB-LAO(GO=CM)        -1.7451            -0.8138

UKB-CGO(GO=1)
UKB-CGO(GO=CM)

--------------------------------------------------------
Used positions of the gauge origin:
CM = center of mass of the LiH molecule
0 = (0.000,0.000,0.000)
1 = (11.000000,-4.500000,5.750000)


