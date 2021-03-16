#!/bin/bash
for dft in LDA BLYP B3LYP CAMB3LYP PBE PBE0 BP86
do
cat <<EOF > H2O_dft.inp
**DIRAC
.TITLE
H2O
.WAVE FUNCTIONS
.ANALYZE
**HAMILTONIAN
.X2C
.DFT
${dft}
**INTEGRALS
*TWOINT
.SCREEN
1.0D-12
*READIN
.UNCONT
**WAVE FUNCTIONS
.SCF
*SCF
.CLOSED SHELL
10
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
*END OF
EOF
pam --inp=H2O_dft  --mol=H2O --outcmo --suffix=${dft}
for occ in 1.0 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9
do
cat <<EOF > H2O_1s_dft.inp
**DIRAC
.TITLE
H2O
.WAVE FUNCTIONS
.ANALYZE
**HAMILTONIAN
.X2C
.DFT
${dft}
**INTEGRALS
*TWOINT
.SCREEN
1.0D-12
*READIN
.UNCONT
**WAVE FUNCTIONS
.SCF
.REORDER
2..5,1
*SCF
.CLOSED SHELL
8
.OPEN SHELL
1
${occ}/2
.OVLSEL
.NODYNSEL
.FOCC
**ANALYZE
.MULPOP
*MULPOP
.VECPOP
1..oo
*END OF
EOF
pam --inp=H2O_1s  --mol=H2O --incmo --suffix=${dft}_occ${occ}
done
done
exit 0
