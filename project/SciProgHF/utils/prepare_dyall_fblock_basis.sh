#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_edited

block=5f
zeta=dz

val="5f correlating,6d/7s correlating,5f polarizing,6d polarizing"
core="5f correlating,6d/7s correlating,5f polarizing,6d polarizing,5d correlating"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=5f
zeta=tz

val="5f correlating,6d/7s correlating g,5f polarizing,6d polarizing"
core="5f correlating,6d/7s correlating g,5f polarizing,6d polarizing,5d correlating"

echo $val
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  



block=5f
zeta=qz


val="5f correlating,6d/7s correlating,6s6p correlating,5f polarizing,6d/7s polarizing"
core="5f correlating,6d/7s correlating,6s6p correlating,5f polarizing,6d/7s polarizing,5d correlating"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  



#
# 4f elements
#
block=4f
zeta=dz

val="4f correlating,5d/6s correlating,4f dipole polarization,5d dipole polarization"
core="4f correlating,5d/6s correlating,4f dipole polarization,5d dipole polarization,4d correlating"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas


python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas


block=4f
zeta=tz

val="4f correlating,5d/6s correlating,5s/5p correlating,4f dipole polarization"
core="4f correlating,5d/6s correlating,5s/5p correlating,4f dipole polarization,4d correlating"

echo $val
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas


python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas



block=4f
zeta=qz


val="4f correlating,5d/6s correlating,4f dipole polarization"
core="4f correlating,5d/6s correlating,4f dipole polarization,4d correlating"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas


python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas

