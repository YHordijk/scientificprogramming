#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_edited


zeta=dz
block=4d
val="valence correlating and dipole polarization"

core="valence correlating and dipole polarization,3d correlating"
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=5d

core="valence correlating and dipole polarization,4f core correlating"
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=6d

core="valence correlating and dipole polarization,5f correlating"
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


#
block=4d
val="valence correlating,4d dipole polarization,4s4p correlating"
core="valence correlating,4d dipole polarization,4s4p correlating,3d correlating"
for zeta in tz qz; do
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done

block=6d
val="valence correlating,6d dipole polarization"
core="valence correlating,6d dipole polarization,5f correlating"
for zeta in tz qz; do
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done

block=5d
val="valence correlating,valence dipole polarization"
core="valence correlating,valence dipole polarization,4f core correlating"
zeta=tz
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


block=5d
val="valence correlating,valence dipole polarization"
core="valence correlating,valence dipole polarization,outer core (5s5p) correlating g,4f core correlating'"
zeta=qz
echo "processing "$block" "$basis" "$zeta" "$val" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
echo "processing "$block" "$basis" "$zeta" "$core" set"
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

