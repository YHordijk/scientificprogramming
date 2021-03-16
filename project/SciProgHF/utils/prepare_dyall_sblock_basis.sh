#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_edited


#
# 4s
#
block=4s

zeta=dz

val="correlating,4s4p polarizing"
core="correlating,4s4p polarizing"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                    --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                    --polarization_type="$val" \
                    --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                    --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                    --polarization_type="$core" \
                    --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


val="4s4p correlating,3s3p correlating,3s3p polarizing"
core="4s4p correlating,3s3p correlating,3s3p polarizing,2s2p correlating"

for zeta in tz qz; do
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas
done

#
# 5s
#
block=5s

zeta=dz

val="correlating,4s4p polarizing"
core="correlating,4s4p polarizing"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                    --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                    --polarization_type="$val" \
                    --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                    --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                    --polarization_type="$core" \
                    --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas


val="5s5p correlating,4s4p correlating,4s4p polarizing"
core="5s5p correlating,4s4p correlating,4s4p polarizing,3d correlating"

for zeta in tz qz; do
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done

#
# 6s
# 
block=6s

zeta=dz

val="correlating,5s5p polarizing"
core="correlating,5s5p polarizing"

python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                    --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                    --polarization_type="$val" \
                    --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas
python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                    --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                    --polarization_type="$core" \
                    --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas


val="6s6p correlating,5s5p correlating,5s5p polarizing"
core="6s6p correlating,5s5p correlating,5s5p polarizing,4d correlating"

for zeta in tz qz; do
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done

#
# 7s
#
block=7s

zeta=dz
val="correlating,4s4p polarizing"
core="correlating,4s4p polarizing"
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas

zeta=tz

val="7s7p correlating,6s6p correlating,6s6p polarizing"
core="7s7p correlating,6s6p correlating,6s6p polarizing,5d correlating"

for zeta in tz qz; do
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$val" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  
   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                       --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                       --polarization_type="$core" \
                       --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  
done
