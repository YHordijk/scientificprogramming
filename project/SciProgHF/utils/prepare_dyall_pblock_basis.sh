#!/bin/sh

#
# $Id$
#
# shell script to prepare intermediate files for assembling dyall basis dirac/dalton basis set files
#

ARCHIVE_DIR=../basis/reference_edited

for block in  4p 5p 6p 7p; do

#note, for the double-zeta, add the polarizing f function by hand on the core.bas files?
   zeta=dz

   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  

   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating','polarizing f' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  


   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating','diffuse' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  

   python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                          --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                          --polarization_type='correlating','diffuse','polarizing f' \
                          --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  

   
   for zeta in tz qz; do

      if [ $block != "7p" ]; then
         python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence correlating' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  


         python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence correlating','core correlating' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

         python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence correlating','diffuse' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  


         python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence correlating','core correlating','diffuse' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  


      else

         python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence 7s7p correlating' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.val.bas  

         python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence 7s7p correlating','diffuse' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-val.bas  

         if [ $zeta == "tz" ] ; then

            python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence 7s7p correlating','outer core 6d correlating' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

            python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                   --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                   --polarization_type='valence 7s7p correlating','outer core 6d correlating','diffuse' \
                                   --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  
         else

            python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                --polarization_type='valence 7s7p correlating','core 6d correlating' \
                                --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.core.bas  

            python basis_set.py --file=$ARCHIVE_DIR/$block"_"$zeta"_archive" \
                                   --print_polarization --print_scfbas --uncontracted --equal_exponent_threshold=0.05 \
                                   --polarization_type='valence 7s7p correlating','core 6d correlating','diffuse' \
                                   --output=dyall.$block.$zeta.bas > dyall.$block.$zeta.aug-core.bas  
         fi
      fi
   done
done
