#!/bin/bash

if [ $# -eq 0 ]
then
  echo "give the number of the test you want to check"
else

INP_file=cc_uf6.inp
MOL_file=UF6.xyz

if [ "$1" -eq 1 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme1
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 2 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme2
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 3 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme3
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

if [ "$1" -eq 4 ] || [ "$1" -gt 8 ]
then
  FOLD=scheme4
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

INP_file=ccd_uf6.inp
if [ "$1" -eq 5 ] || [ "$1" -gt 8 ]
then
  FOLD=ccd
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

MOL_file=CHFClBr.xyz
INP_file=exacc_chfclbr.inp
if [ "$1" -eq 6 ] || [ "$1" -gt 8 ]
then
  FOLD=chfclbr
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

INP_file=cc_laf3.inp
MOL_file=LaF3.xyz
if [ "$1" -eq 6 ] || [ "$1" -gt 8 ]
then
  FOLD=script_exatensor
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

INP_file=cc_sf6.inp
MOL_file=SF6.xyz
if [ "$1" -eq 7 ] || [ "$1" -gt 8 ]
then
  FOLD=script_talsh
  echo ""
  echo " compare results in ${FOLD}  "
  echo "======================================================"
  cd $FOLD
    fname="${INP_file%.*}_${MOL_file%.*}.out"
    grep -n -A 8 -a "Final results from EXACORR" ${fname} | tail -n 5
    echo "-----------------------------"
    grep -n -A 8 -a "Final results from EXACORR" ../result/${fname} | tail -n 5
    echo ""
    grep -n -A 1 -a "Final L1amp" ${fname}
    echo "-----------------------------"
    grep -n -A 1 -a "Final L1amp" ../result/${fname}
    echo ""
    grep -n -A 2 -a "Dipole length: X" ${fname} | tail -n 3
    echo "-----------------------------"
    grep -n -A 2 -a "Dipole length: X" ../result/${fname} | tail -n 3
  cd ..
  echo "======================================================"
fi

fi
