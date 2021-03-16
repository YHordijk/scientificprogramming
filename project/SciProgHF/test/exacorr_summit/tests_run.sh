#!/bin/bash

# script to run tests on summit
echo "Don't forget to keep ~/.diracrc up to date"

# path for dirac.x and exacorr.x
export DIRAC_PATH=/gpfs/alpine/world-shared/chp109/dirac_install

if [ $# -eq 0 ]
then
  echo "give the number of the test you want to run"
else

pam_script=job_pam.pbs

export INP_file=cc_uf6.inp
export MOL_file=UF6.xyz
if [ "$1" -eq 1 ] || [ "$1" -gt 8 ]
then
  fold=scheme1
  cp ${pam_script} ${fold}
  cd ${fold}
    bsub ${pam_script}
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 2 ] || [ "$1" -gt 8 ]
then
  fold=scheme2
  cp ${pam_script} ${fold}
  cd ${fold}
    bsub ${pam_script}
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 3 ] || [ "$1" -gt 8 ]
then
  fold=scheme3
  cp ${pam_script} ${fold}
  cd ${fold}
    bsub ${pam_script}
  cd ..
  echo "${fold} submited"
fi

if [ "$1" -eq 4 ] || [ "$1" -gt 8 ]
then
  fold=scheme4
  cp ${pam_script} ${fold}
  cd ${fold}
    bsub ${pam_script}
  cd ..
  echo "${fold} submited"
fi

export INP_file=ccd_uf6.inp
if [ "$1" -eq 5 ] || [ "$1" -gt 8 ]
then 
  fold=ccd
  cp ${pam_script} ${fold}
  cd ${fold}
    bsub ${pam_script} 
  cd ..
  echo "${fold} submited"
fi

export MOL_file=CHFClBr.xyz
export INP_file=exacc_chfclbr.inp
if [ "$1" -eq 6 ] || [ "$1" -gt 8 ]
then
  fold=chfclbr
  cp ${pam_script} ${fold}
  cd ${fold}
    bsub ${pam_script}
  cd ..
  echo "${fold} submited"
fi

#tests not using pam

if [ "$1" -eq 7 ] || [ "$1" -gt 8 ]
then
  fold=script_exatensor
  cd ${fold}
    bsub job_exatensor.pbs
  cd .. 
  echo "${fold} submited"
fi

if [ "$1" -eq 8 ] || [ "$1" -gt 8 ]
then
  fold=script_talsh
  cd ${fold}
    bsub job_talsh.pbs
  cd ..
  echo "${fold} submited"
fi

fi
