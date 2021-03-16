#!/bin/bash

cp DFCOEF_sym DFCOEF

export DIRAC_MPI_COMMAND="mpiexec -n 2 -hosts localhost"

/Users/johannpototschnig/Documents/DIRAC/frozen_open_shell/build-omp/pam --inp=cc_ru.inp --mol=Ru.mol --put="DFCOEF" --get="DFCOEF" --mw=600

mv DFCOEF DFCOEF_sym
