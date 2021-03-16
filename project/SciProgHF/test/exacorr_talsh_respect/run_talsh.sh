#!/bin/bash

OLD_PATH=$(pwd)

echo $OLD_PATH

mkdir /tmp/respect_try
cd /tmp/respect_try

#ExaTensor specific:
export QF_PATH=~/DIRAC/caar/build-omp              #full path to ExaTENSOR root directory
export QF_NUM_PROCS=1             #total number of MPI processes
export QF_PROCS_PER_NODE=1        #number of MPI processes per logical node (logical nodes are created by node resource isolation)
export QF_CORES_PER_PROCESS=1     #number of physical CPU cores per MPI process (no less than 1)
export QF_MEM_PER_PROCESS=1576    #host RAM memory limit per MPI process in MB
export QF_NVMEM_PER_PROCESS=0     #non-volatile memory limit per MPI process in MB
export QF_HOST_BUFFER_SIZE=1024   #host buffer size per MPI process in MB (must be less than QF_MEM_PER_PROCESS)
export QF_GPUS_PER_PROCESS=0      #number of discrete NVIDIA GPU's per MPI process (optional)
export QF_MICS_PER_PROCESS=0      #number of discrete Intel Xeon Phi's per MPI process (optional)
export QF_AMDS_PER_PROCESS=0      #number of discrete AMD GPU's per MPI process (optional)
export QF_NUM_THREADS=4           #initial number of CPU threads per MPI process (irrelevant, keep it 8)

#OpenMP generic:
export OMP_NUM_THREADS=$QF_NUM_THREADS #initial number of OpenMP threads per MPI process
export OMP_DYNAMIC=false               #no OpenMP dynamic threading
export OMP_NESTED=true                 #OpenMP nested parallelism is mandatory
export OMP_MAX_ACTIVE_LEVELS=3         #max number of OpenMP nesting levels (at least 3)
export OMP_THREAD_LIMIT=256            #max total number of OpenMP threads per process
export OMP_WAIT_POLICY=PASSIVE         #idle thread behavior
#export OMP_STACKSIZE=200M             #stack size per thread
#export OMP_DISPLAY_ENV=VERBOSE        #display OpenMP environment variables
#export GOMP_DEBUG=1                   #GNU OpenMP debugging
#export LOMP_DEBUG=1                   #IBM XL OpenMP debugging

#OpenMP thread binding:
export OMP_PLACES_DEFAULT=threads                                      #default thread binding to CPU logical cores
export OMP_PROC_BIND="close,spread,spread" #nest1: Functional threads (DSVU)
                                           #nest2: TAVP-WRK:Dispatcher spawns coarse-grain Executors
                                           #nest3: TAVP-WRK:Dispatcher:Executor spawns execution threads


rm core.* *.tmp *.log *.x RSD_MOS *.out

cp $QF_PATH/exacorr.x .
cp $OLD_PATH/exacorr.inp .
cp $OLD_PATH/RSD_MOS .

./exacorr.x


cd $OLDPATH

#ulimit -s unlimited


