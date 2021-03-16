:orphan:
 

Installation on `HPC cluster <http://services-numeriques.unistra.fr/hpc.html>`_ in Strasbourg
=============================================================================================

Configuration and compilation
-----------------------------

In order to find the proper software needed for configuration and compilation the following line is 
added to `.bashrc` ::

  source $HOME/modules

where the file `modules` read (as of June 2013)::

  module delete mpi/openmpi-1.4.i11
  module delete compilers/intel11
  module delete libs/mkl11 
  module load mpi/openmpi-1.6.i13
  module load compilers/intel13
  module load libs/mkl13 

To properly set environmental variables for the Intel MKL math library  we follow instructions 
`here <../mkl.html>`_ and add to `.bashrc` the following lines::

  export MKL_NUM_THREADS=1
  export MKL_DYNAMIC="FALSE"
  export OMP_NUM_THREADS=1

Finally we configure and build the binary::

  $ ./setup --fc=mpif90 --cc=mpicc
  $ cd build
  $ make

Setting up pam
--------------

In your top directory you may configure the pam run script by adding a file ``.diracrc`` with contents::

  --scratch=/scratch  
  --noarch
  --debugger=/usr/bin/gdb
  --basis="$HOME/Dirac/basis:$HOME/Dirac/basis_dalton"
  --profiler=/usr/bin/gprof

Example run script
------------------

Here is a simple example of a run script::

  #! /bin/bash
  #SBATCH -p pri2008
  #SBATCH -A qossaue
  #SBATCH -t 00:10:00
  #SBATCH -n 24   # xx Coeurs 
  #SBATCH --mail-type=END
  #SBATCH --mail-user=trond.saue@irsamc.ups-tlse.fr
  . ~/.bashrc
  scontrol show hostname > machinelist
  pam --inp=HF --mol=H2O --mpi=$SLURM_NPROCS --machfile=machinelist
  exit 0


