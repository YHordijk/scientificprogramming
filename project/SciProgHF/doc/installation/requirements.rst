:orphan:


Requirements
============

DIRAC requires CMake, Fortran 90, C, and C++ compilers and a Python environment
for the setup and run scripts. The program is designed to run on a unix-like
operating system, and uses MPI for parallel calculations.  Optionally DIRAC can
(and should) be linked to platform specific BLAS and LAPACK libraries.


Minimum requirements
--------------------

* CMake 3.0 - 3.11
* Python 2.6 - 2.7
* GFortran/GCC/G++ 4.9 - 5.2 or Intel compilers 13.0 - 15.0.6

If your distribution does not provide a recent CMake, you can install it
without root permissions within few minutes. Download the binary distribution
from https://cmake.org/download/, extract it, and set the environment variables
PATH and LD_LIBRARY_PATH to point to the bin/ and lib/ subdirectories of the
extracted path, respectively. Verify that your new CMake is recognized with
``cmake --version``.


Requirements for building with the parallel library ExaTensor
-------------------------------------------------------------

* CMake 3.0 - or later
* Python 3.6 - or later
* GFortran/GCC/G++ 8 (9 has a bug and cannot be used) or Intel compilers 18 or later

At runtime you need to set environment variables to inform the library about the optimal 
parallel setup in both the installation as well as the execution stage. This is explained 
in detail at the `ExaTensor website <https://github.com/ORNL-QCI/ExaTENSOR>`_. 
The setup script of DIRAC should usually take care of the installation step, but the appropriate set up at
runtime should be carefully considered by the user to obtain optimal performance. Examples of setup 
for different architectures are given `here <http://github.com/ORNL-QCI/ExaTENSOR/blob/master/run.sh>`_.

Supported, but optional dependencies
------------------------------------

* OpenMPI 1.6.2 - 1.8.5 (MPI parallelization)
* IntelMPI 4.1 - 5.0 (MPI parallelization)
* Boost 1.54.0 - 1.60.0 (PCMSolver; installed along DIRAC if below 1.54.0)
* HDF5 1.8.14 - 1.8.16 (DMRG code)
* Alps 2.2 (DMRG code)


Additional dependencies due to PCMSolver
----------------------------------------

It is optionally possible to enable the polarizable continuum model functionality.
The functionality is provided using the external module
`PCMSolver <http://pcmsolver.github.io/pcmsolver-doc/>`_. The module requires additional
dependencies:

* the zlib compression library, version 1.2 or higher;
* the Boost set of libraries, version 1.54 or higher.

The module performs checks on the system to find an installation of Boost that has
the suitable version. The search is directed towards conventional system directories.
If the search fails, the module will build a copy of Boost on the spot.
If your system has Boost installed in a nonstandard location, e.g. /boost/location, you have to direct the
search by adding the following flags to the ``setup`` script::

  $ ./setup -DBOOST_INCLUDEDIR=/boost/location/include -DBOOST_LIBRARYDIR=/boost/location/lib

By default, compilation of the module is enabled. CMake will check that the additional
dependencies are met. In case they are not, compilation of PCMSolver will be disabled.
CMake will print out which dependencies were not satisfied.
Passing the option ``-DENABLE_PCMSOLVER=OFF`` to the ``setup`` script will disable compilation
of the module and skip detection of additional dependencies.
