:orphan:
 

Basic installation
==================

General
-------

DIRAC is configured using `CMake <http://www.cmake.org/>`_ , typically via the ``setup`` script,
and subsequently compiled using make (or gmake).
The ``setup`` script is a useful front-end to CMake.
You need python to run ``setup``. To see all options, run::

  $ ./setup --help


Sequential build
----------------

The default installation proceeds through three commands::

  $ ./setup [--flags]
  $ cd build
  $ make

The ``setup`` script creates the directory "build" and
calls CMake with appropriate environment variables and flags.
By default CMake builds out of source. This means that all object files and the
final binary are generated outside of the source directory. Typically the build
directory is called "build", but you can change the name of the build directory
(e.g. "build_gfortran")::

  $ ./setup [--flags] build_gfortran
  $ cd build_gfortran
  $ make

You can compile the code on all available cores::

  $ make -j

We strongly recommend that the installation is followed by testing your installation::

  make test

For more information about testing, see `here <testing.html>`_.


Parallel build
--------------

Building DIRAC for parallel runs using `MPI <http://www.mpi-forum.org/docs/docs.html>`_ 
is in principle a very simple modification of the above procedure::

  $ ./setup [--flags] --mpi
  $ cd build
  $ make

Once again we recommend `testing the code <testing.html>`_ afterwards.


Typical examples
----------------

In order to get familiar with the configuration setup, let us demonstrate
some typical configuration scenarios.

Configure for parallel compilation using MPI (make sure to properly export MPI
paths)::

  $ ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx

There is a shortcut for it::

  $ ./setup --mpi

Configure for sequential compilation using ifort/icc/icpc and link against parallel mkl::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --mkl=parallel

Configure for sequential compilation using gfortran/gcc/g++::

  $ ./setup --fc=gfortran --cc=gcc --cxx=g++

You get the idea. The configuration is usually good at detecting math libraries
automatically, provided you export the proper environment variable ``MATH_ROOT``,
see :ref:`linking_to_math`.


What to do if CMake is not available or too old?
------------------------------------------------

If it is your machine and you have an Ubuntu or Debian-based distribution::

  $ sudo apt-get install cmake

On Fedora::

  $ sudo yum install cmake

Similar mechanisms exist for other distributions or
operating systems. Please consult Google.

If it is a cluster, please ask the Administrator to install/upgrade CMake.

If it is a cluster, but you prefer to install it yourself (it's easy):

1. `Download <http://www.cmake.org/cmake/resources/software.html>`_ the latest pre-compiled CMake tarball
2. Extract the tarball
3. Set correct PATH variable
