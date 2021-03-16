:orphan:
 

.. _linking_to_math:

Linking to math libraries
=========================

General
-------

The DIRAC program uses libraries like BLAS or LAPACK for which a standard (but slow) 
implementation is found among the DIRAC source codes. 
To achieve optimal performance one should use optimized libraries at the linking stage 
and turn to the included versions (builtin math) only when no optimized version is available.

DIRAC requires BLAS and LAPACK libraries. Typically you will want to link to
external math (BLAS and LAPACK) libraries, for instance provided by MKL or
Atlas.

By default the CMake configuration script will automatically detect these libraries::

  $ ./setup --blas=auto --lapack=auto          # this is the default

if you define MATH_ROOT, for instance::

  $ export MATH_ROOT=/opt/intel/mkl

Do not use full path MATH_ROOT='/opt/intel/mkl/lib/ia32'. CMake will append the
correct paths depending on the processor and the default integer type.  If the
MKL libraries that you want to use reside in
/opt/intel/mkl/10.0.3.020/lib/em64t, then MATH_ROOT is defined as::

  $ export MATH_ROOT=/opt/intel/mkl/10.0.3.020

Then::

  $ ./setup [--flags]
  $ cd build
  $ make


Intel/MKL
---------

If you compile with Intel compilers and have the MKL library available, you
should use the --mkl flag which will automatically link to the MKL libraries
(in this case you do not have to set MATH_ROOT).
You have to specify whether you want to use the sequential or parallel
(threaded) MKL version. For a parallel DIRAC runs you should probably link to
the sequential MKL::

  $ ./setup --fc=mpif90 --cc=mpicc --cxx=mpicxx --mkl=sequential

For a sequential compilation you may want to link to the parallel MKL::

  $ ./setup --fc=ifort --cc=icc --cxx=icpc --mkl=parallel

The more general solution is to link to the parallel MKL and control the number
of threads using MKL environment variables.


Explicitly specifying BLAS and LAPACK libraries
-----------------------------------------------

If automatic detection of math libraries fails for whatever reason, you can
always call the libraries explicitly like here::

  $ ./setup --blas=/usr/lib/libblas.so --lapack=/usr/lib/liblapack.so

Alternatively you can use the --explicit-libs option. But in this case you should
disable BLAS/LAPACK detection::

  $ ./setup --blas=none --lapack=none --explicit-libs="-L/usr/lib -lblas -llapack"


Builtin BLAS and LAPACK implementation
--------------------------------------

If no external BLAS and LAPACK libraries are available, you can use the (slow) builtin
BLAS/LAPACK implementation. However note that these are not optimized and you will sacrifice
performance. This should be the last resort if nothing else is available::

  $ ./setup --blas=builtin --lapack=builtin


Linking to Atlas libraries on SUSE or Fedora
--------------------------------------------

DIRAC may not be able to detect Atlas BLAS/LAPACK
libraries on SUSE or Fedora because of nonstandard suffix::

  $ cd /usr/lib64/atlas
  $ ls -l

  lrwxrwxrwx. 1 root root   17 Oct 30 14:41 libatlas.so.3 -> ./libatlas.so.3.0
  -rwxr-xr-x. 1 root root 5.0M Sep  2  2011 libatlas.so.3.0
  lrwxrwxrwx. 1 root root   17 Oct 30 14:41 libcblas.so.3 -> ./libcblas.so.3.0
  -rwxr-xr-x. 1 root root 128K Sep  2  2011 libcblas.so.3.0
  lrwxrwxrwx. 1 root root   19 Oct 30 14:41 libclapack.so.3 -> ./libclapack.so.3.0
  -rwxr-xr-x. 1 root root  96K Sep  2  2011 libclapack.so.3.0
  lrwxrwxrwx. 1 root root   19 Oct 30 14:41 libf77blas.so.3 -> ./libf77blas.so.3.0
  -rwxr-xr-x. 1 root root 117K Sep  2  2011 libf77blas.so.3.0
  lrwxrwxrwx. 1 root root   18 Oct 30 14:41 liblapack.so.3 -> ./liblapack.so.3.0
  -rwxr-xr-x. 1 root root 5.4M Sep  2  2011 liblapack.so.3.0
  lrwxrwxrwx. 1 root root   19 Oct 30 14:41 libptcblas.so.3 -> ./libptcblas.so.3.0
  -rwxr-xr-x. 1 root root 128K Sep  2  2011 libptcblas.so.3.0
  lrwxrwxrwx. 1 root root   21 Oct 30 14:41 libptf77blas.so.3 -> ./libptf77blas.so.3.0
  -rwxr-xr-x. 1 root root 118K Sep  2  2011 libptf77blas.so.3.0

This can be fixed by creating soft links with .so suffix::

  $ su
  $ for file in *3; do ln -s $file ${file%.3}; done
  $ ls

  lrwxrwxrwx. 1 root root      13 Oct 30 17:45 libatlas.so -> libatlas.so.3
  lrwxrwxrwx. 1 root root      17 Oct 30 14:41 libatlas.so.3 -> ./libatlas.so.3.0
  -rwxr-xr-x. 1 root root 5225848 Sep  2  2011 libatlas.so.3.0
  lrwxrwxrwx. 1 root root      13 Oct 30 17:45 libcblas.so -> libcblas.so.3
  lrwxrwxrwx. 1 root root      17 Oct 30 14:41 libcblas.so.3 -> ./libcblas.so.3.0
  -rwxr-xr-x. 1 root root  130824 Sep  2  2011 libcblas.so.3.0
  lrwxrwxrwx. 1 root root      15 Oct 30 17:45 libclapack.so -> libclapack.so.3
  lrwxrwxrwx. 1 root root      19 Oct 30 14:41 libclapack.so.3 -> ./libclapack.so.3.0
  -rwxr-xr-x. 1 root root   97752 Sep  2  2011 libclapack.so.3.0
  lrwxrwxrwx. 1 root root      15 Oct 30 17:45 libf77blas.so -> libf77blas.so.3
  lrwxrwxrwx. 1 root root      19 Oct 30 14:41 libf77blas.so.3 -> ./libf77blas.so.3.0
  -rwxr-xr-x. 1 root root  119760 Sep  2  2011 libf77blas.so.3.0
  lrwxrwxrwx. 1 root root      14 Oct 30 17:45 liblapack.so -> liblapack.so.3
  lrwxrwxrwx. 1 root root      18 Oct 30 14:41 liblapack.so.3 -> ./liblapack.so.3.0
  -rwxr-xr-x. 1 root root 5566480 Sep  2  2011 liblapack.so.3.0
  lrwxrwxrwx. 1 root root      15 Oct 30 17:45 libptcblas.so -> libptcblas.so.3
  lrwxrwxrwx. 1 root root      19 Oct 30 14:41 libptcblas.so.3 -> ./libptcblas.so.3.0
  -rwxr-xr-x. 1 root root  130824 Sep  2  2011 libptcblas.so.3.0
  lrwxrwxrwx. 1 root root      17 Oct 30 17:45 libptf77blas.so -> libptf77blas.so.3
  lrwxrwxrwx. 1 root root      21 Oct 30 14:41 libptf77blas.so.3 -> ./libptf77blas.so.3.0
  -rwxr-xr-x. 1 root root  119824 Sep  2  2011 libptf77blas.so.3.0

Finally set::

  $ export MATH_ROOT=/usr/lib64/atlas

Now DIRAC will correctly detect them.


Linking to Atlas libraries on Linux Mint
----------------------------------------

DIRAC may not be able to detect Atlas BLAS/Lapack
libraries on Linux Mint because of some nonstandard (?) suffix.

After having installed them::

  $ sudo apt-get install liblapack3 libatlas3-base

the libraries were not detected by DIRAC. I had to do this::

  $ cd /usr/lib/atlas-base
  $ su
  $ for file in *3; do ln -s $file ${file%.3}; done

and this::

  $ cd /usr/lib
  $ sudo ln -s liblapack.so.3 liblapack.so

With this the automatic detection of math libraries works.


Linking to the ACML library
---------------------------

The Core Math Library (ACML) is distributed for most used Fortran compilers,
both for 4-byte and 8-byte integers, as dynamic and static, see the `ACML
web-page
<http://developer.amd.com/tools-and-sdks/cpu-development/amd-core-math-library-acml/>`_.

The user has to point the MATH_ROOT variable to the proper library directory.
For instance, here one wants to use 64-bit integer version of ACML, in
connection with Intel compilers::

  $ export MATH_ROOT='/home/milias/bin/acml5.3.1_ifort_int64/ifort64_int64/lib'

The complete variety of compiler- and integer-specific ACML's clones can be
found on the `dedicated web-page
<http://developer.amd.com/tools-and-sdks/cpu-development/amd-core-math-library-acml/acml-downloads-resources//>`_.

We report, however, that the gfortran generated ACML library can give these linking errors::

 Linking Fortran executable dirac.x
 /home/milias/bin/acml5.3.1_gnu_int64/gfortran64_int64/lib/libacml.so: undefined reference to `_gfortran_transfer_integer_write@GFORTRAN_1.4'
 /home/milias/bin/acml5.3.1_gnu_int64/gfortran64_int64/lib/libacml.so: undefined reference to `_gfortran_transfer_character_write@GFORTRAN_1.4'
 /home/milias/bin/acml5.3.1_gnu_int64/gfortran64_int64/lib/libacml.so: undefined reference to `_gfortran_transfer_real_write@GFORTRAN_1.4'
 collect2: ld returned 1 exit status

If this happens, we recommend user to prefer some non-gfortran ACML version.
