:orphan:

Introduction
------------

The DIRAC program suite is parallelized using the MPI2 protocol.

See related works involving parallelization implementation, :cite:`Pernpointner2003`, :cite:`Saue1997`,
:cite:`Knecht2008`, :cite:`Knecht2010a`,
:cite:`Fleig2006`, :cite:`Fleig2001`, :cite:`Fleig2003`.
 

Parallel run information
------------------------

My shell environment variables are not forwarded to the compute nodes. What do
I have to do?

Assuming you are using the pam script to submit your parallel
calculation you can forward any environment variable, e.g., LD\_LIBRARY\_PATH,
to the slave-node shell by adding to you pam command line::

  ./pam --mpiarg="-x LD_LIBRARY_PATH"    # for OpenMPI
  ./pam --mpiarg="-envall"               # for IntelMPI and MPICH(2)

In the latter case ALL environment variables will be forwarded in one shot
whereas for **OpenMPI** you may need to provide the most essential environment
variables individually. These may be for a typical parallel DIRAC run::

  LD_LIBRARY_PATH
  BASDIR
  PATH
