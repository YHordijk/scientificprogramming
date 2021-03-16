:orphan:
 
.. _exacc:

=====================================================
Highly parallelised relativistic Coupled Cluster Code
=====================================================

starstar(EXACC)

Highly parallelised relativistic Coupled Cluster Code. 
In this release only computations using the X2C Hamiltonian (with either :ref:`HAMILTONIAN_.X2Cmmf` or :ref:`HAMILTONIAN_.X2C`) are possible. 

This code is based on the math libraries `TAL-SH <https://github.com/DmitryLyakh/TAL_SH>`_ and `ExaTENSOR <https://github.com/ORNL-QCI/ExaTENSOR>`_ by Dmitry Lyakh. The tensors are kept in working memeory, sufficent RAM needs to be available. In order to test memory requirements instructions can be found in the ``exacorr_talsh_memory`` and ``exacorr_exatensor_memory`` tests. TAL-SH runs on a single node which has to have enough memory (:ref:`EXACC_.TALSH_BUFF`). In ExaTENSOR the memory is distributed, so each additional node will contribute its memory to the memory pool accessible by the library. Currently, it is recommended to use enough nodes that the tensors fit, but not substantially more. 

In the current release, if the library runs out of memory the code will not stop but enter a blocked state and calculations will not advance. So carefully control the advancement of your calculations, stopping them if they appear to hang.

**Mandatory keywords**
======================

keyword(OCCUPIED)

Defines occupied orbitals. Specification of list or energy range (see :ref:`orbital_strings`).

::

    .OCCUPIED
    energy -1.0 0.0 0.00001


keyword(VIRTUAL)

Defines virtual orbitals. Specification of list or energy range (see :ref:`orbital_strings`).

::

    .VIRTUAL
    20..30

**Optional keywords**
=====================

keyword(PRINT)

Print level.

*Default:*

::

    .PRINT
     0

keyword(TCONVERG)

Set convergence criteria (CC iterations, Lambda equations) 

*Default:*

::

    .TCONVERG
     1.0D-9

keyword(NCYCLES)

Maximum number of allowed CC iterations to solve the CC and LAMBDA equations. 

*Default:*

::

    .NCYCLES
     30

keyword(EXATENSOR)

This keyword activates the full multinode EXATENSOR library, which is designed for massively
parallel supercomputers. The additional infrastructure needed for parallel
communication makes this implementation inefficient when used for single node runs.
For such purposes the use of only the TALSH library component is recommended, which is designed for one node but
will make use of GPUs (if available and suitable).

*Default:*

::

    Do not use EXATENSOR


keyword(LAMBDA)

Solve Lambda-equations, needs to be activated in order to compute the one particle density matrix and molecular properties.

This calculation generates the file CCDENS, which contains the CC ground-state density
matrix in AO basis. In this release, CCDENS is used by the property module to calculate 
ground-state expectation values.

If saved, CCDENS can be used in a property calculation (see :ref:`PROPERTIES_.RDCCDENS`) 
without the need to invoke this module.


*Default:*

::

    Lambda equations are not solved

keyword(NOTRIPLES)

Deactivates computation of triples energy corrections (useful for ExaTENSOR as the current implementation is not efficient)

*Default:*

::

    Triples are done

keyword(CC2)

Performs a CC2 calculation instead of the default CCSD. Currently supported only for energies.

*Default:*

::

    CC2 is not activated


keyword(MOINT_SCHEME)

Expert option to choose another AO to MO integral transformation scheme. Change at your own risk. 

In TALSH only schemes 3 (default) and 42 (using Cholesky decompostion) are available. 

In ExaTensor schemes 1-4 and 42 are available with 42 using Cholesky decompostion. 
Scheme 4 is default for ExaTensor as it reduces the memory footprint by only keeping part of the AO integrals in memory.
The other methds keep all AO integrals in memeory. 
Scheme 0 prints the memory requirements and attempts to allocate the memory 
without doing the calculation. 

*Default:*

::

    .MOINT_SCHEME
     3

keyword(OCC_BETA)

Can be used to specify a "high-spin" reference determinant with a different number of "barred" occupied orbitals,
than "unbarred" occupied spinors. If .OCC_BETA is specified .OCCUPIED is interpreted as a list of unbarred (alpha) spinors.
NB: alpha and beta are used in a loose sense in relativistic calculations to indicate the (un)barred spinors.

::

    .OCC_BETA
    energy -1.0 0.0 0.00001

keyword(VIR_BETA)

Can be used to specify a different number of "barred" virtual orbitals than "unbarred" occupied spinors. 
If .VIR_BETA is specified .VIRTUAL is interpreted as a list of unbarred (alpha) spinors.
NB: alpha and beta are used in a loose sense in relativistic calculations to indicate the (un)barred spinors.

::

    .VIR_BETA
    20..30

keyword(CCDOUBLES)

Performs a CCD calculation instead of the default CCSD (switch off the contributions of single excitations).

*Default:*

::

    CCDOUBLES is not activated

keyword(EXA_BLOCKSIZE)

Expert option: Number to tune the parallel distribution (branching) of the spinor spaces.

*Default:*

::

    .EXA_BLOCKSIZE
     75

keyword(TALSH_BUFF)

Maximum memory (in gigabytes) used in TALSH, aim at about 80% of available memory on your machine.

*Default:*

::

    .TALSH_BUFF
     50

keyword(CHOLESKY)

Threshold to define the accuracy of the Cholesky decomposition (MOINT scheme 42), resulting
in inaccuracies of the computed energy of this order of magnitude (in Hartree units).

*Default:*

::

    .CHOLESKY
     1.0D-9

keyword(LSHIFT)

Expert option: Level shift of orbital energies, ignored for values smaller 0.

*Default:*

::

    .LSHIFT
     0.0D0



