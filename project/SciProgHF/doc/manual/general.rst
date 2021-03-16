:orphan:
 

starstar(GENERAL)


keyword(DIRECT)

Direct evaluation of two-electron integrals for Fock matrices (all two-electron
integrals for other uses, e.g. CI, CCSD, MCSCF, are always evaluated directly,
i.e. never read from disk).

The default is to evaluate LL, SL, and SS integrals directly (1 = evaluate
directly; 0 = do not evaluate directly)::

  .DIRECT
   1 1 1

keyword(SPHTRA)

Transformation to spherical harmonics embedded in the transformation to
orthonormal basis; totally symmetric contributions are deleted.

The default is a spherical transformation of large and small components,
respectively (1 = on; 0 = off)::

  .SPHTRA
   1 1

The transformation of the large components is a standard transformation from
Cartesian functions to spherical harmonics. The transformation of the small
components is modified, however, in accordance with the restricted kinetic
balance relation.


keyword(PCMOUT)

Write MO coefficients to the formatted file DFPCMO.  This is useful for porting
coefficients between machines with different binary structure. For reading the 
DFPCMO file, there is no keyword - simply copy this file to the working directory.

keyword(ACMOUT)

Write coefficients in C1 symmetry to the unformatted file DFACMO.


keyword(ACMOIN)

Import coefficients in C1 symmetry from the unformatted file DFACMO to current
symmetry. This assumes that the current symmetry is lower than the symmetry
used for obtaining the original coefficients.


keyword(LOWJACO)

Use Jacobi diagonalization in the Löwdin orthogonalization procedure
(subroutine LOWGEN). This is much slower than the default Householder 
method but does not mix AOs in the case of block-diagonal overlap matrix. 

keyword(DOJACO)

Use the Jacobi method for matrix diagonalization (currently limited to real
matrices). The default `Householder <http://www.maths.lancs.ac.uk/~gilbert/m306c/node21.html>`_  method is generally more efficient, but may mix
degenerate eigenvectors of different symmetries. 

keyword(QJACO)

Employ pure Jacobi diagonalization of quaternion matrixes. Properly handles
degenerate eigenvectors. Slower than .DOJACO and exclusive. Experimental
option.


keyword(LINDEP)

Thresholds for linear dependence in large and small components; refer to
the smallest acceptable values of eigenvalues of the overlap matrix.
The default is::

  .LINDEP
   1.0D-6 1.0D-8


keyword(RKBIMP)

Import SCF coefficients calculated using restricted kinetic balance (RKB) and
add the UKB component (all small component). This option is a nice way to get
(unrestricted) magnetic balance in response calculations involving a uniform
magnetic field (e.g. NMR shielding and magnetizability), in particular when
combined with London orbitals, which makes the magnetic balance *atomic*.


keyword(PRJTHR)

RKBIMP projects out the RKB coefficients transformed to orthonormal
basis and then adds the remained, corresponding to the UKB complement.
With the keyword you can set the threshold for projection.
The default is::

  .PRJTHR
  1.0D-5


keyword(PRINT)

General print level. The higher the number is, the more output the user
gets. Option of this type is useful for code debugging.
By default set to::

  .PRINT
   0


keyword(CVALUE)

Speed of light in a.u. Default::

  .CVALUE
   137.0359998


keyword(LOGMEM)

Write out a line in the output for each memory allocation done by DIRAC. This
is mainly useful for programmers or for solving out-of-memory issues.


keyword(CODATA)

Select the fundamental physical constants used througout the code.
The values are taken from different CODATA sets. Current available data sets are
CODATA86, CODATA98, CODATA02, CODATA06, CODATA10, CODATA14 and CODATA18.
Default::

  .CODATA
   CODATA18


keyword(NOSET)

.. warning:: documentation missing


keyword(FAMILY)

.. warning:: documentation missing


keyword(SKIP2E)

.. warning:: documentation missing

