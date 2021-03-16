C
C     Small common block to pass on information used in the
C     intermediate Hamiltonian formalism.
C
C     For the Intermediate Hamiltonian method we need to identify spinors as either
C     element of Pm or Pi space. We use an array that gives this info that has
C     the same structure as the EPS array (that is passed as argument).
      INTEGER MAXSPIM
      PARAMETER (MAXSPIM=5000)
      REAL*8 CIH,ER_IH,AIH,SHIFT_IH
      COMMON /CCIHI/ IPIORB(MAXSPIM),IHSCHEME,NIH
      INTEGER IPIORB,IHSCHEME,NIH
      COMMON /CCIHR/ CIH,ER_IH(2,2),AIH,SHIFT_IH(2,6)
      LOGICAL DOIH
      COMMON /CCIHL/ DOIH
