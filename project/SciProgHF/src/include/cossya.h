!
      PARAMETER ( N_NFSYM   = 16 ) ! Number of fermion irrep's.
      PARAMETER ( N_NREP    = 32 ) ! Number of (fermion+boson)
!                                    irrep's.
      PARAMETER ( N_NOPEN   = 16 ) ! Number of GAS.
!
      LOGICAL KCOSINP,SYA_KTRDMAT
      COMMON /COSSYA0/ KCOSINP,SYA_KTRDMAT
!
      INTEGER SYA_NREP
      INTEGER SYA_PRINT
      INTEGER SYA_INACT(N_NFSYM)
      INTEGER SYA_ACTIV(N_NFSYM)
      INTEGER SYA_NOPEN
      INTEGER SYA_GASO(N_NREP,N_NOPEN)
      INTEGER SYA_IELC(N_NOPEN)
      INTEGER SYA_MAXE(N_NOPEN)
      INTEGER SYA_MINE(N_NOPEN)
      INTEGER SYA_NTRDM(N_NREP)
      INTEGER MSTATE
      COMMON /COSSYA1/ 
     &   SYA_PRINT, SYA_INACT, SYA_ACTIV, SYA_NOPEN, SYA_GASO, 
     &   SYA_IELC,  SYA_MAXE,  SYA_MINE,  SYA_NTRDM, SYA_NREP,
     &   MSTATE

