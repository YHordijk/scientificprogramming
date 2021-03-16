C
C     dcbstex.h -- information for STEX model in prp/
C                  (Ulf Ekstroem, 2005-2006)
C
      INTEGER STEX_MAX_HOLES, STEX_NR_OP
      PARAMETER (STEX_MAX_HOLES = 10,STEX_NR_OP = 6)
c
C     STEX_NHOLES - total number of holes in all irreps.
c     STEX_HOLES(I) - orbital index of orbital I
c     STEX_NCUBE - Nr of virtuals to make cube files from.
C     STEX_IFSYM(I) - fermion symmetry of hole I
C     STEX_NO2E - Neglect coulomb and exchange part of the Hamiltonian.
C     STEX_DIAG - Neglect all off-diagonal elements of the Hamiltonian
C     STEX_NOCOUPLING - Neglect the coupling between excitations from
C                       different orbitals
C     STEX_COFACTOR - Use a separate reference state
C     STEX_BATCHSIZE - Nr of simultaneous fock matrices
C                      to get from TWOFCK
C      
      INTEGER STEX_NHOLES,STEX_HOLES,
     &     STEX_PRINT,STEX_OP_SYM,STEX_OP_TIM,
     &     STEX_OP_FER,STEX_OP_IND,STEX_NVIR,STEX_NVIRT,
     &     STEX_I2VIRX, STEX_N2VIRX, STEX_N2VIRXQ,
     &     STEX_IFSYM, STEX_NPARAM,
     &     STEX_LUINT,STEX_BATCHSIZE,STEX_NCUBE

      LOGICAL STEX_COFACTOR,STEX_NO2E,STEX_DIAG,STEX_NOCOUPLING

      DOUBLE PRECISION STEX_CUTOFF, STEX_SCRFCK

      COMMON /CBISTEX/ STEX_NHOLES,STEX_HOLES(STEX_MAX_HOLES),
     &     STEX_PRINT,STEX_OP_SYM(6),STEX_OP_TIM(6),
     &     STEX_OP_FER(6),STEX_OP_IND(6), STEX_NVIR(2),STEX_NVIRT,
     &     STEX_I2VIRX(2,2), STEX_N2VIRX, STEX_N2VIRXQ,
     &     STEX_IFSYM(STEX_MAX_HOLES), STEX_NPARAM,
     &     STEX_LUINT,STEX_BATCHSIZE,STEX_NCUBE
      COMMON /CBLSTEX/ STEX_COFACTOR,STEX_NO2E,STEX_DIAG,
     &     STEX_NOCOUPLING
      COMMON /CBDSTEX/  STEX_CUTOFF, STEX_SCRFCK
