!
!
!  MP1STPR:  Parameters for calculation of 1st order properties at
!            the MP2 level
!
!
      INTEGER LUMPPR
      PARAMETER (LUMPPR = 92)

      INTEGER NTMATR1, NTMATR2, IOFFTM1(2, 0:7, 2), IOFFTM2(2, 0:7, 2)
      COMMON /MP2PRP/ NTMATR1, NTMATR2, IOFFTM1, IOFFTM2

      CHARACTER*72 MPR_INDSTR(4, 2)
      COMMON /MP2CPR/ MPR_INDSTR
