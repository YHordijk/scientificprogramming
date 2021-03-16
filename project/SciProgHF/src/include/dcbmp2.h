!
! FILE    : ./include/dcbmp2.F
!
!  DCBMP2:  Parameters for MP2 calculation
!
!     MP2_NSTR(OCC/VIR, Fermion irrep)   - number of occ/virt pairs
!     MP2_INDSTR(OCC/VIR, Fermion irrep) - string of active pairs
!
      LOGICAL MP2_MSOUT, MP2_ANTIS, MP2NOP, SRTSHL, CMMNFO
      COMMON /DCLMP2/ MP2_MSOUT, MP2_ANTIS, MP2NOP, SRTSHL, CMMNFO

      CHARACTER*72 MP2_INDSTR(2, 2)
      COMMON /DCCMP2/ MP2_INDSTR

      REAL*8 DMP2_VIRTHR
      INTEGER IPRMP2, ISTRMP2, MP2_INTFLG, MP2_NSTR(2, 2)
      COMMON /DCBMP2/ DMP2_VIRTHR, IPRMP2, ISTRMP2, MP2_INTFLG,         &
     &                MP2_NSTR

      REAL*8 SCRMP2, EMP2, TPRI34, TPRI44
      COMMON /DCRMP2/ SCRMP2, EMP2, TPRI34, TPRI44
