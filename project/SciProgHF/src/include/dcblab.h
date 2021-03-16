!*** Labels ***
!    PLABEL(*,1), IATTR(*,1), IPLAB(*,1) for AO basis
!    PLABEL(*,2), IATTR(*,2), IPLAB(*,2) for SO basis
!    NPLAB(0:2): number of labels total/large comp/small comp
      INTEGER MAXLAB
      PARAMETER (MAXLAB = 12000)

      CHARACTER*12 PLABEL(MAXLAB, 2)
      COMMON /LABNAM/ PLABEL

      INTEGER IATTR(MAXLAB, 2), IPLAB(MXCORB, 2), NPLAB(0:2)
      COMMON /LABELS/ IATTR, IPLAB, NPLAB
