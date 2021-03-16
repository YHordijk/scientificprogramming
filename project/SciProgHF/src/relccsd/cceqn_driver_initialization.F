      SUBROUTINE cceqn_driver_init(MAXIT,FVO,T1,S1,T2,S2,
     &                 BUF1,BUF2,BUF3,CC1,BB1,BB2,
     &                 NBUF2,NBUF3,RMS,FOCKSP,DELS1,DELS2)

                             
      implicit none
C
C     Intermediaries etc.

      REAL*8 FVO(*),T1(*),S1(*),T2(*),S2(*)
      REAL*8 BUF1(*),BUF2(*),BUF3(*)

C     for diis ...

      REAL*8 CC1(*)
      REAL*8 BB1(*),BB2(*)
      REAL*8 RMS
      REAL*8 DIAG

      integer maxit,nbuf2,nbuf3

C---------------Common Blocks--------------------------------------
#include "files.inc"
#include "param.inc"
#include "symm.inc"
#include "eqns.inc"
#include "inpt.inc"
#include "ihm.inc"
#include "results.inc"
#include "complex.inc"
#include "ccpar.inc"
!#if defined (VAR_MPI)
!      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
!#endif
C---------------Local variables--------------------------------------
      LOGICAL TEQ,FOCKSP      
      CHARACTER*2 DELS1
      CHARACTER*4 DELS2
      integer it,itc,nit
C---------------Executable code--------------------------------------
C
C
      IF (EQNS.EQ.'AMPLTD') THEN
         TEQ = .TRUE.            
         FOCKSP = .FALSE.
      ELSEIF (EQNS.EQ.'LAMBDA') THEN
         TEQ = .FALSE.            
         FOCKSP = .FALSE.
      ELSEIF (EQNS.EQ.'FOCKSP') THEN
         TEQ = .TRUE.            
         FOCKSP = .TRUE.
         DELS1 = 'DD'
         IF (ICURFSS.EQ. 2) DELS1 = 'DK'
         IF (ICURFSS.EQ. 3) DELS1 = 'KD'
         IF (ICURFSS.EQ. 4) DELS1 = 'KK'
         IF (ICURFSS.EQ. 5) DELS1 = 'DK'
         IF (ICURFSS.EQ. 6) DELS1 = 'KD'
         DELS2 = 'DDDD'
         IF (ICURFSS.EQ. 2) DELS2 = 'DDDK'
         IF (ICURFSS.EQ. 3) DELS2 = 'KDDD'
         IF (ICURFSS.EQ. 4) DELS2 = 'KDDK'
         IF (ICURFSS.EQ. 5) DELS2 = 'DDKK'
         IF (ICURFSS.EQ. 6) DELS2 = 'KKDD'

C        Initialize S1 and S2 in case we want to do zero
C        iterations (tests). In normal runs S1 and S2 will
C        be initialized in T1EQNS and T2EQNS.

         CALL XCOPY(NDIMT1,A0,0,S1,1)
         CALL XCOPY(NDIMT2,A0,0,S2,1)
         CALL S1S2_ADDFO (FVO,BUF1,BUF2,S1,S2)

      ELSE
         CALL QUIT('Wrong keyword in call to CCEQNS')
      ENDIF

C     Zero amplitudes if the NOCCX options are specified

      IF (TEQ.AND.NOCCS) CALL XCOPY (NDIMT1,A0,0,T1,1)
      IF (TEQ.AND.NOCCD) CALL XCOPY (NDIMT2,A0,0,T2,1)

C     Initialize the work arrays to zero.

      CALL XCOPY (NDIMT2,A0,0,BUF1,1)
      CALL XCOPY (NBUF2,A0,0,BUF2,1)
      CALL XCOPY (NBUF3,A0,0,BUF3,1)

C     Initialize DIIS for this sector, first zeroing amplitudes from higher sectors.

      IT = -1
      CALL DELFCK ('VO',DELS1,T1)
      CALL DELINT ('VVOO',DELS2,T2,0,0,0)
      CALL DIISX(NDIMT1,NDIMT2,T1,T2,
     &           IT,CC1,BB1,ITC,BB2,RMS)

      ECCSD = 0.D0
      NIT = 0

      CALL T1DIAG(DIAG,T1,BUF1)
      CALL CCENG (T1,T2,FVO,BUF1,BUF2,ECCSD,ECCSDIM)
      IF (IPRNT.GE.3) CALL PAIREN (T1,T2,FVO,BUF1,BUF2)

      IF (IPRNT.GE.1) THEN
         WRITE(IW,6003)
         WRITE(IW,6001) NIT,ECCSD,RMS,DIAG
      ENDIF
                     
 6001 FORMAT(2X,I3,5X,2F20.15,F10.5)
 6003 FORMAT(//3X,'NIT',10X,'ENERGY',20X,'RMS',7X,'T1-DIAGN')

      END
