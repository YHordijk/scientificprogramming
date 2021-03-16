!      Copyright (c) 2019 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://www.diracprogram.org

C
C
#ifdef OLD_LOG
/* Comdeck geoana_log */
950505-hjaaj:
 RADIUS: bug fix - made function double precision instead of real
         (real values would not always be read correctly by geoan1)
942508-kr: As P.Taylor pointed out, IPRINT never used, and thus removed from
CBIANA.
900207-hjaaj: check for non-existent atoms for .ANGLES and .DIHEDR in GEOANA.
Idea: make angle output as in zmat???
#endif
C  /* Deck anainp */
      SUBROUTINE ANAINP(WORD)
C
C  5-Jul-1985 Hans Jorgen Aa. Jensen
C
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
      PARAMETER (NTABLE = 4)
      PARAMETER (MAXANG = 20)
      LOGICAL NEWDEF
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
#include "abainf.h"
      LOGICAL SKIP
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
      DATA TABLE /'.SKIP  ', '.XXXXXX', '.ANGLES', '.DIHEDR'/
      DATA MANG/0/, MDIHED/0/
C
      CALL ANAINI
C
      NEWDEF = (WORD .EQ. '*GEOANA')
      ICHANG = 0
      IF (NEWDEF) THEN
         WORD1 = WORD
  100    CONTINUE
            READ (LUCMD, '(A7)') WORD
            CALL UPCASE(WORD)
            PROMPT = WORD(1:1)
            IF (PROMPT .EQ. '!' .OR. PROMPT .EQ. '#') THEN
               GO TO 100
            ELSE IF (PROMPT .EQ. '.') THEN
               ICHANG = ICHANG + 1
               DO 200 I = 1, NTABLE
                  IF (TABLE(I) .EQ. WORD) THEN
                     GO TO (1,2,3,4), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in ANAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in ANAINP.')
    1          CONTINUE
                  SKIP = .TRUE.
               GO TO 100
    2          CONTINUE
               GO TO 100
    3          CONTINUE
                  READ (LUCMD,*) NANG
                  MANG = MIN(MAXANG,NANG)
                  DO 310 I = 1,MANG
                     READ(LUCMD,*) (IANG(J,I),J=1,3)
  310             CONTINUE
                  MANG = NANG - MANG
                  DO 320 I = 1,MANG
                     READ(LUCMD,'()')
  320             CONTINUE
               GO TO 100
    4          CONTINUE
                  READ (LUCMD,*) NDIHED
                  MDIHED = MIN(MAXANG,NDIHED)
                  DO 410 I = 1,MDIHED
                     READ(LUCMD,*) (IDIHED(J,I),J=1,4)
  410             CONTINUE
                  MDIHED = NDIHED - MDIHED
                  DO 420 I = 1,MDIHED
                     READ(LUCMD,'()')
  420             CONTINUE
               GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     *            '" not recognized in ANAINP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ANAINP.')
            END IF
      END IF
  300 CONTINUE
      IF (ICHANG .GT. 0) THEN
         CALL HEADER('Changes of defaults for GEOANA:',0)
         IF (SKIP) THEN
            WRITE (LUPRI,'(A)') ' GEOANA skipped in this run.'
         ELSE
            IF (NANG .GT. 0) THEN
               WRITE (LUPRI,'(/A/)')
     *            ' Following angles will be calculated:'
               DO 1310 I = 1,NANG
                  WRITE (LUPRI,'(I10,A,4I5)') I,' : ',(IANG(J,I),J=1,3)
 1310          CONTINUE
               IF (MANG .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3,A)') ' The last',MANG,
     *               ' angles specified go beyond current maximum',
     *               ' and will not be printed.'
               END IF
            END IF
            IF (NDIHED .GT. 0) THEN
               WRITE (LUPRI,'(/A/)')
     *            ' Following dihedral angles will be calculated:'
               DO 1410 I = 1,NDIHED
                  WRITE (LUPRI,'(I10,A,4I5)')I,' : ',(IDIHED(J,I),J=1,4)
 1410          CONTINUE
               IF (MDIHED .GT. 0) THEN
                  WRITE (LUPRI,'(/A,I3,A)') ' The last',MDIHED,
     *               ' dihedral angles specified go beyond current',
     *               ' maximum and will not be printed.'
               END IF
            END IF
         END IF
         WRITE (LUPRI,'(/)')
      END IF
      RETURN
      END
C  /* Deck anaini */
      SUBROUTINE ANAINI
C
C     Initialize /CBIANA/
C
#include "implicit.h"
      PARAMETER (MAXANG = 20)
      LOGICAL SKIP
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
      NANG   = 0
      NDIHED = 0
      SKIP   = .FALSE.
      RETURN
      END
C  /* Deck geoana */
      SUBROUTINE GEOANA(COORD,PRINT,DIF,NBONDS,PUNCH,WORK,LWORK)
C
C 30-Jun-1985 Hans Jorgen Aa. Jensen
C Modified for symmetry 25-Sep-1989 tuh
C Modified for differential geomtries 18-Oct-1989 tuh
C
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      LOGICAL PRINT, DIF, PUNCH
      DIMENSION WORK(LWORK)
#include "nuclei.h"
C
      CALL QENTER('GEOANA')
      KFREE = 1
      LFREE = LWORK
      CALL MEMGET('REAL',KVEC,3*NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KDIST2,NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('LOGI',KBOND ,NUCDEP*NUCDEP,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KCHRG ,NUCDEP,       WORK,KFREE,LFREE)
      CALL MEMGET('INTE',IPAIR,2*NUCDEP*NUCDEP,WORK,KFREE,LFREE)
C
      CALL GEOAN1(COORD,PRINT,DIF,NBONDS,PUNCH,WORK(KVEC),
     &            WORK(KDIST2),WORK(KBOND),WORK(KCHRG),WORK(IPAIR))
C
      CALL MEMREL('GEOANA',WORK,1,1,KFREE,LFREE)
      CALL QEXIT('GEOANA')
      RETURN
      END
C  /* Deck geoan1 */
      SUBROUTINE GEOAN1(COORD,PRINT,DIF,NBONDS,PUNCH,VEC,DIST2,
     &                  BONDED,ICHARG,IPAIRS)
C
C Modified for more selective printing of bonded atoms,
C     Jan-1995 Hanne Heiberg
C     Added punching atom bonds for Gamess graphic output, K.Ruud-95
C
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
      PARAMETER (MAXANG = 20)
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0)
      LOGICAL SKIP, DIF, PRINT, PUNCH
      COMMON /CBIANA/ IANG(3,MAXANG),IDIHED(4,MAXANG),NANG,NDIHED,
     *                SKIP
C
#include "nuclei.h"
#include "symmet.h"
#include "dgroup.h"
C
      PARAMETER (LUPUNCH = 4)
      DIMENSION COORD(3,*),
     *          DIST(MXCENT*(MXCENT+1)/2),
     *          ANGLE(MAXANG), DIHED(MAXANG), ICHARG(NUCDEP)
      DIMENSION VEC(3,NUCDEP,NUCDEP), DIST2(NUCDEP,NUCDEP), 
     &          IPAIRS(2,NUCDEP*NUCDEP)
      LOGICAL   BONDED(NUCDEP,NUCDEP)
      CHARACTER*4 NUCNAM(4)
      SAVE DIST, ANGLE, DIHED
#include "facang.h"
      ARCCOS(ARG) = FACANG*ACOS(ARG)
#include "ibtfun.h"
      IPAK(I,J) = MAX(I,J)*(MAX(I,J)-1)/2 + MIN(I,J)
C
      NBONDS = 0
      IF (NUCDEP .EQ. 1) RETURN
C
C     set up bond vectors
C
      IATOMA = 0
      DO 100 ICENTA = 1, NUCIND
         DO 110 IA = 0, MAXOPR
            IF (IBTAND(IA,ISTBNU(ICENTA)) .EQ. 0) THEN
               IATOMA = IATOMA + 1
               CXA = PT(IBTAND(ISYMAX(1,1),IA))*COORD(1,ICENTA)
               CYA = PT(IBTAND(ISYMAX(2,1),IA))*COORD(2,ICENTA)
               CZA = PT(IBTAND(ISYMAX(3,1),IA))*COORD(3,ICENTA)
C
               ICHARG(IATOMA) = IZATOM(ICENTA)
C
               IATOMB = 0
               DO 200 ICENTB = 1, NUCIND
                  DO 210 IB = 0, MAXOPR
                     IF (IBTAND(IB,ISTBNU(ICENTB)) .EQ. 0) THEN
                        IATOMB = IATOMB + 1
                        IF (IATOMB .GT. IATOMA) GO TO 110
C                       ... next IATOMA, only IATOMB .le. IATOMA needed
                          CXB=PT(IBTAND(ISYMAX(1,1),IB))*COORD(1,ICENTB)
                          CYB=PT(IBTAND(ISYMAX(2,1),IB))*COORD(2,ICENTB)
                          CZB=PT(IBTAND(ISYMAX(3,1),IB))*COORD(3,ICENTB)
                          VEC(1,IATOMB,IATOMA) = XTANG*(CXA - CXB)
                          VEC(2,IATOMB,IATOMA) = XTANG*(CYA - CYB)
                          VEC(3,IATOMB,IATOMA) = XTANG*(CZA - CZB)
                          VEC(1,IATOMA,IATOMB) = -VEC(1,IATOMB,IATOMA)
                          VEC(2,IATOMA,IATOMB) = -VEC(2,IATOMB,IATOMA)
                          VEC(3,IATOMA,IATOMB) = -VEC(3,IATOMB,IATOMA)
                     END IF
  210             CONTINUE
  200          CONTINUE
            END IF
  110    CONTINUE
  100 CONTINUE
C
C     Set up distance matrix in Angstrom
C
      IJ = 0
      DO 400 I = 1,NUCDEP
         DO 300 J = 1,I
            IJ = IJ + 1
            DISTAN = VEC(1,J,I)*VEC(1,J,I) + VEC(2,J,I)*VEC(2,J,I)
     *             + VEC(3,J,I)*VEC(3,J,I)
            DISTAN = SQRT(DISTAN)
            DIST2(I,J) = DISTAN
            DIST2(J,I) = DISTAN
            IF (.NOT.DIF) THEN
               DIST(IJ) = DISTAN
            ELSE
               DIST(IJ) = DISTAN - DIST(IJ)
            END IF
  300    CONTINUE
  400 CONTINUE
C
      IF (PRINT) THEN
         IF (DIF) THEN
            CALL HEADER
     *        ('Differential interatomic separations (in Angstroms):',2)
         ELSE
            CALL HEADER('Interatomic separations (in Angstroms):',2)
         END IF
         CALL PRIDIS(NAMDPX,DIST,NUCDEP)
      END IF
C
      DO 10, J= 1,NUCDEP
         RADJ = RADIUS(ICHARG(J))
         DO 20, I= 1, J-1
            RADI = RADIUS(ICHARG(I))
            IF (RADI .LT. D0 .OR. RADJ .LT. D0) THEN
               WRITE(LUPRI,*) 'RADII FOR ATOMS WITH EITHER/OR ATOMIC'
     &              //' NUMBER ', ICHARG(I), ' OR ', ICHARG(J),
     &              ' ARE UNAVAILABLE '
               BONDED(I,J) = .FALSE.
               BONDED(J,I) = .FALSE.
            ELSE IF (DIST2(I,J) .LT. (1.2D0 * (RADI + RADJ))) THEN
               NBONDS = NBONDS + 1
               IPAIRS(1,NBONDS) = I
               IPAIRS(2,NBONDS) = J
               BONDED(I,J) = .TRUE.
               BONDED(J,I) = .TRUE.
            ELSE
               BONDED(I,J) = .FALSE.
               BONDED(J,I) = .FALSE.
            END IF
 20      CONTINUE
         BONDED(J,J) = .FALSE.
 10   CONTINUE
C
      IF (PRINT .AND. .NOT. DIF) THEN
         CALL HEADER('Bond distances (angstroms):',1)
         WRITE (LUPRI,'(14X,A/14X,A)')
     $     '    atom 1     atom 2                           distance',
     $     '    ------     ------                           --------'
         DO 500, I = 1, NUCDEP
            DO 510, J = 1, I-1
               IF (BONDED(I,J)) THEN
                  NUCNAM(1) = NAMDPX(3*I)(1:4)
                  NUCNAM(2) = NAMDPX(3*J)(1:4)
                  WRITE(LUPRI,'(2X,A,4X,A4,7X,A4,22X,F13.6)')
     $                 'bond distance:',
     &                  NUCNAM(1), NUCNAM(2), DIST2(I,J)
               END IF
 510        CONTINUE
 500     CONTINUE
C
      IF (NUCDEP .GT. 2 .AND. NANG .LE. 0) THEN
         CALL HEADER('Bond angles (degrees):',1)
         WRITE (LUPRI,'(14X,A/14X,A)')
     $   '    atom 1     atom 2     atom 3                   angle',
     $   '    ------     ------     ------                   -----'
C
         IJK = 0
         DO 40, I= 1,NUCDEP
            DO 50, J= 1,NUCDEP
               DO 60, K= 1, J - 1
                  IF (BONDED(I,J) .AND. BONDED(I,K)) THEN
                     IJK = IJK + 1
                     NUCNAM(1) = NAMDPX(3*J)(1:4)
                     NUCNAM(2) = NAMDPX(3*I)(1:4)
                     NUCNAM(3) = NAMDPX(3*K)(1:4)
                     ANG = VEC(1,I,J) * VEC(1,I,K)
     *                   + VEC(2,I,J) * VEC(2,I,K)
     *                   + VEC(3,I,J) * VEC(3,I,K)
                     ANG = (ANG / DIST2(I,J)) / DIST2(I,K)
                     IF (ABS(ANG) .GT. D1) ANG = SIGN(D1,ANG)
                     ANG = ARCCOS(ANG)
                     WRITE(LUPRI,'(2X,A,7X, A4,7X, A4,7X,A4,11X,F13.3)')
     *                  'bond angle:',NUCNAM(1),NUCNAM(2),NUCNAM(3),ANG
                  END IF
 60            CONTINUE
 50         CONTINUE
 40      CONTINUE
         IF (IJK .EQ. 0) WRITE(LUPRI,'(5X,A)') 'No angles found'
      END IF
      END IF
C
C     Punch bonding information in Gamess output format on unit LUPUNCH
C
      IF (PUNCH) THEN
         IF(NBONDS.LE.6) THEN
            WRITE(LUPUNCH,8010) (IPAIRS(1,I),IPAIRS(2,I),I=1,NBONDS)
         ELSE
            WRITE(LUPUNCH,8020) (IPAIRS(1,I),IPAIRS(2,I),I=1,6)
            WRITE(LUPUNCH,8030) (IPAIRS(1,I),IPAIRS(2,I),I=7,NBONDS)
         END IF
      END IF
C         
      IF (NANG .GT. 0) THEN
         IF (PRINT) THEN
            CALL HEADER('Angles according to input list:',2)
            WRITE (LUPRI,'(A/A)')
     *       '    atom 1     atom 2     atom 3         angle (degrees)',
     *       '    ------     ------     ------         ---------------'
         END IF
         DO 1000 I = 1,NANG
            I1 = IANG(1,I)
            I2 = IANG(2,I)
            I3 = IANG(3,I)
            IMX = MAX(I1,I2,I3)
            IF (IMX .GT. NUCDEP) THEN
               IF (PRINT) WRITE (LUPRI,'(/A/)')
     &            ' *GEOANA input error for .ANGLES: non-existent atom'
               GO TO 1000
            END IF
            NUCNAM(1) = NAMDPX(3*I1)(1:4)
            NUCNAM(2) = NAMDPX(3*I2)(1:4)
            NUCNAM(3) = NAMDPX(3*I3)(1:4)
            IF (I1 .NE. I2 .AND. I2 .NE. I3) THEN
               ANG = VEC(1,I2,I1)*VEC(1,I2,I3)
     *             + VEC(2,I2,I1)*VEC(2,I2,I3)
     *             + VEC(3,I2,I1)*VEC(3,I2,I3)
               ANG = ANG / DIST2(I2,I1) / DIST2(I2,I3)
               IF (ABS(ANG) .GT. D1) ANG = SIGN(D1,ANG)
               ANG = ARCCOS(ANG)
               IF (.NOT.DIF) THEN
                  ANGLE(I) = ANG
               ELSE
                  ANGLE(I) = ANG - ANGLE(I)
               END IF
               IF (PRINT) WRITE (LUPRI,'(5X,A4,7X,A4,7X,A4,F20.3)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),ANGLE(I)
            ELSE
               IF (PRINT) WRITE (LUPRI,'(5X,A4,7X,A4,7X,A4,10X,A)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),'undefined'
            END IF
 1000    CONTINUE
      END IF
C
      IF (NDIHED .GT. 0) THEN
          IF (PRINT) WRITE (LUPRI,'(//A,A/A,A)')
     *       '    atom 1     atom 2     atom 3     atom 4',
     *       '    dihedral angle (degrees)',
     *       '    ------     ------     ------     ------',
     *       '    ------------------------'
         DO 2000 I = 1,NDIHED
            I1 = IDIHED(1,I)
            I2 = IDIHED(2,I)
            I3 = IDIHED(3,I)
            I4 = IDIHED(4,I)
            IMX = MAX(I1,I2,I3,I4)
            IF (IMX .GT. NUCDEP) THEN
               IF (PRINT) WRITE (LUPRI,'(/A/)')
     &            ' *GEOANA input error for .DIHEDR: non-existent atom'
               GO TO 2000
            END IF
            NUCNAM(1) = NAMDPX(3*I1)(1:4)
            NUCNAM(2) = NAMDPX(3*I2)(1:4)
            NUCNAM(3) = NAMDPX(3*I3)(1:4)
            NUCNAM(4) = NAMDPX(3*I4)(1:4)
            X1 = VEC(2,I2,I1)*VEC(3,I2,I3) - VEC(2,I2,I3)*VEC(3,I2,I1)
            X2 = VEC(3,I2,I1)*VEC(1,I2,I3) - VEC(3,I2,I3)*VEC(1,I2,I1)
            X3 = VEC(1,I2,I1)*VEC(2,I2,I3) - VEC(1,I2,I3)*VEC(2,I2,I1)
            Y1 = VEC(2,I3,I2)*VEC(3,I3,I4) - VEC(2,I3,I4)*VEC(3,I3,I2)
            Y2 = VEC(3,I3,I2)*VEC(1,I3,I4) - VEC(3,I3,I4)*VEC(1,I3,I2)
            Y3 = VEC(1,I3,I2)*VEC(2,I3,I4) - VEC(1,I3,I4)*VEC(2,I3,I2)
            Z1 = X2*Y3 - X3*Y2
            Z2 = X3*Y1 - X1*Y3
            Z3 = X1*Y2 - X2*Y1
            SENSE = Z1*VEC(1,I2,I3) + Z2*VEC(2,I2,I3) + Z3*VEC(3,I2,I3)
            SENSE = SIGN(D1,SENSE)
            ANG = X1*Y1 + X2*Y2 + X3*Y3
            DDD = (X1*X1 + X2*X2 + X3*X3) * (Y1*Y1 + Y2*Y2 + Y3*Y3)
            IF (DDD .GT. 1.D-10) THEN
               ANG = ANG / SQRT(DDD)
               IF (ABS(ANG) .GT. D1) ANG = SIGN(D1,ANG)
               ANG = SENSE*ARCCOS(ANG)
               IF (.NOT.DIF) THEN
                  DIHED(I) = ANG
               ELSE
                  DIHED(I) = ANG - DIHED(I)
               END IF
               IF (PRINT) WRITE(LUPRI,'(5X,A4,7X,A4,7X,A4,7X,A4,F20.3)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),NUCNAM(4),DIHED(I)
            ELSE
               IF (PRINT) WRITE(LUPRI,'(5X,A4,7X,A4,7X,A4,7X,A4,10X,A)')
     *            NUCNAM(1),NUCNAM(2),NUCNAM(3),NUCNAM(4),'undefined'
            END IF
 2000    CONTINUE
      END IF
C
      IF (PRINT) WRITE (LUPRI,'(/)')
      RETURN
 8010 FORMAT('BONDATOMS ',6(I4,I4,2X))
 8020 FORMAT('BONDATOMS ',6(I4,I4,2X),' >')
 8030 FORMAT(7(I4,I4,2X),:,' >')
      END
C  /* Deck pridis */
      SUBROUTINE PRIDIS (NAMEX,DISMAT,NROW)
C
C 30-Jun-1985 Hans Jorgen Aa. Jensen
C (based on OUTPAK by Nelson H.F. Beebe)
C
C Print bond distance matrix (or other matrix over atoms)
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO=0.D00, KCOL=6)
      INTEGER BEGIN
      CHARACTER*6 NAMEX(*)
      DIMENSION DISMAT(*)
C
      LAST = MIN(NROW,KCOL)
      BEGIN = 1
 1050 NCOL = 1
      WRITE (LUPRI,1000) (NAMEX(3*I)(1:4),I = BEGIN,LAST)
      WRITE (LUPRI,'()')
      DO 40 K = BEGIN,NROW
         KTOTAL = (K*(K-1))/2 + BEGIN - 1
         WRITE (LUPRI,2000) ' ',NAMEX(3*K)(1:4),
     *      (DISMAT(KTOTAL+J),J = 1,NCOL)
         IF (K .LT. (BEGIN+KCOL-1)) NCOL = NCOL + 1
   40 CONTINUE
      WRITE (LUPRI,'()')
      LAST = MIN(LAST+KCOL,NROW)
      BEGIN = BEGIN+NCOL
      IF (BEGIN.LE.NROW) GO TO 1050
      WRITE (LUPRI,'()')
      RETURN
 1000 FORMAT (8X,6(4X,A4,4X),(4X,A4,4X))
 2000 FORMAT (A1,2X,A4,6F12.6)
      END
C  /* Deck radius */
      FUNCTION RADIUS(NCHARGE)
#include "implicit.h"
C     Returns atomic radius in Angstrom;
C     based on covalent radii and metallic radii. 
C     Returns radius of MAXCHR when data not available.
C     Revision March 17 1999 - Trond Saue
      PARAMETER ( DPH = 1.0D-02,MAXCHR=83)
      DIMENSION RAD(MAXCHR)
      DATA (RAD(I), I = 1, MAXCHR)/
     & 30.D0,  155.D0,  160.D0,  110.D0,   90.D0,
     & 80.D0,   70.D0,   68.D0,   65.D0,  154.D0,
     &190.D0,  160.D0,  140.D0,  110.D0,  110.D0,
     &105.D0,  105.D0,  190.D0,  238.D0,  200.D0,
     &165.D0,  145.D0,  135.D0,  130.D0,  125.D0,
     &125.D0,  125.D0,  125.D0,  125.D0,  140.D0,
     &140.D0,  130.D0,  120.D0,  120.D0,  120.D0,
     &200.D0,  255.D0,  215.D0,  180.D0,  160.D0,
     &145.D0,  140.D0,  135.D0,  130.D0,  130.D0,
     &135.D0,  140.D0,  155.D0,  160.D0,  160.D0,
     &140.D0,  140.D0,  140.D0,  220.D0,  270.D0,
     &220.D0,  185.D0,  180.D0,  180.D0,  180.D0,
     &180.D0,  180.D0,  200.D0,  180.D0,  175.D0,
     &175.D0,  175.D0,  175.D0,  170.D0,  170.D0,
     &170.D0,  155.D0,  145.D0,  140.D0,  135.D0,
     &135.D0,  135.D0,  135.D0,  145.D0,  155.D0,
     &170.D0,  175.D0,  170.D0/

C
      IF(NCHARGE.GT.MAXCHR) THEN
        RADIUS = DPH*RAD(MAXCHR)
      ELSEIF(NCHARGE.LE.0) THEN
        RADIUS = -1.0D0
      ELSE
        RADIUS = DPH*RAD(NCHARGE)
      ENDIF
      RETURN
      END
C  /* Deck vdwrad */
      FUNCTION VDWRAD(NCHARGE)
#include "implicit.h"
C     based on van der Waals radii. 
C     returns radius of MAXCHR when data not available
C     Revised by T.Saue March 17 1999
      PARAMETER ( DPH = 1.0D-02,MAXCHR=24)
      DIMENSION RAD(MAXCHR)
      DATA (RAD(I), I = 1, MAXCHR)/
     &110.,  220.,  122.,   63.,  155.,
     &155.,  140.,  135.,  130.,  154.,
     &190.,  160.,  140.,  110.,  202.,
     &220.,  150.,  150.,  220.,  188.,
     &181.,  175.,  277.,  239./
C
      IF(NCHARGE.GT.MAXCHR) THEN
        VDWRAD = DPH*RAD(MAXCHR)
      ELSE
        VDWRAD = DPH*RAD(NCHARGE)
      ENDIF
      RETURN
      END
