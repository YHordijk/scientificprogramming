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

! FILE    : hersym.F

! ::: HERSYM_DEBUG - Debug print if > 0 + IPRSYM set to HERSYM_DEBUG
#define HERSYM_DEBUG -1
! ::: IPRTHR - Print level threshold for verbose output
!#define IPRTHR 3 ! Miro: I had to deactivate it due to runtime checks !
!
!========================================================================
!950208-Vebjoern Bakken
!SYMADD: Completely new module to determine molecular symmetry, find the
!        symmetry elements, remove symmetry-dependent atoms and finally
!        return a string containing the full point group.
!========================================================================
C  /* Deck symadd */
      SUBROUTINE SYMADD(WORK,LWORK,NSYMOP,KATOM,KASYM,
     &                  POINT_GROUP,TOLLRN_input,IPRSYM_input,
     &                  internal_file)
C*****************************************************************************
C
C     This module determines molecular symmetry automatically
C
C     vebjorn - 950204 - completely new module
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "mxsymm.h"

! infpar.h : MYTID, MPARID
! optinf.h : VRML_SYM
#include "infpar.h"
#include "optinf.h"

      DIMENSION WORK(LWORK)
      CHARACTER KASYM(3,3)*1, POINT_GROUP*15
      LOGICAL :: internal_file

      CALL QENTER('SYMADD')

      IPRTHR = 3 ! Miro: for runtime check it must be defined as variable
      IPRSYM = MAX(HERSYM_DEBUG, IPRSYM_input)

C     Set control variables in "mxsymm.h"
      IF (TOLLRN_INPUT .LE. 0.0D0) THEN
         TOLLRN = 5.0D-6
      ELSE
         TOLLRN = TOLLRN_INPUT
      END IF
      TOLLRN_2 = TOLLRN**2
      ZERTOL = 1.0D-12
      MAXAXS = 100
      MAXMIR = 25
CSK   only MPARID needs to print the information
      IF( MYTID .eq. MPARID .OR. IPRSYM .GT. 0)THEN
        CALL HEADER('SYMADD: Detection of molecular symmetry',-1)
        WRITE(LUPRI,'(A,1P,G10.2)') ' Symmetry test threshold:', TOLLRN
      END IF

      KFRSAV = 1
      KFREE  = KFRSAV
      LFREE  = LWORK
      CALL MEMGET2('REAL','ATMARR',KCOOR ,6*(MXCENT+1),WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DRTAXS',KROT  ,6*(MAXAXS+1),WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DIRAXS',KIROT ,6*(MAXAXS+1),WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DMRPLN',KMIR  ,6*(MAXMIR+1),WORK,KFREE,LFREE)
      LENTMP = MAX(2*MXCENT,MXCENT+MAXAXS+MAXMIR)
      CALL MEMGET2('REAL','TEMP ',KTEMP ,7*(LENTMP+1),WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','TEMP2',KTEMP2,6*(MXCENT+1),WORK,KFREE,LFREE)
      CALL FIND_SYMMETRY(WORK(KCOOR),WORK(KROT),WORK(KIROT),WORK(KMIR),
     &      WORK(KTEMP),WORK(KTEMP2),NSYMOP,
     &      KATOM,KASYM,LENTMP,POINT_GROUP,VRML_SYM,
     &      internal_file,WORK,KFREE,LFREE)
      CALL MEMREL('SYMADD',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('SYMADD')
      RETURN
      END
C  /* Deck FIND_SYMMETRY */
      SUBROUTINE FIND_SYMMETRY(ATMARR,DRTAXS,DIRAXS,DMRPLN,
     &            TMPARR,TMP2AR,NSYMOP,
     &            KATOM,KASYM,LENTMP,POINT_GROUP,VRML_SYM,
     &            internal_file,WORK,KFREE,LFREE)
C
#include "implicit.h"
#include "mxsymm.h"
#include "mxcent.h"
#include "nuclei.h"
#include "priunit.h"
#include "cbisol.h"
#include "cbihr1.h"
#include "molinp.h"
#include "numder.h"
! infpar.h: MPARID, MYTID
#include "infpar.h"
C ***********************************************************
C *   ATMARR - Atomic coordinates and more
C *   POINT_GROUP  - Full point group
C *   DRTAXS - Rotational axes
C *   DIRAXS - Improper rotational axes
C *   DMRPLN - Mirror planes
C *   TMPARR + TMP2AR - Temporary arrays
C ***********************************************************
      DIMENSION ATMARR(6,0:MXCENT), DRTAXS(6,0:MAXAXS)
      DIMENSION DIRAXS(6,0:MAXAXS), DMRPLN(6,0:MAXMIR)
      DIMENSION TMPARR(7,0:LENTMP), WORK(*)
      DIMENSION TMP2AR(6,0:MXCENT)
      LOGICAL   VRML_SYM, LINMOL
      LOGICAL :: internal_file
      INTEGER, ALLOCATABLE :: IBUF(:)
      PARAMETER (D0 = 0.0D0)
      CHARACTER*1 KASYM(3,3)
      CHARACTER SYMM*9, TMPLN*80
      CHARACTER POINT_GROUP*15, TEXT*43

      IPRTHR = 3 ! Miro: for runtime check it must be defined as variable
C
C     Just to make sure, we zero all temporary arrays
C
      CALL DZERO(ATMARR,6*(MXCENT + 1))
      CALL DZERO(DRTAXS,6*(MAXAXS + 1))
      CALL DZERO(DIRAXS,6*(MAXAXS + 1))
      CALL DZERO(DMRPLN,6*(MAXMIR + 1))
      CALL DZERO(TMPARR,7*(LENTMP + 1))
      CALL DZERO(TMP2AR,6*(MXCENT + 1))
C     IF (SOLVNT .OR. DORLM) THEN
C       N_ATOMS = NUCIND - 1
C     ELSE
        N_ATOMS = NUCIND
C     END IF
      DO I = 1, N_ATOMS
         ATMARR(1,I) = CORD(1,I)
         ATMARR(2,I) = CORD(2,I)
         ATMARR(3,I) = CORD(3,I)
Chjaaj   ATMARR(4,I) = CHARGE(I)
         ATMARR(4,I) = IZATOM(I)
C        ... use true nuclear charge, not effective charge in CHARGE(I)
C            (As different atoms can have same effective charge with
C             ECP, using CHARGE for sorting is not reliable for ECP)
C            /hjaaj Mar 2004
         ATMARR(5,I) = I
         ATMARR(6,I) = ISOTOP(I)
      END DO
      ATMARR(1,0) = N_ATOMS
      IF ( MYTID .eq. MPARID .and. IPRSYM .gt. IPRTHR) THEN
         TEXT = 'Original Coordinates'
         CALL PRNARR(TEXT, MXCENT, ATMARR)
      END IF
         IF (.NOT. NOMOVE) THEN
C           Place molecular center of mass at the origin and
C           with principal axes of inertia along coordinate axes
            CALL CENTER_MOLECULE(ATMARR)
         END IF
C     Sort atoms after atom number and distance from origo
CRF      We have to sort the atoms AFTER we center the molecule
C        Else errors occur when non-standard isotopes are used
      CALL SRTATM(MXCENT, ATMARR, .TRUE.)
      IF ( MYTID .eq. MPARID .and. IPRSYM .gt. IPRTHR) THEN
         TEXT = 'Coordinates after sort'
         CALL PRNARR(TEXT, MXCENT, ATMARR)
      END IF

C     Determine molecular point group
      CALL FIND_PGROUP(MXCENT, ATMARR, DRTAXS, DIRAXS, DMRPLN,
     &                 TMPARR, POINT_GROUP)

      IF (POINT_GROUP(1:4).EQ.'C(1)') THEN
         !  if moving allowed, rotate molecule without symmetry so
         !  principle axes of inertia follow the coordinate axes
         IF (.NOT. NOMOVE) CALL PRINCIPLE_AXIS_MOLECULE(ATMARR)

         IF (VRML_SYM) WRITE(LUPRI,'(/A)')
     &   ' INFO: VRML file with symmetry info not created'//
     &   ' because molecule has C(1) symmetry.'

      ELSE IF (VRML_SYM) THEN
C        ::: If requested, a VRML representation of :::
C        ::: the symmetry elements is constructed.  :::
         CALL MKVRSY(ATMARR,DRTAXS,MAXAXS,
     &               DMRPLN,MAXMIR,WORK(KFREE),LFREE)
      END IF
C
      IF ( MYTID .eq. MPARID  .OR. IPRSYM .GT. 0) WRITE(LUPRI,'(/2A)')
     &      ' Symmetry point group found: ', POINT_GROUP
      LENTMP = NINT(ATMARR(1,0) + DRTAXS(1,0) + DMRPLN(1,0))
      IF (.NOT. NOMOVE) THEN
         CALL TURN_MOLECULE(MXCENT,ATMARR,DRTAXS,DIRAXS,DMRPLN,
     &                      TMPARR,LENTMP)
         IF (MYTID .eq. MPARID .and. iprsym .gt. iprthr) then
            TEXT = 'Centered and Rotated in original atom order'
            KWORK = KFREE
            allocate(IBUF(N_ATOMS))
            CALL DCOPY(N_ATOMS,ATMARR(5,1),MXCENT+1,TMPARR,1) ! atmarr(5,1) has original atom order
            CALL INDEXX(N_ATOMS,TMPARR,IBUF)
            CALL PRNARR2(text, mxcent, atmarr, ibuf) ! prints atoms in original order
            deallocate(IBUF)
         END IF
      ELSE IF ( MYTID .EQ. MPARID .OR. IPRSYM .GT. 0) THEN
         WRITE (LUPRI,'(/A/)') ' Molecule was not centered nor'//
     &      ' rotated as requested in input.'
      END IF
      IF (IPRSYM .GE. MIN(IPRTHR,1)) THEN
      IF (.NOT. LINMOL(MXCENT, ATMARR)) THEN
         IF (NINT(DRTAXS(1,0)) .EQ. 0) THEN
            WRITE(LUPRI,'(/A)') ' No rotational axes were found.'
         ELSE
           TEXT = 'Rotational Axes'
           CALL PRNARR(TEXT, MAXAXS, DRTAXS)
         END IF
         IF (NINT(DIRAXS(1,0)) .LT. 0) THEN
            WRITE(LUPRI,'(/A)') ' Due to number of improper '
     &        // 'rotational axes, they were not determined.'
         ELSE IF (NINT(DIRAXS(1,0)) .EQ. 0) THEN
            WRITE(LUPRI,'(/A)')
     &        ' No unique improper rotational axes were found.'
         ELSE
            TEXT = 'Improper Rotational Axes'
            CALL PRNARR(TEXT, MAXAXS, DIRAXS)
         END IF
         IF (NINT(DMRPLN(1,0)) .EQ. 0) THEN
            WRITE(LUPRI,'(/A)') ' No mirror planes were found.'
         ELSE
            CALL PRNMIR(MAXMIR, DMRPLN)
         END IF
      END IF
      END IF

C ::: Find symmetry operations which can be used by Dalton/Dirac
C ::: and remove symmetry-dependent atoms
      CALL SYM_REDUCE(MXCENT, ATMARR,TMPARR, TMP2AR, SYMM,
     &      internal_file)
      IF (IPRSYM .GE. MIN(3,IPRTHR)) THEN
         TEXT = 'Symmetry Independent Centres'
         CALL PRNARR(TEXT, MXCENT, ATMARR)
      END IF

C ::: The new symmetry is added. :::
      NSYMOP = 0
      DO 200 I = 1, 3
         DO 210 J = 1, 3
            KASYM(J,I) = SYMM((I-1)*3+J:(I-1)*3+J)
 210     CONTINUE
         IF (SYMM((I-1)*3+1:(I-1)*3+3) .NE. '   ') NSYMOP=NSYMOP+1
 200  CONTINUE
C ::: The input file is modified to include symmetry. :::
C ::: If there is any Aangstroem mark, it is removed. :::
      if (internal_file) then
        TMPLN = MLINE(NMLAU)
        WRITE(TMPLN(10:10), '(I1)') NSYMOP
        WRITE(TMPLN(11:20), '(A9,A1)') SYMM, ' '
      end if
      IF (MYTID .eq. MPARID .OR. IPRSYM .GE. 1) THEN
         IF (NSYMOP .EQ. 0) THEN
            WRITE(LUPRI, '(/A)') ' No symmetry elements were found.'
         ELSE IF (NSYMOP .EQ. 1) THEN
            WRITE(LUPRI, '(/2A)')
     &         ' The following symmetry element was found:   ',SYMM
         ELSE
            WRITE(LUPRI, '(/2A)')
     &         ' The following symmetry elements were found: ',SYMM
         END IF
      END IF

C ::: The reduced coordinates are moved back. :::
      CALL SRTOLD(MXCENT, ATMARR)
      IF (IPRSYM .GE. IPRTHR) THEN
         CALL PRNARR('Symmetry Independent Centres in original order',
     &      MXCENT, ATMARR)
      END IF
      NCTOT = NINT(ATMARR(1,0))
      IF (NCTOT .LE. 0) CALL QUIT ('No atoms defined')

      DO 250 I = 1, NCTOT
CMI    I don't know why is there this assignement but do this
CMI    only if no moving of coordinates required
         IF (.NOT.NOMOVE) THEN
           CORD(1:3,I) = ATMARR(1:3,I)
C        ... IOLD .ge. I after CALL SRTOLD
           IOLD = NINT(ATMARR(5,I))
         ELSE
           IOLD = I
         ENDIF
         IATOMTYP(I) = IATOMTYP(IOLD)
         CHARGE(I) = CHARGE(IOLD)
         IZATOM(I) = IZATOM(IOLD)
         ISOTOP(I) = ISOTOP(IOLD)
         GNUEXP(I) = GNUEXP(IOLD)
         NOORBT(I) = NOORBT(IOLD)         
#ifndef PRG_DIRAC
         RSPH(I)   = RSPH(IOLD)
         IDXSPH(I) = IDXSPH(IOLD)
         MULBSI(I) = MULBSI(IOLD)
#endif
         NCLINE(I) = NCLINE(IOLD)
         NAMN(I)   = NAMN(IOLD)
         NAMEX(3*I)   = NAMN(I)//' z'
         NAMEX(3*I-1) = NAMN(I)//' y'
         NAMEX(3*I-2) = NAMN(I)//' x'
 250  CONTINUE
C ::: Then the different elements are counted, and arrays updated. :::
      NONT(1:KATOM) = 0
      NONTYP = 1
      NONT(NONTYP) = 1
      DO I = 2, NCTOT
         IF (IATOMTYP(I) .NE. IATOMTYP(I-1)) THEN
            NONTYP = NONTYP + 1
         END IF
         NONT(NONTYP) = NONT(NONTYP) + 1
      END DO

      DO I = NCTOT + 1, KATOM
         NAMN(I) = '    '
         NAMEX(3*I)   = '      '
         NAMEX(3*I-1) = '      '
         NAMEX(3*I-2) = '      '
         CORD(1:3,I) = D0
         CHARGE(I) = D0
         IZATOM(I) = 0
         ISOTOP(I) = 1
      END DO
C ::: We update the coordinates of the input file :::
      if (internal_file) then
         DO 400 I = 1, NCTOT
            NCI   = NCLINE(I)
            TMPLN = MLINE(NCI)
            IF (ABS(CORD(1,I)).ge.10000.0D0 .OR.
     &          ABS(CORD(2,I)).ge.10000.0D0 .OR.
     &          ABS(CORD(3,I)).ge.10000.0D0) THEN
               WRITE(TMPLN(5:65),'(3G20.12,1X)')
     &            CORD(1,I),CORD(2,I),CORD(3,I)
            ELSE
               WRITE(TMPLN(5:65),'(3F20.14,1X)')
     &            CORD(1,I),CORD(2,I),CORD(3,I)
            END IF
C
            IF (ISOTOP(I) .GT. 1) THEN
              QMASS = DISOTP(IZATOM(I),ISOTOP(I),'MASS')
              WRITE(TMPLN(66:76),'(A8,I3)') 'Isotope=', NINT(QMASS)
            END IF
            MLINE(NCI) = TMPLN
 400     CONTINUE
C ::: Finally the number of atoms are updated. :::
         I = 1
         DO K = 1, NONTYP
            NCI   = NCLINE(I) - 1
           !TMPLN = MLINE(NCI)
           !CALL UPCASE(TMPLN)

               WRITE(MLINE(NCI)(12:15),'(I4)') NONT(K)
            I = I + NONT(K)
         END DO
      end if ! internal_file
C     IF (SOLVNT .OR. DORLM) THEN
C        NUCIND = NUCIND + 1
C        NCNTCV = NUCIND
C        NAMN(NUCIND) = 'cav '
C        NCLINE(NUCIND) = 0
C        NAMEX(3*NUCIND)     = NAMN(NUCIND)//' z'
C        NAMEX(3*NUCIND - 1) = NAMN(NUCIND)//' y'
C        NAMEX(3*NUCIND - 2) = NAMN(NUCIND)//' x'
C        CORD(1,NUCIND) = D0
C        CORD(2,NUCIND) = D0
C        CORD(3,NUCIND) = D0
C        CHARGE(NUCIND) = D0
C     END IF
      NUCIND = NCTOT
      RETURN
      END  ! subroutine FIND_SYMMETRY
C  /* Deck prnarr */
      SUBROUTINE PRNARR(TEXT, MXM, ARR)
C ::::::::::::::::::::::::::::::::::::::::::::::
C ::                                          ::
C ::       Prints all elements of array       ::
C ::(VB, part of hersym.F)                    ::
C ::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) TEXT
      DIMENSION ARR(6,0:MXM)
      I = LEN(TEXT)
 10   CONTINUE
         IF ((I .GE. 1) .AND. (TEXT(I:I) .EQ. ' ')) THEN
            I = I - 1
            GOTO 10
         END IF
      IF (I .GT. 0) WRITE(LUPRI,'(/1X,A/1X,80A1)') TEXT(1:I),('-',J=1,I)
      DO I=1,NINT(ARR(1,0))
         WRITE (LUPRI,11) NINT(ARR(4,I)),(ARR(J,I), J=1,3),
     &                    NINT(ARR(6,I))
      END DO
   11 FORMAT (I8, ' : ', 3F15.8,'  Isotope',I3)
      RETURN
      END
C  /* Deck prnmir */
      SUBROUTINE PRNMIR(MXM, MIRPLN)
C ::::::::::::::::::::::::::::::::::::::::::::::
C ::                                          ::
C ::       Prints all mirror planes           ::
C ::       (based on PRNARR)                  ::
C ::                                          ::
C ::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
      REAL*8    MIRPLN(6,0:MXM)

      WRITE(LUPRI,'(6(/A))')
     &   ' Mirror planes',
     &   ' -----------------------------',
     &   ' #1 = Type: 2=h, 1=v, 0=other',
     &   ' #2 = Number of atoms in plane',
     &   ' #3-5 = normal vector of plane',
     &   ' -----------------------------'
      DO I=1,NINT(MIRPLN(1,0))
         IT = NINT(MIRPLN(4,I))
         N  = MOD(IT,100000) 
         IT = IT / 100000
         WRITE (LUPRI,11) IT,N,(MIRPLN(J,I), J=1,3)
      END DO
   11 FORMAT (I2,I6, ' : ', 3F15.8)
      RETURN
      END
C  /* Deck prnarr2 */
      SUBROUTINE PRNARR2(TEXT, MXM, ARR, IARR)
C ::::::::::::::::::::::::::::::::::::::::::::::
C ::                                          ::
C ::       Prints all elements of array       ::
C ::       according to ordering in IARR      ::
C ::       (by T.Saue, based on PRNARR)       ::
C ::                                          ::
C ::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
      CHARACTER*(*) TEXT
      CHARACTER*40 UNDER
      DIMENSION ARR(6,0:MXM),IARR(MXM)
      I = LEN(TEXT)
 10   CONTINUE
         IF ((I .GE. 1) .AND. (TEXT(I:I) .EQ. ' ')) THEN
            I = I - 1
            GOTO 10
         END IF
      IF (I .GT. 0) WRITE(LUPRI,'(/1X,A/1X,80A1)') TEXT(1:I),('-',J=1,I)
      DO 200 I=1,NINT(ARR(1,0))
         II = IARR(I)
         WRITE (LUPRI,11) NINT(ARR(4,II)),(ARR(J,II), J=1,3),
     &                    NINT(ARR(6,II))
  11     FORMAT (I5, T10,3F15.8, T60,I3)
 200  CONTINUE
      RETURN
      END
C  /* Deck vecang */
      FUNCTION VECANG(V1, V2)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                 ::
C ::       Finds the angle between two vectors       ::
C ::                                                 ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "pi.h"
#include "mxsymm.h"
         DIMENSION V1(3), V2(3)
         TEMP = (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))*
     &          (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
         IF (TEMP .GT. ZERTOL) THEN
            TEMP = (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))/SQRT(TEMP)
C ::: Parallel vectors might yield a value here slightly greater than one :::
C ::: ACOS is undefined for these values, so we have to round  them off,  :::
C ::: removing the "excess", accumulated numerical error.                 :::
C ::: Also, ACOS is extremely sensitive around +/-1, so we round off      :::
C ::: numbers close to these values
            IF ((1.0D0-ABS(TEMP)) .LT. ZERTOL) TEMP = ANINT(TEMP)
            VECANG = ACOS(TEMP)
         ELSE
            VECANG = 0.5D0*PI
         END IF
         RETURN
      END
C  /* Deck add2ar */
      SUBROUTINE ADD2AR(MXM, ARR, IORDER, VEC)
C ::::::::::::::::::::::::::::::::::
C ::                              ::
C ::       Add Vec to array       ::
C ::                              ::
C ::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
      DIMENSION ARR(6,0:MXM), VEC(3)
      IF (ABS(VEC(1)) .LT. TOLLRN) VEC(1) = 0.0D0
      IF (ABS(VEC(2)) .LT. TOLLRN) VEC(2) = 0.0D0
      IF (ABS(VEC(3)) .LT. TOLLRN) VEC(3) = 0.0D0
      VCNORM = SQRT(VEC(1)*VEC(1)+VEC(2)*VEC(2)+VEC(3)*VEC(3))
      IF (VCNORM .LT. TOLLRN) RETURN
      VEC(1) = VEC(1)/VCNORM
      VEC(2) = VEC(2)/VCNORM
      VEC(3) = VEC(3)/VCNORM
      Nvecs  = NINT(ARR(1,0))
C
Chjaaj Sep08: remove duplicates already here, to include as many axes
Chjaaj        as possible before Nvecs .ge. MXM.
      DO I = 1, Nvecs
         IF (IORDER .GT. ARR(1,I)) CYCLE ! want to find highest order axis
         J = 0
         IF (ABS(ARR(1,I) - VEC(1)) .LT. TOLLRN) J = J + 1
         IF (ABS(ARR(2,I) - VEC(2)) .LT. TOLLRN) J = J + 1
         IF (ABS(ARR(3,I) - VEC(3)) .LT. TOLLRN) J = J + 1
         IF (J .EQ. 3) THEN
Cdbg        WRITE(LUPRI,*) 'ADD2AR info   : duplicate axis',I,N+1
            RETURN
         END IF
      END DO
C
      IF (Nvecs .GE. MXM) THEN
Cdbg     WRITE(LUPRI,*) 'ADD2AR warning: axis not added',Nvecs+1,MXM
Cdbg     WRITE(LUPRI,*) 'ADD2AR warning: the axis',VEC
         RETURN
      ELSE
Cdbg     WRITE(LUPRI,*) 'ADD2AR info   : axis added',Nvecs+1,MXM
Cdbg     WRITE(LUPRI,*) 'ADD2AR info   : the axis',VEC
      END IF
C ::: The vector is placed in its correct position in the array, :::
C ::: so no sorting is necessary.                                :::
      I = 1
  100 CONTINUE
         IF ((NINT(ARR(4,I)) .GT. IORDER) .AND. (I .LE. Nvecs)) THEN
            I = I + 1
            GOTO 100
         END IF
      DO J = Nvecs  , I, -1
         ARR(1,J+1) = ARR(1,J)
         ARR(2,J+1) = ARR(2,J)
         ARR(3,J+1) = ARR(3,J)
         ARR(4,J+1) = ARR(4,J)
         ARR(5,J+1) = ARR(5,J)
      END DO
      ARR(1,I) = VEC(1)
      ARR(2,I) = VEC(2)
      ARR(3,I) = VEC(3)
      ARR(4,I) = IORDER
      ARR(1,0) = Nvecs + 1
C
      RETURN
      END
C  /* Deck dldpax */
      SUBROUTINE DLDPAX(MXM, AXARR)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                             ::
C ::       Deletes all duplicate axes (parallel) in array        ::
C ::                                                             ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
      DIMENSION AXARR(6,0:MXM), V1(3), V2(3)

      N_ATOMS = NINT(AXARR(1,0))
      I = 1
 100  CONTINUE
            V1(1) = AXARR(1,I)
            V1(2) = AXARR(2,I)
            V1(3) = AXARR(3,I)
            J = I + 1
 150        CONTINUE
            IF (J .LE. N_ATOMS) THEN
C ::: Parallel axes yields a zero vector as cross product.  :::
C ::: Since the array is sorted, the last one will have     :::
C ::: equal or lesser order, and should be removed.         :::
               V2(1) =  V1(2)*AXARR(3,J) - AXARR(2,J)*V1(3)
               V2(2) = -V1(1)*AXARR(3,J) + AXARR(1,J)*V1(3)
               V2(3) =  V1(1)*AXARR(2,J) - AXARR(1,J)*V1(2)
               V2NORM2 = V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
               IF (V2NORM2 .LT. TOLLRN_2) THEN
                  DO 200 I = 1, N_ATOMS - J
                     AXARR(1,J+I-1) = AXARR(1,J+I)
                     AXARR(2,J+I-1) = AXARR(2,J+I)
                     AXARR(3,J+I-1) = AXARR(3,J+I)
                     AXARR(4,J+I-1) = AXARR(4,J+I)
 200              CONTINUE
                  N_ATOMS = N_ATOMS - 1
               ELSE
                  J = J + 1
               END IF
               GOTO 150
            END IF
         I = I + 1
      IF (I .LT. N_ATOMS) GOTO 100

      AXARR(1,0) = N_ATOMS
      RETURN
      END
C  /* Deck dldpat */
      SUBROUTINE DLDPAT(MXM, ATMARR)
C :::::::::::::::::::::::::::::::::::::::::::::
C ::                                         ::
C ::       Deletes all duplicate atoms       ::
C ::                                         ::
C :::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), V1(3), V2(3)
         N_ATOMS = NINT(ATMARR(1,0))
         CALL SRTATM(MXM, ATMARR, .TRUE.)
         I = 1
 100     CONTINUE
         IF (I .LT. N_ATOMS) THEN
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            J = I + 1
 150        CONTINUE
            IF ((J .LE. N_ATOMS) .AND.
     &        (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J)))) THEN
C ::: Subtracting two position yields a zero vector for duplicate atoms. :::
               V2(1) = ATMARR(1,J) - V1(1)
               V2(2) = ATMARR(2,J) - V1(2)
               V2(3) = ATMARR(3,J) - V1(3)
               IF (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3) .LT.
     &                                           TOLLRN_2) THEN
                  DO 200 K = 1, N_ATOMS - J
                     ATMARR(1,J+K-1) = ATMARR(1,J+K)
                     ATMARR(2,J+K-1) = ATMARR(2,J+K)
                     ATMARR(3,J+K-1) = ATMARR(3,J+K)
                     ATMARR(4,J+K-1) = ATMARR(4,J+K)
                     ATMARR(5,J+K-1) = ATMARR(5,J+K)
                     ATMARR(6,J+K-1) = ATMARR(6,J+K)
 200              CONTINUE
                  N_ATOMS = N_ATOMS - 1
               ELSE
                  J = J + 1
               END IF
               GOTO 150
            END IF
            I = I + 1
            GOTO 100
         END IF
         ATMARR(1,0) = 1.0D0*N_ATOMS
         RETURN
      END
C  /* Deck srtatm */
      SUBROUTINE SRTATM(IUPPER, ARR, USENRM)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                      ::
C ::       Sorts atoms after atom number, isotope number,                 ::
C ::       and distance from origo.                                       ::
C ::       Rotational axes are sorted after order, the norm does not      ::
C ::       matter for these (USENRM turns on and off the norm criteria).  ::
C ::       The method used is bubble sorting.                             ::
C ::                                                                      ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
         DIMENSION ARR(6,0:IUPPER)
         LOGICAL USENRM, SORTED
C
         N_ATOMS = NINT(ARR(1,0))
         I = 1
         SORTED = .FALSE.
 100     CONTINUE
         IF ((I .LT. N_ATOMS) .AND. (.NOT. SORTED)) THEN
            SORTED = .TRUE.
            DO 200 J = 1, N_ATOMS-I
               IF (USENRM) THEN
C              ... sort atoms
C              test 1) higher charge before lower charge
                  IF ( ARR(4,J+1) .GT. ARR(4,J) )  GOTO 300
                  IF ( ARR(4,J+1) .LT. ARR(4,J) )  GOTO 310
C              test 2) for same charge: higher isotope before lower isotope
                  IF ( ARR(6,J+1) .GT. ARR(6,J) ) GOTO 300
                  IF ( ARR(6,J+1) .LT. ARR(6,J) ) GOTO 310
                  RATM1 = ARR(1,J  )**2 + ARR(2,J  )**2 + ARR(3,J  )**2
                  RATM2 = ARR(1,J+1)**2 + ARR(2,J+1)**2 + ARR(3,J+1)**2
C              test 3) shorter distance to origo before longer distance
                  IF ( RATM1 .GT. RATM2 ) GOTO 300
                  GOTO 310 ! J and J+1 ar in correct order
C
  300             CONTINUE
C                 ... J and J+1 are NOT in correct order, switch them
                     DO K = 1, 6
                        TEMP       = ARR(K,J)
                        ARR(K,J)   = ARR(K,J+1)
                        ARR(K,J+1) = TEMP
                     END DO
                     SORTED = .FALSE.
  310             CONTINUE
               ELSE
C                 ... sort axes
                  IF (ARR(4,J+1) .GT. ARR(4,J)) THEN
C                 ... J and J+1 are NOT in correct order, switch them
                     DO K = 1, 5
                        TEMP       = ARR(K,J)
                        ARR(K,J)   = ARR(K,J+1)
                        ARR(K,J+1) = TEMP
                     END DO
                     SORTED = .FALSE.
                  END IF
               END IF
 200        CONTINUE
            I = I + 1
            GOTO 100
         END IF
         RETURN
      END
C  /* Deck srtold */
      SUBROUTINE SRTOLD(IUPPER, ARR)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                               ::
C ::       Sorts atoms after original order. The method used       ::
C ::       is bubble sorting.                                      ::
C ::                                                               ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
         DIMENSION ARR(6,0:IUPPER)
         LOGICAL SORTED
         N_ATOMS = NINT(ARR(1,0))
         I = 1
         SORTED = .FALSE.
 100     CONTINUE
         IF ((I .LT. N_ATOMS) .AND. (.NOT. SORTED)) THEN
            SORTED = .TRUE.
            DO 200 J = 1, N_ATOMS-I
               IF (NINT(ARR(5,J)) .GT. NINT(ARR(5,J+1))) THEN
                  DO 300 K = 1, 6
                     TEMP = ARR(K,J)
                     ARR(K,J) = ARR(K,J+1)
                     ARR(K,J+1) = TEMP
 300              CONTINUE
                  SORTED = .FALSE.
               END IF
 200        CONTINUE
            I = I + 1
            GOTO 100
         END IF
         RETURN
      END
C ---
      SUBROUTINE CENTER_MOLECULE(ATMARR)
C :::::::::::::::::::::::::::::::::::::::::::
C ::                                       ::
C ::       Centres molecule in space       ::
C ::                                       ::
C :::::::::::::::::::::::::::::::::::::::::::
C
#include "implicit.h"
#include "priunit.h"
      DIMENSION ATMARR(6,0:*), CMXYZ(3)
      DIMENSION CINRTP(6),VPCMOM(3,3), TMP(3), ITMP(3)

! infpar.h : MYTID, MPARID
#include "infpar.h"

      N_ATOMS = NINT(ATMARR(1,0))

#if HERSYM_DEBUG > 0
         WRITE (LUPRI,*) 'CENTER_MOLECULE: initial coord.s etc.'
         CALL OUTPUT(ATMARR(1,0),1,6,1,N_ATOMS+1,6,N_ATOMS+1,1,LUPRI)
         CALL PRNARR('CENTER_MOLECULE: ATMARR start', N_ATOMS, ATMARR)
#endif

      IF (N_ATOMS .EQ. 1) THEN ! atom, just move to (0,0,0)
         ATMARR(1,1) = 0.0D0
         ATMARR(2,1) = 0.0D0
         ATMARR(3,1) = 0.0D0
         RETURN
      END IF

      CMXYZ(1:3) = 0.0D0
C ::: Finds the three coordinates of the center of mass,
C ::: then we calculate the inertia tensor around the center of mass,
C ::: then the molecule is moved so that the center of mass coincides with the origin
C ::: and the molecule is rotated such that the principal moment of inertia axes
C ::: are along the coordinate axes.
      TOTMAS = 0.0D0
      DO I = 1,N_ATOMS
         IZ = NINT(ATMARR(4,I))
      IF (IZ .LE. 0) CYCLE
C... hjaaj exclude point charges, MM atoms and floating orbitals
C       in the centering of the molecule.
         AMASS  = DISOTP(IZ,NINT(ATMARR(6,I)),'MASS')
         TOTMAS = TOTMAS + AMASS
         CMXYZ(1) = CMXYZ(1) + ATMARR(1,I)*AMASS
         CMXYZ(2) = CMXYZ(2) + ATMARR(2,I)*AMASS
         CMXYZ(3) = CMXYZ(3) + ATMARR(3,I)*AMASS
      END DO
      CMXYZ(1) = CMXYZ(1)/TOTMAS
      CMXYZ(2) = CMXYZ(2)/TOTMAS
      CMXYZ(3) = CMXYZ(3)/TOTMAS
C     CMXYZ(1:3) now contains center-of-mass coordinates

#if HERSYM_DEBUG > 0
      write(lupri,*) 'Total mass ',TOTMAS
      write(lupri,*) 'Center of mass',CMXYZ(1:3)
#endif

C
C     Now we are ready to center molecule:
C
      DO I=1,N_ATOMS
         ATMARR(1,I) = ATMARR(1,I) - CMXYZ(1)
         ATMARR(2,I) = ATMARR(2,I) - CMXYZ(2)
         ATMARR(3,I) = ATMARR(3,I) - CMXYZ(3)
      END DO

#if HERSYM_DEBUG > 0
      write (lupri,*) 'CENTER_MOLECULE: transformed coord.s'
      CALL OUTPUT(ATMARR(1,0),1,6,1,N_ATOMS+1,6,N_ATOMS+1,1,LUPRI)
#endif

      IF (MYTID .EQ. MPARID) WRITE(LUPRI,'(/A)')
     &   '    The molecule has been centered at center of mass'
#if HERSYM_DEBUG > 0
      CALL PRNARR('CENTER_MOLECULE: coord.s now', N_ATOMS, ATMARR)
#endif
      RETURN
      END ! subroutine CENTER_MOLECULE
C ---
      SUBROUTINE PRINCIPLE_AXIS_MOLECULE(ATMARR)
C :::::::::::::::::::::::::::::::::::::::::::
C ::                                       ::
C ::  Rotate molecule to principle axes    ::
C ::  of inertia                           ::
C ::                                       ::
C :::::::::::::::::::::::::::::::::::::::::::
C
C   Revision Aug 2011 HJAaJ (merged some of my old code from 2002):
C    - rotate to axes corresponding to
C      principal moment of mass inertia
C    - do it twice to minimize round-off errors
C
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
      DIMENSION ATMARR(6,0:*), CMXYZ(3)
      DIMENSION CINRTP(6),VPCMOM(3,3), TMP(3), ITMP(3)

! infpar.h: MPARID, MYTID
#include "infpar.h"

      N_ATOMS = NINT(ATMARR(1,0))

#if HERSYM_DEBUG > 0
         WRITE (LUPRI,*) 'PRINCIPLE_AXIS_MOLECULE: initial coord.s etc.'
         CALL OUTPUT(ATMARR(1,0),1,6,1,N_ATOMS+1,6,N_ATOMS+1,1,LUPRI)
         CALL PRNARR('PRINCIPLE_AXIS_MOLECULE: ATMARR start',
     &      N_ATOMS, ATMARR)
#endif

      IF (N_ATOMS .EQ. 1) THEN ! atom, just move to (0,0,0)
         ATMARR(1,1) = 0.0D0
         ATMARR(2,1) = 0.0D0
         ATMARR(3,1) = 0.0D0
         GO TO 100
      END IF

      ITURN = 0
  10  ITURN = ITURN + 1
         CMXYZ(1:3) = 0.0D0
C ::: Finds the three coordinates of the center of mass,
C ::: then we calculate the inertia tensor around the center of mass,
C ::: then the molecule is moved so that the center of mass coincides with the origin
C ::: and the molecule is rotated such that the principal moment of inertia axes
C ::: are along the coordinate axes.
         TOTMAS = 0.0D0
         DO I = 1,N_ATOMS
            IZ = NINT(ATMARR(4,I))
         IF (IZ .LE. 0) CYCLE
C... hjaaj exclude point charges, MM atoms and floating orbitals
C          in the centering of the molecule.
            AMASS  = DISOTP(IZ,NINT(ATMARR(6,I)),'MASS')
            TOTMAS = TOTMAS + AMASS
            CMXYZ(1) = CMXYZ(1) + ATMARR(1,I)*AMASS
            CMXYZ(2) = CMXYZ(2) + ATMARR(2,I)*AMASS
            CMXYZ(3) = CMXYZ(3) + ATMARR(3,I)*AMASS
         END DO
         CMXYZ(1) = CMXYZ(1)/TOTMAS
         CMXYZ(2) = CMXYZ(2)/TOTMAS
         CMXYZ(3) = CMXYZ(3)/TOTMAS
C        CMXYZ(1:3) now contains center-of-mass coordinates

#if HERSYM_DEBUG > 0
         write(lupri,*) 'Total mass ',TOTMAS
         write(lupri,*) 'Center of mass',CMXYZ(1:3)
#endif

         CINRTP(1:6) = 0.0D0
         DO I=1,N_ATOMS
            IZ = NINT(ATMARR(4,I))
         IF (IZ .LE. 0) CYCLE
            AMASS  = DISOTP(IZ,NINT(ATMARR(6,I)),'MASS')
            XI = ATMARR(1,I) - CMXYZ(1)
            YI = ATMARR(2,I) - CMXYZ(2)
            ZI = ATMARR(3,I) - CMXYZ(3)
            CINRTP(1) = CINRTP(1) + AMASS*(YI*YI + ZI*ZI)
            CINRTP(2) = CINRTP(2) - AMASS* XI*YI
            CINRTP(3) = CINRTP(3) + AMASS*(XI*XI + ZI*ZI)
            CINRTP(4) = CINRTP(4) - AMASS* XI*ZI
            CINRTP(5) = CINRTP(5) - AMASS* YI*ZI
            CINRTP(6) = CINRTP(6) + AMASS*(XI*XI + YI*YI)
         END DO
#if HERSYM_DEBUG > 0
         write(lupri,*) 'mass inertia matrix'
         call outpak(cinrtp,3,1,lupri)
#endif
         CALL DUNIT(VPCMOM,3)
         CALL JACO(CINRTP,VPCMOM,3,3,3,TMP,ITMP)
#if HERSYM_DEBUG > 0
         call outpak(cinrtp,3,1,lupri)
#endif
         CINRTP(2) = CINRTP(3)
         CINRTP(3) = CINRTP(6)
         IF (CINRTP(1) - CINRTP(3) .LT. TOLLRN) THEN
         ! hjaaj Feb 2016: do not rotate input coordinates for T(d), O(h) and other
         ! high symmetry point groups where small changes can give big rotations
         ! because of (near-)degeneracy
            CALL DUNIT(VPCMOM,3)
         END IF
C        Check if molecule is mostly oblate or mostly prolate

         TMP(1:3) = CINRTP(1:3)
         CALL ORDER(DUMMY,TMP,3,0)
         IF (TMP(1)/TMP(2) .LT. 0.8D0*TMP(2)/TMP(3)) THEN ! prolate type
C        ... definitely more prolate than oblate, maybe (close to) linear,
C            sort after decreasing value to get molecular axis along z-axis
            IF (CINRTP(1) .LT. CINRTP(2)) THEN ! switch 1 and 2
               CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
               CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
               VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
            END IF
            IF (CINRTP(2) .LT. CINRTP(3)) THEN ! switch 2 and 3
               CALL DSWAP(3,VPCMOM(1,2),1,VPCMOM(1,3),1)
               CALL DSWAP(1,CINRTP(2),1,CINRTP(3),1)
               VPCMOM(1:3,2) = -VPCMOM(1:3,2)  ! to not change chirality!!!
               IF (CINRTP(1) .LT. CINRTP(2)) THEN ! switch 1 and 2
                  CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
                  CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
                  VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
               END IF
            END IF
         ELSE ! oblate type
C        ... Sorting after increasing value - will put oblate molecule
C            primarily in the x-y plane.
            IF (CINRTP(1) .GT. CINRTP(2)) THEN ! switch 1 and 2
               CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
               CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
               VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
            END IF
            IF (CINRTP(2) .GT. CINRTP(3)) THEN ! switch 2 and 3
               CALL DSWAP(3,VPCMOM(1,2),1,VPCMOM(1,3),1)
               CALL DSWAP(1,CINRTP(2),1,CINRTP(3),1)
               VPCMOM(1:3,2) = -VPCMOM(1:3,2)  ! to not change chirality!!!
               IF (CINRTP(1) .GT. CINRTP(2)) THEN ! switch 1 and 2
                  CALL DSWAP(3,VPCMOM(1,1),1,VPCMOM(1,2),1)
                  CALL DSWAP(1,CINRTP(1),1,CINRTP(2),1)
                  VPCMOM(1:3,1) = -VPCMOM(1:3,1)  ! to not change chirality!!!
               END IF
            END IF
         END IF

#if HERSYM_DEBUG > 0
         write (lupri,*) 'PCMOM, VPCMOM turn',ITURN
         call output(CINRTP,1,1,1,3,1,3,1,lupri)
         call output(VPCMOM,1,3,1,3,3,3,1,LUPRI)
#endif
C
C        Now we are ready to center and rotate molecule:
C        Note:  VPCMOM describes how to rotate coordinate axes
C               thus molecule should be rotated with VPCMOM(inverse)
C               which is same as VPC(transposed)
C
         DO I=1,N_ATOMS
            XI = ATMARR(1,I) - CMXYZ(1)
            YI = ATMARR(2,I) - CMXYZ(2)
            ZI = ATMARR(3,I) - CMXYZ(3)
            ATMARR(1,I) = XI*VPCMOM(1,1)+YI*VPCMOM(2,1)+ZI*VPCMOM(3,1)
            ATMARR(2,I) = XI*VPCMOM(1,2)+YI*VPCMOM(2,2)+ZI*VPCMOM(3,2)
            ATMARR(3,I) = XI*VPCMOM(1,3)+YI*VPCMOM(2,3)+ZI*VPCMOM(3,3)
         END DO

#if HERSYM_DEBUG > 0
         write (lupri,*) 'PRINCIPLE_AXIS_MOLECULE: transformed coord.s'
         CALL OUTPUT(ATMARR(1,0),1,6,1,N_ATOMS+1,6,N_ATOMS+1,1,LUPRI)
#endif

      IF (ITURN .LE. 1) GO TO 10

 100  IF (MYTID .EQ. MPARID) WRITE(LUPRI,'(/A/A)')
     &   '    The molecule has been centered at center of mass and',
     &   '    rotated so principal axes of inertia are along'//
     &      ' the coordinate axes.'
#if HERSYM_DEBUG > 0
      CALL PRNARR('PRINCIPLE_AXIS_MOLECULE: coord.s now',
     &   N_ATOMS, ATMARR)
#endif
      RETURN
      END ! subroutine PRINCIPLE_AXIS_MOLECULE
C  /* Deck rotmol */
      SUBROUTINE ROTMOL(MXM, ATMARR, AXVEC, DEG)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                         ::
C ::       Rotate molecule in space Deg degrees around       ::
C ::           an arbitrary axis defined by AxVec            ::
C ::                                                         ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), AXVEC(3), V1(3),V2(3),V3(3),V4(3)
         AXNRM=SQRT(AXVEC(1)*AXVEC(1) + AXVEC(2)*AXVEC(2) +
     &                                  AXVEC(3)*AXVEC(3))
         AXVEC(1) = AXVEC(1)/AXNRM
         AXVEC(2) = AXVEC(2)/AXNRM
         AXVEC(3) = AXVEC(3)/AXNRM
         CSD = COS(DEG)
         SND = SIN(DEG)
         DO 100 I=1,NINT(ATMARR(1,0))
            V1(1)=ATMARR(1,I)
            V1(2)=ATMARR(2,I)
            V1(3)=ATMARR(3,I)
            DOT = AXVEC(1)*V1(1)+AXVEC(2)*V1(2)+AXVEC(3)*V1(3)
            VCLEN = ABS(DOT) - SQRT(V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3))
C ::: Decomposes position to two vectors, one parallel(V2) and    :::
C ::: one perpendicular(V3) to the rotational axis after checking :::
C ::: wether the rotational axis passes through the atom          :::
            IF (ABS(VCLEN) .GT. ZERTOL) THEN
               V2(1) = (DOT*AXVEC(1))
               V3(1) = V1(1)-(DOT*AXVEC(1))
               V2(2) = (DOT*AXVEC(2))
               V3(2) = V1(2)-(DOT*AXVEC(2))
               V2(3) = (DOT*AXVEC(3))
               V3(3) = V1(3)-(DOT*AXVEC(3))
C ::: Finds another vector(V4) in the plane affected by rotation :::
               V4(1) =  AXVEC(2)*V3(3) - V3(2)*AXVEC(3)
               V4(2) = -AXVEC(1)*V3(3) + V3(1)*AXVEC(3)
               V4(3) =  AXVEC(1)*V3(2) - V3(1)*AXVEC(2)
               VCLEN = SQRT(V4(1)*V4(1)+V4(2)*V4(2)+V4(3)*V4(3))
               V4(1)= V4(1)/VCLEN
               V4(2)= V4(2)/VCLEN
               V4(3)= V4(3)/VCLEN
C ::: Stores the norm of the perpendicular component :::
C ::: (which is to be rotated).                      :::
               VCLEN = SQRT(V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
               V3(1)= V3(1)/VCLEN
               V3(2)= V3(2)/VCLEN
               V3(3)= V3(3)/VCLEN
C ::: V3 and V4 are now normalized and spans a plan perpendicular to the :::
C ::: rotational axis. A normalized vector is then rotated the desired   :::
C ::: amount, before it is multiplied with the stored norm, L.           :::
C ::: Finally the rotated coordinates are given back to the atom. :::
               ATMARR(1,I) = V2(1) + VCLEN*(CSD*V3(1)+SND*V4(1))
               ATMARR(2,I) = V2(2) + VCLEN*(CSD*V3(2)+SND*V4(2))
               ATMARR(3,I) = V2(3) + VCLEN*(CSD*V3(3)+SND*V4(3))
            END IF
 100     CONTINUE
         RETURN
      END
C  /* Deck turn_molecule */
      SUBROUTINE TURN_MOLECULE(MXM, ATMARR, DRTAXS, DIRAXS, DMRPLN,
     &                         TMPARR, LENTMP)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                    ::
C ::       Turns molecule to sensible orientation       ::
C ::                                                    ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
      DIMENSION ATMARR(6,0:MXM), DRTAXS(6,0:MAXAXS)
      DIMENSION DIRAXS(6,0:MAXAXS), DMRPLN (6,0:MAXMIR)
      DIMENSION V1(3), V2(3), V3(3), TMPARR(6,0:LENTMP)
      LOGICAL   LINMOL

      IPRTHR = 3 ! Miro: for runtime check it must be defined as variable
C
      N_ROTAXS = NINT(DRTAXS(1,0))
      N_MIRROR = NINT(DMRPLN(1,0))
      N_ATOMS  = NINT(ATMARR(1,0))

      IF (IPRSYM .GT. IPRTHR+1) THEN
         WRITE(LUPRI,*) 'TURN_MOLECULE: N_ROTAXS, N_MIRROR, N_ATOMS =',
     &      N_ROTAXS, N_MIRROR, N_ATOMS
      END IF

      IF (N_ROTAXS .LE. 0 . AND. .NOT. LINMOL(MXM, ATMARR)) RETURN ! RETURN, no turning if no rotation axes (C1, Cs, or Ci symmetry)
C
      CALL DZERO(TMPARR,6*(LENTMP+1))
C ::: First we make a new array of vectors, consisting of mirror :::
C ::: planes, even-ordered rotational axes and atom positions.   :::
      ITURN = 0
c     NTOT  = N_ATOMS + N_ROTAXS + N_MIRROR
C        ... we count instead
      NTOT  = 0
C 1) rotational axes
chj Nov 2007: include both even and odd orders
      IF (N_ROTAXS .EQ. 7 .AND. NINT(DRTAXS(4,1)) .EQ. 3) THEN
      ! but for T(d) we do not want the highest rot. axis along z-axis
         XXX = 2.0D0 ! select only even-ordered rotations
      ELSE
         XXX = 1.0D0 ! select all rotation
      END IF
      DO I = 1, N_ROTAXS
      IF (ABS(MOD(DRTAXS(4,I),XXX)) .LT. 1.0D-8) THEN
         NTOT = NTOT + 1
         DO K = 1, 4
            TMPARR(K,NTOT) = DRTAXS(K,I)
         END DO
chj      TMPARR(4,NTOT) = 1.0D0
chj: no, we want to save the order to find highest rotation axis for z-direction
      END IF
      END DO
C 2) mirror planes
      IF (N_ROTAXS .NE. 1 .OR. N_MIRROR .NE. 2)
     &THEN
      ! do not use the mirror planes for C2v molecules (e.g. H2O)
      ! we want to be sure that the molecule ends up in the yz-plane
      DO I = 1, N_MIRROR
         NTOT = NTOT + 1
         DO K = 1, 3
            TMPARR(K,NTOT) = DMRPLN(K,I)
         END DO
         TMPARR(4,NTOT) = 1.0D0
      END DO
      END IF
C 3) and all the nuclei/point charges
      DO I = 1, N_ATOMS
         NTOT = NTOT + 1
         DO K = 1, 3
            TMPARR(K,NTOT) = ATMARR(K,I)
         END DO
         TMPARR(4,NTOT) = 0.0D0
         TMPARR(6,NTOT) = 1.0D0
      END DO
C
      TMPARR(1,0) = NTOT
C
C     Find first non-zero vector, place it along z-axis
C     (if no symmetry axes, it will become the coordinates of
C     the first of the nuclei with highest charge).
C
      I = 1
      V1(1) = TMPARR(1,1)
      V1(2) = TMPARR(2,1)
      V1(3) = TMPARR(3,1)
 150  CONTINUE
         IF (((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)) .LT. TOLLRN_2)
     &                        .AND. (I .LT. NTOT)) THEN
            I = I + 1
            V1(1) = TMPARR(1,I)
            V1(2) = TMPARR(2,I)
            V1(3) = TMPARR(3,I)
            GOTO 150
         END IF
         IF ((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)).GT.TOLLRN_2)THEN
C           V2 = target direction (along z-axis)
            V2(1) = 0.0D0
            V2(2) = 0.0D0
            V2(3) = 1.0D0
C           V3 = axis to rotate around
            V3(1) =  V1(2)
            V3(2) = -V1(1)
            V3(3) = 0.0D0
            IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                   .GT. TOLLRN_2) THEN
               DEG = VECANG(V1, V2)
               CALL ROTMOL(MXM   , ATMARR, V3, DEG)
               CALL ROTMOL(MAXAXS, DRTAXS, V3, DEG)
               CALL ROTMOL(MAXAXS, DIRAXS, V3, DEG)
               CALL ROTMOL(MAXMIR, DMRPLN, V3, DEG)
               CALL ROTMOL(NTOT  , TMPARR, V3, DEG)
            END IF
            ITURN = 1
         ELSE
            ITURN = 2
         END IF
C
C        Find second non-zero vector, rotate it so x-coordinate zero
C        (place planar C2v molecule in yz-plane, following :
C           R. S. Mulliken, J. Chem. Phys. 23, 1997 (1955);
C           http://dx.doi.org/10.1063/1.1740655 )
C        First look for a vector with zero z-coordinate
C
         JJ = I + 1
 200     CONTINUE
         IF ((ITURN .LT. 2) .AND. (I .LE. NTOT)) THEN
            IBEF = I
 250        CONTINUE
            IF (I .GT. NTOT) THEN
               I = IBEF
            ELSE IF (ABS(TMPARR(3,I)) .GT. 1.0D-8) THEN
               I = I + 1
               GOTO 250
            END IF
            V1(1) = TMPARR(1,I)
            V1(2) = TMPARR(2,I)
            V1(3) = 0.0D0
            IF ((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .GT. TOLLRN_2)
     &          .AND. (ABS(TMPARR(3,I)) .LT. TOLLRN)) THEN
C              V2 = target direction (along y-axis)
               V2(1) = 0.0D0
               V2(2) = 1.0D0
               V2(3) = 0.0D0
C              V3 = axis to rotate around
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                                        .GT. TOLLRN_2) THEN
                  DEG = VECANG(V1, V2)
                  CALL ROTMOL(MXM   , ATMARR, V3, DEG)
                  CALL ROTMOL(MAXAXS, DRTAXS, V3, DEG)
                  CALL ROTMOL(MAXAXS, DIRAXS, V3, DEG)
                  CALL ROTMOL(MAXMIR, DMRPLN, V3, DEG)
               END IF
               ITURN = 2
            END IF
            I = I + 1
            GOTO 200
         END IF
         IF (ITURN .EQ. 1) THEN
C           did not find a second non-zero vector with zero z-coordinate
C           - zero z-coordinate in temporary array and try again,
C             i.e. do not require true z-coordinate to be zero this time.
            DO 300 K = 1, NTOT
               TMPARR(3,K) = 0.0D0
 300        CONTINUE
            ITURN = -1
            I = JJ
            GOTO 200
         END IF
C
      END
C  /* Deck linmol */
      LOGICAL FUNCTION LINMOL(MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                              ::
C ::       Checks if the molecule is linear       ::
C ::                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), V1(3), V2(3), V3(3)
         N_ATOMS = NINT(ATMARR(1,0))
         LINMOL = .TRUE.
         IF (N_ATOMS .GT. 2) THEN
C ::: Builds first vector from bond between the two first atoms :::
            V1(1) = ATMARR(1,2)-ATMARR(1,1)
            V1(2) = ATMARR(2,2)-ATMARR(2,1)
            V1(3) = ATMARR(3,2)-ATMARR(3,1)
            DO I = 3, N_ATOMS
C ::: If molecule is to be linear, all subsequent bonds must be parallel  :::
C ::: with the first one, and hence give a zero vector as cross product.  :::
               V2(1) = ATMARR(1,I)-ATMARR(1,1)
               V2(2) = ATMARR(2,I)-ATMARR(2,1)
               V2(3) = ATMARR(3,I)-ATMARR(3,1)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                   .GT. TOLLRN_2) LINMOL = .FALSE.
            END DO
         END IF
         RETURN
      END
C  /* Deck plnmol */
      LOGICAL FUNCTION PLNMOL(MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                              ::
C ::       Checks if the molecule is planar       ::
C ::                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), V1(3), V2(3), V3(3)
         LOGICAL LINMOL
         N_ATOMS = NINT(ATMARR(1,0))
         PLNMOL = .TRUE.
C ::: With three atoms or less, we can be sure it is planar. :::
C ::: If it is linear, it is also planar of course.          :::
         IF ((N_ATOMS .GT. 3) .AND. (.NOT. LINMOL(MXM, ATMARR))) THEN
C ::: We then have to find two bonds that span the possible molecular :::
C ::: plane. Since we have taken care of linear molecules, two such   :::
C ::: vectors do indeed exist. The cross product, V3, of these two    :::
C ::: vectors will be perpendicular to the plane.                     :::
            V1(1) = ATMARR(1,1) - ATMARR(1,2)
            V2(1) = ATMARR(1,2) - ATMARR(1,3)
            V1(2) = ATMARR(2,1) - ATMARR(2,2)
            V2(2) = ATMARR(2,2) - ATMARR(2,3)
            V1(3) = ATMARR(3,1) - ATMARR(3,2)
            V2(3) = ATMARR(3,2) - ATMARR(3,3)
            V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
            V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
            V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
            I = 3
 100        CONTINUE
            IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3) .LT. TOLLRN_2)
     &          .AND. (I .LE. N_ATOMS)) THEN
               V2(1) = ATMARR(1,2) - ATMARR(1,I)
               V2(2) = ATMARR(2,2) - ATMARR(2,I)
               V2(3) = ATMARR(3,2) - ATMARR(3,I)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               I = I + 1
               GOTO 100
            END IF
C ::: We the go through the rest of the bonds, verifing that they    :::
C ::: are perpendicular to the perpendicular vector, V3, that is,    :::
C ::: that they lie in the same plane. The test is done by checking  :::
C ::: that the dot product is very close to zero.                    :::
            I = 2
 200        CONTINUE
            IF ((I .LT. N_ATOMS) .AND. PLNMOL) THEN
               V1(1) = ATMARR(1,I) - ATMARR(1,I+1)
               V1(2) = ATMARR(2,I) - ATMARR(2,I+1)
               V1(3) = ATMARR(3,I) - ATMARR(3,I+1)
               PLNMOL = (ABS(V1(1)*V3(1)+V1(2)*V3(2)+V1(3)*V3(3))
     &                                                  .LT. TOLLRN)
               I = I + 1
               GOTO 200
            END IF
         END IF
         RETURN
      END
C  /* Deck invmol */
      LOGICAL FUNCTION INVMOL(MXM, ATMARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                              ::
C ::       Checks if the molecule has a center of inversion       ::
C ::                                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
#include "mxcent.h"
         DIMENSION ATMARR(6,0:MXM), V1(3)
         DIMENSION IUSED(MXCENT)
         N_ATOMS = NINT(ATMARR(1,0))
         IUSED(1:MXCENT) = 0
C ::: The routine depends on the atoms being sorted and centered.    :::
C ::: The matrix 'Used' keeps track of which atoms have been paired. :::
C ::: First, atom(s) lying in origo are dropped.                     :::
         DO 100 I=1,N_ATOMS
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .LT. TOLLRN_2) THEN
               IUSED(I) = 1
            ELSE
               IUSED(I) = 0
            END IF
 100     CONTINUE
C ::: We then proceed with the rest of the atoms, trying to pair them. :::
         I = 1
 200     CONTINUE
         IF (I .LT. N_ATOMS) THEN
            IF (IUSED(I) .EQ. 0) THEN
               J = I + 1
 250           CONTINUE
                  IF ((IUSED(I) .EQ. 0) .AND. (J .LE. N_ATOMS)) THEN
                     IF (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J)) .AND.
     &                   NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J))) THEN
                        V1(1) = ATMARR(1,I) + ATMARR(1,J)
                        V1(2) = ATMARR(2,I) + ATMARR(2,J)
                        V1(3) = ATMARR(3,I) + ATMARR(3,J)
C ::: The sum of two anti-parallel vectors of equal length is the   :::
C ::: zero vector. We use this to find each atoms "partner".        :::
                        IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .LT.
     &                       TOLLRN_2) THEN
                           IUSED(I) = 1
                           IUSED(J) = 1
                        END IF
                        J = J + 1
                        GOTO 250
                     END IF
                  END IF
               END IF
            I = I + 1
            GOTO 200
         END IF
C ::: Finally, we have to check that all the atoms of the molecule :::
C ::: have been inverted into another atom of the same kind.       :::
         ISUM = 0
         DO 300 I=1,N_ATOMS
            ISUM = ISUM + IUSED(I)
 300     CONTINUE
         INVMOL = (ISUM .EQ. N_ATOMS)
         RETURN
      END
C  /* Deck mirmol */
      LOGICAL FUNCTION MIRMOL(MXM,ORGARR,ATMARR,VEC,N_ATOMS_IN_PLANE)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                   ::
C ::       Checks if the mol. has mirror plane defined by vector       ::
C ::                                                                   ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
Chj Oct 2007: new parameter N_ATOMS_IN_PLANE
C             used for ROTMOL to make sure the convention is followed
C             that as many atoms as possible are placed in the x-z plane.
C             (After z-axis has been chosen as one of the axes of highest
C              order)
C
#include "implicit.h"
#include "mxcent.h"
#include "mxsymm.h"
         DIMENSION ORGARR(6,0:MXM), ATMARR(6,0:MXM), VEC(3)
         DIMENSION VROT(3), V1(3), IUSED(MXCENT)
         N_ATOMS = NINT(ORGARR(1,0))
C ::: Instead of working on the original coordinates, we make a copy.   :::
C ::: By doing so, no numerical error is built up by repeated rotation. :::
         DO 100 I = 0,N_ATOMS
            ATMARR(1,I) = ORGARR(1,I)
            ATMARR(2,I) = ORGARR(2,I)
            ATMARR(3,I) = ORGARR(3,I)
            ATMARR(4,I) = ORGARR(4,I)
            ATMARR(6,I) = ORGARR(6,I)
 100     CONTINUE
C ::: The molecule is rotated so that the possible mirror plane lies :::
C ::: in the XY-plane.                                               :::
         V1(1) = 0.0D0
         V1(2) = 0.0D0
         V1(3) = 1.0D0
         VROT(1) =  V1(2)*VEC(3) - VEC(2)*V1(3)
         VROT(2) = -V1(1)*VEC(3) + VEC(1)*V1(3)
         VROT(3) =  V1(1)*VEC(2) - VEC(1)*V1(2)
         IF (VROT(1)*VROT(1)+VROT(2)*VROT(2)+VROT(3)*VROT(3) .GT.
     &                                            ZERTOL*ZERTOL) THEN
            DEG = VECANG(V1, VEC)
            CALL ROTMOL(MXM, ATMARR, VROT, -DEG)
         END IF
C ::: All atoms lying in the XY-plane are dropped. :::
         N_ATOMS_IN_PLANE = 0
         DO 200 I=1,N_ATOMS
            IF (ATMARR(3,I)*ATMARR(3,I).LT.TOLLRN_2) THEN
               N_ATOMS_IN_PLANE = N_ATOMS_IN_PLANE + 1
               IUSED(I) = 1
            ELSE
               IUSED(I) = 0
            END IF
 200     CONTINUE
C ::: Then we go through the rest, trying to find each atoms mirror image :::
         I = 1
 300     CONTINUE
         IF (I .LT. N_ATOMS) THEN
            IF (IUSED(I) .EQ. 0) THEN
               J = I + 1
 350           CONTINUE
!radovan: test dft_dirac_ao_eval goes out of bounds
!Array reference out of bounds for array 'atmarr', upper bound of dimension 2 exceeded (3 > 2)
               IF ((IUSED(I) .EQ. 0) .AND. (J .LE. N_ATOMS) .AND.
     &             (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J))) .AND.
     &             (NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J)))) THEN
C ::: The position of the possible mirror atom is investigated.      :::
C ::: First, we subtract the XY-comp. of the first atom. If the two  :::
C ::: project into the same point in the XY-plane, Atom #2 only has  :::
C ::: its Z-comp. left. We add the Z-comp. of #1, and the sum should :::
C ::: be a zero vector if they are each others mirror image.         :::
                  V1(1) = ATMARR(1,J) - ATMARR(1,I)
                  V1(2) = ATMARR(2,J) - ATMARR(2,I)
                  V1(3) = ATMARR(3,J) + ATMARR(3,I)
                  IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3) .LT.
     &                                    10.D0*TOLLRN_2) THEN
                     IUSED(I) = 1
                     IUSED(J) = 1
                  END IF
                  J = J + 1
                  GOTO 350
               END IF
            END IF
            I = I + 1
            GOTO 300
         END IF
C ::: Finally, we have to check that all atoms have been used. :::
         ISUM = 0
         DO 400 I = 1, N_ATOMS
            ISUM = ISUM + IUSED(I)
 400     CONTINUE
         MIRMOL = (ISUM .EQ. N_ATOMS)
         RETURN
      END
C  /* Deck rotaxs */
      LOGICAL FUNCTION ROTAXS(MXM,ORGARR,ATMARR,IORDER,ROTVEC,PROPER)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                ::
C :: PROPER true: Check if the molecule has a proper rotational     ::
C ::       axis of order IORDER around the vector ROTVEC            ::
C ::                                                                ::
C :: PROPER false: Check if the molecule has an improper rotational ::
C ::       axis of order 2*IORDER around the vector ROTVEC          ::
C ::                                                                ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "pi.h"
#include "mxcent.h"
#include "mxsymm.h"
      DIMENSION ORGARR(6,0:MXM), ATMARR(6,0:MXM), ROTVEC(3)
      LOGICAL   PROPER, EQL
      DIMENSION VROT(3), V1(3), V2(3)
      DIMENSION IUSED(MXCENT)

      IF (IORDER .LE. 1) thEN
         write (lupri,*) 'ERROR ROTAXS: PROPER, IORDER',proper,iorder
         call quit('ROTAXS ERROR: iorder < 2')
      END IF
C ::: If we are looking for improper rotational axes, we have :::
C ::: to multiply the order by a factor of two.               :::
      IF (PROPER) THEN
         IORDR = IORDER
      ELSE
         IORDR=2*IORDER
      END IF
      CON2PO = 2.0D0*PI/DBLE(IORDR)
      TOLDEG = MAX(TOLLRN,1.D-4)

C ::: The original matrix is copied because we modify ATMARR below in ROTMOL. :::
      N_ATOMS = NINT(ORGARR(1,0))
      ATMARR(1,0) = N_ATOMS
      DO 100 I=1,N_ATOMS
         ATMARR(1,I) = ORGARR(1,I)
         ATMARR(2,I) = ORGARR(2,I)
         ATMARR(3,I) = ORGARR(3,I)
         ATMARR(4,I) = ORGARR(4,I)
         ATMARR(6,I) = ORGARR(6,I)
 100  CONTINUE
C
      TEMP=SQRT(ROTVEC(1)*ROTVEC(1)+ROTVEC(2)*ROTVEC(2)+
     &          ROTVEC(3)*ROTVEC(3))
      ROTVEC(1)=ROTVEC(1)/TEMP
      ROTVEC(2)=ROTVEC(2)/TEMP
      ROTVEC(3)=ROTVEC(3)/TEMP
C ::: The molcule is rotated so that the possible rotational axis :::
C ::: lies along the Z-axis                                       :::
      V1(1) = 0.0D0
      V1(2) = 0.0D0
      V1(3) = 1.0D0
      VROT(1) =  V1(2)*ROTVEC(3) - ROTVEC(2)*V1(3)
      VROT(2) = -V1(1)*ROTVEC(3) + ROTVEC(1)*V1(3)
      VROT(3) =  V1(1)*ROTVEC(2) - ROTVEC(1)*V1(2)
      IF (VROT(1)*VROT(1)+VROT(2)*VROT(2)+VROT(3)*VROT(3) .GT.
     &    TOLLRN_2) THEN
         DEG = VECANG(ROTVEC, V1)
         CALL ROTMOL(MXM, ATMARR, VROT, -DEG)
      END IF

      ROTAXS = .TRUE.

C ::: Initialize IUSED array to all atoms not used :::
      IUSED(1:N_ATOMS) = 0
      IF (PROPER) THEN
C ::: All atoms lying on the Z-axis are dropped for proper rotations:::
         DO I = 1,N_ATOMS
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            IF (V1(1)*V1(1)+V1(2)*V1(2).LT.TOLLRN_2) IUSED(I)=1
         END DO
      ELSE
C ::: An atom lying in (0,0,0) is dropped for improper rotations:::
         DO I = 1,N_ATOMS
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            IF (V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3).LT.TOLLRN_2)
     &         IUSED(I)=1
         END DO
      END IF

C ::: Then we go through the rest :::
      DO I = 1, N_ATOMS
         IF (IUSED(I) .EQ. 1) CYCLE

C ::: We take the XY-component. :::
         V1(1) = ATMARR(1,I)
         V1(2) = ATMARR(2,I)
         V1(3) = 0.0D0
         V1NRM2 = V1(1)*V1(1)+V1(2)*V1(2)
         IORD = 0 ! IORD counts how many atoms related by the proper/improper rotation to atom I
         TOLLDS = 2.0D0*TOLLRN*SQRT(V1NRM2) + TOLLRN_2
         DO J = I + 1, N_ATOMS
            IF (IUSED(J) .EQ. 1) CYCLE

C ::: Check that atoms I and J are equivalent; i.e. same charge and same isotope number :::
            EQL = (NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J)) .AND.
     &             NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J)))
            IF (.NOT. EQL) CYCLE

            V2(1) = ATMARR(1,J)
            V2(2) = ATMARR(2,J)
            V2(3) = 0.0D0

C ::: Check that the two atoms have same distance to z-axis :::
            V2NRM2 = V2(1)*V2(1)+V2(2)*V2(2)
            IF ( ABS(V1NRM2-V2NRM2) .GT. TOLLDS) CYCLE   ! | r^2 - (r+d)^2 | = | 2*d*r + d^2 |, thus TOLLDS

C ::: Check that the angle between them is correct based on the rotational order :::
            DEG = VECANG(V1, V2)
            IJ_ORD = NINT(DEG/CON2PO)
            DEG = DEG - IJ_ORD*CON2PO
            IF ( ABS(DEG) .GT. TOLDEG ) CYCLE

C ::: Finally check if the two atoms have matching Z-coordinates  :::
C ::: (This test is last because we need to know IJ_ORD)          :::
            IJ_ORD = MOD(IJ_ORD,2)
            IF (PROPER .OR. IJ_ORD .EQ. 0) THEN ! all even order (S_2n)^(2a) = (C_n)^a, i.e. a proper rotation
               Z_TEST = ABS(ATMARR(3,J) - ATMARR(3,I)) ! Z_I = +Z_J ?
            ELSE
               Z_TEST = ABS(ATMARR(3,J) + ATMARR(3,I)) ! Z_I = -Z_J ?
            END IF
            IF (Z_TEST .GT. TOLLRN) CYCLE

C ::: If we get here, then atom J is related to atom I by the specified symmetry element :::
            IUSED(J) = 1
            IORD = IORD + 1

            IF (IORD .EQ. IORDR-1) EXIT ! OK, finished for atom I, found IORDR-1 related atoms
         END DO  !  J = I + 1, N_ATOMS

C ::: We also have to check that the symmetry is complete :::
         IF (IORD .EQ. IORDR-1) THEN
            IUSED(I) = 1 ! not really needed, but in case someone wants to add debug print of IUSED
         ELSE
            ROTAXS = .FALSE.
            EXIT
         END IF
      END DO ! I = 1, N_ATOMS
      RETURN
      END ! logical function ROTAXS
C  /* Deck findax */
      SUBROUTINE FINDAX(MXM, ATMARR, DRTAXS, TMPARR, PROPER)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                       ::
C ::       Finds all rotational axes in the molecule.      ::
C ::       PROPER determines if we are looking for         ::
C ::       proper or improper rotations.                   ::
C ::                                                       ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxsymm.h"
#include "mxcent.h"
      DIMENSION ATMARR(6,0:MXM), DRTAXS(6,0:MAXAXS)
      DIMENSION TMPARR(7,0:MXM)
      DIMENSION V1(3), V2(3), V3(3), VSTORE(3)
      DIMENSION ICNT(2), IUSED(MXCENT)
      LOGICAL   ROTAXS, PROPER

      DRTAXS(1,0) = 0.0D0
      N_ATOMS = NINT(ATMARR(1,0))
      IF (PROPER) THEN
         MINORD = 2
      ELSE
!        MINORD = 1
         MINORD = 2
         ! NB! minimum n tested for S_n is the double of this order, i.e. S_4.
         ! hjaaj Feb 2016: ! for improper rotations order 2 is inversion
         ! for _any_ rotation axes, so MINORD=2 would add a very
         ! big number of different S_2 axes, which would not be deleted
         ! by DLDPAX (and are redundant because inversion S_2 has
         ! already been detected elsewhere).
         ! Old code had MINORD = 1 corresponding to minimum S_2, but:
         ! IORD=1 (S_1) are reflection planes, they are found in FIND_MIR.
         ! IORD=2 (S_2) is inversion, it is found in INVMOL
      END IF
      VSTORE(1:3) = 0.0D0
      IUSED(1:N_ATOMS) = 0
#if HERSYM_DEBUG > 10
      write(lupri,*) 'FINDAX entry: PROPER,MINORD',PROPER,MINORD
#endif
      I = 1
 100  CONTINUE
      IF (I .LE. N_ATOMS) THEN
C ::: "I" will be the lowest index in a set of equivalent atoms. :::
C ::: We need to find the rest
         J = I+1
         IORD = 1
         V1(1) = ATMARR(1,I)
         V1(2) = ATMARR(2,I)
         V1(3) = ATMARR(3,I)
         V1NRM2 = V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
         TOLLDS = 2.0D0*TOLLRN*SQRT(V1NRM2) + TOLLRN_2   ! | r^2 - (r+d)^2 | = | 2*d*r + d^2 |, thus TOLLDS
C ::: The number of equivalent atoms is counted in IORD. :::
 110  CONTINUE
      IF (J .LE. N_ATOMS) THEN
         V2(1) = ATMARR(1,J)
         V2(2) = ATMARR(2,J)
         V2(3) = ATMARR(3,J)
         V2NRM2 = V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
         IF ((NINT(ATMARR(4,I)) .EQ. NINT(ATMARR(4,J))) .AND. ! same charge
     &       (NINT(ATMARR(6,I)) .EQ. NINT(ATMARR(6,J))) .AND. ! same isotope
     &       (ABS( V1NRM2 - V2NRM2 ) .LT. TOLLDS )) THEN      ! same distance to origo
            IORD = IORD + 1
            J = J + 1
            GOTO 110
         END IF
      END IF
#if HERSYM_DEBUG > 10
         write(lupri,*) 'FINDAX: atom I, IORD',I,IORD
         write(lupri,*) 'FINDAX: K range',J-IORD,J-1
#endif
C ::: Number of equivalent atoms must be two or higher to be of interest. :::
         IF (IORD .GT. 1) THEN
C ::: We "add" all equivalent atoms to yield a vector. If this vector is :::
C ::: different from the zero vector, it is the only possible rotational :::
C ::: axis. If it adds up to the zero vector, the atoms defines a        :::
C ::: linear, planar or highly symmetrical structure (it has to be       :::
C ::: symmetrical around origo).                                         :::
            V1(1) = 0.0D0
            V1(2) = 0.0D0
            V1(3) = 0.0D0
            DO 120 K=I,I+IORD-1
               V1(1) = V1(1) + ATMARR(1,K)
               V1(2) = V1(2) + ATMARR(2,K)
               V1(3) = V1(3) + ATMARR(3,K)
 120        CONTINUE
            V1NRM2 = V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
            ! write(lupri,*) 'FINDAX: V1, V1NRM2',V1,V1NRM2
            IF (V1NRM2 .LT. TOLLRN_2) THEN
C ::: All but the linear case has order higher than two. :::
               IF (IORD .GT. 2) THEN
                  DO 133 J=I,I+IORD-1
                     V1(1) = ATMARR(1,J)
                     V1(2) = ATMARR(2,J)
                     V1(3) = ATMARR(3,J)
C ::: We let each possible pair of atoms define a plane, and examine :::
C ::: the normal vector of this plane as a possible rotational axis. :::
                  DO 130 K=J,I+IORD-1 ! only check each pair once
                     V2(1)= ATMARR(1,K)
                     V2(2)= ATMARR(2,K)
                     V2(3)= ATMARR(3,K)
                     V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
                     V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
                     V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
                     IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                  .GT. TOLLRN_2) THEN
                        DO L=IORD,MINORD,-1
                           IF (ROTAXS(MXM,ATMARR,TMPARR,L,V3,PROPER))
     &                     THEN
                              ! write(lupri,*) 'FINDAX-a found axis',V3,L
                              CALL ADD2AR(MAXAXS, DRTAXS, L, V3)
                              CALL DLDPAX(MAXAXS, DRTAXS)
                           END IF
                        END DO
                     END IF
 130              CONTINUE
C ::: We also examine all axes defined by ONE atom. :::
!     hjaaj: at most of order IORD-1, as one atom is on the axis
!            meaning we have max IORD-1 equivalent atoms not on the axis.
                  IF ( V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3)
     &               .GT. TOLLRN_2) THEN
                     DO L = IORD - 1, MINORD, -1
                     IF (ROTAXS(MXM,ATMARR,TMPARR,L,V1,PROPER)) THEN
                        ! write(lupri,*) 'FINDAX-b found axis',V1,L
                        CALL ADD2AR(MAXAXS, DRTAXS, L, V1)
                        CALL DLDPAX(MAXAXS, DRTAXS)
                     END IF
                     END DO
                  END IF
 133              CONTINUE
C ::: If the order is even, and the symmetry high, we might have missed    :::
C ::: axes lying between the atoms. We therefore test all possibilities.   :::
C ::: These tests are so thorough, they only have to be performed once.    :::
                  IF ((MOD(IORD,2).EQ.0).AND.(IUSED(I).EQ.0)) THEN
                     V1(1) = ATMARR(1,I)
                     V1(2) = ATMARR(2,I)
                     V1(3) = ATMARR(3,I)
C ::: We are going to add from 2 to 3 vectors. The matrix ICnt keeps :::
C ::: track of indexes, so that all combinations are tried.          :::
                     DO 160 K = 1, MIN(IORD/2-1,2)
                        ICNT(1) = 1
                        ICNT(2) = 1
C                       DO WHILE (ICNT(K) .LE. IORD)
 166                    CONTINUE
                        IF (ICNT(K) .LE. IORD) THEN
                           V2(1) = V1(1)
                           V2(2) = V1(2)
                           V2(3) = V1(3)
                           DO L = 1, K
                              M = I-1+ICNT(L)
                              V2(1) = V2(1) + ATMARR(1,M)
                              V2(2) = V2(2) + ATMARR(2,M)
                              V2(3) = V2(3) + ATMARR(3,M)
                           END DO ! L = 1,K
C ::: We check for all possible orders of rotation from Ord/2 and down :::
                           IF((V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3))
     &                        .GT. TOLLRN_2) THEN
                              DO M = IORD/2, MINORD, -1
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V2,PROPER))
     &                        THEN
                              ! write(lupri,*) 'FINDAX-c found axis',V2,M
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V2)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                              END DO
                           END IF
                           ICNT(1) = ICNT(1) + 1
                           IF ((K.GT.1).AND.(ICNT(1).GT.IORD)) THEN
                              ICNT(1) = 1
                              ICNT(2) = ICNT(2) + 1
                           END IF
                           GOTO 166
                        END IF
C                       ... END DO WHILE
 160                 CONTINUE
C ::: Molecules with high symmetry and an inversion center, might have  :::
C ::: a central axis that is difficult to get at. To find it, we add    :::
C ::: all atoms on one side of the origo. All these have positive dot   :::
C ::: products between themselves.                                      :::
                     V1(1) = 0.0D0
                     V1(2) = 0.0D0
                     V1(3) = 0.0D0
                     DO 169 ICRDAX = 1, 3
                        V1(ICRDAX) = 1.0D0
                        V3(1) = 0.0D0
                        V3(2) = 0.0D0
                        V3(3) = 0.0D0
                        DO L = I, I+IORD-1
                           V2(1) = ATMARR(1,L)
                           V2(2) = ATMARR(2,L)
                           V2(3) = ATMARR(3,L)
                           IF (V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3)
     &                        .GT. 0.0D0) THEN
                              V3(1) = V3(1) + V2(1)
                              V3(2) = V3(2) + V2(2)
                              V3(3) = V3(3) + V2(3)
                           END IF
                        END DO
                        IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                      .GT. TOLLRN_2) THEN
                           DO M = IORD/2, MINORD, -1
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V3,PROPER))
     &                           THEN
                              ! write(lupri,*) 'FINDAX-d found axis',V3,M
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V3)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                           END DO
                        END IF
                        V1(ICRDAX) = 0.0D0
 169                 CONTINUE
C ::: To be sure to get all the rotational axes in high order symmetry :::
C ::: species, all sums and differences of rotational axes are tested  :::
C ::: as well. Maximum order is the highest order among the axes.      :::
                     NAX = NINT(DRTAXS(1,0))
                     MXORDR = NINT(DRTAXS(4,1))
                     IF (NAX .GE. 2) THEN
                        K = 1
                        ICHANG = 0
 175                    CONTINUE
                        IF (K .LT. NAX) THEN
                           L = K + 1
                           V1(1) = DRTAXS(1,K)
                           V1(2) = DRTAXS(2,K)
                           V1(3) = DRTAXS(3,K)
 177                       CONTINUE
                           IF ((L .LE. NAX) .AND. (ICHANG .EQ. 0))THEN
                             V2(1) = V1(1) + DRTAXS(1,L)
                             V2(2) = V1(2) + DRTAXS(2,L)
                             V2(3) = V1(3) + DRTAXS(3,L)
                             V3(1) = V1(1) - DRTAXS(1,L)
                             V3(2) = V1(2) - DRTAXS(2,L)
                             V3(3) = V1(3) - DRTAXS(3,L)
                             IF (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3)
     &                                     .LT. ZERTOL) V2(1) = 1.0D0
                             IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                                     .LT. ZERTOL) V3(1) = 1.0D0
                             DO M = MXORDR, MINORD, -1
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V2,PROPER))
     &                           THEN
                              ! write(lupri,*) 'FINDAX-e found axis',V2,M
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V2)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                              IF (ROTAXS(MXM,ATMARR,TMPARR,M,V3,PROPER))
     &                           THEN
                              ! write(lupri,*) 'FINDAX-f found axis',V3,M
                                 CALL ADD2AR(MAXAXS, DRTAXS, M, V3)
                                 CALL DLDPAX(MAXAXS, DRTAXS)
                              END IF
                             END DO
                             IF(NAX .LT. NINT(DRTAXS(1,0))) ICHANG = 1
                             L = L + 1
                             GOTO 177
                           END IF
                           K = K + 1
                           IF (ICHANG .NE. 0) THEN
                              K = 1
                              NAX = NINT(DRTAXS(1,0))
                              MXORDR = NINT(DRTAXS(4,1))
                              ICHANG = 0
                           END IF
                           GOTO 175
                        END IF
                     END IF
                  END IF
C ::: Then we move on to the linear case :::
               ELSE
                  ! write(lupri,*) 'FINDAX: linear case'
                  V1(1) = ATMARR(1,I)
                  V1(2) = ATMARR(2,I)
                  V1(3) = ATMARR(3,I)
C ::: We check for an axis of infinite (more precisely 99.) order. :::
                  IF (ROTAXS(MXM, ATMARR, TMPARR, 99, V1, PROPER)) THEN
                     CALL ADD2AR(MAXAXS, DRTAXS, 99, V1)
                     CALL DLDPAX(MAXAXS, DRTAXS)
C ::: Otherwise we look for all other possible orders along the line. :::
                  ELSE
                     DO L = 2, 8
                        IF (ROTAXS(MXM,ATMARR,TMPARR,L,V1,PROPER))
     &                     THEN
                           CALL ADD2AR(MAXAXS, DRTAXS, L, V1)
                           CALL DLDPAX(MAXAXS, DRTAXS)
                        END IF
                     END DO
C ::: Another possibility is a C2 axis perpendicular to the line. :::
C ::: We need another vector to specify the axis, so we store it. :::
C ::: If a vector already is stored, we check if the two span a   :::
C ::: plane. If they do, a C2 axis perp. to the plane is tested.  :::
                     IF (VSTORE(1)*VSTORE(1)+VSTORE(2)*VSTORE(2)+
     &                   VSTORE(3)*VSTORE(3) .LT. TOLLRN_2) THEN
                        VSTORE(1) = V1(1)
                        VSTORE(2) = V1(2)
                        VSTORE(3) = V1(3)
                     ELSE
                        V3(1) =  V1(2)*VSTORE(3) - VSTORE(2)*V1(3)
                        V3(2) = -V1(1)*VSTORE(3) + VSTORE(1)*V1(3)
                        V3(3) =  V1(1)*VSTORE(2) - VSTORE(1)*V1(2)
                        IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                      .GT. TOLLRN_2) THEN
                           IF (ROTAXS(MXM,ATMARR,TMPARR,2,V3,PROPER))
     &                        THEN
                              CALL ADD2AR(MAXAXS, DRTAXS, 2, V3)
                              CALL DLDPAX(MAXAXS, DRTAXS)
                           END IF
                        END IF
                     END IF
                  END IF
               END IF
C ::: The easiest case is when the atoms define the only possible axis. :::
            ELSE
               IF (ROTAXS(MXM, ATMARR, TMPARR, IORD, V1, PROPER))
     &            THEN
                  ! write(lupri,*) 'FINDAX-g found axis',V1,IORD
                  CALL ADD2AR(MAXAXS, DRTAXS, IORD, V1)
                  CALL DLDPAX(MAXAXS, DRTAXS)
               END IF
            END IF
         END IF ! IORD .GT. 1
C ::: Finally the atoms are marked as used. :::
         DO M = I, I+IORD-1
            IUSED(M) = 1
         END DO
CRF   Set I to first atom not in this equivalent set and repeat
         I = I + IORD
         GOTO 100
      END IF
      CALL DLDPAX(MAXAXS, DRTAXS)
      CALL SRTATM(MAXAXS, DRTAXS, .FALSE.) ! sort so highest order first
      RETURN
      END ! subroutine FINDAX
C  /* Deck fndima */
      SUBROUTINE FNDIMA(MXM, ATMARR, AXIMAR, AXARR, TMPARR)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                                ::
C ::       Finds all improper rotational axes in the molecule       ::
C ::                                                                ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), AXIMAR(6,0:MAXAXS)
         DIMENSION AXARR(6,0:MAXAXS), TMPARR(7,0:MXM)
         DIMENSION V1(3), V2(3)
         N_ATOMS = NINT(ATMARR(1,0))
C ::: And we call upon FindAx to find all improper axes for this molecule. :::
         CALL FINDAX(MXM, ATMARR, AXIMAR, TMPARR, .FALSE.)
C ::: We have to multiply all orders by a factor of two to get S_n not S_2n :::
         DO 100 I = 1, NINT(AXIMAR(1,0))
            AXIMAR(4,I) = 2.0D0*AXIMAR(4,I)
 100     CONTINUE
C ::: Finally we remove all improper rotational axes with the same order :::
C ::: and parallellity to proper axes.                                   :::
         DO 200 I = 1, NINT(AXARR(1,0))
            V1(1) = AXARR(1,I)
            V1(2) = AXARR(2,I)
            V1(3) = AXARR(3,I)
            ORD = NINT(AXARR(4,I))
            J = 1
 210        CONTINUE
            IF (J .LE. NINT(AXIMAR(1,0))) THEN
               V2(1) =  V1(2)*AXIMAR(3,J) - AXIMAR(2,J)*V1(3)
               V2(2) = -V1(1)*AXIMAR(3,J) + AXIMAR(1,J)*V1(3)
               V2(3) =  V1(1)*AXIMAR(2,J) - AXIMAR(1,J)*V1(2)
               IF ((ORD .GE. NINT(AXIMAR(4,J))) .AND. (V2(1)*V2(1)+
     &               V2(2)*V2(2)+V2(3)*V2(3) .LT. TOLLRN_2)) THEN
                  DO 205 K = 1, NINT(AXIMAR(1,0)) - J
                     AXIMAR(1,J+K-1) = AXIMAR(1,J+K)
                     AXIMAR(2,J+K-1) = AXIMAR(2,J+K)
                     AXIMAR(3,J+K-1) = AXIMAR(3,J+K)
                     AXIMAR(4,J+K-1) = AXIMAR(4,J+K)
 205              CONTINUE
                  AXIMAR(1,0) = AXIMAR(1,0) - 1.0D0
               ELSE
                  J = J + 1
               END IF
            GOTO 210
            END IF
 200     CONTINUE
         RETURN
      END
C  /* Deck FIND_MIR */
      SUBROUTINE FIND_MIR(MXM, ATMARR, DMRPLN, AXARR, TMPARR)
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                     ::
C ::       Finds all mirror planes in the molecule       ::
C ::                                                     ::
C :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
         DIMENSION ATMARR(6,0:MXM), DMRPLN(6,0:MAXMIR)
         DIMENSION AXARR(6,0:MAXAXS), TMPARR(7,0:MXM)
         DIMENSION V1(3), V2(3), V3(3), VSTORE(3)
         LOGICAL   MIRMOL, PLNMOL, LINMOL
         DMRPLN(1,0) = 0.0D0
         N_ATOMS = NINT(ATMARR(1,0))
         VSTORE(1:3) = 0.0D0
C ::: First we check all planes perpendicular to one atom pos. :::
         DO 100 I = 1, N_ATOMS
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            IF ((V1(1)*V1(1)+V1(2)*V1(2)+V1(3)*V1(3).GT.TOLLRN_2)
     &         .AND. (MIRMOL(MXM, ATMARR, TMPARR, V1, N2))) THEN
               CALL ADD2AR(MAXMIR, DMRPLN, 0, V1)
               CALL DLDPAX(MAXMIR, DMRPLN)
            END IF
 100     CONTINUE
C ::: Next thing, we go through all pairs of equal atoms, examining      :::
C ::: possible mirror planes through and between them. For the linear    :::
C ::: cases, a vector is stored. If one already is stored, we got enough :::
C ::: information to test a mirror plane.                                :::
         DO 200 I = 1, N_ATOMS-1
            V1(1) = ATMARR(1,I)
            V1(2) = ATMARR(2,I)
            V1(3) = ATMARR(3,I)
            J = I + 1
            IATMNR = NINT(ATMARR(4,I))
            IATMMS = NINT(ATMARR(6,I))
 210        CONTINUE
            IF ( NINT(ATMARR(4,J)).EQ.IATMNR .AND.
     &           NINT(ATMARR(6,J)).EQ.IATMMS ) THEN
               V2(1) = ATMARR(1,J)
               V2(2) = ATMARR(2,J)
               V2(3) = ATMARR(3,J)
C ::: Checking the cross product, tells us whether the atoms lie on a line. :::
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)).GT. TOLLRN_2)
     &         THEN
C ::: The cross product defines the plane both atoms lie in. :::
                  IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                     CALL ADD2AR(MAXMIR, DMRPLN, 0, V3)
                     CALL DLDPAX(MAXMIR, DMRPLN)
                  END IF
C ::: The difference between the two vectors, define the plane between :::
C ::: the atoms.                                                       :::
                  V3(1) = V1(1) - V2(1)
                  V3(2) = V1(2) - V2(2)
                  V3(3) = V1(3) - V2(3)
                  IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                     CALL ADD2AR(MAXMIR, DMRPLN, 0, V3)
                     CALL DLDPAX(MAXMIR, DMRPLN)
                  END IF
               ELSE
C ::: The linear case calls for another vector to specify a plane. :::
                  IF ((VSTORE(1)*VSTORE(1)+VSTORE(2)*VSTORE(2)+
     &                 VSTORE(3)*VSTORE(3)) .LT. TOLLRN_2) THEN
                     VSTORE(1) = V1(1)
                     VSTORE(2) = V1(2)
                     VSTORE(3) = V1(3)
                  ELSE
                     V3(1) =  V1(2)*VSTORE(3) - VSTORE(2)*V1(3)
                     V3(2) = -V1(1)*VSTORE(3) + VSTORE(1)*V1(3)
                     V3(3) =  V1(1)*VSTORE(2) - VSTORE(1)*V1(2)
                     IF ((V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
     &                   .GT. TOLLRN_2) THEN
                        IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                           CALL ADD2AR(MAXMIR, DMRPLN, 0, V3)
                           CALL DLDPAX(MAXMIR, DMRPLN)
                        END IF
                     END IF
                  END IF
               END IF
               J = J + 1
               IF (J.LE.N_ATOMS) GOTO 210
            END IF
 200     CONTINUE
C ::: If the molecule is planar, not linear, and has no other symmetry :::
C ::: element than this plane, we have not found it. This must be      :::
C ::: taken care of.                                                   :::
         IF ((N_ATOMS .GT. 2) .AND. PLNMOL(MXM, ATMARR) .AND.
     &       (.NOT. LINMOL(MXM, ATMARR))) THEN
            V1(1) = ATMARR(1,1) - ATMARR(1,2)
            V2(1) = ATMARR(1,2) - ATMARR(1,3)
            V1(2) = ATMARR(2,1) - ATMARR(2,2)
            V2(2) = ATMARR(2,2) - ATMARR(2,3)
            V1(3) = ATMARR(3,1) - ATMARR(3,2)
            V2(3) = ATMARR(3,2) - ATMARR(3,3)
C ::: Take cross product of the two bonds. :::
            V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
            V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
            V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
            I = 4
 250        CONTINUE
            IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3) .LT. TOLLRN_2) THEN
               V2(1) = ATMARR(1,2) - ATMARR(1,I)
               V2(2) = ATMARR(2,2) - ATMARR(2,I)
               V2(3) = ATMARR(3,2) - ATMARR(3,I)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
               I = I + 1
               IF (I.LE.N_ATOMS) GOTO 250
            END IF
            IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3) .GT. TOLLRN_2) THEN
               IF (MIRMOL(MXM, ATMARR, TMPARR, V3, N2)) THEN
                  CALL ADD2AR(MAXMIR, DMRPLN, 0, V3)
                  CALL DLDPAX(MAXMIR, DMRPLN)
               END IF
            END IF
         END IF
C ::: After finding all planes, we would like to classify them, provided :::
C ::: we have found rotational axes. The axis with the highest order is  :::
C ::: used as reference, and the code is as follows:                     :::
C ::: Horizontal - 2   Vertical - 1   Other/Undecided - 0                :::
         IF (NINT(AXARR(1,0)) .GT. 0) THEN
            V1(1) = AXARR(1,1)
            V1(2) = AXARR(2,1)
            V1(3) = AXARR(3,1)
            DO 300 I = 1, NINT(DMRPLN(1,0))
               V2(1) = DMRPLN(1,I)
               V2(2) = DMRPLN(2,I)
               V2(3) = DMRPLN(3,I)
               V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
               V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
               V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
C ::: First the horizontal. :::
               IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3)
     &                                    .LT. TOLLRN_2) THEN
                  DMRPLN(4,I) = 2.0D0
C ::: Then the vertical. :::
               ELSE IF (ABS(V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))
     &               .LT. TOLLRN) THEN
                  DMRPLN(4,I) = 1.0D0
C ::: The rest must be diagonal/undecided. :::
               ELSE
                  DMRPLN(4,I) = 0.0D0
               END IF
 300        CONTINUE
            CALL SRTATM(MAXMIR, DMRPLN, .FALSE.)
         END IF
         RETURN
      END
C  /* Deck FIND_PGROUP */
      SUBROUTINE FIND_PGROUP(MXM, ATMARR, DRTAXS, DIRAXS, DMRPLN,
     &                       TMPARR, POINT_GROUP)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                              ::
C ::       Determines which point group molecule belongs to       ::
C ::                                                              ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "mxsymm.h"
#include "priunit.h"
      DIMENSION ATMARR(6,0:MXM), DRTAXS(6,0:MAXAXS)
      DIMENSION DIRAXS(6,0:MAXAXS), DMRPLN(6,0:MAXMIR)
      DIMENSION TMPARR(7,0:MXM)
      DIMENSION V1(3), V2(3)
      CHARACTER*15 POINT_GROUP
      LOGICAL   LINMOL, INVMOL
 111  FORMAT(I1)
 222  FORMAT(I2)
C
C ::: Find proper rotation axes
         CALL FINDAX(MXM, ATMARR, DRTAXS, TMPARR, .TRUE.)
C ::: If the number of rotational axes exceeds 10, the number :::
C ::: of improper rotational axes is likely to be very high   :::
C ::: and we choose to skip their determination.              :::
      IF (NINT(DRTAXS(1,0)) .LE. 10) THEN
         CALL FNDIMA(MXM, ATMARR, DIRAXS, DRTAXS, TMPARR)
      ELSE
         DIRAXS(1,0) = -1.0D0
      END IF
      CALL FIND_MIR(MXM, ATMARR, DMRPLN, DRTAXS, TMPARR)
      IF (LINMOL(MXM, ATMARR)) THEN
         IF (INVMOL(MXM, ATMARR)) THEN
            POINT_GROUP = 'D(oo,h)'
         ELSE
            POINT_GROUP = 'C(oo,v)'
         END IF
      ELSE ! not linear
C ::: We count rotational axes with order > 2 :::
         ITEMP = 0
         DO 100 I = 1, NINT(DRTAXS(1,0))
            IF (NINT(DRTAXS(4,I)) .GT. 2) ITEMP = ITEMP + 1
 100     CONTINUE
         IF (ITEMP .GE. 2) THEN
            IF (.NOT. INVMOL(MXM, ATMARR)) THEN
               POINT_GROUP = 'T(d)'
            ELSE IF (NINT(DRTAXS(4,1)) .EQ. 5) THEN
               POINT_GROUP = 'I(h)'
            ELSE
               POINT_GROUP = 'O(h)'
            END IF
         ELSE IF (NINT(DRTAXS(1,0)) .EQ. 0) THEN
            ! no rotation axis
            IF (NINT(DMRPLN(1,0)) .GT. 0) THEN
               POINT_GROUP = 'C(s)'
            ELSE IF (INVMOL(MXM, ATMARR)) THEN
               POINT_GROUP = 'C(i)'
            ELSE
               POINT_GROUP = 'C(1)'
            END IF
         ELSE ! we have at least one rotation axis
            V1(1) = DRTAXS(1,1)
            V1(2) = DRTAXS(2,1)
            V1(3) = DRTAXS(3,1)
C ::: Number of perpendicular C2 axes is counted. :::
            ITEMP = 0
            DO 200 I = 2, NINT(DRTAXS(1,0))
               V2(1) = DRTAXS(1,I)
               V2(2) = DRTAXS(2,I)
               V2(3) = DRTAXS(3,I)
               IF ((NINT(DRTAXS(4,I)).EQ.2) .AND. (ABS(V1(1)*V2(1)+
     &           V1(2)*V2(2)+V1(3)*V2(3)).LT.TOLLRN)) ITEMP = ITEMP+1
 200        CONTINUE
            IF (ITEMP .GE. NINT(DRTAXS(4,1))) THEN
               ! perpendicular C2 axes found
C              ... any horizontal mirror plane is first in DMRPLN
               IF ((NINT(DMRPLN(1,0)) .GT. 0) .AND.
     &             (NINT(DMRPLN(4,1)) .EQ. 2)) THEN
                  IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                     POINT_GROUP = 'D(nh)'
                     WRITE(POINT_GROUP(3:3),111) NINT(DRTAXS(4,1))
                  ELSE
                     POINT_GROUP = 'D(nnh)'
                     WRITE(POINT_GROUP(3:4),222) NINT(DRTAXS(4,1))
                  END IF
               ELSE
C ::: Number of perpendicular mirror planes is counted. :::
                  ITEMP = 0
                  DO 300 I = 1, NINT(DMRPLN(1,0))
                     IF (NINT(DMRPLN(4,I)) .EQ. 1) ITEMP = ITEMP + 1
 300              CONTINUE
                  IF (ITEMP .GE. NINT(DRTAXS(4,1))) THEN
                     IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                        POINT_GROUP = 'D(nd)'
                        WRITE(POINT_GROUP(3:3),111) NINT(DRTAXS(4,1))
                     ELSE
                        POINT_GROUP = 'D(nnd)'
                        WRITE(POINT_GROUP(3:4),222) NINT(DRTAXS(4,1))
                     END IF
                  ELSE
                     IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                        POINT_GROUP = 'D(n)'
                        WRITE(POINT_GROUP(3:3),111) NINT(DRTAXS(4,1))
                     ELSE
                        POINT_GROUP = 'D(nn)'
                        WRITE(POINT_GROUP(3:4),222) NINT(DRTAXS(4,1))
                     END IF
                  END IF
               END IF
            ELSE IF (NINT(DMRPLN(1,0)) .GT. 0) THEN
             ! no perpendicular C2 axes, but we have mirror planes
             IF (NINT(DMRPLN(4,1)) .EQ. 2) THEN
               ! we have a horizontal mirror plane
               IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                  POINT_GROUP = 'C(nh)'
                  WRITE(POINT_GROUP(3:3),111) NINT(DRTAXS(4,1))
               ELSE
                  POINT_GROUP = 'C(nnh)'
                  WRITE(POINT_GROUP(3:4),222) NINT(DRTAXS(4,1))
               END IF
             ELSE ! only vertical mirror planes
C ::: Number of perpendicular mirror planes is counted. :::
               ITEMP = 0
               DO 400 I = 1, NINT(DMRPLN(1,0))
                  IF (NINT(DMRPLN(4,I)) .EQ. 1) ITEMP = ITEMP + 1
 400           CONTINUE
               IF (ITEMP .GE. NINT(DRTAXS(1,0))) THEN
                  IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                     POINT_GROUP = 'C(nv)'
                     WRITE(POINT_GROUP(3:3),111) NINT(DRTAXS(4,1))
                  ELSE
                     POINT_GROUP = 'C(nnv)'
                     WRITE(POINT_GROUP(3:4),222) NINT(DRTAXS(4,1))
                  END IF
               ELSE
                  WRITE(LUPRI,'(//A)')
     &               'FIND_PGROUP ERROR, cannot find point group'
                  CALL QUIT(
     &               'FIND_PGROUP ERROR, cannot find point group')
               END IF
             END IF
            ELSE ! no mirror planes
               IF (NINT(DIRAXS(1,0)) .LT. 0)
     &            CALL FNDIMA(MXM, ATMARR, DIRAXS, DRTAXS, TMPARR)
               IF ((NINT(DIRAXS(1,0)) .GE. 1) .AND.
     &             (NINT(DIRAXS(4,1)).GE.2*NINT(DRTAXS(4,1))))THEN
                  ITEMP = 2*NINT(DRTAXS(4,1))
                  IF (ITEMP .LT. 10) THEN
                     POINT_GROUP = 'S(n)'
                     WRITE(POINT_GROUP(3:3),111) ITEMP
                  ELSE
                     POINT_GROUP = 'S(nn)'
                     WRITE(POINT_GROUP(3:4),222) ITEMP
                  END IF
               ELSE
                  IF (NINT(DRTAXS(4,1)) .LT. 10) THEN
                     POINT_GROUP = 'C(n)'
                     WRITE(POINT_GROUP(3:3),111) NINT(DRTAXS(4,1))
                  ELSE
                     POINT_GROUP = 'C(nn)'
                     WRITE(POINT_GROUP(3:4),222) NINT(DRTAXS(4,1))
                  END IF
               END IF
            END IF ! if perpendicular C2 - else if mirror planes
         END IF
      END IF
      RETURN
      END
C  /* Deck tstpar */
      SUBROUTINE TSTPAR(MXM,ARR,TARR,IATNO1,IATNO2,ICODE,VEC,GENELM)
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                                    ::
C ::       Test pair of atoms for given operation       ::
C ::                                                    ::
C ::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
         DIMENSION ARR(7, 0:MXM), TARR(6, 0:MXM)
         DIMENSION TMPARR(6, 0:2), VEC(3)
         LOGICAL   GENELM, MIRMOL, ROTAXS, INVMOL
         GENELM = .FALSE.
         TMPARR(1,0) = 2.0D0
         DO 100 I = 1,6
            TMPARR(I,1) = ARR(I, IATNO1)
            TMPARR(I,2) = ARR(I, IATNO2)
 100     CONTINUE
C ::: The icode determines what should be tested     :::
C ::: 1 - mirror plane   2 - C2 axis   3 - inversion :::
         IF( ((ICODE .EQ. 1) .AND. (MIRMOL(MXM,TMPARR,TARR,VEC,N2)))
     &   .or.((ICODE .EQ. 2) .AND. (ROTAXS(2,TMPARR,TARR,2,
     &                                         VEC,.TRUE.)))
     &   .or.((ICODE .EQ. 3) .AND. (INVMOL(2,TMPARR)))) THEN
C             GENELM = ((NINT(ARR(6,IATNO1)) .EQ. 0) .OR.
C     &                 (NINT(ARR(6,IATNO2)) .EQ. 0))
            TEMP = 0.0D0
            DO 110 I = 1, 3
               TEMP = TEMP + TMPARR(I,1) - TMPARR(I,2)
 110        CONTINUE
            IF (TEMP .GT. 0.0D0) THEN
               GENELM = (NINT(ARR(7,IATNO2)) .LT. 1)
               ARR(7, IATNO2) = ANINT(ARR(7,IATNO2)+1.0D0)
            ELSE
               GENELM = (NINT(ARR(7,IATNO1)) .LT. 1)
               ARR(7, IATNO1) = ANINT(ARR(7,IATNO1)+1.0D0)
            END IF
         END IF
         RETURN
      END
C  /* Deck sym_reduce */
      SUBROUTINE SYM_REDUCE(MXM, ATMARR, EXTARR, TMPARR, SYMELM,
     &                  internal_file)
C :::::::::::::::::::::::::::::::::::::::::::::::::
C ::                                             ::
C ::       Removes symmetry-dependent atoms      ::
C ::                                             ::
C :::::::::::::::::::::::::::::::::::::::::::::::::
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "mxsymm.h"
#include "molinp.h"
         DIMENSION ATMARR(6, 0:MXM), EXTARR(7, 0:MXM)
         DIMENSION TMPARR(6, 0:MXM), V1(3), V2(3), V3(3)
         CHARACTER STRING*9, SYMELM*9
         LOGICAL   MIRMOL, ROTAXS, INVMOL, NECESS, GENELM
         LOGICAL :: internal_file
         N_ATOMS = NINT(ATMARR(1,0))
         EXTARR(1,0) = N_ATOMS
         DO 100 I = 1, N_ATOMS
            DO 110 J = 1, 6
               EXTARR(J,I) = ATMARR(J,I)
 110        CONTINUE
            EXTARR(7,I) = 0.0D0
 100     CONTINUE
         IELMNT = 0
         SYMELM = '         '
         STRING = 'XYZYZXZXY'
         V1(1) = 0.0D0
         V1(2) = 0.0D0
         V1(3) = 0.0D0
         DO 190 I = 1, 3
            V1(I) = 1.0D0
            IF (MIRMOL(MXM, ATMARR, TMPARR, V1, N2)) THEN
               SYMELM(IELMNT*3+1:IELMNT*3+1) = STRING(I:I)
               IELMNT = IELMNT + 1
              J = 1
 130           CONTINUE
               IF (J .LT. N_ATOMS) THEN
                  IF ((ABS(EXTARR(I,J)) .GT. TOLLRN) .AND.
     &                 (NINT(EXTARR(7,J)) .EQ. 0)) THEN
                     K = J + 1
 150                 CONTINUE
                     IF ((K .LE. N_ATOMS) .AND.
     &                (NINT(EXTARR(4,K)) .EQ. NINT(EXTARR(4,J))) .AND.
     &                (NINT(EXTARR(6,K)) .EQ. NINT(EXTARR(6,J)))) THEN
                       CALL TSTPAR(MXM,EXTARR,TMPARR,J,K,1,V1,GENELM)
                        K = K + 1
                        GOTO 150
                     END IF
                  END IF
                  J = J + 1
                  GOTO 130
               END IF
            END IF
            V1(I) = 0.0D0
 190     CONTINUE
         NROT = 0
         DO 290 I = 1, 3
            V1(I) = 1.0D0
            IF( (IELMNT.LT.3) .AND.
     &          ROTAXS(MXM,ATMARR,TMPARR,2,V1,.TRUE.) ) THEN
               SYMELM(IELMNT*3+1:IELMNT*3+2) = STRING(2+2*I:3+2*I)
               IELMNT = IELMNT + 1
               NROT = NROT + 1
               J = 1
               NECESS = .FALSE.
 250           CONTINUE
               IF (J .LT. N_ATOMS) THEN
                  V2(1) = EXTARR(1,J)
                  V2(2) = EXTARR(2,J)
                  V2(3) = EXTARR(3,J)
                  V3(1) =  V1(2)*V2(3) - V2(2)*V1(3)
                  V3(2) = -V1(1)*V2(3) + V2(1)*V1(3)
                  V3(3) =  V1(1)*V2(2) - V2(1)*V1(2)
                  IF (V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3) .GT. TOLLRN_2)
     &               THEN
                     K = J + 1
 260                 CONTINUE
                     IF ((K .LE. N_ATOMS) .AND.
     &                (NINT(EXTARR(4,K)) .EQ. NINT(EXTARR(4,J))) .AND.
     &                (NINT(EXTARR(6,K)) .EQ. NINT(EXTARR(6,J)))) THEN
                        CALL TSTPAR(MXM,EXTARR,TMPARR,J,K,2,V1,GENELM)
                        IF (GENELM) NECESS = .TRUE.
                        K = K + 1
                        GOTO 260
                     END IF
                  END IF
                  J = J + 1
                  GOTO 250
               END IF
               IF (.NOT. NECESS) THEN
                  IELMNT = IELMNT - 1
                  SYMELM(IELMNT*3+1:IELMNT*3+2) = '  '
               END IF
C ::: Three rotations (D2) causes a problem, only two rotations are     :::
C ::: necessary. The last is removed here (extremely dirty solution!!!) :::
C ::: 12-1999 Torgeir: Expanded hack to make sure D2 symmetry is ok, and:::
C ::: make use of as much symmetry as possible.                         :::
               IF ((NROT .EQ. 3) .AND. (IELMNT .GT. 1)) THEN
                  CALL D2SYMM(ATMARR,N_ATOMS,MXM,SYMELM)
               END IF
C ::: Some molecules causes a problem, with two mirror planes (A & B) :::
C ::: and  a redundant rotation (AB), the last is removed             :::
C ::: (another dirty solution!!!)                                     :::
               IF ((IELMNT .EQ. 3) .AND. (NROT .EQ. 1)) THEN
                  IF ((SYMELM(1:1)//SYMELM(4:4)).EQ.SYMELM(7:8))THEN
                     IELMNT = IELMNT - 1
                     SYMELM(IELMNT*3+1:IELMNT*3+2) = '  '
                  END IF
               END IF
            END IF
            V1(I) = 0.0D0
 290     CONTINUE
         IF (IELMNT .EQ. 0) THEN
            IF (INVMOL(MXM, ATMARR)) THEN
               SYMELM(IELMNT*3+1:IELMNT*3+3) = 'XYZ'
               IELMNT = IELMNT + 1
               J = 1
 350           CONTINUE
               IF (J .LT. N_ATOMS) THEN
                  V2(1) = EXTARR(1,J)
                  V2(2) = EXTARR(2,J)
                  V2(3) = EXTARR(3,J)
                  IF ( (V2(1)*V2(1)+V2(2)*V2(2)+V2(3)*V2(3).GT.TOLLRN_2)
     &                .AND. (NINT(EXTARR(7,J)).EQ.0) ) THEN
                     K = J + 1
 360                 CONTINUE
                     IF ((K .LE. N_ATOMS) .AND.
     &                (NINT(EXTARR(4,K)) .EQ. NINT(EXTARR(4,J))) .AND.
     &                (NINT(EXTARR(6,K)) .EQ. NINT(EXTARR(6,J)))) THEN
                        CALL TSTPAR(MXM,EXTARR,TMPARR,J,K,3,V1,GENELM)
                        K = K + 1
                        GOTO 360
                     END IF
                  END IF
                  J = J + 1
                  GOTO 350
               END IF
            END IF
         END IF
C ::: All symmetry dependant atoms are removed. :::
         I = 1
         DO 600 J = 1, N_ATOMS
            IF (NINT(EXTARR(7,J)) .EQ. 0) THEN
               DO 610 K = 1, 6
                  ATMARR(K,I) = EXTARR(K,J)
 610           CONTINUE
               DO 615 K = 1, 3
                  IF (ABS(ATMARR(K,I)).LT.TOLLRN) ATMARR(K,I) = 0.0D0
 615           CONTINUE
               I = I + 1
            ELSE
               if (internal_file) then
               NC = NCLINE(NINT(EXTARR(5,J)))
               DO 620 K = NC, NMLINE - 1
                  MLINE(K) = MLINE(K+1)
 620           CONTINUE
               DO 640 K = 1, N_ATOMS
                  IF (NCLINE(NINT(EXTARR(5,K))) .GE. NC)
     &           NCLINE(NINT(EXTARR(5,K)))=NCLINE(NINT(EXTARR(5,K)))-1
 640           CONTINUE
               NMLINE = NMLINE - 1
               end if
            END IF
 600     CONTINUE
         ATMARR(1,0) = 1.0D0*(I - 1)
         RETURN
      END
C   /*Deck d2symm*/
      SUBROUTINE D2SYMM(ATMARR,NATOM,MXM,SYMELM)
C
C       ***************************************************************
C       *** This subroutine checks whether there is problems, if we ***
C       *** have D2-symmetry. These problems arise from atoms that  ***
C       *** are degenerate with respect to some rotations, and not  ***
C       *** degenerate with respect to others.                      ***
C       ***************************************************************
C
#include "implicit.h"
#include "priunit.h"
C
        PARAMETER (THRS = 1.0D-10)
        LOGICAL   ONAXS, HAVAXS
        CHARACTER*(*) SYMELM
        CHARACTER*2 ROTAXIS(3)
C
        DIMENSION ATMARR(6, 0:MXM), IATAXS(3), ONAXS(3), HAVAXS(3)
C
        ROTAXIS(1) = 'YZ'
        ROTAXIS(2) = 'XZ'
        ROTAXIS(3) = 'XY'
C
        DO I = 1, 3
           ONAXS (I) = .FALSE.
           HAVAXS(I) = .FALSE.
           IATAXS(I) = 0
        END DO
C
C       ********************************************************
C       *** Find out if there are pairs of atoms on the axis ***
C       ********************************************************
C
        NPRAXS = 0
        DO IATOM = 1, NATOM
           DIST = SQRT(ATMARR(1,IATOM)**2 + ATMARR(2,IATOM)**2
     &               + ATMARR(3,IATOM)**2)
           DO ICART = 1, 3
              IF ((ABS(DIST-ABS(ATMARR(ICART,IATOM))) .LT. THRS)
     &                                    .AND. (DIST .GE. THRS)) THEN
                 IATAXS(ICART) = IATAXS(ICART) + 1
                 IF (IATAXS(ICART) .EQ. 2) THEN
                    ONAXS(ICART) = .TRUE.
                    NPRAXS = NPRAXS + 1
                 END IF
              END IF
           END DO
        END DO
C
C       ***********************************************************
C       ***We figure out which of the three situations we have. ***
C       ***********************************************************
C
        IF (NPRAXS.EQ.0) THEN
C
C       ************************************************************
C       *** No pairs of atoms on any axis. We can treat it in D2 ***
C       ***             symmetry without modification.           ***
C       ************************************************************
C
        ELSE IF (NPRAXS.EQ.3) THEN
C
C       **********************************************************
C       *** Pairs of atoms on all three axis. We have to treat ***
C       *** it in C2 symmetry.                                 ***
C       **********************************************************
C
           SYMELM(2*3+1:2*3+2) = '  '
           SYMELM(  3+1:  3+2) = '  '
        ELSE
C
C       ****************************************************************
C       *** 1 axis with atom-pairs on it. We have to be careful with ***
C       *** which axes we use to generate the dependent atoms.       ***
C       ****************************************************************
C
C          *** Make sure the axes, that have atoms on it, are used. ***
C
           SYMELM(1:9) = '         '
           IEL = 0
           DO IAXS = 3, 1, -1
              IF (ONAXS(IAXS)) THEN
                 SYMELM(IEL*3+1:IEL*3+2) = ROTAXIS(IAXS)
                 HAVAXS(IAXS) = .TRUE.
                 IEL = IEL + 1
              END IF
           END DO
C
C          *** If necessary the second axis is arbitrarily chosen ***
C
           DO IAXS = 3, 1, -1
              IF ((IEL.LT.2) .AND. (.NOT.HAVAXS(IAXS))) THEN
                 SYMELM(IEL*3+1:IEL*3+2) = ROTAXIS(IAXS)
                 HAVAXS(IAXS) = .TRUE.
                 IEL = IEL + 1
              END IF
           END DO
        END IF
C
        RETURN
        END
! --- end of hersym.F ---
