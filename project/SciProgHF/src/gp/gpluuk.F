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
C  /* Deck inpck */
      INTEGER FUNCTION INPCK (NI,INTS)
C
C     Packs small positive integers into a compound integer
C     L. Visscher
C
C     August 1997
C
#include "implicit.h"
#if defined (INT_STAR8)
      PARAMETER (NBITS=64)
#else
      PARAMETER (NBITS=32)
#endif
C
      DIMENSION INTS(NI)
C
      NIBIT = NBITS / NI
C
      IWORD = 0
C
      NFACT = 1
      MFACT = 2**NIBIT
C
      DO I = 1, NI
         IWORD = IWORD + INTS(I) * NFACT
         NFACT = NFACT * MFACT
      ENDDO
C
      INPCK = IWORD
C
      RETURN
C
      END
C  /* Deck iunpck */
      SUBROUTINE IUNPCK (INT,NI,INTS)
C
C     Unpacks small integers from a compound integer
C     L. Visscher
C
C     August 1997
C
#include "implicit.h"
#if defined (INT_STAR8)
      PARAMETER (NBITS=64)
#else
      PARAMETER (NBITS=32)
#endif
C
      DIMENSION INTS(NI)
C
      IWORD = INT
      NIBIT = NBITS / NI
C
      NFACT = 2**((NI-1)*NIBIT)
      MFACT = 2**NIBIT
C
      DO I = NI, 1, -1
         INTS(I) = IWORD / NFACT
         IWORD = MOD(IWORD,NFACT)
         NFACT = NFACT / MFACT
      ENDDO
C
      RETURN
C
      END
C  /* Deck indexx */
      SUBROUTINE INDEXX(N,ARR,INDX)
C***********************************************************************
C
C     This subroutine goes through array ARR and returns the ordering
C     index array INDX giving indices from smallest to largest
C
C     The array ARR is not touched.
C
C***********************************************************************
      INTEGER N,INDX(N),M,NSTACK
      REAL*8 ARR(N)
      PARAMETER (M=7,NSTACK=50)
      INTEGER I,INDXT,IR,ITEMP,J,JSTACK,K,L,ISTACK(NSTACK)
      REAL*8 A
      DO J=1,N
         INDX(J)=J
      END DO
      JSTACK=0
      L=1
      IR=N
 1    IF(IR-L.GE.M)THEN
         K=(L+IR)/2
         ITEMP=INDX(K)
         INDX(K)=INDX(L+1)
         INDX(L+1)=ITEMP
         IF(ARR(INDX(L+1)).GT.ARR(INDX(IR)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ARR(INDX(L)).GT.ARR(INDX(IR)))THEN
            ITEMP=INDX(L)
            INDX(L)=INDX(IR)
            INDX(IR)=ITEMP
         ENDIF
         IF(ARR(INDX(L+1)).GT.ARR(INDX(L)))THEN
            ITEMP=INDX(L+1)
            INDX(L+1)=INDX(L)
            INDX(L)=ITEMP
         ENDIF
         I=L+1
         J=IR
         INDXT=INDX(L)
         A=ARR(INDXT)
 3       CONTINUE
         I=I+1
         IF(ARR(INDX(I)).LT.A)GOTO 3
 4       CONTINUE
         J=J-1
         IF(ARR(INDX(J)).GT.A)GOTO 4
         IF(J.LT.I)GOTO 5
         ITEMP=INDX(I)
         INDX(I)=INDX(J)
         INDX(J)=ITEMP
         GOTO 3
 5       INDX(L)=INDX(J)
         INDX(J)=INDXT
         JSTACK=JSTACK+2
         IF(JSTACK.GT.NSTACK)
     $        CALL QUIT('NSTACK TOO SMALL IN INDEXX')
         IF(IR-I+1.GE.J-L)THEN
            ISTACK(JSTACK)=IR
            ISTACK(JSTACK-1)=I
            IR=J-1
         ELSE
            ISTACK(JSTACK)=J-1
            ISTACK(JSTACK-1)=L
            L=I
         ENDIF
      ELSE
         DO J=L+1,IR
            INDXT=INDX(J)
            A=ARR(INDXT)
            DO I=J-1,1,-1
               IF(ARR(INDX(I)).LE.A)GOTO 2
               INDX(I+1)=INDX(I)
            END DO
            I=0
 2          INDX(I+1)=INDXT
         END DO
         IF(JSTACK.EQ.0)RETURN
         IR=ISTACK(JSTACK)
         L=ISTACK(JSTACK-1)
         JSTACK=JSTACK-2
      ENDIF
      GOTO 1
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck IRECLEN*/
      FUNCTION IRECLEN (NREAL,NINTE,NLOGI)
C
C     Returns the record length for a record consisting of the specified
C     number of variables.
C
C     Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      IMPLICIT NONE
#include "priunit.h"
      INTEGER    NREAL,NINTE,NLOGI,IRECLEN, IRECL
      INTEGER*8  IRECL8
C
#if  defined (INT_STAR8)
      IRECL8 = 8 * NREAL + 8 * NINTE + 8 * NLOGI
#else
      IRECL8 = 8 * NREAL + 4 * NINTE + 4 * NLOGI
#endif

      IRECL  = IRECL8
      IF (NREAL .LT. 0 .OR. NINTE .LT. 0 .OR. 
     &    NLOGI .LT. 0 .OR. IRECL .LT. 0) THEN
         WRITE(LUPRI,'(//A,5(/A,I20))')
     &   'IRECLEN error: negative record length',
     &   '   NREAL',NREAL,
     &   '   NINTE',NINTE,
     &   '   NLOGI',NLOGI,
     &   '   record length (64 bit integers)',IRECL8,
     &   '   record length (32 bit integers)',IRECL
         CALL QUIT(
     &   'IRECLEN error: negative record length')
      END IF
      IF (IRECL .NE. IRECL8) THEN
         WRITE(LUPRI,'(//A,5(/A,I20))')
     &   'IRECLEN error: record length too big for 32 bit integers',
     &   '   NREAL',NREAL,
     &   '   NINTE',NINTE,
     &   '   NLOGI',NLOGI,
     &   '   record length (64 bit integers)',IRECL8,
     &   '   record length (32 bit integers)',IRECL
         CALL QUIT(
     &   'IRECLEN error: record length too big for 32 bit integers')
      END IF
C
      IRECLEN = IRECL
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLANKL (N,A)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Circular shift of character array such that the blanks are to the
C     left side of the array.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N
      CHARACTER*1 A(N)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IS = 0
      DO I = N, 1, -1
         IF (A(I).NE.' ') GOTO 2
         IS = IS + 1
      ENDDO
C
    2 IF (IS.EQ.0) RETURN
C
      DO I = N-IS, 1, -1
            A(I+IS) = A(I)
      ENDDO
      DO I = 1, IS
         A(I) = ' '
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BLANKR (N,A)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Circular shift of character array such that the blanks are to the
C     right side of the array.
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N
      CHARACTER*1 A(N)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
C---------------Executable code--------------------------------------
C
      IS = 0
Caspg
C for some reason N was being corrupted, at least for some
C compilers, so we save its value in ntt
Caspg
      NTT = N
      DO I = 1, N
         IF (A(I).NE.' ') GOTO 2
         IS = IS + 1
      ENDDO
C
    2 IF (IS.EQ.0) RETURN
C
Caspg write (*,*) 'before loop2: i,n,is',i,n,is
      DO I = 1, N-IS
C           write (*,*) 'a(i):',a(i),': <- a(i+is):',a(i+is),':'
            A(I) = A(I+IS)
Caspg       write (*,*) 'inside loop2: i, a(i)',i,a(i)
      ENDDO

Caspg write (*,*) 'before crashing loop: i,ntt,n,is',i,n,ntt,is
      DO I = NTT-IS+1, NTT
C        write (*,*) 'inside crashing loop: i,n,ntt,is',i,n,ntt,is
C        write (*,*) 'a(i) before conversion to whitespace :',a(i),':'
         A(I) = ' '
      ENDDO
C
Caspg write (*,*) 'final line',a
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CONJUGA (N,X,INCX)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C---------------Description--------------------------------------------
C
C     Take complex conjugate of vector X
C
C---------------Routines called----------------------------------------
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher
C
C---------------Calling variables--------------------------------------
C
      INTEGER N,INCX
      COMPLEX*16 X(*)
C
C---------------Common Blocks--------------------------------------
C
C---------------Local variables--------------------------------------
C
      INTEGER II
C
C---------------Executable code--------------------------------------
C
      II = 1
      DO I = 1, N
         X(II) = DCONJG(X(II))
C64B     X(II) =  CONJG(X(II))
         II = II + INCX
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine CpuUsed (TIMCPU)
C     Little wrapper for GETTIM
      Real*8 TIMCPU,TIMWAL
      CALL GETTIM (TIMCPU,TIMWAL)
      Return
      End
C
      Real*8 Function CpuDelta()
C
      Real*8 Sec1,Sec0
      Save Sec0
      Data Sec0 /0.d0/
C
      Call CPUused(Sec1)
      CPUdelta = Sec1 - Sec0
      Sec0 = Sec1
C
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine WallUsed (TIMWAL)
C     Another little wrapper for GETTIM
      Real*8 TIMCPU,TIMWAL
      CALL GETTIM (TIMCPU,TIMWAL)
      Return
      End
C
      Real*8 Function WallDelta()
C
      Real*8 Sec1,Sec0
      Save Sec0
      Data Sec0 /0.d0/
C
      Call WallUsed(Sec1)
      WallDelta = Sec1 - Sec0
      Sec0 = Sec1
C
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      Subroutine DayTime (Datum,Tijd)
      Character Datum*10,Tijd*8

#if defined (SYS_CRAY) || defined (SYS_T3D)
      Call Date(Datum)
      Call Clock(Tijd)
#else
#if defined (SYS_HPUX) || defined (SYS_LINUX) || defined (SYS_DARWIN) \
   || defined (SYS_DEC) || defined (SYS_IRIX)  || defined (SYS_AIX) \
   || defined (SYS_SUN) || defined (SYS_SX)
      CHARACTER FDATE*24, String*24
C
      STRING = FDATE()
C     format from fdate(): 'Fri Sep 13 00:00:00 1986'
      DATUM = ' ' // STRING(9:10) // STRING(4:8) // STRING(23:24)
C     format of DATUM: ' 13 Sep 86'
      TIJD  = STRING(12:19)
C     format of TIJD: '00:00:00'
#else
      Datum = '          '
      Tijd  = '        '
#endif
#endif
      Return
      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      INTEGER FUNCTION INTOWP(N)
#if  defined (INT_STAR8)
      INTOWP = N
#else
      INTOWP=2*N
#endif
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRT_EV (LUPRI,N_IRPS,N_EV,REPNA,E_OFF,ECCIM,IRCW,EV,
     &                   ISIZE,INDX,IRPS,IND,esels,indsel)
C
#include "implicit.h"
#include "codata.h"
C
C---------------Description--------------------------------------------
C
C     Print list of eigenvalues from different symmetries, sorted on
C     increasing energy value.
C
C     An offset is added to this value so that e.g. the constant
C     CCSD nergy can be added to the Fockspace eigenvalues.
C
C     The routine is based on an older routine in GOSCIP but changed
C     to get rid of the common blocks.
C
C     Parameters :
C     - LUPRI   : Output unit
C     - N_IRPS  : Number of irreps
C     - N_EV    : Number of eigenvalues in each irrep
C     - REPNA   : Name of each irrep
C     - E_OFF   : Energy offset (real part)
C     - ECCIM   : Energy offset (imaginary part, for calculations with CAP)
C     - IRCW    : 1 (real values) or 2 (complex energy values)
C     - EV      : Eigenvalues
C     - ISIZE   : size of INDX & IRPS arrays (checked inside)
C     - INDX    : Work array of the same length as EV
C     - IRPS    : Work array of the same length as EV
C     - IND     : Work array of length N_IRPS
C     - esels   : energy of selected (FSCC) state (for numerical gradient runs)
C     - indsel  : index of selected (FSCC) state (for numerical gradient runs)
C
C----------------------------------------------------------------------
C
C   Called from FS_ANALYSIS1
C
C---------------Routines called----------------------------------------
C
C     INDEXX
C
C---------------Last modified------------------------------------------
C
C     Author : Luuk Visscher, august 2003.
C
C     Miro Ilias, July 2007 -  total imaginary eigenenergies are printed
C      if CARITH=T (IRCW,ECCIM parameters)
C
C---------------Calling variables--------------------------------------
C
      DIMENSION N_EV(*),INDX(*)
      DIMENSION EV(IRCW,*), esels(ircw)
      DIMENSION IRPS(*),IND(*)
      CHARACTER*4 REPNA(*)
      REAL*8   , ALLOCATABLE :: EVI(:)
      INTEGER :: IERR, ISIZE, indsel
      LOGICAL DODBG
C
C---------------Local variables--------------------------------------
C
      PARAMETER (DEGEN=1.D-9)
C
C---------------Executable code----------------
C
C     Compute the total number of eigenvalues and fill
C     the array with the identication
C
      N_TOTAL = 0
      DO IRP = 1, N_IRPS
         CALL ICOPY(N_EV(IRP),IRP,0,IRPS(N_TOTAL+1),1)
         N_TOTAL = N_TOTAL + N_EV(IRP)
      ENDDO
C
      IF (N_TOTAL.EQ.0) RETURN
C
C     Sort them from low to high
C
      DODBG = .FALSE.

      IF (DODBG) THEN
      write(LUPRI,*) 'before sorting'
      DO I =1, N_TOTAL
       IF (IRCW.EQ.1) THEN
         write(LUPRI,*) 're:',EV(1,I)
       ELSE IF (IRCW.EQ.2) THEN
         write(LUPRI,*) 're:',EV(1,I),' im:',EV(2,I)
       ENDIF
      ENDDO
      write(LUPRI,*)
      ENDIF

      call FLSHFO(LUPRI)

      IF (IRCW.EQ.2) THEN
        ALLOCATE(EVI(N_TOTAL),STAT=IERR)
        IF (IERR.NE.0) THEN
          CALL QUIT('PRT_EV: EVI wrong allocation !')
        ENDIF
        DO I=1,N_TOTAL
         EVI(I) = EV(1,I)
        ENDDO
        CALL INDEXX(N_TOTAL,EVI,INDX)
        DEALLOCATE(EVI,STAT=IERR)
        IF (IERR.NE.0) THEN
          CALL QUIT('PRT_EV: EVI wrong deallocation !')
        ENDIF
      ELSE IF (IRCW.EQ.1) THEN
        CALL INDEXX(N_TOTAL,EV(1,1),INDX)
      ELSE
        write(LUPRI,*) 'PRT_EV: wrong parameter IRCW=',IRCW
        CALL QUIT('PRT_EV...wrong parameter !!! ')
      ENDIF

      IF (DODBG) THEN
      write(LUPRI,*) 'after sorting'
CMI ... own printout
      IF (IRCW.EQ.2) THEN
       DO I=1,N_TOTAL
         write(LUPRI,*) 're:',EV(1,INDX(I)),' im:',EV(2,INDX(I))
       ENDDO
      ELSE IF (IRCW.EQ.1) THEN
       DO I=1,N_TOTAL
         write(LUPRI,*) 're:',EV(1,INDX(I))
       ENDDO
      ENDIF
      ENDIF
C
      WRITE(LUPRI,'(/A,F14.10,A,F17.9,A)')
     & ' ( 1 au =',XTEV,' eV  / ',XTKAYS,' cm-1)' 

      IF (IRCW.EQ.1) THEN
        WRITE (LUPRI,1005)
      ELSE IF (IRCW.EQ.2) THEN
        WRITE (LUPRI,1007)
      ELSE
        WRITE(LUPRI,*) 'PRT_EV: wrong IRCW=',IRCW
        CALL QUIT('PRT_EV: wrong IRCW value !')
      ENDIF
C
      ILEVEL=0
      EW=0.0D0
      EWIM = 0.0D0
      EWF=EV(1,INDX(1))
      NEW=0
      AVER=0.0D0
      AVERIM=0.0D0

!     save energy for selected state (if any)
      if(indsel > 0)then
                      esels(1)    = EV(1,INDX(indsel)) + E_OFF
        IF(IRCW.EQ.2) esels(ircw) = EV(2,INDX(indsel)) + ECCIM
      end if

      DO I=1,N_TOTAL
        E = EV(1,INDX(I))
        IF (IRCW.EQ.2) EIM = EV(2,INDX(I))

        IF (DODBG) THEN
         write(lupri,*) 'I,NEW=',I,NEW,' E=',E,' E-EW=',E-EW
         IF (IRCW.EQ.2) write(lupri,*) 'EIM=',EIM,' EIM-EWIM=',EIM-EWIM
        ENDIF

        AVER=AVER+E
        IF (IRCW.EQ.2) AVERIM=AVERIM+EIM

       IF      (IRCW.EQ.1) THEN
        IF ( ABS(E-EW).GT.DEGEN) THEN
          IF (NEW.NE.0) THEN
             WRITE(LUPRI,1010) ILEVEL,EW-EWF,EW,EW+E_OFF,NEW
          ENDIF
          ILEVEL=ILEVEL+1
          NEW=1
          EW=E
        ELSE
          NEW=NEW+1
          EW=((NEW-1)*EW+E)/DBLE(NEW)
        ENDIF  

       ELSE IF (IRCW.EQ.2) THEN

        IF ( ABS(E-EW).GT.DEGEN.OR.ABS(EIM-EWIM).GT.DEGEN ) THEN  
          IF (NEW.NE.0) THEN
            WRITE(LUPRI,1015) ILEVEL,EW-EWF,EW,EW+E_OFF,
     &      EWIM+ECCIM, NEW
          ENDIF
          ILEVEL=ILEVEL+1
          NEW=1
          EW=E
          EWIM = EIM
        ELSE
          NEW=NEW+1
          EW=((NEW-1)*EW+E)/DBLE(NEW)
          EWIM=((NEW-1)*EWIM+EIM)/DBLE(NEW)
        ENDIF

       ELSE
         call quit('prt_ev blind branch !')
       ENDIF

      ENDDO
      IF (NEW.NE.0) THEN
        IF (IRCW.EQ.1) THEN
          WRITE(LUPRI,1010) ILEVEL,EW-EWF,EW,EW+E_OFF,NEW
        ELSE IF (IRCW.EQ.2) THEN
          WRITE(LUPRI,1015) ILEVEL,EW-EWF,EW,EW+E_OFF,
     &    EWIM+ECCIM,NEW
C    &    EV(2,INDX(I)),NEW
        ENDIF
      ENDIF

      WRITE(LUPRI,1020) AVER/DBLE(N_TOTAL)+E_OFF

CMI  ... printout with symmetry classification ...real energies only...
      DO IRP1 = 1, N_IRPS, 8
      IRP2 = MIN0(N_IRPS,IRP1+7)
      WRITE(LUPRI,1025) (REPNA(I),I=IRP1,IRP2)
      ILEVEL=0
      EW=0.0D0
      EWF=EV(1,INDX(1))
      NEW=0
      DO JRP = IRP1, IRP2
         IND(JRP) = 0
      END DO

      DO 140 I=1,N_TOTAL

!MI: check index range 
        IF (I.LT.1.OR.I.GT.ISIZE) THEN
          WRITE(LUPRI,*) "I,ISIZE=",I,ISIZE,"...INDX(I)=",INDX(I)
          WRITE(LUPRI,*) "N_TOTAL=",N_TOTAL 
          CALL FLSHFO(LUPRI)
          CALL QUIT('Index I in in INDX(I) out of bounds !')
        ENDIF
 
        E=EV(1,INDX(I))

        IF (ABS(E-EW).GT.DEGEN) THEN
          IF (NEW.NE.0)
     &      WRITE(LUPRI,1030) ILEVEL,
     &                     XTEV*(EW-EWF),XTKAYS*(EW-EWF),
     &                    (IND(JRP),JRP=IRP1,IRP2)
          DO JRP = IRP1, IRP2
            IND(JRP) = 0
          ENDDO
          ILEVEL = ILEVEL + 1

          IND(IRPS(INDX(I))) = 1

          NEW=1
          EW=E
        ELSE
          NEW=NEW+1
          IND(IRPS(INDX(I))) = IND(IRPS(INDX(I))) + 1
          EW=((NEW-1)*EW+E)/DBLE(NEW)
        ENDIF

 140    CONTINUE
      IF (NEW.NE.0) THEN

         WRITE(LUPRI,1030) ILEVEL,
     &   XTEV*(EW-EWF),XTKAYS*(EW-EWF),
     &   (IND(JRP),JRP=IRP1,IRP2)

      ENDIF

      ENDDO

 1000 FORMAT(F17.8)
 1008 FORMAT(17X,I8,2X,16(A1),F16.8)
 1005 FORMAT(//' Energy eigenvalues in atomic units'
     +            //' Level   Rel eigenvalue     ',
     +            'Abs eigenvalue      Total Energy    Degeneracy'/)
 1010 FORMAT(I5,F17.10,F20.12,'  ',F20.12,' (',I4,' * )')
 1007 FORMAT(//' Energy eigenvalues in atomic units'
     & //' Level   Rel eigenvalue     Abs eigenvalue     ',
     & ' Total Energy       Total imag.ener.  Degeneracy'/)
 1015 FORMAT(I5,F17.10,F20.12,'  ',F20.12,' ',F20.15,' (',I4,' * )')
 1020 FORMAT(//' Total average: ',F20.10)
 1025 FORMAT(//' Relative real eigenvalues in other units;'/' Symmetry',
     +       ' Classification in the Abelian subgroup',
     +       //' Level  eigenvalue (eV) ',
     +      ' Eigenvalue (cm-1)',2X,8(A4,'|')/)
 1030 FORMAT(I5,F19.9,F18.6,2X,8(I4,'|'))

      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE PRT_EV_VEC (LUPRI,N_IRPS,N_EV,REPNA,E_OFF,EV,
     &                       INDX,IRPS,IND,VEC,LBL_SEC1,LBL_SEC2,
     &                       IRSEC1,IRSEC2,IRSEC12,REPNAO,NACTV,
     &                       TSHOLD)
C
#include "implicit.h"
#include "codata.h"
C
C---------------Description--------------------------------------------
C
C     Print list of eigenvalues from different symmetries, sorted on
C     increasing energy value.
C
C     An offset is added to this value so that e.g. the constant
C     CCSD energy can be added to the Fockspace eigenvalues.
C
C     The routine is based on an older routine in GOSCIP but changed
C     to get rid of the common blocks.
C
C     Parameters :
C     - LUPRI   : Output unit
C     - N_IRPS  : Number of irrep
C     - N_EV    : Number of eigenvalues in each irrep
C     - REPNA   : Name of each irrep after adding particles
C     - REPNAO  : Name of each irrep before adding particles
C     - E_OFF   : Energy offset
C     - EV      : Eigenvalues
C     - VEC     : Eigenvectors
C     - NACTV   : Number of active virtual for each irrep
C     - LBL_SEC1: Sector 1 spinor for eigenvectors in each irrep
C     - LBL_SEC2: Sector 2 spinor for eigenvectors in each irrep
C     - IRSEC1  : Sector 1 irrep for eigenvectors in each irrep
C     - IRSEC2  : Sector 2 irrep for eigenvectors in each irrep
C     - IRSEC12 : Bosono mult sector 1  sector 2
C     - INDX    : Work array of the same length as EV
C     - IRPS    : Work array of the same length as EV
C     - IND     : Work array of length N_IRPS
C
C---------------Routines called----------------------------------------
C
C     INDEXX
C
C---------------Last modified------------------------------------------
C
C     Author : Ivan Infante, september 2005
C
C---------------Calling variables--------------------------------------
C
      DIMENSION N_EV(*),INDX(*), INDXV(600)
      DIMENSION EV(*), VEC(*), VEC_TMP(600*600)
      DIMENSION LBL_SEC1(*), LBL_SEC2(*)
      DIMENSION IRSEC1(*), IRSEC2(*), IRSEC12(*)
      DIMENSION LSEC1(600*600), LSEC2(600*600)
      DIMENSION IRPSEC1(600*600), IRPSEC2(600*600), IRPSEC12(600*600)
      DIMENSION N_TOT_EV(64), N_TOT_VEC(64), NACTV(*)
      DIMENSION IRPS(*),IND(*)
      CHARACTER*4 REPNA(*), REPNAO(*)
      INTEGER TSHOLD
      REAL*8 E_OFF
C
C---------------Local variables--------------------------------------
C
      PARAMETER (DEGEN=1.D-9)
C
C---------------Executable code----------------
C
C     Compute the total number of eigenvalues and fill
C     the array with the identication
C     Select total number of eigenvectors to use like offset later
C
      N_TOTAL = 0
      N_TOT_EV(1) = 0
      N_TOT_VEC(1) = 0
      DO IRP = 1, N_IRPS
         CALL ICOPY(N_EV(IRP),IRP,0,IRPS(N_TOTAL+1),1)
         N_TOTAL = N_TOTAL + N_EV(IRP)
         N_TOT_EV(IRP+1) = N_TOT_EV(IRP) + N_EV(IRP)
         N_TOT_VEC(IRP+1) = N_TOT_VEC(IRP) + N_EV(IRP)* N_EV(IRP)
      ENDDO
C
      IF (N_TOTAL.EQ.0) RETURN
C
C
C
C     Sort them from low to high
C
      CALL INDEXX(N_TOTAL,EV,INDX)
C
CMI ... print out conversion factors (from codata)
      WRITE(LUPRI,'(/A,F14.10,A,F17.9,A)')
     & ' ( 1 au =',XTEV,' eV  / ',XTKAYS,' cm-1)' 


      WRITE (LUPRI,1005)
      ILEVEL=0
      EW=0.0D0
      EWF=EV(INDX(1))
      NEW=0
      AVER=0.0D0
      DO I=1,N_TOTAL
        E = EV(INDX(I))
        AVER=AVER+E
        IF (ABS(E-EW).GT.DEGEN) THEN
          IF (NEW.NE.0)
     &      WRITE(LUPRI,1010) ILEVEL,EW-EWF,EW,EW+E_OFF,NEW
          ILEVEL=ILEVEL+1
          NEW=1
          EW=E
        ELSE
          NEW=NEW+1
          EW=((NEW-1)*EW+E)/DBLE(NEW)
        ENDIF
      ENDDO
      IF (NEW.NE.0)
     &WRITE(LUPRI,1010) ILEVEL,EW-EWF,EW,EW+E_OFF,NEW
      WRITE(LUPRI,1020) AVER/DBLE(N_TOTAL)+E_OFF
C
      DO IRP1 = 1, N_IRPS, 8
        IRP2 = MIN0(N_IRPS,IRP1+7)
        WRITE(LUPRI,1025) (REPNA(I),I=IRP1,IRP2)
        ILEVEL=0
        EW=0.0D0
        EWF=EV(INDX(1))
        NEW=0
        DO JRP = IRP1, IRP2
           IND(JRP) = 0
        END DO
        DO 140 I=1,N_TOTAL
          E=EV(INDX(I))
          IVEC = (INDX(I)-N_TOT_EV(IRPS(INDX(I)))-1)*N_EV(IRPS(INDX(I)))
     &            + N_TOT_VEC(IRPS(INDX(I))) + 1
          IF (ABS(E-EW).GT.DEGEN) THEN
             IF (NEW.NE.0) THEN
               WRITE(LUPRI,1030) ILEVEL,
     &                       XTEV*(EW-EWF),XTKAYS*(EW-EWF),
     &                      (IND(JRP),JRP=IRP1,IRP2)
             ENDIF
             DO JRP = IRP1, IRP2
                IND(JRP) = 0
             ENDDO
             ILEVEL = ILEVEL + 1
             IND(IRPS(INDX(I))) = 1
             NEW=1
             EW=E
          ELSE
             NEW=NEW+1
             IND(IRPS(INDX(I))) = IND(IRPS(INDX(I))) + 1
             EW=((NEW-1)*EW+E)/DBLE(NEW)
          ENDIF
 140     CONTINUE
         IF (NEW.NE.0) THEN
             WRITE(LUPRI,1030) ILEVEL,
     &                XTEV*(EW-EWF),XTKAYS*(EW-EWF),
     &                (IND(JRP),JRP=IRP1,IRP2)
         ENDIF
      ENDDO
C
C     New part (I.I.):
C     Print composition of each excited state
C     PS: It should be generalized for other sectors
C
      ILEVEL=1
      EWF=EV(INDX(1))
      LBL = 0
      WRITE(LUPRI,1080)
      DO IRP = 1, N_IRPS
         DO I = 1, NACTV(IRP)
            LBL = LBL + 1
            WRITE(LUPRI,1085) LBL, REPNAO(IRP)
         ENDDO
      ENDDO
      DO I=1,N_TOTAL
        E = EV(INDX(I))
        IVEC = (INDX(I)-N_TOT_EV(IRPS(INDX(I)))-1)*N_EV(IRPS(INDX(I)))
     &          + N_TOT_VEC(IRPS(INDX(I))) + 1
        WRITE(LUPRI,1055) ILEVEL,E-EWF,XTEV*(E-EWF),XTKAYS*(E-EWF)
        WRITE(LUPRI,1090) E_OFF + E
        CALL DCOPY(N_EV(IRPS(INDX(I))),VEC(IVEC),1,VEC_TMP,1)
        CALL ICOPY(N_EV(IRPS(INDX(I))),LBL_SEC1(IVEC),1,LSEC1,1)
        CALL ICOPY(N_EV(IRPS(INDX(I))),LBL_SEC2(IVEC),1,LSEC2,1)
        CALL ICOPY(N_EV(IRPS(INDX(I))),IRSEC1(IVEC),1,IRPSEC1,1)
        CALL ICOPY(N_EV(IRPS(INDX(I))),IRSEC2(IVEC),1,IRPSEC2,1)
        CALL ICOPY(N_EV(IRPS(INDX(I))),IRSEC12(IVEC),1,IRPSEC12,1)
        WRITE(LUPRI,1060)
        DO J = 1, N_EV(IRPS(INDX(I)))
          VEC_TMP(J) =  VEC_TMP(J)*VEC_TMP(J)*100
          IF (VEC_TMP(J).GT.TSHOLD) THEN
           WRITE(LUPRI,1065) LSEC2(J),REPNAO(IRPSEC2(J)),
     &                       LSEC1(J),REPNAO(IRPSEC1(J)),
     &                       REPNA(IRPSEC12(J)),VEC_TMP(J)
           ENDIF
        ENDDO
        ILEVEL=ILEVEL+1
      ENDDO
C
 1000 FORMAT(F17.8)
 1008 FORMAT(17X,I8,2X,16(A1),F16.8)
 1005 FORMAT(//' Energy eigenvalues in atomic units'
     +            //' Level   Rel eigenvalue     ',
     +            'Abs eigenvalue      Total Energy    Degeneracy'/)
 1010 FORMAT(I5,F17.12,F20.12,'  ',F20.12,' (',I4,' * )')
 1015 FORMAT(I5,F17.12,F19.9,F18.6)
 1020 FORMAT(//' Total average: ',F20.10)
 1025 FORMAT(//' Relative eigenvalues in other units;'/' Symmetry',
     +       ' Classification in the Abelian subgroup',
     +       //' Level  eigenvalue (eV) ',
     +      ' Eigenvalue (cm-1)',2X,8A4/)
 1030 FORMAT(I5,F19.9,F18.6,2X,8I4)
 1035 FORMAT(/'Real Vector:',/(1x,6f12.7))
 1055 FORMAT(//'Energy ',I5,':',F12.9,' (au)',3X,F12.5,' (eV)',3X,F19.5
     +       ,' (cm-1)')
 1060 FORMAT(/,6x,'Composition:'/,11x,'Spinor 2',2x,'Rep',6x,'Spinor 1',
     +       2x,'Rep',5x,'Tot Rep',3x,'weight (%)')
 1065 FORMAT(10x,I5,4x,A4,7x,I5,4x,A4,5x,A4,5x,F10.5)
 1080 FORMAT(/,2x,'Composition of all excited states in (0h,2p) sector'
     +        /'  Active Spinor',5x,'Abelian'/)
 1085 FORMAT(4x,I5,9x,A6)
 1090 FORMAT('Total CCSD energy + affinity: ',F20.12)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTIME (TIMERT,FUNCT,NAME_IN)
C
C     Initial routine written by M. Pernpointner
C     Documented and revised by L. Visscher
C
C     Input parameters
C         - TIMERT : Type of timer (Used in MPE calls
C                                   0:Do not log in MPE
C                                   1:Initialization
C                                   2:Wait,
C                                   3:Send/receive task,
C                                   4:Send/receive Data,
C                                   5:Perform computation)
C         - FUNCT  : TCLEAR=0 : Set the timer to zero
C                    TSTART=1 : Start adding time to this timer
C                    TSTOP=2  : Stop the timer
C                    TPRNTI=3 : Write contents of the timer to output
C                    TPRNTA=4 : Write contents of all timers to output
C                    TINIT=-1 : Initialize all timers
C         - NAME   : Name of the timer (max 30 characters)
C
      use interface_to_mpi
      IMPLICIT INTEGER (A-Z)
      CHARACTER*(*) NAME_IN
      CHARACTER*30  NAME, T_REPORT
      PARAMETER ( D0 = 0.0D0 )
#include "timtbl.h"
#include "priunit.h"
#include "infpar.h"
      REAL*8 A,MY_T,CPU_AV,CPU_MI,CPU_MA
C
      NAME = NAME_IN   ! extend with blanks if NAME_IN shorter than 30 chars
C
C     Check first for operations that involve all timers
C
      IF (FUNCT.EQ.TINIT) THEN
C
C        Initialization of all timers
C
         DO I = 1, MXTIMF
          CPUTF(I)=D0
          CPUTFT(I)=D0
          CPUTT(I)=TIMERF_NAME
         ENDDO
         LAST_MPE_ACTION = 0
         LAST_MPE_STATE  = 0
         RETURN
      ELSEIF (FUNCT.EQ.TPRNTA) THEN
C
C        Print out of all timers
C
        IF (MYTID.EQ.MPARID) WRITE(LUPRI,220) NAME
#if ! defined (VAR_MPI)
C       Serial code : no need to communicate with slaves
        DO I=1,MXTIMF
          IF(CPUTFT(I).NE.D0) WRITE(LUPRI,222) CPUTT(I),CPUTFT(I)
        ENDDO
#else
C
C       We may have the parallel version running in serial mode,
C       handle this first (same code as serial)
C
        IF (NUMNOD.EQ.0) THEN
          DO I=1,MXTIMF
            IF(CPUTFT(I).NE.D0) WRITE(LUPRI,222) CPUTT(I),CPUTFT(I)
          ENDDO
        ELSE
C
C        The parallel version runs in parallel. We need to communicate
C        to gather all data on the master.
C
C        Run over all master timers and collect data from slaves
C        NB : This also gets the inactive timers, but the loop structure
C        is easier this way
C
         IF (MYTID.EQ.MPARID)  WRITE(LUPRI,330)
C
         DO I = 1, MXTIMF
C          Tell the slaves to look for timer I
           T_REPORT = CPUTT(I)
           call interface_mpi_BCAST(T_REPORT,30,
     &                              MPARID,global_communicator)
C          Find the timer (we do this code on master also for simplicity)
           TIMERN = -1
           DO J = 1, MXTIMF
              IF (T_REPORT.EQ.CPUTT(J)) TIMERN = J
           ENDDO
           IF (TIMERN.EQ.-1) THEN
C             Timer not used, give zero values
              MY_T = D0
           ELSE
              MY_T = CPUTFT(TIMERN)
           ENDIF
C          NB : The following computations are only correct on the master,
C          but since we report only there there is no need for an ALLREDUCE call.
C          Compute the average time spent by the slaves
           IF (MYTID.EQ.MPARID) THEN
              A = D0
           ELSE
              A = MY_T
           ENDIF
           call interface_mpi_REDUCE(A,CPU_AV,1,op_MPI_SUM,MPARID,
     &                               global_communicator)
           CPU_AV = (CPU_AV)/DBLE(NUMNOD)
C          Compute the minimum time
           IF (MYTID.EQ.MPARID) A = 1.D8
           call interface_mpi_REDUCE(A,CPU_MI,1,
     &                      op_MPI_MIN,MPARID,global_communicator)
C          Compute the maximum time
           IF (MYTID.EQ.MPARID) A = -1.D8
           call interface_mpi_REDUCE(A,CPU_MA,1,
     &                      op_MPI_MAX,MPARID,global_communicator)
C
C          Write the result to the output
C
           IF (MYTID.EQ.MPARID.AND.(CPU_AV.NE.D0.OR.MY_T.NE.D0))
     &        WRITE(LUPRI,340) T_REPORT,MY_T,CPU_MI,CPU_MA,CPU_AV
C
         ENDDO
C
C        The slaves may also have timers that master does not have, let
C        the primus inter pares report.
C
         DO I = 1, MXTIMF
           SPARTACUS = MPARID+1
           T_REPORT = CPUTT(I)
           call interface_mpi_BCAST(T_REPORT,30,
     &                     SPARTACUS,global_communicator)
           TIMERN = -1
           DO J = 1, MXTIMF
              IF (T_REPORT.EQ.CPUTT(J)) TIMERN = J
           ENDDO
           IF (TIMERN.EQ.-1) THEN
C             Timer not used, give zero values
              MY_T = D0
           ELSE
              MY_T = CPUTFT(TIMERN)
           ENDIF
           IF (MYTID.EQ.MPARID) THEN
              A = D0
           ELSE
              A = MY_T
           ENDIF
C          Compute the average time spent by the slaves
           call interface_mpi_REDUCE(A,CPU_AV,1,
     &                      op_MPI_SUM,MPARID,global_communicator)
           CPU_AV = (CPU_AV)/DBLE(NUMNOD)
C          Compute the minimum time
           IF (MYTID.EQ.MPARID) A = 1.D8
           call interface_mpi_REDUCE(A,CPU_MI,1,
     &                      op_MPI_MIN,MPARID,global_communicator)
C          Compute the maximum time
           IF (MYTID.EQ.MPARID) A = -1.D8
           call interface_mpi_REDUCE(A,CPU_MA,1,
     &                      op_MPI_MAX,MPARID,global_communicator)
C
C          Write the result to the output
C          The master has to have TIMERN -1, otherwise this timer
C          has already been reported in the previous loop.
C
           IF (MYTID.EQ.MPARID.AND.TIMERN.EQ.-1.AND.CPU_AV.NE.D0)
     &        WRITE(LUPRI,340) T_REPORT,MY_T,CPU_MI,CPU_MA,CPU_AV
        ENDDO
      ENDIF ! end of 1 / multiple processor clause
      IF (MYTID.EQ.MPARID) WRITE(LUPRI,350)
#endif
      RETURN
      ENDIF ! end of all timers IF-THEN-ELSE clause
C
C     We come here only if we work on a single timer, find out which one
C
      TIMERN = -1
      DO I = 1, MXTIMF
         IF (NAME.EQ.CPUTT(I)) TIMERN = I
      ENDDO
C
C     We may have a new timer, look for free slots
C
      IF (TIMERN.EQ.-1) THEN
         TIMERF = 0
         DO I = MXTIMF,1,-1
            IF (CPUTT(I).EQ.TIMERF_NAME) TIMERF = I
         ENDDO
C
         IF (TIMERF.NE.0) THEN
C           Use the first free slot and initialize
            CPUTT(TIMERF)=NAME
            CPUTFT(TIMERF)=D0
            TIMERN = TIMERF
         ELSE
C           Too bad, we ran out of slots
            WRITE(LUPRI,130) MYTID,NAME
            RETURN
         ENDIF
      ENDIF
C
C     Perform the desired action
C
      IF(FUNCT.EQ.TCLEAR) THEN
          CPUTF(TIMERN)=D0
          CPUTFT(TIMERN)=D0
          CPUTT(TIMERN)=TIMERF_NAME
      ELSE IF (FUNCT.EQ.TSTART) THEN
        CALL CPUUSED(A)
        CPUTF(TIMERN)=A
        CPUTT(TIMERN)=NAME
#if defined (MPE)
        IF (TIMERT.NE.0) THEN
C          We may need to end the previous state
           IF (TIMERT.NE.LAST_MPE_STATE.AND.LAST_MPE_ACTION.EQ.1) THEN
              I_MPE = 2 * LAST_MPE_STATE
              CALL MPE_LOG_EVENT(I_MPE,0,NAME)
           ENDIF
C          Announce the new state if necessary
           IF (TIMERT.NE.LAST_MPE_STATE.OR.LAST_MPE_ACTION.NE.1) THEN
              I_MPE = 2 * TIMERT - 1
              CALL MPE_LOG_EVENT(I_MPE,0,NAME)
              LAST_MPE_STATE  = TIMERT
              LAST_MPE_ACTION = 1
           ENDIF
        ENDIF
#endif
      ELSE IF (FUNCT.EQ.TSTOP) THEN
        CALL CPUUSED(A)
        CPUTF(TIMERN)=A-CPUTF(TIMERN)
        CPUTFT(TIMERN)=CPUTFT(TIMERN)+CPUTF(TIMERN)
#if defined (MPE)
        IF (TIMERT.NE.0) THEN
C          The state to be ended may be interrupted and ended by a start
C          Only do something if everything matches
           IF (TIMERT.EQ.LAST_MPE_STATE.AND.LAST_MPE_ACTION.EQ.1) THEN
              I_MPE = 2 * TIMERT
              CALL MPE_LOG_EVENT(I_MPE,0,NAME)
              LAST_MPE_STATE  = TIMERT
              LAST_MPE_ACTION = 2
           ENDIF
        ENDIF
#endif
      ELSE IF (FUNCT.EQ.TPRNTI) THEN
        IF (MYTID.EQ.MPARID) THEN
           IF(CPUTF(TIMERN).NE.D0) THEN
             WRITE(LUPRI,222) CPUTT(TIMERN),CPUTFT(TIMERN)
           ELSE
             WRITE(LUPRI,224) TIMERN
           ENDIF
        ENDIF
      ELSE
        WRITE(LUPRI,135)
      ENDIF
C
 130  FORMAT(/' **WARNING** No free timer slots for proces',I4/,
     &        ' ** no action taken')
 135  FORMAT(/' **WARNING** Unknown timer function'/,
     &        ' ** no action taken')
 220  FORMAT(/'------ Timing report (in CPU seconds) of module ',A/)
 222  FORMAT(' Time in ',A30,T40,F12.3,' seconds')
 223  FORMAT(T37,'-------------',/,T25,'Sum: ',F12.3,' seconds',
     &      /T37,'-------------')
 224  FORMAT(' Timer ',I3,' not used for measurements')
 330  FORMAT(//,T37,'Master',T46,'--------- Slaves ------------'//,
     &          T2,'Timer',T37,'Master',T46,'Minimum',T56,'Maximum',
     &          T68,'Average'/)
 340  FORMAT(T2,A30,T33,F10.1,T43,F10.1,T53,F10.1,T63,F10.1)
 350  FORMAT(/'------ End of timing report ------'/)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE XTIMRD (ELTIM,NAME)
C
C  purpose:  read out a specific timer previously set by XTIME
C
      IMPLICIT INTEGER (A-Z)
      CHARACTER*30 NAME
      REAL*8 ELTIM
      PARAMETER ( D0 = 0.0D0 )
#include "timtbl.h"
#include "priunit.h"
      TIMERN = -1
      DO I = 1, MXTIMF
         IF (NAME.EQ.CPUTT(I)) TIMERN = I
      ENDDO
      IF(TIMERN.LT.1.OR.TIMERN.GT.MXTIMF) THEN
        ELTIM=D0
      ELSE
        ELTIM=CPUTFT(TIMERN)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LFNAM */
      CHARACTER*10 FUNCTION LFNAM (FILNAM)
C
C     Add the node identification to the file name for shared disk systems
C     Written by Luuk Visscher, February 1998
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "dcbgen.h"
#include "infpar.h"
      CHARACTER*5  FILNAM
      CHARACTER*10 LFILNAM
      NODE = MYTID
C     ... MYTID will be zero if this is a sequential version of Dirac
C
 1000 FORMAT (A5,I5)
C
      WRITE (LFILNAM,1000) FILNAM,NODE
      DO I = 6, 10
         IF (LFILNAM(I:I).EQ.' ') LFILNAM(I:I) = 'X'
      ENDDO
C
      LFNAM = LFILNAM
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_INFO (NBUF,NSIZE)
C
C     Report the maximum number of buffer and record size that
C     will be used in the sort
C
C     Calling parameters
C
C     NBUF   : Maximum number of buffers (hardwired in the common block)
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
C
      NBUF = N_BUFFER
      NSIZE = NGBFSZ
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_INIT (I_FILE)
C
C     Initialize a buffer file to transpose a large matrix that is generated in parts
C     Currently used in MOLTRA for the halftransformed two-electron
C     integrals but it may be useful at other places as well.
C
C     Calling parameters
C
C     I_FILE   : Bucket sort indentifier (we may have more than one active)
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
      CHARACTER*10 FNODE,LFNAM
      CHARACTER*5  NAME
C
C     The record size needs to be NGBFSZ real + (2*NGBFSZ + 2) integer
C     because we also need to store the labels
C
      IRECL = IRECLEN(NGBFSZ+1,2*NGBFSZ+2,0)
C
C     The opening of a file is a bit tricky in FORTRAN, pray that the unit
C     is not yet in use....
C
      WRITE (NAME,'(A3,I2.2)') 'BUF',I_FILE
      FNODE = LFNAM(NAME)
      LGFIL = I_BF_UNIT+I_FILE

#if !defined (VAR_PFS)
      OPEN (UNIT=LGFIL,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
#endif
C
C     Initialize the data that is kept in the common block
C
      I_BF_REC(I_FILE) = 0
      BF_NAME(I_FILE) = FNODE
      CALL IZERO(LGREC(1,I_FILE),N_BUFFER)
      CALL IZERO(LGBUF(1,I_FILE),N_BUFFER)
      N_WRITTEN = 0
      N_DELETED = 0
C
#if defined (VAR_MPI)
      DO I = 1, MAX_SEND_BUFFERS
         IREQ(I) = df_MPI_REQUEST_NULL
      ENDDO
#endif
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_OUT (I_FILE,RGBUF,IGBUF,IROWS,
     &                            A,NROWS,ICOL,TRESHOLD)
C
C     Scatter row of matrix A over the buffers, flush them if necessary
C
C     Calling parameters
C
C     I_FILE   : File unit
C     RGBUF    : Real data in buffers
C     IGBUF    : Integer data data in buffers
C     IROWS    : Pointer from the rows to buffer
C     A        : Column of matrix that is to be buffered away
C     NROWS    : Number of rows in A
C     ICOL     : Column index of A
C     TRESHOLD : Treshold for writing the matrix element
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
      REAL*8 RGBUF(NGBFSZ,*), A(NROWS)
      DIMENSION IGBUF(2,NGBFSZ,*)
      DIMENSION IROWS(3,NROWS)
C
C     Loop through the row and check if the data need be stored
C
      DO I = 1, NROWS
         IBUF = IROWS(1,I)
         IF (IBUF.GT.0.AND.ABS(A(I)).GT.TRESHOLD) THEN
            N_WRITTEN = N_WRITTEN + 1
            N = LGBUF(IBUF,I_FILE) + 1
            RGBUF(N,IBUF) = A(I)
            IGBUF(1,N,IBUF) = ICOL
            IGBUF(2,N,IBUF) = IROWS(2,I)
            IF (N.GE.NGBFSZ) THEN
               INODE = IROWS(3,I)
#if defined (VAR_PFS)
               CALL GLOBAL_FILE_WRITE (
     &         RGBUF(1,IBUF),IGBUF(1,1,IBUF),N,MYTID,INODE,IBUF)
#else
#if defined (VAR_MPI)
               IF (INODE.NE.MYTID) THEN
                   CALL BUFFER_FILE_SEND (I_FILE,IBUF,RGBUF(1,IBUF),
     &              IGBUF(1,1,IBUF),N,INODE)
               ELSE
#endif
                   CALL BUFFER_FILE_WRITE (I_FILE,IBUF,RGBUF(1,IBUF),
     &              IGBUF(1,1,IBUF),N)
#if defined (VAR_MPI)
               ENDIF
#endif
#endif
               N = 0
            ENDIF
            LGBUF(IBUF,I_FILE) = N
         ELSE
            N_DELETED = N_DELETED + 1
         ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_WRITE (I_FILE,IBUF,RGBUF,IGBUF,N)
C
C     Write local integrals
C
C     Calling parameters
C
C     I_FILE   : File unit
C     IBUF     : Buffer number
C     RGBUF    : One buffer with real data
C     IGBUF    : One buffer with integer data
C     N        : Actual size of the buffer
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
      REAL*8 RGBUF(NGBFSZ)
      DIMENSION IGBUF(2,NGBFSZ)
C
C     Get file unit number and file pointer
C
      LGFIL    = I_BF_UNIT+I_FILE
      IREC     = I_BF_REC(I_FILE) + 1
      WRITE (LGFIL,REC=IREC) LGREC(IBUF,I_FILE),N,RGBUF,IGBUF
C
C     Update file pointers
C
      LGREC(IBUF,I_FILE) = IREC
      I_BF_REC(I_FILE) = IREC
      records_in_buffer(ibuf) = records_in_buffer(ibuf)+1      
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_FLUSH (I_FILE,RGBUF,IGBUF,IROWS,NROWS)
C
C     Flush the buffers
C
C     Calling parameters
C
C     LGFIL    : File unit
C     RBUF     : Real data in buffers
C     IBUF     : Integer data data in buffers
C     IROWS    : Pointer from rows to buffer and node
C     NROWS    : Number of rows
C
C     Written by Luuk Visscher, August 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
      REAL*8 RGBUF(NGBFSZ,*)
      DIMENSION IGBUF(2,NGBFSZ,*)
      DIMENSION IROWS(3,NROWS)
C
C     Loop through the rows and check whether the data is to be flushed
C
      N=0
      DO I = 1, NROWS
         IBUF = IROWS(1,I)
         IF (IBUF.NE.0) N = LGBUF(IBUF,I_FILE)
         IF (IBUF.GT.0.AND.N.GT.0) THEN
            INODE = IROWS(3,I)
#if defined (VAR_PFS)
             CALL GLOBAL_FILE_WRITE (
     &        RGBUF(1,IBUF),IGBUF(1,1,IBUF),N,MYTID,INODE,IBUF)
#else
#if defined (VAR_MPI)
            IF (INODE.NE.MYTID) THEN
                CALL BUFFER_FILE_SEND (I_FILE,IBUF,RGBUF(1,IBUF),
     &                                 IGBUF(1,1,IBUF),N,INODE)
            ELSE
#endif
                CALL BUFFER_FILE_WRITE (I_FILE,IBUF,RGBUF(1,IBUF),
     &                                 IGBUF(1,1,IBUF),N)
#if defined (VAR_MPI)
            ENDIF
#endif
#endif
            N = 0
         ENDIF
         IF (IBUF.NE.0) LGBUF(IBUF,I_FILE) = N
      ENDDO
C
      IF((N_WRITTEN+N_DELETED).EQ.0) THEN
        WRITE(LUPRI,'(A,I5,A)') 
     &    ' Node',MYTID,' slept during first half transformation.'
      ELSE
        PERC = N_WRITTEN+N_DELETED
        PERC = 100.D0/PERC
        PERC = PERC*N_WRITTEN
#if defined (INT_STAR8)
        GIGABYTES = N_WRITTEN * 24.0
#else
        GIGABYTES = N_WRITTEN * 16.0
#endif
        GIGABYTES = GIGABYTES / (1024.0**3)
        WRITE (LUPRI,1000) MYTID,N_WRITTEN,PERC, GIGABYTES
      ENDIF
 1000 FORMAT (" Node",I5," finished first half transformation",I16,
     &        " HT integrals written (",F6.2,"%, ",F9.2," GB)") 
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_REMOVE (I_FILE)
C
C     Write local integrals
C
C     Calling parameters
C
C     I_FILE   : File unit
C
C     Delete the buffer file (free some disk space)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
C
C     Get file unit number and file pointer
C
      LGFIL  = I_BF_UNIT+I_FILE
      CLOSE (LGFIL,STATUS='DELETE')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     The following routines are only needed in parallel runs, they take
C     care of the communication between the nodes
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (VAR_MPI)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_FREE (IFREE)
C
C     Check the status of buffers that were sent via a nonblocking send
C
C     Calling parameters
C
C     IFREE    : Contains on output the first free send buffer or zero
C                if none is free.
C
C     Written by Luuk Visscher, October 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "bucket_sort.h"
      LOGICAL READY
      integer :: status_array(df_mpi_status_size)
C
      IFREE = 0
      READY = .FALSE.
      DO I = 1, MAX_SEND_BUFFERS
         CALL interface_MPI_TEST(IREQ(I), READY, status_array)
         IF (READY) THEN
            IFREE = I
            RETURN
         ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_SEND (I_FILE,IBUF,RGBUF,IGBUF,N,INODE)
C
C     Send buffer IBUF to node INODE
C
C     Calling parameters
C
C     I_FILE   : File unit
C     IBUF     : Buffer number
C     RGBUF    : One buffer with real data
C     IGBUF    : One buffer with integer data
C     N        : Actual size of the buffer
C     INODE    : Node where the buffer should go
C
C     Written by Luuk Visscher, October 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
      REAL*8 RGBUF(NGBFSZ)
!     INTEGER*4 IGBUF(2,NGBFSZ) ! should be default int as far as i can see
      INTEGER   IGBUF(2,NGBFSZ)
C
C     Loop until we find a free send buffer IFREE
    1 CONTINUE
         CALL BUFFER_FILE_FREE(IFREE)
C        Check for incoming buffers
         CALL BUFFER_FILE_RECV()
         IF (IFREE.EQ.0) GOTO 1
C     Pack all data in the integer array ISGBUF
      ISGBUF(1,IFREE) = N
      ISGBUF(2,IFREE) = IBUF
      ISGBUF(3,IFREE) = I_FILE
      ISGBUF(4,IFREE) = IFREE
      CALL ICOPY (2*N,IGBUF,1,ISGBUF(5,IFREE),1)
      CALL DCOPY (N,RGBUF,1,ISGBUF(5+2*N,IFREE),1)
      call interface_mpi_issend_i2_work_f77(ISGBUF(1,IFREE),4*NGBFSZ+4,
     &               INODE,95,global_communicator,IREQ(IFREE))
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_RECV ()
C
C     Receive and write buffers that were sent to this node
C
C     Calling parameters
C
C     Written by Luuk Visscher, October 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
      INTEGER ISTAT(df_MPI_STATUS_SIZE)
      LOGICAL NEWMAIL
      INTEGER IHBUF(4*NGBFSZ+4)
C
C     Check whether there are any messages waiting
C
   1  call interface_mpi_iprobe(df_MPI_ANY_SOURCE,95,NEWMAIL,
     &                          global_communicator,istat)
      IF (NEWMAIL) THEN
C        Find out who was sending this
         INODE = ISTAT(df_mpi_source)
C        Receive the buffer
         call interface_mpi_RECV(IHBUF,4*NGBFSZ+4,INODE,95,
     &                 global_communicator)
C        Decode the information and determine where the buffer should go
         N      = IHBUF(1)
         IBUF   = IHBUF(2)
         I_FILE = IHBUF(3)
         IFREE  = IHBUF(4)
C        Write the buffers
         CALL BUFFER_FILE_WRITE (I_FILE,IBUF,IHBUF(5+2*N),IHBUF(5),N)
C        There may be more messages still waiting, go back to the start
         NEWMAIL = .FALSE.
         GOTO 1
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE BUFFER_FILE_END ()
C
C     Monitor status of our own pending messages and
C     stay in receiving mode until master tells us to stop
C
C     Calling parameters
C
C     Written by Luuk Visscher, October 2004
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
      INTEGER ISTAT(df_MPI_STATUS_SIZE)
      INTEGER ISTATS(df_MPI_STATUS_SIZE,MAX_SEND_BUFFERS)
      LOGICAL NEWMAIL
      INTEGER IEND
      INTEGER IBUFMPI(5)
      LOGICAL READY, ALLRECV
C
C     The actual value of IEND does not matter, we need to send
C     something. Likewise for IBUFMPI, only the third element is used
C
      CALL IZERO (IBUFMPI,5)
      IBUFMPI(3) = -1
      IEND = 1
C
C     The code for the slaves
C
      IF (MYTID.NE.MPARID) THEN
C
C     Start loop that will only end when master sends the termination
C     signal.
C
      ALLRECV = .FALSE.
   1  CONTINUE
C     Verify status of our own messages
      IF (.NOT.ALLRECV) THEN
         CALL interface_MPI_TESTALL(MAX_SEND_BUFFERS,IREQ, ALLRECV,
     &                              ISTATS)
C        Tell master in case all our messages were received
         IF (ALLRECV) CALL interface_MPI_ISEND(IEND,1,MPARID,98,
     &               global_communicator,IREQ(1))
      ENDIF
C     Check for end signal
      call interface_mpi_iprobe(mparid,99,NEWMAIL,global_communicator,
     &                          istat)
      IF (NEWMAIL) THEN
C        Receive the end signal message and return
         CALL interface_MPI_RECV(IEND,1,MPARID,99,
     &                           global_communicator)
         RETURN
      ELSE
C        There are still buffers waiting out there, see if we get some
         CALL BUFFER_FILE_RECV ()
         GOTO 1
      ENDIF
C
C     The code for the master
C
      ELSE
C
      NREADY = 0
    2 CONTINUE
C
C     Initially the slaves are still waiting for a new integral task,
C     tell them that they should do something else
C
      call interface_mpi_iprobe(df_mpi_any_source,20,NEWMAIL,
     &                          global_communicator,istat)
      IF (NEWMAIL) THEN
C        Find out who was sending this
         INODE = ISTAT(df_mpi_source)
C        Receive this signal
         CALL interface_MPI_RECV(IEND,1,INODE,20,
     &                           global_communicator)
         CALL interface_MPI_SEND(IBUFMPI,5,INODE,30,
     &                           global_communicator)
      ENDIF
C
C     Check for ready signals from the slaves
C
      call interface_mpi_iprobe(df_mpi_any_source,98,NEWMAIL,
     &                          global_communicator,istat)
      IF (NEWMAIL) THEN
C        Find out who was sending this
         INODE = ISTAT(df_mpi_source)
C        Receive this signal
         CALL interface_MPI_RECV(IEND,1,INODE,98,
     &                           global_communicator)
         NREADY = NREADY + 1
      ENDIF
C
C     Send out end signal when all nodes are ready
C
      IF (NREADY.EQ.NODES) THEN
         DO INODE = 1, NODES
            CALL interface_MPI_SEND(IEND,1,INODE,99,
     &                              global_communicator)
         ENDDO
         RETURN
      ELSE
C        There are still buffers waiting out there, see if we get some
         CALL BUFFER_FILE_RECV ()
         GOTO 2
      ENDIF
C
C     End of master-slave branching
C
      ENDIF
C
      RETURN
      END
#endif

c 
      SUBROUTINE BUFFER_FILE_READ (I_FILE,IBUF,A,NROW,ICOL,do_finesort)
C
C     Read back the sorted integrals and put them at the right position
C
C     Written by Luuk Visscher, August 2004
c
C        adapted by Andre Gomes,   March  2009
c
c        - difference now is that the records related to this buffer are sorted
c          and we can just read those related to the column we want. the information
c          of which columns are in which records is in the metadata file generated
c          by the fine sorting routines before the 2HT actually starts 
C
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
#include "infpar.h"
       REAL*8 RGBUF(NGBFSZ),A(NROW)
       DIMENSION IGBUF(2,NGBFSZ)
       logical do_finesort

       LGFIL         = I_BF_UNIT+I_FILE
       irec_initial  = LGREC(IBUF,I_FILE)
       irec_final    = 0

       if (do_finesort) then
          call set_record_range(IBUF,ICOL,irec_initial,irec_final)
       endif
 
       CALL DZERO (A,NROW)
 
       irec = irec_initial

    1  if ((irec.eq.0).or.(irec.lt.irec_final)) then
          return
       endif
C
C     Read a batch and the link to the previous record.
C
      READ (LGFIL,REC=IREC) IRECN,N,RGBUF,IGBUF
      IREC = IRECN
C
C     Pick out the integrals that we need.
C
      DO I = 1, N
         IF (IGBUF(2,I).EQ.ICOL) THEN
            A(IGBUF(1,I)) = RGBUF(I)
         else
         ENDIF
      ENDDO
C
C     We read in reverse order, if we have the first block we're done
C
      GOTO 1
C
      RETURN
      END


      subroutine set_record_range(IBUF,ICOL,irec_initial,irec_final)
#include "implicit.h"
#include "priunit.h"
#include "bucket_sort.h"
         integer, intent(in)  :: ibuf, icol
         integer, intent(out) :: irec_initial, irec_final 
         integer :: meta(3), imdrec, imdrecn

c getting metadata for this buffer, column combination in the metadata file

         LMFIL  = I_MD_UNIT+1
         imdrec = irec_md_lastwritten(ibuf)

         if (imdrec.eq.0.or.irec_initial.eq.0) then
            return
         endif

 11      if (imdrec.eq.0) goto 12
            read(LMFIL,REC=imdrec) imdrecn,iactive,kl,meta
            if (meta(1).eq.icol) goto 12

            imdrec = imdrecn
            goto 11
 12      continue

         irec_initial = meta(2)
         irec_final   = meta(3)

      end subroutine set_record_range

