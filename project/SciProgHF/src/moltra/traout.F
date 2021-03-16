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
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTI */
      SUBROUTINE MDINTI(IPRINT)
C
C     Initialize MDINT file
C     Written by Luuk Visscher, March 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 
      Use moltra_labeling

#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "infpar.h"
      CHARACTER*8 CDATE,TIMEX,DATEX*10
      CHARACTER*10 FNODE,LFNAM
C
C     Define a unique name for MDCINT. For backwards compatibility we name 
C     the file  on the master node simply MDCINT
C
      IF (MYTID.EQ.0) THEN
         FNODE = 'MDCINT    '
      ELSE
         FNODE = LFNAM('MDCINT')
      ENDIF
C
C     The date is kept in a 10-character variable in RELCCSD
C
      CALL GETDAT (CDATE,TIMEX)
      DATEX = ' '//CDATE//' '
C
      OPEN (LUMLF2,FILE=FNODE,FORM='UNFORMATTED')
      WRITE (LUMLF2) DATEX,TIMEX,
     & moltra_nkr,(moltra_isp(i),moltra_isp(-i),i=1,moltra_nkr)
      NTKL = 0
      N2EFIL = 0
C
      IF (TRA_ASCII) THEN
         IF (MYTID.EQ.0) THEN
            WRITE (LUASCII,*) "*** 2-electron integrals ***"
         ELSE
            CALL QUIT ('ASCII output not possible in parallel runs')
         END IF
      END IF
        
      IF (IPRINT.GE.5) WRITE (LUPRI,'(//3X,A)') 'Integrals written :'
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTF */
      SUBROUTINE MDINTF
C
C     Finish MDINT file
C     Written by Luuk Visscher, March 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dcbgen.h"
#include "infpar.h"
      WRITE (LUMLF2) 0,0,0
      CLOSE (LUMLF2,STATUS='KEEP')
      IF (TRA_ASCII) THEN
         WRITE (LUASCII,*) '*** End of File ***'
         CLOSE (LUASCII,STATUS='KEEP')
      ENDIF
C
      IF (MYTID.EQ.0) THEN
         WRITE (LUPRI,1000)
         IF (TRA_ASCII) WRITE (LUPRI,1001)
      END IF
      CALL FLSHFO(LUPRI)
 1000 FORMAT (' - Binary  file MDCINT was written.')
 1001 FORMAT (' - ASCII   file MO_integrals.txt was written.')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTN */
      SUBROUTINE MDINTN (NZKL)
C
C     Split MDINT file if necessary
C     Written by Luuk Visscher, Spring 2000
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
      CHARACTER*5 NEXTFIL
      CHARACTER*10 FNODE,LFNAM
C
      NTKL = NTKL + NZKL
C     Every integral takes at most 32 bytes on 4 byte integer machines
C     Make a new file if we have written more than 60 million integrals.
C     This solves the 2Gb problem on Linux.
      IF (NTKL.GT.60 000 000) THEN
         N2EFIL = N2EFIL + 1
         IF (N2EFIL.LT.10) THEN
            WRITE (NEXTFIL,'(A4,I1)') 'MDC0',N2EFIL
         ELSEIF (N2EFIL.LT.100) THEN
            WRITE (NEXTFIL,'(A3,I2)') 'MDC',N2EFIL
         ELSE
            CALL QUIT('*** ERROR in MDINTN: exceeding 100 integral'//
     &                ' files! ***')
         ENDIF
         FNODE = LFNAM(NEXTFIL)
         WRITE (LUMLF2) 0,N2EFIL,0
         WRITE (LUMLF2) FNODE
         CLOSE (LUMLF2,STATUS='KEEP')
         OPEN (LUMLF2,FILE=FNODE,FORM='UNFORMATTED')
         NTKL = 0
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDINTB */
      SUBROUTINE MDINTB(IPRINT,RCW,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &                  FMAT,NFM,INDK,INDL,RKL)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     Write block of integrals to MDINT file. This is the old MOLFDIR
C     format and should work with both RELCCSD and DIRRCI.
C
C     26-6-2002 : Slightly modified to avoid writing of imaginary zeroes
C     in case of real groups.
C
C     Called from SYMFINT.
C
C     Written by Luuk Visscher, March 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

      Use moltra_labeling

#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      INTEGER RCW
      DIMENSION NSTR3(2),NSTR4(2)
      DIMENSION FMAT(NFM,RCW),INDK(*),INDL(*)
      REAL*8 RKL(*)

      CALL QENTER('MDINTB')

      IF (IPRINT.GE.5) THEN
      CALL HEADER('Output from MDINTB routine',-1)
      SELECT CASE (ICLASS)
       CASE(1); write(LUPRI,'(A)') 'Class 1 : (i    j    | k    l   )'
       CASE(2); write(LUPRI,'(A)') 'Class 2 : (i    j    | kbar lbar)'
       CASE(3); write(LUPRI,'(A)') 'Class 3 : (i    jbar | k    lbar)'
       CASE(4); write(LUPRI,'(A)') 'Class 4 : (i    jbar | kbar l   )'
       CASE(9); write(LUPRI,'(A)') 'Class 9 : (i    j    | kbar l   )'
       CASE(10); write(LUPRI,'(A)') 'Class 10 : (i    j    | k    lbar)'
       CASE(11); write(LUPRI,'(A)') 'Class 11 : (ibar j    | k    l   )'
       CASE(12); write(LUPRI,'(A)') 'Class 12 : (i    jbar | k    l   )'
       CASE DEFAULT
         write(LUPRI,'(A,I5)')'WARNING Unknown integral class',ICLASS
      END SELECT
C     ...      
      ENDIF
C
      NZKL = 0
      IKR1 = IKR
      JKR1 = JKR
      KL = 0
      LL = 0
      DO IREPL = 1, NFSYM
         IREPK = MOD(IREPL+IREPIJ,2) + 1

         IF (IPRINT.GE.5) THEN
          WRITE(LUPRI,'(2X,A,I3,A,I3)') 'IREPK=',IREPK,' IREPL=',IREPL
          WRITE(LUPRI,'(2X,A,I3,A,I3)') ' NSTR3(IREPK)=',NSTR3(IREPK), 
     &                                  ' NSTR4(IREPL)=',NSTR4(IREPL)
         ENDIF
         
         DO L = 1, NSTR4(IREPL)
            LL = LL + 1
            DO K = 1, NSTR3(IREPK)
               IF (IREPK.EQ.1) THEN
                  KK = K
               ELSE
                  KK = K + NSTR3(1)
               ENDIF
               KL = KL + 1
               IF (RCW.EQ.1) THEN
                  IF (DABS(FMAT(KL,1)).GT.THROUT) THEN
                     NZKL = NZKL + 1
                     INDK(NZKL) = KK
                     INDL(NZKL) = LL
                     RKL(NZKL) = FMAT(KL,1)
                  ENDIF
               ELSE
                  IF (DABS(FMAT(KL,1)).GT.THROUT.OR.
     &                DABS(FMAT(KL,2)).GT.THROUT) THEN
                     NZKL = NZKL + 1
                     INDK(NZKL) = KK
                     INDL(NZKL) = LL
                     RKL(2*NZKL-1) = FMAT(KL,1)
                     RKL(2*NZKL)   = FMAT(KL,2)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      IF (ICLASS.EQ.1) THEN
         CONTINUE
      ELSEIF (ICLASS.EQ.2) THEN
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.3) THEN
         JKR1 = - JKR1
         DO INZ = 1, NZKL
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.4) THEN
         JKR1 = - JKR1
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.9) THEN
         DO INZ = 1, NZKL
            INDK(INZ) = - INDK(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.10) THEN
         DO INZ = 1, NZKL
            INDL(INZ) = - INDL(INZ)
         ENDDO
      ELSEIF (ICLASS.EQ.11) THEN
         IKR1 = - IKR1
      ELSEIF (ICLASS.EQ.12) THEN
         JKR1 = - JKR1
      ELSE
         CALL QUIT ('Class not defined in MDINTB')
      ENDIF
C
C     ... write the integral batch to the file 
      WRITE (LUMLF2) IKR1,JKR1,NZKL,
     &               (INDK(INZ),INDL(INZ),INZ=1,NZKL),
     &               (RKL(INZ),INZ=1,NZKL*RCW)

C     Write out integrals in ASCII format

      IF (TRA_ASCII) THEN
        i = Kramer_to_SpinorIndex(ikr1)
        j = Kramer_to_SpinorIndex(jkr1)
        do inz = 1, nzkl
          k = Kramer_to_SpinorIndex(indk(inz))
          l = Kramer_to_SpinorIndex(indl(inz))
          IF (NZ.EQ.1) THEN
          write (LUASCII,'(E30.16,4I6)') 
     &          (RKL(RCW*(INZ-1)+IW),IW=1,RCW),
     &          i,j,k,l
          ELSE
          write (LUASCII,'(2E30.16,4I6)') 
     &          (RKL(RCW*(INZ-1)+IW),IW=1,RCW),
     &          i,j,k,l
          ENDIF
        end do
      END IF
C
C     Write out integrals to output if desired. 
C   ------------------------------------------------------------------
      IF (IPRINT.GE.25) THEN
         WRITE(LUPRI,'(2x,a,3i4)')
     &  'MOLTRA: integrals batch saved;  IKR1,JKR1,NZKL:',
     &   IKR1,JKR1,NZKL
         DO INZ = 1, NZKL
           WRITE (LUPRI,'(I4,A,I4,A,A,I4,A,I4,A,2f30.16)')
     &     INZ,'./',NZKL,' >',
     &     ' INDK=',INDK(INZ),' INDL=',INDL(INZ),
     &     ' intgrl:',(RKL(RCW*(INZ-1)+I),I=1,RCW)
         ENDDO
      ENDIF
C
      CALL QEXIT('MDINTB')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck LINZERO */
      SUBROUTINE LINZERO(ICLASS,IREPIJB,IREPKLB,NKL,FMAT)
C
C     Zero integrals that should be zero due to symmetry
C     This means that the MJ values stored in IREPIJB and IREPKLB should match.
C
C     Written by Luuk Visscher, November 2003
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
      DIMENSION FMAT(NKL)
      DIMENSION IREPIJB(2),IREPKLB(2,NKL)
C
C     To determine the MJ value we need to know the classes :
C     Class 1 : (i    j    | k    l   )
C     Class 2 : (i    j    | kbar lbar)
C     Class 3 : (i    jbar | k    lbar)
C     Class 4 : (ibar j    | k    lbar)
C
      IF (ICLASS.LE.2) THEN
         MIJ = - IREPIJB(1) + IREPIJB(2)
      ELSEIF (ICLASS.EQ.3) THEN
         MIJ = - IREPIJB(1) - IREPIJB(2)
      ELSE
         MIJ =   IREPIJB(1) + IREPIJB(2)
      ENDIF
C
      DO I = 1, NKL
         IF (ICLASS.EQ.1) THEN
            MKL = - IREPKLB(1,I) + IREPKLB(2,I)
         ELSEIF (ICLASS.EQ.2) THEN
            MKL =   IREPKLB(1,I) - IREPKLB(2,I)
         ELSE
            MKL = - IREPKLB(1,I) - IREPKLB(2,I)
         ENDIF
!        write(lupri,'(2x,a,2i6)') '(linzero): MIJ, MKL are:', MIJ,MKL
!        IF (MIJ+MKL/=0)
!    &   write(lupri,'(2x,a,i6,a,1p,d10.2)')
!    &   '(linzero): removed integral(',iclass,')',FMAT(I)
!        IF (MIJ+MKL==0)
!    &   write(lupri,'(2x,a,i6,a,1p,d10.2)')
!    &   '(linzero): kept integral(',iclass,')',FMAT(I)
         IF (MIJ+MKL.NE.0) FMAT(I) = 0.D0
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SYMFINT */
      SUBROUTINE SYMFINT(IPRINT,IREPIJ,IREPIJB,IREPKLB,IJ,IKR,JKR,
     &           NSTR3,NSTR4,NKL,FMAT)
C
C     Symmetrize the final integrals and write them to file
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dcbham.h"
#include "dgroup.h"
      INTEGER NSTR3(2),NSTR4(2)
      DIMENSION FMAT(NKL,NZ,NZ,*)
      DIMENSION IREPIJB(2),IREPKLB(2,NKL)
      INTEGER RCW
      Real(8), Allocatable :: FTMAT(:,:),RKL(:,:)
      Integer, Allocatable :: INDK(:),INDL(:)
C
C     Set parameter to indicate whether the integrals are complex or real.
C
      IF (SPINFR.OR.NZ.EQ.1) THEN
         RCW = 1
      ELSE
         RCW = 2
      ENDIF
C
      allocate(FTMAT(NKL,RCW))
      allocate(  RKL(NKL,RCW))
      allocate(INDK(NKL))
      allocate(INDL(NKL))
C
      IF (SPINFR) THEN
C
C     Special code that takes out the quaternion phase factors
C     The resulting integrals have full 8-fold permutational
C     symmetry and can be fed into non-relativistic codes.
C     Only one class of integrals (all-unbarred) suffices for
C     this purpose.
C
      CALL DZERO(FTMAT,RCW*NKL)
      IREPABR = 0
      DO IREPAB = 0, NBSYM-1
         IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
            IREPABR = IREPABR + 1
C
C           Zero out the spin-orbit coupling for electron 1
C
            IF (IREPAB.NE.IREPIJB(1)) GOTO 100
            DO I = 1, NKL
               IF (IREPAB.NE.IREPKLB(1,I)) THEN
C
C                  Zero out the spin-orbit coupling for electron 2
C
                   FMAT(I,1,1,IREPABR) = 0.D0
               ELSE
C
C                  Take out quaternion phase factor by multiplying
C                  with the compound phase factor that was introduced
C                  in the consecutive multiplications
C 
                   FMAT(I,1,1,IREPABR) = 
     &             FMAT(I,1,1,IREPABR) *
     &             IREPIJB(2)*IREPKLB(2,I)
                   
               ENDIF
            ENDDO
            FACTR = 1.D0
            CALL DAXPY (NKL,FACTR,
     &                  FMAT(1,1,1,IREPABR),1,
     &                  FTMAT,1)
 100        CONTINUE
         ENDIF
      ENDDO
C
C     Only one class remains in the spinfree case
C
      ICLASS = 1
      CALL MDINTB(IPRINT,RCW,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &            FTMAT,NKL,INDK,INDL,RKL)
C
      ELSE ! not spinfree
C
C     Real and complex classes
C
      DO ICLASS = 1, 4
         CALL DZERO(FTMAT,RCW*NKL)
         IREPABR = 0
         DO IREPAB = 0, NBSYM-1
            IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
               IREPABR = IREPABR + 1
               DO IZ1 = 1, NZ
                  IQ1 = IPQTOQ(IZ1,IREPAB)
                  DO IZ2 = 1, NZ
                     IQ2 = IPQTOQ(IZ2,IREPAB)
                     IFAC = INTCL(ICLASS,IQ1,IQ2)
                     IF (IFAC.NE.0) THEN
                        FACTR = IFAC
                        CALL DAXPY (NKL,FACTR,
     &                              FMAT(1,IZ2,IZ1,IREPABR),1,
     &                              FTMAT(:,1),1)
                     ENDIF
                     IF (NZ.GT.1) THEN
                        IFAC=INTCL(ICLASS+4,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                           FACTR = IFAC
                           CALL DAXPY (NKL,FACTR,
     &                                 FMAT(1,IZ2,IZ1,IREPABR),1,
     &                                 FTMAT(:,2),1)
                        ENDIF
                     ENDIF
                  ENDDO
               ENDDO
            ENDIF
         ENDDO
C
C        Take out integrals that should be zero due to the linear symmetry (can be non-zero due to noise)
C
         IF (LINEAR)
     &   CALL LINZERO (ICLASS,IREPIJB,IREPKLB,NKL,FTMAT)

         CALL MDINTB(IPRINT,RCW,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &               FTMAT,NKL,INDK,INDL,RKL)
      ENDDO
C
C     Quaternion classes
C
      IF (NZ.EQ.4) THEN
         DO ICLASS = 9, 12
            CALL DZERO(FTMAT,RCW*NKL)
            IREPABR = 0
            DO IREPAB = 0, NBSYM-1
               IF (IREPIJ.EQ.JBTOF(IREPAB,1)) THEN
                  IREPABR = IREPABR + 1
                  DO IZ1 = 1, NZ
                     IQ1 = IPQTOQ(IZ1,IREPAB)
                     DO IZ2 = 1, NZ
                        IQ2 = IPQTOQ(IZ2,IREPAB)
                        IFAC = INTCL(ICLASS,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                            FACTR = IFAC
                            CALL DAXPY (NKL,FACTR,
     &                                  FMAT(1,IZ2,IZ1,IREPABR),1,
     &                                  FTMAT(:,1),1)
                        ENDIF
                        IFAC=INTCL(ICLASS+4,IQ1,IQ2)
                        IF (IFAC.NE.0) THEN
                           FACTR = IFAC
                           CALL DAXPY (NKL,FACTR,
     &                                 FMAT(1,IZ2,IZ1,IREPABR),1,
     &                                 FTMAT(:,2),1)
                        ENDIF
                     ENDDO
                  ENDDO
               ENDIF
            ENDDO
            CALL MDINTB(IPRINT,RCW,ICLASS,IREPIJ,NSTR3,NSTR4,IKR,JKR,
     &                  FTMAT,NKL,INDK,INDL,RKL)
         ENDDO
      ENDIF
C
      ENDIF ! spinfree
C
      deallocate(FTMAT)
      deallocate(RKL)
      deallocate(INDK)
      deallocate(INDL)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MKINDXKR*/
      SUBROUTINE MKINDXKR (NSTR,TRIAN,INDXKR)
C
C     Make an array that gives the Kramers pairs from the ij-value
C
C     Written by Luuk Visscher, April 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      LOGICAL TRIAN
      INTEGER NSTR(2,0:2,2)
      DIMENSION INDXKR(2,*)
C
      IJ = 0
      DO IREPIJ =1, NFSYM
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
            IF (.NOT.TRIAN.OR.(TRIAN.AND.IREPI.GT.IREPJ)) THEN
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = 1, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    INDXKR(1,IJ) = I
                    INDXKR(2,IJ) = J
                    IF (IREPI.NE.1) INDXKR(1,IJ) = I + NSTR(IREPI-1,0,1)
                    IF (IREPJ.NE.1) INDXKR(2,IJ) = J + NSTR(IREPJ-1,0,1)
                 ENDDO
              ENDDO
            ELSEIF (IREPI.EQ.IREPJ) THEN
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = J, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    INDXKR(1,IJ) = I
                    INDXKR(2,IJ) = J
                    IF (IREPI.NE.1) INDXKR(1,IJ) = I + NSTR(IREPI-1,0,1)
                    IF (IREPJ.NE.1) INDXKR(2,IJ) = J + NSTR(IREPJ-1,0,1)
                 ENDDO
              ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRI */
      SUBROUTINE MDSCRI
C
C     Initialize MDINT scratch files
C     Written by Luuk Visscher, May 1997
C     Modified by HJAaJ Aug 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      LOGICAL   FEXIST
      CHARACTER*10 FNODE,LFNAM,FNAM*5
      Real(8), Allocatable :: FTMAT(:)
C
C     Calculate size of integral batches that will be written
      NQQCLASS = NZ * NZ * NBSYM / NFSYM        ! Number of quaternion units active for each electron
      NKL = NFPCK34(1)                          ! Dimension of MO matrix (gg+uu)
      IF (NFSYM.EQ.2) NKL = MAX(NKL,NFPCK34(2)) ! Dimension of MO matrix (ug+gu)
      LFMAT = NQQCLASS*NKL                      ! Total maximum dimension of bacth

      call alloc(FTMAT,LFMAT)
      CALL DZERO(FTMAT,LFMAT)
C
C     Open direct access scratch file for 4-indexed transformed integrals
C
      DO IREPIJ = 1, NFSYM
         LFMAT = NQQCLASS*NFPCK34(IREPIJ)
         IRECL = IRECLEN(LFMAT,0,0)
         IRECL = MAX(IRECL,IRECLEN(1,0,0))
         WRITE (FNAM,'(A4,I1)') '4IND',IREPIJ
         FNODE = LFNAM(FNAM)

C        Open and delete first if old file, for some systems
C        leaving out this leads to errors in MCSCF where the record
C        length changes as the calculation continues
         INQUIRE(FILE=FNODE,EXIST=FEXIST)
         IF (FEXIST) THEN
            JRECL = IRECLEN(1,0,0) 
            OPEN (LUMLF3+IREPIJ,FILE=FNODE,STATUS='OLD',
     &            ACCESS='DIRECT',RECL=JRECL)
            CLOSE (LUMLF3+IREPIJ,STATUS='DELETE')
         END IF

C        Now open the new file and initialize the records with zeroes
         OPEN (LUMLF3+IREPIJ,FILE=FNODE,ACCESS='DIRECT',RECL=IRECL)
         NIJ = NFPCK12(IREPIJ)
         CALL MDSCRI1(LUMLF3+IREPIJ,NIJ,LFMAT,FTMAT)
      ENDDO
C
      call dealloc(FTMAT)
      END

      SUBROUTINE MDSCRI1(LU,NIJ,LFMAT,FMAT)
C
C     HJAaJ Aug 2001
C
#include "implicit.h"
      DIMENSION FMAT(LFMAT)
      DO IJ = 1, NIJ
         WRITE (LU,REC=IJ) FMAT
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRR */
      SUBROUTINE MDSCRR(IREPIJ,IJ,LFMAT,FMAT)
C
C     Read the final integrals from the scratch file
C
C     Written by Luuk Visscher, July 1997
C     Modified by HJAaJ Aug 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbtra.h"
      DIMENSION FMAT(LFMAT)
C
C     Read integrals.
C
      READ (LUMLF3+IREPIJ,REC=IJ) FMAT
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRU */
      SUBROUTINE MDSCRU(IREPIJ,IJ,LFMAT,FMAT,FTMP)
C
C     Update the scratch file with final integrals
C
C     Written by Luuk Visscher, July 1997
C     Revised HJAaJ Aug 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbtra.h"
      DIMENSION FMAT(LFMAT), FTMP(LFMAT)
C
C     Read integrals from previous batch and add up.
C
      READ (LUMLF3+IREPIJ,REC=IJ) FTMP
      DO I = 1,LFMAT
         FTMP(I) = FTMP(I) + FMAT(I)
      END DO
      WRITE (LUMLF3+IREPIJ,REC=IJ) FTMP
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRW */
      SUBROUTINE MDSCRW(IREPIJ,IJ,LFMAT,FMAT)
C
C     Write the final integrals to scratch file.
C
C     Written by Luuk Visscher, July 1997
C     Modified HJAaJ Aug 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbtra.h"
      DIMENSION FMAT(LFMAT)
C
      WRITE (LUMLF3+IREPIJ,REC=IJ) FMAT
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MDSCRF */
      SUBROUTINE MDSCRF
C
C     Delete MDSCR file
C     Written by Luuk Visscher, May 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbtra.h"
#include "dgroup.h"
C
      DO IREPIJ = 1, NFSYM
         IF (.NOT. NOMDCINT) THEN
C           delete the files
            CLOSE (LUMLF3+IREPIJ,STATUS='DELETE')
         ELSE
C           close the files
            CLOSE (LUMLF3+IREPIJ,STATUS='KEEP')
         END IF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MKINDXB*/
      SUBROUTINE MKINDXB (NSTR,IBE11,IBE21,IBE12,IBE22,TRIAN,INDXB)
C
C    Make an array that gives the compound boson irrep from the ij-value
C    It also gives the sign of the quaternion phase factor that this
C    density obtains in the multiplication of the MO coefficients. This
C    sign is required when we want to backtransform (in spinfree calcs)
C    the integrals such that only class 1 integrals with full 8-fold
C    permutation symmetry remain.
C
C    Written by Luuk Visscher, April 1997
C    Last modification, Februari 2001
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dcbtra.h"
#include "dgroup.h"
      LOGICAL TRIAN
      INTEGER NSTR(2,0:2,2)
      DIMENSION IBE11(*),IBE21(*),IBE12(*),IBE22(*)
      DIMENSION INDXB(2,*)
#include "ibtfun.h"
C
      IJ = 0
      DO IREPIJ =1, NFSYM
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
            IF (.NOT.TRIAN.OR.(TRIAN.AND.IREPI.GT.IREPJ)) THEN
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = 1, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    IF (IREPI.EQ.1) IREPIB = IBE11(I)
                    IF (IREPI.EQ.2) IREPIB = IBE21(I)
                    IF (IREPJ.EQ.1) IREPJB = IBE12(J)
                    IF (IREPJ.EQ.2) IREPJB = IBE22(J)
C                   The extra label may be due to linearity or due to the spinfree option
                    IF (LINEAR) THEN
C
C                      Store the MJ values in the INDXB array
C
                       INDXB(1,IJ) = IREPIB
                       INDXB(2,IJ) = IREPJB
                    ELSE
C
C                      The compound boson irrep can be obtained
C                      by an xor operation.
C
                       IREPIJB = IBTXOR(IREPIB,IREPJB)
                       INDXB(1,IJ) = IREPIJB
C
C                      We first determine which quaternion phase
C                      factors the individual functions receive and
C                      than obtain the compound phase factor by
C                      quaternion multiplication.
C
                       IQI = IPQTOQ(1,IREPIB)
                       IQJ = IPQTOQ(1,IREPJB)
C
C                      The sign of the quaternion multiplication is given by iqphase,
C                      but since IQJ is conjugated (backtransform !) we also have to 
C                      multiply by IQSIGN (-1 for q-imaginary, 1 for real)
C
                       INDXB(2,IJ) = IQPHASE(IQI,IQJ,1)*IQSIGN(IQJ,2,1)
                    ENDIF
                 ENDDO
              ENDDO
            ELSEIF (IREPI.EQ.IREPJ) THEN
C
C             Same procedure for lower trangular gg or uu case..
C
              DO J = 1, NSTR(IREPJ,0,2)
                 DO I = J, NSTR(IREPI,0,1)
                    IJ = IJ + 1
                    IF (IREPI.EQ.1) IREPIB = IBE11(I)
                    IF (IREPI.EQ.2) IREPIB = IBE21(I)
                    IF (IREPJ.EQ.1) IREPJB = IBE12(J)
                    IF (IREPJ.EQ.2) IREPJB = IBE22(J)
                    IF (LINEAR) THEN
                       INDXB(1,IJ) = IREPIB
                       INDXB(2,IJ) = IREPJB
                    ELSE
                       IREPIJB = IBTXOR(IREPIB,IREPJB)
                       INDXB(1,IJ) = IREPIJB
                       IQI = IPQTOQ(1,IREPIB)
                       IQJ = IPQTOQ(1,IREPJB)
                       INDXB(2,IJ) = IQPHASE(IQI,IQJ,1)*IQSIGN(IQJ,2,1)
                    ENDIF
                 ENDDO
              ENDDO
            ENDIF
         ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck INI4INDINFO */
      SUBROUTINE INI4INDINFO
C
C     Initialize 4INDINFO file
C     Authored by Joern Thyssen, 1/11/2000.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "dcbtra.h"
#include "dgroup.h"
#include "priunit.h"
C
      CHARACTER*8 DATEX,TIMEX
C
C     Open file
C
      OPEN (LUMLF2,FILE='4INDINFO',FORM='UNFORMATTED')
C
C     Write date and time
C
      CALL GETDAT (DATEX,TIMEX)
      WRITE (LUMLF2) DATEX,TIMEX
C
C     The rest of the file is written in various places. 
C     The contents is:
C
C     Number of IJ indices for each fermion symmetry
C     Number of KL indices for each fermion symmetry
C     I and J indices
C     Boson irrep for IJ indices
C     Boson irrep for KL indices
C
C
      RETURN
      END
