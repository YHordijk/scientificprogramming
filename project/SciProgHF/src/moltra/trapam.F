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

! define task symbols for CALL DIRAC_PARCTL( task )
#include "dirac_partask.h"
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PAMTRA */
      SUBROUTINE PAMTRA()
C***********************************************************************
C
C     Driver routine for 4-index transformation
C
C     Written by T.Saue, J.Laerdahl and L.Visscher Jan 1997
C
C     Called from main/dirac.F if DOTRA=.true.
C
C***********************************************************************

      Use memory_allocator
      Use moltra_labeling
#ifdef MOD_EXACORR
      Use exacorr_ao_to_mo
#endif

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "maxorb.h"
#include "dcbpsi.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"

! until the KRMCSCF/KRMCOLD issue with IBEIG has been fixed...
#include "dcbham.h"

#include "dcbmp2.h"
      CHARACTER*24 DAYTID
      LOGICAL TOBE,TOBEK,TOBEKK,TRASAM,ACTIVE
      DIMENSION NSTR(2,0:2,4),KVEC(2,5),KQ(2,4),KIBE(2,4),
     &          IR(2,2),NR(2),NQ(2,4),NQT(4),KE(2,4),
     &          NSPC(2,0:2),KQC(2),NSTRT(4),NSPC2(2,0:3)
      Real(8), Allocatable :: CMO(:),EIG(:)
      Integer, Allocatable :: IBEIG(:)
      real(8), allocatable :: WORK(:)
C     Initialize (for debug mode)
      KE=0
      KIBE=0
      KQC=0
C
      CALL QENTER('PAMTRA')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in PAMTRA')
      KFRSAV = KFREE

C
      CALL GETTIM(CPU1,WALL1)
      CALL TRAHI(IPRTRA,3)

C     Define NOMDCINT from **MOLTRA input

      NOMDCINT = NOMDCINT_TRAINP
      ISTRAT_save = ISTRAT
      IF (NOMDCINT .AND. ISTRAT .NE. 4) THEN
         WRITE(LUPRI,'(/A,I0/)')
     &   'INFO, .NOMDCINT is only implemented for strategy 4;'//
     &   ' therefore use stategy 4 instead of ',ISTRAT
         ISTRAT = 4
      END IF

C
C     Define number of double quaternionic classes
C
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Set up index arrays for active orbitals
C     =======================================
C
      IF(NOPAIR) THEN
        DO I = 1,NFSYM
          NR(I)   = NESH(I)
          IR(1,I) = 1
          IR(2,I) = NESH(I)
        ENDDO
      ELSE
        DO I = 1,NFSYM
          NR(I)   =  NORB(I)
          IR(1,I) = -NPSH(I)
          IR(2,I) =  NESH(I)
        ENDDO
      ENDIF
C
C     NSTR(I): Number of Spinors to be TRansformed for index I
C
      DO I = 1,4
        NSTRT(I) = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET2('INTE','VEC tmp',KVEC(IFRP,I),NR(IFRP),
     &       WORK,KFREE,LFREE)
          NSTR(IFRP,0,I) = - 1
          CALL  NUMLST(TRA4_INDSTR(I,IFRP),WORK(KVEC(IFRP,I)),
     &                 NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                 IFRP,NSTR(IFRP,0,I))
          CALL ORBCNT(WORK(KVEC(IFRP,I)),NSTR(IFRP,0,I),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          NSTR(IFRP,0,I) = NSTR(IFRP,1,I) + NSTR(IFRP,2,I)
          NSTRT(I) = NSTRT(I) +  NSTR(IFRP,0,I)
          CALL MEMREL('PAMTRA.1',WORK,KFRSAV,KVEC(IFRP,I),KFREE,LFREE)
          CALL MEMGET2('INTE','VEC',KVEC(IFRP,I),NSTR(IFRP,0,I),
     &                WORK,KFREE,LFREE)
        ENDDO
      ENDDO
C
 
      CALL Make_Kramer_to_SpinorIndex (NFSYM,NSTR)

      KFRSAV = KFREE
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(/A)') '* Orbital ranges for 4-index transformation:'
      DO I = 1, NFSYM
         CALL TRAPRI(4,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),
     &                 WORK(KVEC(I,3)),WORK(KVEC(I,4)),NSTR)
      ENDDO

C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 4 arrays of two matrices each
C     They may share the same memory if the transformation ranges
C     are identical
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1, 4
         NQT(I) = 0
         ICMOQR(1,I) = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP,I) = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NQT(I) = NQT(I) + NQ(IFRP,I)
            IF (IFRP.LT.NFSYM)
     &         ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP,I)
         ENDDO
      ENDDO
C
C     Allocate the memory for the coefficients, check their relations
C
      K = 2
      IF (NFSYM.EQ.1) K=1
      DO I = 1, 4
C
         DO J = I, 1, -1
            IF (TRASAM(WORK(KVEC(1,I)),WORK(KVEC(K,I)),
     &                 WORK(KVEC(1,J)),WORK(KVEC(K,J)),
     &          NSTR(1,0,I),NSTR(1,0,J))) ISAME(I) = J
         ENDDO
C
         IF (ISAME(I).EQ.I) THEN
            IF (I.NE.1) CALL QUIT ('All ranges should be equal')
            CALL MEMGET2('REAL','Qmat',KQ(1,I),NQT(I),WORK,KFREE,LFREE)
         ELSE
            KQ(1,I) = KQ(1,ISAME(I))
         ENDIF
         KQ(2,I) = KQ(1,I) + NQ(1,I)
C
      ENDDO
C
C     Selected eigenvalues
C
      DO I = 1,4
        NDIM = 0
        DO IFRP = 1, NFSYM
          NDIM = NDIM + NSTR(IFRP,0,I)
        ENDDO
        CALL MEMGET2('REAL','EIGVAL',KE(1,I),NDIM,WORK,KFREE,LFREE)
        KE(2,I) = KE(1,I) + NSTR(1,0,I)
      ENDDO
C
C     Allocate memory for the integer array with information about the
C     spinors.
C
      DO I = 1, 4
       DO IFRP = 1, NFSYM
        CALL MEMGET2('INTE','IBEIG',KIBE(IFRP,I),NSTR(IFRP,0,I),
     &     WORK,KFREE,LFREE)
       ENDDO
      ENDDO
C
C     Get all coefficients
C
      Allocate (CMO(NCMOTQ))
      Allocate (EIG(NORBT))
      Allocate (IBEIG(NORBT))
      CALL IZERO(IBEIG,NORBT)
      CALL DZERO(EIG,NORBT)

      SELECT CASE ( MOFILE_TRAINP )
         CASE ( "DFCOEF" )
            TOBEK  = .FALSE.
            TOBEKK = .FALSE.
         CASE ( "KRMCSCF" )
            TOBEK  = .TRUE.
            TOBEKK = .FALSE.
         CASE ( "KRMCOLD" )
            TOBEK  = .FALSE.
            TOBEKK = .TRUE.
         CASE ( "UNKNOWN" )
            INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
            INQUIRE(FILE='KRMCOLD',EXIST=TOBEKK)
         CASE DEFAULT
            WRITE(LUPRI,'(/2A)')
     &      'FATAL ERROR for **MOLTRA input: Unknown .MOFILE name : ',
     &      MOFILE_TRAINP
            CALL QUIT('Unknown .MOFILE from **MOLTRA input')
      END SELECT


      IF (TOBEK) THEN
         CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','PAMTRA')
         JRDMO = -1
         CALL RREADMO(CMO,JRDMO,1,LUKRMC)
         IF (JRDMO .EQ. 0) THEN
            IF (IPRTRA .GE. 0) THEN
               WRITE(LUPRI,'(/A)')
     &              ' (PAMTRA)  Orbitals read from' //
     &              ' label NEWORB on file KRMCSCF'
            END IF
         END IF
         ! stefan - may 2013: the above fix (zeroing) is dirty - fix
         ! it soon!
         WRITE(LUPRI,'(/A)')
     &   ' (PAMTRA)  WARNING!!! zero   ' //
     &   ' boson information and eigenvalues. FIXME!'
         CLOSE (LUKRMC)
         IF(spinfr.or.linear) 
     &   call quit('pamtra: fix spinfr/linear cases with KRMCSCF')
      ELSE IF (TOBEKK) THEN
         CALL OPNFIL(LUKRMC,'KRMCOLD','OLD','PAMTRA')
         JRDMO = -1
         CALL RREADMO(CMO,JRDMO,3,LUKRMC)
         IF (JRDMO .EQ. 0) THEN
           WRITE(LUPRI,'(/A)')
     &     ' (PAMTRA)  Orbitals read from' //
     &     ' label NEWNATOB on file KRMCOLD'
         END IF
         ! stefan - may 2013: the above fix (zeroing) is dirty - fix
         ! it soon!
         WRITE(LUPRI,'(/A)')
     &   ' (PAMTRA)  WARNING!!! zero   ' //
     &   ' boson information and eigenvalues. FIXME!'
         IF(spinfr.or.linear) 
     &   call quit('pamtra: fix spinfr/linear cases with KRMCOLD')
         
         CLOSE (LUKRMC)
      ELSE
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF (.NOT.TOBE) THEN
            WRITE(LUPRI,'(/A)')
     &        'PAMTRA ERROR: MO coefficient file not present'
            CALL QUIT('PAMTRA: No MO coefficient file !')
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMTRA')
            CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)
            CLOSE(LUCOEF,STATUS='KEEP')
         END IF
         IF (IPRTRA .GE. 0) THEN
            WRITE(LUPRI,'(/A)')
     &          ' (PAMTRA)  Orbitals read from DFCOEF'
         END IF
      END IF
C
C     Select the set of MOs(CMO), eigenvalues(EIG) and boson array id(IBEIG) that we need for core transformation
C
      DO I = 1, 4
        DO IFRP = 1, NFSYM
          IF(NSTR(IFRP,0,I).GT.0) THEN
             CALL SELCFS (CMO(1+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,I)),
     &                   NSTR(IFRP,0,I),WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                   NFBAS(IFRP,0),NORB(IFRP))
             CALL SELEIG(EIG(1+IORB(IFRP)),IFRP,WORK(KE(IFRP,I)),
     &                 WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I))
             CALL SELIBEIG(IBEIG,IORB(IFRP),IFRP,
     &                     WORK(KIBE(IFRP,I)),WORK(KVEC(IFRP,I)),
     &                     NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          ENDIF
        ENDDO
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
C
      DeAllocate (CMO)
      DeAllocate (EIG)
      DeAllocate (IBEIG)
C
C     If necessary recanonize (virtual) orbitals
C     Assuming that all sets are equal, we use the first index !!!
C     ============================================================
C
      IF (RCORBS) THEN
C        Overdimension CMO because it is used as Focksize scratch space
         Allocate (CMO(N2BBASXQ))
         CALL VIRCAN(WORK,KFREE,LFREE,IPRTRA,KQ(1,1),KE(1,1),
     &               KIBE(1,1),NDMOQR(1,1,1),ICMOQR(1,1),
     &               CMO)
         DeAllocate (CMO)
      ENDIF
C
C     Set up index arrays for core orbitals
C     =====================================
C
      IC   = 1
      DO I = 1,NFSYM
        NR(I)   = NESH(I)
        IR(1,I) = 1
        IR(2,I) = NESH(I)
      ENDDO
C
      DO IFRP = 1,NFSYM
         CALL MEMGET2('INTE','VEC tmp',KVEC(IFRP,5),NR(IFRP),
     &      WORK,KFREE,LFREE)
         NSPC(IFRP,0) = - 1
C
C        If the user did not specify the core index range explictly we
C        take the occupied orbitals that do not belong to the active
C        set.
C
         IF (INDEX(TRA_CORSTR(IFRP),'not specified').NE.0) THEN
            CALL NUMCOR (NOCC(IFRP),NSTR(IFRP,1,1),WORK(KVEC(IFRP,1)),
     &                   NSPC(IFRP,1),WORK(KVEC(IFRP,5)))
            NSPC(IFRP,2) = 0
         ELSE
            CALL NUMLST(TRA_CORSTR(IFRP),WORK(KVEC(IFRP,5)),
     &                  NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                  IFRP,NSPC(IFRP,0))
            CALL ORBCNT(WORK(KVEC(IFRP,5)),NSPC(IFRP,0),
     &                  NPSH(IFRP),NESH(IFRP),
     &                  NSPC(IFRP,2),NSPC(IFRP,1))
         ENDIF
         IF (INDEX(TRA_CORSTR2(IFRP),'not specified').EQ.0) THEN
             WRITE(LUPRI,'(/A)') '* frozen open-shell added'
             Allocate (IBEIG(NR(IFRP)))
             NSPC2(IFRP,0)=1
             CALL NUMLST(TRA_CORSTR2(IFRP),IBEIG,
     &                  NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &                  IFRP,NSPC2(IFRP,0))
             CALL ORBCNT(IBEIG,NSPC2(IFRP,0),
     &                  NPSH(IFRP),NESH(IFRP),
     &                  NSPC2(IFRP,2),NSPC2(IFRP,1))
             CALL ADDCORE2(NR(IFRP),WORK(KVEC(IFRP,5)),NSPC(IFRP,0),
     &                     IBEIG,NSPC2(IFRP,0))
             CALL ORBCNT(WORK(KVEC(IFRP,5)),NSPC(IFRP,0),
     &                  NPSH(IFRP),NESH(IFRP),
     &                  NSPC(IFRP,2),NSPC(IFRP,1))
             DeAllocate (IBEIG)
             WRITE(LUPRI,*) 'WARNING: Frozen open shell is a continous',
     &                  ' set of ',NSPC2(IFRP,0),'. orbitals'
             NSPC2(IFRP,3)=NCORE2
         ELSE
             NSPC2(IFRP,3)=-1
         ENDIF
         NSPC(IFRP,0) = NSPC(IFRP,1) + NSPC(IFRP,2)
         CALL MEMREL('PAMTRA.4',WORK,KFRSAV,KVEC(IFRP,5),KFREE,LFREE)
         CALL MEMGET2('INTE','VEC',KVEC(IFRP,5),NSPC(IFRP,0),
     &               WORK,KFREE,LFREE)
         NDMOQC(1,IFRP,1) = NFBAS(IFRP,0)
         NDMOQC(2,IFRP,1) = NSPC(IFRP,0)
         NDMOQC(1,IFRP,2) = NFBAS(IFRP,0)
         NDMOQC(2,IFRP,2) = NSPC(IFRP,0)
      ENDDO
      IF (NSPC2(1,3).LT.0) THEN
        NCORE2=0
      ELSE
        I=0
        DO IFRP = 1,NFSYM
          I=I+NSPC2(IFRP,0)*2
        ENDDO
        NCORE2=I-NCORE2
      ENDIF
C####################################################
C
C     number of active occupied
C
      NAOCCT=0
      DO IFRP = 1, NFSYM
         NAOCC(IFRP) = NOCC(IFRP) - NSPC(IFRP,0)
         NAOCCT = NAOCCT + NAOCC(IFRP)
      ENDDO
C        
C     number of inactive occupied
C    
      NIOCCT=0
      DO IFRP = 1, NFSYM
         NIOCC(IFRP) = NOCC(IFRP) - NAOCC(IFRP)
         NIOCCT = NIOCCT + NIOCC(IFRP)
      ENDDO
C
C     number of active virtuals
C     
      NAVIRT=0
      DO IFRP = 1, NFSYM
         NAVIR(IFRP) = NSTR(IFRP,0,1) - NAOCC(IFRP)
         NAVIRT = NAVIRT + NAVIR(IFRP)
      ENDDO
C     
C     number of inactive virtuals
C 
      NIVIRT=0
      DO IFRP = 1, NFSYM 
         NIVIR(IFRP) = NORB(IFRP) - NOCC(IFRP) - NPSH(IFRP) -NAVIR(IFRP)
         NIVIRT = NIVIRT + NIVIR(IFRP)
      ENDDO
C####################################################
C
C     Print section
C     =============
C
      WRITE(LUPRI,'(/A)')
     &     '* Core orbital ranges for 2-index transformation:'
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,5)),0,0,0,NSPC)
      ENDDO
C
C     Calculate dimensions of the coefficient array
C
C     ICMOQC : Index to start of representation in the coefficient array
C
      NQT(1) = 0
      ICMOQC(1,1) = 1
      DO IFRP = 1, NFSYM
         NQ(IFRP,1) = NFBAS(IFRP,0)*NSPC(IFRP,0)*NZ
         NQT(1) = NQT(1) + NQ(IFRP,1)
         IF (IFRP.LT.NFSYM) ICMOQC(IFRP+1,1) = ICMOQC(IFRP,1)+NQ(IFRP,1)
      ENDDO
C
C     Allocate the memory for the coefficients
C
      CALL MEMGET2('REAL','QC mat',KQC(1),NQT(1),WORK,KFREE,LFREE)
      KQC(2) = KQC(1) + NQ(1,1)
C
C     Get all coefficients
C
      Allocate (CMO(NCMOTQ))
      Allocate (EIG(NORBT))
      Allocate (IBEIG(NORBT))
      IF (TOBEK) THEN
         CALL OPNFIL(LUKRMC,'KRMCSCF','OLD','PAMTRA')
         JRDMO = -1
         CALL RREADMO(CMO,JRDMO,1,LUKRMC)
         IF (JRDMO .EQ. 0) THEN
           WRITE(LUPRI,'(/A)')
     &     ' (PAMTRA)  Orbitals read from' //
     &     ' label NEWORB on file KRMCSCF'
         END IF
         CLOSE (LUKRMC)
      ELSE IF (TOBEKK) THEN
         CALL OPNFIL(LUKRMC,'KRMCOLD','OLD','PAMTRA')
         JRDMO = -1
         CALL RREADMO(CMO,JRDMO,3,LUKRMC)
         IF (JRDMO .EQ. 0) THEN
           WRITE(LUPRI,'(/A)')
     &     ' (PAMTRA)  Orbitals read from' //
     &     ' label NEWNATOB on file KRMCOLD'
         END IF
         CLOSE (LUKRMC)
      ELSE
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF (.NOT.TOBE) THEN
            WRITE(LUPRI,'(/A)') 'MO coefficient file not present'
            CALL QUIT('PAMTRA: No MO coefficient file !')
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PAMTRA')
            CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)
         ENDIF
      END IF
C
C     Select the set of MOs that we need for core transformation
C
      DO IFRP = 1, NFSYM
         CALL SELCFS (CMO(1+ICMOQ(IFRP)),IFRP,WORK(KQC(IFRP)),
     &                NSPC(IFRP,0),WORK(KVEC(IFRP,5)),
     &                NSPC(IFRP,2),NSPC(IFRP,1),
     &                NFBAS(IFRP,0),NORB(IFRP))
      ENDDO
C
C     Throw away the complete set, keep only the selected sets
      DeAllocate (CMO)
      DeAllocate (EIG)
      DeAllocate (IBEIG)
C
      CALL FLSHFO(LUPRI)
C

!.s/sya,2007.0907/ SK - 30-11-2009
      IF (PRPTRA .OR. PRPSYA ) CALL TRAPRP(WORK,KFREE,LFREE)
!.q
C
      CALL FLSHFO(LUPRI)
C
      IF (.NOT.NO2IND) THEN
         CALL FLSHFO(LUPRI)
         CALL TRAONE(TOTERG,WORK,KFREE,LFREE,NSTR,NSTRT,NSPC,
     &                  KQ,KQC,KE,KIBE,NSPC2)
      ENDIF
C
      CALL FLSHFO(LUPRI)

      if (tra_exacorr) then
#ifndef MOD_EXACORR
          call quit(
     &    'moltra error: tra_exacorr not included in this version')
#else
          call exa_moltra()
          no4ind = .true. ! skip conventional 4-index transformation
#endif
      end if

C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      IF (.NOT.NO4IND)
     &   CALL PAMTR1(WORK,KFREE,LFREE,IPRTRA,KQ,KE,KIBE,
     &               NDMOQR,ICMOQR,NSTR,.FALSE.,DUMMY,
     &               TRA_ANTIS,ITRA_INTFL4)
C
C      Print timing information
C
      CALL GETTIM(CPU2,WALL2)
      WALL   = WALL2 - WALL1
C
      IMINS  = INT(WALL)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(WALL) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(//A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total wall time used in PAMTRA :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CPU    = CPU2 - CPU1
      IMINS  = INT(CPU)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(CPU) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total CPU  time used in PAMTRA (master only) :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24/)') ' Transformation ended at : ',DAYTID
C
      CALL FLSHFO(LUPRI)

      ISTRAT = ISTRAT_save
      call dealloc(WORK)
      CALL QEXIT('PAMTRA')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck NUMCOR */
      SUBROUTINE NUMCOR (NO,NA,KAVEC,NC,KCVEC)

C***********************************************************************
C
C     Routine to find the core orbitals, given the positions of the
C     active orbitals and the number of occupied orbitals
C
C     Written  L.Visscher May 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      INTEGER KAVEC(NA),KCVEC(NO)
      LOGICAL CORE
C
      NC = 0
      DO I = 1, NO
         CORE = .TRUE.
         DO J = 1, NA
            IF (KAVEC(J).EQ.I) CORE = .FALSE.
         ENDDO
         IF (CORE) THEN
            NC = NC + 1
            KCVEC(NC) = I
         ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck PSIMP2 */
      SUBROUTINE PSIMP2()
C***********************************************************************
C
C     Driver routine for MP2-calculation
C
C     Written by T.Saue, J.Laerdahl and L.Visscher Jan 1997
C
C***********************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "dcbmp2.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
      LOGICAL TOBE,TRANOP
      DIMENSION NSTR(2,0:2,4),KVEC(2,2),KQ(2,4),KIBE(2,4),
     &          NR(2),IR(2,2),NQ(2),KE(2,4)
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4)

      Real(8), Allocatable :: CMO(:),EIG(:)
      Integer, Allocatable :: IBEIG(:)
      real(8), allocatable :: WORK(:)
C
      CALL QENTER('PSIMP2')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in PSIMP2')
      KFRSAV = KFREE
      TRANOP = .false.
      CALL MP2HI(IPRMP2)
C
C     Check to see if coefficients are there
C     ======================================
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
        WRITE(LUPRI,'(/A)')
     &    'PSIMP2 ERROR: MO coefficient file not present'
        CALL QUIT('PSIMP2: No MO coefficient file !')
      ELSE
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PSIMP2')
      ENDIF
C
C     Transfer information from DCBMP2 to DCBTRA
C     ==========================================
C
      TRANOP = NOPAIR
      NOPAIR = MP2NOP
      SCRBUF = SCRTRA
      SCRTRA = SCRMP2
      IPRBUF = IPRTRA
      IPRTRA = IPRMP2
      ISTBUF = ISTRAT
      ISTRAT = ISTRMP2
C     number of double quaternionic classes for strategy 3 and 4
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Set up index arrays for active occupieds
C     ========================================
C
      DO IFRP = 1,NFSYM
        IF(NOPAIR) THEN
          NR(IFRP)   =  NOCC(IFRP)
          IR(1,IFRP) =  1
          IR(2,IFRP) =  NOCC(IFRP)
        ELSE
          CALL QUIT('PSIMP2: e-p not implemented !')
          NR(IFRP)   =  NOCC(IFRP) + NPSH(IFRP)
          IR(1,IFRP) = -NPSH(IFRP)
          IR(2,IFRP) =  NOCC(IFRP)
        ENDIF
C
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,1),NR(IFRP),
     &     WORK,KFREE,LFREE)
        NSTR(IFRP,0,1) = - 1
        CALL NUMLST(MP2_INDSTR(1,IFRP),WORK(KVEC(IFRP,1)),
     &              NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &              IFRP,NSTR(IFRP,0,1))
        CALL ORBCNT(WORK(KVEC(IFRP,1)),NSTR(IFRP,0,1),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,1),NSTR(IFRP,1,1))
        NSTR(IFRP,0,1)   = NSTR(IFRP,1,1) + NSTR(IFRP,2,1)
      ENDDO
C
C     Set up index arrays for active virtuals
C     ========================================
C
      DO IFRP = 1,NFSYM
        IF(NOPAIR) THEN
          NR(IFRP)   = NSSH(IFRP)
          IR(1,IFRP) = NOCC(IFRP) + 1
          IR(2,IFRP) = NESH(IFRP)
        ELSE
          CALL QUIT('PSIMP2: e-p not implemented !')
          NR(IFRP)   =  NSSH(IFRP) + NPSH(IFRP)
          IR(1,IFRP) = -NPSH(IFRP)
          IR(2,IFRP) =  NESH(IFRP)
        ENDIF
C
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,2),NR(IFRP),
     &     WORK,KFREE,LFREE)
        NSTR(IFRP,0,2) = - 1
        CALL  NUMLST(MP2_INDSTR(2,IFRP),WORK(KVEC(IFRP,2)),
     &               NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &               IFRP,NSTR(IFRP,0,2))
        CALL ORBCNT(WORK(KVEC(IFRP,2)),NSTR(IFRP,0,2),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,2),NSTR(IFRP,1,2))
        NSTR(IFRP,0,2)   = NSTR(IFRP,1,2) + NSTR(IFRP,2,2)
      ENDDO
      
C
C     If energy threshold, truncate virtual space
C     ===========================================
C
      IF(DMP2_VIRTHR.NE.DUMMY) CALL VIRTHR(WORK,KFREE,LFREE,
     &     LUCOEF,NSTR(1,0,2),KVEC(1,2),DMP2_VIRTHR)
C
C     Calculate dimensions of coefficient arrays
C     The coefficients are stored as 2 arrays of two matrices each
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
C
C     Selected eigenvalues
C
      DO I = 1,2
        NDIM = 0
        DO IFRP = 1, NFSYM
          NDIM = NDIM + NSTR(IFRP,0,I)
        ENDDO
        CALL MEMGET2('REAL','EIGVAL',KE(1,I),NDIM,WORK,KFREE,LFREE)
        KE(2,I) = KE(1,I) + NSTR(1,0,I)
      ENDDO
C
C     Allocate memory for the integer array with information about the
C     spinors.
C
      DO I = 1, 2
       DO IFRP = 1, NFSYM
        CALL MEMGET2('INTE','IBEIG',KIBE(IFRP,I),NSTR(IFRP,0,I),
     &     WORK,KFREE,LFREE)
       ENDDO
      ENDDO
C
C     Selected coefficients
C
      DO I = 1,2
        NDIM = 0
        DO IFRP = 1, NFSYM
          NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
          NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
          ICMOQR(IFRP,I)   = NDIM + 1
          NQ(IFRP)         = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
          NDIM             = NDIM + NQ(IFRP)
        ENDDO
        CALL MEMGET2('REAL','Q mat',KQ(1,I),NDIM,WORK,KFREE,LFREE)
        KQ(2,I)   = KQ(1,I) + NQ(1)
      ENDDO
C
C     Get all coefficients and eigenvalues
C     ====================================
C
      KBUF = KFREE
      Allocate (CMO(N2BBASXQ))
      Allocate (EIG(NORBT))
      Allocate (IBEIG(NORBT))
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(/A)') 'MO coefficient file not present'
         CALL QUIT('PSIMP2: No MO coefficient file !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PSIMP2')
         CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,ESCF,14)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')
      WRITE(LUPRI,'(A,F18.10,/)')
     &  'Coefficients imported. Total DHF energy:', ESCF
C
C     Select the set that we need
C
      DO I = 1, 2
        ISAME(I)   = I
        DO IFRP = 1, NFSYM
          CALL SELCFS (CMO(1+ICMOQ(IFRP)),
     &         IFRP,WORK(KQ(IFRP,I)),NSTR(IFRP,0,I),
     &         WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &         NFBAS(IFRP,0),NORB(IFRP))
          CALL SELEIG(EIG(1+IORB(IFRP)),IFRP,WORK(KE(IFRP,I)),
     &         WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          CALL SELIBEIG(IBEIG,IORB(IFRP),IFRP,
     &                  WORK(KIBE(IFRP,I)),WORK(KVEC(IFRP,I)),
     &                  NSTR(IFRP,2,I),NSTR(IFRP,1,I))
        ENDDO
      ENDDO
C
C     Copy for indices 3 and 4
C
      DO I = 1,2
        KQ(1,I+2) = KQ(1,I)
        KQ(2,I+2) = KQ(2,I)
        KE(1,I+2) = KE(1,I)
        KE(2,I+2) = KE(2,I)
        KIBE(1,I+2) = KIBE(1,I)
        KIBE(2,I+2) = KIBE(2,I)
        ISAME(I+2) = I
        DO IFRP = 1,NFSYM
          NSTR(IFRP,0,I+2) = NSTR(IFRP,0,I)
          NSTR(IFRP,1,I+2) = NSTR(IFRP,1,I)
          NSTR(IFRP,2,I+2) = NSTR(IFRP,2,I)
          NDMOQR(1,IFRP,I+2) = NDMOQR(1,IFRP,I)
          NDMOQR(2,IFRP,I+2) = NDMOQR(2,IFRP,I)
          ICMOQR(IFRP,I+2)   = ICMOQR(IFRP,I)
        ENDDO
      ENDDO
C
C     If necessary recanonize virtual orbitals
C     ========================================
C
      CALL VIRCAN(WORK,KFREE,LFREE,IPRMP2,KQ(1,2),KE(1,2),
     &            KIBE(1,2),NDMOQR(1,1,2),ICMOQR(1,2),
     &            CMO)
C
C     Throw away the full coefficients and eigenvalues
C     ================================================
C
      DeAllocate (CMO)
      DeAllocate (EIG)
      DeAllocate (IBEIG)
C
C
C     Print section
C     =============
C
      IF (IPRMP2.GE.0) THEN
        WRITE(LUPRI,'(/1X,A/)') '* Active MP2 space:'
        DO I = 1, NFSYM
          CALL MP2PRI(I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),NSTR)
        ENDDO
        IF (IPRMP2.GE.1) THEN
          WRITE(LUPRI,'(/A)') 'Active Eigenvals:'
          DO I = 1, 2
            DO IFRP = 1,NFSYM
              WRITE(LUPRI,'(/2I8)') I, IFRP
              WRITE(LUPRI,'(5D15.8)')
     &             (WORK(KE(IFRP,I)+IXX-1), IXX=1,NSTR(IFRP,0,I))
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      CALL PAMTR1(WORK,KFREE,LFREE,IPRMP2,KQ,KE,KIBE,
     &            NDMOQR,ICMOQR,NSTR,.TRUE.,EMP2,
     &            MP2_ANTIS,MP2_INTFLG)
C
C     Print MP2 energy
C     ----------------
C
      WRITE(LUPRI,1000) 
      WRITE(LUPRI,1001) ESCF
      WRITE(LUPRI,1002) EMP2
      WRITE(LUPRI,1003) ESCF+EMP2
      CALL FLSHFO(LUPRI)
C
C     Recover information from DCBTRA
C
      NOPAIR = TRANOP
      SCRTRA = SCRBUF
      IPRTRA = IPRBUF
      ISTRAT = ISTBUF
C
      CALL MEMREL('PSIMP2',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      call dealloc(WORK)
      CALL QEXIT('PSIMP2')
C
      RETURN
 1000 FORMAT (/' Overview of calculated energies')
 1001 FORMAT('@ SCF energy :',T40,F25.15)
 1002 FORMAT('@ MP2 correlation energy :',T40,F25.15)
 1003 FORMAT('@ Total MP2 energy :',T40,F25.15)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck lagtau */
      SUBROUTINE LAGTAU(XMO,TMATRX1,TMATRX2,WORK,LWORK)
C***********************************************************************
C
C     Driver routine for tau dependent part of the Lagrangian
C
C     Written by Joost van Stralen -  may 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "maxorb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "mp1stpr.h"
      CHARACTER*24 DAYTID
      LOGICAL TOBE,TRASAM,ACTIVE
      DIMENSION TMATRX1(*),TMATRX2(*)
      DIMENSION XMO(*)
      DIMENSION WORK(*),NSTR(2,0:2,4),KVEC(2,5),KQ(2,4),KIBE(2,4),
     &          IR(2,2),NR(2),NQ(2,4),KE(2,4),
     &          NSPC(2,0:2),KQC(2),NSTRT(4)

      Real(8), Allocatable :: CMO(:),EIG(:)
      Integer, Allocatable :: IBEIG(:)
C     
      CALL QENTER('LAGTAU')
#include "memint.h"
      KFRSAV = KFREE
C     
      CALL GETTIM(CPU1,WALL1)
      CALL TRAHI(IPRTRA,5)
C     
C     Define number of double quaternionic classes
C     
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C
C     +-----------------------------------------+
C     | Set up index arrays for active orbitals |
C     +-----------------------------------------+
C     
C     index 1: occupied
C     =================
C
      DO IFRP = 1,NFSYM
        IF(NOPAIR) THEN
          NR(IFRP)   =  NOCC(IFRP)
          IR(1,IFRP) =  1
          IR(2,IFRP) =  NOCC(IFRP)
        ELSE
          NR(IFRP)   =  NOCC(IFRP) + NPSH(IFRP)
          IR(1,IFRP) = -NPSH(IFRP)
          IR(2,IFRP) =  NOCC(IFRP)
        ENDIF
C
        IF(NOCC(IFRP).EQ.0) THEN
          ISTART = 0
        ELSE
          ISTART = 1 + NIOCC(IFRP)
        ENDIF
        IEND = NOCC(IFRP)
        WRITE(MPR_INDSTR(1,IFRP),'(i4,A2,i4)')ISTART,'..',IEND
C
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,1),NR(IFRP),
     &     WORK,KFREE,LFREE)
        NSTR(IFRP,0,1) = - 1
        CALL NUMLST(MPR_INDSTR(1,IFRP),WORK(KVEC(IFRP,1)),
     &              NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &              IFRP,NSTR(IFRP,0,1))
        CALL ORBCNT(WORK(KVEC(IFRP,1)),NSTR(IFRP,0,1),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,1),NSTR(IFRP,1,1))
        NSTR(IFRP,0,1) = NSTR(IFRP,1,1) + NSTR(IFRP,2,1)
      ENDDO
C
C
C     index 2: virtual
C     =================
C
      DO IFRP = 1,NFSYM
        IF(NOPAIR) THEN
          NR(IFRP)   = NSSH(IFRP)
          IR(1,IFRP) = NOCC(IFRP) + 1
          IR(2,IFRP) = NESH(IFRP)
        ELSE
          NR(IFRP)   =  NSSH(IFRP) + NPSH(IFRP)
          IR(1,IFRP) = -NPSH(IFRP)
          IR(2,IFRP) =  NESH(IFRP)
        ENDIF
C
        ISTART = NOCC(IFRP) + 1
        IF(NAVIR(IFRP).EQ.0) THEN
          IEND = NOCC(IFRP) + 1
        ELSE
          IEND = NOCC(IFRP) + NAVIR(IFRP)
        ENDIF
        WRITE(MPR_INDSTR(2,IFRP),'(i4,A2,i4)')ISTART,'..',IEND
C
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,2),NR(IFRP),
     &     WORK,KFREE,LFREE)
        NSTR(IFRP,0,2) = - 1
        CALL NUMLST(MPR_INDSTR(2,IFRP),WORK(KVEC(IFRP,2)),
     &              NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &              IFRP,NSTR(IFRP,0,2))
        CALL ORBCNT(WORK(KVEC(IFRP,2)),NSTR(IFRP,0,2),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,2),NSTR(IFRP,1,2))
        NSTR(IFRP,0,2)   = NSTR(IFRP,1,2) + NSTR(IFRP,2,2)
      ENDDO
C
C
C     index 3: virtual
C     =================
C
      DO IFRP = 1,NFSYM
        IF(NOPAIR) THEN
          NR(IFRP)   = NSSH(IFRP)
          IR(1,IFRP) = NOCC(IFRP) + 1
          IR(2,IFRP) = NESH(IFRP)
        ELSE
          NR(IFRP)   =  NSSH(IFRP) + NPSH(IFRP)
          IR(1,IFRP) = -NPSH(IFRP)
          IR(2,IFRP) =  NESH(IFRP)
        ENDIF
C
        ISTART = NOCC(IFRP) + 1
        IF(NAVIR(IFRP).EQ.0) THEN
          IEND = NOCC(IFRP) + 1
        ELSE
          IEND = NOCC(IFRP) + NAVIR(IFRP)
        ENDIF
        WRITE(MPR_INDSTR(3,IFRP),'(i4,A2,i4)')ISTART,'..',IEND
C
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,3),NR(IFRP),
     &     WORK,KFREE,LFREE)
        NSTR(IFRP,0,3) = - 1
        CALL NUMLST(MPR_INDSTR(3,IFRP),WORK(KVEC(IFRP,3)),
     &              NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &              IFRP,NSTR(IFRP,0,3))
        CALL ORBCNT(WORK(KVEC(IFRP,3)),NSTR(IFRP,0,3),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,3),NSTR(IFRP,1,3))
        NSTR(IFRP,0,3)   = NSTR(IFRP,1,3) + NSTR(IFRP,2,3)
      ENDDO
C     
C       
C     index 4: occupied (Because we don't have a fourth index, because
C                        we do a 3/4 transformation, this is the 3rd
C                        index of the 2nd tau dependent part of the
C                        Lagrangian)
C     ================================================================
C         
      DO IFRP = 1,NFSYM
        IF(NOPAIR) THEN
          NR(IFRP)   =  NOCC(IFRP)
          IR(1,IFRP) =  1
          IR(2,IFRP) =  NOCC(IFRP)
        ELSE
          NR(IFRP)   =  NOCC(IFRP) + NPSH(IFRP)
          IR(1,IFRP) = -NPSH(IFRP)
          IR(2,IFRP) =  NOCC(IFRP)
        ENDIF
C    
        IF(NOCC(IFRP).EQ.0) THEN
          ISTART = 0
        ELSE
          ISTART = 1 + NIOCC(IFRP)
        ENDIF
        IEND = NOCC(IFRP)
        WRITE(MPR_INDSTR(4,IFRP),'(i4,A2,i4)')ISTART,'..',IEND
C
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,4),NR(IFRP),
     &      WORK,KFREE,LFREE)
        NSTR(IFRP,0,4) = - 1
        CALL NUMLST(MPR_INDSTR(4,IFRP),WORK(KVEC(IFRP,4)),
     &              NR(IFRP),IR(1,IFRP),IR(2,IFRP),
     &              IFRP,NSTR(IFRP,0,4))
        CALL ORBCNT(WORK(KVEC(IFRP,4)),NSTR(IFRP,0,4),
     &              NPSH(IFRP),NESH(IFRP),
     &              NSTR(IFRP,2,4),NSTR(IFRP,1,4))
        NSTR(IFRP,0,4)   = NSTR(IFRP,1,4) + NSTR(IFRP,2,4)
      ENDDO
C     
C     
C     Print section
C     =============
C     
      WRITE(LUPRI,'(/A)') '* Orbital ranges for 4-index transformation:'
      DO I = 1, NFSYM
         CALL TRAPRI(4,I,WORK(KVEC(I,1)),WORK(KVEC(I,2)),
     &                 WORK(KVEC(I,3)),WORK(KVEC(I,4)),NSTR)
      ENDDO
C     
C     Calculate dimensions of coefficient arrays 
C     The coefficients are stored as 4 arrays of two matrices each
C     They may share the same memory if the transformation ranges
C     are identical
C
C     NDMOQR : 1st dimension row/column dimension of coefficient array
C              2nd dimension fermion irrep
C              3rd dimension transformation index
C
C     ICMOQR : Index to start of representation in the coefficient array
C
      DO I = 1, 4
        NDIM = 0
         ICMOQR(1,I) = 1
         DO IFRP = 1, NFSYM
            NDMOQR(1,IFRP,I) = NFBAS(IFRP,0)
            NDMOQR(2,IFRP,I) = NSTR(IFRP,0,I)
            NQ(IFRP,I) = NFBAS(IFRP,0)*NSTR(IFRP,0,I)*NZ
            NDIM = NDIM + NQ(IFRP,I)
            IF (IFRP.LT.NFSYM)
     &         ICMOQR(IFRP+1,I) = ICMOQR(IFRP,I) + NQ(IFRP,I)
         ENDDO
        CALL MEMGET2('REAL','Q mat',KQ(1,I),NDIM,WORK,KFREE,LFREE)
        KQ(2,I)   = KQ(1,I) + NQ(1,I)
      ENDDO
C
C     Selected eigenvalues
C     
      DO I = 1,4
        NDIM = 0
        DO IFRP = 1, NFSYM
          NDIM = NDIM + NSTR(IFRP,0,I)
        ENDDO  
        CALL MEMGET2('REAL','EIGVAL',KE(1,I),NDIM,WORK,KFREE,LFREE)
        KE(2,I) = KE(1,I) + NSTR(1,0,I)
      ENDDO
C     
C     Allocate memory for the integer array with information about the
C     spinors.
C        
      DO I = 1, 4
       DO IFRP = 1, NFSYM
        CALL MEMGET2('INTE','IBEIG',KIBE(IFRP,I),NSTR(IFRP,0,I),
     &     WORK,KFREE,LFREE)
       ENDDO
      ENDDO 
C    
C     Get all coefficients
C     
      KBUF = KFREE
      Allocate (CMO(NCMOTQ))
      Allocate (EIG(NORBT))
      Allocate (IBEIG(NORBT))
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(/A)')
     &     'LAGTAU ERROR: MO coefficient file not present'
         CALL QUIT('LAGTAU: No MO coefficient file !')
      ELSE  
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','LAGTAU')
         CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)
         CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF 
C           
C     Select the set that we need
C        
      DO I = 1, 4
        ISAME(I)   = I
        DO IFRP = 1, NFSYM
          IF(NSTR(IFRP,0,I).GT.0) THEN
             CALL SELCFS (CMO(1+ICMOQ(IFRP)),IFRP,WORK(KQ(IFRP,I)),
     &                   NSTR(IFRP,0,I),WORK(KVEC(IFRP,I)),
     &                   NSTR(IFRP,2,I),NSTR(IFRP,1,I),
     &                   NFBAS(IFRP,0),NORB(IFRP))
             CALL SELEIG(EIG(1+IORB(IFRP)),IFRP,WORK(KE(IFRP,I)),
     &                 WORK(KVEC(IFRP,I)),NSTR(IFRP,2,I),NSTR(IFRP,1,I))
             CALL SELIBEIG(IBEIG,IORB(IFRP),IFRP,
     &                     WORK(KIBE(IFRP,I)),WORK(KVEC(IFRP,I)),
     &                     NSTR(IFRP,2,I),NSTR(IFRP,1,I))
          ENDIF
        ENDDO
      ENDDO
C     
C     Throw away the complete set, keep only the selected sets
C     
      DeAllocate (CMO)
      DeAllocate (EIG)
      DeAllocate (IBEIG)
C     
C     If necessary recanonize (virtual) orbitals
C     Assuming that all sets are equal, we use the first index !!!
C     ============================================================
C     
      IF (RCORBS) THEN 
C        Overdimension CMO because it is used as Focksize scratch space
         Allocate (CMO(N2BBASXQ))
         CALL VIRCAN(WORK,KFREE,LFREE,IPRMP2,KQ(1,2),KE(1,2),
     &            KIBE(1,2),NDMOQR(1,1,2),ICMOQR(1,2),CMO)
         DeAllocate (CMO)
      ENDIF
C
C     Call driver for 4-index transformation
C     --------------------------------------
C
      IF (.NOT.NO4IND)
     &   CALL LAGTAU1(XMO,TMATRX1,TMATRX2,WORK,KFREE,LFREE,IPRTRA,KQ,KE,
     &               KIBE,NDMOQR,ICMOQR,NSTR,
     &               TRA_ANTIS,ITRA_INTFL4)
C
C      Print timing information
C
      CALL GETTIM(CPU2,WALL2)
      WALL   = WALL2 - WALL1
C
      IMINS  = INT(WALL)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(WALL) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(//A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total wall time used in LAGTAU :',
     &      IHOURS,':',IMINS,':',ISECS
C
      CPU    = CPU2 - CPU1
      IMINS  = INT(CPU)/60
      IHOURS = IMINS/60
      IMINS  = IMINS - 60*IHOURS
      ISECS  = NINT(CPU) - 3600*IHOURS - 60*IMINS
      WRITE(LUPRI,'(A,I5.2,A,I2.2,A,I2.2)')
     &      ' Total CPU  time used in LAGTAU :',
     &      IHOURS,':',IMINS,':',ISECS
C
      WRITE(LUPRI,'(/A,F20.6)')
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24)') ' Transformation ended at : ',DAYTID
C
      CALL QEXIT('LAGTAU')
      CALL FLSHFO(LUPRI)
C     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MP2PRI */
      SUBROUTINE MP2PRI(I,IOS1,IOS2,NSTR)
C***********************************************************************
C
C     Initial print routine for MP2 module
C
C     Written by T.Saue Jan 19 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION IOS1(*),IOS2(*),NSTR(2,0:2,4)
C
      WRITE(LUPRI,'(1X,A,A3)') '* Fermion ircop ',FREP(I)
C
C     Index 1: Active occupied
C
      IF(NSTR(I,0,1).GT.0) THEN
        WRITE(LUPRI,'(6X,A,I5)')
     &     '- Active occupied orbitals:', NSTR(I,0,1)
        WRITE(LUPRI,1002) (IOS1(K),K=1,NSTR(I,0,1))
      ELSE
        WRITE(LUPRI,'(6X,A)')
     &     '- No active occupied orbitals'
      ENDIF
C
C     Index 2: Active virtuals
C
      IF(NSTR(I,0,2).GT.0) THEN
        WRITE(LUPRI,'(6X,A,I5)')
     &     '- Active virtual orbitals: ', NSTR(I,0,2)
        WRITE(LUPRI,1002) (IOS2(K),K=1,NSTR(I,0,2))
      ELSE
        WRITE(LUPRI,'(6X,A/)')
     &     '- No active virtual orbitals'
      ENDIF
C
 1002 FORMAT (6X,10I5/)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck TRAPRI */
      SUBROUTINE TRAPRI(NINDEX,I,IOS1,IOS2,IOS3,IOS4,NSTR)
C***********************************************************************
C
C     Initial print routine for 4-index transformation module
C
C     Written by T.Saue Jan 19 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION IOS1(*),IOS2(*),IOS3(*),IOS4(*)
      DIMENSION NSTR(2,0:2,NINDEX)
C
      WRITE(LUPRI,'(//A,A3)') ' * Fermion ircop ',FREP(I)
C
C     Index 1
C
      IF(NSTR(I,0,1).GT.0) THEN
        WRITE(LUPRI,1000) 1,NSTR(I,0,1)
        WRITE(LUPRI,1002) (IOS1(K),K=1,NSTR(I,0,1))
      ELSE
        WRITE(LUPRI,1001) 1
      ENDIF
C
C       Index 2
C
      IF (NINDEX.LT.2) RETURN
      IF(NSTR(I,0,2).GT.0) THEN
        WRITE(LUPRI,1000) 2,NSTR(I,0,2)
        WRITE(LUPRI,1002) (IOS2(K),K=1,NSTR(I,0,2))
      ELSE
        WRITE(LUPRI,1001) 2
      ENDIF
C
C       Index 3
C
      IF (NINDEX.LT.3) RETURN
      IF(NSTR(I,0,3).GT.0) THEN
        WRITE(LUPRI,1000) 3,NSTR(I,0,3)
        WRITE(LUPRI,1002) (IOS3(K),K=1,NSTR(I,0,3))
      ELSE
        WRITE(LUPRI,1001) 3
      ENDIF
C
C       Index 4
C
      IF (NINDEX.LT.4) RETURN
      IF(NSTR(I,0,4).GT.0) THEN
        WRITE(LUPRI,1000) 4,NSTR(I,0,4)
        WRITE(LUPRI,1002) (IOS4(K),K=1,NSTR(I,0,4))
      ELSE
        WRITE(LUPRI,1001) 4
      ENDIF
C
 1000 FORMAT (/6X,'Index ',I1,I5,' orbitals')
 1001 FORMAT (/6X,'No orbitals for index ',I1)
 1002 FORMAT (6X,12I5)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck PAMTR1 */
      SUBROUTINE PAMTR1(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &           NDMOQR,ICMOQR,NSTR,LMP2,EMP2,ANTIS,INTFLG)
C***********************************************************************
C
C     Subdriver for 4index transformation.
C
C     KQ(ifrp,n)     - adresses of coefficient set for index n
C     KE(ifrp,2)     - adresses of eigenvalues (for MP2 calculation)
C     KIBE(ifrp,n)   - adresses of single group identification (for spinfree calculation)
C     NSTR(ifrp,0,i) - total number of orbitals for index i
C     NSTR(ifrp,1,i) - number of electronic orbitals for index i
C     NSTR(ifrp,2,i) - number of positronic orbitals for index i
C     LMP2  = .TRUE. - evaluate MP2 energy
C     ANTIS = .TRUE. - anti-symmetrize integrals
C     note TODO: ANTIS is not implemented!
C     INTFLG - flag of what integral types to transform
C
C     Written by Luuk Visscher, Trond Saue and J.Laerdahl Jan 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),NDMOQR(2,2,4),ICMOQR(2,4),
     &          KQ(2,4),KE(2,4),KIBE(2,4),DINTSKP(4,2),WORK(*)
C
#include "dcbibt.h"
C
      KFRSAV = KFREE
      call memchk('Start of PAMTR1',WORK,1)
C
C     Set up symmetry packing information.
C
      CALL MEMGET2('INTE','INDX',KINDX,3*NTBAS(0),WORK,KFREE,LFREE)
      CALL GISPCK(NSPCK,ISPCK,0,0,0,IDUM,IDUM,IDUM,IDUM,IDUM,IDUM)
C.....Create pointer array INDX:
C     INDX(1,INDA) : Position of function in shell
C     INDX(2,INDA) : Irreducible representation of function
C     INDX(3,INDA) : Position of function within this particular irrep
C.....MAXINSH is the maximum numbers of functions in a shell
      MXINSH = NINSH(0,-1,WORK(KINDX),IPRINT)
C
C     Define the MS integral types that we are going to make
C
      IF (NZ.EQ.1) THEN
         NCLASS = 4
      ELSEIF (NZ.EQ.2) THEN
         NCLASS = 8
      ELSE
         NCLASS = 16
      ENDIF
C
      CALL DEFCL (NCLASS)
C
C
C     We have four possible strategies :
C
C     1) In-core test version : half- and fully-transformed integrals
C        should fit in memory
C     2) Loop over restricted range of MS-labels such that all transf.
C        integrals fit in memory
C        Loop over distributions transform 4 indices(rest. ranges first)
C        Accumulate integrals in memory, calculate MP2 energy
C     3) Loop over distributions, transform 2 indices and write to disk
C        Loop over transformed indices and transform last 2 indices
C     4) Like 3, but 2-index transformed integrals written to disk in
C        batches; parallel scheme.
C
C     In all cases we transform in mulliken notation
C
      CALL QENTER('PAMTR1')
C
C
C     Prepare for screening
C
      IF(SCRTRA .GT. 0.0d0) THEN
        CALL PR4SCR(KGAB,KDRIJ,WORK,KFREE,LFREE,DINTSKP,
     &              WORK(KQ(1,1)),WORK(KQ(1,2)),
     &              WORK(KQ(1,3)),WORK(KQ(1,4)),
     &              NDMOQR,ICMOQR,NSTR,ISAME,IPRINT)
      ELSE
        KDRIJ = KFREE
        KGAB  = KFREE
      ENDIF
C
      IF (ANTIS) THEN
         CALL QUIT('PAMTR1: .ANTIS is not implemented!')
      ENDIF
      IF (ISTRAT.EQ.1) THEN
          IF (LMP2)
     &    CALL QUIT('PAMTR1: MP2 not implemented for strategy 1')
          CALL TRDR1T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,NSTR,
     &                ANTIS,LMP2,KINDX,KQ,KE,KIBE,DINTSKP)
      ELSEIF (ISTRAT.EQ.2) THEN
          IF (.NOT.LMP2)
     &    CALL QUIT('PAMTR1: Only MP2 implemented for strategy 2')
          CALL TRDR2T(WORK,KFREE,LFREE,IPRINT,INTFLG,NDMOQR,ICMOQR,NSTR,
     &                ANTIS,LMP2,EMP2,KINDX,KQ,KE,KIBE,DINTSKP)
      ELSEIF (ISTRAT.EQ.3) THEN
          IF (LMP2)
     &    CALL QUIT('PAMTR1: MP2 not implemented for strategy 3')
          CALL TRDR3T(WORK,KFREE,LFREE,IPRINT,INTFLG,NSTR,
     &                ANTIS,LMP2,WORK(KINDX),KQ,KE,KIBE,DINTSKP)
      ELSEIF (ISTRAT.EQ.4) THEN
          IF (LMP2)
     &    CALL QUIT('PAMTR1: MP2 not implemented for strategy 4')
          CALL TRDR4T(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                KGAB,KDRIJ,NSTR,
     &                ANTIS,LMP2,KINDX,KQ,KE,KIBE,DINTSKP)
      ELSEIF (ISTRAT.EQ.5) THEN
          CALL QUIT('PAMTR1: strategy 5 can only do lagrangians')
      ELSEIF (ISTRAT.EQ.6) THEN
          CALL TRDR6T(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                KGAB,KDRIJ,NSTR,
     &                ANTIS,LMP2,KINDX,KQ,KE,KIBE,DINTSKP)
      ELSE
        WRITE(LUPRI,'(A)')
     &        'PAMTR1: Strategy ',ISTRAT,' not implemented'
        CALL QUIT('PAMTR1: Unknown strategy !')
      ENDIF
C
      CALL MEMREL('PAMTR1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('PAMTR1')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck lagtau1 */
      SUBROUTINE LAGTAU1(XMO,TMATRX1,TMATRX2,WORK,KFREE,LFREE,IPRINT,KQ,
     &                   KE,KIBE,NDMOQR,ICMOQR,NSTR,ANTIS,INTFLG)
C***********************************************************************
C       
C     Subdriver for 4index transformation.
C       
C     KQ(ifrp,n)     - adresses of coefficient set for index n
C     KE(ifrp,2)     - adresses of eigenvalues (for MP2 calculation)
C     KIBE(ifrp,n)   - adresses of single group identification (for spinfree calculation)
C     NSTR(ifrp,0,i) - total number of orbitals for index i
C     NSTR(ifrp,1,i) - number of electronic orbitals for index i
C     NSTR(ifrp,2,i) - number of positronic orbitals for index i
C     ANTIS = .TRUE. - anti-symmetrize integrals
C     note TODO: ANTIS is not implemented!
C     INTFLG - flag of what integral types to transform
C
C     Written by Joost van Stralen september 2003
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dgroup.h"
C
      LOGICAL ANTIS
      DIMENSION TMATRX1(*),TMATRX2(*)
      DIMENSION XMO(*)
      DIMENSION NSTR(2,0:2,4),NDMOQR(2,2,4),ICMOQR(2,4),
     &          KQ(2,4),KE(2,4),KIBE(2,4),DINTSKP(4,2),WORK(*)
C    
#include "dcbibt.h"
C
      KFRSAV = KFREE
C
C     Set up symmetry packing information.
C
      CALL MEMGET2('INTE','INDX',KINDX,3*NTBAS(0),WORK,KFREE,LFREE)
      CALL GISPCK(NSPCK,ISPCK,0,0,0,IDUM,IDUM,IDUM,IDUM,IDUM,IDUM)
      MXINSH = NINSH(0,-1,WORK(KINDX),IPRINT)
C
C     Define the MS integral types that we are going to make
C
      IF (NZ.EQ.1) THEN
         NCLASS = 4
      ELSEIF (NZ.EQ.2) THEN
         NCLASS = 8
      ELSE
         NCLASS = 16
      ENDIF
C
      CALL DEFCL (NCLASS)
C
      CALL QENTER('LAGTAU1')
C
      IF (ANTIS) THEN
         CALL QUIT('LAGTAU1: .ANTIS is not implemented!')
      ENDIF
C
      CALL TRDR5T(XMO,TMATRX1,TMATRX2,WORK,KFREE,LFREE,IPRINT,INTFLG,
     &            NSTR,ANTIS,LMP2,KINDX,KQ,KE,KIBE,DINTSKP)
C
      CALL MEMREL('LAGTAU1',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      CALL QEXIT('LAGTAU1')
C     
      RETURN
C     
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck MP2HI */
      SUBROUTINE MP2HI(IPRINT)
C
C     Written by Trond Saue and Luuk Visscher Jan 17 1996
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbmp2.h"
C
C     Now say hi and give orbital information
C
      CALL TITLER('MP2 Calculation ','*',128)
      WRITE (LUPRI,'(/10X,A)') 'Written by Jon Laerdahl'//
     &     ', Luuk Visscher & Trond Saue'
      WRITE (LUPRI,'(10X,A//)')
     &     'Odense (Jan. 1997), Oslo (June 1997), Auckland (Sept. 1998)'
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('SETDC2: Distribution of symmetry orbitals '//
     &     'among fermion ircops',-1)
        WRITE(LUPRI,'(7X,3A5)')
     &       'Total',' Lbas',' Sbas'
        WRITE(LUPRI,'(A3,4X,3I5)')
     &        (FREP(I),(NFBAS(I,J),J=0,2),I = 1,NFSYM)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck trahi */
      SUBROUTINE TRAHI(IPRINT,IMODUL)
C
C     Written by Luuk Visscher 05-11-1996.
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dgroup.h"
#include "maxorb.h"
#include "dcbtra.h"
      CHARACTER*24 DAYTID
      LOGICAL TOBE,TOBEK
C
C     Now say hi and give orbital information
C
      CALL TITLER('Transformation to Molecular Spinor Basis','*',116)
      WRITE (LUPRI,'(/10X,A)') 'Written by Luuk Visscher'//
     &     ', Jon Laerdahl & Trond Saue'
      WRITE (LUPRI,'(10X,A//)') 'Odense, 1997'
      IF (IMODUL.EQ.1) THEN
         CALL TITLER('Transformation of property integrals','*',116)
      ELSEIF (IMODUL.EQ.2) THEN
         CALL TITLER('Transformation of core Fock matrix','*',116)
      ELSEIF (IMODUL.EQ.3) THEN
         CALL TITLER('Transformation of 2-electron integrals','*',116)
      ELSEIF (IMODUL.EQ.4) THEN
         CALL TITLER('Transformation of MCSCF 2-electron integrals',
     &        '*',116)
      ELSEIF (IMODUL.EQ.5) THEN
         CALL TITLER('Formation of Lagrangian','*',116)
      ENDIF
C
      IF (SMLV1C)
     & WRITE(LUPRI,'(/A)') ' WARNING : Incomplete projection on 1-',
     & ' center integrals with SMLV1C : results may be inconsistent'
C
      CALL GTINFO(DAYTID)
      WRITE(LUPRI,'(/A,A24)') ' Transformation started at : ',DAYTID
C
      IF(IPRINT.GE.3) THEN
        CALL HEADER('SETDC2: Distribution of symmetry orbitals '//
     &     'among fermion ircops',-1)
 1      WRITE(LUPRI,'(7X,3A5)')
     &       'Total',' Lbas',' Sbas'
        WRITE(LUPRI,'(A3,4X,3I5)')
     &        (FREP(I),(NFBAS(I,J),J=0,2),I = 1,NFSYM)
      ENDIF
C
      IF (IMODUL .NE. 4 ) THEN
         INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
         INQUIRE(FILE='KRMCOLD',EXIST=TOBEK)
         IF (.NOT.TOBEK) THEN
            INQUIRE(FILE='DFCOEF',EXIST=TOBE)
            IF (.NOT.TOBE) THEN
               WRITE(LUPRI,'(/A)') 'MO coefficient file not present'
               CALL QUIT('TRAHI: No MO coefficient file !')
            ELSE
               CALL OPNFIL(LUCOEF,'DFCOEF','OLD','TRAHI')
               CALL REACMO(LUCOEF,'DFCOEF',DUM,DUM,IDUM,TOTERG,1)
               CLOSE(LUCOEF,STATUS='KEEP')
            ENDIF
         END IF
      END IF
C
      CALL FLSHFO(LUPRI)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck MS4IN1 */
      SUBROUTINE MS4IN1(WORK,KFREE,LFREE,IPRINT,ICS,ICF,TRIAN,
     &                  NSTR1,NSTR2,NDMOQR,ICMOQR,IJPAIRS,Q1,Q2,
     &                  GMAT,HMAT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher January 1997.
C
C     PURPOSE : Do 4-index transformation to molecular spinor basis
C               First index-pair transformation.
C
C     Input :
C
C     - IPRINT       Print flag
C     - ICS          First class of integrals 1 : (LL|XX), 2 : (SS|XX), 3 : (SL|XX)
C     - ICL          Last class of integrals
C     - TRIAN        Store the transformed integrals as a lower triangle
C     - NSTR1        Number of active spinors for index 1
C     - NSTR2        Number of active spinors for index 2
C     - NFPCK12      Number of spinor pairs for each compound symmetry
C     - IFPCK12      Pointers to spinor pairs
C     - NDMOQR       Dimensions of the coefficient array
C     - ICMOQR       Pointers to coefficients
C     - IJPAIRS      Number of blocks for each boson symmetry
C     - Q1           Coefficients for index 1
C     - Q2           Coefficients for index 2
C     - GMAT         Symmetry packed scalar integrals
C     - HMAT         Zero matrix
C
C     Output :
C
C     - HMAT         Symmetry packed half-transformed integrals
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0, D1 = 1.0D0)
C
      LOGICAL   TRIAN
      INTEGER   NSTR1(2),NSTR2(2)
      DIMENSION WORK(*)
      DIMENSION Q1(*),Q2(*),GMAT(*),HMAT(*)
      DIMENSION IJPAIRS(0:7)
      DIMENSION NDMOQR(2,2,2),ICMOQR(2,2)
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbtra.h"
#include "dcbibt.h"
C
      CALL QENTER('MS4IN1')
C
      KFRSAV = KFREE
C
C     Allocate memory for temporary storage of unsorted block of
C     half-transformed integrals. Allocate absolute maximum instead
C     of the exact maximum size.
C
      IF (NFSYM .GT. 1) THEN
         LHTMP = MAX(NFPCK12(1),NFPCK12(2))*NZ
      ELSE
         LHTMP = 0
      END IF
      CALL MEMGET2('REAL','HTMP',KHTMP,LHTMP,WORK,KFREE,LFREE)
      IF (TRIAN) THEN
         IF (NFSYM .EQ. 1) THEN
            NI = NSTR1(1)
            NJ = NSTR2(1)
         ELSE
            NI = MAX(NSTR1(1),NSTR1(2))
            NJ = MAX(NSTR2(1),NSTR2(2))
         END IF
         CALL MEMGET2('REAL','BMAT',KBMAT,NI*NJ*NZ,WORK,KFREE,LFREE)
      ELSE
         CALL MEMGET2('REAL','BMAT',KBMAT,0,WORK,KFREE,LFREE)
      END IF
C
      IXYZ = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C
C     Transform lefthand indices to molecular spinor basis
C
      IOFF = 1
      JOFF = 1
C
C     Loop over the boson symmetries of the densities
C
      DO IREPPQ = 0, NBSYM-1
C
C        The parity of the spinor product to which this density
C        contributes need to be established. Since the large and
C        small component have opposite parity the assignment is
C        different in case of Gaunt integrals.
C
         IREPIJ = JBTOF(IREPPQ,1)
         IF (ICS.EQ.3) IREPIJ = JBTOF(IREPPQ,2)
C
         IF (NFSYM.GT.1) CALL DZERO (WORK(KHTMP),NFPCK12(IREPIJ)*NZ)
C
C        Loop over the distributions with this boson symmetry in this
C        shell-block batch
C
         DO IJ = 1, IJPAIRS(IREPPQ)
C
C           Loop over the compound component (LL, SS, SL) of the lefthand
C           densities
C
            IF (ICS.EQ.2) IOFF = IOFF + NSPCK(IREPPQ,1)
            DO IC = ICS, ICF
               CALL ICTYPES (IC,IC1,IC2)
C
C              Loop over the boson symmetries of the second
C              (untransformed) index. Get the fermion symmetry (parity)
C              of that index and get the boson and fermion symmetry
C              of the first index.
C
               DO IREPQ = 0, NBSYM-1
                  IREPP = IBTXOR(IREPPQ,IREPQ)
                  IREPI =  JBTOF(IREPP,IC1)
                  IREPJ =  JBTOF(IREPQ,IC2)
                  KOFF = ICMOQR(IREPI,1) + IBBAS(IREPP,IC1) 
     &                 - IBAS(IREPI)
                  LOFF = ICMOQR(IREPJ,2) + IBBAS(IREPQ,IC2) 
     &                 - IBAS(IREPJ)
                  NP = NBBAS(IREPP,IC1)
                  NQ = NBBAS(IREPQ,IC2)
                  NI = NSTR1(IREPI)
                  NJ = NSTR2(IREPJ)
                  NRQ1 = NDMOQR(1,IREPI,1)
                  NCQ1 = NDMOQR(2,IREPI,1)
                  NRQ2 = NDMOQR(1,IREPJ,2)
                  NCQ2 = NDMOQR(2,IREPJ,2)
                  IF (IC1.EQ.1) THEN
                     IREPPI = IREPP
                  ELSE
                     IREPPI = IBTXOR(IXYZ,IREPP)
                  ENDIF
                  IF (IC2.EQ.1) THEN
                     IREPQI = IREPQ
                  ELSE
                     IREPQI = IBTXOR(IXYZ,IREPQ)
                  ENDIF
C
C                 Do the actual transformation.
C
                  IF (.NOT.TRIAN.OR.(TRIAN.AND.IREPI.GT.IREPJ)) THEN
C
C                 Square transformation
C
                     IF ((NP*NQ.NE.0).AND.(NI*NJ.NE.0)) THEN
                     IF (NFSYM.GT.1) THEN
                       CALL QTRANS('AOMO','S',D0,NP,NQ,NI,NJ,
     &                    GMAT(IOFF),NP,NQ,1,IPQTOQ(1,0),
     &                    WORK(KHTMP),NI,NJ,NZ,IPQTOQ(1,IREPPQ),
     &                    Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                    Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                    WORK(KFREE),LFREE,IPRINT)
C
C                      Order the integrals with NKL as first and IZ
C                      as second index
C
                       DO IZ = 1, NZ
                          JOFF1 = JOFF + (IZ-1)*NFPCK12(IREPIJ) +
     &                                   IFPCK12(IREPI,IREPJ)
                          JOFF2 = KHTMP + (IZ-1)*NI*NJ
                          CALL DAXPY(NI*NJ,D1,WORK(JOFF2),1,
     &                               HMAT(JOFF1),1)
                       ENDDO
                     ELSE
                       CALL QTRANS('AOMO','S',D1,NP,NQ,NI,NJ,
     &                    GMAT(IOFF),NP,NQ,1,IPQTOQ(1,0),
     &                    HMAT(JOFF),NI,NJ,NZ,IPQTOQ(1,IREPPQ),
     &                    Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                    Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                    WORK(KFREE),LFREE,IPRINT)
                     ENDIF
                     ENDIF
C
                  ELSEIF (TRIAN.AND.IREPI.EQ.IREPJ) THEN
C
C                 Triangular storage after square transformation
C
                     IF ((NP*NQ.NE.0).AND.(NI*NJ.NE.0)) THEN
                       CALL QTRANS('AOMO','S',D0,NP,NQ,NI,NJ,
     &                           GMAT(IOFF),NP,NQ,1,IPQTOQ(1,0),
     &                           WORK(KBMAT),NI,NJ,NZ,IPQTOQ(1,IREPPQ),
     &                           Q1(KOFF),NRQ1,NCQ1,NZ,IPQTOQ(1,IREPPI),
     &                           Q2(LOFF),NRQ2,NCQ2,NZ,IPQTOQ(1,IREPQI),
     &                           WORK(KFREE),LFREE,IPRINT)
                       NII = NI*(NI+1)/2
                       IF (NFSYM.GT.1) THEN
                         CALL QTRIAN('PACK',D0,NI,NZ,WORK(KBMAT),NI,NJ,
     &                               WORK(KHTMP),NII)
C
C                        Order the integrals with NKL as first and IZ
C                        as second index
C
                         DO IZ = 1, NZ
                            JOFF1 = JOFF + (IZ-1)*NFPCK12(IREPIJ) +
     &                                     IFPCK12(IREPI,IREPJ)
                            JOFF2 = KHTMP + (IZ-1)*NII
                            CALL DAXPY(NII,D1,WORK(JOFF2),1,
     &                                 HMAT(JOFF1),1)
                         ENDDO
                       ELSE
                         CALL QTRIAN('PACK',D1,NI,NZ,WORK(KBMAT),NI,NJ,
     &                               HMAT(JOFF),NII)
                       ENDIF
                     ENDIF
                  ENDIF
                  IOFF = IOFF + NP * NQ
               ENDDO
C              ... END DO IREPQ
            ENDDO
C           ... END DO IC
            JOFF = JOFF + NFPCK12(IREPIJ)*NZ
         ENDDO
C        ... END DO IJ
      ENDDO
C     ... END DO IREPPQ
C
      CALL MEMREL('MS4IN1',WORK,1,KFRSAV,KFREE,LFREE)
      CALL QEXIT('MS4IN1')
C
      RETURN
      END
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck caldis */
      SUBROUTINE CALDIS(I2TYP,GMAT,INDX,INDXAB,NODV,NOPV,
     &                  GABRAO,DRIJ,DINTSKP,SCRTRA,WORK,LWORK,IPRINT)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Luuk Visscher 22-10-1996
C
C     PURPOSE : Driver routine for the calculation of the two-electron
C               distributions (**|cd).
C               The distributions are stored as full squares with
C               symmetry reduction.
C
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C
#include "comdis.h"
#include "maxaqn.h"
#include "mxcent.h"
      DIMENSION GMAT(*),INDX(3,*),INDXAB(*),GABRAO(*),DRIJ(*)
      DIMENSION WORK(LWORK),DINTSKP(*)
#include "cbihr2.h"
#include "dorps.h"
#include "ccom.h"
#include "dcborb.h"
#include "dgroup.h"
C
      LOGICAL NODV, NOPV, NOCONT
C
C
C-------------------------------------------------------------
C     Setup information for the two-electron integralroutines.
C-------------------------------------------------------------
C
      NODV = .FALSE.
      NOPV = .FALSE.
      IF (NASHT .EQ. 0) NODV = .TRUE.
      IF (NASHT .LT. 2) NOPV = .TRUE.
      IPRALL =  0
      DO 100 I = 0,7
         DOREPS(I) = .TRUE.
  100 CONTINUE
      DO 110 I = 1,MXCENT
         DOCOOR(1,I) = .TRUE.
         DOCOOR(2,I) = .TRUE.
         DOCOOR(3,I) = .TRUE.
  110 CONTINUE
      NOCONT = .FALSE.
C
      MAXDIF = 0
      ITYPE  = 4
C
C     Initialize buffer file if necessary
C
      IF (NPASS.GT.1) CALL INITGBF(LGFIL,NGBFSZ,NPASS,IRECG,
     &                             GMAT(KLGREC),GMAT(KLGBUF))
C
      CALL TWOINT(WORK,LWORK,DUMMY,DUMMY,NFSYM,IDUMMY,IDUMMY,GMAT,INDX,
     &            INDXAB,ITYPE,MAXDIF,0,NODV,NOPV,NOCONT,
     &            TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,IPRNTD,RTNTWO,
     &            IDUMMY,I2TYP,IDUMMY,SCRTRA,GABRAO,DRIJ,DUMMY,DINTSKP,
     &            .TRUE.,.false.,DUMMY,DUMMY)
C
C     Finish buffer file if necessary
C
      IF (NPASS.GT.1) CALL FINAGBF (LGFIL,NPASS,IRECG,GMAT(KLGREC),
     &                              GMAT(KLGBUF),NGBFSZ,GMAT(KRGBUF),
     &                              GMAT(KIGBUF))
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SHRNGE */
      SUBROUTINE SHRNGE (IC,INTFLG,I2TYP,IASTRT,IBSTRT,IASMAX,IBSMAX)
C
C     Written by Luuk Visscher, january 1997
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#include "implicit.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dcbbas.h"
C
      LOGICAL L2TYP(0:4),LBIT
C
      L2TYP(0) = INTFLG.EQ.0
      L2TYP(1) = LBIT(INTFLG,1)
      L2TYP(2) = LBIT(INTFLG,2)
      L2TYP(3) = LBIT(INTFLG,3)
      L2TYP(4) = LBIT(INTFLG,4)
C
      IF (IC.EQ.1) THEN
         IASTRT = 1
         IBSTRT = 1
         IASMAX = NLRGBL
         IBSMAX = NLRGBL
         IF (L2TYP(0)) THEN
C           Non-relativistic option
            IASTRT = 1
            IBSTRT = 1
            IASMAX = MAXSHL
            IBSMAX = MAXSHL
            I2TYP = 0
         ELSEIF (L2TYP(1).AND.L2TYP(2)) THEN
C           (LL|LL) + (LL|SS) - integrals
            I2TYP = 12
         ELSEIF (L2TYP(1)) THEN
C           (LL|LL) - integrals
            I2TYP = 1
         ELSEIF (L2TYP(2)) THEN
C           (LL|SS) - integrals
            I2TYP = 2
         ELSE
C           No valid class of integrals for this value of IC
            IASTRT = 1
            IBSTRT = 1
            IASMAX = 0
            IBSMAX = 0
            I2TYP = -1
         ENDIF
      ELSEIF (IC.EQ.2) THEN
         IASTRT = NLRGBL + 1
         IBSTRT = NLRGBL + 1
         IASMAX = MAXSHL
         IBSMAX = MAXSHL
         IF (L2TYP(0)) THEN
C           Non-relativistic option
            I2TYP = -1
         ELSEIF (L2TYP(2).AND.L2TYP(3)) THEN
C           (SS|LL) + (SS|SS) - integrals
            I2TYP = 12
         ELSEIF (L2TYP(2)) THEN
C           (SS|LL) - integrals
            I2TYP = 1
         ELSEIF (L2TYP(3)) THEN
            I2TYP = 2
C           (SS|SS) - integrals
         ELSE
C           No valid class of integrals for this value of IC
            IASTRT = 1
            IBSTRT = 1
            IASMAX = 0
            IBSMAX = 0
            I2TYP = -1
         ENDIF
      ELSEIF (IC.EQ.3) THEN
         IASTRT = NLRGBL + 1
         IBSTRT = 1
         IASMAX = MAXSHL
         IBSMAX = NLRGBL
         IF (L2TYP(4)) THEN
C           Gaunt integrals
            I2TYP = 4
         ELSE
C           No valid class of integrals for this value of IC
            IASTRT = 1
            IBSTRT = 1
            IASMAX = 0
            IBSMAX = 0
            I2TYP = -1
         ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gispck */
      SUBROUTINE GISPCK(NSPCK,ISPCK,ISHELA,ISHELBS,ISHELBE,
     &                  NBBASA,NBBASB,IBBASA,IBBASB,IBASA,IBASB)
C*****************************************************************************
C
C     Make arrays needed for symmetry packing in 4-index transformation.
C     Luuk Visscher, january 1997
C     Revised May 2001 for range in ISHELB /jth+hjaaj
C
C*****************************************************************************
#include "implicit.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
C
      DIMENSION NSPCK(0:7,0:3),ISPCK(0:7,0:7,3)
      DIMENSION NBBASA(0:7,0:2),NBBASB(0:7,0:2)
      DIMENSION IBBASA(0:7,0:2),IBBASB(0:7,0:2)
      DIMENSION IBASA(2),IBASB(2)
C
#include "ibtfun.h"
C
C     Find out to which component these shells belong
C
      ICA = 0
      IF (ISHELA.NE.0) THEN
         IF (ISHELA.LE.NLRGBL) THEN
            ICA = 1
         ELSE
            ICA = 2
         ENDIF
C
C        We generate the basis function pointer arrays.
C
         IBASI = 0
         DO IFSYM = 1, NFSYM
            IBASA(IFSYM) = IBASI
            DO IREPA = 0, NBSYM-1
               IF (IFSYM.EQ.JBTOF(IREPA,ICA)) THEN
                  NBA = NINSH(ISHELA,IREPA,INDX,0)
                  NBBASA(IREPA,0)     = NBA
                  NBBASA(IREPA,ICA)   = NBA
                  NBBASA(IREPA,3-ICA) = 0
                  IBBASA(IREPA,ICA)    = IBASI
                  IBASI = IBASI + NBA
               ENDIF
            ENDDO
         ENDDO
         IF (ICA.EQ.1) THEN
            CALL ICOPY (NBSYM,IBASI,0,IBBASA(0,3-ICA),1)
         ELSE
            CALL ICOPY (NBSYM,0,0,IBBASA(0,3-ICA),1)
         ENDIF
      ENDIF
C
      ICB = 0
      IF (ISHELBS.NE.0) THEN
         IF (ISHELBS.GT.ISHELBE)
     &      CALL QUIT('GISPCK: ISHELBS .gt. ISHELBE')
         IF (ISHELBE.LE.NLRGBL) THEN
            ICB = 1
         ELSE IF (ISHELBS.GT.NLRGBL) THEN
            ICB = 2
         ELSE
            CALL QUIT('GISPCK: ISHELB both L and S type!')
         ENDIF
C
C        We generate the basis function pointer arrays.
C
         IBASI = 0
         DO IFSYM = 1, NFSYM
            IBASB(IFSYM) = IBASI
            DO IREPB = 0, NBSYM-1
               IF (IFSYM.EQ.JBTOF(IREPB,ICB)) THEN
                  NBB = 0
                  DO ISHELB = ISHELBS,ISHELBE
                     NBB = NBB + NINSH(ISHELB,IREPB,INDX,0)
                  END DO
                  NBBASB(IREPB,0)     = NBB
                  NBBASB(IREPB,ICB)   = NBB
                  NBBASB(IREPB,3-ICB) = 0
                  IBBASB(IREPB,ICB)    = IBASI
                  IBASI = IBASI + NBB
               ENDIF
            ENDDO
         ENDDO
         IF (ICB.EQ.1) THEN
            CALL ICOPY (NBSYM,IBASI,0,IBBASB(0,3-ICB),1)
         ELSE
            CALL ICOPY (NBSYM,0,0,IBBASB(0,3-ICB),1)
         ENDIF
      ENDIF
C
C     We generate the boson symmetry packing, based on square matrices.
C
      CALL IZERO(NSPCK,32)
      CALL IZERO(ISPCK,192)
C
      DO IREPAB = 0, NBSYM-1
C
C        Dimensions for the Coulomb type integrals are
C        in IC=1 (LL), 2 (SS); the ones for the Gaunt
C        are in IC=3 (SL), they are not added to the total
C        in IC=0.
C
         DO IC = 1, 2
            DO IREPB = 0, NBSYM-1
               IREPA = IBTXOR(IREPAB,IREPB)
               ISPCK(IREPA,IREPB,IC) = NSPCK(IREPAB,0)
               IF (ISHELA.EQ.0) THEN
                  NBA = NBBAS(IREPA,IC)
               ELSE
                  NBA = NBBASA(IREPA,IC)
               ENDIF
               IF (ISHELBS.EQ.0) THEN
                  NBB = NBBAS(IREPB,IC)
               ELSE
                  NBB = NBBASB(IREPB,IC)
               ENDIF
               NBAB = NBA * NBB
               NSPCK(IREPAB,0)  = NSPCK(IREPAB,0) + NBAB
               NSPCK(IREPAB,IC) = NSPCK(IREPAB,IC) + NBAB
            ENDDO
         ENDDO
C
C        The Gaunt type integrals
C
            DO IREPB = 0, NBSYM-1
               IREPA = IBTXOR(IREPAB,IREPB)
               ISPCK(IREPA,IREPB,3) = NSPCK(IREPAB,3)
               IF (ISHELA.EQ.0) THEN
                  NBA = NBBAS(IREPA,2)
               ELSE
                  NBA = NBBASA(IREPA,2)
               ENDIF
               IF (ISHELBS.EQ.0) THEN
                  NBB = NBBAS(IREPB,1)
               ELSE
                  NBB = NBBASB(IREPB,1)
               ENDIF
               NBAB = NBA * NBB
               NSPCK(IREPAB,3) = NSPCK(IREPAB,3) + NBAB
            ENDDO
C
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ninsh */
      FUNCTION NINSH(ISHELA,IREPA,INDX,IPRINT)
C
C     Counts number of functions in a shell for a given irrep.
C     When called with IREPA = -1 it gives the total number for a shell.
C     When called with ISHELA = 0 and IREPA = -1, it loops over all
C     shells to create a pointer array from the complete list to the
C     position in the shell and the symmetry. The function value is than
C     the maximum number of functions in a shell.
C
C     INDX(1,INDA) : Position of function in shell
C     INDX(2,INDA) : Irreducible representation of function
C     INDX(3,INDA) : Position of function within this particular irrep
C
C     Luuk Visscher, 24-1-97
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
#include "dcbtra.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
C
      DIMENSION INDX(3,*),IIREP(0:7),IIIREP(0:7)
C
#include "ibtfun.h"
C
!     hjaaj: disabled print of header
!     as no other print in this routine!
!     IF (IPRINT .GT. 6) CALL HEADER('Function NINSH',-1)
C
      IF (ISHELA.NE.0) THEN
         ISHELS = ISHELA
         ISHELF = ISHELA
      ELSE
         ISHELS = 1
         ISHELF = MAXSHL
      ENDIF
C
      IF (IREPA.NE.-1) THEN
         IREPS = IREPA
         IREPF = IREPA
      ELSE
         IREPS = 0
         IREPF = MAXREP
      ENDIF
C
      MXINSH = 0
      II = 0
      DO IREPA1 = IREPS, IREPF
         IIIREP(IREPA1) = 0
      ENDDO
C
      DO ISHELA1 = ISHELS, ISHELF
C
         I = 0
         DO IREPA1 = IREPS, IREPF
            IIREP(IREPA1) = 0
         ENDDO
         NHKTA1 = NHKTSH(ISHELA1)
         KHKTA1 = KHKTSH(ISHELA1)
         MULA1  = ISTBSH(ISHELA1)
         NORBA1 = NORBSH(ISHELA1)
         NSTRA1 = IORBSB(IORBSH(ISHELA1,1))
C
C        Loop over components
C
         DO NA = 1,KHKTA1
            NSTRNA = NSTRA1 + NA
            ITYNA  = ISYMAO(NHKTA1,NA)
C
C           Loop over symmetries
C
            DO IREPA1 = IREPS, IREPF
               IF (IBTAND(MULA1,IBTXOR(IREPA1,ITYNA)) .EQ. 0) THEN
C
C                 Loop over contracted orbitals
C
                  DO IA = 1,NORBA1
                     INDA = IPTSYM(NSTRNA + KHKTA1*(IA-1),IREPA1)
                     I = I + 1
                     IIREP(IREPA1) = IIREP(IREPA1) + 1
                     IF (ISHELA.EQ.0) THEN
                       INDX(1,INDA) = I
                       INDX(2,INDA) = IREPA1
                       INDX(3,INDA) = IIREP(IREPA1)
                     ENDIF
                  ENDDO
               END IF
            ENDDO
         ENDDO
         II = II + I
         DO IREPA1 = IREPS, IREPF
            IIIREP(IREPA1) = IIIREP(IREPA1) + IIREP(IREPA1)
         ENDDO
         MXINSH = MAX(MXINSH,I)
      ENDDO
C
      IF (IREPA.NE.-1) THEN
         NINSH = IIIREP(IREPA)
      ELSE
         IF (ISHELA.NE.0) THEN
             NINSH = II
         ELSE
            IF (NTBAS(0).NE.II) THEN
               write(lupri, *)'ERROR in NINSH; ntbas(0), II',ntbas(0),II
               CALL QUIT (' Counting error in NINSH')
            END IF
            NINSH = MXINSH
         ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck ndistrn */
      FUNCTION NDISTRN(I2TYP,IPRINT,MKINAB,MKIJPA,MKIAB,ABSPOI,INDX,
     &                 INDXAB,IJPAIRS,IABIND,NWDS,IJPASS)
C
C     Calculate number of distributions in a shell pair
C     Set up pointer array INDXAB
C
C     Luuk Visscher 24-01-97
C
C     Input :
C
C     - ISHLAB : Compound label of ISHELA and ISHELB
C     - I2TYP  : Integral class
C     - IPRINT : Print Flag
C     - MKINAB : .TRUE.   Make pointer array
C                .FALSE.  Only calculate number of distributions
C     - MKIJPA : .TRUE.   Make pointer array for each pass
C                .FALSE.  Only calculate number of passes
C     - MKIAB  : .TRUE.   Make another pointer array from distribution to pair indices
C     - ABSPOI : .TRUE.   Absolute pointers in INDXAB(,,4)
C                .FALSE.  Pointers within pass in INDXAB(,,4)
C     - INDX   : Pointer from complete list to position in the shell
C     - NWDS   : Maximum number of words for symmetry packed int batch
C
C     Output :
C
C     - INDXAB :  Pointer array for symmetry packing of distributions
C     - IJPAIRS : Number of distributions for each boson symmetry
C     - NWDS :    Number of words required to store symmetry packed ints
C     - NPASS :   Number of passes for this distribution set (kept in
C                 disbuf)
C     - IJPASS  : Like IJPAIRS, but for each batch separately
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "blocks.h"
      LOGICAL SHAEQB1,DCMPAB,MKINAB,MKIJPA,MKIAB,DIAGAB1,ABSPOI
      DIMENSION INDX(3,*),INDXAB(NINSHA,NINSHB,5),IJPAIRS(0:7)
      DIMENSION IJPASS(0:7,*)
      DIMENSION IABIND(*)
      DIMENSION IND(2)
#include "twocom.h"
#include "dgroup.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "symmet.h"
#include "dcbtra.h"
C
#include "ibtfun.h"
C
!     hjaaj: disabled print of header
!     as no other print in this routine!
!     IF (IPRINT .GT. 5) CALL HEADER('Function NDISTRN',-1)
      IF (MKINAB.AND.(.NOT.ABSPOI).AND.(.NOT.MKIJPA))
     &    CALL QUIT ('Conflicting options in NDISTRN')
C
      IF (MKINAB) CALL IZERO(INDXAB,NINSHA*NINSHB*5)
      CALL IZERO(IJPAIRS,8)
      IF (MKIJPA) CALL IZERO(IJPASS,NPASS*8)
      MWDS  = 0
      MPASS = NPASS
      NPASS = 1
C
      IF (I2TYP.EQ.1) THEN
         IC = 1
      ELSEIF (I2TYP.EQ.4) THEN
         IC = 3
      ELSE
         IC = 0
      ENDIF
C
      IOFF = 0
      IJ = 0
C
      ISHELA1 = ISHLA
      ISHELB1 = ISHLB
      NHKTA1 = NHKTSH(ISHELA1)
      KHKTA1 = KHKTSH(ISHELA1)
      MULA1  = ISTBSH(ISHELA1)
      NORBA1 = NORBSH(ISHELA1)
      NSTRA1 = IORBSB(IORBSH(ISHELA1,1))
      NHKTB1 = NHKTSH(ISHELB1)
      KHKTB1 = KHKTSH(ISHELB1)
      MULB1  = ISTBSH(ISHELB1)
      NORBB1 = NORBSH(ISHELB1)
      NSTRB1 = IORBSB(IORBSH(ISHELB1,1))
      SHAEQB1 = ISHELA1 .EQ. ISHELB1
      DIAGAB1 = SHAEQB1 .AND. .NOT.BIGVEC
C
      DO IREPAB = 0, NBSYM-1
C
C     Loop over components
C
         DO NA = 1,KHKTA1
            NSTRNA = NSTRA1 + NA
            ITYNA  = ISYMAO(NHKTA1,NA)
            KHKTBB = KHKTB1
            IF (DIAGAB1) KHKTBB = NA
            DO NB = 1,KHKTBB
               NSTRNB = NSTRB1 + NB
               ITYNB  = ISYMAO(NHKTB1,NB)
               DCMPAB = SHAEQB1 .AND. NA .EQ. NB
C
C              Loop over symmetries
C
               DO IREPA = 0, MAXREP
               IF (IBTAND(MULA1,IBTXOR(IREPA,ITYNA)) .EQ. 0) THEN
                  DO IREPB = 0, MAXREP
                  IF (IBTAND(MULB1,IBTXOR(IREPB,ITYNB)) .EQ. 0) THEN
C
C                 This is the sorting step : the compound symmetry
C                 should match the one in the highest loop
C
                     IF (IREPAB.EQ.IBTXOR(IREPA,IREPB)) THEN
C
C                    Loop over contracted orbitals
C
                        DO IA = 1,NORBA1
C
C                          INDA  is the position in the complete list
C                          INDAR is the position in the shell
C
                           INDA = IPTSYM(NSTRNA + KHKTA1*(IA-1),IREPA)
                           IND(1) = INDA
                           INDAR = INDX(1,INDA)
                           NORBBB = NORBB1
                           IF (DCMPAB) NORBBB = IA
                           IF (DCMPAB.AND.IREPA.LT.IREPB) NORBBB = IA-1
                           DO IB = 1,NORBBB
                              INDB = IPTSYM(NSTRNB+KHKTB1*(IB-1),IREPB)
                              INDBR = INDX(1,INDB)
                              IND(2) = INDB
                              INDAB = INPCK(2,IND)
C
C                             Check whether we should start a new batch
C
                              IF (IOFF+NSPCK(IREPAB,IC).GT.NWDS) THEN
C
                                 IF (IOFF.EQ.0) THEN
                                    write(LUPRI,'(/A,2(/A,I10)/A,2I5)')
     &                                 'FATAL ERROR in NDISTRN',
     &                                 'maximum size of gmat',NWDS,
     &                                 'size of block',NSPCK(IREPAB,IC),
     &                                 'ic, irepab=',ic,irepab
                                    CALL QUIT(
     &                                 'Maximum size of GMAT too small')
                                 END IF
C
C                                If this happens one should either
C                                decrease the shell size or increase the
C                                total memory
C
                                 MWDS = MAX(IOFF,MWDS)
                                 NPASS = NPASS + 1
                                 IOFF = 0
                              ENDIF
C
                              IJ = IJ + 1
                              IJPAIRS(IREPAB) = IJPAIRS(IREPAB) + 1
                              IF (MKIJPA) IJPASS(IREPAB,NPASS) =
     &                                    IJPASS(IREPAB,NPASS) + 1
C
                              IF (MKINAB) THEN
                                 INDXAB(INDAR,INDBR,1) = INDAB
                                 INDXAB(INDAR,INDBR,2) = IOFF
                                 INDXAB(INDAR,INDBR,3) = IREPAB
                                 IF (ABSPOI) THEN
                                    INDXAB(INDAR,INDBR,4) =
     &                              IJPAIRS(IREPAB)
                                 ELSE
                                    INDXAB(INDAR,INDBR,4) =
     &                              IJPASS(IREPAB,NPASS)
                                 ENDIF
                                 INDXAB(INDAR,INDBR,5) = NPASS
                              ENDIF
C
                              IF (MKIAB) IABIND(IJ) = INDAB
C
                              IOFF = IOFF + NSPCK(IREPAB,IC)
                            ENDDO
                        ENDDO
                     ENDIF
                  ENDIF
                  ENDDO
               ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
      NWDS = MAX(IOFF,MWDS)
      NDISTRN = IJ
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck TRASAM */
      LOGICAL FUNCTION TRASAM(IOS11,IOS21,IOS12,IOS22,NSTR1,NSTR2)
C***********************************************************************
C
C     Check relations between indices that are to be transformed
C
C     Written by L.Visscher January 1997
C
C***********************************************************************
******
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION IOS11(*),IOS12(*),IOS21(*),IOS22(*)
      DIMENSION NSTR1(2,0:2),NSTR2(2,0:2)
C
C     Check 1-2
C
      IF (NSTR1(1,0).EQ.NSTR2(1,0)) THEN
         ISAME = 1
         DO IS = 1, NSTR1(1,0)
            IF (IOS12(IS).NE.IOS11(IS)) ISAME = 2
         ENDDO
         IF (NFSYM.GE.2) THEN
            IF (NSTR1(2,0).EQ.NSTR2(2,0)) THEN
               DO IS = 1, NSTR1(2,0)
                  IF (IOS22(IS).NE.IOS21(IS)) ISAME = 2
               ENDDO
            ELSE
               ISAME = 2
            ENDIF
         ENDIF
      ELSE
        ISAME = 2
      ENDIF
C
      TRASAM = ISAME.EQ.1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck INDSAM */
      SUBROUTINE INDSAM(IOS1,IOS2,N1,N2,IBUF,NBUF,IP,NSAM)
C***********************************************************************
C
C   Count the number of elements in IOS2 identical to those of
C   IOS1.
C
C   Written by T.Saue Feb 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IOS1(*),IOS2(*),IBUF(*)
C
      CALL IZERO(IBUF,NBUF)
      DO I = 1,N1
        IBUF(IP*IOS1(I)) = 1
      ENDDO
C
      DO I = 1,N2
        NSAM = NSAM + IBUF(IP*IOS2(I))
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DEFCL*/
      SUBROUTINE DEFCL (NCLASS)
C
C     Define integral classes
C
C     Luuk Visscher
C
C     Input :  Integral class, quaternion units
C
C     Output : Contribution to class
C
C     Classes are defined such that the complete list can be generated
C     without using permutation symmetry, except for the quaternion
C     case where we assume particle-particle symmetry.
C
C     Real groups       : classes 1-4 are used
C     Complex groups    : classes 1-8 are used
C     Quaternion groups : classes 1-16 are used
C
C     Class 1  : Real part of integrals      (i    j    |k    l   )
C     Class 2  : Real part of integrals      (i    j    |kbar lbar)
C     Class 3  : Real part of integrals      (i    jbar |k    lbar)
C     Class 4  : Real part of integrals      (i    jbar |kbar l   )
C     Class 5  : Imaginary part of integrals (i    j    |k    l   )
C     Class 6  : Imaginary part of integrals (i    j    |kbar lbar)
C     Class 7  : Imaginary part of integrals (i    jbar |k    lbar)
C     Class 8  : Imaginary part of integrals (i    jbar |kbar l   )
C     Class 9  : Real part of integrals      (i    j    |kbar l   )
C     Class 10 : Real part of integrals      (i    j    |k    lbar)
C     Class 11 : Real part of integrals      (ibar j    |k    l   )
C     Class 12 : Real part of integrals      (i    jbar |k    l   )
C     Class 13 : Imaginary part of integrals (i    j    |kbar l   )
C     Class 14 : Imaginary part of integrals (i    j    |k    lbar)
C     Class 15 : Imaginary part of integrals (ibar j    |k    l   )
C     Class 16 : Imaginary part of integrals (i    jbar |k    l   )
C
#include "implicit.h"
#include "dcbtra.h"
#include "dgroup.h"
C
      IF (NCLASS.GT.16) CALL QUIT ('Only 16 integral classes defined')
C
      CALL IZERO (ICLARR,16*NCLASS)
C
C     Real classes
C
      ICLARR(1,1,1) =  1
      ICLARR(2,2,1) = -1
      ICLARR(1,1,2) =  1
      ICLARR(2,2,2) =  1
      ICLARR(3,3,3) =  1
      ICLARR(4,4,3) = -1
      ICLARR(3,3,4) = -1
      ICLARR(4,4,4) = -1
C
      IF (NCLASS.LE.4) RETURN
C
C     Imaginary classes
C
      ICLARR(1,2,5) =  1
      ICLARR(2,1,5) =  1
      ICLARR(1,2,6) = -1
      ICLARR(2,1,6) =  1
      ICLARR(4,3,7) =  1
      ICLARR(3,4,7) =  1
      ICLARR(4,3,8) = -1
      ICLARR(3,4,8) =  1
C
      IF (NCLASS.LE.8) RETURN
C
C     Quaternion classes
C
      ICLARR(1,3,9)  = -1
      ICLARR(2,4,9)  = -1
      ICLARR(1,3,10) =  1
      ICLARR(2,4,10) = -1
      ICLARR(3,1,11) = -1
      ICLARR(4,2,11) = -1
      ICLARR(3,1,12) =  1
      ICLARR(4,2,12) = -1

      ICLARR(2,3,13) = -1
      ICLARR(1,4,13) =  1
      ICLARR(1,4,14) =  1
      ICLARR(2,3,14) =  1

      ICLARR(3,2,15) = -1
      ICLARR(4,1,15) =  1
      ICLARR(4,1,16) =  1
      ICLARR(3,2,16) =  1
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck INTCL*/
      FUNCTION INTCL (ICLASS,IQ1,IQ2)
C
C     Give contribution of the quaternion combination to the class
C
C     Luuk Visscher
C
#include "implicit.h"
#include "dcbtra.h"
C
      INTCL = ICLARR(IQ1,IQ2,ICLASS)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck QTRIAN*/
      SUBROUTINE QTRIAN (TASK,FACTOR,N,NZ,ASQ,NROW,NCOL,ATRI,NDIM)
C
C     Pack or unpack symmetric quaternion matrix to lower triangular
C     form. Add result to existing matrix multiplied by factor
C
C     Luuk Visscher
C
#include "implicit.h"
C
      CHARACTER*4 TASK
      DIMENSION ASQ(NROW,NCOL,NZ),ATRI(NDIM,NZ)
C
      IF (TASK.EQ.'PACK') THEN
         DO IZ = 1, NZ
            IJ = 0
            DO J = 1, N
               DO I = J, N
                  IJ = IJ + 1
                  ATRI(IJ,IZ) = ATRI(IJ,IZ) * FACTOR + ASQ(I,J,IZ)
               ENDDO
            ENDDO
         ENDDO
      ELSE
         DO IZ = 1, NZ
            IJ = 0
            DO J = 1, N
               DO I = J, N
                  IJ = IJ + 1
                  ASQ(I,J,IZ) = ASQ(I,J,IZ) * FACTOR + ATRI(IJ,IZ)
                  ASQ(J,I,IZ) = ASQ(J,I,IZ) * FACTOR + ATRI(IJ,IZ)
               ENDDO
            ENDDO
         ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mkdrij */
      SUBROUTINE MKDRIJ(DRIJ,SAME,IPRINT,
     &                  Q1,NDMOQ1,ICMOQ1,NSTR1,
     &                  Q2,NDMOQ2,ICMOQ2,NSTR2,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     This routine does the following:
C       1) Transform vector sets 1 and 2 to AO-basis
C       2) For given AO-indices p and q find the largest
C          product (c_1*c_2)
C       3) Reduce the corresponding matrix to block indices
C
C       Q1/Q2       - vector sets 1 and 2
C       NDMOQn(2,ifrp) : row/column dimensions of coefficient n
C       ICMOQn(ifrp)   : offsets for coefficient n
C       NSTRn(ifrp,0:2): number of vectors
C         0 - total
C         1 - electronic
C         2 - positronic
C       SAME=.TRUE. - the vector sets are identical
C
C     Written by T.Saue - Aug 19 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "blocks.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      LOGICAL SAME
      DIMENSION DRIJ(NSYMBL,NSYMBL,*),WORK(*),
     &          Q1(*),NDMOQ1(2,2),ICMOQ1(2),NSTR1(2,0:2),
     &          Q2(*),NDMOQ2(2,2),ICMOQ2(2),NSTR2(2,0:2),
     &          IV(2,2)
C
C     Initialize
C
      KFRSAV = KFREE
      NDIM   = NSYMBL*NSYMBL*NFSYM
      CALL DZERO(DRIJ,NDIM)
C
C     Make index array from AOs to blocks
C
      CALL MEMGET2('INTE','INRAO',KINRAO,NTBAS(0),WORK,KFREE,LFREE)
      CALL PTRDAO(WORK(KINRAO))
C
C     Get maximal elements for vector set 1 in AO-basis
C     =================================================
C
      NRV1 = 0
      DO IFRP = 1,NFSYM
        IV(IFRP,1) = NRV1
        NRV1       = NRV1 + NSTR1(IFRP,0)*NSYMBL
      ENDDO
      CALL MEMGET2('REAL','NRV1',KNRV1,NRV1,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KNRV1),NRV1)
      CALL GETCFB(WORK(KNRV1),WORK(KINRAO),Q1,NDMOQ1,ICMOQ1,NSTR1,
     &            IPRINT,WORK,KFREE,LFREE)
C
C     If not same as vector set 1,
C     get maximal elements for vector set 2 in AO-basis
C     =================================================
C
      IF(SAME) THEN
        KNRV2 = KNRV1
        NRV2  = NRV1
        DO IFRP = 1,NFSYM
          IV(IFRP,2) = IV(IFRP,1)
        ENDDO
      ELSE
        NRV2 = 0
        DO IFRP = 1,NFSYM
          IV(IFRP,2) = NRV2
          NRV2       = NRV2 +  NSTR2(IFRP,0)*NSYMBL
        ENDDO
        CALL MEMGET2('REAL','NRV2',KNRV2,NRV2,WORK,KFREE,LFREE)
        CALL DZERO(WORK(KNRV2),NRV2)
        CALL GETCFB(WORK(KNRV2),WORK(KINRAO),Q2,NDMOQ2,ICMOQ2,NSTR2,
     &            IPRINT,WORK,KFREE,LFREE)
      ENDIF
C
C     Construct the density screening matrix
C
      DO IREPIJ = 1,NFSYM
        DO IREPJ = 1,NFSYM
          IREPI = MOD(IREPJ+IREPIJ,2) + 1
          JOFF = IV(IREPJ,2) + KNRV2  - 1
          DO JC = 1,NSTR2(IREPJ,0)
            DO JR = 1,NSYMBL
              IOFF = IV(IREPI,1) + KNRV1 - 1
              DO IC = 1,NSTR1(IREPI,0)
                DO IR = 1,NSYMBL
                  A = WORK(JOFF+JR)*WORK(IOFF+IR)
                  DRIJ(IR,JR,IREPIJ) = MAX(DRIJ(IR,JR,IREPIJ),A)
                ENDDO
                IOFF = IOFF + NSYMBL
              ENDDO
            ENDDO
            JOFF = JOFF + NSYMBL
          ENDDO
        ENDDO
      ENDDO
C
C     Symmetrize the density screening matrix (to allow for use
C     of permutation symmetry in the integral generation step)
C
      DO IREPIJ = 1, NFSYM
         DO JR = 1, NSYMBL
            DO IR = 1, JR-1
               A = MAX(DRIJ(IR,JR,IREPIJ),DRIJ(JR,IR,IREPIJ))
               DRIJ(IR,JR,IREPIJ) = A
               DRIJ(JR,IR,IREPIJ) = A
            ENDDO
         ENDDO
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        DO I = 1,NFSYM
          WRITE(LUPRI,'(A,A3)') '*MKDRIJ: Screening matrix ',FREP(I)
          CALL OUTPUT(DRIJ(1,1,I),1,NSYMBL,1,NSYMBL,
     &                NSYMBL,NSYMBL,1,LUPRI)
        ENDDO
      ENDIF
      CALL MEMREL('MKDRIJ.AO',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck getcao */
      SUBROUTINE GETCAO(CAO,Q,IFRP,NVEC,NSEL,ISEL,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Get vectors in AO-basis
C
C     Written by T.Saue August 20 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION CAO(*),Q(*),WORK(*)
      IF(NSEL.EQ.NVEC) THEN
C
C       The whole vector set is selected
C
        IF(NBSYM.EQ.1) THEN
          CALL QEXPAND(CAO,Q(ISEL),NSEL,IFRP,NFBAS(IFRP,0))
        ELSE
          NCBF = NFBAS(IFRP,0)*NSEL*4
          CALL MEMGET2('REAL','CBF',KCBF,NCBF,WORK,KFREE,LFREE)
          CALL QEXPAND(WORK(KCBF),Q(ISEL),NSEL,IFRP,NFBAS(IFRP,0))
          CALL WTSOAO(WORK(KCBF),CAO,IFRP,NFBAS(IFRP,0),
     &                NTBAS(0),NSEL*4,IPRINT)
          CALL MEMREL('GETCAO.AO',WORK,KCBF,KCBF,KFREE,LFREE)
        ENDIF
      ELSE
C
C       Only some vectors are selected
C
        NCBF = NTBAS(0)*NSEL*4
        CALL MEMGET2('REAL','CBF',KCBF,NCBF,WORK,KFREE,LFREE)
        NDSEL = NSEL*NFBAS(IFRP,0)
        NDVEC = NVEC*NFBAS(IFRP,0)
        IF(NBSYM.EQ.1) THEN
          IOFF = ISEL
          JOFF = KCBF
          DO IZ = 1,NZ
            CALL DCOPY(NDSEL,Q(IOFF),1,WORK(JOFF),1)
            IOFF = IOFF + NDVEC
            JOFF = JOFF + NDSEL
          ENDDO
          CALL QEXPAND(CAO,WORK(KCBF),NSEL,IFRP,NFBAS(IFRP,0))
        ELSE
          IOFF = ISEL
          JOFF = 1
          DO IZ = 1,NZ
            CALL DCOPY(NDSEL,Q(IOFF),1,CAO(JOFF),1)
            IOFF = IOFF + NDVEC
            JOFF = JOFF + NDSEL
          ENDDO
C.........expand coefficients to quaternion format
          CALL QEXPAND(WORK(KCBF),CAO,NSEL,IFRP,NFBAS(IFRP,0))
C.........transform to AO-basis
          CALL WTSOAO(WORK(KCBF),CAO,IFRP,NFBAS(IFRP,0),NTBAS(0),
     &                NSEL*4,IPRINT)
        ENDIF
        CALL MEMREL('GETCAO.AO',WORK,KCBF,KCBF,KFREE,LFREE)
      ENDIF
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A,I5)')
     &   '* GETCAO: Coefficients in AO-basis. Fermion ircop: ',IFRP
            CALL PRQMAT(CAO,NFBAS(IFRP,0),NSEL,NFBAS(IFRP,0),NSEL,
     &                  4,IQDEF,LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck getcfb */
      SUBROUTINE GETCFB(QB,IND,Q,NDMOQR,ICMOQR,NSTR,IPRINT,
     &                  WORK,KFREE,LFREE)
C***********************************************************************
C
C     For vector set Q return vector set in AO-basis compressed
C     on block indices
C       Q      - coefficient
C       NDMOQR(2,ifrp) : row/column dimensions of coefficient set
C       ICMOQR(ifrp)   : offsets for coefficient set
C       NSTR(ifrp,0:2): number of vectors
C         0 - total
C         1 - electronic
C         2 - positronic
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1 = 1.0D0)
C
#include "maxorb.h"
#include "aovec.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "blocks.h"
      DIMENSION QB(*),Q(*),NDMOQR(2,2),ICMOQR(2),NSTR(2,0:2),
     &          IND(*),WORK(*)
C
      KFRSAV = KFREE
      KOFF = 1
      DO 10 IFRP = 1,NFSYM
        IF(NDMOQR(2,IFRP).EQ.0) GOTO 10
C
C       1. Transform to AO-basis
C
        NCAO  = NTBAS(0)*NSTR(IFRP,0)
        NCAOQ = NCAO*4
        CALL MEMGET2('REAL','CAO',KCAO,NCAOQ,WORK,KFREE,LFREE)
        CALL GETCAO(WORK(KCAO),Q,IFRP,NDMOQR(2,IFRP),
     &              NSTR(IFRP,0),ICMOQR(IFRP),
     &              IPRINT,WORK,KFREE,LFREE)
C
C       2. Square all elements, sum up and take square root
C
        KEND = KCAO + NCAOQ - 1
        DO I = KCAO,KEND
          WORK(I) = WORK(I)*WORK(I)
        ENDDO
        DO I = 1,3
          IOFF = KCAO + NCAO * I
          CALL DAXPY(NCAO,D1,WORK(IOFF),1,WORK(KCAO),1)
        ENDDO
        KEND = KCAO + NCAO - 1
        DO I = KCAO,KEND
          WORK(I) = SQRT(WORK(I))
        ENDDO
C
C       3. Reduce vector set to block indices..
C
        CALL GATVCC(0,NSTR(IFRP,0),WORK(KCAO),NTBAS(0),
     &              QB(KOFF),IND,NSYMBL)
        IF(IPRINT.GE.3) THEN
          WRITE(LUPRI,'(A,A3)')
     &      '* GETCFB: Vector set on block indices. Fermion ircop ',
     &      FREP(IFRP)
          CALL OUTPUT(QB(KOFF),1,NSYMBL,1,NSTR(IFRP,0),
     &                NSYMBL,NSTR(IFRP,0),1,LUPRI)
        ENDIF
        KOFF = KOFF + NSYMBL*NSTR(IFRP,0)
        CALL MEMREL('GETCFB.AO',WORK,KCAO,KCAO,KFREE,LFREE)
 10     CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck pr4scr */
      SUBROUTINE PR4SCR(KGAB,KDRIJ,WORK,KFREE,LFREE,DINTSKP,
     &                  Q1,Q2,Q3,Q4,NDMOQR,ICMOQR,NSTR,
     &                  ISAME,IPRINT)
C***********************************************************************
C
C     Prepare for 4-index screening
C     NSYMBL
C      - set in PAOVE2, called by PAOVEC
C
C     Written by T.Saue Sep 8 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "dgroup.h"
#include "blocks.h"
C
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4),NSTR(2,0:2,4),ISAME(4),
     &          Q1(*),Q2(*),Q3(*),Q4(*),WORK(*),DINTSKP(2,4)
C
C     Initialize
C
C     NSYMBL is the total number of blocks, 
C     degeneracy of centers into account
      N2GAB   = NSYMBL*NSYMBL
      NDRIJ   = N2GAB*NFSYM
      CALL DZERO(DINTSKP,8)
C
C     Get Cauchy-Schwartz integrals
C
      CALL MEMGET2('REAL','GAB',KGAB,N2GAB,WORK,KFREE,LFREE)
      IJOB   = 0
      ITYPE  = 0
      IGTYP  = 1
      MAXDF = 0
      CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDF,
     &     IPRINT,WORK(KGAB),WORK(KFREE),LFREE)
C
C     Get screening matrices
C
      CALL PR4SC1(KDRIJ,WORK,KFREE,LFREE,WORK(KGAB),
     &                  Q1,Q2,Q3,Q4,NDMOQR,ICMOQR,NSTR,
     &                  ISAME,IPRINT)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck pr4sc1 */
      SUBROUTINE PR4SC1(KDRIJ,WORK,KFREE,LFREE,GABAO,
     &                  Q1,Q2,Q3,Q4,NDMOQR,ICMOQR,NSTR,
     &                  ISAME,IPRINT)
C***********************************************************************
C
C     Prepare for 4-index screening
C
C     Written by T.Saue Sep 8 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "dgroup.h"
#include "blocks.h"
C
      LOGICAL SAME
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4),NSTR(2,0:2,4),ISAME(4),
     &          Q1(*),Q2(*),Q3(*),Q4(*),WORK(*)
      DIMENSION GABAO(*)
C
C     Initialize
C
      N2GAB   = NSYMBL*NSYMBL
      NDRIJ   = N2GAB*NFSYM
      IF(IPRINT.GE.3) THEN
        DO I = 1,NFSYM
          WRITE(LUPRI,'(A,A3)') '*MKDRIJ: GAB-matrix ',FREP(I)
          CALL OUTPUT(GABAO,1,NSYMBL,1,NSYMBL,
     &                NSYMBL,NSYMBL,1,LUPRI)
        ENDDO
      ENDIF
C
C     Get screening matrices
C
      CALL MEMGET2('REAL','DRIJ',KDRIJ,NDRIJ*2,WORK,KFREE,LFREE)
C
C     Screening matrix for electron 1 ; vectors 1 and 2 !
C
      IOFF = KDRIJ
      SAME = ISAME(1).EQ.ISAME(2)
      CALL MKDRIJ(WORK(IOFF),SAME,IPRINT,
     &                  Q1,NDMOQR(1,1,1),ICMOQR(1,1),NSTR(1,0,1),
     &                  Q2,NDMOQR(1,1,2),ICMOQR(1,2),NSTR(1,0,2),
     &                  WORK,KFREE,LFREE)
C
C     Screening matrix for electron 2 ; vectors 3 and 4 !
C     Check for sameness relations
C
      IOFF = KDRIJ + NDRIJ
      IF(ISAME(1).EQ.ISAME(3).AND.ISAME(2).EQ.ISAME(4)) THEN
        CALL DCOPY(NDRIJ,WORK(KDRIJ),1,WORK(IOFF),1)
      ELSE
        SAME = ISAME(3).EQ.ISAME(4)
        CALL MKDRIJ(WORK(IOFF),SAME,IPRINT,
     &                  Q3,NDMOQR(1,1,3),ICMOQR(1,3),NSTR(1,0,3),
     &                  Q4,NDMOQR(1,1,3),ICMOQR(1,4),NSTR(1,0,4),
     &                  WORK,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck pr4sc2 */
      SUBROUTINE PR4SC2(KDRIJ,WORK,KFREE,LFREE,GABAO,
     &                  Q1,Q2,Q3,Q4,NDMOQR,ICMOQR,NSTR,
     &                  ISAME,IPRINT)
C***********************************************************************
C
C     Prepare for 4-index screening for integrals that are to be anti-symmetrized
C     Fixes old screening bug in MP2
C
C     Written by L.Visscher 2013
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "dgroup.h"
#include "blocks.h"
C
      LOGICAL SAME
      DIMENSION NDMOQR(2,2,4),ICMOQR(2,4),NSTR(2,0:2,4),ISAME(4),
     &          Q1(*),Q2(*),Q3(*),Q4(*),WORK(*)
      DIMENSION GABAO(*)
C
C     Initialize
C
      N2GAB   = NSYMBL*NSYMBL
      NDRIJ   = N2GAB*NFSYM
      CALL MEMGET2('REAL','DRIJ',KDRIJ,NDRIJ*2,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DRIJE',KDRIJE,NDRIJ*2,WORK,KFREE,LFREE)
C
C     Get screening matrices for direct contribution
C
C     Screening matrix for electron 1 ; vectors 1 and 2 !
C
      IOFF = KDRIJ
      SAME = ISAME(1).EQ.ISAME(2)
      CALL MKDRIJ(WORK(IOFF),SAME,IPRINT,
     &                  Q1,NDMOQR(1,1,1),ICMOQR(1,1),NSTR(1,0,1),
     &                  Q2,NDMOQR(1,1,2),ICMOQR(1,2),NSTR(1,0,2),
     &                  WORK,KFREE,LFREE)
C
C     Screening matrix for electron 2 ; vectors 3 and 4 !
C     Check for sameness relations
C
      IOFF = KDRIJ + NDRIJ
      IF(ISAME(3).EQ.ISAME(1).AND.ISAME(4).EQ.ISAME(2)) THEN
        CALL DCOPY(NDRIJ,WORK(KDRIJ),1,WORK(IOFF),1)
      ELSE
        SAME = ISAME(3).EQ.ISAME(4)
        CALL MKDRIJ(WORK(IOFF),SAME,IPRINT,
     &                  Q3,NDMOQR(1,1,3),ICMOQR(1,3),NSTR(1,0,3),
     &                  Q4,NDMOQR(1,1,4),ICMOQR(1,4),NSTR(1,0,4),
     &                  WORK,KFREE,LFREE)
      ENDIF
C
C     Get screening matrices for exchange contribution
C
C
C     Screening matrix for electron 1 ; vectors 1 and 4 !
C
      IOFF = KDRIJE
      SAME = ISAME(1).EQ.ISAME(4)
      CALL MKDRIJ(WORK(IOFF),SAME,IPRINT,
     &                  Q1,NDMOQR(1,1,1),ICMOQR(1,1),NSTR(1,0,1),
     &                  Q4,NDMOQR(1,1,4),ICMOQR(1,4),NSTR(1,0,4),
     &                  WORK,KFREE,LFREE)
C
C     Screening matrix for electron 2 ; vectors 3 and 2 !
C     Check for sameness relations
C
      IOFF = KDRIJE + NDRIJ
      IF(ISAME(3).EQ.ISAME(1).AND.ISAME(2).EQ.ISAME(4)) THEN
        CALL DCOPY(NDRIJ,WORK(KDRIJE),1,WORK(IOFF),1)
      ELSE
        SAME = ISAME(3).EQ.ISAME(2)
        CALL MKDRIJ(WORK(IOFF),SAME,IPRINT,
     &                  Q3,NDMOQR(1,1,3),ICMOQR(1,3),NSTR(1,0,3),
     &                  Q2,NDMOQR(1,1,2),ICMOQR(1,2),NSTR(1,0,2),
     &                  WORK,KFREE,LFREE)
      ENDIF
C
C     Take maximum of the direct and exchange contributions
C
      DO I = 1, 2 * NDRIJ
         WORK(KDRIJ+I-1) = MAX(WORK(KDRIJ+I-1),WORK(KDRIJE+I-1))
      ENDDO
C
C     Clean up
C
      CALL MEMREL('PR4SC2',WORK,KDRIJE,KDRIJE,KFREE,LFREE)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck St4scr */
      SUBROUTINE ST4SCR(DINTSKP)
C***********************************************************************
C
C     Screening statistics for 4-index transformations
C
C     Written by T.Saue Sep 8 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (DC = 100.0D0,D1 = 1.0D0,D0 = 0.0D0)
C
      DIMENSION DINTSKP(2,4)
C
      DDEL = D0
      DTOT = D0
      WRITE(LUPRI,'(A)') '* Screening statistics:'
      IF(DINTSKP(1,1).GT.D0) THEN
        DDEL  = DDEL + DINTSKP(2,1)
        DTOT  = DTOT + DINTSKP(1,1)
        DSTEP = DC*(DINTSKP(2,1)/DINTSKP(1,1))
        WRITE(LUPRI,500) '(LL|LL)ints : ',DSTEP,'%'
      ENDIF
      IF(DINTSKP(1,2).GT.D0) THEN
        DDEL  = DDEL + DINTSKP(2,2)
        DTOT  = DTOT + DINTSKP(1,2)
        DSTEP = DC*(DINTSKP(2,2)/DINTSKP(1,2))
        WRITE(LUPRI,500) '(SS|LL)ints : ',DSTEP,'%'
      ENDIF
      IF(DINTSKP(1,3).GT.D0) THEN
        DDEL  = DDEL + DINTSKP(2,3)
        DTOT  = DTOT + DINTSKP(1,3)
        DSTEP = DC*(DINTSKP(2,3)/DINTSKP(1,3))
        WRITE(LUPRI,500) '(SS|SS)ints : ',DSTEP,'%'
      ENDIF
      IF(DINTSKP(1,4).GT.D0) THEN
        DDEL  = DDEL + DINTSKP(2,4)
        DTOT  = DTOT + DINTSKP(1,4)
        DSTEP = DC*(DINTSKP(2,4)/DINTSKP(1,4))
        WRITE(LUPRI,500) '(SL|SL)ints : ',DSTEP,'%'
      ENDIF
      IF(DTOT.GT.D0) THEN
        DSTEP = DC*(DDEL/DTOT)
        WRITE(LUPRI,500) 'Total       : ',DSTEP,'%'
      ENDIF
  500 FORMAT(3X,A,F6.2,A1)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck pck2in */
      SUBROUTINE PCK2IN(NSTR,TRIAN,IPRINT)
C***********************************************************************
C
C     Define the packing of the 2-index transformed integrals
C
C     Written by T.Saue Sep 16 , based on existing code by Luuk
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbtra.h"
#include "dgroup.h"
      LOGICAL TRIAN(2)
      DIMENSION NSTR(2,0:2,4)
C
      NFPCK12T = 0
      NFPCK34T = 0
      DO IREPIJ =1, NFSYM
         NFPCK12(IREPIJ) = 0
         NFPCK34(IREPIJ) = 0
         DO IREPJ = 1, NFSYM
            IREPI = MOD(IREPJ+IREPIJ,2) + 1
C
C           Packing can be triangular (i=>j) or square.
C
            IF (.NOT.TRIAN(1).OR.(TRIAN(1).AND.IREPI.GT.IREPJ)) THEN
               IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
               NFPCK12(IREPIJ)      = NFPCK12(IREPIJ) +
     &                        NSTR(IREPI,0,1)*NSTR(IREPJ,0,2)
            ELSEIF (TRIAN(1).AND.IREPI.EQ.IREPJ) THEN
               IFPCK12(IREPI,IREPJ) = NFPCK12(IREPIJ)
               NFPCK12(IREPIJ)      = NFPCK12(IREPIJ) +
     &                        NSTR(IREPI,0,1)*(NSTR(IREPI,0,1)+1)/2
            ELSE
               IFPCK12(IREPI,IREPJ) = - 1
            ENDIF
            IF (.NOT.TRIAN(2).OR.(TRIAN(2).AND.IREPI.GT.IREPJ)) THEN
               IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
               NFPCK34(IREPIJ)      = NFPCK34(IREPIJ) +
     &                        NSTR(IREPI,0,3)*NSTR(IREPJ,0,4)
            ELSEIF (TRIAN(2).AND.IREPI.EQ.IREPJ) THEN
               IFPCK34(IREPI,IREPJ) = NFPCK34(IREPIJ)
               NFPCK34(IREPIJ)      = NFPCK34(IREPIJ) +
     &                        NSTR(IREPI,0,3)*(NSTR(IREPI,0,3)+1)/2
            ELSE
               IFPCK34(IREPI,IREPJ) = - 1
            ENDIF
         ENDDO
         NFPCK12T = NFPCK12T + NFPCK12(IREPIJ)
         NFPCK34T = NFPCK34T + NFPCK34(IREPIJ)
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         WRITE(LUPRI,*)
         WRITE(LUPRI,'(A,2I8)') 'NFPCK12',(NFPCK12(I), I=1,NFSYM)
         WRITE(LUPRI,'(A,2I8)') 'NFPCK34',(NFPCK34(I), I=1,NFSYM)
         WRITE(LUPRI,'(A,4I8)') 'IFPCK12',
     &        ((IFPCK12(I,J), I=1,NFSYM), J=1,NFSYM)
         WRITE(LUPRI,'(A,4I8)') 'IFPCK34',
     &        ((IFPCK34(I,J), I=1,NFSYM), J=1,NFSYM)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck pck3in */
      SUBROUTINE PCK3IN(ITQMT,NTQMT,NSTR,IPRINT)
C***********************************************************************
C
C     Define the packing of the 3-index transformed integrals
C
C     Written by T.Saue Sep 16 , based on existing code by Luuk
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbtra.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION ITQMT(0:7,0:7,2),NSTR(2,0:2,4)
C
#include "ibtfun.h"
C
C     Initialize
C
      NTQMT = 0
      CALL IZERO(ITQMT,128)
C
      DO IREPRS = 0, NBSYM-1
        IREPIJ = JBTOF(IREPRS,1)             ! Get parity of rs
        NIJ    = NFPCK12(IREPIJ)             ! Number of ij pairs
        DO IREPS = 0, NBSYM-1
          DO 11 IC = 1,2
            IREPR = IBTXOR(IREPRS,IREPS)
            NR = NBBAS(IREPR,IC)             ! Number of r functs
            ITQMT(IREPRS,IREPS,IC) = NTQMT   ! Pointers
            IREPL = JBTOF(IREPS,IC)          ! Parity that s contributes to
            NS = NBBAS(IREPS,IC)             ! Number of s functs
            IF (IREPIJ.EQ.1) THEN
              IREPK=IREPL
            ELSE
              IREPK=3-IREPL
            ENDIF
C              WRITE(LUPRI,'(A,I10,10I6)')
C    &              'NTQMT',NTQMT,NZ,NS,NIJ,NSTR(IREPK,0,3),IREPK,
C    &              IREPIJ,IREPS,IREPRS,NR,IREPR
            IF (NR.EQ.0) GOTO 11
            NTQMT = NTQMT + NIJ*NSTR(IREPK,0,3)*NS*NZ*NZ
11        CONTINUE
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck vircan */
      SUBROUTINE VIRCAN(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,NDMOQR,
     &                  ICMOQR,CMO)
C***********************************************************************
C
C     Check for canonical virtuals; if not, recanonize
C
C     Written by T.Saue Sep 30 1998
C
C***********************************************************************
      use interface_to_mpi
#include "implicit.h"
#include "priunit.h"
#include "infpar.h"
#include "dummy.h"
#include "aovec.h"
#include "maxorb.h"
C
#include "dcbgen.h"
#include "dcbtra.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"
#include "blocks.h"
      LOGICAL RECAN
      DIMENSION NDMOQR(2,2),ICMOQR(2),KQ(2),KE(2),KIBE(2),
     &          CMO(*),WORK(*)
      DIMENSION NFMO(2,2),IFMO(2)
C
      CALL QENTER('VIRCAN')
      KFRSAV = KFREE
C
C     Scan through active virtual eigenvalues to see if
C     recanonization is necessary
C
      RECAN = .FALSE.
      DO IFRP = 1,NFSYM
        DO J = 0,NDMOQR(2,IFRP)-1
          RECAN = RECAN.OR.(WORK(KE(IFRP)+J).EQ.DUMMY)
        ENDDO
      ENDDO
C
      RECAN = RECAN.OR.RCORBS
C
C     Recanonization
C     ==============
C     Note that closed shells are assumed !
C
      IF(RECAN) THEN
        CALL HEADER('Recanonization of virtual orbitals',-1)
C
C       Get hold of the slaves if we run in parallel.
C        ( ITASK = 1 for Fock matrices )
C
        IF (PARCAL) CALL DIRAC_PARCTL( HERFCK_PAR )
C
        NELM = 0
        DO IFRP = 1,NFSYM
          IFMO(IFRP) = NELM
          NFMO(1,IFRP) = NDMOQR(2,IFRP)
          NFMO(2,IFRP) = NDMOQR(2,IFRP)
          NELM = NELM + NDMOQR(2,IFRP)*NDMOQR(2,IFRP)*NZ
        ENDDO
        CALL MEMGET2('REAL','FMO',KFMO,NELM    ,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','FAO',KFAO,N2BBASXQ,WORK,KFREE,LFREE)
        CALL VIRCA1(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,NDMOQR,ICMOQR,
     &       INTDEF,CMO,WORK(KFMO),NFMO,IFMO,WORK(KFAO))
        CALL MEMREL('VIRCAN',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C       Release the slaves if we run in parallel.
C        ( ITASK = -1 )
C
        IF (PARCAL) CALL DIRAC_PARCTL( EXIT_NODEMENU )
      END IF
C
      CALL QEXIT('VIRCAN')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck virca1 */
      SUBROUTINE VIRCA1(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,NDMOQR,
     &              ICMOQR,INTDHF,CMO,FMO,NFMO,IFMO,FAO)
C***********************************************************************
C
C     Recanonize virtual orbitals...
C     NB!! For the moment this routine assumes closed shell !!!!
C     Written by T.Saue Jan 11 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0=0.0D0)
C
#include "dcbgen.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbham.h"
      DIMENSION NDMOQR(2,2),ICMOQR(2),KQ(2),KE(2),KIBE(2),
     &          CMO(*),FAO(*),FMO(*),WORK(*)
      DIMENSION NFMO(2,2),IFMO(2)
C
C     Get Fock-matrix in active orbitals
C     ==================================
C
      CALL GTVFCK(FMO,NFMO,IFMO,INTDHF,CMO,FAO,
     &            WORK(KQ(1)),NDMOQR,ICMOQR,WORK(KQ(1)),NDMOQR,ICMOQR,
     &            IPRINT,WORK,KFREE,LFREE)
C
C     Diagonalize Fock-matrix
C     =======================
C
      DO 30 I = 1,NFSYM
        IF(NFMO(1,I).EQ.0) GOTO 30
        IF (.NOT.SPINFR) THEN
           CALL QDIAG(NZ,NFMO(1,I),FMO(IFMO(I)+1),
     &             NFMO(1,I),NFMO(2,I),
     &             WORK(KE(I)),1,
     &             FAO(IFMO(I)+1),NFMO(1,I),NFMO(2,I),
     &             WORK(KFREE),LFREE,IERR)
        ELSE
           NFMO2 = NFMO(1,I)*NFMO(1,I)
           CALL DZERO(FAO(1+IFMO(I)),NFMO2*NZ)
           NBO = NFMO(1,I)
           CALL RSJACO(NBO,NBO,NBO,FMO(IFMO(I)+1),
     &                 WORK(KE(I)),1,1,0,FAO(IFMO(I)+1))
           CALL MEMGET2('INTE','IBEIG2',KIBEIG2,NBO,WORK,KFREE,LFREE)
           CALL IDBOS (NBO,WORK(KIBE(I)),WORK(KIBEIG2),FAO(IFMO(I)+1))
           CALL MEMREL('VIRCA1',WORK,KIBEIG2,KIBEIG2,KFREE,LFREE)
        ENDIF
   30 CONTINUE
C
C     Backtransform coefficients
C     ==========================
C
      DO 40 I = 1,NFSYM
        IF(NFMO(1,I).EQ.0) GOTO 40
        CALL BCKTRA(CMO(ICMOQR(I)),NDMOQR(1,I),NDMOQR(2,I),
     &              FAO(IFMO(I)+1),NDMOQR(2,I),NDMOQR(2,I),
     &              NDMOQR(2,I),NZ,
     &              NDMOQR(2,I),1,NDMOQR(1,I),
     &              WORK(KQ(I)),NDMOQR(1,I),NDMOQR(2,I),NZ,
     &              IPRINT)
        NELM = NDMOQR(1,I)*NDMOQR(2,I)*NZ
        CALL DCOPY(NELM,CMO(ICMOQR(I)),1,WORK(KQ(I)),1)
 40   CONTINUE
C
C     Print section
C
      IF(IPRINT.GE.1) THEN
        WRITE(LUPRI,'(A)')
     &   '* Eigenvalues of recanonized orbitals:'
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A3)') '* Fermion ircop ',FREP(IFRP)
          WRITE(LUPRI,'(5D15.8)')
     &             (WORK(KE(IFRP)+IXX), IXX=0,(NFMO(1,IFRP)-1))
        ENDDO
        IF(IPRINT.GE.5) THEN
        WRITE(LUPRI,'(A)')
     &   '* Coefficients of recanonized orbitals:'
          DO I = 1,NFSYM
            WRITE(LUPRI,'(A,A3)') '* Fermion ircop ',FREP(I)
            CALL PRQMAT(WORK(KQ(I)),NDMOQR(1,I),NDMOQR(2,I),
     &                  NDMOQR(1,I),NDMOQR(2,I),NZ,IPQTOQ(1,0),LUPRI)
          ENDDO
        ENDIF
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck eigtst */
      SUBROUTINE EIGTST(EIG,NVEC,IVEC,VAL,TESTOK)
C***********************************************************************
C
C     Written by T.Saue Sep 29 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL TESTOK
      DIMENSION EIG(*),IVEC(*)
C
      DO I = 1,NVEC
        TESTOK = TESTOK.AND.(EIG(IVEC(I)).NE.VAL)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck eigthr */
      SUBROUTINE EIGTHR(EIG,NVEC,NSEL,IVEC,IBUF,THRESH)
C***********************************************************************
C
C     Written by T.Saue Sep 29 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION EIG(*),IVEC(*),IBUF(*)
C
      DO I = 1,NVEC
      IF(EIG(IVEC(I)).LE.THRESH) THEN
        NSEL       = NSEL + 1
        IBUF(NSEL) = IVEC(I)
      ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck virthr */
      SUBROUTINE VIRTHR(WORK,KFREE,LFREE,LUCOEF,NSTR,KVEC,
     &                  THRESH)
C***********************************************************************
C
C     Truncate virtual space using eigenvalue threshold THRESH
C
C     Written by T.Saue Sep 30 1998
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
      DIMENSION WORK(*),NSTR(2,0:2),KVEC(*)
      KFRSAV = KFREE
C
C     Read and check eigenvalues
C
      CALL MEMGET2('REAL','EIG',KEIG,NORBT,WORK,KFREE,LFREE)
      CALL REACMO(LUCOEF,'DFCOEF',DUM,WORK(KEIG),DUM,DUM,4)
      DO IFRP = 1,NFSYM
        NSEL = 0
        NDEL = 0
        CALL MEMGET2('INTE','IBUF',KBUF,NSTR(IFRP,1),WORK,KFREE,LFREE)
        IOFF = NPSH(IFRP) + IORB(IFRP)
        CALL EIGTHR(WORK(KEIG+IOFF),NSTR(IFRP,1),NSEL,
     &              WORK(KVEC(IFRP)),WORK(KBUF),THRESH)
        IF(NSEL.LT.NSTR(IFRP,1)) THEN
          NDEL = NSTR(IFRP,1) - NSEL
          NSTR(IFRP,1) = NSEL
          NSTR(IFRP,0) = NSTR(IFRP,1) + NSTR(IFRP,2)
          CALL ICOPY(NSEL,WORK(KBUF),1,WORK(KVEC(IFRP)),1)
        ENDIF
        WRITE(LUPRI,'(A,A3,A,I5)')
     &      '* Fermion ircop ',FREP(IFRP),'- Virtuals deleted: ',NDEL
      ENDDO
      CALL MEMREL('VIRTHR',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GTVFCK */
      SUBROUTINE GTVFCK(FMO,NFMO,IFMO,INTDHF,CMO,FAO,
     &                  Q1,NCMO1,ICMO1,Q2,NCMO2,ICMO2,
     &                  IPRINT,WORK,KFREE,LFREE)
C*****************************************************************************
C
C     Get a block of the MO Fock-matrix; CMO is destroyed in the process.
C
C     Written by T.Saue, Jan 12 1999
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "aovec.h"
#include "maxorb.h"
      PARAMETER(D1 = 1.0D0,D0 = 0.0D0)
C
#include "dcbgen.h"
#include "dcbfir.h"
#include "dcbdhf.h"
#include "blocks.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "cbihr2.h"
CTOBS: Note the use of dcbdhf ! Should be generalized !!!
      LOGICAL TOBE, SAVEFLAGS(4)
      DIMENSION FMO(*),NFMO(2,2),IFMO(2),FAO(*),CMO(*),Q1(*),Q2(*),
     &          NCMO1(2,2),NCMO2(2,2),ICMO1(2),ICMO2(2),
     &          WORK(*)
C
      CALL QENTER('GTVFCK')
      KFRSAV = KFREE
      LUBUF = 22
C
C     1. Get two-electron Fock matrix in AO-basis
C     ===========================================
C
      INQUIRE(FILE='DFFCK2',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK2','OLD','GTVFCK')
        CALL REAFCK(LUBUF,CMO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GETFCK: No 2-Fock found. Regenerating.'

        call SaveTaskDistribFlags(saveflags)
       call SetTaskDistribFlags((/ .TRUE. , .TRUE. , .TRUE. , .TRUE. /))
        call SetIntTaskArrayDimension(NPOS,PARCAL)
        if (NPOS.GT.0) THEN
           CALL MEMGET2('INTE','IPOS',KPOS,NPOS,WORK,KFREE,LFREE)
        else
           KPOS = KFREE
        endif

C
C       Closed shell Hartree - Fock
C       ===========================
CTROND: This must be generalized !!! TODO
        NFMAT = 1
C       Totally symmetric operator
        ISYMOP(1) = 1
C       Fock matrix type
        IFCKOP(1) = 1
C       Hermitian operator
        IHRMOP(1) = 1
        INQUIRE(FILE='DFDENS',EXIST=TOBE)
        IF(TOBE) THEN
          CALL OPNFIL(LUBUF,'DFDENS','OLD','GTVFCK')
          CALL READNS(LUBUF,FAO)
          CLOSE(LUBUF,STATUS='KEEP')
        ELSE
          CALL DENSTY(FAO,CMO,IPRINT)
        ENDIF
        CALL TWOFCK(ISYMOP,IHRMOP,IFCKOP,CMO,FAO,NFMAT,
     &       WORK(KPOS),INTDHF,IPRINT,WORK(KFREE),LFREE)

        IF(PARCAL) call SetTaskDistribFlags(saveflags)

      ENDIF
C
C     3a. Solvent: get solvent contribution
C     =====================================
C
      IF(SOLVEN) THEN
C        CALL DENSTY(FAO,CMO,IPRINT)
        NF = 1
        CALL SOLFCK(CMO,FAO,NF,ESOLVE,ESOLVN,WORK(KFREE),LFREE,IPRSOL)
      ENDIF
C
C     3. Get one-electron Fock matrix
C     ===============================
C
      INQUIRE(FILE='DFFCK1',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUBUF,'DFFCK1','OLD','GTVFCK')
        CALL REAFCK(LUBUF,FAO,.TRUE.,1)
        CLOSE(LUBUF,STATUS='KEEP')
      ELSE
        WRITE(LUPRI,'(A)')
     &    '*** WARNING *** GETFCK: No 1-Fock found. Regenerating.'
        CALL ONEFCK(FAO,IPRINT,WORK(KFREE),LFREE)
      ENDIF
C
C     4. Add 1- and 2-electron matrix and transform to MO-basis
C     =========================================================
C
      CALL DAXPY(N2BBASXQ,D1,CMO,1,FAO,1)
C
C     Make the matrix spinfree if desired.
C
      IF (SPINFR) CALL SPFAO (CMO,FAO,WORK(KFREE),LFREE)
C
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTVFCK:Total Fock matrix in AO-basis',-1)
        CALL PRQMAT(FAO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &              IPQTOQ(1,0),LUPRI)
      ENDIF

      DO 10 I = 1,NFSYM
        N2FMO = NFMO(1,I)*NFMO(2,I)
        IF(N2FMO.EQ.0) GOTO 10
        CALL QTRANS('AOMO','S',D0,
     &              NFBAS(I,0),NFBAS(I,0),NFMO(1,I),NFMO(2,I),
     &              FAO(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &              FMO(IFMO(I)+1),NFMO(1,I),NFMO(2,I),NZ,IPQTOQ(1,0),
     &              Q1(ICMO1(I)),NCMO1(1,I),NCMO1(2,I),NZ,IPQTOQ(1,0),
     &              Q2(ICMO2(I)),NCMO2(1,I),NCMO2(2,I),NZ,IPQTOQ(1,0),
     &              WORK(KFREE),LFREE,IPRINT)
      IF(IPRINT.GE.5) THEN
        CALL HEADER('GTVFC1:MO Fock matrix '//FREP(I),-1)
        CALL PRQMAT(FMO(IFMO(I)+1),NFMO(1,I),NFMO(2,I),
     &              NFMO(1,I),NFMO(2,I),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
 10   CONTINUE
      CALL MEMREL('GTVFCK',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      CALL QEXIT('GTVFCK')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck RESOLV */
      SUBROUTINE RESOLV()
C***********************************************************************
C
C     Driver routine for 4-index transformation for
C     resolution of electronic states within open-shell manifold
C
C     Based on PAMTRA
C
C     Written by T.Saue Jan 18 1999
C     Modified by Miro for DFT-COSCI, August 2016
C
C***********************************************************************

      Use memory_allocator
      Use moltra_labeling

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "infpar.h"
C
#include "maxorb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbres.h"
#include "dcbgen.h"
#include "dcbdhf.h"
      CHARACTER DAYTID*24,SECTID*12,CPUTID*12
      LOGICAL TOBE,TRANOP,SKPONE,SKPTWO
      DIMENSION NSTR(2,0:2,4),NSTRT(2),NSPC(2,0:2),NSPC2(2,0:3),
     &          KVEC(2,2),KQ(2,4),KE(2,4),KIBE(2,4),NQ(2),KQC(2)
CTROND: to fix !
      DIMENSION MINE(16),MAXE(16),NGASO(256)
      real(8), allocatable :: WORK(:)
C     Initialize (for debug mode)
      KE=0
      KIBE=0
      KQC=0
C
      CALL QENTER('RESOLV')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in RESOLV')
      KFRSAV = KFREE
      CALL TITLER('Resolution of open-shell states','*',120)
C
C     Determine number of active electrons
C
      IF (NAELEC.EQ.0) THEN
         WRITE(LUPRI,'(A)')
     &   '*WARNING: No open-shell electrons, no states to resolve!'
         GOTO 9999
      END IF
      IND = 1
      NGAS_DC = NOPEN
      NSYMRP = 2*NFSYM
      DO IOPEN = 1,NOPEN
        ITEMP = 0
        DO IFRP = 1,NFSYM
          NGASO(IND) = NACSH(IFRP,IOPEN)
          NGASO(IND+1) = NGASO(IND)
          IND = IND + 2
          ITEMP = ITEMP + NACSH(IFRP,IOPEN)      
        ENDDO
        MINE(IOPEN) = ITEMP
        MAXE(IOPEN) = MINE(IOPEN)
      ENDDO
C
C     Check for necessary files:
C       MRCONEE - one-electron integrals and other information
C       MDCINT  - two-electron integrals
C
      INQUIRE(FILE='MRCONEE',EXIST=SKPONE)
      IF(SKPONE) THEN
        WRITE(LUPRI,'(A)')
     & '*WARNING: MRCONEE found. Skipping 2-index transformation.'
      ENDIF
      INQUIRE(FILE='MDCINT',EXIST=SKPTWO)
      IF(SKPTWO) THEN
        WRITE(LUPRI,'(A)')
     & '*WARNING: MDCINT found. Skipping 4-index transformation.'
      ENDIF
      IF(SKPONE.AND.SKPTWO) GOTO 10
C
C     Transfer information from DCBRES to DCBTRA
C     ==========================================
C
      TRANOP = NOPAIR
      NOPAIR = .TRUE.
      SCRBUF = SCRTRA
      SCRTRA = SCRRES
      IPRBUF = IPRTRA
      IPRTRA = IPRRES
      ISTBUF = ISTRAT
      ISTRAT = ISTRES
      THRBUF = THROUT
      THROUT = RESOUT
C     number of double quaternionic classes for strategy 3 and 4
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Set up index arrays for active orbitals;
C     =======================================
C
      NQT  = 0
      ISAME(1) = 1
      NSTRT(1) = 0
      DO IFRP = 1,NFSYM
        CALL MEMGET2('INTE','VEC',KVEC(IFRP,1),NASH(IFRP),
     &     WORK,KFREE,LFREE)
        CALL INDARR(WORK(KVEC(IFRP,1)),NISH(IFRP),NASH(IFRP),IPRINT)
        NSTR(IFRP,1,1) = NASH(IFRP)
        NSTR(IFRP,2,1) = 0
        NSTR(IFRP,0,1) = NSTR(IFRP,1,1) + NSTR(IFRP,2,1)
        NDMOQR(1,IFRP,1) = NFBAS(IFRP,0)
        NDMOQR(2,IFRP,1) = NASH(IFRP)
        ICMOQR(IFRP,1) = NQT + 1
        NSTRT(1) = NSTRT(1) +  NSTR(IFRP,0,1)
C
        NQ(IFRP) = NFBAS(IFRP,0)*NASH(IFRP)*NZ
        NQT      = NQT + NQ(IFRP)
        CALL MEMGET2('INTE','IBEIG',KIBE(IFRP,1),NSTR(IFRP,0,1),
     &     WORK,KFREE,LFREE)
      ENDDO
      CALL MEMGET2('REAL','KQ',KQ(1,1),NQT,WORK,KFREE,LFREE)
      KQ(2,1) = KQ(1,1) + NQ(1)
      CALL MEMGET2('REAL','KE',KE(1,1),NSTRT(1),WORK,KFREE,LFREE)
      KE(2,1) = KE(1,1) + NSTR(1,0,1)
C
C     Copy for indices 2..4
C
      DO I = 2,4
        ISAME(I) = 1
        DO IFRP = 1,NFSYM
          KQ(IFRP,I)     = KQ(IFRP,1)
          KE(IFRP,I)     = KE(IFRP,1)
          KIBE(IFRP,I)   = KIBE(IFRP,1)
          NSTR(IFRP,0,I) = NSTR(IFRP,0,1)
          NSTR(IFRP,1,I) = NSTR(IFRP,1,1)
          NSTR(IFRP,2,I) = NSTR(IFRP,2,1)
          NDMOQR(1,IFRP,I) = NDMOQR(1,IFRP,1)
          NDMOQR(2,IFRP,I) = NDMOQR(2,IFRP,1)
          ICMOQR(IFRP,I)   = ICMOQR(IFRP,1)
        ENDDO
      ENDDO
C
      CALL Make_Kramer_to_SpinorIndex (NFSYM,NSTR)

      KTEMP = KFREE
      IF(.NOT.SKPONE) THEN
C
C       Set up index arrays for core orbitals
C       =======================================
        NQT  = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET2('INTE','KVEC',KVEC(IFRP,2),NISH(IFRP),
     &       WORK,KFREE,LFREE)
          CALL INDARR(WORK(KVEC(IFRP,2)),0,NISH(IFRP),IPRINT)
          NSPC(IFRP,1) = NISH(IFRP)
          NSPC(IFRP,2) = 0
          NSPC(IFRP,0) = NSPC(IFRP,1) + NSPC(IFRP,2)
          NDMOQC(1,IFRP,1) = NFBAS(IFRP,0)
          NDMOQC(2,IFRP,1) = NISH(IFRP)
          ICMOQC(IFRP,1) = NQT + 1
          NQ(IFRP) = NFBAS(IFRP,0)*NISH(IFRP)*NZ
          NQT      = NQT + NQ(IFRP)
          KE(IFRP,2) = KE(IFRP,1)
        ENDDO
        CALL MEMGET2('REAL','KQC',KQC(1),NQT,WORK,KFREE,LFREE)
        KQC(2) = KQC(1) + NQ(1)
        NSTRT(2) = NSTRT(1)
      ENDIF
C
      WRITE(LUPRI,'(A,I5)') '- Number of active electrons:',NAELEC
      WRITE(LUPRI,'(A)')    '- Active orbitals:'
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,1)),0,0,0,NSTR)
      ENDDO
      WRITE(LUPRI,'(A)')    '- Core orbitals:'
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,2)),0,0,0,NSPC)
      ENDDO
C
C     Allocate buffer spaced for coefficients AND core Fock matrix
C
!Miro: extra space (KBUF1) for full set of MO coeff for DFT-COSCI
      CALL MEMGET2('REAL','KBUF',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','KBUF1',KBUF1,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','KEIG',KEIG,NORBT ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','KIBEIG',KIBEIG,NORBT,WORK,KFREE,LFREE)
C
C     Get all coefficients and eigenvalues
C     ====================================
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
         WRITE(LUPRI,'(/A)')
     &      'RESOLV ERROR: MO coefficient file not present'
         CALL QUIT('RESOLV: No MO coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RESOLV')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KBUF),WORK(KEIG),WORK(KIBEIG),
     &               TOTERG,14)
        CALL REACMO(LUCOEF,'DFCOEF',WORK(KBUF1),WORK(KEIG),WORK(KIBEIG),
     &               TOTERG,14)
         CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF
C
C     Select active orbitals
C     ======================
C
      DO IFRP = 1, NFSYM
        CALL SELCFS (WORK(KBUF+ICMOQ(IFRP)),
     &       IFRP,WORK(KQ(IFRP,1)),NSTR(IFRP,0,1),WORK(KVEC(IFRP,1)),
     &       NSTR(IFRP,2,1),NSTR(IFRP,1,1),
     &       NFBAS(IFRP,0),NORB(IFRP))
        CALL SELEIG (WORK(KEIG+IORB(IFRP)),IFRP,WORK(KE(IFRP,1)),
     &               WORK(KVEC(IFRP,1)),NSTR(IFRP,2,1),NSTR(IFRP,1,1))
        CALL SELIBEIG(WORK(KIBEIG),IORB(IFRP),IFRP,
     &                WORK(KIBE(IFRP,1)),WORK(KVEC(IFRP,1)),
     &                NSTR(IFRP,2,1),NSTR(IFRP,1,1))
      ENDDO
      IF(.NOT.SKPONE) THEN
C
C       Select core orbitals
C       ======================
C
        DO IFRP = 1, NFSYM
          CALL SELCFS (WORK(KBUF+ICMOQ(IFRP)),
     &         IFRP,WORK(KQC(IFRP)),NSPC(IFRP,0),WORK(KVEC(IFRP,2)),
     &         NSPC(IFRP,2),NSPC(IFRP,1),
     &         NFBAS(IFRP,0),NORB(IFRP))
        ENDDO
C
C       2-index transformation
C       ======================
C
        CALL TRAHI(IPRRES,2)
C
C       Get hold of the slaves if we run in parallel.
C        ( ITASK = 1 for Fock matrices )
C
        IF (PARCAL) CALL DIRAC_PARCTL( HERFCK_PAR )
C
C       don't use frozen open shells
        NSPC2=-1
C       Transform core-fock matrix to active MO basis
C       ---------------------------------------------
C Miro/Aug 2016 - introduced DFT-COSCI, core-fock mtx contains xc contrib 

        CALL GTCFCK(WORK(KBUF),IPRRES,WORK(KQ(1,1)),WORK(KQ(1,2)),
     &            WORK(KQC(1)),WORK(KBUF1),NDMOQR,ICMOQR,NSTR(1,0,1),
     &             NDMOQC,ICMOQC,NSPC(1,0),NSPC2,
     &             INTRES,WORK,KFREE,LFREE)
C 
C       Write MOLFDIR file MRCONEE
C       --------------------------
C
        CALL MEMGET2('REAL','KORBE',KORBE,NSTRT(1),WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','KIBORB',KIBORB,NSTRT(1),WORK,KFREE,LFREE)
C
        CALL MRCONE(WORK(KBUF),WORK(KORBE),WORK(KE(1,1)),WORK(KE(2,1)),
     &              WORK(KIBORB),WORK(KIBE(1,1)),WORK(KIBE(2,1)),
     &              NSPC,NSTR,NSTRT,IPRRES,TOTERG,
     &              WORK,KFREE,LFREE)
C
C       Release the slaves if we run in parallel.
C        ( ITASK = -1 )
C
        IF (PARCAL) CALL DIRAC_PARCTL( EXIT_NODEMENU )
      ENDIF
C
C     Release buffer memory
C     ---------------------
C
      CALL MEMREL('RESOLV.buf',WORK,KFRSAV,KTEMP,KFREE,LFREE)
C
C     4-index transformation
C     ======================
C
      IF(.NOT.SKPTWO) THEN
        CALL GETTIM(CPU1,WALL1)
        CALL TRAHI(IPRRES,3)
C
C       Call driver for 4-index transformation
C       --------------------------------------
C
        CALL PAMTR1(WORK,KFREE,LFREE,IPRRES,KQ,KE,KIBE,
     &              NDMOQR,ICMOQR,NSTR,.FALSE.,DUMMY,
     &              .FALSE.,INTRES)
C
C        Print timing information
C
        CALL GETTIM(CPU2,WALL2)
        WALL   = WALL2 - WALL1
        CPU    = CPU2 - CPU1
        CPUTID = SECTID(WALL)
        WRITE(LUPRI,'(//A,A12)')
     &      ' Total wall time used in RESOLV :',CPUTID
        CPUTID = SECTID(CPU)
        WRITE(LUPRI,'(A,A12)')
     &      ' Total CPU  time used in RESOLV :',CPUTID
C
        CALL GTINFO(DAYTID)
        WRITE(LUPRI,'(/A,A24)') ' Transformation ended at : ',DAYTID
      ENDIF
C
C     Release all memory
C     ------------------
C
      CALL MEMREL('RESOLV',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C     Do full CI in open-shell manifold
C     =================================
C
 10   CONTINUE
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      CALL ASNUNIT(5,'-Yon',IERR)
      IF (IERR.NE.0)  WRITE (IW,'(/A,I4/)') ' Assign error ',IERR
#endif
      OPEN (5,FILE='DIRAC.INP')
      CALL GOSRES(WORK(KFREE),LFREE)
      CLOSE(5,STATUS='KEEP')
C
C     Restore DCBTRA information
C
      NOPAIR = TRANOP
      SCRTRA = SCRBUF
      IPRTRA = IPRBUF
      ISTRAT = ISTBUF
      THROUT = THRBUF
C
 9999 CONTINUE
      call dealloc(WORK)
      CALL QEXIT('RESOLV')
C
      RETURN
      END
C/* Deck SY_COSCI */
      SUBROUTINE SY_COSCI()
C***********************************************************************
C
C     Driver routine for 4-index transformation for
C     resolution of electronic states within open-shell manifold
C
C     Based on RESOLV
C     Written by T.Saue Jan 18 1999
!     Modified by sya, Jan 23, 2007
!       SYA separated Transformation routine from COS_CI.
!
!     Called by : PAMPSI (main/dirac.F)
C
C***********************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "infpar.h"
C
#include "maxorb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbres.h"
#include "dcbgen.h"
#include "dcbdhf.h"
!
!     ...SYA_MAXE
#include "cossya.h"
      CHARACTER DAYTID*24,SECTID*12,CPUTID*12
      LOGICAL   TOBE,TRANOP,SKPONE,SKPTWO
      DIMENSION NSTR(2,0:2,4),NSTRT(2),NSPC(2,0:2),
     &          KVEC(2,2),KQ(2,4),KE(2,4),KIBE(2,4),NQ(2),KQC(2)
CTROND: to fix !
      DIMENSION MINE(16),MAXE(16),NGASO(256)
      real(8), allocatable :: WORK(:)
!----------------------------------------------------------------------
C
!@m   CALL QENTER('RESOLV')
      CALL QENTER('SY_COSCI')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in SY_COSCI')

! .. definitions to satisfy ifort runtime check
      TRANOP = .false.
      SCRBUF = 0.0D0
      THRBUF = 0.0D0
      IPRBUF = 0
      ISTBUF = 0
      KFRSAV = KFREE
!@m   CALL TITLER('Resolution of open-shell states','*',120)
      CALL TITLER('Stand-alone COS_CI','*',120)
C
C     Determine number of active electrons
C
      IF (NAELEC.EQ.0) THEN
!@m      WRITE(LUPRI,'(A)')
!@m  &   '*WARNING: No open-shell electrons, no states to resolve!'
        WRITE(LUPRI,'(A)')
     &    '*WARNING: No open-shell electrons, no states to sy_cosci!'
!.s/sya,2007.06.14
!#       GOTO 9999
        IF( NAELEC .EQ. 0 ) THEN
          IACEL = 0
          DO I = 1, SYA_NOPEN
            IACEL = IACEL + SYA_IELC(I)
          END DO
          IF( IACEL .NE. 0 ) THEN
            NAELEC = IACEL
            WRITE(LUPRI,"(
     &        '*WARNING: NALEC has been calculated from SYA_IELC: ',
     &        I5)") IACEL
          ELSE
            DO I = 1, SYA_NOPEN
              IACEL = MAX(IACEL,SYA_MAXE(I))
            END DO
            NAELEC = IACEL
            WRITE(LUPRI,"(
     &        '*WARNING: NALEC has been calculated from SYA_MAXE: ',
     &        I5)") IACEL
          END IF
        END IF
!.q
      END IF
      IND = 1
      NGAS = NOPEN
      NSYMRP = 2*NFSYM
      DO IOPEN = 1,NOPEN
        ITEMP = 0
        DO IFRP = 1,NFSYM
          NGASO(IND) = NACSH(IFRP,IOPEN)
          NGASO(IND+1) = NGASO(IND)
          IND = IND + 2
          ITEMP = ITEMP + NACSH(IFRP,IOPEN)
        ENDDO
        MINE(IOPEN) = ITEMP
        MAXE(IOPEN) = MINE(IOPEN)
      ENDDO
C
C     Check for necessary files:
C       MRCONEE - one-electron integrals and other information
C       MDCINT  - two-electron integrals
C
      INQUIRE(FILE='MRCONEE',EXIST=SKPONE)
      IF(SKPONE) THEN
        WRITE(LUPRI,'(A)')
     & '*WARNING: MRCONEE found. Skipping 2-index transformation.'
      ENDIF
      INQUIRE(FILE='MDCINT',EXIST=SKPTWO)
      IF(SKPTWO) THEN
        WRITE(LUPRI,'(A)')
     & '*WARNING: MDCINT found. Skipping 4-index transformation.'
      ENDIF
      IF(SKPONE.AND.SKPTWO) GOTO 10
C
C     Transfer information from DCBRES to DCBTRA
C     ==========================================
C
      TRANOP = NOPAIR
      NOPAIR = .TRUE.
      SCRBUF = SCRTRA
      SCRTRA = SCRRES
      IPRBUF = IPRTRA
      IPRTRA = IPRRES
      ISTBUF = ISTRAT
      ISTRAT = ISTRES
      THRBUF = THROUT
      THROUT = RESOUT
C     number of double quaternionic classes for strategy 3 and 4
      NQQCLASS = NZ * NZ * NBSYM / NFSYM
C
C     Set up index arrays for active orbitals;
C     =======================================
C
      NQT  = 0
      ISAME(1) = 1
      NSTRT(1) = 0
      DO IFRP = 1,NFSYM
        CALL MEMGET2('INTE','KVEC',KVEC(IFRP,1),NASH(IFRP),
     &     WORK,KFREE,LFREE)
        CALL INDARR(WORK(KVEC(IFRP,1)),NISH(IFRP),NASH(IFRP),IPRINT)
        NSTR(IFRP,1,1) = NASH(IFRP)
        NSTR(IFRP,2,1) = 0
        NSTR(IFRP,0,1) = NSTR(IFRP,1,1) + NSTR(IFRP,2,1)
        NDMOQR(1,IFRP,1) = NFBAS(IFRP,0)
        NDMOQR(2,IFRP,1) = NASH(IFRP)
        ICMOQR(IFRP,1) = NQT + 1
        NSTRT(1) = NSTRT(1) +  NSTR(IFRP,0,1)
C
        NQ(IFRP) = NFBAS(IFRP,0)*NASH(IFRP)*NZ
        NQT      = NQT + NQ(IFRP)
        CALL MEMGET2('INTE','KIBE',KIBE(IFRP,1),NSTR(IFRP,0,1),
     &     WORK,KFREE,LFREE)
      ENDDO
      CALL MEMGET2('REAL','KQ',KQ(1,1),NQT,WORK,KFREE,LFREE)
      KQ(2,1) = KQ(1,1) + NQ(1)
      CALL MEMGET2('REAL','KE',KE(1,1),NSTRT(1),WORK,KFREE,LFREE)
      KE(2,1) = KE(1,1) + NSTR(1,0,1)
C
C     Copy for indices 2..4
C
      DO I = 2,4
        ISAME(I) = 1
        DO IFRP = 1,NFSYM
          KQ(IFRP,I)     = KQ(IFRP,1)
          KE(IFRP,I)     = KE(IFRP,1)
          KIBE(IFRP,I)   = KIBE(IFRP,1)
          NSTR(IFRP,0,I) = NSTR(IFRP,0,1)
          NSTR(IFRP,1,I) = NSTR(IFRP,1,1)
          NSTR(IFRP,2,I) = NSTR(IFRP,2,1)
          NDMOQR(1,IFRP,I) = NDMOQR(1,IFRP,1)
          NDMOQR(2,IFRP,I) = NDMOQR(2,IFRP,1)
          ICMOQR(IFRP,I)   = ICMOQR(IFRP,1)
        ENDDO
      ENDDO
C
      KTEMP = KFREE
      IF(.NOT.SKPONE) THEN
C
C       Set up index arrays for core orbitals
C       =======================================
        NQT  = 0
        DO IFRP = 1,NFSYM
          CALL MEMGET2('INTE','KVEC',KVEC(IFRP,2),NISH(IFRP),
     &       WORK,KFREE,LFREE)
          CALL INDARR(WORK(KVEC(IFRP,2)),0,NISH(IFRP),IPRINT)
          NSPC(IFRP,1) = NISH(IFRP)
          NSPC(IFRP,2) = 0
          NSPC(IFRP,0) = NSPC(IFRP,1) + NSPC(IFRP,2)
          NDMOQC(1,IFRP,1) = NFBAS(IFRP,0)
          NDMOQC(2,IFRP,1) = NISH(IFRP)
          ICMOQC(IFRP,1) = NQT + 1
          NQ(IFRP) = NFBAS(IFRP,0)*NISH(IFRP)*NZ
          NQT      = NQT + NQ(IFRP)
          KE(IFRP,2) = KE(IFRP,1)
        ENDDO
        CALL MEMGET2('REAL','KQC',KQC(1),NQT,WORK,KFREE,LFREE)
        KQC(2) = KQC(1) + NQ(1)
        NSTRT(2) = NSTRT(1)
      ENDIF
C
      WRITE(LUPRI,'(A,I5)') '- Number of active electrons:',NAELEC
      WRITE(LUPRI,'(A)')    '- Active orbitals:'
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,1)),0,0,0,NSTR)
      ENDDO
      WRITE(LUPRI,'(A)')    '- Core orbitals:'
      DO I = 1, NFSYM
         CALL TRAPRI(1,I,WORK(KVEC(I,2)),0,0,0,NSPC)
      ENDDO
C
C     Allocate buffer spaced for coefficients AND core Fock matrix
C
      CALL MEMGET2('REAL','KBUF',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','KEIG',KEIG,NORBT ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','KIBEIG',KIBEIG,NORBT,WORK,KFREE,LFREE)
C
C     Get all coefficients and eigenvalues
C     ====================================
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF (.NOT.TOBE) THEN
!@m      WRITE(LUPRI,'(/A)')
!@m  &      'RESOLV ERROR: MO coefficient file not present'
         WRITE(LUPRI,'(/A)')
     &      'SY_COSCI ERROR: MO coefficient file not present'
         CALL QUIT('RESOLV: No MO coefficients !')
      ELSE
         CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RESOLV')
         CALL REACMO(LUCOEF,'DFCOEF',WORK(KBUF),WORK(KEIG),WORK(KIBEIG),
     &               TOTERG,14)
         CLOSE(LUCOEF,STATUS='KEEP')
      ENDIF
C
C     Select active orbitals
C     ======================
C
!@d   DO IFRP = 1, NFSYM
!@d     CALL SELCFS (WORK(KBUF+ICMOQ(IFRP)),
!@d  &       IFRP,WORK(KQ(IFRP,1)),NSTR(IFRP,0,1),WORK(KVEC(IFRP,1)),
!@d  &       NSTR(IFRP,2,1),NSTR(IFRP,1,1),
!@d  &       NFBAS(IFRP,0),NORB(IFRP))
!@d     CALL SELEIG (WORK(KEIG+IORB(IFRP)),IFRP,WORK(KE(IFRP,1)),
!@d  &               WORK(KVEC(IFRP,1)),NSTR(IFRP,2,1),NSTR(IFRP,1,1))
!@d     CALL SELIBEIG(WORK(KIBEIG),IORB(IFRP),IFRP,
!@d  &                WORK(KIBE(IFRP,1)),WORK(KVEC(IFRP,1)),
!@d  &                NSTR(IFRP,2,1),NSTR(IFRP,1,1))
!@d   ENDDO
!@d   IF(.NOT.SKPONE) THEN
C
C       Select core orbitals
C       ======================
C
!@d     DO IFRP = 1, NFSYM
!@d       CALL SELCFS (WORK(KBUF+ICMOQ(IFRP)),
!@d  &         IFRP,WORK(KQC(IFRP)),NSPC(IFRP,0),WORK(KVEC(IFRP,2)),
!@d  &         NSPC(IFRP,2),NSPC(IFRP,1),
!@d  &         NFBAS(IFRP,0),NORB(IFRP))
!@d     ENDDO
C
C       2-index transformation
C       ======================
C
!@d     CALL TRAHI(IPRRES,2)
C
C       Get hold of the slaves if we run in parallel.
C        ( ITASK = 1 for Fock matrices )
C
!@d     IF (PARCAL) CALL DIRAC_PARCTL( HERFCK_PAR )
C
C       Transform core-fock matrix to active MO basis
C       ---------------------------------------------
C
!@d     CALL GTCFCK(WORK(KBUF),IPRRES,WORK(KQ(1,1)),WORK(KQ(1,2)),
!@d  &              WORK(KQC(1)),NDMOQR,ICMOQR,NSTR(1,0,1),
!@d  &              NDMOQC,ICMOQC,NSPC(1,0),
!@d  &              INTRES,WORK,KFREE,LFREE)
C
C       Write MOLFDIR file MRCONEE
C       --------------------------
C
!@d     CALL MEMGET('REAL',KORBE,NSTRT(1),WORK,KFREE,LFREE)
!@d     CALL MEMGET('INTE',KIBORB,NSTRT(1),WORK,KFREE,LFREE)
C
!@d     CALL MRCONE(WORK(KBUF),WORK(KORBE),WORK(KE(1,1)),WORK(KE(2,1)),
!@d  &              WORK(KIBORB),WORK(KIBE(1,1)),WORK(KIBE(2,1)),
!@d  &              NSPC,NSTR,NSTRT,IPRRES,
!@d  &              WORK,KFREE,LFREE)
C
C       Release the slaves if we run in parallel.
C        ( ITASK = -1 )
C
!@d     IF (PARCAL) CALL DIRAC_PARCTL( EXIT_NODEMENU )
!@d   ENDIF
C
C     Release buffer memory
C     ---------------------
C
!@d   CALL MEMREL('RESOLV.buf',WORK,KFRSAV,KTEMP,KFREE,LFREE)
C
C     4-index transformation
C     ======================
C
!@d   IF(.NOT.SKPTWO) THEN
!@d     CALL GETTIM(CPU1,WALL1)
!@d     CALL TRAHI(IPRRES,3)
C
C       Call driver for 4-index transformation
C       --------------------------------------
C
!@d     CALL PAMTR1(WORK,KFREE,LFREE,IPRRES,KQ,KE,KIBE,
!@d  &              NDMOQR,ICMOQR,NSTR,.FALSE.,DUMMY,
!@d  &              .TRUE.,.FALSE.,INTRES)
C
C        Print timing information
C
!@d     CALL GETTIM(CPU2,WALL2)
!@d     WALL   = WALL2 - WALL1
!@d     CPU    = CPU2 - CPU1
!@d     CPUTID = SECTID(WALL)
!@m     WRITE(LUPRI,'(//A,A12)')
!@m  &      ' Total wall time used in RESOLV :',CPUTID
!@d     WRITE(LUPRI,'(//A,A12)')
!@d  &      ' Total wall time used in SY_COSCI :',CPUTID
!@d     CPUTID = SECTID(CPU)
!@m     WRITE(LUPRI,'(A,A12)')
!@m  &      ' Total CPU  time used in RESOLV :',CPUTID
!@d     WRITE(LUPRI,'(A,A12)')
!@d  &      ' Total CPU  time used in SY_COSCI :',CPUTID
C
!@d     CALL GTINFO(DAYTID)
!@d     WRITE(LUPRI,'(/A,A24)') ' Transformation ended at : ',DAYTID
!@d   ENDIF
C
C     Release all memory
C     ------------------
C
!@m   CALL MEMREL('RESOLV',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
!@d   CALL MEMREL('SY_COSC',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
C     Do full CI in open-shell manifold
C     =================================
C
 10   CONTINUE
#if defined (SYS_CRAY) || defined (SYS_T3D)
!@d   IERR = 0
      CALL ASNUNIT(5,'-Yon',IERR)
      IF (IERR.NE.0)  WRITE (IW,'(/A,I4/)') ' Assign error ',IERR
#endif
      OPEN (5,FILE='DIRAC.INP')
!.s/sya,2007.02.04
!#    CALL GOSRES(WORK(KFREE),LFREE)
      CALL SYA_GOSRES(WORK(KFREE),LFREE)
!.q
      CLOSE(5,STATUS='KEEP')
C
C     Restore DCBTRA information
C
      NOPAIR = TRANOP
      SCRTRA = SCRBUF
      IPRTRA = IPRBUF
      ISTRAT = ISTBUF
      THROUT = THRBUF
C
!@m 9999 CALL QEXIT('RESOLV')
 9999 CONTINUE
      call dealloc(WORK)
      CALL QEXIT('SY_COSCI')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck indarr */
      SUBROUTINE INDARR(IVEC,IOFF,NVEC,IPRINT)
C***********************************************************************
C
C     Make index array
C
C     Written by T.Saue Jan 18 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
      DIMENSION IVEC(*)
      DO I = 1,NVEC
        IVEC(I) = IOFF + I
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck idbos */
      SUBROUTINE IDBOS(NBO,IBEIG,IBEIG2,VEC)
C***********************************************************************
C
C     Identify the boson irreps of recanonized virtual orbitals.
C     Note that the input array IBEIG is overwritten by the new data.
C     Written by L.Visscher may 1999
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
      PARAMETER(D0=0.0D0)
C
      DIMENSION IBEIG(NBO),IBEIG2(NBO),VEC(NBO,NBO)
C
      DO I = 1, NBO
         IBEIG2(I)= 0
         DO J = 1, NBO
            IF (VEC(J,I).NE.D0) THEN
               IBEIG2(I) = IBEIG(J)
            ENDIF
         ENDDO
      ENDDO
C
C     Consistency check
C
      DO I = 1, NBO
         DO J = 1, NBO
            IF (VEC(J,I).NE.D0) THEN
               IF (IBEIG2(I).NE.IBEIG(J))
     &         CALL QUIT('Spinfree canonization error')
            ENDIF
         ENDDO
      ENDDO
C
      DO I = 1, NBO
         IBEIG(I) = IBEIG2(I)
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C Start of parallel section : activitate compilation only when necessary
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
#if defined (VAR_MPI)
C/* Deck TRAPARI*/
      SUBROUTINE TRAPARI(WORK,KFREE,LFREE,KQ,KE,KIBE,
     &                   KINDX,KGAB,KDRIJ,INTFLG,NSTR,ANTIS,LMP2)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Luuk Visscher (based on HERNOD) August 1997
C     Polish by T.Saue Sep 8 1998
C
C    *******************************************************************
C    *   Initializes the nodes for the parallel 4-index transformation *
C    *******************************************************************
C
#include "implicit.h"
      PARAMETER(D0=0.0D0)
C
#include "infpar.h"
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "dcbgen.h"
C
#include "dorps.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dgroup.h"
#include "cbihr2.h"
C
      LOGICAL ANTIS,LMP2
      DIMENSION NSTR(2,0:2,4),KQ(2,4),KE(2,4),KIBE(2,4),
     &          WORK(*)
      LOGICAL NODV,NOPV,NOCONT,RETUR,NEWGEO,FINISH
C
      CALL QENTER('TRAPARI')
C
      KFRSAV = KFREE
C
      FINISH = .FALSE.
C
      IF (MYTID.EQ.0) THEN
C
C        **********************************************************
C        ***** Send HERMIT initialization data to the other nodes *
C        **********************************************************
C
C----------------------------------------------------------------
C        Setup information for the two-electron integralroutines.
C----------------------------------------------------------------
C
         NODV = .FALSE.
         NOPV = .FALSE.
         IF (NASHT .EQ. 0) NODV = .TRUE.
         IF (NASHT .LT. 2) NOPV = .TRUE.
         IPRALL =  0
         DO 100 I = 0,7
            DOREPS(I) = .TRUE.
  100    CONTINUE
         DO 110 I = 1,MXCENT
            DOCOOR(1,I) = .TRUE.
            DOCOOR(2,I) = .TRUE.
            DOCOOR(3,I) = .TRUE.
  110    CONTINUE
         NOCONT = .FALSE.
C
         MAXDIF = 0
         ITYPE  = 4
         IF(SCRTRA.GT.D0) THEN
           NDMAT = 2*NFSYM
         ELSE
           NDMAT = 0
         ENDIF
C
C        I2TYP will defined at each node
         CALL SDINIT(DUMMY,NDMAT,IDUMMY,IDUMMY,ITYPE,MAXDIF,IDUMMY,
     &               NODV,NOPV,NOCONT,TKTIME,RTNTWO,IDUMMY,
     &               IDUMMY,SCRTRA,WORK(KGAB),WORK(KDRIJ),DUMMY)
C
      ELSE
C
C        **********************************************************
C        ***** Receive HERMIT initialization data from master *****
C        **********************************************************
C
         CALL RVINIT(WORK,KFREE,LFREE,NDMAT,NFMAT,ITYPE,
     &               IATOM,MAXDIF,MYTID,NODV,NOPV,NOCONT,RETUR,
     &               TKTIME,NEWGEO,FINISH,I2TYP,ICEDIF,SCRTRA,
     &               KFMAT,KDMAT,KIFC,KIRD,KGAB,KDRIJ,KDMRSO)
      ENDIF
C
C        *****************************************************
C        ***** Send or receive DIRAC initialization data *****
C        *****************************************************
         CALL TRAPRD(WORK,KFREE,LFREE,KQ,KE,KIBE,KINDX,
     &               WORK(KGAB),WORK(KDRIJ),INTFLG,NSTR,ANTIS,LMP2)
C
C
C     Do not call MEMREL since the variable size arrays are represented
C     as pointers when the routine is called from a slave
C
#if defined (VAR_PFS)
      CALL GLOBAL_FILE_OPEN_WRITE
#endif
      CALL QEXIT('TRAPARI')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck traprd */
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE TRAPRD (WORK,KFREE,LFREE,KQ,KE,KIBE,KINDX,
     &                   GABRAO,DRIJ,INTFLG,NSTR,ANTIS,LMP2)
C
C     Communicate DIRAC common blocks for parallel 4-index transform.
C     Communicate arguments received from main program
C
C     Common blocks initialized : dcborb, dcbbas, dgroup, dcbtri,
C                                 dcbtr3, dcbmp2, comdis and symmet
C
C     Written by Luuk Visscher, august 1997.
C
      use interface_to_mpi
#include "implicit.h"
      PARAMETER(D0=0.0D0)
C
#include "infpar.h"
C
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbmp2.h"
#include "dcbmpt.h"
#include "dcbbas.h"
#include "symmet.h"
#include "dgroup.h"
#include "parint.h"
#include "blocks.h"
#include "twosta.h"
C
      LOGICAL ANTIS,LMP2,NSAME
      DIMENSION NSTR(2,0:2,4),KQ(2,4),KE(2,4),KIBE(2,4)
      DIMENSION DRIJ(*),GABRAO(*)
      DIMENSION WORK(*)
C
C
C     ********************************
C     ****** Send COMMON blocks ******
C     ********************************
C
C
C     This is highly unauthorized!!!
C     No F77 standard says that elements in a common block are
C     consequtive, and, in fact, with -Ofast on SGI the compiler
C     may pad common blocks with junk. /jth Feb 12 2001
C
C     This routine should be improved, the problem is that the slaves
C     do not read input so do not get their common blocks initialized.
C     A rigorous solution is to transfer the input file to all nodes
C     and have the slaves go though it as well. This requires some
C     non-trivial changes of the parallelization setup. LV 18-5-01
C
C     dcborb.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCBORB,I2_DCBORB)
      call interface_MPI_BCAST(I1_DCBORB,NUMELM,MPARID,
     &                         global_communicator)
C
C     dcbbas.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCBBAS,I2_DCBBAS)
      call interface_MPI_BCAST(I1_DCBBAS,NUMELM,MPARID,
     &                         global_communicator)
C
C     dcbham.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DCLHAM,I2_DCLHAM)
      call interface_mpi_bcast_l0(QED,NUMELM,MPARID,
     &                                global_communicator)
C
C     dgroup.h
C     --------
C
      NUMELM = ICOMMSIZE(I1_DGROUP,I2_DGROUP)
      call interface_MPI_BCAST(I1_DGROUP,NUMELM,MPARID,
     &                         global_communicator)
      call interface_mpi_bcast_l0(LINEAR,1,MPARID,
     &                                global_communicator)
C
C     dcbtra.h
C     --------
C
C     integers...
C
      NUMELM = NDCBTRI
      call interface_MPI_BCAST(IPRTRA,NUMELM,MPARID,global_communicator)
C
C     reals...
C
      NUMELM = NDCBTRR
      call interface_MPI_BCAST(THROUT,NUMELM,MPARID,global_communicator)
C
C     logicals...
C
      NUMELM = NDCBTRL
      call interface_mpi_bcast_l0(TRA_ANTIS,NUMELM,MPARID,
     &                                global_communicator)
C
C     comdis.h
C     --------
C
      NUMELM = NCOMDIS
      call interface_MPI_BCAST(NSPCK,NUMELM,MPARID,global_communicator)
C
C     symmet.h
C     --------
C
C     integers...
C
      NUMELM = ICOMMSIZE(I1_SYMMTI,I2_SYMMTI)
      call interface_MPI_BCAST(I1_SYMMTI,NUMELM,MPARID,
     &                         global_communicator)
C
C     reals...
C
      NUMELM = NSYMMTR
      call interface_MPI_BCAST(FMULT,NUMELM,MPARID,global_communicator)
C
C     The rest
C     --------
C
      NUMELM = 1
      CALL interface_mpi_bcast_l0(TESTLS,NUMELM,MPARID,
     &                                      global_communicator)
      NUMELM = 1
      CALL interface_mpi_bcast_l0(TESTSL,NUMELM,MPARID,
     &                                      global_communicator)
      NUMELM = 1
      CALL interface_mpi_bcast_l0(MP2ORG,NUMELM,MPARID,
     &                                      global_communicator)
      NUMELM = 1
      CALL interface_mpi_bcast_l0(GTHSHL,NUMELM,MPARID,
     &                                      global_communicator)
      NUMELM = 1
      CALL interface_mpi_bcast_l0(CMMNFO,NUMELM,MPARID,
     &                                      global_communicator)
      NUMELM = 1
      CALL interface_MPI_BCAST(TPRI34,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 1
      CALL interface_MPI_BCAST(TPRI44,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 1
      CALL interface_MPI_BCAST(INTFLG,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 1
      CALL interface_mpi_bcast_l0(ANTIS,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 1
      CALL interface_mpi_bcast_l0(LMP2,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 24
      CALL interface_MPI_BCAST(NSTR,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 16
      CALL interface_MPI_BCAST(NDMOQR,NUMELM,MPARID,
     &               global_communicator)
      NUMELM = 8
      CALL interface_MPI_BCAST(ICMOQR,NUMELM,MPARID,
     &               global_communicator)
C
C
C
C     HERMIT memory estimates
C     -----------------------
C
      CALL interface_MPI_BCAST(MXHERMEM,4,MPARID,global_communicator)
C
C
      IF (MYTID.EQ.MPARID) THEN
C
C     *************************
C     ****** MASTER send ******
C     *************************
C
C       Symmetry packing information
        CALL interface_mpi_bcast_i1_work_f77(WORK(KINDX),NTBAS(0)*3,
     &                                       MPARID,global_communicator)
C       Coefficients; use sameness relations to compress
        DO I = 1,4
          NSAME = .TRUE.
          DO K = (I-1),1,-1
          IF(ISAME(I).EQ.ISAME(K)) THEN
            NSAME = .FALSE.
            JSAME = K
          ENDIF
          ENDDO
          IF(NSAME) THEN
            NSIZE = 0
            DO J = 1, NFSYM
              NSIZE = NSIZE + NDMOQR(1,J,I)*NDMOQR(2,J,I)*NZ
            ENDDO
            CALL interface_mpi_bcast_r1_work_f77(WORK(KQ(1,I)),NSIZE,
     &                 MPARID,global_communicator)
          ENDIF
        ENDDO
C       Send orbital energies
        DO I = 1,4
          NSIZE = 0
          DO J = 1,NFSYM
            NSIZE = NSIZE + NSTR(J,0,I)
          ENDDO
          CALL interface_mpi_bcast_r1_work_f77(WORK(KE(1,I)),NSIZE,
     &                   MPARID,global_communicator)
        ENDDO
C       Send orbital (single group) info
        DO I = 1,4
          DO J = 1,NFSYM
            CALL interface_mpi_bcast_i1_work_f77(WORK(KIBE(J,I)),
     &                     NSTR(J,0,I),MPARID,global_communicator)
          ENDDO
        ENDDO
C       If screening send GAB-matrix and screening matrix
        IF(SCRTRA.GT.D0) THEN
          N2GAB   = NSYMBL*NSYMBL
          NDRIJ   = N2GAB*NFSYM*2
          CALL interface_mpi_bcast_r1_work_f77(GABRAO,N2GAB,
     &                   MPARID,global_communicator)
          CALL interface_mpi_bcast_r1_work_f77(DRIJ,NDRIJ,
     &                   MPARID,global_communicator)
        ENDIF
      ELSE
C
C     **************************
C     ****** SLAVE receive *****
C     **************************
C
C       Symmetry packing information
        CALL MEMGET2('INTE','KINDX',KINDX,3*NTBAS(0),WORK,KFREE,LFREE)
        CALL interface_mpi_bcast_i1_work_f77(WORK(KINDX),NTBAS(0)*3,
     &                                       MPARID,global_communicator)
C       Coefficients; use sameness relations to compress
        DO I = 1,4
          NSAME = .TRUE.
          DO K = (I-1),1,-1
          IF(ISAME(I).EQ.ISAME(K)) THEN
            NSAME = .FALSE.
            JSAME = K
          ENDIF
          ENDDO
          IF(NSAME) THEN
            NSIZE = 0
            DO J = 1, NFSYM
              NSIZE = NSIZE + NDMOQR(1,J,I)*NDMOQR(2,J,I)*NZ
            ENDDO
            CALL MEMGET2('REAL','KQ',KQ(1,I),NSIZE,WORK,KFREE,LFREE)
            KQ(2,I) = KQ(1,I) + NDMOQR(1,1,I)*NDMOQR(2,1,I)*NZ
            CALL interface_mpi_bcast_r1_work_f77(WORK(KQ(1,I)),NSIZE,
     &              MPARID,global_communicator)
          ELSE
            KQ(1,I) = KQ(1,JSAME)
            KQ(2,I) = KQ(2,JSAME)
          ENDIF
        ENDDO
C       Receive orbital energies
        DO I = 1,4
          NSIZE = 0
          DO J = 1,NFSYM
            NSIZE = NSIZE + NSTR(J,0,I)
          ENDDO
          CALL MEMGET2('REAL','KE',KE(1,I),NSIZE,WORK,KFREE,LFREE)
          KE(2,I) = KE(1,I) + NSTR(1,0,I)
          CALL interface_mpi_bcast_r1_work_f77(WORK(KE(1,I)),NSIZE,
     &                   MPARID,global_communicator)
        ENDDO
C       Receive orbital (single group) info
        KIBE = 1 ! initialize these memory locations to avoid bound check errors if nfsym = 1
        DO I = 1,4
          DO J = 1,NFSYM
            CALL MEMGET2('REAL','KIBE',KIBE(J,I),NSTR(J,0,I),
     &         WORK,KFREE,LFREE)
            CALL interface_mpi_bcast_i1_work_f77(WORK(KIBE(J,I)),
     &                                       NSTR(J,0,I),
     &                                       MPARID,global_communicator)
          ENDDO
        ENDDO
C       If screening receive GAB-matrix and screening matrix
        IF(SCRTRA.GT.D0) THEN
          N2GAB   = NSYMBL*NSYMBL
          NDRIJ   = N2GAB*NFSYM*2
          CALL interface_mpi_bcast_r1_work_f77(GABRAO,N2GAB,MPARID,
     &                                     global_communicator)
          CALL interface_mpi_bcast_r1_work_f77(DRIJ,NDRIJ,MPARID,
     &                                     global_communicator)
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRANOD*/
      SUBROUTINE TRANOD()
C
C     Written by Luuk Visscher, August 1997
C     Extensions by T.Saue Sep 1998
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      use memory_allocator
      use interface_to_mpi

#include "implicit.h"
#include "priunit.h"
C
#include "infpar.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dcbtr3.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "aovec.h"
#include "maxorb.h"
#include "blocks.h"
      LOGICAL   ANTIS,LMP2
      INTEGER   NSTR(2,0:2,4),KQ(2,4),KE(2,4),KIBE(2,4)
      REAL(8)   DINTSKP(8),TIM1(7)
      real(8), allocatable :: WORK(:)
C
      CALL QENTER('TRANOD')
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in TRANOD')
      KFRSAV = KFREE
      CALL DZERO(DINTSKP,8)

C
C     Initialize the common blocks.
C
      PARHER = .TRUE.
      IPRINT = -1
      CALL TRAPARI(WORK,KFREE,LFREE,KQ,KE,KIBE,
     &             KINDX,KGAB,KDRIJ,
     &             INTFLG,NSTR,ANTIS,LMP2)
C
      IF (ISTRAT.EQ.1) THEN
         WRITE(LUPRI,'(A)')
     &     'TRANOD: Strategy 1 not implemented in parallel.'
         CALL QUIT('TRANOD: No parallel strategy 1 (yet)!')
      ELSEIF (ISTRAT.EQ.2) THEN
         IF (.NOT. MP2ORG) THEN
            CALL TR2NOD(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &           KINDX,WORK(KGAB),WORK(KDRIJ),INTFLG,NSTR,ANTIS,LMP2,
     &           DINTSKP,TIM1)
         ELSE
            CALL TR2PAR(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &           KINDX,WORK(KGAB),WORK(KDRIJ),INTFLG,NSTR,ANTIS,
     &           LMP2,EMP2,DINTSKP,TIM1)
         END IF
      ELSEIF (ISTRAT.EQ.3) THEN
         WRITE(LUPRI,'(A)')
     &     'TRANOD: Strategy 3 not implemented in parallel.'
         CALL QUIT('TRANOD: No parallel strategy 3 (yet)!')
      ELSEIF (ISTRAT.EQ.4) THEN
          CALL TR4PAR(WORK,KFREE,LFREE,IPRINT,KQ,KE,KIBE,
     &             KINDX,WORK(KGAB),WORK(KDRIJ),INTFLG,NSTR,ANTIS,
     &             LMP2,EMP2,DINTSKP)
      ELSEIF (ISTRAT.EQ.5) THEN
         WRITE(LUPRI,'(A)')
     &     'TRANOD: Strategy 5 not implemented in parallel.'
         CALL QUIT('TRANOD: No parallel strategy 5 (yet)!')
      ELSEIF (ISTRAT.EQ.6) THEN
          CALL TRDR6T(WORK,KFREE,LFREE,IPRINT,INTFLG,
     &                KGAB,KDRIJ,NSTR,
     &                ANTIS,LMP2,KINDX,KQ,KE,KIBE,DINTSKP)
      ELSE
         WRITE(LUPRI,'(A,I10,A)')
     &     'TRANOD: Strategy ',ISTRAT,' not implemented'
         CALL QUIT('TRANOD: Unknown strategy !')
      ENDIF
C
C     We check whether all the work is done and return to the
C     main slave menu.
C
      CALL interface_MPI_BCAST(NTEST,1,MPARID,global_communicator)
      IF (NTEST.NE.-1) CALL QUIT('TRANOD: Received wrong NTEST message')
C
#if defined (VAR_PFS)
      IF( ISTRAT .ne. 6 )THEN
C
C       close global file (opened in TRAPARI)
C
        CALL GLOBAL_FILE_CLOSE_READ
C
      END IF
#endif
C
      CALL MEMREL('TRANOD',WORK,1,KFRSAV,KFREE,LFREE)
      call dealloc(WORK)
      PARHER = .FALSE.
      CALL QEXIT('TRANOD')
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRARES*/
      SUBROUTINE TRARES(WORK,KFREE,LFREE,DINTSKP,LMP2,EMP2,TIM1,WALLBEF)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
C     Written by Trond Saue Sep 15 1998
C
C     Slave --> Master:
C       - Screening information
C       - MP2 results
C***********************************************************************
      use interface_to_mpi
#include "implicit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
C
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "dcbgen.h"
C
#include "dorps.h"
#include "dcborb.h"
#include "dcbtra.h"
#include "dgroup.h"
C
      LOGICAL LMP2
      DIMENSION DINTSKP(8),TIM1(7,*),WORK(*),IBUFMPI(5)
C
      KFRSAV = KFREE
      IF (MYTID.EQ.0) THEN
C
C     **************************
C     ***** Master receive *****
C     **************************
C
        IF(LMP2) EBUF = D0
        CALL MEMGET2('REAL','KBUF',KBUF,8,WORK,KFREE,LFREE)
        IF(ISTRAT.EQ.4) THEN
C         For historical reasons we need to set the third variable in the buf to -1 to
C         signal the slave to stop computing & transforming integrals.
          CALL IZERO (IBUFMPI,5)
          IBUFMPI(3) = -1
        ENDIF
        DO I = 1,NODES
 100      CONTINUE
#if defined (MPE)
          CALL MPE_LOG_EVENT(3, 0, "start wait")
#endif
          call interface_mpi_RECV(ITEST,1,df_MPI_ANY_SOURCE,20,
     &                            global_communicator,ISTAT)
#if defined (MPE)
          CALL MPE_LOG_EVENT(4, 0, "end wait")
#endif
          IF (ITEST .NE. 0) THEN
C         ... this node sends the wrong message, abort
             CALL QUIT('TRARES: wrong code sent by slave')
             GO TO 100
          END IF
          NWHO = ISTAT(df_MPI_SOURCE)
#if defined (MPE)
          CALL MPE_LOG_EVENT(7, 0, "start rcvint")
#endif
          CALL GETTIM(CPUST,WALLST)
          WRITE(LUPRI,'(3X,(A,I4),A,F14.2,A)')
     &           '- Process ',NWHO,
     &           ' finished at ',WALLST-WALLBEF,' seconds after start'
          CALL FLSHFO(LUPRI)
          IF(ISTRAT.EQ.4) THEN
            call interface_mpi_SEND(IBUFMPI,5,NWHO,30,
     &                    global_communicator)
          ELSEIF(ISTRAT.EQ.2) THEN
            ITEST = -1
            call interface_mpi_SEND(ITEST,1,NWHO,30,
     &                    global_communicator)
          ELSE
            WRITE(6,*) 'Inside TRARES in scheme ..',ISTRAT
            CALL QUIT('TRARES not adapted for this strategy !')
          ENDIF
          call interface_mpi_RECV(WORK(KBUF),8,
     &          NWHO,40,global_communicator)
          CALL DAXPY(8,D1,WORK(KBUF),1,DINTSKP,1)
          IF(LMP2) THEN
            call interface_mpi_RECV(EBUF,1,
     &            NWHO,40,global_communicator)
            EMP2 = EMP2 + EBUF
            call interface_mpi_RECV(TIM1(1,NWHO),7,
     &            NWHO,40,global_communicator)
          ENDIF
#if defined (MPE)
          CALL MPE_LOG_EVENT(8, 0, "end rcvint")
#endif
        ENDDO
        CALL MEMREL('TRARES',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
      ELSE
C
C     ***********************
C     ***** Slave send  *****
C     ***********************
C
        call interface_mpi_SEND(DINTSKP,8,
     &          MPARID,40,global_communicator)
        IF(LMP2) THEN
          call interface_mpi_SEND(EMP2,1,
     &          MPARID,40,global_communicator)
          call interface_mpi_send_r1_work_f77(TIM1,7,
     &          MPARID,40,global_communicator)
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C& End of parallel section : activitate compilation only when necessary
#endif  /* ifdef VAR_MPI */
C/* Deck MXMEMHER */
      FUNCTION MXMEMHER(IC,I2TYP)
C
C HJAaJ Aug 2001: return maximum memory needed in Hermit for the
C                 type of integrals specified with IC, I2TYP
C
#include "implicit.h"
#include "maxaqn.h"
#include "twosta.h"
C
C     Memory requirements for HERMIT:
C
C     IC = 1: (LL|**) integrals
C     IC = 2: (SS|**) integrals
C     IC = 3: (SL|**) integrals
C
C     I2TYP =  1: (**|LL) integrals
C     I2TYP =  2: (**|SS) integrals
C     I2TYP = 12: (**|LL) + (**|SS) integrals
C     I2TYP =  4: (**|SL) integrals
C
      IF ( I2TYP .EQ. 0 ) THEN
C        ...(GG|GG) integrals
         MWHER = MAX(MXHERMEM(1),MXHERMEM(2),MXHERMEM(3))
      ELSE IF ( I2TYP .EQ. 1 ) THEN
         IF ( IC .EQ. 1 ) THEN
C           ...(LL|LL) integrals
            MWHER = MXHERMEM(1)
         ELSE
C           ...(SS|LL) integrals
            MWHER = MXHERMEM(2)
         END IF
      ELSE IF ( I2TYP .EQ. 2 ) THEN
         IF ( IC .EQ. 1 ) THEN
C           ...(LL|SS) integrals
            MWHER = MXHERMEM(2)
         ELSE
C           ...(SS|SS) integrals
            MWHER = MXHERMEM(3)
         END IF
      ELSE IF ( I2TYP .EQ. 12 ) THEN
         IF ( IC .EQ. 1 ) THEN
C           ...(LL|LL) + (LL|SS) integrals
            MWHER = MAX(MXHERMEM(1),MXHERMEM(2))
         ELSE
C           ...(SS|LL) + (SS|SS) integrals
            MWHER = MAX(MXHERMEM(2),MXHERMEM(3))
         END IF
      ELSE IF ( I2TYP .EQ. 4 ) THEN
C        ...(SL|SL) integrals
Cluuk    MWHER = MXHERMEM(4)
Cluuk   not yet implemented use value for SS
         MWHER = MXHERMEM(3)
      ELSE
         CALL QUIT('MXHERMEM: unknown I2TYP')
      END IF
      MXMEMHER = MWHER
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck GOSRES */
      SUBROUTINE GOSRES(WORK,LWORK)
C***********************************************************************
C
C     Driver routine setting up correctly
C     the resolution of electronic states within 
C     open-shell manifold
C
C     Written on a rainy day in Nagoya 2006....
C     Trond Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (NORDER=16)
      CHARACTER*4 REPA(4*NORDER)
#include "dcbtra.h"
#include "dcbdhf.h"
      DIMENSION WORK(LWORK)
#include "memint.h"
C
      CALL OPNFIL(LUMLF1,'MRCONEE','OLD','GOSRES')
      READ(LUMLF1) NORB
C.....Skip symmetry information for parent group
      READ(LUMLF1) 
C.....Read symmetry information for Abelian subgroup
      READ(LUMLF1) NREP,(REPA(IRP),IRP=1,2*NREP)
      CLOSE(LUMLF1,STATUS='KEEP')

C.....Allocate memory for orbital symmetry information
      NDIM = NREP*NOPEN
      CALL MEMGET2('INTE','KGASO',KGASO,NDIM, WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','KIELC',KIELC,NOPEN,WORK,KFREE,LFREE)
C      
      CALL MEMGET2('INTE','KIPA',KIRPA,NORB,WORK,KFREE,LFREE)
      CALL GOSRE1(REPA,NREP,WORK(KIRPA),
     &            WORK(KGASO),WORK(KIELC),WORK(KFREE),LFREE)
      CALL MEMREL('GOSRES',WORK,1,KIRPA,KFREE,LFREE)
C.....Call GOSCIP module
      CALL MEMGET2('INTE','KIGAS',KIGAS,NOPEN,WORK,KFREE,LFREE)
      CALL PAMGOS(WORK(KFREE),LFREE,NAELEC_DHF,NOPEN,WORK(KIGAS),
     &            WORK(KIELC),WORK(KIELC),WORK(KGASO))
      CALL MEMREL('PAMGOS',WORK,1,KFREE,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/*   Deck GOSRE1 */
      SUBROUTINE GOSRE1(REPA,NREP,IRPAMO,NGASO,IELC,WORK,LWORK)
C***********************************************************************
C
C     Driver routine setting up correctly
C     the resolution of electronic states within 
C     open-shell manifold
C
C     Written on a rainy day in Nagoya 2006....
C     Trond Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2=2.0D0)
      CHARACTER*4 REPA(*)
#include "dcbtra.h"
#include "dgroup.h"
#include "dcbdhf.h"
      DIMENSION IRPAMO(*),NGASO(NREP,*),IELC(*),WORK(LWORK)
C
      CALL OPNFIL(LUMLF1,'MRCONEE','OLD','GOSRES')
      READ (LUMLF1) NORB
      READ (LUMLF1) 
      READ (LUMLF1) NREP,(REPA(IRP),IRP=1,2*NREP)
      READ (LUMLF1) 
      READ (LUMLF1) (IDUM,IRPAMO(IMO),DUM,IMO=1,NORB)
      CLOSE(LUMLF1,STATUS='KEEP')
C
      NDIM = NREP*NOPEN
      CALL IZERO(NGASO,NDIM)
      CALL IZERO(IELC,NOPEN)
      IORB = 0
      DO IFRP = 1,NFSYM
        DO IOPEN = 1,NOPEN
          IELC(IOPEN) = IELC(IOPEN) + NACSH(IFRP,IOPEN)
          DO I = 1,NACSH(IFRP,IOPEN)
            IORB = IORB + 1
C...........Add orbital to GAS
            IREP = IRPAMO(IORB)
            NGASO(IREP,IOPEN) = NGASO(IREP,IOPEN) + 1
C...........Add Kramers partner to GAS
            IREP = IRPAMO(IORB+NASH_DHF(IFRP))
            NGASO(IREP,IOPEN) = NGASO(IREP,IOPEN) + 1
          ENDDO
        ENDDO
        IORB = IORB + NASH_DHF(IFRP)
      ENDDO
C
      DO IOPEN = 1,NOPEN
        IELC(IOPEN) = NINT(DF(IOPEN)*D2*IELC(IOPEN))
      ENDDO   
C
      WRITE(LUPRI,'(A)') 'GASRES: Set up the following information:'
      WRITE(LUPRI,'(5X,32(1X,A4))') (REPA(I),I=1,NREP)
      DO IOPEN = 1,NOPEN
        WRITE(LUPRI,'(32I5)') 
     &    IELC(IOPEN),(NGASO(IREP,IOPEN),IREP = 1,NREP)
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SYA_GOSRES */
      SUBROUTINE SYA_GOSRES(WORK,LWORK)
C***********************************************************************
C
C     Driver routine setting up correctly
C     the resolution of electronic states within 
C     open-shell manifold
C
C     Written on a rainy day in Nagoya 2006....
C     Trond Saue
!
!     Modified by sya, Jan 31, 2007
!
!     Called by: SYA_COSCI
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (NORDER=16)
      CHARACTER*4 REPA(4*NORDER)
#include "dcbtra.h"
#include "dcbdhf.h"
!.s/sya,2007.02.02
#include "cossya.h"
!.q  
      DIMENSION WORK(LWORK)
#include "memint.h"
C----------------------------------------------------------------------
      CALL OPNFIL(LUMLF1,'MRCONEE','OLD','GOSRES')
      READ(LUMLF1) NORB
C.....Skip symmetry information for parent group
      READ(LUMLF1) 
C.....Read symmetry information for Abelian subgroup
      READ(LUMLF1) NREP,(REPA(IRP),IRP=1,2*NREP)
      CLOSE(LUMLF1,STATUS='KEEP')
!.s/sya,2007.09.21/
!     ...set NOPEN from COSCI input
      IF( NOPEN .EQ. 0 ) THEN
        IF( SYA_NOPEN .GT. 0 ) THEN
          WRITE(LUPRI,'(A)') 
     &      'SYA_GOSRES: NOPEN being replaced by SYA_NOPEN.'
          NOPEN = SYA_NOPEN
        ELSE
          WRITE(LUPRI,'(A,2I5)') 
     &      'GOSRES: Inconsitent data.',NOPEN,SYA_NOPEN
          STOP
        END IF
      END IF
!.q
C.....Allocate memory for orbital symmetry information
      NDIM = NREP*NOPEN
      CALL MEMGET2('INTE','KGASO',KGASO,NDIM, WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','KIELC',KIELC,NOPEN,WORK,KFREE,LFREE)
!.s/sya,2007.01.31
      CALL MEMGET2('INTE','KIELC2',KIELC2,NOPEN,WORK,KFREE,LFREE)
!.q
C      
      CALL MEMGET2('INTE','KIRPA',KIRPA,NORB,WORK,KFREE,LFREE)
!.s/sya,2007.02.04
!#    CALL GOSRE1(REPA,NREP,WORK(KIRPA),
!#   &            WORK(KGASO),WORK(KIELC),WORK(KFREE),LFREE)
      CALL SYA_GOSRE1(REPA,NREP,WORK(KIRPA),
     &            WORK(KGASO),WORK(KIELC),WORK(KFREE),LFREE)
!.q
      CALL MEMREL('GOSRES',WORK,1,KIRPA,KFREE,LFREE)
C.....Call GOSCIP module
      CALL MEMGET2('INTE','KIGAS',KIGAS,NOPEN,WORK,KFREE,LFREE)
!.s/sya,2007.01.31
!#    CALL PAMGOS(WORK(KFREE),LFREE,NAELEC_DHF,NOPEN,WORK(KIGAS),
!#   &            WORK(KIELC),WORK(KIELC),WORK(KGASO))
      CALL SY_PAMGOS(WORK(KFREE),LFREE,NAELEC_DHF,NOPEN,WORK(KIGAS),
     &      WORK(KIELC),WORK(KIELC2),WORK(KGASO))
!.q
      CALL MEMREL('GOSRES',WORK,1,KFREE,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/*   Deck SYA_GOSRE1 */
      SUBROUTINE SYA_GOSRE1(REPA,NREP,IRPAMO,NGASO,IELC,WORK,LWORK)
C***********************************************************************
C
C     Driver routine setting up correctly
C     the resolution of electronic states within 
C     open-shell manifold
C
C     Written on a rainy day in Nagoya 2006....
C     Trond Saue
!
!     Last-update: sya, 2007.02.01, ulp
!
!     Called by : SYA_GOSRES (moltra/trapam.F)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2=2.0D0)
      CHARACTER*4 REPA(*)
#include "dcbtra.h"
#include "dgroup.h"
#include "dcbdhf.h"
!.s/sya,2007.02.01
#include "cossya.h"
!.q
      DIMENSION IRPAMO(*),NGASO(NREP,*),IELC(*),WORK(LWORK)
!----------------------------------------------------------------------
! Comment by sya, 2007.02.02.
!     (REPA(IRP),IRP=1,2*NREP)
!       latter part from boson irrep.
!----------------------------------------------------------------------
      CALL OPNFIL(LUMLF1,'MRCONEE','OLD','GOSRES')
      READ (LUMLF1) NORB
      READ (LUMLF1) 
      READ (LUMLF1) NREP,(REPA(IRP),IRP=1,2*NREP)
      READ (LUMLF1) 
      READ (LUMLF1) (IDUM,IRPAMO(IMO),DUM,IMO=1,NORB)
      CLOSE(LUMLF1,STATUS='KEEP')
C
      NDIM = NREP*NOPEN
      CALL IZERO(NGASO,NDIM)
      CALL IZERO(IELC,NOPEN)
      IORB = 0
      DO IFRP = 1,NFSYM
        DO IOPEN = 1,NOPEN
          IELC(IOPEN) = IELC(IOPEN) + NACSH(IFRP,IOPEN)
          DO I = 1,NACSH(IFRP,IOPEN)
            IORB = IORB + 1
C...........Add orbital to GAS
            IREP = IRPAMO(IORB)
            NGASO(IREP,IOPEN) = NGASO(IREP,IOPEN) + 1
C...........Add Kramers partner to GAS
            IREP = IRPAMO(IORB+NASH_DHF(IFRP))
            NGASO(IREP,IOPEN) = NGASO(IREP,IOPEN) + 1
          ENDDO
        ENDDO
        IORB = IORB + NASH_DHF(IFRP)
      ENDDO
C
      DO IOPEN = 1,NOPEN
        IELC(IOPEN) = NINT(DF(IOPEN)*D2*IELC(IOPEN))
      ENDDO   
!
!.s/sya,2007.02.01
!
!     Recover data defined by COSCI input data.
!     Overwrite the data from DHF by COSCI input data.
!
      IF( KCOSINP ) THEN
         DO I = 1, SYA_NOPEN
            IELC(I) = SYA_IELC(I)
            DO J = 1, NREP
               NGASO(J,I) = SYA_GASO(J,I)
            END DO
         END DO
         WRITE(LUPRI,'(A)') 
     &      'GOSRES: Recovering the data defined by COSCI input.'
      END IF
!.q
C
      WRITE(LUPRI,'(A)') 'GOSRES: Set up the following information:'
      WRITE(LUPRI,'(5X,32(1X,A4))') (REPA(I),I=1,NREP)
      DO IOPEN = 1,NOPEN
        WRITE(LUPRI,'(32I5)') 
     &    IELC(IOPEN),(NGASO(IREP,IOPEN),IREP = 1,NREP)
      ENDDO
!
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck TRAPAM */
      SUBROUTINE ADDCORE2(NLIST,TLIST,TNUM,
     &                     SLIST,SNUM)
C***********************************************************************
C
C     Written by J.Pototschnig, Fall 2019
C
C***********************************************************************
C     add SLIST to TLIST
#include "priunit.h"
      INTEGER NLIST,TNUM,SNUM,I
      INTEGER TLIST(NLIST),SLIST(NLIST)

      DO I=1,SNUM
        TLIST(TNUM+I)=SLIST(I)
      ENDDO
      TNUM=TNUM+SNUM
      WRITE(LUPRI,'(6X,12I5)') (SLIST(I),I=1,SNUM)

      RETURN
      END
