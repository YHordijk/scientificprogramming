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

***********************************************************************

      SUBROUTINE GASDIAS_REL(IMK2NUM,NMS2VAL,MS2VAL,NAEL,
     &           IASTR,NBEL,IBSTR,
     &           NORB,DIAG,NSMST,NIRREP,
     &           H1A,H1B,XA,XB,SCR,RJ,RK,
     &           RK2,NSSOA,NSSOB,LUDIA,ECORE,
     &           PLSIGN,PSSIGN,IPRNT,NTOOB,ICISTR,RJKAA,IHAM12,
     &           IBLTP,NBLOCK,IBLKFO,IASPGPTP,IBSPGPTP,
     &           ISPINFREE,NPARBLOCK,ILASTBLK,IDIA_OFFSET)
*
* Calculate determinant diagonal
* Turbo-ras version
*
* Driven by IBLKFO, May 97
*
* ========================
* General symmetry version
* ========================
*
* Jeppe Olsen, July 1995, GAS version
*
* Fixed bug in first exchange term.    Timo Fleig, 07.Aug.98
* Fixed bug in third integral term.    Timo and Jeppe, 25.Aug.98
* Correct sym-type ordering for all integrals implemented.
*                                      Timo Fleig, fall 2000
*
* IHAM12 = 1 => only one-body part
*        = 2 =>      one+two-body part
*
* separate alpha and beta 1-electron integrals  May 1998
      use interface_to_mpi
      IMPLICIT REAL*8           (A-H,O-Z)
#include "ipoist8.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*
*.General input
      DIMENSION NSSOA(NSMST,*),NSSOB(NSMST,*)
      DIMENSION H1A(NORB), H1B(NORB), NPARBLOCK(*)
*. Specific input
      DIMENSION IBLTP(*),IBLKFO(8,NBLOCK)
      dimension MS2VAL(NMS2VAL)
*. Scratch
      DIMENSION RJ(NTOOB,NTOOB),RK(NTOOB,NTOOB)
      DIMENSION RK2(NTOOB,NTOOB)
*
      DIMENSION XA(NORB),XB(NORB),SCR(2*NORB)
      DIMENSION IASTR(NAEL,*),IBSTR(NBEL,*)
      DIMENSION RJKAA(*)
#if defined (VAR_MPI2)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IDIA_OFFSET
#endif
*. Output
      DIMENSION DIAG(*)
*
      call qenter('GASDIAS')
*
!     NTEST =  5000 ! debug
      NTEST =  0 
      NTEST = MAX(NTEST,IPRNT/1000)
*
      IF(PSSIGN.EQ.-1.0D0) THEN
         XADD = 1000000.0
      ELSE
         XADD = 0.0D0
      END IF
*
      IF (NTEST.GE.20) THEN
        if (IMK2NUM.eq.1) then
          WRITE(6,*)
     &    ' Diagonal one electron integrals, unbarred'
          CALL WRTMAT(H1A,1,NORB,1,NORB)
          WRITE(6,*)
     &    ' Diagonal one electron integrals,   barred'
          CALL WRTMAT(H1B,1,NORB,1,NORB)
          WRITE(6,*) ' Core energy ', ECORE
          IF(IHAM12.EQ.2) THEN
            WRITE(6,*) ' Coulomb and exchange integrals '
            CALL WRTMAT(RJ,NTOOB,NTOOB,NTOOB,NTOOB)
            WRITE(6,*)
            CALL WRTMAT(RK,NTOOB,NTOOB,NTOOB,NTOOB)
            if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
              WRITE(6,*)
              CALL WRTMAT(RK2,NTOOB,NTOOB,NTOOB,NTOOB)
            end if
          END IF
        end if
*
!       WRITE(6,*) ' TTSS for Blocks '
!       DO IBLOCK = 1, NBLOCK
!         WRITE(6,'(10X,4I3,2I8)') (IBLKFO(II,IBLOCK),II=1,4)
!       END DO

!       WRITE(6,*) ' Number of alpha strings , type 1 '
!       CALL IWRTMA(NSSOA,1,NSMST,1,NSMST)
!       WRITE(6,*) ' Number of beta strings , type 1 '
!       CALL IWRTMA(NSSOB,1,NSMST,1,NSMST)
*
      END IF
*
*  Diagonal elements according to Timo's formula
*
*   DIAG(IDET)
*          = sum_I (NIU*HII + NIB*HI_I_)
*            + sum_KM 1/2 * NKU*NMU [(KK|MM)     - (KM|MK)    ]
*                   + 1/2 * NKB*NMB [(K_K_|M_M_) - (K_M_|M_K_)]
*                   +       NKU*NMB [(KK|M_M_)   - (KM_|M_K)  ]
*
      IDET = 0
      ITDET = 0
*
      DO IBLK = 1, NBLOCK
#if defined (VAR_MPI2)
        JBLOCK = IBLK + ILASTBLK
        IF( NPARBLOCK( JBLOCK ) .ne. MYPROC ) GOTO 2425 
#endif
        if (NTEST.ge.50) then
          write(6,*) '-----------------------------------'
          write(6,'(A14,I8,A13)') 
     &               ' Block number ',IBLK,' in diagonal '
          write(6,*) '-----------------------------------'
        end if
*
        IATP = IBLKFO(1,IBLK)
        IBTP = IBLKFO(2,IBLK)
        IASM = IBLKFO(3,IBLK)
        IBSM = IBLKFO(4,IBLK)
        if (NTEST.ge.10)
     &    WRITE(6,'(A21,4I4)') 
     &    ' IATP IBTP IASM IBSM ',IATP,IBTP,IASM,IBSM
*
        IF(IBLTP(IASM).EQ.2) THEN
          IREST1 = 1
        ELSE
          IREST1 = 0
        END IF
*. Obtain alpha strings of sym IASM and type IATP
        CALL GETSTR_TOTSM_SPGP_REL(IASPGPTP,IATP,IASM,NAEL,NASTR1,IASTR)
*. Obtain beta strings of sym IBSM and type IBTP
        CALL GETSTR_TOTSM_SPGP_REL(IBSPGPTP,IBTP,IBSM,NBEL,NBSTR1,IBSTR)
*
        IBSTRT = 1
        IBSTOP =  NSSOB(IBSM,IBTP)
        DO IB = IBSTRT,IBSTOP
C?        WRITE(6,*)  ' IB = ', IB
          IBREL = IB - IBSTRT + 1
*
          IF(IREST1.EQ.1.AND.IATP.EQ.IBTP) THEN
            IASTRT =  IB
          ELSE
            IASTRT = 1
          END IF
          IASTOP = NSSOA(IASM,IATP)
C?        WRITE(6,*)' IASTRT,IASTOP', IASTRT,IASTOP
*
          DO IA = IASTRT,IASTOP
C?          WRITE(6,*)  ' IA = ', IA
            IDET = IDET + 1
            ITDET = ITDET + 1
*
            X = ECORE
*
            DO IEL = 1, NAEL
* Occupation of alpha spinors
              X = X + H1A(IASTR(IEL,IA))
              if (NTEST.ge.1000) then
                 write(6,'(A,I3,A,2X,1F18.10)') 
     &               'Adding H1A(',IASTR(IEL,IA),')',
     &                       H1A(IASTR(IEL,IA))
              end if
            END DO
*
            DO IEL = 1, NBEL
* Occupation of beta spinors
              X = X + H1B(IBSTR(IEL,IB))
              if (NTEST.ge.1000) then
                 write(6,'(A,I3,A,2X,1F18.10)') 
     &               'Adding H1B(',IBSTR(IEL,IB),')',
     &                       H1B(IBSTR(IEL,IB))
              end if
            END DO
*
            DO IEL = 1, NAEL
              DO JEL = 1, NAEL
                X = X + 0.5*
     &            RJ(IASTR(IEL,IA),IASTR(JEL,IA))
     &                - 0.5*
     &            RK(IASTR(IEL,IA),IASTR(JEL,IA))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding 0.5 J(',
     &                          IASTR(IEL,IA),
     &                          IASTR(JEL,IA),
     &                          ') = ',
     &            0.5*RJ(IASTR(IEL,IA),IASTR(JEL,IA)),
     &                          '- 0.5 K (',
     &                          IASTR(IEL,IA),
     &                          IASTR(JEL,IA),') =',
     &            0.5*RK(IASTR(IEL,IA),IASTR(JEL,IA))
                  end if
              END DO
            END DO
*
            DO IEL = 1, NBEL
              DO JEL = 1, NBEL
               X = X + 0.5*
     &            RJ(IBSTR(IEL,IB),IBSTR(JEL,IB))
     &               - 0.5*
     &            RK(IBSTR(IEL,IB),IBSTR(JEL,IB))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding 0.5 J(',
     &                          IBSTR(IEL,IB),
     &                          IBSTR(JEL,IB),
     &                          ') = ',
     &            0.5*RJ(IBSTR(IEL,IB),IBSTR(JEL,IB)),
     &                          '- 0.5 K (',
     &                          IBSTR(IEL,IB),
     &                          IBSTR(JEL,IB),') =',
     &            0.5*RK(IBSTR(IEL,IB),IBSTR(JEL,IB))
                  end if
              END DO
            END DO
*
            DO IEL = 1, NAEL
              DO JEL = 1, NBEL
                X = X + 
     &            RJ(IASTR(IEL,IA),IBSTR(JEL,IB))
                if (IDENSI.ge.1.or.ISPINFREE.eq.0) then
                  X = X - 
     &            RK2(IASTR(IEL,IA),IBSTR(JEL,IB))
                  if (NTEST.ge.1000) then
                     write(6,'(2X,A,2I3,A,1F18.10,3X,A,2I3,A,1F18.10)') 
     &                          'Adding J(',
     &                          IASTR(IEL,IA),
     &                          -IBSTR(JEL,IB),
     &                          ') = ',
     &            RJ(IASTR(IEL,IA),IBSTR(JEL,IB)),
     &                          '- K2 (',
     &                          IASTR(IEL,IA),
     &                          -IBSTR(JEL,IB),') =',
     &            RK2(IASTR(IEL,IA),IBSTR(JEL,IB))
                  end if
                end if
              END DO
            END DO
*
            if (NTEST.ge.100) then
               write(6,*) 'Element ITDET and value X:'
               write(6,*) ITDET,X
               write(6,*) 'offset IDET:'
               write(6,*) IDET
               WRITE(6,*) ' alpha and beta string : '
               CALL IWRTMA(IASTR(1,IA),1,NAEL,1,NAEL)
               CALL IWRTMA(IBSTR(1,IB),1,NBEL,1,NBEL)
               WRITE(6,*)
            end if
            DIAG(IDET) = X
            IF(IB.EQ.IA) DIAG(IDET) = DIAG(IDET) + XADD
          END DO
*         ^ End of loop over alpha strings
        END DO
*       ^ End of loop over betastrings
*. Yet a RAS block of the diagonal has been constructed
2345    continue
        IF(NTEST.GE.100) THEN
          if(idet .gt. 0)then
            write(6,*) ' number of diagonal elements to disc ',IDET
            CALL WRTMAT(DIAG,1,IDET,1,IDET)
          end if
        END IF
#if defined (VAR_MPI2)
        call interface_mpi_FILE_WRITE_AT_r(IDIA,IDIA_OFFSET,DIAG,IDET,
     &                         ISTAT)
*       new offset
        IDIA_OFFSET = IDIA_OFFSET + IDET
#endif
        IF( .NOT. NO_CDIAF )THEN
          CALL ITODS(IDET,1,-1,LUDIA)
          CALL TODSC(DIAG,IDET,-1,LUDIA)
        END IF
        IDET = 0
2425    CONTINUE
      END DO ! loop over TTSS blocks

      if (NTEST.ge.500) then
         WRITE(6,*)
     &   '  Total number of diagonal elements:',ITDET
      end if
      call qexit('GASDIAS')
 
      END
***********************************************************************

      SUBROUTINE GASDIAT_REL(DIAG,LUDIA,ECORE,ICISTR,IHAM12,
     &           IBLTP,NBLOCK,IBLKFO,T_BUFF,T_BUFF_D,
     &           IDENSI,NIRREP,ISPINFREE,ENVIRO,SIGDEN_ROUTE,IRECO,
     &           NPARBLOCK,CPIC_PASS_D)
*
* CI diagonal in SD basis for state with symmetry ISM in internal
* space ISPC
*
* GAS version, Winter of 95
*
* Driven by table of TTS blocks, May97
*
* For DBG-RELCI, March21 98
* Modifed fetch of 1 electron integrals, May 98
* Modified for parallel runs, S. Knecht, May 2007
*
      use memory_allocator
      use mospinor_info
      use symmetry_setup_krci
      use interface_to_mpi
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
* =====
*.Input
* =====
*
*./ORBINP/ : NACOB used
*
#include "mxpdim.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "wrkspc.inc"
#include "cprnt.inc"
#include "cgas.inc"
#include "gasstr.inc"
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*
      LOGICAL CPIC_PASS_D
      DIMENSION IBLTP(*), T_BUFF_D(*)
      DIMENSION IBLKFO(8,NBLOCK), NPARBLOCK(*)
#if defined (VAR_MPI2)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IDIA_OFFSET
#endif
      integer, allocatable :: astring(:,:)
      integer, allocatable :: bstring(:,:)
*
* ======
*.Output
* ======
      DIMENSION DIAG(*)
*
      CALL QENTER('CIDIA')
*
      NTEST = 000
      NTEST = MAX(NTEST,IPRDIA/1000)
*. Largest Number of electrons around
      MAXNEL = 0
      DO IMK2 = 1, NMS2VAL
        IATP = IST_FOR_DT(1,IMK2)
        IBTP = IST_FOR_DT(2,IMK2)
        MAXNEL = MAX(MAXNEL,NELEC(IATP),NELEC(IBTP))
      END DO
*. Local memory
      IDUM = 0
      call memmar(KDUM,  IDUM,    'MARK  ',IDUM,'GASDIA')
*
CTF  Total number of orbitals and active orbitals (spinors, actually)
CTF  from BOTH sets, barred and unbarred !
      NTOOBTOT = NTOOB + NTOOB2
      NACOBTOT = NACOB + NACOB2
*
*
      call memmar(KLJ   ,NTOOB**2,'ADDL  ',2,'KLJ   ')
      call memmar(KLK   ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLK2  ,NTOOB**2,'ADDL  ',2,'KLK   ')
      call memmar(KLSCR2,2*NTOOBTOT**2,'ADDL  ',2,'KLSC2 ')
      call memmar(KLXA  ,NACOBTOT,   'ADDL  ',2,'KLXA  ')
      call memmar(KLXB  ,NACOBTOT,   'ADDL  ',2,'KLXB  ')
      call memmar(KLSCR ,2*NACOBTOT, 'ADDL  ',2,'KLSCR ')
C     CALL MEMMAR(KLH1D ,NACOBTOT, 'ADDL  ',2,'KLH1D ')
      call memmar(KLH1DA,NTOOB, 'ADDL  ',2,'KLH1DA')
      call memmar(KLH1DB,NTOOB, 'ADDL  ',2,'KLH1DB')
*. Space for blocks of strings
      call alloc(astring,MAXNEL,MXNSTR)
      call alloc(bstring,MAXNEL,MXNSTR)
      KLSVST = 1
      call memmar(KLRJKA,MXNSTR, 'ADDL  ',2,'KLRJKA')
!
!     initialize diagonal 1e-integral arrays
      call dzero(work(klh1da),ntoob)
      call dzero(work(klh1db),ntoob)
*
      if (IRECO.eq.1) then
        if( CPIC_PASS_D )then
C
          IOFF_T_D = 1
C
          CALL DCOPY(NTOOB,T_BUFF_D,1,WORK(KLH1DA),1)
          IOFF_T_D = IOFF_T_D + NTOOB
          CALL DCOPY(NTOOB,T_BUFF_D(IOFF_T_D),1,WORK(KLH1DB),1)
          IOFF_T_D = IOFF_T_D + NTOOB
          CALL DCOPY(NTOOB**2,T_BUFF_D(IOFF_T_D),1,WORK(KLJ),1)
          IOFF_T_D = IOFF_T_D + NTOOB**2
          CALL DCOPY(NTOOB**2,T_BUFF_D(IOFF_T_D),1,WORK(KLK),1)
          IOFF_T_D = IOFF_T_D + NTOOB**2
          CALL DCOPY(NTOOB**2,T_BUFF_D(IOFF_T_D),1,WORK(KLK2),1)
C
CCC       WRITE(LUWRT,*) ' T_BUFF_D:'
CCC       CALL WRTMATMN8(T_BUFF_D,1,MY_T_LEN_D,1,MY_T_LEN_D,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLH1DA):'
CCC       CALL WRTMATMN(WORK(KLH1DA),1,NTOOB,1,NTOOB,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLH1DB):'
CCC       CALL WRTMATMN(WORK(KLH1DB),1,NTOOB,1, NTOOB,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLJ):'
CCC       CALL WRTMATMN(WORK(KLJ),1, NTOOB**2,1, NTOOB**2,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLK):'
CCC       CALL WRTMATMN(WORK(KLK),1, NTOOB**2,1, NTOOB**2,LUWRT)
CCC       WRITE(LUWRT,*) ' WORK(KLK2):'
CCC       CALL WRTMATMN(WORK(KLK2),1, NTOOB**2,1, NTOOB**2,LUWRT)
C
        else
          if (IHAM12.eq.1.or.IHAM12.eq.2) then
            call geth1dia4(WORK(KLH1DA),WORK(KLH1DB),T_BUFF)
          end if
*         obtain 2e integrals needed for diagonal
          if (IHAM12.EQ.2.or.IHAM12.eq.3) then
            call get_dia_2e(WORK(KLJ),WORK(KLK),WORK(KLK2),T_BUFF,
     &                      IDENSI,ISPINFREE,NGAS,IPRDIA)
          end if
        end if
C       ^  IT_SHL?
      end if
C
#if defined (VAR_MPI2)
      IDIA_OFFSET = 0
*
*     calculate offset IDIA_OFFSET for LUDIA file
*
      IDIA_OFFSET = IDIA_OFFSET + MY_DIA_OFF + 
     &              MY_VEC1_IOFF  * ( IRECO - 1 )
CSK      WRITE(LUWRT,*)' Start writing on LUDIA at',IDIA_OFFSET
*
#endif
*
      DO IMK2 = 1, NMS2VAL
        if (NTEST.ge.10) then
          write(LUWRT,*) '=============================='
          write(LUWRT,*) 'Treating MK space in diagonal '
          write(LUWRT,'(A15,I4)') '         MK2 = ',MS2VAL(IMK2)
          write(LUWRT,*) '=============================='
        end if
        IF(IMK2.EQ.1) THEN
          IBOFF = 1
          ILASTBLK = 0
        ELSE
          IBOFF = IBOFF + NBLK_MS2(IMK2-1)
          ILASTBLK = ILASTBLK + NBLK_MS2(IMK2-1)
        END IF
        NNBLK = NBLK_MS2(IMK2)
CSK        WRITE(LUWRT,*)' IBOFF, NNBLK, ILASTBLK',IBOFF, NNBLK,ILASTBLK
*
        IATP = IST_FOR_DT(1,IMK2)
        IBTP = IST_FOR_DT(2,IMK2)
*
        NAEL = NELEC(IATP)
        NBEL = NELEC(IBTP)
*
        NOCTPA = NOCTYP(IATP)
        NOCTPB = NOCTYP(IBTP)
*
        IASPGPTP = IATP
        IBSPGPTP = IBTP
*
*. Offsets for alpha and beta supergroups
        IOCTPA = IBSPGPFTP(IATP)
        IOCTPB = IBSPGPFTP(IBTP)
        IF(NTEST.GE.10) THEN
          WRITE(LUWRT,*) ' ================'
          WRITE(LUWRT,*) ' GASDIA speaking '
          WRITE(LUWRT,*) ' ================'
          WRITE(LUWRT,*) ' IATP IBTP NAEL NBEL ',IATP,IBTP,NAEL,NBEL
          write(LUWRT,*) ' NOCTPA NOCTPB  : ', NOCTPA,NOCTPB
          write(LUWRT,*) ' IOCTPA IOCTPB  : ', IOCTPA,IOCTPB
        END IF
*
** Info on block structure of internal state
*
**. Diagonal of one-body integrals and coulomb and exchange integrals
*
*
        CALL GASDIAS_REL(IMK2,NMS2VAL,MS2VAL,NAEL,astring,NBEL,
     &                   bstring,NTOOB,DIAG,NSMST,NIRREP,
     &                   WORK(KLH1DA),WORK(KLH1DB),
     &                   WORK(KLXA),WORK(KLXB),WORK(KLSCR),WORK(KLJ),
     &                   WORK(KLK),WORK(KLK2),
     &                   WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &                   LUDIA,ECORE,PLSIGN,PSSIGN,IPRDIA,NTOOB,ICISTR,
     &                   WORK(KLRJKA),IHAM12,IBLTP,NNBLK,
     &                   IBLKFO(1,IBOFF),IASPGPTP,IBSPGPTP,ISPINFREE,
     &                   NPARBLOCK,ILASTBLK,IDIA_OFFSET)
*
      END DO
*     ^End of loop over MK2 spaces
*
*.Flush local memory
      call dealloc(astring)
      call dealloc(bstring)
      call memmar(KDUM,  IDUM,    'FLUSM ',IDUM,'GASDIA')
      CALL QEXIT('CIDIA')
*
      END
***********************************************************************

      SUBROUTINE GET_DIA_2E(UUUUC,UUUUE,UBBUE,T_BUFF,IDENSI,ISPINFREE,
     &                      NGAS,IPRNT)
*
* Obtain 2electron integrals needed for
* constructing the determinant diagonal
*
*===============================================
* Jeppe Olsen, Timo Fleig
*  SIGDEN route implemented, 
*     Timo Fleig, August 2001
*===============================================
*
* UUUUC(I,J) = ( I I  ! J  J)
* UUUUE(I,J) = ( I J  ! J  I)
* UBBUE(I,J) = ( I J_ ! J_ I)
*
* ISPINFREE == 0 --> relativistic run (including SOC)
* ISPINFREE == 1 --> relativistic run (no        SOC)
*
      use symmetry_setup_krci
      use mospinor_info
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
      DIMENSION UUUUC(NTOOB,NTOOB),UUUUE(NTOOB,NTOOB)
      DIMENSION UBBUE(NTOOB,NTOOB)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "glbbas.inc"
#include "ctcc.inc"
#include "parluci.h"
#include "integrals_off.inc"
*
      dimension INDGATH(MXPORB),IKPERM(4)
      DIMENSION T_BUFF(*)
*     Input list of integrals (or in SIGDEN on WORK(KT_CC))
*     Input list of integrals (or in SIGDEN in "shared memory"
*     mode on T_BUFF (for "node master" == N_MASTER)
*
!     NTESTL = 120 ! debug
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT/1000)
*
      if(NTEST.ge.2) then
        write(6,*) ' //////////////////////////////////// '
        write(6,*) '  Entering get_dia_2e                 '
        write(6,*) ' //////////////////////////////////// '
      end if
*
      NINT12 = N1ELINT + N2ELINT
*
* UUUUC and UUUUE obtained from integral list 3
*
* ( I I ! J J )
*
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal Coulomb ints. (II|JJ):'
      end if
*
      IKFAC = 1
      IKPERM(1) = 1
      IKPERM(2) = 2
      IKPERM(3) = 3
      IKPERM(4) = 4
      if(ISPINFREE.eq.0) then
* Storage mode for integrals contains "-" sign for beta/barred indices.
        IKFAC = -1
* Use Kramers symmetry on second index pair.
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 4
        IKPERM(4) = 3
      end if  
      IGSOFF = 0
      JGSOFF = 0
      ISMOFF = 0
      JSMOFF = 0
      do IGS = 1,NGAS,1
        if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
        do JGS = 1,NGAS,1
          if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
          do ISM = NSMOB/2+1, NSMOB
            if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
            do JSM = NSMOB/2+1, NSMOB
              if (JSM.gt.NSMOB/2+1) JSMOFF = JSMOFF + NGSOB(JSM-1,JGS)
              NI = NGSOB(ISM,IGS)
              NJ = NGSOB(JSM,JGS)
              IOFF = IGSOFF + ISMOFF
              JOFF = JGSOFF + JSMOFF
              if (NTEST.ge.10) then
                write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
              end if
              do I = 1+IOFF,NI+IOFF,1
                do J = 1+JOFF,NJ+JOFF,1
                  call geth2dia_sigden(UUUUC(I,J),
     &                                 WORK(KT_CC),T_BUFF,
     &                                 I,J,
     &                                 'COUL',
     &                                 I2DIATP,IKFAC,IKPERM,
     &                                 IREOTS,IREOTS,
     &                                 imosp_luci2dirac1,
     &                                 imosp_luci2dirac1,
     &                                 WORK(KMDIASV),              
     &                                 N1DIAINT,NDIAINT,NTEST)
                end do
              end do
            end do
            JSMOFF = 0
          end do
          ISMOFF = 0
        end do
        JGSOFF = 0
      end do
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,'(/a )') '  Two-electron integrals for CI diagonal'
      WRITE(6,'( a )') '  --------------------------------------'
      WRITE(6,'(/a/)') '            (I I | J J) matrix '
      CALL WRTMAT(UUUUC,NTOOB,NTOOB,NTOOB,NTOOB)
#endif
!#undef LUCI_DEBUG
*
*
* ( I J ! J I )
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal exchange ints. (IJ|JI):'
      end if
*
      IKFAC = 1
      IKPERM(1) = 1
      IKPERM(2) = 2
      IKPERM(3) = 3
      IKPERM(4) = 4
      if(ISPINFREE.eq.0) then
        IKFAC = -1
* Use Kramers symmetry on second index pair.
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 4
        IKPERM(4) = 3
      end if  
      IGSOFF = 0
      JGSOFF = 0
      ISMOFF = 0
      JSMOFF = 0
      do IGS = 1,NGAS,1
        if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
        do JGS = 1,NGAS,1
          if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT(JGS-1)
          do ISM = NSMOB/2+1, NSMOB
            if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
            do JSM = NSMOB/2+1, NSMOB
              if (JSM.gt.NSMOB/2+1) JSMOFF = JSMOFF + NGSOB(JSM-1,JGS)
              NI = NGSOB(ISM,IGS)
              NJ = NGSOB(JSM,JGS)
              IOFF = IGSOFF + ISMOFF
              JOFF = JGSOFF + JSMOFF
              if (NTEST.ge.10) then
                write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
              end if
              do I = 1+IOFF,NI+IOFF,1
                do J = 1+JOFF,NJ+JOFF,1
                  call geth2dia_sigden(UUUUE(I,J),
     &                                 WORK(KT_CC),T_BUFF,
     &                                 I,J,
     &                                 'EXC1',
     &                                 I2DIATP,IKFAC,IKPERM,
     &                                 IREOTS,IREOTS,
     &                                 imosp_luci2dirac1,
     &                                 imosp_luci2dirac1,
     &                                 WORK(KMDIASV),              
     &                                 N1DIAINT,NDIAINT,NTEST)
                end do
              end do
            end do
            JSMOFF = 0
          end do
          ISMOFF = 0
        end do
        JGSOFF = 0
      end do
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,'(/a )') '  Two-electron integrals for CI diagonal'
      WRITE(6,'( a )') '  --------------------------------------'
      WRITE(6,'(/a/)') '            (I J | J I) matrix '
      CALL WRTMAT(UUUUE,NTOOB,NTOOB,NTOOB,NTOOB)
#endif
!#undef LUCI_DEBUG
*
* ( I J_ ! J_ I)
      if (NTEST.ge.2) then
        write(6,*)
        write(6,*) ' Getting diagonal exchange ints. (IJ_|J_I):'
      end if
*
*( obtained from list 3 )
      if(IDENSI.ge.1.or.ISPINFREE.eq.0)then
        IKFAC     = -1
        IKPERM(1) = 1
        IKPERM(2) = 2
        IKPERM(3) = 3
        IKPERM(4) = 4

        IGSOFF    = 0
        JGSOFF    = 0
        ISMOFF    = 0
        JSMOFF    = 0
        do IGS = 1,NGAS,1
          if (IGS.gt.1) IGSOFF = IGSOFF + NGSOBT(IGS-1)
          do JGS = 1,NGAS,1
            if (JGS.gt.1) JGSOFF = JGSOFF + NGSOBT2(JGS-1)
            do ISM = NSMOB/2+1, NSMOB
              if (ISM.gt.NSMOB/2+1) ISMOFF = ISMOFF + NGSOB(ISM-1,IGS)
              do JSM = NSMOB/2+1, NSMOB
                if (JSM.gt.NSMOB/2+1) 
     &                     JSMOFF = JSMOFF + NGSOB2(JSM-1,JGS)
                NI = NGSOB(ISM,IGS)
                NJ = NGSOB2(JSM,JGS)
                IOFF = IGSOFF + ISMOFF
                JOFF = JGSOFF + JSMOFF
                if (NTEST.ge.10) then
                  write(6,*) 'IGS,JGS,ISM,JSM ',IGS,JGS,ISM,JSM
                  write(6,*) ' Index range :',1+IOFF,'...',NI+IOFF
                  write(6,*) ' Index range :',1+JOFF,'...',NJ+JOFF
                end if
                do I = 1+IOFF,NI+IOFF,1
                  do J = 1+JOFF,NJ+JOFF,1
                    call geth2dia_sigden(UBBUE(I,J),
     &                                 WORK(KT_CC),T_BUFF,
     &                                 I,J,
     &                                 'EXC2',
     &                                 I2DIATP,IKFAC,IKPERM,
     &                                 IREOTS,IREOTS2,
     &                                 imosp_luci2dirac1,
     &                                 imosp_luci2dirac2,
     &                                 WORK(KMDIASV),              
     &                                 N1DIAINT,NDIAINT,NTEST)
                  end do
                end do
              end do
              JSMOFF = 0
            end do
            ISMOFF = 0
          end do
          JGSOFF = 0
        end do
      end if ! relativistic run or scalar-relativistic run
*
!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(6,'(/a )') '  Two-electron integrals for CI diagonal'
      WRITE(6,'( a )') '  --------------------------------------'
      if(ISPINFREE.eq.0)then
        WRITE(6,'(/a/)') '            (I J_| J_I) matrix '
        CALL WRTMAT(UBBUE,NTOOB,NTOOB,NTOOB,NTOOB)
      end if
#endif
!#undef LUCI_DEBUG

      END
***********************************************************************

      SUBROUTINE GETH1DIA4(H1DIA_U,H1DIA_B,T_BUFF)
*
* Fetch unbarred and barred diagonal one-electron integrals
*
      use symmetry_setup_krci
      use mospinor_info
      use memory_allocator
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
*. General input
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "ctcc.inc"
#include "parluci.h"
#include "integrals_off.inc"
*
*     Input list of integrals (or in SIGDEN on WORK(KT_CC))
*     Input list of integrals (or in SIGDEN in "shared memory"
*     mode on T_BUFF (for "node master")  
*. Output
      DIMENSION H1DIA_U(*),H1DIA_B(*),T_BUFF(*)
!     scratch
      real(8), allocatable :: h1dia_local(:)
      integer, allocatable :: ngsobx(:,:)
      integer, allocatable :: ngsobtx(:)
      integer, allocatable :: ireotsx(:)
      integer, allocatable :: imosp_luci2diracx(:)

!     NTESTL = 0100 ! debug
      NTESTL = 0000
      NTEST = MAX(NTESTL,IPRDIA/1000)
*
      if(NTEST.gt.1)then
        write(LUWRT,*) '----------------------------------------'
        write(LUWRT,*) ' GETH1DIA4 speaking '
        write(LUWRT,*) '----------------------------------------'
      end if
*
      call alloc(ngsobtx,         mxpngas, id ='ngsobtx_h1dia')
      call alloc(ngsobx, mxndgirr,mxpngas, id ='ngsobx_h1dia')
      call alloc(ireotsx        ,mxporb,   id ='ireotsx_h1dia')
      call alloc(imosp_luci2diracx, mxporb,id ='luci2dirac   ')
      call alloc(h1dia_local   ,ntoob**2,  id ='h1dia_local')
      ngsobtx           = 0
      ngsobx            = 0
      ireotsx           = 0
      imosp_luci2diracx = 0
      
      do ispin = 1, 2
 
        if(ispin .eq. 1)then
          call icopy(         mxpngas,ngsobt ,1,ngsobtx,1)
          call icopy(mxndgirr*mxpngas,ngsob  ,1,ngsobx,1)
          call icopy(mxporb,ireots ,1,ireotsx,1)
          call icopy(mxporb,imosp_luci2dirac1 ,1,imosp_luci2diracx,1)
        else
          call icopy(         mxpngas,ngsobt2,1,ngsobtx,1)
          call icopy(mxndgirr*mxpngas,ngsob2 ,1,ngsobx,1)
          call icopy(mxporb     ,ireots2,1,ireotsx,1)
          call icopy(mxporb,imosp_luci2dirac2,1,imosp_luci2diracx,1)
        end if

!       GAS-SYMM ordering (GAS is outer loop!)

!       initialize
        call dzero(h1dia_local,ntoob**2)
        IOFFSD  = 1
        ISMOFF  = 0
        IGASOFF = 0
        ISTART  = 0
        IEND    = 0
        NI      = 0
        do IGAS = 1, NGAS
          if(IGAS.gt.1) IGASOFF = IGASOFF + NGSOBTx(IGAS-1)
!         print *,'IGASOFF ==> ',IGASOFF 
          ISMOFF  = 0
          do ISMOB = nirr_dg/2+1, nirr_dg
            if(ISMOB.gt.nirr_dg/2+1) 
     &      ISMOFF = ISMOFF + NGSOBx(ISMOB-1,IGAS)
            NI     = NGSOBx(ISMOB,IGAS)
            ISTART = IGASOFF + ISMOFF + 1
            IEND   = IGASOFF + ISMOFF + NI
            if (NTEST.ge.50) then
              write(LUWRT,'(A,5I4)') '  IGAS,ISMOB,ISTART,IEND',
     &                                  IGAS,ISMOB,ISTART,IEND
            end if
!           check for occupation of current space
            if(MXNELGS(IGAS).ne.0)then
!             check for spinors in current space and symmetry irrep
              if(NGSOBx(ISMOB,IGAS).gt.0)then
                call geth1dia_sigden(ispin,h1dia_local,WORK(KT_CC),
     &                               T_BUFF,WORK(KMDIASV),N1ELINT,
     &                               N1DIAINT,NGSOBx(ISMOB,IGAS),
     &                               ISTART,IEND,ireotsx,
     &                               imosp_luci2diracx,NTEST)
                if(ispin .eq. 1)then
                  call dcopy(NGSOBx(ISMOB,IGAS),h1dia_local,1,
     &                       H1DIA_U(IOFFSD),1)
                else
                  call dcopy(NGSOBx(ISMOB,IGAS),h1dia_local,1,
     &                       H1DIA_B(IOFFSD),1)
                end if
              endif
            end if
            IOFFSD = IOFFSD + NGSOBx(ISMOB,IGAS)
          end do
        end do
      end do ! unbarred/barred (ispin loop)
 
      call dealloc(h1dia_local)
      call dealloc(imosp_luci2diracx)
      call dealloc(ireotsx)
      call dealloc(NGSOBx)
      call dealloc(NGSOBTx)


!#define LUCI_DEBUG
#ifdef LUCI_DEBUG
      WRITE(LUWRT,'(/a)') ' diagonal integrals'
      WRITE(LUWRT,'( a)') ' ------------------'
      WRITE(LUWRT,'(/a/)') ' unbarred spinors'
      CALL WRTMATMN(H1DIA_U,1,NTOOB,1,NTOOB,LUWRT)
      WRITE(LUWRT,'(/a/)') '   barred spinors'
      CALL WRTMATMN(H1DIA_B,1,NTOOB,1,NTOOB,LUWRT)
#endif
!#undef LUCI_DEBUG

      END
***********************************************************************

      subroutine geth1dia_sigden(IUB,XOUT,XINT,T_BUFF,
     &                           MDIASV,N1E,N1DIAINT,NUMINT,
     &                           ISTART,IEND,IREOTS,
     &                           imosp_luci2diracx,IPRNT)
!
! Get diagonal one-electron integrals in SIGDEN formalism.
!
!
#include "implicit.inc"
#include "mxpdim.inc"
#include "integrals_off.inc"
#include "parluci.h"
*
      DIMENSION IREOTS(*),MDIASV(4,*), imosp_luci2diracx(*)
      DIMENSION XINT(N1E),XOUT(*),T_BUFF(*)
*
!     NTESTL = 100 ! debug
      NTESTL = 000
      NTEST = max(NTESTL,IPRNT/1000)
*
      is_zero = 0
      ICT = 1
      do ISEA = ISTART, IEND
*
* Rename to SYMM-TYPE ordering on external list:
        ISEA_REO = imosp_luci2diracx(IREOTS(ISEA))
        do INT=1,N1DIAINT,1
          if (MDIASV(IUB,INT).eq.ISEA_REO.and.MDIASV(IUB+2,INT).eq.
     &        ISEA_REO.and.MDIASV(IUB+1,INT).eq.is_zero) then
C
            if (NTEST.ge.100) then
              write(LUWRT,'(A,I5,A,I5,A,1F18.12)') 
     &                'Writing integral no. ',IPOSDIA(INT),'to ',ICT,
     &                '  Integral ',XINT(IPOSDIA(INT))
            end if
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
C           write(LUWRT,*) ' fetch integral from pos',IPOSDIA(INT)
            XOUT(ICT) = T_BUFF(IPOSDIA(INT))
#else
            if( lowsrt_ijkl )then
              XOUT(ICT) = T_BUFF(IPOSDIA(INT))
            else
              XOUT(ICT) = XINT(IPOSDIA(INT))
            end if
#endif
            ICT = ICT + 1
            goto 10
          end if
        end do
        write(LUWRT,*) 'No integral found in geth1dia_sigden.'
        write(LUWRT,*) 'Index pair ',ISEA_REO,ISEA_REO
        write(LUWRT,*) 'original index pair ',ISEA,ISEA
        call quit('*** error in geth1dia_sigden: No integral found ***')
  10    continue
      end do
      if ((ICT-1).ne.NUMINT) then
        write(LUWRT,*) 'Number of processed integrals: ',ICT-1
        write(LUWRT,*) 'Number of diagonal elements  : ',NUMINT
        write(LUWRT,*) 'Problem in geth1dia_sigden.'
        call quit( '*** error in geth1dia_sigden: # of processed diag'//
     &  'onal 1e- integrals < total # of diagonal 1e- integrals.***')
      end if

      end
***********************************************************************

      subroutine geth2dia_sigden(XOUT,XINT,T_BUFF,ISEA,JSEA,TYPE,
     &                           I2DIATP,IKFAC,IKPERM,
     &                           IREOTS,IREOTS2,
     &                           imosp_luci2diracx1,imosp_luci2diracx2,
     &                           MDIASV,N1DIAINT,NDIAINT,IPRNT)
*
*     Get diagonal two-electron integrals in SIGDEN formalism.
*     Timo Fleig, August 2001
*
*     Modified to account for integrals over true spinors, i.e.
*     explicit reference to unbarred (+) and barred (-) indices (IKFAC).
*
*     Revised routine - Timo Fleig, February 2002
*
*     "shared memory" MPI implementation by Stefan Knecht - Nov. 2007
*     linear symmetry generalization by Stefan Knecht - March 2011
*
#include "implicit.inc"
#include "mxpdim.inc"
#include "integrals_off.inc"
#include "parluci.h"
*
      character*4 TYPE
*
      dimension IKPERM(4),IREOTS(*),IREOTS2(*),MDIASV(4,*)
      dimension imosp_luci2diracx1(*), imosp_luci2diracx2(*)
      dimension XINT(NINT12), T_BUFF(*)
*
      NTESTL = 000
      NTEST  = max(NTESTL,IPRNT/1000)
*
      IERROR = 1
      do INT = N1DIAINT+1,NDIAINT,1
        if (TYPE.eq.'COUL') then
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(2),INT)
          K = MDIASV(IKPERM(3),INT) * IKFAC
          L = MDIASV(IKPERM(4),INT) * IKFAC
        else if (TYPE.eq.'EXC1') then
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(4),INT) * IKFAC
          K = MDIASV(IKPERM(2),INT)
          L = MDIASV(IKPERM(3),INT) * IKFAC
        else if (TYPE.eq.'EXC2') then
          I = MDIASV(IKPERM(1),INT)
          J = MDIASV(IKPERM(4),INT)
          K = MDIASV(IKPERM(2),INT) * IKFAC
          L = MDIASV(IKPERM(3),INT) * IKFAC
        end if
*
* Rename: External ST ordering, internal TS ordering:
        if (I.eq.imosp_luci2diracx1(IREOTS(ISEA)).and.
     &      J.eq.imosp_luci2diracx1(IREOTS(ISEA)).and.
     &      K.eq.imosp_luci2diracx2(IREOTS2(JSEA)).and.
     &      L.eq.imosp_luci2diracx2(IREOTS2(JSEA)))then
C
          if (NTEST.ge.100) then
            write(LUWRT,'(A,1X,A,1I4,3X,4I3,1F18.10)') 
     &            'Diagonal int type ',TYPE,
     &             IPOSDIA(INT),I,J,K,L,
     &             XINT(IPOSDIA(INT))
          end if
#if defined (ARCH32BIT) && defined (SYS_AIX) && defined (VAR_MPI2)
          XOUT = T_BUFF(IPOSDIA(INT))
#else
          if( lowsrt_ijkl )then 
            XOUT = T_BUFF(IPOSDIA(INT))
          else
            XOUT = XINT(IPOSDIA(INT))
          end if
#endif
          IERROR = 0
* and exit.
          goto 10
        end if
      end do
      if (IERROR.eq.1) then
        write(6,*) 'No integral found in geth2dia_sigden.'
        write(6,*) 'Type of integral = ',TYPE
        write(6,*) 'Index pair ',IREOTS(ISEA),IREOTS2(JSEA)
        write(6,*) 'Indices ',I,J,K,L
        call quit(' *** error in geth2dia_sigden: no 2e- diagonal int'//
     &            'egral found. ***')
      end if
  10  continue

      end
