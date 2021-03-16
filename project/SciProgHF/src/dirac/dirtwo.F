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
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twofck */
      SUBROUTINE TWOFCK(ISYMOP,IHRMOP,IFCKOP,FMAT,DMAT,NFMAT,
     &                  NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for the construction of the two-electron Fock matrix 
C     
C
C     Written by T.Saue - January 1995
C     
C     The reduced density matrices are generated within this routine,
C     whereas Cauchy-Schwarz integrals must be provided externally
C
C*****************************************************************************
      use quaternion_algebra
#ifdef MOD_SRDFT
      use srdft_cfg ! for MCSCF-srDFT calculation
#endif
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "dummy.h"
      PARAMETER(D0=0.0D0, D1=1.0D0)
C
#include "dcbind.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "ccom.h"
#include "dcborb.h"
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),
     &          DMAT(N2BBASX,NZ,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),
     &          IFCKOP(NFMAT),NPOS(*),WORK(LWORK)
C
      CALL QENTER('TWOFCK')
#include "memint.h"
      iprint_save = iprint
!     iprint = 5
C
      IF(IPRINT.GE.5) THEN
        CALL TITLER('Output from TWOFCK','*',103)
      ENDIF
C 
C     In approximate density ZORA only the LL part of the density matrix
C     is to be used
C
      IF (ZORA.AND..NOT.ZORA4) THEN
         DO I = 1,NFMAT
            CALL KEEPLL(DMAT(1,1,I),NZ,'SRT')
         ENDDO
      ENDIF
C
C     Transform density matrix to unsorted basis
C     and insert half-phases
C     ==========================================
C
      DO I = 1,NFMAT
        IREP = ISYMOP(I)-1
        IF(IPRINT.GE.5) THEN
          CALL HEADER('TWOFCK: Unsorted DMAT',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(DMAT(1,1,I),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NTBAS(0),NZ,
     &                IPQTOQ(1,IREP),LUPRI)
        ENDIF
C
C       Insert half-phases
C       ==================
C
        IF(NZ.LT.4) THEN
          DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            CALL Q2BPHASE('D',IQ,1,DMAT(1,IZ,I))
          ENDDO
          IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted DMAT with half-phases inserted',-1)
            WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
            CALL PRQMAT(DMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NZ,
     &                  IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
        IF(.NOT.MDIRAC)CALL BSTOBU(DMAT(1,1,I),NZ,WORK,LWORK)
      ENDDO
C
C     ===========================
      IF (ONECAP.AND.SMLV1C) THEN
C     ===========================
C
C        Adapt density matrix in the case that we are going to neglect
C        multicenter integrals in the VSS part. Save original density
C        matrix to file as it should be unchanged once we leave the
C        routine. Use Fock matrix to store projector.
C        TODO : reduce memory use for DMAT in AO basis
C
         CALL MEMGET2 ('REAL','BUF1',KBUF1,4*N2BBASX,WORK,KFREE,LFREE)
         CALL OPNFIL(LUPMAT,'DFPMAT','OLD','TWOFCK')
         CALL READT (LUPMAT,N2BBASX,FMAT)
         CALL WRITT (LUPMAT,N2BBASXQ*NFMAT,DMAT)
         CLOSE(LUPMAT,STATUS='KEEP')
         DO I = 1,NFMAT
            IREP = ISYMOP(I)-1
            CALL PR1CEX3 (NTBAS(0),IREP,DMAT(1,1,I),
     &                    FMAT,WORK(KBUF1),IPRINT)
         ENDDO
         CALL MEMREL('TWOFCK',WORK,KBUF1,KBUF1,KFREE,LFREE)
C     =====
      ENDIF
C     =====
C
C     Initialize and build Fock matrix
C     ================================
C
      CALL DZERO(FMAT,N2BBASXQ*NFMAT)
      IF (HFXMU.NE.D0) THEN
C
C       Calculate attenuated exchange (hfxatt*erf(hfxmu*r12)/r12) 
C       NOTE that this extra part has to be calculated first,
C       since there is a scaling of Fock matrices in the SO Fock build
C       The scaling is now turned off for HFXMU.NE.0
C
        HFXTMP = HFXFAC
        HFXFAC = HFXATT 
C       ...note that when HFXMU.NE.0, Coulomb contributions are
C       turned off !
C
C       Manu: note also that for MCSCF-srDFT BOTH Coulomb and exchange 
C             contributions should be calculated since long-range and
C             short-range coulomb terms cannot be recombined
C             straightforwardly due to the active density matrix
C             contribution !!!
C
C       Trond:  Gaunt-integrals should be turned off !
C
        CALL FOCKBUILD(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &         IFCKOP,NPOS,INTFLG,IPRINT,WORK,LWORK)
        HFXFAC = HFXTMP
      END IF
C
C        Manu: if no need for long-range coulomb integrals (because they 
C              can be recombined with the short-range ones or we not need them) 
C              then proceed as usual (compute full range integrals). If needed 
C              (for the MCSCF-srDFT for example), they have already been calculated. 
C              The call to FOCKBUILD is then skipped.  
C
#ifdef MOD_SRDFT
      if (.not.srdft_cfg_lrcoulomb_int) then
#endif
         HFXTMP = HFXMU
         HFXMU = D0
         CALL FOCKBUILD(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &                  IFCKOP,NPOS,INTFLG,IPRINT,WORK,LWORK)
         HFXMU = HFXTMP
#ifdef MOD_SRDFT
      end if 
#endif
C
C     ===========================
      IF (ONECAP.AND.SMLV1C) THEN
C     ===========================
C
C        Correct for neglected multi-center contributions via a 1-center
C        expansion that approximates the missing integrals. Restore 
C        density matrix after using it as scratch space.
C        TODO : Reduce memory usage 
C
         CALL MEMGET2 ('REAL','BUF1',KBUF1,4*N2BBASX,WORK,KFREE,LFREE)
         CALL OPNFIL(LUPMAT,'DFPMAT','OLD','TWOFCK')
C
C        CAUTION: Here we read the projector in the scratch array DMAT !
C
         CALL READT(LUPMAT,N2BBASX,DMAT)
         DO I = 1,NFMAT
            IREP = ISYMOP(I)-1
            CALL PR1CEX2 (NTBAS(0),IREP,FMAT,DMAT,WORK(KBUF1),IPRINT)
         ENDDO
C
C        Now we do read the original (full) DMAT !
C
         CALL READT(LUPMAT,N2BBASXQ*NFMAT,DMAT)
         CLOSE(LUPMAT,STATUS='KEEP')
         CALL MEMREL('TWOFCK',WORK,KBUF1,KBUF1,KFREE,LFREE)
C     =====
      ENDIF
C     =====
C
C     Insert half-phases and symmetrize
C     =================================
C
      DO I = 1,NFMAT
C
C       Reindex to sorted basis
C       =======================
C
        IF(IPRINT.GE.5) THEN
          CALL HEADER('Unsorted FMAT before half-phases inserted',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &               NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
        CALL BUTOBS(FMAT(1,1,I),NZ,WORK,LWORK)
        CALL BUTOBS(DMAT(1,1,I),NZ,WORK,LWORK)
        IREP = ISYMOP(I)-1
        IF(NZ.LT.4) THEN
          DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IREP)
            CALL Q2BPHASE('F',IQ,1,FMAT(1,IZ,I))
            CALL Q2BPHASE('D',IQ,1,DMAT(1,IZ,I))
          ENDDO
          IF(IPRINT.GE.5) THEN
            CALL HEADER('Unsorted FMAT with half-phases inserted',-1)
            WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
            CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZ,IPQTOQ(1,IREP),LUPRI)
          ENDIF
        ENDIF
C
C       Symmetrize FMAT matrix
C       ======================
C       IH = 0   Non-symmetric
C       IH = 1   Symmetric
C       IH = 2   Anti-symmetric
C
        IF(IPRINT.GE.5) THEN
          CALL HEADER('unSorted FMAT matrix',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
        DO IZ = 1,NZ
          IQ = IPQTOQ(IZ,IREP)
          IH = IHQMAT(IQ,IHRMOP(I)) 
          IF    (IH.EQ.1) THEN
            CALL FULMAT('S',NTBAS(0),NTBAS(0),FMAT(1,IZ,I))
          ELSEIF(IH.EQ.2) THEN
            CALL FULMAT('A',NTBAS(0),NTBAS(0),FMAT(1,IZ,I))
          ENDIF
        ENDDO
        IF(IPRINT.GE.4) THEN
          CALL HEADER('Sorted FMAT matrix',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(FMAT(1,1,I),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
        ENDIF
      ENDDO
      CALL QEXIT('TWOFCK')
      iprint = iprint_save
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Twoprp */
      SUBROUTINE TWOPRP(WORK,LWORK)
C*****************************************************************************
C
C     Prepare for two-electron integrals
C
C     Written by T.Saue Sep 30 1996
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "cbihr2.h"
#include "dcbind.h"
#include "dcbbas.h"
#include "hrunit.h"
      DIMENSION WORK(LWORK)
      LUINTA = 11
        LUSUPM = 17
        LUDASP = 24
        LUSOL  = 49
        CALL PAOVEC(WORK,LWORK,0,IPRTWO)
C
C     Define bit packing
C     ==================
C
#if defined (SYS_CRAY) || defined (SYS_T3D) || defined (INT_STAR8)
      NIBUF = 1
      IF (NTBAS(0) .LE. 255) THEN
         NBITS = 8
      ELSE
         NBITS = 16
      END IF
#else
      IF (NTBAS(0) .LE. 255) THEN
         NIBUF = 1
         NBITS = 8
      ELSE
         NIBUF = 2
         NBITS = 16
      END IF
#endif
      NBIT1   = 2**NBITS
      NBIT2   = 2**(2*NBITS)
      IBIT1   = NBIT1 - 1
      IBIT2   = NBIT2 - 1
C
C     Calculate integrals to be stored on disk
C     ========================================
C
      CALL TWOGEN(WORK,LWORK)
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twosrt */
      SUBROUTINE TWOSRT(ITYP,NCHUNK,IREPA,WORK,LWORK,IPRINT)
C*****************************************************************************
C
C	This files provides sorting of XX-integrals based on index relations
C		Class D,Ca and Cb sorted in buffers of size N2BBASX
C		Class Ba and Bb placed in matrices of dimension NDIM2
C		Class A placed in array of length NDIM
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
#include "dcbbas.h"
#include "dcbind.h"
#include "dgroup.h"
      CHARACTER FILNAM*6
      DIMENSION	WORK(LWORK),IREPA(N2BBASX),LUINT(0:8)
      CALL QENTER('TWOSRT')
#include "memint.h"
      LU2HER    = 13
      DO I = 0,8
        LUINT(I) = LU2HER + 1 + I
      ENDDO
      GO TO (11,12,13,14), ITYP
   11 CONTINUE
C
C     *****************************
C     ***** LLLL - integrals ******
C     *****************************
C
C
C       **************************
C	***** Coarse sorting *****
C       **************************
C
C       File assignment:
C         LUINTA - file of sorted LL-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted LL-integrals (HERMIT-format, deleted in DO1SRT)
C         LU2BUF - buffer file of UNSORTED LL-integrals, other classes
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from LL1SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER  ,'DFTWLL','OLD'    ,'LL1SRT')
        CALL OPNFIL(LUINT(0),'DFLLSA','UNKNOWN','LL1SRT')
        DO I = 1,NBSYM
          FILNAM = 'DFLLBF'
          WRITE(FILNAM(6:6),'(I1)') I-1
          CALL OPNFIL(LUINT(I),FILNAM,'UNKNOWN','LL1SRT')          
        ENDDO
C
        MBUF = NBSYM*2
        CALL MEMGET2('REAL','TWOMT.LLc',KTWOMT,N2BBASX*MBUF      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','INDMT.LLc',KINDMT,N2BBASX*MBUF*NIBUF,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BUF.LLc'  ,KBUF  ,NCHUNK            ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBUF.LLc' ,KIBUF ,NCHUNK*NIBUF      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','CHUNK.LLc',KIND  ,NCHUNK*2          ,
     &       WORK,KFREE,LFREE)
C
        CALL XX1SRT(LUINT,LU2HER,WORK(KTWOMT),WORK(KINDMT),N2BBASX,
     &              WORK(KBUF),WORK(KIBUF),WORK(KIND),IREPA,
     &              NCHUNK,IPRINT)
C
        CALL MEMREL('LL1SRT',WORK,KWORK,KWORK,KFREE,LFREE)
        CLOSE(LU2HER,STATUS='DELETE')
C
C       *************************
C       ***** Fine sorting ******
C       *************************
C
C
        NNDIM = NTBAS(1)*(NTBAS(1)-1)/2
        CALL MEMGET2('REAL','TWOMT.LLf',KTWOMT,N2BBASX        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','INDMT.LLf',KINDMT,N2BBASX*NIBUF  ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IND.LLf'  ,KIND  ,N2BBASX*4      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','CMAT.LLf' ,KCMAT ,2*N2BBASX      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','ICMAT.LLf',KICMAT,2*NIBUF*N2BBASX,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BMAT.LLf' ,KBMAT ,4*NNDIM        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBMAT.LLf',KIBMAT,4*NNDIM        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','AVEC.LLf' ,KAVEC ,NTBAS(1)       ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IAVEC.LLf',KIAVEC,NTBAS(1)*NIBUF ,
     &       WORK,KFREE,LFREE)
C
        CALL OPNFIL(LU2HER,'DFLLSB','UNKNOWN','LL2SRT')
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from LL2SRT',-1)
        ENDIF
        CALL XX2SRT(LUINT(0),LU2HER,LUINT(1),N2BBASX,
     &           NTBAS(1),NNDIM,WORK(KTWOMT),WORK(KINDMT),WORK(KIND),
     &           WORK(KCMAT),WORK(KICMAT),WORK(KBMAT),WORK(KIBMAT),
     &           WORK(KAVEC),WORK(KIAVEC),IPRINT)
        CLOSE(LUINT(0),STATUS='KEEP')
        CLOSE(LU2HER  ,STATUS='KEEP')
        DO I = 1,NBSYM
          CLOSE(LUINT(I),STATUS = 'DELETE')
        ENDDO
C
        CALL MEMREL('LL2SRT',WORK,KWORK,KWORK,KFREE,LFREE)
      GOTO 10
   12 CONTINUE
C
C     *****************************
C     ***** LLSS - integrals ******
C     *****************************
C
C
C     **************************
C	***** Coarse sorting *****
C     **************************
C
C       File assignment:
C         LUINTA - file of sorted LL-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted LL-integrals (HERMIT-format, deleted in DO1SRT)
C         LU2BUF - buffer file of UNSORTED LL-integrals, other classes
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SL1SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER  ,'DFTWSL','OLD'    ,'SL1SRT')
        CALL OPNFIL(LUINT(0),'DFSLSA','UNKNOWN','SL1SRT')
        DO I = 1,NBSYM
          FILNAM = 'DFSLBF'
          WRITE(FILNAM(6:6),'(I1)') I-1
          CALL OPNFIL(LUINT(I),FILNAM,'UNKNOWN','SL1SRT')          
        ENDDO
C
        MBUF = NBSYM*2
        CALL MEMGET2('REAL','TWOMT.SLc',KTWOMT,N2BBASX*MBUF      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','INDMT.SLc',KINDMT,N2BBASX*MBUF*NIBUF,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BUF.SLc'  ,KBUF  ,NCHUNK            ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBUF.SLc' ,KIBUF ,NCHUNK*NIBUF      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IND.SLc'  ,KIND  ,NCHUNK*2          ,
     &       WORK,KFREE,LFREE)
C
        CALL XY1SRT(LUINT,LU2HER,WORK(KTWOMT),WORK(KINDMT),N2BBASX,
     &              WORK(KBUF),WORK(KIBUF),WORK(KIND),IREPA,
     &              NCHUNK,IPRINT)
C
        CALL MEMREL('SL1SRT',WORK,KWORK,KWORK,KFREE,LFREE)
        CLOSE(LU2HER,STATUS='DELETE')
C
C       *************************
C       ***** Fine sorting ******
C       *************************
C
C       File assignment:
C         LUINTA - file of sorted LL-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted LL-integrals (HERMIT-format)
C         LU2BUF - buffer file of SORTED integrals, other classes
C
C	Memory assignment
C
        NSLDIM = NTBAS(1)*NTBAS(2)
        CALL MEMGET2('REAL','BUF.SLf'  ,KBUF  ,N2BBASX        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBUF.SLf' ,KIBUF ,N2BBASX*NIBUF  ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IND.SLf'  ,KIND  ,N2BBASX*4      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','AMAT.SLf' ,KAMAT ,NSLDIM         ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IAMAT.SLf',KIAMAT,NSLDIM         ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BCMAT.SLf',KBCMAT,2*N2BBASX      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBCMT.SLf',KIBCMT,2*N2BBASX*NIBUF,
     &       WORK,KFREE,LFREE)
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SL2SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER,'DFSLSB','UNKNOWN','SL2SRT')
        CALL XY2SRT(LUINT(0),LU2HER,LUINT(1),N2BBASX,
     &              NSLDIM,WORK(KBUF),WORK(KIBUF),WORK(KIND),
     &              WORK(KAMAT),WORK(KIAMAT),
     &              WORK(KBCMAT),WORK(KIBCMT),IPRINT)
        CLOSE(LUINT(0),STATUS='KEEP')
        CLOSE(LU2HER  ,STATUS='KEEP')
        DO I = 1,NBSYM
          CLOSE(LUINT(I),STATUS = 'DELETE')
        ENDDO
        CALL MEMREL('SL2SRT',WORK,KWORK,KWORK,KFREE,LFREE)
      GOTO 10
   13 CONTINUE
C
C     *****************************
C     ***** SSSS - integrals ******
C     *****************************
C
C       **************************
C	***** Coarse sorting *****
C       **************************
C
C       File assignment:
C         LUINTA - file of sorted SS-integrals , all indices unequal (class D)
C         LUINTB - file of unsorted SS-integrals (HERMIT-format, deleted in DO1SRT)
C         LU2BUF - buffer file of UNSORTED SS-integrals, other classes
C
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SS1SRT',-1)
        ENDIF
        CALL OPNFIL(LU2HER  ,'DFTWSS','OLD'    ,'SS1SRT')
        CALL OPNFIL(LUINT(0),'DFSSSA','UNKNOWN','SS1SRT')
        DO I = 1,NBSYM
          FILNAM = 'DFSSBF'
          WRITE(FILNAM(6:6),'(I1)') I-1
          CALL OPNFIL(LUINT(I),FILNAM,'UNKNOWN','SS1SRT')          
        ENDDO
C
        MBUF = NBSYM*2
        CALL MEMGET2('REAL','TWOMT.SSc',KTWOMT,N2BBASX*MBUF      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','INDMT.SSc',KINDMT,N2BBASX*MBUF*NIBUF,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BUF.SSc'  ,KBUF  ,NCHUNK            ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBUF.SSc' ,KIBUF ,NCHUNK*NIBUF      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IND.SSc'  ,KIND  ,NCHUNK*2          ,
     &       WORK,KFREE,LFREE)
        CALL XX1SRT(LUINT,LU2HER,WORK(KTWOMT),WORK(KINDMT),N2BBASX,
     &              WORK(KBUF),WORK(KIBUF),WORK(KIND),IREPA,
     &              NCHUNK,IPRINT)
        CALL MEMREL('SS1SRT',WORK,KWORK,KBUF,KFREE,LFREE)
        CLOSE(LU2HER,STATUS='DELETE')
C
C       *************************
C       ***** Fine sorting ******
C       *************************
C
C
        NNDIM = NTBAS(2)*(NTBAS(2)-1)/2
        CALL MEMGET2('REAL','TWOMT.SSf',KTWOMT,N2BBASX        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','INDMT.SSf',KINDMT,N2BBASX*NIBUF  ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IND.SSf'  ,KIND  ,N2BBASX*4      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','CMAT.SSf' ,KCMAT ,2*N2BBASX      ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','ICMAT.SSf',KICMAT,2*NIBUF*N2BBASX,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','BMAT.SSf' ,KBMAT ,4*NNDIM        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IBMAT.SSf',KIBMAT,4*NNDIM        ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','AVEC.SSf' ,KAVEC ,NTBAS(2)       ,
     &       WORK,KFREE,LFREE)
        CALL MEMGET2('INTE','IAVEC.SSf',KIAVEC,NTBAS(2)*NIBUF ,
     &       WORK,KFREE,LFREE)
C
        CALL OPNFIL(LU2HER,'DFSSSB','UNKNOWN','SS2SRT')
        IF(IPRINT.GE.3) THEN
          CALL HEADER('Output from SS2SRT',-1)
        ENDIF
        CALL XX2SRT(LUINT(0),LU2HER,LUINT(1),N2BBASX,
     &           NTBAS(2),NNDIM,WORK(KTWOMT),WORK(KINDMT),WORK(KIND),
     &           WORK(KCMAT),WORK(KICMAT),WORK(KBMAT),WORK(KIBMAT),
     &           WORK(KAVEC),WORK(KIAVEC),IPRINT)
        CLOSE(LUINT(0),STATUS='KEEP')
        CLOSE(LU2HER  ,STATUS='KEEP')
        DO I = 1,NBSYM
          CLOSE(LUINT(I),STATUS = 'DELETE')
        ENDDO
C
        CALL MEMREL('SS2SRT',WORK,KWORK,KWORK,KFREE,LFREE)
      GOTO 10
   14 CONTINUE
C
C     Gaunt integrals
C
      GOTO 10
   10 CONTINUE
C
      CALL QEXIT('TWOSRT')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xx1srt */
      SUBROUTINE XX1SRT(LTFIL,LUTWXX,TWOMAT,INDMAT,MAXBUF,
     &                  BUF,IBUF,IND,IREPA,NCHUNK,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Coarse sorting of 2-electron XX - integrals
C
C     METHOD:
C
C     FILE ASSIGNMENT:
C      LUTWXX   - 'AOTWOXX'    - file of unsorted XX-integrals
C      LTFIL(0) - 'DF2XXA'     - file of sorted XX-integrals (D,Cb)
C      LTFIL(1) - 'DF2XXC'     - buffer file of XX-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size N2BBASX) to be read
C       1 - XX, class D
C       2 - XX-buffer
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C
C************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION LTFIL(0:8),ICOUNT(0:7,0:1),NINTS(0:7,0:1)
      DIMENSION IBUF(NCHUNK,NIBUF),BUF(NCHUNK),IND(2,NCHUNK)
      DIMENSION TWOMAT(MAXBUF,0:MAXREP,0:1),
     &          INDMAT(MAXBUF,NIBUF,0:MAXREP,0:1),
     &          IREPA(*)
#include "ibtfun.h"
C
      CALL IZERO(NINTS  ,16)
      CALL IZERO(ICOUNT,16)
      REWIND LUTWXX
      CALL MOLLAB('BASTWOEL',LUTWXX,LUPRI)
C*** DO WHILE - loop: read a chunk of integrals, perform  sorting
      MINT = 0
   10 CONTINUE
      READ(LUTWXX,IOSTAT = IOS) BUF,IBUF,NUT
      IF(NUT.LE.0) GOTO 30
      MINT = MINT + NUT
      IF (NIBUF .EQ. 1) THEN
C
C       Unpack indices and extract information about
C         1. Combined symmetry of IJ
C         2. Whether any indices are equal
C
        DO INT = 1,NUT
          I     = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
          J     = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
          K     = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
          L     = IBTAND(       IBUF(INT,1),         IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - K)*(I - J)*(J - K)*(J - L)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      ELSE
C
C       Unpack indices and extract information about
C         1. Whether any indices are equal
C         2. Combined symmetry of IJ
C
        DO INT = 1,NUT
          I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
          J = IBTAND(       IBUF(INT,1),       IBIT1)
          K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
          L = IBTAND(       IBUF(INT,2),       IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - K)*(I - J)*(J - K)*(J - L)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
          INDMAT(INDEX,2,IREP,IT) = IBUF(INT,2)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      END IF
      GOTO 10
C
C*** END DO WHILE - loop: reading of integrals from LUTWXX
C
   30 CONTINUE
C*** Empty matrices:
      DO IT = 0,1
        DO IREP = 0,MAXREP
          IF(NINTS(IREP,IT).GT.0) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),NINTS(IREP,IT),IREP,MAXBUF,NIBUF,
     &               INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            NINTS(IREP,IT) = NINTS(IREP,IT) + ICOUNT(IREP,IT)*MAXBUF
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
          ELSE
            NINTS(IREP,IT) = ICOUNT(1,1)*MAXBUF
          ENDIF
        ENDDO
      ENDDO       
      DO IUNIT = 0,MAXREP+1
        WRITE(LTFIL(IUNIT)) -1,0
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - coarse sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class D integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Other classes    :'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,0),IREP=0,MAXREP)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xy1srt */
      SUBROUTINE XY1SRT(LTFIL,LUTWXX,TWOMAT,INDMAT,MAXBUF,
     &                  BUF,IBUF,IND,IREPA,NCHUNK,IPRINT)
C*****************************************************************************
C
C     PURPOSE: Coarse sorting of 2-electron XY - integrals
C
C     METHOD:
C
C     FILE ASSIGNMENT:
C      LUTWXX   - 'AOTWOXX'    - file of unsorted XX-integrals
C      LTFIL(0) - 'DF2XXA'     - file of sorted XX-integrals (D,Cb)
C      LTFIL(1) - 'DF2XXC'     - buffer file of XX-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size N2BBASX) to be read
C       1 - XX, class D
C       2 - XX-buffer
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C
C************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION LTFIL(0:8),ICOUNT(0:7,0:1),NINTS(0:7,0:1)
      DIMENSION IBUF(NCHUNK,NIBUF),BUF(NCHUNK),IND(2,NCHUNK)
      DIMENSION TWOMAT(MAXBUF,0:MAXREP,0:1),
     &          INDMAT(MAXBUF,NIBUF,0:MAXREP,0:1),
     &          IREPA(*)
#include "ibtfun.h"
C
      CALL IZERO(NINTS  ,16)
      CALL IZERO(ICOUNT,16)
      REWIND LUTWXX
      CALL MOLLAB('BASTWOEL',LUTWXX,LUPRI)
C*** DO WHILE - loop: read a chunk of integrals, perform  sorting
      MINT = 0
   10 CONTINUE
      READ(LUTWXX,IOSTAT = IOS) BUF,IBUF,NUT
      IF(NUT.LE.0) GOTO 30
      MINT = MINT + NUT
      IF (NIBUF .EQ. 1) THEN
C
C       Unpack indices and extract information about
C         1. Combined symmetry of IJ
C         2. Whether any indices are equal
C
        DO INT = 1,NUT
          I     = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
          J     = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
          K     = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
          L     = IBTAND(       IBUF(INT,1),         IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - J)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      ELSE
C
C       Unpack indices and extract information about
C         1. Whether any indices are equal
C         2. Combined symmetry of IJ
C
        DO INT = 1,NUT
          I = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
          J = IBTAND(       IBUF(INT,1),       IBIT1)
          K = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
          L = IBTAND(       IBUF(INT,2),       IBIT1)
          IND(1,INT) = IBTXOR(IREPA(I),IREPA(J))
          IT         = IABS((I - J)*(K - L))
          IND(2,INT) = MIN(1,IT)
        ENDDO
C 
C       Place integrals in buffer
C
        DO INT = 1,NUT
          IREP  = IND(1,INT)
          IT    = IND(2,INT)
          NINTS(IREP,IT) = NINTS(IREP,IT) + 1
          INDEX         = NINTS(IREP,IT)
          TWOMAT(INDEX,IREP,IT)   = BUF(INT)
          INDMAT(INDEX,1,IREP,IT) = IBUF(INT,1)
          INDMAT(INDEX,2,IREP,IT) = IBUF(INT,2)
C*** Test to see whether matrix is full:
          IF (INDEX.EQ.MAXBUF) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),INDEX,IREP,MAXBUF,NIBUF,
     &                 INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
            NINTS  (IREP,IT) = 0
          ENDIF
        ENDDO
      END IF
      GOTO 10
C
C*** END DO WHILE - loop: reading of integrals from LUTWXX
C
   30 CONTINUE
C*** Empty matrices:
      DO IT = 0,1
        DO IREP = 0,MAXREP
          IF(NINTS(IREP,IT).GT.0) THEN
            IUNIT = (1-IT)*(1+IREP)
            CALL WRINT(LTFIL(IUNIT),NINTS(IREP,IT),IREP,MAXBUF,NIBUF,
     &               INDMAT(1,1,IREP,IT),TWOMAT(1,IREP,IT))
            NINTS(IREP,IT) = NINTS(IREP,IT) + ICOUNT(IREP,IT)*MAXBUF
            ICOUNT(IREP,IT) = ICOUNT(IREP,IT) + 1
          ELSE
            NINTS(IREP,IT) = ICOUNT(1,1)*MAXBUF
          ENDIF
        ENDDO
      ENDDO       
      DO IUNIT = 0,MAXREP+1
        WRITE(LTFIL(IUNIT)) -1,0
      ENDDO
C
C     Print section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - coarse sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class D integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Other classes    :'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NINTS(IREP,0),IREP=0,MAXREP)
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xx2srt */
      SUBROUTINE XX2SRT(LUINTA,LUINTB,LUBUF,MAXBUF,
     &                  NDIM,NNDIM,BUF,IBUF,IND,CMAT,ICMAT,
     &                  BMAT,IBMAT,AVEC,IAVEC,IPRINT)
C*****************************************************************************
C
C     PURPOSE:
C        Fine sorting of 2-electron XX-integrals:
C
C     METHOD:
C        Twoelectron XX-integrals IJKL can be classified as
C        follows according to indices:
C
C        *Group A: All indices equal
C            Class  1: I = J = K = L
C        *Group B: Three indices equal
C            Class a1: I > J = K = L
C            Class a2: I = J = K > L
C            Class b1: I = J > K = L
C            Class b2: I = K > J = L
C        *Group C: two indices equal
C            Class a1: I = K > J > L
C            Class a2: I > J = K > L
C            Class a3: I > K > J = L
C            Class b1: I > K = L > J
C            Class b2: I > J > K = L
C            Class b3: I = J > K > L
C        *Group D : all indices different(sorted in XX1SRT):
C            Class  1: I > K > J > L
C            Class  2: I > K > L > J
C            Class  3: I > J > K > L
C
C        In this routine the integrals are sorted, that is
C        distributed on different files, according to how they
C        contribute to the twoelectron Fock matrix
C
C     FILE ASSIGNMENT:
C      LUINTA - 'DF2XXA'     - file of sorted XX-integrals (D,Cb)
C      LUINTB - 'DF2XXB'     - file of sorted XX-integrals (Ca,Bb,Ba,A)
C      LU2BUF - 'DF2XXC'     - buffer file of XX-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size N2BBASX) to be read
C       1 - Input: XX-buffer.  Output: XX, class Ca
C       2 - XX, class Cb
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C     Written by T.Saue Feb 27 1995
C     Last revision: Feb 27 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION BUF(MAXBUF),IBUF(MAXBUF,NIBUF),IND(4,MAXBUF),
     &          CMAT(MAXBUF,2),ICMAT(MAXBUF,NIBUF,2),
     &          BMAT(NNDIM,4)  ,IBMAT(NNDIM,4),
     &          AVEC(NDIM),IAVEC(NDIM),LUBUF(0:7)
      DIMENSION NC(0:7,2),ICOUNT(0:7,2),NB(0:7,4),NBT(4)
#include "ibtfun.h"
C
C     Initialize
C
      REWIND LUINTB
      CALL IZERO(ICOUNT,16)
      CALL IZERO(NC,16)
      CALL IZERO(NB,32)
      CALL IZERO(NBT,4)
      NA = 0
C
C     Loop over irreps
C      
      DO IREP = 0,MAXREP
        LUINTC = LUBUF(IREP)
        REWIND LUINTC
C
C       Loop over batches in buffer
C
   10   CONTINUE
        CALL REAINT(LUINTC,NINTS,IRP,MAXBUF,NIBUF,IBUF,BUF)
        IF(NINTS.EQ.-1) GOTO 30
C
C       Unpack indices
C
        IF (NIBUF .EQ. 1) THEN
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
            IND(2,INT) = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,1),         IBIT1)
          ENDDO
        ELSE
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
            IND(2,INT) = IBTAND(       IBUF(INT,1),       IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,2),       IBIT1)
          ENDDO
        END IF
C
C       Begin sorting
C
        DO 50 INT = 1,NINTS
          I = IND(1,INT)
          J = IND(2,INT)
          K = IND(3,INT)
          L = IND(4,INT)
          IJ = I - J
          JL = J - L
          IF((IJ*JL).NE.0) THEN
C
C***  Group C: Class a1,a2,b1,b2
C
            KL = K - L
            IF(KL.EQ.0) THEN
C
C*** Group C:  Class b1,b2 need not be distinguished
C
              IND1 = NBIT1*I + J
              IND2 = K
              NC(IREP,2) = NC(IREP,2) + 1
              CMAT(NC(IREP,2),2) = BUF(INT)
              IF(NIBUF.EQ.1) THEN
                ICMAT(NC(IREP,2),1,2) = IND1*NBIT2+IND2
              ELSE
                ICMAT(NC(IREP,2),1,2) = IND1
                ICMAT(NC(IREP,2),2,2) = IND2
              ENDIF
              IF(NC(IREP,2).EQ.MAXBUF) THEN
                CALL WRINT(LUINTB,NC(IREP,2),IREP,MAXBUF,NIBUF,
     &                     ICMAT(1,1,2),CMAT(1,2))
                ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
                NC(IREP,2) = 0
              ENDIF
            ELSE
C
C***  Group C: Class a1,a2
C
              IF(I.EQ.K) THEN
C**** Group C: Class a1:
                IND1 = NBIT1*J + L
                IND2 = I
              ELSE
C**** Group C: Class a2:
                IND1 = NBIT1*I + L
                IND2 = J
              ENDIF
              NC(IREP,1) = NC(IREP,1) + 1
              CMAT(NC(IREP,1),1) = BUF(INT)
              IF(NIBUF.EQ.1) THEN
                ICMAT(NC(IREP,1),1,1) = IND1*NBIT2+IND2
              ELSE
                ICMAT(NC(IREP,1),1,1) = IND1
                ICMAT(NC(IREP,1),2,1) = IND2
              ENDIF
              IF(NC(IREP,1).EQ.MAXBUF) THEN
                CALL WRINT(LUINTA,NC(IREP,1),IREP,MAXBUF,NIBUF,
     &                     ICMAT(1,1,1),CMAT(1,1))
                ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
                NC(IREP,1) = 0
              ENDIF
            ENDIF
          ELSE
            IK = I - K
            KL = K - L
            IF((IK*KL).GT.0) THEN
C
C*** Group C : Class a3,b3
C
              IF(IJ.EQ.0) THEN
C
C*** Group C: Class b3
C
                NC(IREP,2)         = NC(IREP,2) + 1
                CMAT(NC(IREP,2),2) = BUF(INT)
                IND1               = NBIT1*K + L
                IND2               = I
                IF(NIBUF.EQ.1) THEN
                  ICMAT(NC(IREP,2),1,2) = IND1*NBIT2+IND2
                ELSE
                  ICMAT(NC(IREP,2),1,2) = IND1
                  ICMAT(NC(IREP,2),2,2) = IND2
                ENDIF
                IF(NC(IREP,2).EQ.MAXBUF) THEN
                  CALL WRINT(LUINTB,NC(IREP,2),IREP,MAXBUF,NIBUF,
     &                       ICMAT(1,1,2),CMAT(1,2))
                  ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
                  NC(IREP,2) = 0
                ENDIF
              ELSE
C
C*** Group C: Class a3
C
                NC(IREP,1)         = NC(IREP,1) + 1
                CMAT(NC(IREP,1),1) = BUF(INT)
                IND1          = NBIT1*I + K
                IND2          = J
                IF(NIBUF.EQ.1) THEN
                  ICMAT(NC(IREP,1),1,1) = IND1*NBIT2+IND2
                ELSE
                  ICMAT(NC(IREP,1),1,1) = IND1
                  ICMAT(NC(IREP,1),2,1) = IND2
                ENDIF
                IF(NC(IREP,1).EQ.MAXBUF) THEN
                    CALL WRINT(LUINTA,NC(IREP,1),IREP,MAXBUF,NIBUF,
     &                       ICMAT(1,1,1),CMAT(1,1))
                  ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
                  NC(IREP,1) = 0
                ENDIF
              ENDIF
            ELSEIF(K.GT.J) THEN
C
C*** Group B - b2
C
              NB(IREP,1)      = NB(IREP,1) + 1
              NBT(1)          = NBT(1) + 1
              BMAT(NBT(1),1)  = BUF(INT)
              IBMAT(NBT(1),1) = NBIT1*J+K
            ELSEIF(K.LT.J) THEN
C
C*** Group B - b1
C
              NB(IREP,4)      = NB(IREP,4) + 1
              NBT(4)          = NBT(4) + 1
              BMAT(NBT(4), 4) = BUF(INT)
              IBMAT(NBT(4),4) = NBIT1*J+K
            ELSEIF(I.GT.J) THEN
C
C*** Group B - a1
C
              NB(IREP,2)      = NB(IREP,2)+1
              NBT(2)          = NBT(2)+1
              BMAT(NBT(2), 2) = BUF(INT)
              IBMAT(NBT(2),2) = NBIT1*I+J
            ELSEIF(I.GT.L) THEN
C
C*** Group B - a2
C
              NB(IREP,3)      = NB(IREP,3) + 1
              NBT(3)          = NBT(3) + 1
              BMAT(NBT(3), 3) = BUF(INT)
              IBMAT(NBT(3),3) = NBIT1*L + I
            ELSE
C
C*** Group A:
C
              NA        = NA + 1
              AVEC(NA)  = BUF(INT)
              IAVEC(NA) = I
            ENDIF
          ENDIF
 50       CONTINUE
        GOTO 10
   30   CONTINUE
C*** Empty matrices:
        IF(NC(IREP,1).GT.0) THEN
          CALL WRINT(LUINTA,NC(IREP,1),IREP,MAXBUF,NIBUF,
     &               ICMAT(1,1,1),CMAT(1,1))
          NC(IREP,1)  = NC(IREP,1) + ICOUNT(IREP,1)*MAXBUF
          ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
        ELSE
          NC(IREP,1) = ICOUNT(IREP,1)*MAXBUF
        ENDIF
        IF(NC(IREP,2).GT.0) THEN
          CALL WRINT(LUINTB,NC(IREP,2),IREP,MAXBUF,NIBUF,
     &               ICMAT(1,1,2),CMAT(1,2))
          NC(IREP,2) = NC(IREP,2) + ICOUNT(IREP,2)*MAXBUF
          ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
        ELSE
          NC(IREP,2) = ICOUNT(IREP,2)*MAXBUF
        ENDIF
      ENDDO
C
      WRITE(LUINTA) -1,0
      WRITE(LUINTB) -1,0
      IOFF = 1
      DO IREP = 0,MAXREP
        CALL WRINT(LUINTB,NB(IREP,1),IREP,NNDIM,1,
     &                      IBMAT(IOFF,1),BMAT(IOFF,1))
        IOFF = IOFF + NB(IREP,1)
      ENDDO
C     Class Ba1, Ba2 and Bb1 has totally symmetric densities !
      DO I = 2,4
        CALL WRINT(LUINTB,NB(0,I),0,NNDIM,1,
     &                      IBMAT(1,I),BMAT(1,I))
      ENDDO
      CALL WRINT(LUINTB,NA,0,NDIM,1,IAVEC,AVEC)
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - fine sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class Ca integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NC(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Cb integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NC(IREP,2),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Ba1 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,2),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Ba2 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,3),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Bb1 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,4),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class Bb2 integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP,1),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class A   integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       'Irrep ',REP(0),NA
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrint */
      SUBROUTINE WRINT(IUNIT,NINTS,IREP,NDIM,NIBUF,IBUF,BUF)
C*****************************************************************************
C
C     Write integral and indices to unformatted file
C
C     Written by T.Saue, Feb 27 1995
C     Last revision: Feb 27 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BUF(NDIM),IBUF(NDIM,NIBUF)
C
      WRITE(IUNIT) NINTS,IREP
      IF(NINTS.LE.0) RETURN
      CALL WRITI(IUNIT,NINTS,IBUF(1,1))
      IF(NIBUF.EQ.2) CALL WRITI(IUNIT,NINTS,IBUF(1,2))
      CALL WRITT(IUNIT,NINTS,BUF)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reaint */
      SUBROUTINE REAINT(IUNIT,NINTS,IREP,NDIM,NIBUF,IBUF,BUF)
C*****************************************************************************
C
C     Read integral and indices from unformatted file
C
C     Written by T.Saue, Feb 27 1995
C     Last revision: Feb 27 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION BUF(NDIM),IBUF(NDIM,NIBUF)
C
      READ(IUNIT) NINTS,IREP
      IF(NINTS.LE.0) RETURN
      CALL READI(IUNIT,NINTS,IBUF(1,1))
      IF(NIBUF.EQ.2) CALL READI(IUNIT,NINTS,IBUF(1,2))
      CALL READT(IUNIT,NINTS,BUF)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xy2srt */
      SUBROUTINE XY2SRT(LUINTA,LUINTB,LUBUF,MAXBUF,
     &                  NSLDIM,BUF,IBUF,IND,AMAT,IAMAT,
     &                  BCMAT,IBCMAT,IPRINT)
C*****************************************************************************
C
C     PURPOSE:
C        Fine sorting of 2-electron XY-integrals:
C
C     METHOD:
C        Two-electron XY-integrals IJKL can be classified as
C        follows according to indices:
C
C        *Group A: I = J,  K = L
C        *Group B: I > J,  K = L
C        *Group C: I = J,  K > L
C        *Group D: I > J,  K > L
C
C        In this routine the integrals are sorted, that is
C        distributed on different files, according to how they
C        contribute to the twoelectron Fock matrix
C
C     FILE ASSIGNMENT:
C      LUINTA - 'DF2XYA'     - file of sorted XY-integrals (D,B)
C      LUINTB - 'DF2XYB'     - file of sorted XY-integrals (C,A)
C      LU2BUF - 'DF2XYC'     - buffer file of XY-integrals (C,B,A)
C
C     ICOUNT: number of buffers (size MAXBUF) to be read
C       1 - Input: XY-buffer.  Output: XY, class Ca
C       2 - XY, class Cb
C
C     LIBRARY ROUTINES:
C         This routine employs machine-dependent bitoperations
C         which are defined in statement function statements
C
C     Written by T.Saue Feb 27 1995
C     Last revision : Feb 27 1995
C
C************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "pgroup.h"
#include "dcbind.h"
      DIMENSION BUF(MAXBUF),IBUF(MAXBUF,NIBUF),IND(4,MAXBUF),
     &          AMAT(NSLDIM),IAMAT(NSLDIM),
     &          BCMAT(MAXBUF,2),IBCMAT(MAXBUF,NIBUF,2),
     &          ICOUNT(0:7,2),LUBUF(0:7),NB(0:7),NC(0:7)
#include "ibtfun.h"
      REWIND LUINTB
      CALL IZERO(ICOUNT,16)
      CALL IZERO(NC,8)
      CALL IZERO(NB,8)
      NA = 0
C
C     Loop over irreps
C      
      DO IREP = 0,MAXREP
        LUINTC = LUBUF(IREP)
        REWIND LUINTC
C
C       Loop over batches in buffer
C
   10   CONTINUE
        CALL REAINT(LUINTC,NINTS,IRP,MAXBUF,NIBUF,IBUF,BUF)
        IF(NINTS.EQ.-1) GOTO 30
C
C       Unpack indices
C
        IF (NIBUF .EQ. 1) THEN
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),3*NBITS),IBIT1)
            IND(2,INT) = IBTAND(IBTSHR(IBUF(INT,1),2*NBITS),IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,1),  NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,1),         IBIT1)
          ENDDO
        ELSE
          DO INT = 1,NINTS
            IND(1,INT) = IBTAND(IBTSHR(IBUF(INT,1),NBITS),IBIT1)
            IND(2,INT) = IBTAND(       IBUF(INT,1),       IBIT1)
            IND(3,INT) = IBTAND(IBTSHR(IBUF(INT,2),NBITS),IBIT1)
            IND(4,INT) = IBTAND(       IBUF(INT,2),       IBIT1)
          ENDDO
        END IF
C
C       Begin sorting
C
        DO INT = 1,NINTS
          I = IND(1,INT)
          J = IND(2,INT)
          K = IND(3,INT)
          L = IND(4,INT)
          IF(I.GT.J) THEN
C
C*** Group B : IJKK
C
            NB(IREP)           = NB(IREP) + 1
            BCMAT(NB(IREP),1)  = BUF(INT)
            IND1         = NBIT1*I + J
            IND2         = K
            IF(NIBUF.EQ.1) THEN
              IBCMAT(NB(IREP),1,1) = IND1*NBIT2+IND2
            ELSE
              IBCMAT(NB(IREP),1,1) = IND1
              IBCMAT(NB(IREP),2,1) = IND2
            ENDIF
            IF(NB(IREP).EQ.MAXBUF) THEN
              CALL WRINT(LUINTA,NB(IREP),IREP,MAXBUF,NIBUF,
     &                   IBCMAT(1,1,1),BCMAT(1,1))
              ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
              NB(IREP) = 0
            ENDIF
          ELSEIF(K.GT.L) THEN
C
C*** Group C : IIKL
C
            IND1         = NBIT1*K + L
            IND2         = I
            NC(IREP)           = NC(IREP) + 1
            BCMAT(NC(IREP),2)  = BUF(INT)
            IF(NIBUF.EQ.1) THEN
              IBCMAT(NC(IREP),1,2) = IND1*NBIT2+IND2
            ELSE
              IBCMAT(NC(IREP),1,2) = IND1
              IBCMAT(NC(IREP),2,2) = IND2
            ENDIF
            IF(NC(IREP).EQ.MAXBUF) THEN
              CALL WRINT(LUINTB,NC(IREP),IREP,MAXBUF,NIBUF,
     &                   IBCMAT(1,1,2),BCMAT(1,2))
              ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
              NC(IREP) = 0
            ENDIF
          ELSE
C
C*** Group A : IIKK
C
            NA        = NA + 1
            AMAT(NA)  = BUF(INT)
            IAMAT(NA) = NBIT1*I + K
          ENDIF
        ENDDO
        GOTO 10
   30   CONTINUE
C*** Empty matrices:
        IF(NB(IREP).GT.0) THEN
          CALL WRINT(LUINTA,NB(IREP),IREP,MAXBUF,NIBUF,
     &               IBCMAT(1,1,1),BCMAT(1,1))
          NB(IREP) = NB(IREP) + ICOUNT(IREP,1)*MAXBUF
          ICOUNT(IREP,1) = ICOUNT(IREP,1) + 1
        ELSE
          NB(IREP) = ICOUNT(IREP,1)*MAXBUF
        ENDIF
        IF(NC(IREP).GT.0) THEN
          CALL WRINT(LUINTB,NC(IREP),IREP,MAXBUF,NIBUF,
     &              IBCMAT(1,1,2),BCMAT(1,2))
          NC(IREP) = NC(IREP) + ICOUNT(IREP,2)*MAXBUF
          ICOUNT(IREP,2) = ICOUNT(IREP,2) + 1
        ELSE
          NC(IREP) = ICOUNT(IREP,2)*MAXBUF
        ENDIF
      ENDDO
C
      WRITE(LUINTA) -1,0
      WRITE(LUINTB) -1,0
      CALL WRINT(LUINTB,NA,0,NSLDIM,1,IAMAT,AMAT)
C
C     Output section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(1X,A)') '* Integral statistics - fine sort:'
        WRITE(LUPRI,'(1X,A)') '  - Class C integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NC(IREP),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)') '  - Class B integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       ('Irrep ',REP(IREP),NB(IREP),IREP=0,MAXREP)
        WRITE(LUPRI,'(1X,A)')  ' - Class A integrals:'
        WRITE(LUPRI,'(8X,A,A3,I16)') 
     &       'Irrep ',REP(0),NA
      ENDIF
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twogen */
      SUBROUTINE TWOGEN(WORK,LWORK)
C*****************************************************************************
C
C     Do all necessary preparations for two-electron integrals
C
C     Written by T.Saue May 18 1995
C     Last revision : Sep 30 1996 - tsaue
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "cbihr2.h"
#include "dcbbas.h"
#include "dcbgen.h"
      DIMENSION WORK(LWORK)
#include "memint.h"
C
C     Memory allocation
C
      CALL MEMGET2('INTE','IREPA',KIREPA,NTBAS(0),WORK,KFREE,LFREE)
C
      CALL TWOGE1(WORK(KIREPA),WORK(KFREE),LFREE)
C
C     Memory deallocation
C
      CALL MEMREL('TWOGEN',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck twoge1 */
      SUBROUTINE TWOGE1(IREPA,WORK,LWORK)
C*****************************************************************************
C
C*****************************************************************************
#include "implicit.h"
      INTEGER IREPA(*), LWORK
      REAL(8) WORK(LWORK)
#include "priunit.h"
#include "dummy.h"
#include "maxaqn.h"
C
#include "ccom.h"
#include "dcbgen.h"
#include "dcbham.h"
#include "cbihr2.h"
#include "hrunit.h"
      LOGICAL LBIT,TOBE,REGO
      real(8) :: TIMSTR(2), TIMEND(2)

C     Undifferentiated integrals written on file
      INTTYP = 0
      NCHUNK = 600
C
C     Make array of boson irrep for each basis function
C
      CALL RPGEN(IREPA)
C
C     LL - integrals
C
 10   CONTINUE
      IF(.NOT.LBIT(IDFLAG,1).AND.LBIT(INTGEN,1)) THEN
C
C       Check for restart
C
cjkp+hjaaj: do not allow restart; we may have a new geometry
c       INQUIRE(FILE='DFLLSA',EXIST=TOBE)
c       REGO = TOBE
c       INQUIRE(FILE='DFLLSB',EXIST=TOBE)
c       REGO = TOBE.AND.REGO
c       IF(REGO) THEN
c         WRITE(LUPRI,'(A)')
c    &     '*** INFO *** : TWOGEN - LL-integrals from restart'
c         GOTO 20
c       ENDIF
cjkp+hjaaj end
C
C       Generate integrals
C
        I2TYP = 1
CTROND        CALL OPNFIL(LUINTA,'DFTWLL','UNKNOWN','TWOGE1')          
        CALL GPOPEN(LUINTA,'DFTWLL','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
        CALL TIMER2('START ',TIMSTR,TIMEND)
        CALL TWODRV(INTTYP,I2TYP,WORK,LWORK,IPRTWO)
        CALL TIMER2('LL-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
C
C     SL - integrals
C
 20   CONTINUE
      IF(.NOT.LBIT(IDFLAG,2).AND.LBIT(INTGEN,2)) THEN
C
C       Check for restart
C
        INQUIRE(FILE='DFSLSA',EXIST=TOBE)
        REGO = TOBE
        INQUIRE(FILE='DFSLSB',EXIST=TOBE)
        REGO = TOBE.AND.REGO
C
C       IF one-center model 2,3,4 is used for geometry optimization,
C       do not generate integrals since these are one-center
C       and therefore unchanged during minimization.   /jkp
C
        IF(REGO) THEN
          IF ((ONECAP.AND.INTV1C.NE.1).AND.OPTIMI) THEN
            WRITE(LUPRI,'(/2A)')
     &      '** One-center LS-integrals unchanged in this geometry -->',
     &      '   Read from disk.'
            GOTO 30
cjkp+hjaaj: do not allow restart; we may have a new geometry
c         ELSE 
c            WRITE(LUPRI,'(A)')
c    &      '*** INFO *** : TWOGEN - SL-integrals from restart'
          ENDIF
        ENDIF
C
C       Generate integrals
C
        I2TYP  = 2
        IF(IPRTWO.GE.2)
     +     WRITE(LUPRI,'(1P,A,D9.2)') 'SL-int - threshhold:',THRS
        CALL GPOPEN(LUINTA,'DFTWSL','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
        CALL TIMER2('START ',TIMSTR,TIMEND)
        CALL TWODRV(INTTYP,I2TYP,WORK,LWORK,IPRTWO)
        CALL TIMER2('SL-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
C
C     SS - integrals
C
 30   CONTINUE
      IF(.NOT.LBIT(IDFLAG,3).AND.LBIT(INTGEN,3)) THEN
C
C       Check for restart
C
C        INQUIRE(FILE='DFSSTA',EXIST=TOBE)
C        REGO = TOBE
C        INQUIRE(FILE='DFSSTB',EXIST=TOBE)
C        REGO = TOBE.AND.REGO
        INQUIRE(FILE='DFSSSA',EXIST=TOBE)
        REGO = TOBE
        INQUIRE(FILE='DFSSSB',EXIST=TOBE)
        REGO = TOBE.AND.REGO
C
C       IF one-center model 2,3,4 is used for geometry optimization,
C       do not generate integrals since these are one-center
C       and therefore unchanged during minimization.   /jkp
C
        IF(REGO) THEN
          IF ((ONECAP.AND.INTV1C.NE.1).AND.OPTIMI) THEN
            WRITE(LUPRI,'(/2A)')
     &      '** One-center SS-integrals unchanged in this geometry -->',
     &      '   Read from disk.'
            GOTO 40
cjkp+hjaaj: do not allow restart; we may have a new geometry
c         ELSE 
c           WRITE(LUPRI,'(A)')
c    &      '*** INFO *** : TWOGEN - SS-integrals from restart'
          ENDIF
        ENDIF
C
C       Generate integrals
C
        I2TYP  = 3
        IF(IPRTWO.GE.2)
     +     WRITE(LUPRI,'(A,E9.2)') 'SS-int - threshhold:',THRS
        CALL GPOPEN(LUINTA,'DFTWSS','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
        CALL TIMER2('START ',TIMSTR,TIMEND)
        CALL TWODRV(INTTYP,I2TYP,WORK,LWORK,IPRTWO)
        CALL TIMER2('SS-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
C
C     Gaunt - integrals
C
 40   CONTINUE
      IF(.NOT.LBIT(IDFLAG,4).AND.LBIT(INTGEN,4)) THEN
C
C       Generate integrals
C
        I2TYP  = 4
        IF(IPRTWO.GE.2)
     +     WRITE(LUPRI,'(A,E9.2)') 'GT-int - threshhold:',THRS
        CALL GPOPEN(LUINTA,'DFTWGT','UNKNOWN',' ',' ',IDUMMY,.FALSE.)
        CALL TIMER2('START ',TIMSTR,TIMEND)
        CALL TWODRV(INTTYP,I2TYP,WORK,LWORK,IPRTWO)
        CALL TIMER2('GT-int',TIMSTR,TIMEND)
        CALL FLSHFO(LUPRI)
        CALL TWOSRT(I2TYP,NCHUNK,IREPA,WORK,LWORK,IPRTWO)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xx2mat */
      SUBROUTINE XX2MAT(LUINTA,LUINTB,NDIM,NNDIM,
     +  FMAT,DMAT,NDMAT,IFCTYP,IREPDM,BUF,IBUF,IND,IPRINT)
C*****************************************************************************
C
C     Insert LL/SS-integrals in Fock matrix
C
C     Written by T.Saue Mar 3 1995
C     Last revision : tsaue - Mar 3 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
      DIMENSION IBUF(*),BUF(*),FMAT(*),DMAT(*),IND(*),IFCTYP(*),
     &          IREPDM(*)
C
C     *********************************************
C     ***** Process integrals on file LUINTA ******
C     *********************************************
C
 10   CONTINUE
      REWIND LUINTA
C
C     Class D integrals:
C     ==================
C        I > K > J > L
C        I > K > L > L
C        I > J > K > L
C
      CALL DOIJKL(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class Ca integrals:
C     ===================
C        I = K > J > L
C        I > J = K > L
C        I > K > J = L
C
      CALL DOIJIL(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     *********************************************
C     ***** Process integrals on file LUINTB ******
C     *********************************************
C
 20   CONTINUE
      REWIND LUINTB
C
C     Class Cb integrals:
C     ===================
C        I > K = L > J
C        I > J > K = L
C        I = J > K > L
C
      CALL DOIJKK(LUINTB,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class Bb2 integrals:
C     ====================
C        I = K > J = L
C
      DO IREP = 1,NBSYM
        CALL DOIJIJ(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
      ENDDO
C
C     Class Ba1 integrals:
C     ====================
C       I > J = K = L
C
      CALL DOIJJJ(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
C
C     Class Ba2 integrals:
C     ====================
C        I = J = K > L
C
      CALL DOIJJJ(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
C
C     Class Bb1 integrals:
C     ====================
C        I = J > K = L
C
      CALL DOIIKK(LUINTB,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IND,IPRINT)
C
C     Class A
C     =======
C
      CALL DOIIII(LUINTB,NDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &            IREPDM,IBUF,BUF,IPRINT)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck xy2mat */
      SUBROUTINE XY2MAT(LUINTA,LUINTB,
     +  FMAT,DMAT,NDMAT,IFCTYP,IREPDM,BUF,IBUF,IND,IPRINT)
C*****************************************************************************
C
C     Insert LL/SS-integrals in Fock matrix
C
C     Written by T.Saue Mar 3 1995
C     Last revision : tsaue - Mar 3 1995
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
      DIMENSION IBUF(*),BUF(*),FMAT(*),DMAT(*),IND(*),
     &          IFCTYP(*),IREPDM(*)
C
C     *********************************************
C     ***** Process integrals on file LUINTA ******
C     *********************************************
C
      REWIND LUINTA
C
C     Class D integrals:
C     ==================
C       I > J, K > L
C
      CALL DOIJKL(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class B integrals:
C     ==================
C       I > J, K = L
C
      CALL DOIJKK(LUINTA,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     *********************************************
C     ***** Process integrals on file LUINTB ******
C     *********************************************
C
      REWIND LUINTB
C
C     Class C integrals:
C     ==================
C        I = J, K > L
C
      CALL DOIJKK(LUINTB,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &            IBUF,BUF,IND,IPRINT)
C
C     Class A integrals:
C     ==================
C        I = J, K = L
C
      NSLDIM = NTBAS(1)*NTBAS(2)
      CALL DOIIKK(LUINTB,NSLDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &              IREPDM,IBUF,BUF,IND,IPRINT)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijkl */
      SUBROUTINE DOIJKL(IUNIT,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &                  IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals with all indices unequal
C
C     Written by T.Saue Mar 5 1995
C     Last revision: Jan 7 1998 - jth (Linux cheat)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
C
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,*),DMAT(NBAST,NBAST,*) 
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),*),DMAT(NTBAS(0),NTBAS(0),*) 
#endif
      DIMENSION BUF(N2BBASX),IBUF(N2BBASX,NIBUF),IND(4,N2BBASX),
     &          IFCTYP(NDMAT),IREPDM(NDMAT)
      INTEGER A,B,C,D
#include "ibtfun.h"
C
   10   CONTINUE
        CALL REAINT(IUNIT,NBUF,IREPAB,N2BBASX,NIBUF,IBUF,BUF)
        IF(NBUF.EQ.-1) RETURN
C
C       Unpack indices
C       ==============
C
        IF (NIBUF .EQ. 1) THEN
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),3*NBITS),IBIT1)
            IND(2,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),2*NBITS),IBIT1)
            IND(3,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),  NBITS),IBIT1)
            IND(4,INTS) = IBTAND(       IBUF(INTS,1),         IBIT1)
          ENDDO
        ELSE
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),NBITS),IBIT1)
            IND(2,INTS) = IBTAND(       IBUF(INTS,1),       IBIT1)
            IND(3,INTS) = IBTAND(IBTSHR(IBUF(INTS,2),NBITS),IBIT1)
            IND(4,INTS) = IBTAND(       IBUF(INTS,2),       IBIT1)
          ENDDO
        END IF
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJKL: Integrals'
        WRITE(LUPRI,'(4I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),IND(3,I),IND(4,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,D,I)
              FMAT(C,D,I) = FMAT(C,D,I) + DINT*DMAT(A,B,I)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,D,I)
              FMAT(C,D,I) = FMAT(C,D,I) + DINT*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKL: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              GCD  = D2*BUF(INTS)*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GCD
              FMAT(B,A,I) = FMAT(B,A,I) + GCD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I) = FMAT(C,D,I) + GAB
              FMAT(D,C,I) = FMAT(D,C,I) + GAB
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              D = IND(4,INTS)
              GCD  = D2*BUF(INTS)*(DMAT(C,D,I) + DMAT(D,C,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GCD
              FMAT(B,A,I) = FMAT(B,A,I) + GCD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,D,I) = FMAT(C,D,I) + GAB
              FMAT(D,C,I) = FMAT(D,C,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKL: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJKL ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJKL: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      GOTO 10
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijil */
      SUBROUTINE DOIJIL(IUNIT,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &                  IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJIL
C     Indices are bit packed with the degenerate index as index 3.
C
C     Written by T.Saue Mar 4 1995
C     Last revision Mar 4 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D2 = 2.0D0,D4 = 4.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT) 
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT) 
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(N2BBASX),IBUF(N2BBASX,NIBUF),IND(3,N2BBASX)
      INTEGER A,B,D
#include "ibtfun.h"
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJIL: IFCTYP:'
        WRITE(LUPRI,'(5X,I5)') (IFCTYP(I),I=1,NDMAT)
      ENDIF        
   10   CONTINUE
        CALL REAINT(IUNIT,NBUF,IREPAB,N2BBASX,NIBUF,IBUF,BUF)
        IF(NBUF.EQ.-1) RETURN
C
C       Unpack indices
C       ==============
C
        IF (NIBUF .EQ. 1) THEN
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),3*NBITS),IBIT1)
            IND(2,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),2*NBITS),IBIT1)
            IND(3,INTS) = IBTAND(       IBUF(INTS,1),         IBIT1)
          ENDDO
        ELSE
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),NBITS),IBIT1)
            IND(2,INTS) = IBTAND(       IBUF(INTS,1),       IBIT1)
            IND(3,INTS) = IBUF(INTS,2)
          ENDDO
        END IF
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJIL: Integrals'
        WRITE(LUPRI,'(3I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),IND(3,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(A,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) + DINT*DMAT(A,B,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              DINT = D4*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(A,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) + DINT*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJIL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJIL: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              GAD  = D2*BUF(INTS)*(DMAT(A,D,I) + DMAT(D,A,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GAD
              FMAT(B,A,I) = FMAT(B,A,I) + GAD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(A,D,I) = FMAT(A,D,I) + GAB
              FMAT(D,A,I) = FMAT(D,A,I) + GAB
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(D,B,I)
              FMAT(D,A,I) = FMAT(D,A,I) - BUF(INTS)*DMAT(A,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(D,A,I)
              FMAT(D,B,I) = FMAT(D,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,D,I)
              FMAT(A,D,I) = FMAT(A,D,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(A,D,I)
              FMAT(B,D,I) = FMAT(B,D,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(3,INTS)
              B = IND(1,INTS)
              D = IND(2,INTS)
              GAD  = D2*BUF(INTS)*(DMAT(A,D,I) + DMAT(D,A,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GAD
              FMAT(B,A,I) = FMAT(B,A,I) + GAD
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(A,D,I) = FMAT(A,D,I) + GAB
              FMAT(D,A,I) = FMAT(D,A,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJIL ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJIL: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJIL ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJIL: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      GOTO 10
CC
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijkk */
      SUBROUTINE DOIJKK(IUNIT,FMAT,DMAT,NDMAT,IFCTYP,IREPDM,
     &                  IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJKK
C     Indices are bit packed with the degenerate index as index 3.
C
C     Written by T.Saue Mar 4 1995
C     Last revision Mar 4 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(N2BBASX),IBUF(N2BBASX,NIBUF),IND(3,N2BBASX)
      INTEGER A,B,C
#include "ibtfun.h"
C
 10     CONTINUE
        CALL REAINT(IUNIT,NBUF,IREPAB,N2BBASX,NIBUF,IBUF,BUF)
        IF(NBUF.EQ.-1) RETURN
C
C       Unpack indices
C       ==============
C
        IF (NIBUF .EQ. 1) THEN
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),3*NBITS),IBIT1)
            IND(2,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),2*NBITS),IBIT1)
            IND(3,INTS) = IBTAND(       IBUF(INTS,1),         IBIT1)
          ENDDO
        ELSE
          DO INTS = 1,NBUF
            IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS,1),NBITS),IBIT1)
            IND(2,INTS) = IBTAND(       IBUF(INTS,1),       IBIT1)
            IND(3,INTS) = IBUF(INTS,2)
          ENDDO
        END IF
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIJKK: Integrals'
        WRITE(LUPRI,'(3I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),IND(3,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (1/4) * (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              DINT = D2*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,C,I)
              FMAT(C,C,I) = FMAT(C,C,I) + DINT*DMAT(A,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (1/4) (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              DINT = D2*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(C,C,I)
              FMAT(C,C,I) = FMAT(C,C,I) + DINT*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKK: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              GCC  = D2*BUF(INTS)*DMAT(C,C,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GCC
              FMAT(B,A,I) = FMAT(B,A,I) + GCC
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,C,I) = FMAT(C,C,I) + GAB
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              FMAT(C,A,I) = FMAT(C,A,I) - BUF(INTS)*DMAT(C,B,I)
              FMAT(C,B,I) = FMAT(C,B,I) - BUF(INTS)*DMAT(C,A,I)
              FMAT(A,C,I) = FMAT(A,C,I) - BUF(INTS)*DMAT(B,C,I)
              FMAT(B,C,I) = FMAT(B,C,I) - BUF(INTS)*DMAT(A,C,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = (1/8) * FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              C = IND(3,INTS)
              GCC  = D2*BUF(INTS)*DMAT(C,C,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GCC
              FMAT(B,A,I) = FMAT(B,A,I) + GCC
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              FMAT(C,C,I) = FMAT(C,C,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJKK: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJKK ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJKK: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
      GOTO 10
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijjj */
      SUBROUTINE DOIJJJ(IUNIT,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJJJ
C
C     Written by T.Saue Mar 5 1995
C     Last revision: tsaue Mar 5 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#include "dgroup.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(NNDIM),IBUF(NNDIM),IND(2,NNDIM)
      INTEGER A,B
#include "ibtfun.h"
        CALL REAINT(IUNIT,NBUF,IREPAB,NNDIM,1,IBUF,BUF)
C
C       Unpack indices
C       ==============
C
        DO INTS = 1,NBUF
          IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS),NBITS),IBIT1)
          IND(2,INTS) = IBTAND(       IBUF(INTS),       IBIT1)
        ENDDO
C
C       Print integrals
C       ===============
C
        IF(IPRINT.GE.6) THEN
          WRITE(LUPRI,'(A)') '* DOIJJJ: Integrals'
          IF(NBUF.GT.0) WRITE(LUPRI,'(2I5,3X,F12.6)')
     &       (IND(1,I),IND(2,I),BUF(I),I=1,NBUF)
        ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + BUF(INTS)*DMAT(B,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,B,I) = FMAT(B,B,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(B,B,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
           ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              DINT = D2*BUF(INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + DINT*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + DINT*DMAT(B,A,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJJJ ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJJJ: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              GBB  = BUF(INTS)*DMAT(B,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GBB
              FMAT(B,A,I) = FMAT(B,A,I) + GBB
              FMAT(B,B,I) = FMAT(B,B,I) 
     &                    + BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,B,I) = FMAT(B,B,I) 
     &                    - BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              GBB  = BUF(INTS)*DMAT(B,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) - GBB
              FMAT(B,A,I) = FMAT(B,A,I) - GBB
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,B,I) = FMAT(B,B,I) 
     &                    + D2*BUF(INTS)*(DMAT(A,B,I) + DMAT(B,A,I))
              GBB  = D2*BUF(INTS)*DMAT(B,B,I)
              FMAT(A,B,I) = FMAT(A,B,I) + GBB
              FMAT(B,A,I) = FMAT(B,A,I) + GBB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIJJJ ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIJJJ: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJJJ ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIJJJ: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doiikk */
      SUBROUTINE DOIIKK(IUNIT,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IIKK
C
C     Written by T.Saue Mar 5 1995
C     Last revision: tsaue Mar 5 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(NNDIM),IBUF(NNDIM),IND(2,NNDIM)
      INTEGER A,B
#include "ibtfun.h"
C
      CALL REAINT(IUNIT,NBUF,IREPAB,NNDIM,1,IBUF,BUF)
C
C     Unpack indices
C     ==============
C
      DO INTS = 1,NBUF
        IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS),NBITS),IBIT1)
        IND(2,INTS) = IBTAND(       IBUF(INTS),       IBIT1)
      ENDDO
C
C     Print integrals
C     ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIIKK: Integrals'
        WRITE(LUPRI,'(2I5,3X,F12.6)')
     &     (IND(1,I),IND(2,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
       DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + BUF(INTS)*DMAT(A,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIKK: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + D2*BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + D2*BUF(INTS)*DMAT(A,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + D2*BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + D2*BUF(INTS)*DMAT(A,A,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIKK ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIKK: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIIKK ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIIKK: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doijij */
      SUBROUTINE DOIJIJ(IUNIT,NNDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IBUF,BUF,IND,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IJIJ
C
C     Written by T.Saue Mar 5 1995
C     Sep 3 1996 - tsaue   Generalize to IFCTYP
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0,D3 = 3.0D0,DP5 = 0.5D0)
#include "dcbbas.h"
#include "dcbind.h"
#include "dgroup.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION IFCTYP(NDMAT),IREPDM(NDMAT),
     &          BUF(NNDIM),IBUF(NNDIM),IND(2,NNDIM)
      INTEGER A,B
#include "ibtfun.h"
C
      CALL REAINT(IUNIT,NBUF,IREPAB,NNDIM,1,IBUF,BUF)
C
C     Unpack indices
C     ==============
C
      DO INTS = 1,NBUF
        IND(1,INTS) = IBTAND(IBTSHR(IBUF(INTS),NBITS),IBIT1)
        IND(2,INTS) = IBTAND(       IBUF(INTS),       IBIT1)
      ENDDO
C
C     Print integrals
C     ===============
C
      IF(IPRINT.GE.6.AND.NBUF.GT.0) THEN
        WRITE(LUPRI,'(A)') '* DOIJIJ: Integrals'
        WRITE(LUPRI,'(2I5,3X,F12.6)')
     &   (IND(1,I),IND(2,I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C     IFCTYP = XY
C       X indicates symmetry about diagonal
C         X = 0 No symmetry
C         X = 1 Symmetric
C         X = 2 Anti-symmetric
C       Y indicates contributions
C         Y = 1 Coulomb
C         Y = 2 Exchange
C         Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C       Symmetric singlet Fock matrix
C       =============================
C       F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + D3*BUF(INTS)*DMAT(A,B,I)
              EINT = -DP5*BUF(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + EINT*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + EINT*DMAT(A,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(B,A,I)
              EINT = -DP5*BUF(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + EINT*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) + EINT*DMAT(A,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) + D4*BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJIJ ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
            CALL QUIT(
     &      'ERROR in DOIJIJ: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              GAB  = BUF(INTS)*DMAT(A,B,I)
              GBA  = BUF(INTS)*DMAT(B,A,I)
              FMAT(A,B,I) = FMAT(A,B,I) + D2*GAB + GBA
              FMAT(B,A,I) = FMAT(A,B,I) + D2*GBA + GAB
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              FMAT(A,B,I) = FMAT(A,B,I) - BUF(INTS)*DMAT(B,A,I)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(B,B,I)
              FMAT(B,B,I) = FMAT(B,B,I) - BUF(INTS)*DMAT(A,A,I)
              FMAT(B,A,I) = FMAT(B,A,I) - BUF(INTS)*DMAT(A,B,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(1,INTS)
              B = IND(2,INTS)
              GAB  = D2*BUF(INTS)*(DMAT(A,B,I)+DMAT(B,A,I))
              FMAT(A,B,I) = FMAT(A,B,I) + GAB
              FMAT(B,A,I) = FMAT(A,B,I) + GAB
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIJIJ ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &     'ERROR in DOIJIJ: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &    'DOIJIJ ERROR, specified IFCTYP not implemented yet',
     &    '              specified IFCTYP was',IFCTYP(I),
     &    '              for F,D matrix no.  ',I
          CALL QUIT(
     &    'ERROR in DOIJIJ: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck doiiii */
      SUBROUTINE DOIIII(IUNIT,NDIM,FMAT,DMAT,NDMAT,IFCTYP,
     &                  IREPDM,IND,BUF,IPRINT)
C*****************************************************************************
C
C     Process integrals of the form IIII
C
C     Written by T.Saue Mar 5 1995
C     Sep 3 1996 - tsaue   Genrazlized to IFCTYP
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER (D4 = 4.0D0,D2 = 2.0D0,DP5=0.5D0)
#include "dcbbas.h"
#include "dcbind.h"
#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION FMAT(NBAST,NBAST,NDMAT),
     &          DMAT(NBAST,NBAST,NDMAT)
#else
      DIMENSION FMAT(NTBAS(0),NTBAS(0),NDMAT),
     &          DMAT(NTBAS(0),NTBAS(0),NDMAT)
#endif
      DIMENSION BUF(NDIM),IND(NDIM),IFCTYP(NDMAT),
     &          IREPDM(NDMAT)
      INTEGER A
      CALL REAINT(IUNIT,NBUF,IREPAB,NDIM,1,IND,BUF)
C
C       Print integrals
C       ===============
C
      IF(IPRINT.GE.6) THEN
        WRITE(LUPRI,'(A)') '* DOIIII: Integrals'
        WRITE(LUPRI,'(I5,3X,F12.6)')
     &     (IND(I),BUF(I),I=1,NBUF)
      ENDIF        
C
C     **********************************
C     ***** Contract Fock matrices *****
C     **********************************
C       IFCTYP = XY
C         X indicates symmetry about diagonal
C           X = 0 No symmetry
C           X = 1 Symmetric
C           X = 2 Anti-symmetric
C         Y indicates contributions
C           Y = 1 Coulomb
C           Y = 2 Exchange
C           Y = 3 Coulomb + Exchange
C
      DO 400 I = 1,NDMAT
        IY = MOD(IFCTYP(I),10)
        IX = (IFCTYP(I)-IY)/10
        IC = MOD(IY,2)
        IE = IY - IC
        IF(IREPAB.NE.IREPDM(I)) IC = 0
        IY = IE + IC
        IF(IY.EQ.0) GOTO 400
        IF    (IX.EQ.1.OR.IX.EQ.2.OR.IX.EQ.3) THEN
C
C         Symmetric singlet Fock matrix
C         =============================
C         F(i,j) = (FMAT(i,j) + FMAT(j,i))
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + DP5*BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         Antisymmetric singlet Fock matrix OR
C         symmetric triplet Fock matrix OR
C         antisymmetric triplet Fock matrix
C         =========================================
C         F(i,j) = (FMAT(i,j) +/- FMAT(j,i))
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - DP5*BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         Coulomb contributions only
C         ==========================
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(A,A,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIII ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIII: specified IFCTYP not implemented.')
          ENDIF
        ELSEIF(IX.EQ.0) THEN
C
C         General singlet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          IF(IY.EQ.3) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General triplet case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.2) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) - BUF(INTS)*DMAT(A,A,I)
            ENDDO
C
C         General Coulomb case - no permutational symmetry
C         ================================================
C         F(i,j) = FMAT(i,j)
C
          ELSEIF(IY.EQ.1) THEN
            DO INTS = 1,NBUF
              A = IND(INTS)
              FMAT(A,A,I) = FMAT(A,A,I) + D2*BUF(INTS)*DMAT(A,A,I)
            ENDDO
          ELSE
            WRITE (LUPRI,'(/A,2(/A,I10))')
     &        'DOIIII ERROR, specified IFCTYP not implemented yet',
     &        '              specified IFCTYP was',IFCTYP(I),
     &        '              for F,D matrix no.  ',I
            CALL QUIT(
     &        'ERROR in DOIIII: specified IFCTYP not implemented.')
          ENDIF
        ELSE
          WRITE (LUPRI,'(/A,2(/A,I10))')
     &      'DOIIII ERROR, specified IFCTYP not implemented yet',
     &      '              specified IFCTYP was',IFCTYP(I),
     &      '              for F,D matrix no.  ',I
          CALL QUIT(
     &      'ERROR in DOIIII: specified IFCTYP not implemented.')
        ENDIF
 400  CONTINUE
C
      END
CC&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck fastjk */
      SUBROUTINE FASTJK(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,IFCKOP,
     &              NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct construction of skeleton Fock-matrix
C
C     Written by T.Saue, S.Reine and A.Gomes June 15 2006
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
      PARAMETER(D0 = 0.0D0,DP5 = 0.5D0,D2 = 2.0D0,D10 = 10.0D0)
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),DMAT(N2BBASX,NZ,NFMAT),
     &        ISYMOP(*),IHRMOP(*),IFCKOP(*),
     &        NPOS(*),WORK(LWORK), IFCKOP_TEMP(NFMAT)
C
!     DIM = N2BBASX*NZC1*NFMAT
C       Conventional (4-center integrals) AO Coulomb + exchange
        CALL DIRAOF(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &              IFCKOP,NPOS,INTFLG,
     &              IPRINT,WORK,LWORK)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck diraof */
      SUBROUTINE DIRAOF(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,IFCKOP,
     &                  NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations in AO-basis
C
C     Written by T.Saue July 2 1995
C     tsaue -    Sep 2 1996 Generalized to NFMAT 
C                           and general operator symmetry IRPOP
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
      LOGICAL LBIT
#include "maxorb.h"
#include "aovec.h"
#include "dcbham.h"
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),DMAT(N2BBASX,NZ,NFMAT),
     &          ISYMOP(*),IHRMOP(*),IFCKOP(*),
     &          NPOS(*),WORK(LWORK)
#include "dgroup.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "blocks.h"
#include "cbihr2.h"
#include "ibtfun.h"
C
      CALL QENTER('DIRAOF')
#include "memint.h"
C
C     Memory allocation
      NDMAT   = NFMAT*NZC1
      N2GAB   = NSYMBL*NSYMBL
      NDMRAO  = N2GAB*NDMAT
      CALL MEMGET2('INTE','IFC'  ,KIFC  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IRD'  ,KIRD  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','GAB'  ,KGAB  ,N2GAB  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DMRAO',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
C
C
C     Get G(ab) integrals from AOPROPER
C     =================================
C
      IF (SCRFCK .GT. D0) THEN
C
         IJOB   = 0
         ITYPE  = 0
         IGTYP  = 0
C        IGTYP 0 gives just the LL integrals
C        IGTYP 1 will add the SS type integrals
         IF(LBIT(INTFLG,2).OR.LBIT(INTFLG,3)) IGTYP = 1
C        IGTYP 2 will add also the SL (Gaunt) type integrals
         IF(LBIT(INTFLG,4)) IGTYP = 2
         MAXDIF = 0
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDIF,
     &        IPRINT,WORK(KGAB),WORK(KFREE),LFREE)
      ENDIF
C
      IF(NBSYM.GT.1) THEN
        IF(MDIRAC) THEN
          CALL QUIT('DIRAOF: Dont know how to handle symmetry yet..')
        ENDIF
        CALL MEMGET2('REAL','FAO',KFAO,N2BBASX*NDMAT,WORK,KFREE,LFREE)
        CALL MEMGET2('REAL','DAO',KDAO,N2BBASX*NDMAT,WORK,KFREE,LFREE)
C
C       Transform density matrix from QO to AO
C       ======================================
C
        IOFF = 0
        DO I = 1,NFMAT
          IREP = ISYMOP(I)-1
          IPAR = JBTOF(IREP,1)
          DO IZ = 1,NZC1
            IREPD = IRQMAT(IZ,IREP)
            IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)            
            IPQ   = IQTOPQ(IQ,IREP)
            IF(IPQ.GT.NZ) CALL QUIT('DIRAOF: IPQ.GT.NZ !')
            IF(IPQ.LE.0 ) CALL QUIT('DIRAOF: IPQ.LE.0 !')
            CALL DTSOAO(DMAT(1,IPQ,I),WORK(KDAO+IOFF),
     &                 NTBAS(0),IREPD,IPRTWO)
            IOFF = IOFF + N2BBASX
          ENDDO
          IF(IPRINT.GE.5) THEN
            JOFF = KDAO + (I-1)*4*N2BBASX
            WRITE(LUPRI,'(A,I5)') '** Matrix no.',I
            CALL HEADER('Backtransformed hph. density matrix',-1)
            CALL PRQMAT(WORK(JOFF),NTBAS(0),NTBAS(0),NTBAS(0),
     &                  NTBAS(0),NZC1,IQMULT,LUPRI)
          ENDIF
        ENDDO
C
C       Construct Fock matrix
C       ======================
C
        CALL DZERO(WORK(KFAO),N2BBASX*NDMAT)
        CALL DIRAOF_1(WORK(KFAO),WORK(KDAO),NFMAT,
     &              WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &              WORK(KDMRAO),WORK(KGAB),NPOS,INTFLG,
     &              IPRINT,WORK(KFREE),LFREE)
        IF(IPRINT.GE.5) THEN
          WRITE(LUPRI,'(A)') 'After DIRAOF_1'
          CALL PRQMAT(WORK(KFAO),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NTBAS(0),NZC1,IQMULT,LUPRI)
        ENDIF
C
C       Transform Fock matrix from AO to QO
C       ===================================
C
        IOFF = 0
        DO I = 1,NFMAT
          IREP = ISYMOP(I)-1
          IPAR = JBTOF(IREP,1)
          DO IZ = 1,NZC1
            IREPD = IRQMAT(IZ,IREP)
            IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)            
            IPQ   = IQTOPQ(IQ,IREP)
            CALL MTAOSO(WORK(KFAO+IOFF),FMAT(1,IPQ,I),
     &                 NTBAS(0),IREPD,IPRTWO)
            IOFF = IOFF + N2BBASX
          ENDDO
        ENDDO
      ELSE
C
C       If there is no symmetry involved the density matrix DMAT
C       is sent to DIRAOF_1 module:
C
        CALL DIRAOF_1(FMAT,DMAT,NFMAT,
     &              WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &              WORK(KDMRAO),WORK(KGAB),NPOS,INTFLG,
     &              IPRINT,WORK(KFREE),LFREE)
      ENDIF
      CALL MEMREL('DIRAOF.DIRAOF_1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('DIRAOF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck DIRAOF_1 */
      SUBROUTINE DIRAOF_1(FAO,DAO,NFMAT,
     &                  IFCTYP,IREPDM,ISYMOP,IHRMOP,IFCKOP,
     &                  DMRAO,GABRAO,NPOS,INTFLG,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations in AO-basis
C
C        FAO - Fock matrix in AO-basis
C        DAO - density matrix in AO-basis
C        IFCTYP - array indicating Fock-matrix type, see routine SETFCK
C        ISYMOP - array indicating spatial symmetry, see routine SETFCK
C        IHRMOP - array indicating Hermiticity, see routine SETFCK
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,DP5 = 0.5D0,D2 = 2.0D0,D10 = 10.0D0)
#include "dummy.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "mxcent.h"
C
C Used from COMMON blocks:
C   DCBGEN: PARCAL
C   DCBBAS: NTBAS(0),N2BBASX
C   CCOM  : THRS,NHTYP,KHK
C   SHELLS: KMAX
C   symmet.h : MAXREP
C
#include "dcbgen.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "ccom.h"
#include "shells.h"
#include "blocks.h"
#include "symmet.h"
#if defined (VAR_MPI)
#include "infpar.h"
#endif
      LOGICAL PERTUR, EXPECT, UNDIFF, SPNORB, DIRFCK, DISTRI, SQ12EL,
     &        LONDON, SUSCEP, DDFOCK, ADISTR, LBIT, FIRST
      CHARACTER SECTID*12, TOTTID*12, WALLTID*12, CPUTID*12,
     &          AOTEXT*9
      DIMENSION FAO(N2BBASX*4,NFMAT),DAO(N2BBASX*4,NFMAT),
     &          IFCTYP(4,NFMAT),IREPDM(4,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),
     &          NPOS(*),DMRAO(*),GABRAO(*),WORK(LWORK),DNTSKP(8)
C
      CALL QENTER('DIRAOF_1')
#include "memint.h"
C
Ctec Is FIRST needed???
      FIRST = .TRUE.
      NDMAT = NFMAT*NZC1
      CALL SETFCK(IFCTYP,IREPDM,NFMAT,NZC1,ISYMOP,IHRMOP,IFCKOP,IPRINT)
      CALL DZERO(DNTSKP,8)

      IF (SCRFCK.GT.D0) THEN
C
C     Make reduced density matrix for screening
C     ==========================================
C     (Reduced GAB is already present in GABRAO)
C
         CALL MKDRAO(DAO,DMRAO,NDMAT,WORK,LWORK,IPRINT)
      ENDIF
C
C     ****************************
C     **** Process integrals *****
C     ****************************
C
      IRNTYP = 3
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      KPOS   = 1
      THRSBUF = THRS
      IF (PARCAL) THEN
         CALL SCRSTA('MPIHeader',DUM,DUM)
      ELSE
         CALL SCRSTA('Header',DUM,DUM)
      END IF
C
C....> Start looping over the four possible integral classes
C
      DO I2TYP = 1, 4
C
C....> Set the parameters that are special for each class
C
      IF (I2TYP.EQ.1) THEN
         MTOTTK = NLRGBL*(NLRGBL+1)/2
         THRS   = THRSBUF
         AOTEXT = 'AOfock:LL'
      ELSEIF (I2TYP.EQ.2) THEN
         MTOTTK = NLRGBL*(NLRGBL+1)/2
         THRS   = THRSBUF*THRFAC(1)
         AOTEXT = 'AOfock:SL'
      ELSEIF (I2TYP.EQ.3) THEN
         MTOTTK = NSMLBL*(NSMLBL+1)/2
         THRS   = THRSBUF*THRFAC(2)
         AOTEXT = 'AOfock:SS'
      ELSEIF (I2TYP.EQ.4) THEN
         MTOTTK = NSMLBL*NLRGBL
         THRS   = THRSBUF*THRFAC(1)
         AOTEXT = 'AOfock:GT'
C        Inform the Fock matrix builder that we're going to do Gaunt
         DO IFMAT = 1, NFMAT
C           The real matrix can only have an exchange contribution
            IX = IFCTYP(1,IFMAT)/10
            IY = MOD(IFCTYP(1,IFMAT),10)
            IC = MOD(IY,2)
            IE = (IY-IC)/2
            IFCTYP(1,IFMAT)=10*IX+4+IE
C           For the imaginary parts we'll just add 6 to signal that it is Gaunt
            IFCTYP(2,IFMAT)=IFCTYP(2,IFMAT)+6
            IFCTYP(3,IFMAT)=IFCTYP(3,IFMAT)+6
            IFCTYP(4,IFMAT)=IFCTYP(4,IFMAT)+6
         ENDDO
      ENDIF
C
C....> Check whether the class is active
C
      IF(LBIT(INTFLG,I2TYP).AND.LBIT(IDFLAG,I2TYP)) THEN
C
#if defined (VAR_MPI)
            IF (PARCAL) THEN
C
C             Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
              CALL DIRAC_PARCTL( HERFCK_PAR )
C
              CALL GETTIM(CPU1,WALL1)
              CALL HER_PARDRV(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,
     &             IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &             TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &             ICEDIF,SCRFCK,GABRAO,DMRAO,DUM,DNTSKP)
C
              CALL GETTIM(CPU2,WALL2)
              CPU    = CPU2 - CPU1
              WALL   = WALL2 - WALL1
              TOTWAL = TOTWAL + WALL
              CALL SCRSTA(AOTEXT,DNTSKP,WALL)
              IF (IPRINT .GT. 0) THEN
                 TOTTID = SECTID(TOTWAL)
                 WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &                TOTTID
                 IF (IPRINT .GT. 2) THEN
                    CPUTID = SECTID(CPU)
                    WALLTID = SECTID(WALL)
                    WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &                CPUTID,'/',WALLTID
                 END IF
              END IF
           ELSE
#endif
              CPU = SECOND()
              CALL TWOINT(WORK,LWORK,FAO,DAO,NDMAT,IDUM,IFCTYP,DUM,
     &             IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &             .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &             IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &             GABRAO,DMRAO,DUM,DNTSKP,.TRUE.,.false.,IDUM,DUM)
              CPU = SECOND() - CPU
              CALL SCRSTA(AOTEXT,DNTSKP,CPU)
#if defined (VAR_MPI)
           ENDIF
#endif
      ENDIF
C
C     Restore IFCTYP to the old value (even though we'll not need
C     it anymore)
      IF (I2TYP.EQ.4) THEN
         DO IFMAT = 1, NFMAT
            IFCTYP(1,IFMAT)=10*(IFCTYP(1,IFMAT)/10)+IY
            IFCTYP(2,IFMAT)=IFCTYP(2,IFMAT)-6
            IFCTYP(3,IFMAT)=IFCTYP(3,IFMAT)-6
            IFCTYP(4,IFMAT)=IFCTYP(4,IFMAT)-6
         ENDDO
      ENDIF
C
C---> Update the offset
C
      KPOS = KPOS + MTOTTK
C
C.... End of loop over integral classes
C
         CALL FLSHFO(LUPRI)
      ENDDO ! I2TYP = 1, 4
C
C     Symmetrize Fock matrix
C     ======================
C
!     CALL MEMGET2('REAL','SKLTON',KSKLTN,NMAT*N2BASX,
!    &   WORK,KFREE,LFREE)
!     save memory by reusing DAO for SKLTON matrix
      CALL MEMGET2('REAL','FAC',KFAC,NTBAS(0)*(MAXREP+1),
     &   WORK,KFREE,LFREE)
      NINDEX = KMAX*(MAXREP + 1)*KHK(NHTYP)
      CALL MEMGET2('INTE','IINDEX',KINDEX,NINDEX,
     &   WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IPOINT',KPOINT,NTBAS(0)*(MAXREP+1),
     &   WORK,KFREE,LFREE)
      IXYZ = 0
!     CALL SKLFCK_1(FAO,WORK(KSKLTN),WORK(KFAC),WORK(KINDEX),
      CALL SKLFCK_1(FAO,DAO         ,WORK(KFAC),WORK(KINDEX),
     &               WORK(KPOINT),DUM,DUM,DUM,IDUM,IXYZ,
     &               NDMAT,IREPDM,IFCTYP,IPRINT)
      CALL MEMREL('SKLFCK_1',WORK,KWORK,KWORK,KFREE,LFREE)
      IF(IPRINT.GE.5) THEN
        DO I = 1,NFMAT
          CALL HEADER('Symmetrized skeleton FOCK',-1)
          WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I,' of',NFMAT
          CALL PRQMAT(FAO(1,I),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
     &                NZC1,IQMULT,LUPRI)
        ENDDO
      ENDIF
      CALL QEXIT('DIRAOF_1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck conscf */
      SUBROUTINE CONSCF(FSO,DSO,NFMAT,
     &                  IFCKOP,ISYMOP,IHRMOP,INTFLG,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for conventional SCF
C
C     Written by T.Saue July 2 1995
C     Sep 2 1996 tsaue - Allocate buffer IND for unpacking of indices
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C Used from COMMON blocks:
C   DHFIOU: LUINT
C   INDPCK: NIBUF
C
#include "dcbind.h"
#include "dcbgen.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION FSO(*),DSO(*),
     &          IFCKOP(NFMAT),ISYMOP(NFMAT),IHRMOP(NFMAT),
     &          WORK(LWORK)
C
#include "memint.h"
C
C     Memory allocation
C
      NDMAT = NZ*NFMAT
      CALL MEMGET2('INTE','IFC' ,KIFC  ,NDMAT        ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IRD' ,KIRD  ,NDMAT        ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','BUF' ,KBUF  ,N2BBASX      ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IBUF',KIBUF ,N2BBASX*NIBUF,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IND' ,KIND  ,N2BBASX*4    ,WORK,KFREE,LFREE)
      CALL CONSC1(FSO,DSO,NFMAT,
     &            WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &            INTFLG,WORK(KBUF),WORK(KIBUF),
     &            WORK(KIND),IPRINT,WORK(KFREE),LFREE)
C
C     Memory deallocation
C
      CALL MEMREL('CONSCF',WORK,KWORK,KWORK,KFREE,LFREE)
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck CONSC1 */
      SUBROUTINE CONSC1(FSO,DSO,NFMAT,
     &                  IFCTYP,IREPDM,ISYMOP,IHRMOP,IFCKOP,
     &                  INTFLG,BUF,IBUF,IND,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Conventional SCF
C
C     Written by T.Saue Feb 1997
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbind.h"
#include "dcbgen.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      LOGICAL LBIT
      DIMENSION FSO(N2BBASX,*),DSO(N2BBASX,*),
     &          IFCTYP(NZ,NFMAT),IREPDM(NZ,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),
     &          BUF(*),IBUF(*),IND(*),
     &          WORK(LWORK)
C
      LUINTA = 14
      LUINTB = 15
      NDMAT  = NZ*NFMAT
      CALL SETFCK(IFCTYP,IREPDM,NFMAT,NZ,ISYMOP,IHRMOP,IFCKOP,IPRINT)
C
C     LL-integrals
C     ============

      IF(LBIT(INTFLG,1).AND.(.NOT.LBIT(IDFLAG,1))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: LL-integrals:'
        CALL TIMER('START ',TIMSTR,TIMEND)
        NNDIM = NTBAS(1)*(NTBAS(1)-1)/2
        CALL OPNFIL(LUINTA,'DFLLSA','OLD','LL2SMT')
        CALL OPNFIL(LUINTB,'DFLLSB','OLD','LL2SMT')
        CALL XX2MAT(LUINTA,LUINTB,NTBAS(1),NNDIM,
     &       FSO,DSO,NDMAT,IFCTYP,IREPDM,
     &       BUF,IBUF,IND,IPRINT)
        CLOSE(LUINTA,STATUS='KEEP')
        CLOSE(LUINTB,STATUS='KEEP')
        CALL TIMER('CONSCF(LL)',TIMSTR,TIMEND)
      ENDIF
C
C     SL-integrals
C     ============
C
      IF(LBIT(INTFLG,2).AND.(.NOT.LBIT(IDFLAG,2))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: SL-integrals:'
        CALL TIMER('START ',TIMSTR,TIMEND)
        CALL OPNFIL(LUINTA,'DFSLSA','OLD','SL2SMT')
        CALL OPNFIL(LUINTB,'DFSLSB','OLD','SL2SMT')
        CALL XY2MAT(LUINTA,LUINTB,FSO,DSO,NDMAT,IFCTYP,IREPDM,
     &              BUF,IBUF,IND,IPRINT)
        CLOSE(LUINTA,STATUS='KEEP')
        CLOSE(LUINTB,STATUS='KEEP')
        CALL TIMER('CONSCF(SL)',TIMSTR,TIMEND)
      ENDIF
C
C     SS-integrals
C     ============
C
      IF(LBIT(INTFLG,3).AND.(.NOT.LBIT(IDFLAG,3))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: SS-integrals:'
        CALL TIMER('START ',TIMSTR,TIMEND)
        NNDIM = NTBAS(2)*(NTBAS(2)-1)/2
        CALL OPNFIL(LUINTA,'DFSSSA','OLD','SS2TMT')
        CALL OPNFIL(LUINTB,'DFSSSB','OLD','SS2SMT')
        CALL XX2MAT(LUINTA,LUINTB,NTBAS(2),NNDIM,
     &      FSO,DSO,NDMAT,IFCTYP,IREPDM,
     &      BUF,IBUF,IND,IPRINT)
        CLOSE(LUINTA,STATUS='KEEP')
        CLOSE(LUINTB,STATUS='KEEP')
        CALL TIMER('CONSCF(SS)',TIMSTR,TIMEND)
      ENDIF
C
C     Gaunt-integrals
C     ===============
C     to be programmed
C
C     aspg, 2006-05-14: added if for gaunt integrals, at present
C     just makes the program stop if conventional scf is asked
C     together with gaunt integrals 
C
      IF(LBIT(INTFLG,4).AND.(.NOT.LBIT(IDFLAG,4))) THEN
        WRITE(LUPRI,'(A)') 'CONSCF: Gaunt-integrals not yet implemented'
        CALL QUIT('*** ERROR in CONSC1 ***') 
      ENDIF
C
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C*/Deck setfck */
      SUBROUTINE SETFCK(IFCTYP,IREPDM,NFMAT,MZ,ISYMOP,IHRMOP,IFCKOP,
     &                  IPRINT)
C*****************************************************************************
C
C     This routine will assign IFCTYP, indicating Fock matrices
C     to construct based on information about
C       hermiticiy of matrices     IHRMOP 0 - no symmetry
C                                         1 - Hermitian
C                                        -1 - anti-Hermitian
C       point group symmetry       ISYMOP 
C       Fock matrix type           IFCKOP 1, 2, or 3 - see IPFCK below
C
C     IFCTYP = XY
C       X - symmetry about the diagonal  - based on IPHRM
C         0 - G: (General) no symmetry
C         1 - S: Symmetric
C         2 - A: Anti-symmetric
C         3 - M: (Mixed) some blocks symmetric, other blocks anti-symmetric
C       Y - two-electron contributions   - based on IPFCK
C         1 - Coulomb
C         2 - Exchange
C         3 - Coulomb+Exchange
C
CMI    IREPDM is set up as well.
C
C     For the moment we will assume closed-shell Fock-matrices
C
C     Written by T.Saue Sep 4 1996
C
C*****************************************************************************
      use quaternion_algebra
#ifdef MOD_SRDFT
      use srdft_cfg ! for MCSCF-srDFT calculations
#endif
#include "implicit.h"
#include "priunit.h"    
      PARAMETER (D0 = 0.0D0)
C
#include "pgroup.h"
#include "dgroup.h"
#include "dcbham.h"
C     External variables
      DIMENSION IFCTYP(MZ,NFMAT),IREPDM(MZ,NFMAT),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT)
C     Local variables
      SAVE IPFCK,HTYP
c     DIMENSION IPFCK(4,3),IRFCK(4),IQFCK(4)
      DIMENSION IPFCK(4,0:3),IRFCK(4),IQFCK(4)
      CHARACTER HTYP(0:3)*1
      DATA HTYP /'G','S','A','M'/
c     DATA IPFCK/3,2,2,2, 1,0,0,0, 2,2,2,2/
      DATA IPFCK/0,0,0,0, 3,2,2,2, 1,0,0,0, 2,2,2,2/
#include "ibtfun.h"
C
      IF (IPRINT.GE.2) THEN
        CALL HEADER('Output from SETFCK',-1)
        WRITE(LUPRI,'(A,F10.4)') 
     &      '* Exchange factor: ',HFXFAC,
     &      '* Mu factor      : ',HFXMU,
     &      '* Attenuation    : ',HFXATT
        IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(3X,A,I2,A,I2)')
     &   '* Entering parameters:  NFMAT=',NFMAT,' MZ= ',MZ
         WRITE(LUPRI,*) 'ISYMOP(1:NFMAT) >',(ISYMOP(I),I=1,NFMAT)
         WRITE(LUPRI,*) 'IHRMOP(1:NFMAT) >',(IHRMOP(I),I=1,NFMAT)
         WRITE(LUPRI,*) 'IFCKOP(1:NFMAT) >',(IFCKOP(I),I=1,NFMAT)
        ENDIF
      ENDIF
      CALL IZERO(IFCTYP,MZ*NFMAT)
      DO I = 1,NFMAT
        IREP = ISYMOP(I)-1
C       IREP is always 0 (totally symmetric) for Hartre-Fock
C       (defined in SETDHF(dirac/dirset.F) and GETFC1(prp/pamset.F))
C       and it is the irrep of the perturbation operator in response
C       (defined in SETGMT(prp/pamrvc.F)). /hjaaj Apr 2002
        IHRM = IHRMOP(I)
        IPAR = JBTOF(IREP,1)
        IFCK = IFCKOP(I)
C
C       First determine contributions
C
        DO IZ = 1,NZC1
          IRFCK(IZ) = IRQMAT(IZ,IREP)
          IY        = IPFCK(IZ,IFCK)
          IC        = MOD(IY,2)
          IE        = IY - IC
          IH        = IHQMAT(IZ,IHRM)
C         An anti-symmetric matrix gives no Coulomb contributions
          IF (IH.EQ.2) IC = 0
C         No exchange contributions if HFXFAC = 0
          IF (HFXFAC.EQ.D0) IE = 0
C
C         No Coulomb contributions if ... 
C
C         ...   HFXMU.NE.0 and if this is NOT an MCSCF-srDFT calculation
C
C         Indeed, when building the long-range active Fock operator, 
C         we need both coulomb and exchange long-range integrals. 
C
          if (hfxmu /= 0.0d0) then
#ifdef MOD_SRDFT
             if (.not. srdft_cfg_lrcoulomb_int) then
#endif
                ic = 0
#ifdef MOD_SRDFT
             end if
#endif
          end if
          IY   = IC + IE
          IF (MZ.LT.4) THEN
C           ... we are using .SOFOCK or non-direct (CONSCF)
            IQFCK(IZ)     = IQMULT(1,JQBAS(IRFCK(IZ),IPAR),IZ)
            IPQ           = IQTOPQ(IQFCK(IZ),IREP)
C
            IYOLD = MOD(IFCTYP(IPQ,I),10)
            IXOLD = IFCTYP(IPQ,I)/10
            IFCTYP(IPQ,I) = 10*IBTOR(IXOLD,IH) + IBTOR(IYOLD,IY)
C           ... hjaaj apr 2002: with IBTOR we are sure we include
C           all the requested types. If some blocks are symmetric
C           and some are antisymmetric we end with the new type
C           of IX = 3 (IFCTYP = 10*IX + IY), which is in a way a
C           special case of IX = 0 (unsymmetric). /hjaaj Apr 2002
C
            IREPDM(IPQ,I) = IREP
C           ... MZ.lt.4 means that SOFOCK has been requested, and
C           the density matrix will be of mixed symmetry. We
C           use IREPDM to tell FCKOUT(abacus/her2out.F) which
C           irrep gives symmetric matrix elements and thus
C           Coulomb contributions. /hjaaj Apr 2002
C
          ELSE
C           ... MZ.eq.4: we are using DIRAOF (not .SOFOCK)
            IQFCK(IZ)     = IZ
            IFCTYP(IZ,I)  = 10*IH+IY
            IREPDM(IZ,I) =  IRFCK(IZ)
          ENDIF
        ENDDO
C
C     Print section
C
#if defined (DEBUG_DIRAC_FCKBUILD_DF)
        IF(IPRINT.GE.0) THEN
#else
        IF(IPRINT.GE.2) THEN
#endif
          WRITE(LUPRI,'(A,I5,3X,A,3X,A,A3)') 
     &      '* Fock matrix no.',I,' - ',
     &      'Operator irrep: ',REP(IREP)
          WRITE(LUPRI,'(A,4(2X,A3))')
     &      'Matrix symmetry:',(REP(IRFCK(IZ)),IZ=1,NZC1)
          IF(MZ.LT.4) THEN
            WRITE(LUPRI,'(A,4(3X,A1,1X))')
     &        'Packing        :',(QUNIT(IQFCK(IZ)),IZ=1,NZC1)
            WRITE(LUPRI,'(A,4(3X,A1,1X))')
     &        'Eff.   packing :',
     &        (QUNIT(IQTOPQ(IQFCK(IZ),IREP)),IZ=1,NZC1)
          ENDIF
            WRITE(LUPRI,'(A)') 
     &      '    Mat.  IFC  SYM  Coulomb Exchange Irrep'
          DO IZ = 1,MZ
            IY = MOD(IFCTYP(IZ,I),10)
            IH = (IFCTYP(IZ,I)-IY)/10
            IC = MOD(IY,2)
            IE = IY - IC
            IQ = IZ
            IF(MZ.LT.4) IQ = IPQTOQ(IZ,IREP)
            WRITE(LUPRI,'(4X,I1,3A1,I5,3X,A1,6X,I1,7X,I1,6X,A3)')
     &        IZ,'(',QUNIT(IQ),')',
     &        IFCTYP(IZ,I),HTYP(IH),IC,IE/2,REP(IREPDM(IZ,I))
          ENDDO
          CALL FLSHFO(LUPRI)
        ENDIF
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirsof */
      SUBROUTINE DIRSOF(FMAT,DSO,NFMAT,ISYMOP,IHRMOP,IFCKOP,
     &                  NPOS,INTFLG,IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations in SO-basis
C
C     Written by T.Saue August 28 1996
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER ( D0 = 0.0D0 )
C
#include "maxorb.h"
#include "aovec.h"
      LOGICAL LBIT
      DIMENSION FMAT(N2BBASX,NZ,NFMAT),DSO(N2BBASX,NZ,NFMAT),
     &          IFCKOP(NFMAT),NPOS(*),
     &          IHRMOP(NFMAT),ISYMOP(NFMAT),WORK(LWORK)
#include "dgroup.h"
#include "cbihr2.h"
#include "dcbbas.h"
#include "blocks.h"
C
      CALL QENTER('DIRSOF')
#include "memint.h"
C
C     Memory allocation
C
      NDMAT  = NFMAT*NZ
      N2GAB  = NSYMBL*NSYMBL
      NDMRAO = N2GAB*NDMAT
      NDMRSO = MAXSHL*MAXSHL*NDMAT
      CALL MEMGET2('INTE','IFC'  ,KIFC  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET2('INTE','IRD'  ,KIRD  ,NDMAT  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','GAB'  ,KGAB  ,N2GAB  ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DMRAO',KDMRAO,NDMRAO ,WORK,KFREE,LFREE)
      CALL MEMGET2('REAL','DMRSO',KDMRSO,NDMRSO ,WORK,KFREE,LFREE)
C
C     Get G(ab) integrals from AOPROPER
C     =================================
C
      IF (SCRFCK .GT. D0) THEN
         IJOB   = 0
         ITYPE  = 0
         MAXDIF = 0
         IGTYP  = 1
C        IGTYP 1 gives LL + SS type integrals
         IF(LBIT(INTFLG,2).OR.LBIT(INTFLG,3)) IGTYP = 1
C        IGTYP 2 will add also the SL (Gaunt) type integrals
         IF(LBIT(INTFLG,4)) IGTYP = 2
         CALL GETGAB(IJOB,ITYPE,IGTYP,MAXDIF,
     &        IPRINT,WORK(KGAB),WORK(KFREE),LFREE)
C
      END IF
C
      CALL DIRSO1(FMAT,DSO,NFMAT,
     &            WORK(KIFC),WORK(KIRD),ISYMOP,IHRMOP,IFCKOP,
     &            WORK(KDMRAO),WORK(KDMRSO),WORK(KGAB),
     &            NPOS,INTFLG,IPRINT,WORK(KFREE),LFREE)
      CALL MEMREL('DIRSOF.DIRSO1',WORK,KWORK,KWORK,KFREE,LFREE)
C
      CALL QEXIT('DIRSOF')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck dirso1 */
      SUBROUTINE DIRSO1(FSO,DSO,NFMAT,
     &                  IFCTYP,IREPDM,ISYMOP,IHRMOP,IFCKOP,
     &                  DMRAO,DMRSO,GABRAO,NPOS,INTFLG,
     &                  IPRINT,WORK,LWORK)
C*****************************************************************************
C
C     Driver for direct SCF calculations in SO-basis
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D2 = 2.0D0,D4 = 4.0D0,D8 = 8.0D0,D0=0.0D0)
#include "dummy.h"
#include "maxaqn.h"
C
#include "maxorb.h"
#include "aovec.h"
#include "cbihr2.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcbham.h"
#include "ccom.h"
#include "blocks.h"
#if defined (VAR_MPI)
#include "infpar.h"
#endif
      LOGICAl LBIT
      CHARACTER SECTID*12, INTTID*12, TOTTID*12, CPUTID*12, WALLTID*12,
     &          SOTEXT*9
      DIMENSION FSO(N2BBASX,*),DSO(N2BBASX,*),
     &          IFCTYP(*),IREPDM(*),
     &          ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),
     &          DMRAO(*),DMRSO(*),GABRAO(*),
     &          DNTSKP(8),NPOS(*),WORK(LWORK)
C
      CALL QENTER('DIRSO1')
#include "memint.h"
C
C     Initialization
C
      intflg_local = intflg
      if(twocomp .or. twocompbss) intflg_local = 1
      NDMAT = NFMAT*NZ
      CALL DZERO(DNTSKP,8)
      CALL SETFCK(IFCTYP,IREPDM,NFMAT,NZ,ISYMOP,IHRMOP,IFCKOP,IPRINT)
C
C     Make reduced density matrix for screening
C     =========================================
      IF(SCRFCK.GT.D0) THEN
C       in SO-basis
        CALL MKDRSO(DSO,DMRSO,NDMAT,WORK,LWORK,IPRINT)
C       in AO-basis
        CALL MKDRQO(DSO,DMRAO,NFMAT,ISYMOP,WORK,LWORK,IPRINT)
      ENDIF
C
C     ****************************
C     **** Process integrals *****
C     ****************************
C
      IRNTYP = 9
      MAXDIS = 1
      MAXDIF = 0
      JATOM  = 0
      KPOS   = 1
      THRSBUF = THRS
      IF (PARCAL) THEN
         CALL SCRSTA('MPIHeader',DUM,DUM)
      ELSE
         CALL SCRSTA('Header',DUM,DUM)
      END IF
C
C....> Start looping over the four possible integral classes
C
      DO I2TYP = 1, 4
C
C....> Set the parameters that are special for each class
C
      IF (I2TYP.EQ.1) THEN
         MTOTTK = NLRGBL*(NLRGBL+1)/2
         THRS   = THRSBUF
         SOTEXT = 'SOfock:LL'
      ELSEIF (I2TYP.EQ.2) THEN
         MTOTTK = NLRGBL*(NLRGBL+1)/2
         THRS   = THRSBUF*THRFAC(1)
         SOTEXT = 'SOfock:SL'
      ELSEIF (I2TYP.EQ.3) THEN
         MTOTTK = NSMLBL*(NSMLBL+1)/2
         THRS   = THRSBUF*THRFAC(2)
         SOTEXT = 'SOfock:SS'
      ELSEIF (I2TYP.EQ.4) THEN
         MTOTTK = NSMLBL*NLRGBL
         THRS   = THRSBUF*THRFAC(1)
         SOTEXT = 'SOfock:GT'
Cluuk: scale real part by factor 3 and eliminate Coulomb
Cluuk: this code does only work for C1 and Ci, use AOFOCK
        IF (LBIT(intflg_local,4).AND.NZ.LT.4) 
     &     call quit ('Gaunt not implemented for SOFOCK')
         factor = -3.d0
         CALL DSCAL(N2BBASX,factor,DSO,1)
         IFCTYP(1)=IFCTYP(1)+4
      ENDIF
C
C....> Check whether the class is active
C
      IF(LBIT(intflg_local,I2TYP).AND.LBIT(IDFLAG,I2TYP)) THEN
#if defined (VAR_MPI)
        IF (PARCAL) THEN
C
C          Get hold of the slaves ( ITASK = 1 for Fock matrices )
C
           CALL DIRAC_PARCTL( HERFCK_PAR )
C
           CALL GETTIM(CPU1,WALL1)
           CALL HER_PARDRV(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,
     &          IRNTYP,MAXDIF,JATOM,.TRUE.,.TRUE.,.FALSE.,
     &          TKTIME,.FALSE.,FIRST,NPOS(KPOS),MTOTTK,I2TYP,
     &          ICEDIF,SCRFCK,GABRAO,DMRAO,DMRSO,DNTSKP)
C
           CALL GETTIM(CPU2,WALL2)
           CPU    = CPU2 - CPU1
           WALL   = WALL2 - WALL1
           TOTWAL = TOTWAL + WALL
           CALL SCRSTA(SOTEXT,DNTSKP,WALL)
           IF (IPRINT .GT. 1) THEN
              TOTTID = SECTID(TOTWAL)
              WRITE(LUPRI,'(A,A12)')
     &             '>>>> Total wall time used in HER_PARDRV so far   :',
     &             TOTTID
              IF (IPRINT .GT. 5) THEN
                 CPUTID = SECTID(CPU)
                 WALLTID = SECTID(WALL)
                 WRITE(LUPRI,'(A,A12,A1,A12)')
     &        '>>>> CPU/wall  time used in HER_PARDRV last iteration :',
     &         CPUTID,'/',WALLTID
              END IF
           END IF
        ELSE
#endif
           CPU = SECOND()
           CALL TWOINT(WORK,LWORK,FSO,DSO,NDMAT,IREPDM,IFCTYP,DUM,
     &          IDUM,IDUM,IRNTYP,MAXDIF,JATOM,.TRUE.,
     &          .TRUE.,.FALSE.,TKTIME,IPRTWO,IPRNTA,IPRNTB,IPRNTC,
     &          IPRNTD,RTNTWO,IDUM,I2TYP,ICEDIF,SCRFCK,
     &          GABRAO,DMRAO,DMRSO,DNTSKP,.TRUE.,.false.,IDUM,DUM)
           CPU = SECOND() - CPU
           CALL SCRSTA(SOTEXT,DNTSKP,CPU)
#if defined (VAR_MPI)
        ENDIF
#endif
      ENDIF
cluuk correct for scaling
      IF (I2TYP.EQ.4) THEN
         factor = 1.d0/factor
         CALL DSCAL(N2BBASX,factor,DSO,1)
      ENDIF
C
C---> Update the offset
C
      KPOS = KPOS + MTOTTK
C
C.... End of loop over integral classes
C
!      WRITE(LUPRI,'(A,I5,A,I5)') '*** Matrix no.',I2typ,' of',4
!      CALL PRQMAT(FSO,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),
!    &             NZ,IPQTOQ(1,0),LUPRI)
         CALL FLSHFO(LUPRI)
      ENDDO
C
C     Scale matrices
C
      IF (HFXMU.EQ.D0) THEN
        DO I = 1,NDMAT
          IY = MOD(IFCTYP(I),10)
          IX = (IFCTYP(I)-IY)/10
          IF(IX.EQ.0) THEN
            CALL DSCAL(N2BBASX,D2,FSO(1,I),1)
          ELSE
            CALL DSCAL(N2BBASX,D4,FSO(1,I),1)
          ENDIF
        ENDDO
      ENDIF
      CALL QEXIT('DIRSO1')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Mkdrqo */
      SUBROUTINE MKDRQO(DQO,DMRQO,NFMAT,ISYMOP,
     &                  WORK,LWORK,IPRINT)
C*****************************************************************************
C
C     Make reduced density matrices in AO-basis starting from
C     density matrices in QO-basis
C
C     Written by T.Saue Oct 11 1996
C
C*****************************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "aovec.h"
C
#include "dgroup.h"
#include "dcbbas.h"
#include "blocks.h"
      DIMENSION DQO(N2BBASX,NZ,NFMAT),
     &          DMRQO(NSYMBL,NSYMBL,NZ,NFMAT),
     &          ISYMOP(NFMAT),WORK(LWORK)
C
#include "ibtfun.h"
      CALL QENTER('MKDRQO')
#include "memint.h"
C
C     Initialize
C
      NDMAT = NFMAT*NZ
      N2DRAO = NSYMBL*NSYMBL*NDMAT
      CALL DZERO(DMRQO,N2DRAO)
C
C     Make index array from AOs to blocks
C
      CALL MEMGET2('INTE','INRAO',KINRAO,NTBAS(0),WORK,KFREE,LFREE)
      CALL PTRDAO(WORK(KINRAO))
C
C     Make reduced AO matrix
C
      IF(NBSYM.GT.1) THEN
        CALL MEMGET2('REAL','DBF',KDBF,N2BBASX,WORK,KFREE,LFREE)
        DO I = 1,NFMAT
          IREP = ISYMOP(I)-1
          IPAR = JBTOF(IREP,1)
          DO IZ = 1,NZC1
            IREPDM = IRQMAT(IZ,IREP)
            IQ     = IQMULT(1,JQBAS(IREPDM,IPAR),IZ)
            IQP    = IQTOPQ(IQ,IREP)
            IF(IQP.GT.NZ) CALL QUIT('MKDRQO: IQP.GT.NZ !')
            IF(IQP.LE.0 ) CALL QUIT('MKDRQO: IQP.LE.0 !')
            CALL DTSOAO(DQO(1,IQP,I),WORK(KDBF),
     &                 NTBAS(0),IREPDM,IPRINT)
            CALL GATMAT(0,NTBAS(0),NTBAS(0),WORK(KDBF),
     &              DMRQO(1,1,IQP,I),WORK(KINRAO),NSYMBL)
          ENDDO
        ENDDO
      ELSE
        DO I = 1,NFMAT
          DO IZ = 1,NZ
            CALL GATMAT(0,NTBAS(0),NTBAS(0),DQO(1,IZ,I),
     &              DMRQO(1,1,IZ,I),WORK(KINRAO),NSYMBL)
          ENDDO
        ENDDO
      ENDIF
C
C     Print section
C     =============
C
      IF(IPRINT.GE.4) CALL PRDRAO(DMRQO,NDMAT,WORK(KINRAO))
C
      CALL MEMREL('MKDRQO',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('MKDRQO')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Twodrv */
      SUBROUTINE TWODRV(INTTYP,I2TYP,WORK,LWORK,IPRINT)
C***********************************************************************
C
C     General driver for the calculation of two-electron integrals
C     by either the ERI or the HERMIT code
C
C     INTTYP:
C        0 - Undifferentiated integrals written on file
C
C     Written by T. Saue, Nov 2004
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "aovec.h"
#include "cbihr2.h"
#include "ccom.h"
#include "dgroup.h"
C
      LOGICAL NODV,NOPV,NOCONT,RELCAL
      DIMENSION WORK(LWORK)
      CHARACTER CINT(3)*2
      DATA CINT/'LL','SL','SS'/
#include "memint.h"
C.......Undifferentiated integrals written on file
        IF(INTTYP.EQ.0) THEN
          NDMAT  = 1
          THRBUF = THRS
          IF(I2TYP.EQ.2) THRS = THRS*THRFAC(1)
          IF(I2TYP.EQ.3) THRS = THRS*THRFAC(2)
          IF(IPRINT.GE.2) THEN
            WRITE(LUPRI,'(1P,A2,A,D9.2)') 
     &         CINT(I2TYP),'-int - threshhold:',THRS
          ENDIF
        ENDIF
        MAXDIS = 1
        MAXDIF = 0
        JATOM  = 0
        NODV   = .TRUE.
        NOPV   = .TRUE.
        NOCONT = .FALSE.
        RELCAL = .TRUE.
        ISHLA  = 0
        CALL TWOINT(WORK,LWORK,FMAT,DMAT,NDMAT,IREPDM,IFCTYP,GMAT,
     &               INDX,INDXAB,INTTYP,MAXDIF,JATOM,NODV,
     &               NOPV,NOCONT,TKTIME,IPRINT,IPRNTA,IPRNTB,IPRNTC,
     &               IPRNTD,RTNTWO,ISHLA,I2TYP,ICEDIF,SCREEN,
     &               GABRAO,DMRAO,DMRSO,DINTSKP,RELCAL,.false.,
     &               INDEX,TIMES)
        IF(INTTYP.EQ.0) THEN
          THRS = THRBUF 
        ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Fockbuild */
      SUBROUTINE FOCKBUILD(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &                     IFCKOP,NPOS,INTFLG,
     &                     IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Fock matrix builder
C
C     Written by T. Saue Sep 5 2006
C
C***********************************************************************
#include "implicit.h"
      INTEGER   NFMAT, LWORK, INTFLG
      REAL(8)   FMAT(*),DMAT(*),WORK(LWORK)
      INTEGER   ISYMOP(NFMAT),IHRMOP(NFMAT),IFCKOP(NFMAT),NPOS(*)

#include "priunit.h"
#include "maxaqn.h"
#include "dummy.h"
      PARAMETER(D0=0.0D0, D1=1.0D0)
C
#include "dcbgen.h"
#include "cbihr2.h"
#include "dcbham.h"
C
      real(8) :: TIMSTR(2), TIMEND(2)
C
C     *************************************************
C     *****  D i r e c t  c o n s t r u c t i o n *****
C     *****  o f   F o c k   m a t r i c e s      *****
C     *************************************************
C
C     This module has to come before the conventional SCF-module;
C     otherwise FOCK must be copied into buffer; this stems from
C     the symmetrization procedure
C
C aspg, 2006-05-15: added test so that whenever gaunt integrals
C       and at least one other type of integral is asked for in
C       the direct calculation it assumes all other types are also
C       calculated on the fly. otherwise conventional scf is
C       performed. 
C
      IF(IDFLAG.NE.0.AND.IDFLAG.GT.8) THEN
        IF(SOFOCK) THEN
C.........construction of Fock matrix using ERIs in SO-basis
          CALL TIMER2('START ',TIMSTR,TIMEND)
          CALL DIRSOF(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &                IFCKOP,NPOS,INTFLG,
     &                IPRINT,WORK,LWORK)
          CALL TIMER2('SO Fock',TIMSTR,TIMEND)
        ELSE
C.........construction of Fock matrix using ERIs in AO-basis:
C.........skeleton Fock approach
          CALL TIMER2('START ',TIMSTR,TIMEND)

        !DIM = N2BBASX*NZC1*NFMAT
C       Conventional (4-center integrals) AO Coulomb + exchange
        CALL DIRAOF(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
     &              IFCKOP,NPOS,INTFLG,
     &              IPRINT,WORK,LWORK)
!          CALL FASTJK(FMAT,DMAT,NFMAT,ISYMOP,IHRMOP,
!     &                IFCKOP,NPOS,INTFLG,
!     &                IPRINT,WORK,LWORK)
!         CALL TIMER2('fastjk',TIMSTR,TIMEND)
          CALL TIMER2('AO Fock',TIMSTR,TIMEND)
        ENDIF
      ELSE
C
C     ******************************************
C     *****  D i s k - b a s e d    S C F  *****
C     ******************************************
C
        IF(HFXFAC.NE.D1 .AND. HFXFAC.NE.D0) 
     &       CALL QUIT('HFXFAC not implemented in CONSCF !')
C       hjaaj July 2004: HFXFAC .eq. D0 is defined through SETFCK
        CALL CONSCF(FMAT,DMAT,NFMAT,IFCKOP,ISYMOP,IHRMOP,INTFLG,
     &              IPRINT,WORK,LWORK)
      ENDIF
C
      RETURN
      END

Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C
C aspg, 20080411
C
C  Parallel calculations require an integer array describing the
C  slave tasks as argument to TWOFCK.
C
C  below are some utility routines to handle the transfer of such
C  information. they are pretty much what was in the STEX_GET_NPOS
C  subroutine by ulfek but it seems better to decouple the setting
C  of the FIRST* variables (used to indicate that this is the first
C  time the integral code is being called for a given class of
C  integrals - LL,SL,SS,GT) and the settinf of the dimension of the
C  integer array (heferred to in the code as NPOS)
C
C  i have made the names a bit longer than usual, but at least i hope
C it makes it more obvious what is happening.
C

C
C  SaveTaskDistribFlags:
C
C  saves the content of the FIRST* flags into the array
C  saveflags for later use
C
      SUBROUTINE SaveTaskDistribFlags(saveflags)
#include "implicit.h"
#include "dcbfir.h"
      LOGICAL saveflags(4)

      saveflags(1) = FIRST1
      saveflags(2) = FIRST2
      saveflags(3) = FIRST3
      saveflags(4) = FIRST4
      END
C
C  SetTaskDistribFlags:
C
C  sets the FIRST* variables to make TWOFCK fill the array
C  with new content.
C
      SUBROUTINE SetTaskDistribFlags(flags)
#include "implicit.h"
#include "dcbfir.h"
      LOGICAL FLAGS(4)

      FIRST1 = flags(1)
      FIRST2 = flags(2)
      FIRST3 = flags(3)
      FIRST4 = flags(4)
      END
C
C  SetIntTaskArrayDimension: sets TskArrDim
C
C  TskArrDim will give the dimension of the variable "NPOS" in
C  the case of parallel calculations. for serial runs, it will
C  have value zero
C
      SUBROUTINE SetIntTaskArrayDimension(TskArrDim,isParallel)
#include "implicit.h"
#include "priunit.h"
#include "aovec.h"
#include "maxorb.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcbbas.h"
#include "cbihr2.h"
#include "blocks.h"
      logical isParallel
      integer TskArrDim

      TskArrDim = 0
      IF (isParallel) THEN
C  We will always allocate the maximum size, even if some
C  integral classes are inactive (keep the code simple)
         TskArrDim = TskArrDim + NLRGBL*(NLRGBL+1)/2
         TskArrDim = TskArrDim + NLRGBL*(NLRGBL+1)/2
         TskArrDim = TskArrDim + NSMLBL*(NSMLBL+1)/2
         TskArrDim = TskArrDim + NSMLBL*NLRGBL
      ENDIF
      END
C
Ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

