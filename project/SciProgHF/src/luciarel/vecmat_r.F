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

      SUBROUTINE EXPCIV_REL(ISM,ISPCIN,LUIN,
     &                  ISPCUT,LUUT,LBLK,
     &                  LUSCR,NROOT,ICOPY,IDC,NTESTG)
*
* Expand CI vector in CI space ISPCIN to CI vector in ISPCUT
* Input vector is supposed to be on LUIN
* Output vector will be placed on unit LUUT
*. If ICOPY .ne. 0 the output vectors will be copied to LUIN
*
* Storage form is defined by ICISTR
*
* Jeppe Olsen, February 1994
* GAS version August 1995
*
      use symmetry_setup_krci
      IMPLICIT REAL*8(A-H,O-Z)
#include "ipoist8.inc"
#include "mxpdim.inc"
#include "cicisp.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "cgas.inc"
#include "gasstr.inc"

*
      call memmar(KDUM,IDUM,'MARK  ',IDUM,'EXPCIV')
*
      NTESTL = 1
      NTEST = MAX(NTESTG,NTESTL)
*
      IATP = 1
      IBTP = 2
*
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
*
*. Allowed combinations of strings types for input and output
*. spaces
*
      call memmar(KLIABI,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABI')
      call memmar(KLIABU,NOCTPA*NOCTPB,'ADDL  ',1,'KLIABU')
      CALL IAIBCM_GAS(LCMBSPC(ISPCIN),ICMBSPC(1,ISPCIN),IGSOCCX,
     &            NOCTPA,NOCTPB,
     &            ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,
     &            WORK(KLIABI),NTEST)
*
      CALL IAIBCM_GAS(LCMBSPC(ISPCUT),ICMBSPC(1,ISPCUT),IGSOCCX,
     &            NOCTPA,NOCTPB,
     &            ISPGPFTP(1,IOCTPA),ISPGPFTP(1,IOCTPB),
     &            NELFGP,MXPNGAS,NGAS,
     &            WORK(KLIABU),NTEST)
*
* type of each symmetry block ( full, lower diagonal, absent )
*
      call memmar(KLBLIN,NSMST,'ADDL  ',1,'KLBLIN')
      call memmar(KLBLUT,NSMST,'ADDL  ',1,'KLBLUT')
      CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KLBLIN))
      CALL ZBLTP_REL(ISMOST(1,ISM),NSMST,WORK(KLBLUT))

*
*. Number of dets etc per TTS block
*
      call memmar(KLTTSII,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
      call memmar(KLTTSNI,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLIN),NSMST,WORK(KLIABI),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSII),WORK(KLTTSNI),
     &            NCOMBI,0)
*
      call memmar(KLTTSIU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSI')
      call memmar(KLTTSNU,NSMST*NOCTPA*NOCTPB,'ADDL  ',1,'KLTTSN')
        CALL ZOOS(ISMOST(1,ISM),WORK(KLBLUT),NSMST,WORK(KLIABU),
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO2(IBTP)),
     &            NOCTPA,NOCTPB,IDC,WORK(KLTTSIU),WORK(KLTTSNU),
     &            NCOMBU,0)
*
*. Allocate memory for blocks of input and output space in
*  accordance with ICISTR
*
      IF(ICISTR.EQ.1) THEN
        LENGTHI = NCOMBI
        LENGTHU = NCOMBU
      ELSE IF (ICISTR .EQ. 2 ) THEN
        LENGTHI = MXSB
        LENGTHU = MXSB
      ELSE IF (ICISTR.EQ. 3 ) THEN
        LENGTHI = MXSOOB
        LENGTHU = MXSOOB
      END IF
*
      call memmar(KLBLI,LENGTHI,'ADDL  ',2,'KLBLI ')
      call memmar(KLBLU,LENGTHU,'ADDL  ',2,'KLBLU ')
*
*. and now : Let another subroutine complete the taks
*
      CALL REWINE(LUIN,-1)
      CALL REWINE(LUUT,-1)
*
*. Print for testing initial vectors out
*
C?    WRITE(6,*) ' Initial vectors in EXPCIV '
C?    WRITE(6,*) ' Initial vectors in EXPCIV '
C?    WRITE(6,*) ' Initial vectors in EXPCIV '
*
C?    DO IROOT = 1, NROOT
C?    WRITE(6,*) ' Root number ', IROOT
*
C?      IF(ICISTR.EQ.1) THEN
C?          CALL WRTMAT(WORK(KLBLI),1,NCOMBUT,1,NCOMBUT)
C?      ELSE
C?        CALL WRTVCD(WORK(KLBLI),LUIN,0,-1)
C?      END IF
C?    END DO
C?    CALL REWINE(LUIN,-1)
*
      write(6,*) ' ICISTR = ', ICISTR
      DO IROOT = 1, NROOT
*. Input vector should be first vector on file so
        IF(IROOT.EQ.1) THEN
          LLUIN = LUIN
        ELSE
          IF(ICISTR.EQ.1) THEN
            CALL REWINE(LUSCR,-1)
            CALL FRMDSC(WORK(KLBLI),NCOMBI,-1,LUIN,IMZERO,IAMPACK)
            CALL  TODSC(WORK(KLBLI),NCOMBI,-1,LUSCR)
            CALL REWINE(LUSCR,-1)
            LLUIN = LUSCR
          ELSE
*. With the elegance of an elephant
            CALL REWINE(LUSCR,-1)
            CALL REWINE(LUIN,-1)
            DO JROOT = 1, IROOT
              CALL REWINE(LUSCR,-1)
              CALL COPVCD_REL(LUIN,LUSCR,WORK(KLBLI),0,-1)
            END DO
            CALL REWINE(LUSCR,-1)
C?          write(6,*) ' Vector copied to LUSCR '
C?          CALL WRTVCD(WORK(KLBLI),LUSCR,1,-1)
C?          CALL REWINE(LUSCR,-1)
            LLUIN = LUSCR
          END IF
        END IF
*
        ITTSS_ORD = 2
        CALL EXPCIVS_REL(LLUIN,WORK(KLBLI),NCOMBI,WORK(KLTTSNI),
     &       WORK(KLTTSII),
     &       WORK(KLIABI),NOCTPA,NOCTPB,WORK(KLBLIN),
     &       LUUT,WORK(KLBLU),NCOMBU,WORK(KLTTSNU),
     &       WORK(KLTTSIU),WORK(KLIABU),
     &       WORK(KLBLUT),
     &       ICISTR,IDC,NSMST,
     &       LBLK,ITTSS_ORD)
*
      END DO
*
      IF(ICOPY.NE.0) THEN
*. Copy expanded vectors to LUIN
        CALL REWINE(LUIN,-1)
        CALL REWINE(LUUT,-1)
        DO IROOT = 1, NROOT
          IF(ICISTR.EQ.1) THEN
            CALL FRMDSC(WORK(KLBLU),NCOMBU,-1,LUUT,IMZERO,IAMPACK)
            CALL  TODSC(WORK(KLBLU),NCOMBU,-1,LUIN)
          ELSE
            CALL COPVCD_REL(LUUT,LUIN,WORK(KLBLU),0,-1)
          END IF
        END DO
      END IF
*
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
C?    WRITE(6,*) ' OUtput  vectors in EXPCIV '
*
C?    CALL REWINE(LUUT,-1)
C?    DO IROOT = 1, NROOT
C?    WRITE(6,*) ' Root number ', IROOT
*
C?      IF(ICISTR.EQ.1) THEN
C?          CALL WRTMAT(WORK(KLBLU),1,NCOMBUT,1,NCOMBUT)
C?      ELSE
C?        CALL WRTVCD(WORK(KLBLU),LUUT,0,-1)
C?      END IF
C?    END DO
      call memmar(KDUM,IDUM,'FLUSM ',IDUM,'EXPCV2')
*
      END
***********************************************************************

      SUBROUTINE EXPCIVS_REL(LUIN,VECIN,NCOMBIN,NTTSIN,ITTSIN,IABIN,
     &                       NOCTPA,NOCTPB,IBLTPIN,
     &                       LUUT,VECUT,NCOMBUT,NTTSUT,ITTSUT,IABUT,
     &                       IBLTPUT,
     &                       ICISTR,IDC,NSMST,LBLK,ITTSS_ORD)
*
* Obtain those part of vector in cispace UT ,
* that can be obtained from terms in cispace IN
*
* Input vector on LUIN, Output vector in LUUT
* Output vector is supposed on start of vector
*
* LUIN is assumed to be single vector file,
* so rewinding will place vector on start of vector
*
* Both files are assumed on start of vector
*
* Jeppe Olsen, February 1994
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "ipoist8.inc"
*. Input
      INTEGER IABIN(NOCTPA,NOCTPB),IABUT(NOCTPA,NOCTPB)
      INTEGER IBLTPIN(NSMST),IBLTPUT(NSMST)
      INTEGER NTTSIN(NOCTPA,NOCTPB,NSMST),ITTSIN(NOCTPA,NOCTPB,NSMST)
      INTEGER NTTSUT(NOCTPA,NOCTPB,NSMST),ITTSUT(NOCTPA,NOCTPB,NSMST)
      DIMENSION VECIN(*)
C     DIMENSION IATPUI(*),IBTPUI(*)
*. Output
      DIMENSION VECUT(*)
*
      WRITE(6,*) ' EXPCIVS in action'
      WRITE(6,*) '        Number of input  parameters',NCOMBIN
      WRITE(6,*) '        Number of output parameters',NCOMBUT
*
      IF(ICISTR.EQ.1) THEN
        CALL FRMDSC(VECIN,NCOMBIN,-1,LUIN,IMZERO,IAMPACK)
C       WRITE(6,*) ' Input vector '
C       CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
      END IF
*
*. Loop over TTS blocks of output vector
*
      IATPIN = 1
      IBTPIN = 0
      IASMIN = 1
*
      IATPUT = 1
      IBTPUT = 0
      IASMUT = 1
*
 1000 CONTINUE
*. Next output block
        CALL NXTBLK(IATPUT,IBTPUT,IASMUT,NOCTPA,NOCTPB,NSMST,
     &              IBLTPUT,IDC,NONEWUT,IABUT,ITTSS_ORD)
*. Corresponding input TTS block
        JATPIN = IATPUT
        JBTPIN = IBTPUT
        IF(IABIN(JATPIN,JBTPIN).EQ.0) THEN
          IZERO = 1
        ELSE
          IZERO = 0
        END IF
*
        NELMNT = NTTSUT(IATPUT,IBTPUT,IASMUT)
        IF(ICISTR.EQ.1) THEN
          IF(IZERO.EQ.0)
     &    IOFFIN = ITTSIN(JATPIN,JBTPIN,IASMUT)
          IOFFUT = ITTSUT(IATPUT,IBTPUT,IASMUT)
          IF(IZERO.EQ.0) THEN
            CALL COPVEC(VECIN(IOFFIN),VECUT(IOFFUT),NELMNT)
          ELSE
            ZERO = 0.0D0
            CALL SETVEC(VECUT(IOFFUT),ZERO,NELMNT)
          END IF
        END IF
*
        IF(ICISTR.NE.1.AND.NONEWUT.EQ.0) THEN
          IF(IZERO.EQ.1) THEN
C           ZERO = 0.0D0
C           CALL SETVEC(VECUT,ZERO,NELMNT)
            CALL ITODS(NELMNT,1,-1,LUUT)
            CALL ZERORC_REL(-1,LUUT)
C           CALL TODSC(VECUT,NELMNT,-1,LUUT)
          ELSE
*. Obtain input block
  999      CONTINUE
           CALL NXTBLK(IATPIN,IBTPIN,IASMIN,NOCTPA,NOCTPB,NSMST,
     &                 IBLTPIN,IDC,NONEWIN,IABIN,ITTSS_ORD)
           IF(NONEWIN.NE.0) THEN
             CALL REWINE(LUIN,-1)
             IATPIN = 1
             IBTPIN = 0
             IASMIN = 1
             GOTO 999
           ELSE IF (NONEWIN.EQ. 0 ) THEN
             CALL IFRMDS(LENGTH,1,-1,LUIN)
             CALL FRMDSC(VECIN,LENGTH,-1,LUIN,IMZERO,IAMPACK)
             IF(IATPIN.EQ.JATPIN.AND.IBTPIN.EQ.JBTPIN.AND.
     &          IASMIN.EQ.IASMUT) THEN
*. Correct block, save it
                  CALL ITODS(LENGTH,1,-1,LUUT)
                  IF(IMZERO.EQ.1) THEN
                    CALL ZERORC_REL(-1,LUUT)
                  ELSE
                    IF(IAMPACK.EQ.0) THEN
                      CALL TODSC(VECIN,LENGTH,-1,LUUT)
                    ELSE
                      CALL TODSCP(VECIN,LENGTH,-1,LUUT)
                    END IF
                  END IF
             ELSE
               GOTO 999
             END IF
           END IF
          END IF
        END IF
      IF(NONEWUT.EQ.0) GOTO 1000
*. End of loop over output blocks
      IF(ICISTR.EQ.1) THEN
        CALL TODSC(VECUT,NCOMBUT,-1,LUUT)
      ELSE
        CALL ITODS(-1,1,-1,LUUT)
      END IF
*
      NTEST = 0
      IF(NTEST.NE.0) THEN
        WRITE(6,*) ' EXPCIVS Speaking '
        WRITE(6,*) ' ================='
        WRITE(6,*)
        WRITE(6,*) ' ============ '
        WRITE(6,*) ' Input Vector '
        WRITE(6,*) ' ============ '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECIN,1,NCOMBIN,1,NCOMBIN)
        ELSE
          CALL WRTVCD(VECIN,LUIN,1,LBLK)
        END IF
        WRITE(6,*)
        WRITE(6,*) ' =============== '
        WRITE(6,*) ' Output Vector '
        WRITE(6,*) ' =============== '
        WRITE(6,*)
        IF(ICISTR.EQ.1) THEN
          CALL WRTMAT(VECUT,1,NCOMBUT,1,NCOMBUT)
        ELSE
          CALL WRTVCD(VECUT,LUUT,1,LBLK)
        END IF
      END IF
*
      END
***********************************************************************

      SUBROUTINE GSTTBLD(CTT,IATP,IASM,IBTP,IBSM,
     &                  NSASO,NSBSO,PSSIGN,IDC,
     &                  PLSIGN,LUC,SCR,NSMST,ISCALE,SCLFAC)
*
* obtain  determinant block (iatp iasm, ibtp ibsm )
* from vector packed in combination format according to IDC
*
*. If ISCALE = 1, the routine scales and returns the block
*  in determinant normalization, and SCLFAC = 1.0D0
*
* If ISCALE = 0, the routine does not perform any overall
* scaling, and a scale factor is returned in SCLFAC
*
* IF ISCALE = 0, zero blocks are not set explicitly to zero,
* instead  zero is returned in SCLFAC
*
* ISCALE, SCLFAC added May 97
*
* Simplified version working only for vectors on disc
*
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION CTT(*),NSASO(NSMST, *),NSBSO(NSMST, *)
      DIMENSION SCR(*)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN 
        write(6,*) ' GSTTBLD, IATP,IASM,IBTP,IBSM,ISCALE'
        write(6,*)            IATP,IASM,IBTP,IBSM,ISCALE     
        WRITE(6,*) ' LUC = ', LUC
      END IF
* =================
* Read in from disc
* =================
      CALL IFRMDS(LBL,1,-1,LUC)
      IF(NTEST.GE.100) write(6,*) ' LBL = ', LBL
      IF(ISCALE.EQ.1) THEN
        CALL FRMDSC(SCR,LBL,-1,LUC,IMZERO,IAMPACK)
      ELSE
        NO_ZEROING = 1
        CALL FRMDSC2(SCR,LBL,-1,LUC,IMZERO,IAMPACK,NO_ZEROING)
      END IF
*
      IF(IMZERO.EQ.1.AND.ISCALE.EQ.0) THEN
        SCLFAC = 0.0D0
      ELSE
        NAST = NSASO(IASM,IATP)
        NBST = NSBSO(IBSM,IBTP)
        IF(LBL.NE.0) THEN
          CALL SDCMRF(CTT,SCR,2,IATP,IBTP,IASM,IBSM,NAST,NBST,
     &         IDC,PSSIGN,PLSIGN,ISGVST,LDET,LCOMB,ISCALE,SCLFAC)
        ELSE
          SCLFAC = 0.0D0
        END IF
      END IF
*
C?    WRITE(6,*) ' ISCALE and SCLFAC on return in GSTTBL',
C?   &ISCALE,SCLFAC
*
      END
***********************************************************************

      SUBROUTINE GSTTBLD2(CTT,LUC)
*
*     obtain  determinant block from disc
*
      IMPLICIT REAL*8  (A-H,O-Z)
      DIMENSION CTT(*)
*
      NTEST = 000
*
      IF(NTEST.GE.100) THEN 
        WRITE(6,*) ' LUC = ', LUC
      END IF
*     =================
*     Read in from disc
*     =================
      CALL IFRMDS(LBL,1,-1,LUC)
      CALL FRMDSC(CTT,LBL,-1,LUC,IMZERO,IAMPACK)
*
      END
***********************************************************************

      FUNCTION IOFF_BLOCK_MAT(LROW,LCOL,IBLK)
*
* Offset of block IBLK in blocked matrix
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION LROW(*), LCOL(*)
*
      IOFF = 1
      DO JBLK = 1, IBLK-1
        IOFF = IOFF + LROW(JBLK)*LCOL(JBLK)
      END DO
*
      IOFF_BLOCK_MAT = IOFF
C?    WRITE(6,*) ' Offset for block ', IBLK, ' is ', IOFF
*
      END
***********************************************************************

      SUBROUTINE IOFF_SYMBLK_MAT(NSMST_in,NA,NB,ITOTSM,IOFF)
*
* Offset for symmetry blocked matrix. 
* Row-symmetry used as primary index
*
* Jeppe Olsen, Summer of 99
*               July 2000 : IRESTRICT added
* Modified for relativistic case, Timo Fleig, July 2001
*
      use symmetry_setup_krci
#include "implicit.inc"
*. Specific input 
      INTEGER NA(NSMST_in),NB(NSMST_in)
*. Output
      INTEGER IOFF(NSMST_in)
*
      IB = 1
      DO IASM = 1, NSMST_in
        IBSM       = IDBGMULT(INVELM(IASM),ITOTSM)
        IOFF(IASM) = IB
        IB         = IB + NA(IASM)*NB(IBSM)
      END DO
*
      END
***********************************************************************

      SUBROUTINE IOFF_SYMBLK_MAT_PAR(NSMST_in,NA,NB,ITOTSM,IOFF,
     &                               IACTLIST,IOFF_ACT,IVERS)
*
* Offset for symmetry blocked matrix. 
* Row-symmetry used as primary index
*
* Jeppe Olsen, Summer of 99
*               July 2000 : IRESTRICT added
* Modified for relativistic case, Timo Fleig, July 2001
* adaption of sequential routine, Stefan Knecht May 2007
*
* IVERS     = 0 : c-block routine
*           = 1 : s-block routine
      use symmetry_setup_krci
#include "implicit.inc"
*. General input  
#include "parluci.h"
*. Specific input 
      INTEGER NA(NSMST_in),NB(NSMST_in), IACTLIST(*)
*. Output
      INTEGER IOFF(NSMST_in)
*
      IB = 1
      IOFF_SCR = 0
      IOFF_SCR = IOFF_ACT
      IF( IVERS .eq. 0 ) THEN
        DO IASM = 1, NSMST_in
          IF( IACTLIST( IOFF_SCR ) .eq. 0 ) GOTO 100
          IBSM       = IDBGMULT(INVELM(IASM),ITOTSM)
          IOFF(IASM) = IB
          IB         = IB + NA(IASM)*NB(IBSM)
 100      CONTINUE
          IOFF_SCR = IOFF_SCR + 1
        END DO
      ELSE
        DO IASM = 1, NSMST
          IF( IACTLIST( IOFF_SCR ) .ne. MYPROC ) GOTO 200
          IBSM       = IDBGMULT(INVELM(IASM),ITOTSM)
          IOFF(IASM) = IB
          IB         = IB + NA(IASM)*NB(IBSM)
 200      CONTINUE
          IOFF_SCR = IOFF_SCR + 1
        END DO
      END IF
*
      END
***********************************************************************

      subroutine rdvcdc(VE,LU,LBLK,IRILP,IPRT)
*
* Read vector from disc to array
*
* Complex version      Timo Fleig, Dec 04, 2001
*
#include "implicit.inc"
*. Output
      DIMENSION VE(*)
*
      NTESTL = 00
      NTEST = max(NTESTL,IPRT/1000)
*
      CALL REWINE(LU,LBLK)
*
*
*. Loop over blocks of file
*
      if (LBLK.gt.0) then
        LENGTH = LBLK
        call frmdsc(VE,LENGTH,-1,LU,IMZERO,IAMPACK)
        stop 'Take care: Check me first!!.'
      else
        NASS = 0
        IOFF = 1
        do IRI = 1,IRILP,1
 1000     continue
          call ifrmds(LENGTH,1,LBLK,LU)
          if (LENGTH.ge.0) then
            call frmdsc(VE(IOFF),LENGTH,LBLK,LU,IMZERO,IAMPACK)
            if (NTEST.ge.10) then
              write(6,*) ' Record read in: '
              call wrtmat(VE(IOFF),1,LENGTH,1,LENGTH)
            end if
            IOFF = IOFF + LENGTH
          end if
          if (LENGTH.ge.0) then
            NASS = NASS + LENGTH
            goto 1000
          end if
*
          if (NTEST.ge.5) then
            if (IRI.eq.1) then
              write(6,*) ' Read vector in rdvcdc (real):  '
              call wrtmat(VE(1),1,NASS,1,NASS)
            else
              write(6,*) ' Read vector in rdvcdc (imag):  '
              call wrtmat(VE(1+NASS/2),1,NASS/2,1,NASS/2)
            end if
          end if
*
        end do
      end if
*
      end
***********************************************************************

      SUBROUTINE VECSUMC(C,A,B,FACAR,FACAI,FACBR,FACBI,NDIM)
*
* Multiply two complex vectors with scalars and add
*
*C(I) = FACA*A(I) + FACB*B(I)
*
* NDIM is length of real/imag part
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION C(NDIM,2)
*. Input
      DIMENSION A(NDIM,2),B(NDIM,2)
*
*
      DO I = 1, NDIM
*. Real part of C
        C(I,1) = FACAR * A(I,1)  + FACBR * B(I,1)
     &         - FACAI * A(I,2) -  FACBI * B(I,2)
*. Imaginary part of C
        C(I,2) = FACAR * A(I,2) + FACBR * B(I,2)
     &         + FACAI * A(I,1) + FACBI * B(I,1)
      END DO
*
      RETURN
      END
***********************************************************************

      SUBROUTINE WRSVCD_REL(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,NDIM,LUFORM,
     &                      NBLOCK,IBLOCKL,IBLOCKD,IAM_REAL,JROOT
#if defined (VAR_MPI2)
     &                      ,LU1LIST
#endif
     &                      )
*
* Write scattered vector to disc,
* every node writes only it's part.
* Information about the length of a blcok is stored in IBLOCKL.
* Information about which block belongs to which node is stored in
* IBLOCKD.
*
      use interface_to_mpi
      IMPLICIT REAL*8 (A-H,O-Z)
#if defined (VAR_MPI2)
#include "infpar.h"
      INTEGER   ISTAT(df_MPI_STATUS_SIZE)
#endif
#include "parluci.h"
*.Input
      DIMENSION IPLAC(*),VAL(*), IBLOCKL(*), IBLOCKD(*)
*.Scratch
      DIMENSION VEC1(*)
      LOGICAL FOUND_ONE
#if defined (VAR_MPI2)
      INTEGER(KIND=df_MPI_OFFSET_KIND) IOFF_SCR1, IOFF_SCR2
      DIMENSION LU1LIST(*)
#endif
C
      FOUND_ONE = .FALSE.
      IONE = 1
C
      IF(LBLK.GT.0) THEN
C        Write the vector without markers in one block
         CALL DZERO(VEC1,NDIM)
         DO 100 IEFF = 1, NSCAT
           VEC1(IPLAC(IEFF)) = VAL(IEFF)
  100    CONTINUE
         CALL TODSC(VEC1,NDIM,-1,LU)
      ELSE
C       Write the vector with the block format of file LUFORM
        CALL REWINE(LUFORM,-1)
        IBOT = 1
C       Loop over records
 1000   CONTINUE
#if defined (VAR_MPI2)
C       initialize offset
        IOFF_SCR1 = 0
        IOFF_SCR2 = 0
        INT_IOFF1  = 0
        INT_IOFF2  = 0
C       
        IOFF_SCR1 = MY_LU1_OFF  
     &            + MY_VEC2_IOFF * ( JROOT    - 1 ) 
     &            + MY_VEC1_IOFF * ( IAM_REAL - 1 )
        INT_IOFF1  = 1 + MY_ACT_BLK2 * ( JROOT    - 1 ) 
     &                 + MY_ACT_BLK1 * ( IAM_REAL - 1 )
CSK        WRITE(LUWRT,*) 'OFFSET at START',IOFF_SCR1
CSK        WRITE(LUWRT,*) 'INT_IOFF1 at START ',INT_IOFF1
C
        DO IBLK = 1, NBLOCK
          IF( IBLOCKD( IBLK ) .ne. MYPROC ) GOTO 300
CSK          WRITE(LUWRT,*) 'OFFSET new for block',IOFF_SCR1,IBLK -
CSK     &                                          INT_IOFF2, INT_IOFF2
CSK          WRITE(LUWRT,*) 'INT_IOFF1 new for block ',INT_IOFF1,IBLK -
CSK     &                                              INT_IOFF2, INT_IOFF2
          LBL = IBLOCKL( IBLK ) 
C         
#else
        CALL IFRMDS(LBL,1,-1,LUFORM)
        CALL ITODS(LBL,1,-1,LU)
        IF(LBL.GE.0) THEN
#endif
C?        write(6,*) ' IBOT,IBOT + LBL -1 ',IBOT,
C?   &    IBOT + LBL -1
          CALL DZERO(VEC1,LBL)
C         *******************************************************
C         *******************************************************
C                Find the places where a 1.0 should be set
C         *******************************************************
C         *******************************************************
          DO 200 IEFF = 1, NSCAT
            IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1) THEN
               VEC1(IPLAC(IEFF)-IBOT+1) = VAL(IEFF)
               FOUND_ONE = .TRUE.
            END IF
CSK          IF( IPLAC(IEFF).GE.IBOT.AND.IPLAC(IEFF).LE.IBOT+LBL-1)
CSK     &      write(luwrt,*) ' Catch : IPLAC(IEFF) VAL(IEFF) ',
CSK     &      IPLAC(IEFF),VAL(IEFF)
  200     CONTINUE
C         
#if defined (VAR_MPI2)
C
C         new offset
C
          IOFF_SCR1 = IOFF_SCR1 + IOFF_SCR2
          INT_IOFF1 = INT_IOFF1 + INT_IOFF2
C
C
          IF( FOUND_ONE ) THEN
CSK          WRITE(LUWRT,*) 'THIS IS WHAT I WILL WRITE',MYPROC
CSK          CALL WRTMATMN(VEC1,1,LBL,1,LBL,LUWRT)
CSK          WRITE(LUWRT,*) 'OFFSET',IOFF_SCR1
CSK          WRITE(LUWRT,*) 'INT_IOFF1 is ',INT_IOFF1
C
            LU1LIST( INT_IOFF1 ) = IONE
            CALL interface_MPI_FILE_WRITE_AT_r(ILU1,IOFF_SCR1,VEC1,LBL,
     &                             ISTAT)
C
          END IF
C
          FOUND_ONE = .FALSE.
#else
          CALL TODSC(VEC1,LBL,-1,LU)
#endif
C         Skip record on LUFORM
#if !defined (VAR_MPI2)
          CALL FRMDSC(VEC1,LBL,-1,LUFORM,IMZERO,IAMPACK)
          IBOT = IBOT + LBL
          GOTO 1000
        END IF
#else
C         keep track of correct offset
          IOFF_SCR2 = LBL
          INT_IOFF2 = IONE
  300     CONTINUE
C         keep IBOT up-to-date for all nodes
          LBL = IBLOCKL( IBLK )
          IBOT = IBOT + LBL
        END DO
#endif
C
      END IF
C
      END
***********************************************************************

      SUBROUTINE WRTVCDC(SEGMNT,LU,IREW,IRC,LBLK)
C
C PRINT COMPLEX VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        IF( LBLK .GE. 0 ) THEN
          REWIND LU
        ELSE
          CALL REWINE(LU,LBLK)
        END IF
      END IF
C LOOP OVER BLOCKS
C
      IBLK = 0
      DO IRI = 1,IRC,1
        IF(IRI.EQ.1) THEN
           WRITE(6,*) ' Real Part '
        ELSE
           WRITE(6,*) ' Imaginary part '
        END IF
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LU) LBL
        ELSE
          CALL IFRMDS(LBL,1,-1,LU)
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC(SEGMNT,LBL,KBLK,LU,IMZERO,I_AM_PACKED)
           IF(LBL .GT. 0 ) THEN
             WRITE(6,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,6)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      END DO
*. End of loop over real/imaginary parts
      RETURN
      END
***********************************************************************

      SUBROUTINE ZERORC_REL(MBLOCK,IFIL)
*
* A record was known to be identical  zero
*
*
* Write corresponding info to file IFIL
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISCR(2)
* Zero record
      ISCR(1) = 1
!  no packing
      ISCR(2) = 0
*
      CALL ITODS(ISCR,2,2,IFIL)
*
      END
