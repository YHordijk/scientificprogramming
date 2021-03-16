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

       SUBROUTINE MAKE_H2C(VMAT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  On input : 
C  =========
C       VMAT - Lowdin transformation matrix (4c SA-AO -> 4c UKB MO)
C
C On output :
C ===========
C   VMAT - Lowdin transformation matrix is adapted for 2-comp. SCF calculations,
C    or unchanged is on output
C
C    During this routine H2c AO matrix elements as well as
C picture change transformation matrixes are written to the file no=LUBSS 'BSSMAT'
C
C  Written by Miro Ilias & Hans Jorgen Aa Jensen, Odense, 2002, 
C     reprogrammed by Miro Ilias, Strasburg, 2005, 2006.
C
C****************************************************************************
C
C    '  Dobrorec dusa moja Hospodinu, a cele moje vnutro 
C     Jeho svatemu menu !' (Zalm 103,1)
C
C    '  Praise the Lord, my soul ! All my beeing,
C     praise his holy name!'  (Psalm 103,1)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbdhf.h"

      DIMENSION VMAT(*),WORK(LWORK)
      LOGICAL DOBSS,V2COUT,FNDLAB,TWOBUF
      real(8)              :: timstr(2), timend(2)

#include "memint.h"

      CALL QENTER('MAKE_H2C')
C ... generate 2-component Hamiltonian
C
      CALL TIMER2('START ',TIMSTR,TIMEND)
      CALL MAKE_BSS(VMAT,WORK(KFREE),LFREE) 
      CALL TIMER2('make_bss',TIMSTR,TIMEND)

C    ===================================================================================
C     ... decide whether do the picture change transformation of 4c Fock AO matrixes
C    ===================================================================================
C    ... it's using the picture change transformation matrix obtained from the
C    transformation of the Dirac bare nucleus or of the Fock-Dirac matrix
C    (when USEDF=.true.)
!     print *,'TWOCOMPBSS and I2COFK ==> ',TWOCOMPBSS , I2COFK
      IF (DO4C2C.AND..NOT.START2C) THEN

       IF (IPRHAM.GE.0) THEN
        write(LUPRI,'(/2X,A)')
     &  'MAKE_H2C: Going to do the picture change transformations'//
     &  ' of four-component Fock-Dirac AO matrixes'//
     &  ' (renews files DFFCK1, DFFCK2)'
        IF (USEDF) THEN
C       ... in this case we block-diagonalize the Fock-Dirac operator
         write(LUPRI,'(4X,A/)')
     &   '... the picture change transformation matrix comes from'//
     &   ' the block-diagonalization of the FOCK-DIRAC operator.'
        ELSE
         write(LUPRI,'(4X,A/)')
     &   '... the picture change transformation matrix comes from the'//
     &   ' block-diagonalizations of Dirac BARE-NUCLEUS operator.'
        ENDIF
       ENDIF

       ISIZE = NORB(1)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)   
       CALL MEMGET('REAL',KIBEIG,ISIZE,WORK,KFREE,LFREE)
       CALL FDTR2C(VMAT,WORK(KIBEIG),WORK(KFREE),LFREE)
!      CALL MEMCHK('MAKE_H2C.FDTR2C ',WORK,1)

       IF (TRMO4C2C) THEN 
!mi    ... transform 4c MO coeff -> 2c MO coeff using the picture change tranf.matrix
CMI        .... consider NFMAT for KBUF array allocation !!!
          CALL MEMGET('REAL',KCMO,N2BBASXQ,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KBUF,N2BBASXQ,WORK,KFREE,LFREE)
          CALL MEMGET('REAL',KEIG ,NTBAS(0),      WORK,KFREE,LFREE)
          IF (NFMAT.GT.1) 
     &       CALL QUIT('TRMO4C2C: open-shell not yet implemented')
          CALL CMOTR2C(VMAT,WORK(KBUF),WORK(KCMO),WORK(KEIG),
     &                 WORK(KIBEIG),WORK(KFREE),LFREE)
!         CALL MEMCHK('MAKE_H2C.CMOTR2C ',WORK,1)
          CALL MEMREL('MAKE_H2C.CMOTR2C 1',WORK,KWORK,KEIG,KFREE,LFREE) 
          CALL MEMREL('MAKE_H2C.CMOTR2C 2',WORK,KWORK,KBUF,KFREE,LFREE) 
          CALL MEMREL('MAKE_H2C.CMOTR2C 3',WORK,KWORK,KCMO,KFREE,LFREE) 
       ENDIF

       CALL MEMREL('MAKE_H2C.FDTR2C',WORK,KWORK,KIBEIG,KFREE,LFREE)

      ENDIF

C ... store transformed H2c_AO Hamiltonian in pure "LL" form ...
      IF (TWOCOMPBSS.AND.(I2COFK.EQ.1.OR.I2COFK.EQ.2)) THEN
        CALL H2CTOLL(I2COFK,VMAT,WORK(KFREE),LFREE,IPRHAM)
      ENDIF

C ... for calculating properties in pure 2comp. mode do adapt 4comp.operators
      IF (TWOCOMPBSS) THEN
        CALL PROP2BSS(VMAT,IPRHAM,WORK(KFREE),LFREE)
      ENDIF
C ===============================================================================
C ... export the Lowdin matrix adapted for 2c calculations,
C ===============================================================================
      CALL RFBSSMAT('VMAT2C  ',.TRUE.,'MAKE_H2C',
     &               VMAT,N2BBASX,LUBSS,IPRHAM)

C    ... close the BSSMAT file ...
      CLOSE(LUBSS,STATUS='KEEP')

C ============================================================
C  Finally set up number of positronic shells equal to zero
C ============================================================
      NSYM = 4/NZ
      DO I=1,NFSYM
       NPSH(I) = 0
       NORB(I) = NESH(I) + NPSH(I)
C  ... new variable ....
       NTMO(I) = NORB(I)
       DO ISYM = 1, NSYM
         NBORB(ISYM,I,2) = 0
         NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
       ENDDO
      ENDDO
C     ... recalculate necessary off-sets
      CALL SETDC2(0)
C **********************************************************************************
C    ... Number of transf. matrixes is set to 1 - we use only the Lowdin AO->MO
      NZT = 1
C **********************************************************************************
      WRITE(LUPRI,'(/2X,A)') 
     &'MAKE_H2C: Two-component BSS mode - positronic shells deleted !'
C   .... do set up appropriate integrals...
      INTDEF_SAVE = INTDEF
      INTDEF = 1
      INTGEN = INTDEF
      IF (IPRHAM.GE.3) THEN
        WRITE(LUPRI,'(2x,a)') 
     & 'MAKE_H2C: SL,SS,GT integrals not calculated !'
      ENDIF
C     ... set number of components to 1 (=> two-component mode)      
      MC = 1
      SSMTRC_SAVE = SSMTRC
      SSMTRC = 0.0D0
C
      CALL QEXIT('MAKE_H2C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck geth2cao */
       SUBROUTINE GETH2CAO(I2COFK,FOCK,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Purpose: routine reads selected type of the two-component transformed Hamiltonian
C           depending on the entering I2OFCK value
C
C  Called from: ONEFCK (if BSS=.true.)
C
C On input :   
C ==========
C  FOCK - array used inside
C  I2COFK - determines the integrals type (ie what kind of 2-comp. transformed Hamiltonian)
C
C  Note that it the temporary default, BSS is really run in 4-component
C  mode, that is the small components are simply set to zero.
C
C  True 2-component mode means that the small components are fully absent.
C
C * BSS transformation of one-electron problem:
C    1 - 2c matrix from 4c: h2c from h_D(1)_4C
C    2 - 2c matrix from 4c: h2c_CSC(+AMFI if present) from hD(1)_4C
C
C * BSS transformation after 4c SCF, that is using the mean-field potential:
C    3 - 2c one-electron matrix 
C    4 - 2c full Fock matrix
C    5,6 - 2c h2c+u1_so from 2c full Fock matrix
C    7 - true 2-component mode
C
C  IPRINT - print level
C
C On output :     BSSMAT file updated with LL H2c_AO integrals
C ===========     I2COFK set to 7
C          
C Written by Miro Ilias, Strasburg 2006  
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"

      DIMENSION FOCK(*),WORK(LWORK)
      LOGICAL FNDLAB
      CHARACTER*8 HAM2CL
C
C     HAM2CL is the label of the two-comp. hamiltonian matrix to read
C     ==================================================================
C

!     print *,'entered GETH2CAO, reading type  ==> ',I2COFK
C     Select the label according to the I2COFK parameter
      IF      (I2COFK.EQ.1) THEN       
C.... h2c from h_D(1)_4C ... DEFAULT ...    
        HAM2CL = 'H1AO_DK ' 
      ELSE IF (I2COFK.EQ.2) THEN      
C.... h2c_CSC(+AMFI,if present) from hD(1)_4C                                     
        HAM2CL = 'H1AO_DK2'     
C.... these are post-DC-SCF two-component Hamiltonians ...                         
C     =============================================================                       
      ELSE IF (I2COFK.EQ.3) THEN 
C.... h2c(1) from FD4C - one electron terms !                                          
C.... WE MAY WISH TO PERFORM ONLY ONE SCF ITERATION HERE
        HAM2CL = 'FD1AO_DK' 
      ELSE IF (I2COFK.EQ.4) THEN 
C.... full Fock_2c transformed from FD4C,ie from previous 4c SCF  
        HAM2CL = 'FDAO_DK ' 
      ELSE IF (I2COFK.EQ.5.OR.I2COFK.EQ.6) THEN
C.... f2c(=h1+u1_so, spin-free terms are removed from u1) from FD_4C term              
        HAM2CL = 'FD2C_USO' 
      ELSE IF (I2COFK.EQ.7) THEN    
C.... it's H1AO_DK or H1AO_DK2 in LL basis (ie pure two-component modus)
        HAM2CL = 'H2CAO_LL' 
      ELSE     
        WRITE(LUPRI,*) 'H2CAOLL: I2COFK=',I2COFK   
        CALL QUIT(
     &  'ONEFCK: Wrong value of the entering I2COFK!')
      ENDIF  
!     print *,'entered GETH2CAO, reading label ==> ',HAM2CL

      IF (IPRINT.GE.0) THEN
       WRITE(LUPRI,'(/,2X,A,I1,A,A)')
     & 'GETH2CAO: asked for 2c Hamiltonian (param.I2COFK=',I2COFK,'):',
     & HAM2CL
      ENDIF

C
C     Read chosen two-component hamiltonian integrals from file BSSMAT
C     ==================================================================
C
      CALL RFBSSMAT(HAM2CL,.FALSE.,'GETH2CAO',
     &      FOCK,N2BBASXQ,LUBSS,IPRINT)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck make_bss */
       SUBROUTINE MAKE_BSS(VMAT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C On input : 
C =========
C  VMAT - Lowdin transformation matrix (AO to MO) for 4-comp. Fock matrix
C
C On output :
C ===========
C  none, but VMAT - Lowdin transformation matrix is adapted
C  for 2-comp. SCF calculations and is written to the file BSSMAT.
C          
C Written by Miro Ilias & Hans Jorgen Aa Jensen, Odense, 2002, 
C        deeply reprogrammed by Miro Ilias, Strasburg, 2005.
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     use x2c_utils, only: 
!    &    print_x2cmat
      use common_matvec_op
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbdhf.h"

      DIMENSION VMAT(*),WORK(LWORK)
      LOGICAL DK2SO1,ONLYSO1,SAVESO1,SF_BEG,FPB,FNDLAB,
     &    GETPCTM,ADD_AMFI,ONLYTRA,ELIM_SF, WRITEH2C, LINAD
      integer is_defining_h1mat

        CALL QENTER('MAKE_BSS')
#include "memint.h"

#ifdef MOD_OOF
C      CALL TEST_OOF(WORK(KFREE),LFREE)
#endif /* MOD_OOF  */
C
C     Print section
C
      IF (IPRHAM.GE.1) THEN
        CALL HEADER('*** Output from MAKE_BSS ***',-1)
      ENDIF
C ======================================================================
C            extract information about the BSS Hamiltonian 
C ======================================================================
      CALL EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,SAVESO1,ADD_AMFI,
     &     ONLYSO1,ONLYTRA,SF_BEG,ELIM_SF,NOPRTR,BLOCKD,
     &     DO4C2C,USEDF,CONT2C,START2C,I2CHAM,DO2C4C,INI2C,I2COFK,
     &     YREQ1,IC_AMFI,GETPCTM,NOAMFI,IPRHAM)
C  ... do not compare eigenvalues for ONLYSO1 = T, what is testing mode
C  on BSS-SO1 integrals against AMFI SO1 terms
      CMPEIG = CMPEIG.AND..NOT.ONLYSO1
C     ... initialize the important FreeParticleBasis flag!
C     If .false., we use the one-step RKB approach.
      FPB = .NOT.NOPRTR

C  allocate multipurpose arrays 
      CALL MEMGET('REAL',KTBUF,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KTMAT,N2BBASXQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG ,NTBAS(0),WORK,KFREE,LFREE)

      IF (IPRHAM.GE.5) THEN
        CALL HEADER('MAKE_BSS: Entering Lowdin matrix',-1)
        CALL PRQMAT(VMAT,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),1,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF

C     ... choose the file number !
      LUBSS = 40
      ! .. save 4c Lowdin transformation matrix  
      CALL W2BSSMAT(4,.TRUE.,'VMAT4C  ','MAKE_BSS',
     &              VMAT,N2BBASX,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_BSS: 4c Lowdin AO to MO transf. matrix is '//
     & ' written to the BSSMAT file.'
      ENDIF
C ============================================================
C ... get the RKB transformation matrix into WORK(KTMAT)
C ============================================================
!     flag for the adaptation for the linear symmetry
      LINAD = ( DO4C2C.AND.LINEAR ) 
!     defining h1 matrix control variable (in the X2C module) - 
!     set here for consistency to the maximum value (the true value 
!     needs to be checked though)
      is_defining_h1mat = 4
      CALL MODHAM(WORK(KTMAT),VMAT,WORK(KTBUF),WORK(KEIG),.NOT.FPB,
     &            LINAD,is_defining_h1mat,WORK(KFREE),LFREE)
C.....set variables needed in DFDIAG (N2TMT etc.)
      CALL SET_TMT(IPRHAM)
C     1. Write transformation matrices to file BSSMAT
      CALL HANDLE_RKB(VMAT,WORK(KTMAT),WORK(KTBUF),
     &            LINAD,FPB,IPRHAM,WORK(KFREE),LFREE)

C      ... allocate more memmory 
      CALL MEMGET('INTE',KIBEIG,NTBAS(0),WORK,KFREE,LFREE)
      !CALL IZERO(WORK(KIBEIG),NTBAS(0))
      CALL MEMGET('REAL',KFACKIN,NTBAS(1),WORK,KFREE,LFREE)

C ===============================================================
C ===       Prepare the free particle (FP) basis "k"
C ===============================================================
      IF (FPB) THEN
C generate free particle one-electron matrix in RKB ON basis "theta"
C    ... into WORK(KTBUF), reuse VMAT...
C   ... merge into 1 routine ! - get + diagonalize !
C diagonalize the FP matrix in the ON "theta" basis
        CALL MAKE_FPM(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KEIG),
     &             WORK(KIBEIG),WORK(KFREE),LFREE)
      ELSE
C ===============================================================
C ===   We have the RKB basis, prepare the boson irreps info
C ===============================================================
       IF (SPINFR.OR.DO4C2C) THEN
        CALL INI_SF_RKB(WORK(KIBEIG),IPRHAM)
       ENDIF
C       ... we need to zero them as in enters GET_A_FAC routine !
       CALL DZERO(WORK(KEIG),NTBAS(0))
      ENDIF
C NOTE: If SPINFR=.true. we utilize later the WORK(KIBEIG) array with boson info !!!

C get A factors into WORK(KFACKIN), nothing enters...
      CALL GET_A_FAC(WORK(KFACKIN),WORK(KEIG),VMAT,WORK(KTMAT),
     &              WORK(KTBUF),FPB,USEDF,IPRHAM,WORK(KFREE),LFREE)

C ===============================================================
C    ... Do prepare the MO->AO transformation matrix ...
C ===============================================================
C        ...entering WORK(KTMAT) contains the transformation matrix in AO
      CALL MKSAOMO_Le(WORK(KTBUF),WORK(KTMAT),VMAT,
     &         WORK(KFACKIN),FPB,IPRHAM,WORK(KFREE),LFREE)

!     CALL MEMCHK('MAKE_BSS 1',WORK,1)

C get H1 (Dirac bare nucleus) Hamiltonian in the ON basis into WORK(KBUF)
C    ... that is  the RKB  or FP "k" basis; WORK(KTMAT) entering ... !
      CALL GET_H1_ONB(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KIBEIG),
     &                   SF_BEG,FPB,LINAD,SAVESO1,WORK(KFREE),LFREE)

C ===================================================================
C    get H2C Hamiltonian(inf,DKH2,DKH2SO1,DKH1)  in ON ("k"/RKB) basis
C ===================================================================
      IF (IPZ.EQ.9) THEN
C       ... infinite order - either block diagonalization, or YR equation solving
#ifdef MOD_UNRELEASED
        IF (BLOCKD) THEN
         CALL GET_HINF_BD(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KIBEIG),
     &                    GETPCTM,SF_BEG,WORK(KFREE),LFREE,IPRHAM)
        ELSE
#endif
C        ... we need additional integer array...
         ISIZE=NORB(1)
         IF (NFSYM.EQ.2) THEN
          IF (ISIZE.LT.NORB(2)) ISIZE=NORB(2)
         ENDIF

C        ... both for FPM "k" & RKB 
         CALL GET_HINF_R(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KFACKIN),
     &                   WORK(KEIG),WORK(KIBEIG),FPB,YREQ1,
     &               GETPCTM,SF_BEG,WORK(KFREE),LFREE,IPRHAM)

#ifdef MOD_UNRELEASED
        ENDIF
#endif
C =======================================================================
C       ... second order DKH - only in the free-particle basis
C =======================================================================
      ELSE IF (IPZ.EQ.2) THEN
        CALL GET_HDK2_FPM(WORK(KTBUF),WORK(KTMAT),VMAT,WORK(KEIG),
     &                    WORK(KIBEIG),DK2SO1,
     &                    IPRHAM,WORK(KFREE),LFREE)
      ELSE
C       ...the roughest first order, or SO1 integrals; only FP "k" basis 
        CALL GET_HDK1_FPM(WORK(KTBUF),VMAT, WORK(KIBEIG),ONLYSO1,
     &                            IPRHAM,WORK(KFREE),LFREE)
      ENDIF

C  Do  transform H2c RKB/"k" -> H2c "xhi" LL, save final elements into the file
C     ... to be specified later !!!!
      WRITEH2C = .NOT.(USEDF.AND.DO4C2C.AND..NOT.START2C)
      CALL H2CFINAL(WORK(KTBUF),WORK(KTMAT),VMAT,
     &            WORK(KFACKIN),WORK(KIBEIG),WRITEH2C,SAVESO1,
     &      SPINFR,SF_BEG,GETPCTM,FPB,IPRHAM,WORK(KFREE),LFREE)

C      ... after preparing final H2c Hamiltonian do renew the full relatvity !
      IF (SAVESO1.AND.SPINFR) THEN
       write(lupri,'(/,2X,A)') 'MAKE_BSS: SPINFR set to false !'
       SPINFR = .FALSE.
       SF_BEG = .FALSE.
      ENDIF

C   ===============================================================
C     ... prepare picture change transformation matrixes...
C   ===============================================================
      IF (GETPCTM) THEN
        CALL GETPCTMAT(WORK(KTBUF),WORK(KTMAT),VMAT,FPB,USEDF,
     &               IPRHAM,WORK(KFREE),LFREE)
      ENDIF
C
C     Prepare Lowdin matrix on 2c form
C
      CALL LOWD2C(VMAT,WORK(KTBUF),WORK(KTMAT),WORK(KEIG),
     &            WORK(KFREE),LFREE,IPRHAM)
C
C transform H2C Hamitonian from FP "k" basis to the SO-AO  basis,
C prepare the final picture change transformation matrix if wished,
C save integrals (and pctmtx) into the file
      IF(CMPEIG) THEN
C       ... compare eigenvalues of 4c and 2c 1-el.Hamiltonians
        CALL CMP_EIGVAL(VMAT,WORK(KTMAT),WORK(KTBUF),WORK(KEIG),
     &              WORK(KIBEIG),WRITEH2C,SPINFR,SF_BEG,
     &          ONLYSO1,SAVESO1,USEDF,DO4C2C,START2C,
     &             IPRHAM,WORK(KFREE),LFREE)
      ENDIF

C============================================================
C ... prepare the FULL picture change transformation matrix
C============================================================
      IF (DO2C4C.AND.START2C) THEN
C      ... consider the spin-free case !!!
        CALL FULLPCTM(VMAT,WORK(KTMAT),WORK(KTBUF),WORK(KEIG),
     &              IPRHAM,WORK(KFREE),LFREE)
      ENDIF

C     ...      before calling AMFI stuff release ALL
C     ...    allocated memmory (except entering VMAT)
      CALL MEMREL('AMFI',WORK,KWORK,KTBUF,KFREE,LFREE)

      IF (ADD_AMFI) THEN
C    ... transformation cartesian -> spherical, calling the AMFI stuff,
C        fitting AMFI integrals in the DIRAC order, saving to the BSSMAT
       IF (SPINFR.AND.ELIM_SF) SPINFR=.FALSE.
C      ... VMAT is reused as helping array...
        CALL AMFI_STUFF(VMAT,ONLYSO1,ONLYTRA,IC_AMFI,
     &                 WORK(KFREE),LFREE,IPRHAM,IMFCH)
        I2CHAM = 2
      ENDIF

C     ... prepare the screened 2comp. bare nucleus - only for the infinite order !
      IF (BNCRON .AND. .NOT. ONESYS .AND. GETPCTM) THEN
       CALL MAKEBNCR2C(VMAT,WORK(KFREE),LFREE)
      ENDIF

C   ... release the rest of allocated memmory ...
      CALL MEMREL('MAKE_BSS',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('MAKE_BSS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      SUBROUTINE SET_TMT(IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Set up important variables - needed for DFDIAG routine !
C
C  Written by Miro Ilias, Prievidza, December 2006
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbdhf.h"

C      ... NTMO is set up !!!!
C .... this settings are needed (DFDIAG!!!) for the RKB option !!!
C .. NZT is in dcbbas.h, and after MODHAM for RKB is equal to NZ !
      N2TMT   = 0
      N2TMOTQ = 0
      NORBT  = 0
      N2ORBT = 0
      DO IFRP = 1,NFSYM
        NISHMF(IFRP) = NISH(IFRP)
        NOCCMF(IFRP) = NOCC(IFRP)
        I2TMT(IFRP)  = N2TMT
        N2TMT        = N2TMT + NFBAS(IFRP,0)*NORB(IFRP)*NZ
        I2TMOT(IFRP) = N2TMOTQ
        N2TMO(IFRP)  = NTMO(IFRP)*NTMO(IFRP)
        N2TMOTQ      = N2TMOTQ + N2TMO(IFRP)*NZ
        IORB(IFRP)   = NORBT
C       ...  Count total number of orbitals  ...
        NORBT        = NORBT  + NORB(IFRP)
        I2ORBT(IFRP) = N2ORBT*NZ
        N2ORBT       = N2ORBT + NORB(IFRP)*NORB(IFRP)
      ENDDO

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/,6X,A)') '**** SET_TMT: set variables ****'
       DO I = 1, NFSYM
        write(lupri,'(1X,A,I1,A,I1)') 'fermion symm:',I,'/',NFSYM
        write(lupri,'(3X,A,I1,A,I4)') 'NISHMF(',I,')=',NISHMF(I)
        write(lupri,'(3X,A,I1,A,I4)') 'NOCCMF(',I,')=',NOCCMF(I)
        write(lupri,'(3X,A,I1,A,I4)') 'I2TMT(',I,') =',I2TMT(I)
        write(lupri,'(3X,A,I1,A,I4)') 'IORB(',I,')  =',IORB(I)
        write(lupri,'(3X,A,I1,A,I4)') 'I2ORBT(',I,')=',I2ORBT(I)
       ENDDO
       write(lupri,'(4X,A,I3)') 'N2TMT=',N2TMT
      ENDIF

      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck handle_rkb */
      SUBROUTINE HANDLE_RKB(VMAT,TMAT,TBUF,LINAD,FPB,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  After the MODHAM routine take care of matrices
C
C  On input:
C ------------
C    VMAT - (non-linear) SL resorted RKB AO to MO (for FPB = False)
C    TMAT - Non-linsym RKB transformation matrix 4c AO->MO
C    TBUF - RKB MO to LINSYM RKB_MO (LINAD = T)
C
C  Written by Miro Ilias, Prievidza, December 2006
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
      DIMENSION VMAT(*),TMAT(*),TBUF(*),WORK(LWORK)
      LOGICAL LINAD,FPB

      CALL QENTER('HANDLE_RKB')

      IF (IPRINT.GE.5) THEN
        IMAT =  1
        IVMAT = 1
        IBUF = 1
        DO I = 1, NFSYM                                              
        IF(NORB(I).GT.0) THEN 
          CALL HEADER(
     &  'HANDLE_RKB: Non-linsym RKB transformation matrix 4c AO->MO',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion ircop no.',I,'/',NFSYM                          
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),              
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)        

          IF (LINAD) THEN
           CALL HEADER(
     & 'HANDLE_RKB: RKB MO to LINSYM RKB_MO: ',-1)
           WRITE(LUPRI,'(A,I1,A,I1)')
     &     '* Fermion ircop no.',I,'/',NFSYM                          
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),              
     &          NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)        
          ENDIF

          IF (.NOT.FPB) THEN
           CALL HEADER(
     & 'HANDLE_RKB: (non-linear) SL resorted RKB AO to MO  ',-1)
           WRITE(LUPRI,'(A,I1,A,I1)')
     &     '* Fermion ircop no.',I,'/',NFSYM                          
           CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NORB(I),              
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)        
          ENDIF

          IBUF = IBUF  + NORB(I)*NORB(I)
          IMAT = IMAT  + NFBAS(I,0)*NORB(I)*NZ                         
         IVMAT = IVMAT + NFBAS(I,0)*NORB(I)*NZ  

        ENDIF                                                   
        ENDDO                                                         
      ENDIF

      IF (.NOT.FPB) THEN
C      ... save the 'SL' resorted transformation matrix
       CALL W2BSSMAT(4,.TRUE.,'SL_TM_4C',
     & 'HANDLE_RKB',VMAT,N2TMT,LUBSS,IPRINT)
       IF (IPRINT.GE.2) THEN
        write(lupri,*)
     & 'HANDLE_RKB: SL_TM_4C (VMAT) was written to the BSSMAT !!'
       ENDIF
      ENDIF

      CALL W2BSSMAT(4,.TRUE.,'TMAT4C  ','HANDLE_RKB',TMAT,
     &               N2TMT,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     & 'HANDLE_RKB: RKB transf. matrix (TMAT, basic) '//
     & 'has been written to the BSSMAT file.'
      ENDIF

      IF (LINAD) THEN
        ISIZE = NORB(1)*NORB(1)
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)
        CALL W2BSSMAT(4,.TRUE.,'TMAT_LSY','HANDLE_RKB',
     &      TBUF,ISIZE,LUBSS,IPRINT)
        IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,'(/2X,A)')
     &      'HANDLE_RKB: RKB MO to MO_LSYM (WORK(KTBUF)) '//
     &      'was written to the BSSMAT file.'
        ENDIF
      ENDIF

C       ... dummy vars
      CALL W2BSSMAT(6,.TRUE.,'XXXXXXXX','HANDLE_RKB',
     &               VMAT,N2TMT,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'HANDLE_RKB: BSSMAT file - '//
     & ' EOFLABEL added.'
      ENDIF

      CALL QEXIT('HANDLE_RKB')
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck lowd2c */
      SUBROUTINE LOWD2C(VMAT,TBUF,TMAT,EIG,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Reorder the 4c-Lowdin matrix for 2component calculations and write it to the file
C
C  On input: arrays to be utilized inside, and flags
C  
C
C  Written by Miro Ilias, December 2006, Prievidza,
C                        2009, Tel Aviv - reaarangement for linear symmetry
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
      DIMENSION VMAT(*),TBUF(*),TMAT(*),EIG(*),WORK(LWORK)
      DIMENSION NBORBX(4, 2, 0:2),NPSHX(2)
      LOGICAL PURLO, WANT_LINEAR
#include "memint.h"

      CALL QENTER('LOWD2C')

C  ===================================================================
C     ... get AND reorganize the Lowdin transformation matrix VMAT
C  ===================================================================
      CALL RFBSSMAT('VMAT4C  ',.TRUE.,'LOWD2C',
     &               VMAT,N2BBASX,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/2X,A)')
     &  'LOWD2C: VMAT4C was read from BSSMAT.'
      ENDIF

      IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'LOWD2C: The total Lowdin matrix before BSS reordering',-1)
        IVMAT=1
        DO I=1, NFSYM
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IVMAT),
     &         NFBAS(I,0),NFORB(I,0),NFBAS(I,0),NFORB(I,0),
     &         1,IPQTOQ(1,0),LUPRI)
          IVMAT = IVMAT + NFBAS(I,0)*NFORB(I,0)
        ENDDO
      ENDIF

      PURLO = .TRUE.
      IF (PURLO) THEN
C  ...  purify the Lowdin matrix !!!
C  ...  1. remove zeros for LL !
       IVMAT = 1
       DO I=1,NFSYM
        IPOS=IVMAT
        DO J=1,NFORB(I,1)
         CALL DZERO( VMAT(IPOS+NFBAS(I,1)),NFBAS(I,2) )
         IPOS = IPOS + NFBAS(I,0)
        ENDDO
        IPOS = IVMAT + (NFBAS(I,0)*NFORB(I,1))
        CALL DZERO( VMAT(IPOS), (NFORB(I,2)*NFBAS(I,0)) )
        IVMAT = IVMAT + ( NFBAS(I,0)*NFORB(I,0) )
       ENDDO
      ENDIF

C ===================================================================
C  ... reorder the Lowdin matrix for 2component Hamiltonians ...
C ===================================================================
      IF (NFSYM.EQ.2) THEN
C    HJJ:For TMAT: copy the VMAT for the L block ifsym=2 to
C        follow immediately after the L block ifsym=1.
           NDIM2 = NFBAS(2,0)*NESH(2)
           JV2   = 1 + NFBAS(1,0)*NFORB(1,0)
           JT2   = 1 + NFBAS(1,0)*NESH(1)
           CALL DCOPY(NDIM2,VMAT(JV2),1,VMAT(JT2),1)
      ENDIF

      IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &'LOWD2C: Lowdin matrix (symm.blocked) AFTER BSS reordering',-1)
        IVMAT=1
        DO I=1, NFSYM
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IVMAT),
     &          NFBAS(I,0),NFORB(I,0),NFBAS(I,0),NFORB(I,0),
     &          1,IPQTOQ(1,0),LUPRI)
          IVMAT = IVMAT + NFBAS(I,0)*NFORB(I,0)
        ENDDO

        CALL HEADER(
     &'LOWD2C:Lowdin matrix (symm.blocked -Le) AFTER BSS reordering',-1)
        IVMAT=1
        DO I=1, NFSYM
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IVMAT),
     &          NFBAS(I,0),NESH(I),NFBAS(I,0),NESH(I),
     &          1,IPQTOQ(1,0),LUPRI)
          IVMAT = IVMAT + NFBAS(I,0)*NESH(I)
        ENDDO
      ENDIF

C ==========================================================
C .... write the 2c-adapted Lowdin matrix into the file
C .. this is the non-linear symmetry case !!!!
C ==========================================================
      CALL W2BSSMAT(1,.TRUE.,'VMAT2C  ','LOWD2C',
     &               VMAT,N2BBASX,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'LOWD2C: Adapted (non-linsym) Lowdin matrix, VMAT2C,'//
     & ' was written into the BSSMAT file, EOFLABEL was renewed.'
      ENDIF

      WANT_LINEAR=.FALSE.
      IF (LINEAR) THEN 
!     IF (WANT_LINEAR) THEN 
C     ... for the linear symmetry DO ADAPT the Lowdin matrix !
C      ... the proper setting for the 2c mode is needed !
C   ... when we have Fock-Dirac matrix SCF in linear symmetry, we need
C      the H2c in LINSYM_MO basis !!!
      NSYM = 4/NZ     
      DO I=1,NFSYM     
       NPSHX(I) = NPSH(I)
       NPSH(I) = 0    
       NORB(I) = NESH(I) + NPSH(I)     
       DO ISYM = 1, NSYM     
         NBORBX(ISYM,I,2) =  NBORB(ISYM,I,2)     
         NBORB(ISYM,I,2) = 0      
         NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)     
       ENDDO      
      ENDDO
      CALL LINSYM(TMAT,VMAT,TBUF,EIG,WORK(KFREE),LFREE)
      CALL DCOPY(N2BBASX,TMAT,1,VMAT,1)
      DO I=1,NFSYM     
       NPSH(I) = NPSHX(I)
       NORB(I) = NESH(I) + NPSH(I)     
       DO ISYM = 1, NSYM     
         NBORB(ISYM,I,2) =  NBORBX(ISYM,I,2)     
         NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)     
       ENDDO      
      ENDDO  

      IF (IPRINT.GE.6) THEN 
        CALL HEADER(
     &'LOWD2C: Total Lowdin 2c_matrix AFTER LINEAR adaptation',-1)
        IVMAT=1
        DO I=1, NFSYM
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),
     &        NFBAS(I,0),NFORB(I,0),NFBAS(I,0),NFORB(I,0),
     &        1,IPQTOQ(1,0),LUPRI)
          IVMAT = IVMAT + NFBAS(I,0)*NFORB(I,0)
        ENDDO
      ENDIF
C .... write the 2c-adapted lin-sym adapted Lowdin matrix into the file
      CALL W2BSSMAT(1,.TRUE.,'VMAT2C_L','LOWD2C',
     &              VMAT,N2BBASX,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'LOWD2C: 2c LinSym adapted Lowdin matrix, VMAT2C_L,'//
     &  ' was written to the BSSMAT file, EOFLABEL was renewed.'
      ENDIF

      ENDIF

      CALL QEXIT('LOWD2C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      SUBROUTINE MAKEBNCR2C(VMAT,WORK,LWORK)
C****************************************************************************
C
C   PURPOSE: Do the picture change transformation of the SCREENED
C            one-electron Fock matrix obtained from the BNCORR routine.
C
C   Called from: MAKE_BSS when the picture change transformation is available 
C               (we have infinite order)
C
C   On input: VMAT - only array
C
C   On output: VMAT - picture change transformed 2comp. screened bare nucleus
C
C   Written by Miro ILIAS, March 2006/Strasbourg
C
C****************************************************************************
!     use x2c_utils, only:
!    &    print_x2cmat
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbham.h"
      DIMENSION VMAT(*),WORK(LWORK)
      LOGICAL BSS_SAVE,TWOCOMP_SAVE

#include "memint.h"
      CALL QENTER('MAKEBNCR2C')
      
      BSS_SAVE = BSS
      BSS = .FALSE.
      TWOCOMP_SAVE = TWOCOMP
      TWOCOMP = .FALSE.

CMI   ... get Dirac 1el matrix
      CALL ONEFCK(VMAT,IPRHAM,WORK(KFREE),LFREE)
!     call print_x2cmat(VMAT,NTBAS(0),NTBAS(0),nz,IPQTOQ(1,0),
!    &     'bss - 4ch1sbc',6)      
CMI   ... get screened Dirac bare nucleus
      CALL BNCORR(VMAT,WORK(KFREE),LFREE,IPRHAM)
C
!     call print_x2cmat(VMAT,NTBAS(0),NTBAS(0),nz,IPQTOQ(1,0),
!    &     'bss - 4ch1bnc',6)      
C
      BSS = BSS_SAVE
      TWOCOMP = TWOCOMP_SAVE

      ISIZE = 0
      DO IFRP = 1,NFSYM
C.......................[    L+S    ]*[     L     ]
        ISIZE = ISIZE + NFBAS(IFRP,0)*NFBAS(IFRP,1)
      ENDDO
      ISIZE = ISIZE*NZ
C  ... allocate array for the picture change transformation matrix
      CALL MEMGET('REAL',KPCTM,ISIZE,WORK,KFREE,LFREE)
C     ... array for the picture change transformed operator
      CALL MEMGET('REAL',KPCTO,N2BBASXQ,WORK,KFREE,LFREE)

      CALL DOPCTRA(VMAT,WORK(KPCTM),WORK(KPCTO),ISIZE,.TRUE.,
     &             'BNCR_IOTC       ',1,0,1,'S',
     &              WORK(KFREE),LFREE,IPRHAM)

      CALL MEMREL('MAKEBNCR2C',WORK,KPCTM,KPCTM,KFREE,LFREE)
      IF(IPRHAM.GE.8) THEN                                          
       CALL HEADER(
     & 'MAKEBNCR2C: IOTC-2c matrix screened bare nucleus',-1) 
       CALL PRQMAT(VMAT,NTBAS(0),NTBAS(0),                
     &    NTBAS(0),NTBAS(0),1,IPQTOQ(1,0),LUPRI) 
      END IF 
C     ... save the pctra screened bare nucleus into the BSSMAT file ...
      CALL W2BSSMAT(1,.FALSE.,'BNCR2C_4','PREPBNCR',
     &               VMAT,N2BBASXQ,LUBSS,IPRHAM)

      CALL QEXIT('MAKEBNCR2C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck rfbssmat */
      SUBROUTINE RFBSSMAT(LBL,KEEP_OPEN,FROMSUB,
     &                    TMAT,ISIZE,LUBSS,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE: Read the LBL record label from the BSSMAT file into
C           the entering TMAT array
C
C  On input: LBL - label of the record to be read
C            KEEP_OPEN - flag for keeping the BSSMAT file open after
C                        reading
C            FROMSUB -  name of the routine from which RFBSSMAT is called
C            TMAT - array for the record
C            ISIZE - the size of the record
C            LUBSS - file unit number ...
C
C  On output: TMAT with record
C
C  Written by Miro ILIAS, Strasbourg, March 2006
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      DIMENSION TMAT(*)
      LOGICAL FNDLAB,KEEP_OPEN,IS_HERE,IS_OPEN
      CHARACTER*8 LBL
      CHARACTER*(*) FROMSUB
      
      IF (IPRINT.GE.3) THEN
       CALL HEADER('*** Output from the RFBSSMAT routine ***',-1)
       write(lupri,'(2x,a,a)')
     &            '* RFBSSMAT called from routine:',FROMSUB
       write(lupri,'(2x,a,a)')
     &            '* The label to be read from BSSMAT:',LBL
       write(lupri,'(2x,a,i8)')
     &         '* Length of the record on the file:',ISIZE
      ENDIF

      INQUIRE(FILE='BSSMAT',EXIST=IS_HERE,OPENED=IS_OPEN)

      IF (.NOT.IS_HERE) THEN
       WRITE(LUPRI,'(/,2x,a,a,a)') 'RFBSSMAT from ',FROMSUB,
     &   ': BSSMAT file is not present !'
       CALL QUIT('RFBSSMAT: BSSMAT file is not present !')
      ENDIF

      IF (.NOT.IS_OPEN) THEN
       CALL OPNFIL(LUBSS,'BSSMAT','OLD',FROMSUB)
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)') '* Opening the BSSMAT file.'
       ENDIF
      ELSE
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)') '* BSSMAT file was already open. '//
     &                         ' No opening.'
       ENDIF
      ENDIF

      REWIND LUBSS
      IF (FNDLAB(LBL,LUBSS)) THEN
        CALL READT(LUBSS,ISIZE,TMAT)
        IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(2x,a,a,a)')
     &   '* RFBSSMAT: ',LBL,' record  was read from the file BSSMAT'
        ENDIF
      ELSE
         WRITE(LUPRI,'(2x,a,a,a,a)')
     & 'RFBSSMAT from ',FROMSUB,': ',LBL,' record was'//
     & ' NOT FOUND in the file BSSMAT !'
       CALL QUIT('RFBSSMAT: Record NOT FOUND in the file BSSMAT!') 
      ENDIF

      IF (.NOT.KEEP_OPEN) THEN
        CLOSE(LUBSS,STATUS='KEEP')
        IF (IPRINT.GE.3) THEN
         write(lupri,'(2x,a)')  
     &   '* BSSMAT file left closed.'
        ENDIF
      ELSE
        IF (IPRINT.GE.3) THEN
         write(lupri,'(2x,a)')
     &   '* BSSMAT file left open.'
        ENDIF
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck rfbssmat */
      SUBROUTINE W2BSSMAT(IW,KEEP_OPEN,LBL,
     &           FROMSUB,TMAT,ISIZE,LUBSS,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE: WRITE the LBL record (stored in TMAT) to the BSSMAT file 
C           acording to IW, KEEP_OPEN parameters.
C
C  On input:
C -----------
C            IW - integer number (1-4) specifying how to write the record
C
C              IW = 1: find EOFLABEL, write record, renew EOFLABEL
C              IW = 2: find EOFLABEL, write record
C              IW = 3: write record, renew EOFLABEL
C              IW = 4: write record only
C              IW = 5: empty
C              IW = 6: only renew the EOFLABEL
C
C            KEEP_OPEN - if .true., the BSSMAT file is left open
C            LBL - the label of the record
C            FROMSUB - the name of the routine from where W2BSSMAT is called
C            TMAT - data to be written into BSSMAT
C            ISIZE - the size of the record
C            LUBSS - file unit number.
C
C  On output: none
C
C   Written by Miro ILIAS, Strasbourg March, 2006
C   Latest modificatins: MI, Prievidza, Jun 2006
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"

      DIMENSION TMAT(*)
      LOGICAL FNDLAB,IS_HERE,IS_OPEN,KEEP_OPEN
      CHARACTER*8 LBL
      CHARACTER*(*) FROMSUB

      CALL QENTER('W2BSSMAT')

      IF (IPRINT.GE.3) THEN
        CALL HEADER('*** Output from the W2BSSMAT routine ***',-1)
        write(lupri,'(2x,a,a)')  
     &   '* W2BSSMAT was called from routine: ',FROMSUB
        write(lupri,'(2x,a,i2)') 
     &   '* Entering parameter on how to write IW=',IW
        write(lupri,'(2x,a,a)') 
     &   '* The entering label for record: ',LBL
        write(lupri,'(2x,a,i5)')
     &  '* Entering size of record (double words) ISIZE=',ISIZE
        write(lupri,'(2x,a,i3)')
     &  '* Entering logical unit of the file, LUBSS=',LUBSS
        write(lupri,'(2x,a,l6)')
     &  '* Keeping BSSMAT file open after writing the record (T/F):',
     &  KEEP_OPEN
        CALL FLSHFO(LUPRI)
      ENDIF

      IF (ISIZE.LE.0) THEN
       WRITE(LUPRI,*) 'W2BSSMAT from ',FROMSUB,
     &                ': Entering ISIZE=',ISIZE
       CALL QUIT('W2BSSMAT: wrong value of ISIZE ! <= 0 !')
      ENDIF

C     ... find out about the BSSMAT file
      INQUIRE(FILE='BSSMAT',EXIST=IS_HERE,OPENED=IS_OPEN)

      IF (.NOT.IS_OPEN.AND.IS_HERE) THEN
       CALL OPNFIL(LUBSS,'BSSMAT','OLD',FROMSUB)
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)')
     &  '* BSSMAT does exist, opening as an old file.'
       ENDIF
      ELSE IF (.NOT.IS_HERE.AND..NOT.IS_OPEN) THEN
       CALL OPNFIL(LUBSS,'BSSMAT','UNKNOWN',FROMSUB)
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)')
     &  '* BSSMAT does not exist, opening as an UNKNOWN file.'
       ENDIF
      ELSE IF (IS_HERE.AND.IS_OPEN) THEN
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)')
     &  '* BSSMAT file does exist and is created.'
       ENDIF
       CONTINUE
      ELSE
       WRITE(LUPRI,*) 'W2BSSMAT from ',FROMSUB,
     &                ': dead branch !'
       CALL QUIT('W2BSSMAT: dead branch, check routine !')
      ENDIF

      IF (IW.EQ.1.OR.IW.EQ.2) THEN
        REWIND LUBSS                                                    
        IF (FNDLAB('EOFLABEL',LUBSS)) THEN                            
          BACKSPACE(LUBSS)                   
          CALL NEWLAB(LBL,LUBSS,LUPRI)             
          CALL WRITT(LUBSS,ISIZE,TMAT)                   
          IF (IW.EQ.1) THEN
            CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
            IF (IPRINT.GE.3) THEN
             WRITE(LUPRI,'(2x,a,a,a)')
     &      '* Code 1:',LBL,' record '//
     &      'was WRITTEN to the file BSSMAT after the EOFLABEL,'//
     &      ' which was renewed.'
            ENDIF
          ELSE
           IF (IPRINT.GE.3) THEN
            WRITE(LUPRI,'(2x,a,a,a)')
     & '* Code 2:',LBL,' record was WRITTEN to the file BSSMAT'//
     & ' after the EOFLABEL, which was NOT renewed.'
           ENDIF
          ENDIF
        ELSE                                                       
         WRITE(LUPRI,'(2x,a,a,a)') 
     &     'W2BSSMAT called from ',FROMSUB,
     &     ': EOFLABEL not found on BSSMAT !' 
         CALL QUIT( 
     & 'W2BSSMAT: Error, EOFLABEL not found in BSSMAT file!')    
        ENDIF
      ELSE IF (IW.EQ.3.OR.IW.EQ.4.OR.IW.EQ.5)  THEN
        CALL NEWLAB(LBL,LUBSS,LUPRI)             
        CALL WRITT(LUBSS,ISIZE,TMAT)                   
        IF (IW.EQ.3) THEN
          CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
          IF (IPRINT.GE.3) THEN
            WRITE(LUPRI,'(2x,a,a)')
     &      '* Code 3:',LBL,' record '//
     &     'was WRITTEN to file BSSMAT. EOFLABEL was renewed.'
          ENDIF
        ELSE IF (IW.EQ.4.OR.IW.EQ.5) THEN
          IF (IPRINT.GE.3) THEN
            WRITE(LUPRI,'(2x,a,a,a)')
     &     '* Code 4:',LBL,' record '//
     &     'was WRITTEN to the file BSSMAT. EOFLABEL NOT renewed.'
           ENDIF
        ENDIF
      ELSE IF (IW.EQ.6) THEN
       CALL NEWLAB('EOFLABEL',LUBSS,LUPRI)
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)')
     &  '* Code 6: Only the EOFLABEL was added.'
       ENDIF
      ELSE
       WRITE(LUPRI,'(2x,a,a,i2)') 'W2BSSMAT called from ', FROMSUB,
     &                ': Entering code IW=',IW
       CALL QUIT('W2BSSMAT: Wrong value of entering IW ! Must be 1-5 !')
      ENDIF
 
C     ... close the file if not set to keep it open
      IF (.NOT.KEEP_OPEN) THEN
       CLOSE(LUBSS,STATUS='KEEP')
       IF (IPRINT.GE.3) THEN
        write(lupri,*)
     &  '* BSSMAT file was closed.'
       ENDIF
      ELSE
       IF (IPRINT.GE.3) THEN
        write(lupri,'(2x,a)')
     &  '* BSSMAT file was left open.'
       ENDIF
      ENDIF

      CALL FLSHFO(LUPRI)
      CALL QEXIT('W2BSSMAT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck isonbssmat */
      LOGICAL FUNCTION ISONBSSMAT(LBL,FROMSUB,LUBSS,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Written by Miro ILIAS, Strasbourg March, 2006,
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      LOGICAL FNDLAB,IS_HERE,IS_OPEN
      CHARACTER*8 LBL
      CHARACTER*(*) FROMSUB
      CALL QENTER('ISONBSSMAT')

      INQUIRE(FILE='BSSMAT',EXIST=IS_HERE,OPENED=IS_OPEN)

      IF (.NOT.IS_HERE) THEN
       write(lupri,*) 'ISONBSSMAT: BSSMAT file is not present !'
       write(lupri,*) 'ISONBSSMAT was called from ',FROMSUB
       CALL QUIT('ISONBSSMAT: BSSMAT file is not present !')
      ENDIF

      IF (.NOT.IS_OPEN) THEN
       CALL OPNFIL(LUBSS,'BSSMAT','OLD',FROMSUB)
      ENDIF

      REWIND LUBSS                                                    
      IF (FNDLAB(LBL,LUBSS)) THEN                            
         ISONBSSMAT = .TRUE.
      ELSE
         ISONBSSMAT = .FALSE.
      ENDIF
      CLOSE(LUBSS,STATUS='KEEP')

      CALL QEXIT('ISONBSSMAT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck ini_sf_rkb */
      SUBROUTINE INI_SF_RKB(IBEIG,IPRINT)
C****************************************************************************
C
C   Initialize the IBEIG(*) array with boson info for 2-component
C calculations using the RKB (Trond's) approach for the infinite order !
C
C  Written by Miro ILIAS, Strasbourg, 2005
C
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbgen.h"
#include "dgroup.h"
      DIMENSION IBEIG(*)
      LOGICAL FNDLAB

C     ... initialize NORB_SUB(IB,IFRP,0),ID_SUB_BL(IB,IFRP)
      CALL INISUB
C     ... run over fermion symmetries ...
      DO I = 1, NFSYM
       CALL INIBOS(IBEIG(1+IORB(I)),I,.TRUE.,IPRINT)
!      print *,'IBEIG(1+IORB(I)), I',IBEIG(1+IORB(I)),I
      ENDDO

C    ... write the boson irreps (Trond's approach) into the BSSMAT file !
      ISIZE = NORB(1)                                 
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)         
      CALL W2BSSMAT(1,.TRUE.,'IBOS_ONB','INI_SF_RKB',
     &              IBEIG,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
         WRITE(LUPRI,'(2X,A)')
     &   'INI_SF_RKB: IBEIG (IBOS_ONB) written to BSSMAT'
      ENDIF
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck TEST_UNIT_MTX */
      SUBROUTINE TEST_UNIT_MTX(S,N,LR,LC,NZ)
C****************************************************************************
C
C   PURPOSE: Check the enterning matrix - is it unit mtx ?
C
C  Written by Miro ILIAS, Odense,2003
C             
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, THRNULL = 1.00D-09 )

      DIMENSION S(LR,LC,NZ)

      SUM_DE  = SUM_DIAG_EL(S,N,LR,LC,NZ)
      SUM_DEN = SUM_DE/DFLOAT(N)
      SUM_ODE = SUM_OFFDIAG_EL(S,N,LR,LC,NZ)

      WRITE(LUPRI,'(/2X,A,D12.5)') 
     &'Sum of diagonal elements /N  :',SUM_DEN
      WRITE(LUPRI,'(2X,A,D12.5)') 
     &'Sum of off-diagonal elements :',SUM_ODE 

      IF (DABS(SUM_DEN-D1).GT.THRNULL) 
     & WRITE(LUPRI,'(2X,A)') 
     &'WARNING, SUM OF DIAG.ELEM/N DIFFERS TOO MUCH FROM 1!'

      IF (DABS(SUM_ODE-D0).GT.THRNULL) 
     & WRITE(LUPRI,'(2X,A)') 
     &'WARNING, SUM OF OFF_DIAG.ELEM DIFFERS TOO MUCH FROM 0!'
      CALL FLSHFO(LUPRI)  
      RETURN
      END

      REAL*8 FUNCTION SUM_DIAG_EL(S,N,LR,LC,NZ)
#include "implicit.h"
      PARAMETER (  D0 = 0.00D00 )
      DIMENSION S(LR*LC*NZ)

      PX = D0
      DO I = 1, N
       II = I + (I-1)*LR
       PX = PX + S(II)
      ENDDO

      SUM_DIAG_EL = PX

      RETURN
      END

      REAL*8 FUNCTION SUM_OFFDIAG_EL(S,N,LR,LC,NZ)
#include "implicit.h"
      PARAMETER ( D0 = 0.00D00 )
      DIMENSION S(LR*LC*NZ)

      PX = D0
      DO I = 2, N
      DO J = 1, I-1
       IJ = I + (J-1)*LR
       JI = J + (I-1)*LR
       PX = PX + S(IJ) + S(JI) 
      ENDDO
      ENDDO

      SUM_OFFDIAG_EL = PX

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck TEST_DIAG_MTX */
      SUBROUTINE TEST_DIAG_MTX(S,N,LR,LC,NZX)
C****************************************************************************
C
C  PURPOSE: Check the enterning matrix - is it a diagonal matrix ?
C
C  Diagonal matrix should have diagonal elements between 0 and 1
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C             
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00, THRNULL = 1.00D-09 )
      DIMENSION S(LR,LC,NZX)
      
      SUM_DE  = SUM_DIAG_EL(S,N,LR,LC,NZX)
      SUM_DEN = SUM_DE/DFLOAT(N)
      SUM_ODE = SUM_OFFDIAG_EL(S,N,LR,LC,NZX)

      WRITE(LUPRI,'(/2X,A,D12.5)') 
     &'Sum of diagonal elements /N :',SUM_DEN
      WRITE(LUPRI,'(2X,A,D12.5)') 
     &'Sum of off-diagonal elements :',SUM_ODE 
      CALL FLSHFO(LUPRI)

      IF (DABS(SUM_DEN-D1).GT.THRNULL) THEN
        WRITE(LUPRI,'(2X,A)') 
     & 'TEST_DIAG_MTX: WARNING, sum of diagonal elements/N'//
     &' of a unit matrix differs too much from 1!'
        CALL FLSHFO(LUPRI)
      ENDIF

      IF (DABS(SUM_ODE-D0).GT.THRNULL) THEN 
       WRITE(LUPRI,'(2X,A)') 
     &'TEST_DIAG_MTX: WARNING, sum of off-diagonal elements of a unit'//
     &' matrix differs too much from 0!'
       CALL FLSHFO(LUPRI)
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck EXTR_BSS_INFO */
      SUBROUTINE EXTR_BSS_INFO(IBSS,IPZ,SPINFR,DK2SO1,SAVESO1,ADD_AMFI,
     &                    ONLYSO1,ONLYTRA,SF0,ELIM_SF,
     &                    NOFPB,BLOCKD,DO4C2C,USEDF,CONT2C,
     &                    START2C,I2CHAM,
     &                    DO2C4C,INI2C,I2COFK,
     &                    YREQ1,IC_AMFI,GETPCTM,NOAMFI,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C PURPOSE:
C ========
C     Set control variables based on IBSS coding
C     IBSS = axyz
C
C        a - 1 to 6, defines AMFI flags.. - if NOAMFI=.true., zero the "a" !
C
C        x - .gt. 0: define spinfree before (as 4-comp)
C            .eq. 0: define spinfree at end
C            ... x is active with the SPINFREE option & with the y=1 (means only SO1)
C           
C        y - order of spin-orbit term:
C            0 - spinfree, 1 - first order for DK2/BSS, 2+ - full.
C
C        z - order of Douglas-Kroll transformation
C            (0, 1, 2, or 9), where 9 is full(infinite).
C
C     Pick up the x,y,z variables from IBSS = xyz
C
C     Note: in the case of one-step method xyz must be set up for infinite order !
C
C On input: IBSS, SPINFR, IPRINT, ....
C
C On output: IPZ - order of Hamiltonian
C 
C            DK2SO1,SF0
C
C            ONLYSO1 - flagfo BSS-DKH-SO1 integrals to be compared with AMFI DKH-SO1 integrals
C
C            SAVESO1  - flag to decide to save SO1 elements (from H1) into the BSSMAT file
C
C  to be added at the later stage to the final spin-free 2c Hamilt.(BSSsfB,BSSsfE,DKH2sfB,DKH2sfE)
C
C            ONLYTRA - flag for C->S->C transformation of H2c for test
C                       of AMFI routine working in spherical basis
C
C            ELIM_SF
C
C            IC_AMFI - flag for AMFI subroutine
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
C
      LOGICAL SPINFR,DK2SO1,SF0,ONLYSO1,ONLYTRA,ADD_AMFI,ELIM_SF,
     &        NOFPB,BLOCKD,YREQ1,DO4C2C,START2C,DO2C4C,
     &        USEDF,SAVESO1,GETPCTM,CONT2C,NOAMFI

      CALL QENTER('EXTR_BSS_INFO')

      CALL HEADER('*** Output from the EXTR_BSS_INFO routine ***',-1)

C     ... do extract a,x,y,z digits of the entering IBSS=axyz value
      CALL EXTR_DIG(IBSS,IPA,IPX,IPY,IPZ,IPRINT)

      IF (DO2C4C.AND.START2C) THEN
       WRITE(LUPRI,'(/2X,A,I2/)')
     & 'Before the four-component Dirac-Fock'//
     & '  procedure run preliminary'//
     & ' TWO-COMPONENT (infinite order) SCF ! code INI2C=',INI2C
      ENDIF

      IF (DO4C2C) THEN
       WRITE(LUPRI,'(2X,A)')
     & '* After the Dirac-Fock procedure continue'//
     & ' in the TWO-COMPONENT (infinite order) MODE !'
       IF (USEDF) THEN
         WRITE(LUPRI,'(/,2X,A/)')
     & '...transforming the converged four-component'//
     & ' DIRAC-FOCK operator !'
       ENDIF
       IF (CONT2C) THEN
         WRITE(LUPRI,'(/2X,A,I1)')
     & '...continuing in two-component (infinite order) SCF'//
     & ' procedure ... I2CHAM (3,4,5,6) =',I2CHAM
       ENDIF
      ELSE
        WRITE(LUPRI,'(2X,A)')
     & '* Transforming the one-electron DIRAC bare nucleus operator.'
      ENDIF

      WRITE (LUPRI,'(2X,A,$)')
     & '* Order of the Douglas-Kroll-Hess/BSS transformation :'
      IF (IPZ.EQ.9) THEN
C      write(lupri,'(5X,A)') '...infinite  !'
       write(lupri,'(A)') '...infinite  !'
      ELSE IF (IPZ.EQ.2) THEN
C      write(lupri,'(5X,A)') '...second (DKH2) !'
       write(lupri,'(A)') '...second (DKH2) !'
      ELSE IF (IPZ.EQ.1) THEN
C      write(lupri,'(5X,A)') '...first (DKH1) !'
       write(lupri,'(A)') '...first (DKH1) !'
      ELSE IF (IPZ.EQ.0) THEN
       write(lupri,'(/5X,A)')
     & '...only testing of the SO1 approach, then quit ...'
      ELSE
       write(lupri,*) 'IPZ=',IPZ
       write(lupri,*) 'Can be only 9,2,1 or 0!'
       CALL QUIT('MAKE_BSS: EXTR_BSS_INFO: Wrong value of'//
     &           ' IZ (should be 0,1,2,9) !')
      ENDIF

C  ... check if for some finite-order we have the free-particle basis transformation
      IF (IPZ.NE.9.AND.NOFPB) THEN
        WRITE(LUPRI,'(/,2X,A,A,I3)')
     &  'Wrong BSS input: One-step method and NOT the infinite order !',
     &   ' IPZ=',IPZ
        CALL FLSHFO(LUPRI)
        CALL QUIT('Wrong BSS input: One-step method and NOT'
     &  //' the infinite order !')
      ENDIF
      
C      ... treat the NOFPB & BLOCKD
      IF (IPZ.EQ.9) THEN
       IF (NOFPB) THEN
        WRITE(LUPRI,'(2X,A)')
     &  '* Only the RKB basis (one-step method). No preliminary free-'//
     &  'particle basis transformation according to Barysz.'
       ELSE
        WRITE(LUPRI,'(2X,A)')
     & '* Applied the preliminary free-particle ("k") basis'//
     & ' transformation (two-step method). '
       ENDIF
C      ... by default do not use pure numerical block-diagonalizating approach
C       this is still TODO ...
       IF (BLOCKD) THEN
        WRITE(LUPRI,'(2X,A)')
     &  '* The block diagonalization method (QJACOBI) applied !'
       ELSE
        WRITE(LUPRI,'(2X,A,$)')
     &  '* R-equation solving applied :'
        IF (YREQ1) THEN
          WRITE(LUPRI,'(A)') ' YR-equation 1 (HJAaJ) '
        ELSE
          WRITE(LUPRI,'(A)') ' YR-equation 2 (Trond) '
        ENDIF
       ENDIF
      ELSE
C     ... for lower orders do always use the FP "k"- basis !
         NOFPB = .FALSE.
         BLOCKD = .FALSE.
      ENDIF

C      ... remove AMFI contribution when .NOAMFI present
      IF (IPA.GT.0.AND.NOAMFI) THEN
        WRITE(LUPRI,'(2X,A)') 
     &  '* Hard removal of the AMFI contribution due to .NOAMFI flag'
        IPA = 0
      ENDIF

      ELIM_SF=.FALSE.
      IF (SPINFR.AND.IPY.GE.1) THEN 
       WRITE(LUPRI,'(2X,A)') '* "SPINFR" option, no spin-orbit terms!' 
C      ... to secure that no spin-orbit terms appear
       IPY=0
      ENDIF

      DK2SO1 = .FALSE.
         SF0 = .FALSE.
      IF (IPY .EQ. 0) THEN
         IF (.NOT.SPINFR.AND.IPA.EQ.0) CALL QUIT(
     & 'You forgot to specify .SPINFR in the Hamiltonian input!')
         DK2SO1 = .FALSE.
         IF (SPINFR.AND.IPA.GT.0) THEN 
          WRITE(LUPRI,'(/,2x,a)') 'EXTR_BSS_INFO: Wrong .SPINFR in the'
     &    //' Hamiltonian input with the AMFI ... discard it !'
          CALL QUIT(
     &    'Wrong .SPINFR in the Hamiltonian input with the AMFI on!')
         ENDIF

         IF (.NOT.SPINFR.AND.IPA.GT.0) THEN
          SPINFR=.TRUE.
          ELIM_SF=.TRUE.
         ENDIF
      ELSE IF (IPY.EQ.1.AND.IPZ.EQ.2.AND.IPX.EQ.1
     &           .AND.(.NOT.SPINFR)) THEN
C        ... DKH2sfB + SO1
        DK2SO1 = .TRUE.
        WRITE (LUPRI,'(2X,A)')
     &  '* spin-orbit terms only to the first order -'//
     &  ' (within the DKH2 Hamilt. generation)'//
     &  ' - DKH2sfBSO1 Hamiltonian'
      ELSE IF (.NOT.SPINFR) THEN
         DK2SO1 = .FALSE.
         WRITE (LUPRI,'(2X,A)')
     &   '* spin-orbit terms are included'
      END IF

      IF (IPX .GE. 1 .AND. SPINFR) THEN  
         SF0 = .TRUE.
         WRITE (LUPRI,'(2X,A)')
     & '* Elimination of spin-orbit terms from the Dirac bare '//
     & 'nuclei in RKB/"k" basis (spinfree "BEFORE / AT THE BEGINNING")' 
      ELSE IF(IPX .EQ. 0 .AND. SPINFR) THEN
         SF0 = .FALSE.
         IF (SPINFR) WRITE (LUPRI,'(2X,A)')
     & '* Elimination of spin-orbit terms from the 2comp.Ham. '//
     & ' in RKB/"k" basis (spinfree "AFTER / AT THE END")'
      ENDIF

C     ... by default do not include SO1 addition to spin-free Hamiltonians DK2,BSS
      SAVESO1 = .FALSE.
      IF (IPY.EQ.1.AND..NOT.SPINFR.AND..NOT.DK2SO1
     &            .AND.IPZ.GT.1.AND..NOT.NOFPB) THEN
        SAVESO1=.TRUE.
        SPINFR=.TRUE.
        SF0=.FALSE.
        IF (IPX.GE.1) SF0=.TRUE.
        WRITE(LUPRI,'(/2X,A/)')
     & ' Addition of SO1 terms (extracted from H1_ONB)'//
     & ' to the desired spin-free Hamiltonian !'
        IF (IPZ.EQ.2.AND.SF0) THEN
          write(lupri,'(5X,A/)') 'DK2sfB + SO1 Hamiltonian'
        ELSE IF (IPZ.EQ.2.AND..NOT.SF0) THEN
          write(lupri,'(5X,A/)') 'DK2sfE + SO1 Hamiltonian'
        ELSE IF (IPZ.GT.2.AND.SF0) THEN
          write(lupri,'(5X,A/)') 'BSSsfB + SO1 Hamiltonian'
        ELSE IF (IPZ.GT.2.AND..NOT.SF0) THEN
          write(lupri,'(5X,A/)') 'BSSsfE + SO1 Hamiltonian'
        ENDIF
      ENDIF

      IF (ELIM_SF) WRITE (LUPRI,'(A)')                                
     & '* SPIN-FREE flag is valid only for the generation of the BSS'//
     & ' Hamiltonain, later removed '//
     & ' due to the AMFI spin-orbit integrals.'
     
C      .... flag IBSS=010 - activate the ONLYSO1 flag
      ONLYSO1=.FALSE.
      ADD_AMFI=.FALSE.
CMI   IF (IPX.EQ.0.AND.IPY.EQ.1.AND.IPZ.EQ.0) THEN
      IF (IPX.EQ.0.AND.IPY.EQ.1.AND.IPZ.EQ.0.AND..NOT.NOAMFI) THEN
         WRITE (LUPRI,'(/2X,A)')
     &  '* only SO1 terms (for comparison with AMFI, then quit)'
        ONLYSO1=.TRUE.
C       ... AMFI gives only SO1 integrals
        IC_AMFI = 1
        ADD_AMFI=.TRUE.
      ENDIF

      ONLYTRA=.FALSE.
      IF (IPA.EQ.1) THEN
        ONLYTRA=.TRUE.
        ADD_AMFI=.TRUE.
C     ... AMFI calculated SO1+MFSO2 integrals,
C        but these are added to the Hamiltionian with ZERO weight !
        IC_AMFI = 3
        WRITE(LUPRI,'(/2X,A)')
     &'* C->S->C transformation + AMFI run(test), but no adding'//
     &' of AMFI integrals.'  
      ELSE IF (IPA.EQ.2) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 2
        WRITE(LUPRI,'(/2X,A)')
     &'* AMFI: MFSSO2 terms (spin-same orbit only) added.'  
      ELSE IF (IPA.EQ.3) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 3
        WRITE(LUPRI,'(/2X,A)')
     &'* AMFI: MFSO2 terms (both spin-same&spin-other orbit) added.'  
      ELSE IF (IPA.EQ.4) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 4
        WRITE(LUPRI,'(/2X,A)')
     &'* AMFI: SO1+MFSSO2 terms (only spin-same orbit) added.'
      ELSE IF (IPA.EQ.5) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 5
        WRITE(LUPRI,'(/2X,A)')
     &'* AMFI: SO1+MFSO2 terms (spin-same&spin-other orbit) added.'  
      ELSE IF (IPA.EQ.6) THEN
        ADD_AMFI=.TRUE.
        IC_AMFI = 1
        WRITE(LUPRI,'(/2X,A)')
     &'* AMFI: SO1 terms added.'  
      ENDIF

      IF (.NOT.ADD_AMFI) THEN
        WRITE (LUPRI,'(2X,A)')
     &  '* AMFI contribution is NOT included.'
      ENDIF

C ...   treat the GETPCTM flag - for the infinite order always prepare
C     the picture change transformation matrix !
      GETPCTM = .FALSE.
C     IF (.NOT.SPINFR.AND..NOT.SAVESO1.AND.IPZ.EQ.9) THEN
      IF (.NOT.SAVESO1.AND.IPZ.EQ.9) THEN
       GETPCTM = .TRUE. 
       WRITE(LUPRI,'(2X,A)')
     & '* Prepare the picture-change transformation matrixes'
       IF (SPINFR) WRITE(LUPRI,'(4X,A)')
     & '  - special case - spin-free  '
      ELSE
       WRITE(LUPRI,'(2X,A)')
     & '* Do NOT prepare the picture-change transformation matrix !'
      ENDIF

C     ... consistency checks - do not repeat 1-el SO interactions !!!
C =========================================================================
      IF (IPY.GT.0.AND.(IPA.EQ.6.OR.IPA.EQ.5.OR.IPA.EQ.4)) THEN
       write(lupri,*)
       WRITE(LUPRI,*) 'EXTR_BSS_INFO: IPY=',IPY,' IPA=',IPA
       WRITE(LUPRI,*) 'EXTR_BSS_INFO: Doubled inclusion of 1el'//
     & ' SO terms - both from AMFI & from BSS Hamilatonian !!! '
       CALL  QUIT('EXTR_BSS_INFO: Double inclusion of 1el SO terms! ')
      ENDIF

C ... START2C preliminary case - verify other settings
      IF ( (START2C.AND.DO2C4C) .OR. DO4C2C) THEN
C       WRITE(LUPRI,'(/2X,A/)')
C    &  'Starting the preliminary 2comp. BSS-SCF procedure.'
C     ... can be full 999 or 999 and SPINFR
        IF (.NOT.GETPCTM.OR.(GETPCTM.AND.SPINFR.AND..NOT.SF0)) THEN
         write(lupri,*) 'EXTR_BSS_INFO: GETPCTM=',GETPCTM,
     &   ' SPINFR=',SPINFR,' SF0=',SF0
         CALL QUIT(
     &   'EXTR_BSS_INFO: Preliminary BSS-SCF and/or post-DC-BSS-SCF'//
     &   ' is not possible !')
        ENDIF
      ENDIF
C  ... set up I2COFK parameter for the ONEFCK according to the I2CHAM
      I2COFK = 1
      IF (ADD_AMFI) I2COFK=2
      IF (DO4C2C.AND..NOT.START2C) THEN
C     ... pure 2c mode after the DC-SCF ...
       IF (I2CHAM.EQ.3) THEN
         I2COFK = 3
       ELSE IF (I2CHAM.EQ.4) THEN
         I2COFK = 4
       ELSE IF (I2CHAM.EQ.5) THEN
         I2COFK = 5
       ELSE IF (I2CHAM.EQ.6) THEN
         I2COFK = 6
       ELSE
        CONTINUE
!         CALL QUIT(
!     &   'EXTR_BSS_INFO: I2CHAM parameter '//
!     &   ' is wrong (not 3/4/5/6) !' )
       ENDIF
      ENDIF

      WRITE(LUPRI,'(2X,A,I1)')
     &'* I2COFK parameter for H2c specification set to:',I2COFK

      CALL FLSHFO(LUPRI)
      CALL QEXIT('EXTR_BSS_INFO')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck EXTR_DIG */
      SUBROUTINE EXTR_DIG(IBSS,IPA,IPX,IPY,IPZ,IPRINT)
C****************************************************************************
C
C  PURPOSE: For entering integer number of type xyz,
C            extract a,x,y,z into IPA,IPX,IPY,IPZ
C             
C  Written: M Ilias, Odense, aug. 2002             
C
C****************************************************************************
#include "implicit.h"
#include "priunit.h"

      IF (IBSS.LT.1.OR.IBSS.GT.9999)
     & CALL QUIT('EXTR_DIG: The entering number not in range 1-9999! ')
  
C     ... save the IBSS
      IBSS_SAVE = IBSS
      IBSS0     = IBSS

      IF (IBSS.GE.1000) THEN
       IPA  = IBSS/1000
       IBSS0 = MOD(IBSS,1000)
       IBSS = MOD(IBSS0,100)
       IPX  = IBSS0/100
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.100) THEN 
       IBSS = MOD(IBSS,100)
       IPA = 0
       IPX  = IBSS0/100
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.10) THEN
       IPA = 0
       IPX = 0
       IPY = IBSS/10
       IPZ = MOD(IBSS,10)
      ELSE IF (IBSS.GE.0) THEN
       IPA = 0
       IPX = 0
       IPY = 0
       IPZ = IBSS
      ENDIF

      IF (IPRINT.GE.3) THEN
       WRITE (LUPRI,'(/,4X,A,I4,A,4I3)')
     &'EXTR_DIG: IBSS=axyz=',IBSS_SAVE,
     &' ...extracted digits  a,x,y,z:',IPA,IPX,IPY,IPZ
      ENDIF

C     ... restore the IBSS
      IBSS = IBSS_SAVE
C     ... check the extraction procedure
      IBSS1 = (1000*IPA)+(100*IPX)+(10*IPY)+IPZ
      IF (IBSS1.NE.IBSS) THEN 
        WRITE(LUPRI,*) 'EXTR_DIG: input IBSS=axyz',IBSS
        WRITE(LUPRI,*) 'EXTR_DIG: extracted digits:',IPA,IPX,IPY,IPZ
        CALL QUIT('EXTR_DIG: Error in extraction of a,x,y,z from axyz!')
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck spf_hk */
      SUBROUTINE SPF_HK(HK,SPFR,IBS1,IBS2,
     &                 NROW,NCOL,LROW,LCOL,NZ,IPRINT)
C****************************************************************************
C
C  PURPOSE: Eliminate spin-orbit off-diagonal elements of entering matrix
C           in MO basis (SPFR=.true.) or eliminate spin-free elements(SPFR=.false.) 
C
C  On input: HK - quaternion hermitian matrix
C            IBS1 - boson symetries of first index
C            IBS2 - boson symetries of second index
C            SPFR - if true, elim.SO terms, otherwise elim.SF terms
C
C  On output: HK - the same matrix without spin-orbit/spin-free components
C             
C  Written: M Ilias, Odense, aug. 2002             
C           last update MI,2005,Strasbourg             
C
C****************************************************************************
#include "implicit.h"
      PARAMETER ( D0 = 0.00D00 )
#include "priunit.h"

      DIMENSION HK(LROW,LCOL,NZ),IBS1(*),IBS2(*)
      LOGICAL SPFR

      CALL QENTER('SPF_HK')

      IF (NROW.NE.NCOL) THEN
        WRITE(LUPRI,'(2X,A)') 'SPF_HK: Warning, NROW<>NCOL!'
        CALL QUIT('SPF_HK: NROW<>NCOL! How it is possible?')
      ENDIF

      IF (IPRINT. GE. 5) THEN
       write(lupri,*)
       IF (SPFR) THEN
        write(lupri,*)
     &  '..eliminating SPIN-ORBIT TERMS..'
       ELSE
        write(lupri,*)
     &  '..eliminating SPIN-FREE TERMS..'
       ENDIF
       WRITE(LUPRI,'(2X,A)')
     & 'SPF_HK: Entering RKB/"|k>" ON basis boson symmetries:'
       DO I = 1, NROW
        WRITE(LUPRI,'(5X,A,I3,3X,A,I2)')
!     &  'orbital: ',I,'boson symmetry:',IBS1(I)+1
     &  'orbital: ',I,'boson symmetry:',IBS1(I)
       ENDDO
      ENDIF

      IF (SPFR) THEN
C     ... eliminate spin-orbit components...
       IZ = 1
       DO I=1,NROW
       DO J=1,NCOL
        IF (IBS1(I).NE.IBS2(J)) HK(I,J,IZ)=D0
       ENDDO
       ENDDO

       IF (NZ.GT.1) THEN
        DO IZ = 2, NZ
        DO I=1,NROW
        DO J=1,NCOL
         HK(I,J,IZ)=D0
        ENDDO
        ENDDO
        ENDDO
       ENDIF

      ELSE
C     ... eliminate spin-free components
       IZ = 1
       DO I=1,NROW
       DO J=1,NCOL
         IF (IBS1(I).EQ.IBS2(J)) HK(I,J,IZ)=D0
       ENDDO
       ENDDO

      ENDIF

      CALL QEXIT('SPF_HK')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck dopctra */
      SUBROUTINE DOPCTRA(PRPAO,PCTM,PRPAO1,ISIZEPCTM,MAKEPCTRA,
     &          PROPNAME,ISYM,IREP,ITIM,TREV,WORK,LWORK,IPRINT)
C****************************************************************************
C
C  PURPOSE: Do picture change transformation of the entering AO 
C           one-electron matrix of the FOUR-COMPONENT property operator.
C
C           If MAKEPCTRA = .false., do not perform the picture change tranf.,
C           only take the LL block of a four-component property 
C           operator (this is not rigorous, but shows how important is
C           the picture change transformation of property operator)
C
C  Called from: 
C ==============
C    MAKEBNCR2C when picture change of screened bare nucleus is needed
C    PRPMSAO  when one needs to perform the picture change transformation of 
C             the entering property operator
C
C  On input: PRPAO - 4component operator in AO basis
C            PRPAO1 - array for 2comp. operator
C            INDXPR - index of the property operator
C            MAKEPCTRA - perform or not picture change transf.
C            to be transformed (if -1, th's diagonal SCF start from BNCORR)
C            IREP, ITIM - boson and time-rev.symmetry of the operator
C            TREV - symmetric (S) or antisymmetric (A) operator
C
C  On output: PRPAO - picture change transformed operator 
C             matrix (LL block) of AO basis, or only LL block
C
C  Written: M Ilias, Prievidza, June 2004
C  Last changes: MI, Strasbourg, October 2005    
C                MI, Prievidza, June 2006, August 2006  
C      
C****************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D0 = 0.00D00, D1 = 1.00D00 )
#include "dcbgen.h"
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION WORK(LWORK),PRPAO(N2BBASX,NZ),
     &        PRPAO1(N2BBASX,NZ),PCTM(ISIZEPCTM)
      LOGICAL MAKEPCTRA
      CHARACTER PROPNAME*16,TREV*1

      CALL QENTER('DOPCTRA')
#include "memint.h"

      IF (IPRINT.GE.5) THEN
        CALL HEADER('Output from DOPCTRA',-1)
        write(lupri,'(2X,A,A)')  '         Operator name: ',PROPNAME 
        write(lupri,'(2X,A,I2)') 'time reversal symmetry: ',ITIM
        write(lupri,'(2X,A,I2)') '      fermion rep IREP: ',IREP
        write(lupri,'(2X,A,I2)') '   boson symmetry ISYM: ',ISYM
        write(lupri,*)
      ENDIF

      CALL MEMGET('REAL',Kxxxx,1,WORK,KFREE,LFREE)
C  ... keep only LL blocks of the four-component property operator
      IF (.NOT.MAKEPCTRA) THEN
        CALL ONLY_LL(PRPAO,PROPNAME,IREP,ITIM,IPRINT)
      ELSE
C ...  read (load) the picture change transformation matrix:
       CALL RFBSSMAT('U_PICTRM',.FALSE.,'DOPCTRA',PCTM,
     &               ISIZEPCTM,LUBSS,IPRINT)

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/A,A16,/,10X,A,I2,A,I1,A,I1,A,A1)')
     &'DOPCTRA: Total entering FOUR-COMPONENT SAAO-matrix of property ',
     &PROPNAME,' time rev.symm ITIM=',ITIM,' IREP=',IREP,
     &' number of matrixes NZ=',NZ,' TREV=',TREV
        CALL PRQMAT(PRPAO(1,1),NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,IREP),LUPRI)
        CALL FLSHFO(LUPRI)
      ENDIF

CMI... this cleaning is necessary in order to keep PRPAO1 without
C     waste....
      CALL DZERO(PRPAO1,N2BBASXQ)

      IOPSY = JBTOF(IREP,1)
      DO I1 = 1, NFSYM
         I2 = MOD(I1+IOPSY,2) + 1
C ----------------------------------------------------------
C     Do the picture change transformation X_2c =U+.X_4c.U
C   NOTE: It works only as the symmetric transformation !!!
C ----------------------------------------------------------
        IF (IPRINT.GE.6) THEN
         WRITE(LUPRI,'(/2X,A,I1,I2,A,I1)')
     &   'DOPCTRA: Fermion coreps (left,right/total) I1,I2 / NFSYM: ',
     &   I1,I2,'/',NFSYM
         WRITE(LUPRI,'(2X,A,I3,A,I3,A,I3,A,I3)')
     &   '...sizes of matrixes before and after the picture'//
     &   ' change transformation: ',
     &   NFBAS(I1,0),' x ',NFBAS(I2,0),' -->',
     &   NFBAS(I1,1),' x ',NFBAS(I2,1)
         CALL FLSHFO(LUPRI)
        ENDIF

        IF (NFBAS(I1,1).GT.0.AND.NFBAS(I2,1).GT.0) THEN

C       ... find offsets for picture change transformation matrix
        IF (I1.EQ.1) THEN
           IPCTM1=1
        ELSE IF (I1.EQ.2) THEN
           IPCTM1 = 1 + NFBAS(1,0)*NFBAS(1,1)*NZ
        ENDIF

        IF (I2.EQ.1) THEN
           IPCTM2=1
        ELSE IF (I2.EQ.2) THEN
           IPCTM2 = 1 + NFBAS(1,0)*NFBAS(1,1)*NZ
        ENDIF

         CALL QTRANS('AOMO',TREV,D0,
     &     NFBAS(I1,0),NFBAS(I2,0),
     &     NFBAS(I1,1),NFBAS(I2,1),
     &     PRPAO(1+I2BASX(I1,I2),1),
     &     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
     &     PRPAO1(1+I2BASX(I1,I2),1),
     &     NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),
     &     PCTM(IPCTM1),
     &     NFBAS(I1,0),NFBAS(I1,1),NZ,IPQTOQ(1,0),
     &     PCTM(IPCTM2),
     &     NFBAS(I2,0),NFBAS(I2,1),NZ,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)

C .... control output
       IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/A,1X,A16)')
     & 'DOPCTRA: Symm. blocked FOUR-COMPONENT AO-matrix of property',
     &  PROPNAME

        WRITE(LUPRI,'(3X,A,I1,I2,A,I1)')
     &  '*** Fermion coreps (left,right/total) I1,I2 / NFSYM: ',
     &  I1,I2,'/',NFSYM

        CALL PRQMAT(PRPAO(1+I2BASX(I1,I2),1),NFBAS(I1,0),NFBAS(I2,0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)

         WRITE(LUPRI,'(/A)')
     & 'DOPCTRA:SYMMETRY BLOCKED picture change transformation matrix: '
         WRITE(LUPRI,'(3X,A,I1,I2,A,I1)')
     &      '*** Fermion coreps I1,I2/NFSYM ',I1,I2,'/',NFSYM
         IF (I1.EQ.I2) THEN
          write(lupri,*) 'acting on left=right'
          CALL PRQMAT(PCTM(IPCTM1),NFBAS(I1,0),NFBAS(I1,1),
     &        NFBAS(I1,0),NFBAS(I1,1),NZ,IPQTOQ(1,0),LUPRI)
         ELSE
          write(lupri,*) 'matrix on left'
          CALL PRQMAT(PCTM(IPCTM1),NFBAS(I1,0),NFBAS(I1,1),
     &        NFBAS(I1,0),NFBAS(I1,1),NZ,IPQTOQ(1,0),LUPRI)
          write(lupri,*) 'matrix on right'
          CALL PRQMAT(PCTM(IPCTM2),NFBAS(I2,0),NFBAS(I2,1),
     &        NFBAS(I2,0),NFBAS(I2,1),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         WRITE(LUPRI,'(/A,1X,A16)')
     & 'DOPCTRA: Symm. blocked picture change TRANSFORMED'//
     & ' TWO-COMPONENT AO-matrix of property',
     &  PROPNAME
         WRITE(LUPRI,'(5X,A,I1,I2,A,I1)')
     &    '*** Fermion coreps (left,right/total): ',I1,I2,'/',NFSYM
         CALL PRQMAT(PRPAO1(1+I2BASX(I1,I2),1),NFBAS(I1,1),NFBAS(I2,1),
     &             NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
         CALL FLSHFO(LUPRI)

         ENDIF      
        ENDIF      
      ENDDO 

      IF (IPRINT.GE.5) THEN
       WRITE(LUPRI,'(/A,1X,A16)')
     & 'DOPCTRA: Total picture change TRANSFORMED TWO-COMPONENT'//
     & ' AO-matrix of transformed property', PROPNAME
       WRITE(LUPRI,'(2X,A,I3,A,I2)')
     & '2-comp. matrix has to be in the upper LL block of size ',
     & NTBAS(1),' x ',NTBAS(1)
       CALL PRQMAT(PRPAO1,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &          IPQTOQ(1,IREP),LUPRI)
       CALL FLSHFO(LUPRI)
      ENDIF

      CALL DCOPY(N2BBASXQ,PRPAO1,1,PRPAO,1)

C ... check the pict change transformed operator...
      CALL CHECK_LLOP(IREP,IOPSY,PROPNAME,PRPAO,IPRINT)

      ENDIF
      CALL QEXIT('DOPCTRA')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck only_ll */
      SUBROUTINE ONLY_LL(PRPAO,PROPNAME,IREP,ITIM,IPRINT)
C**************************************************************************
C
C    PURPOSE: Keep only LL block of entering property matrix PRPAO.
C
C   Written by Miro ILIAS, Strasbourg, 2006
C   Last modifications:  MI, June 2006, Prievidza
C
C**************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D0 = 0.00D00, D1 = 1.00D00, DEB_THRS=1.0D-7 )
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION PRPAO(N2BBASX,NZ)
      CHARACTER PROPNAME*16

      CALL QENTER('ONLY_LL')

      IOPSY = JBTOF(IREP,1)
      DO I1 = 1, NFSYM
         I2 = MOD(I1+IOPSY,2) + 1

       IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/A,1X,A16,/,10X,A,I2,A,I1)')
     &'ONLY_LL: Symm. blocked FOUR-COMPONENT AO-matrix of property "',
     &  PROPNAME,'" time reversal symm:',ITIM,
     &  ' number of matrixes NZ: ',NZ
       WRITE(LUPRI,'(/5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I1,'/',NFSYM
       CALL PRQMAT(PRPAO(1+I2BASX(I1,I2),1),NFBAS(I1,0),NFBAS(I2,0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
       ENDIF

C      ... delet all other blocks (SL,LS,SS) except LL
       DO IZ = 1, NZ
C          ... first LS block
         ILS1 = 1 + I2BASX(I1,I2) + NFBAS(I1,1) 
         DO J = 1, NFBAS(I2,1)
           ILS = ILS1 + ((J-1)*NTBAS(0)) 
           CALL DZERO(PRPAO(ILS,IZ),NFBAS(I1,2))
         ENDDO
C          ... then SS,SL blocks
         ILS1 = 1 + I2BASX(I1,I2) + (NTBAS(0)*NFBAS(I2,1)) 
         DO J = 1, NFBAS(I2,2)
           ILS = ILS1 + ((J-1)*NTBAS(0)) 
           CALL DZERO(PRPAO(ILS,IZ),NFBAS(I1,0))
         ENDDO
       ENDDO

       IF (IPRINT.GE.5) THEN
        WRITE(LUPRI,'(/A,1X,A16,/,10X,A,I2,A,I1)')
     &'ONLY_LL: LL block of FOUR-COMPONENT AO-matrix of property',
     &  PROPNAME,' time reversal symm:',ITIM,
     &  ' number of matrixes NZ=',NZ
       WRITE(LUPRI,'(/5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I1,'/',NFSYM
       CALL PRQMAT(PRPAO(1+I2BASX(I1,I2),1),NFBAS(I1,0),NFBAS(I2,0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)
       ENDIF
      ENDDO

      CALL CHECK_LLOP(IREP,IOPSY,PROPNAME,PRPAO,IPRINT)

      CALL QEXIT('ONLY_LL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck CHECK_LLOP */
      SUBROUTINE CHECK_LLOP(IREP,IOPSY,PROPNAME,PRPAO,IPRINT)
C**************************************************************************
C
C   PURPOSE:  Check whether the (picture change transformed) operator 
C             matrix PRPAO contains only LL block
C
C    Written: M Ilias, Prievidza, June 2004
C    Last updates: MI, Strasbourg, Oct.2005
C
C**************************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D0 = 0.00D00, D1 = 1.00D00, DEB_THRS=1.0D-7 )
#include "dgroup.h"
#include "dcbbas.h"

      DIMENSION PRPAO(NTBAS(0),NTBAS(0),NZ)
      CHARACTER PROPNAME*16

      CALL QENTER('CHECK_LLOP')

      IF (IPRINT.GE.6) THEN
        WRITE(LUPRI,'(/A,1X,A16)')
     & '*** CHECK_LLOP: Total AO-matrix of 2comp. (transf.) property ',
     &   PROPNAME 
       CALL PRQMAT(PRPAO(1,1,1),NTBAS(0),NTBAS(0),
     &      NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)

      ENDIF

      DOCCUPB=D0
      DEMPTYB=D0
      IOCCUPB=0
      IEMPTYB=0

      DO I1 = 1, NFSYM
         I2 = MOD(I1+IOPSY,2) + 1
C ... filled matrix blocke LL
         DO IZ = 1, NZ
         DO I = 1, NFBAS(I1,0)
         DO J = 1, NFBAS(I2,0)
          IF (I.LE.NFBAS(I1,1).AND.J.LE.NFBAS(I2,1)) THEN
           IOCCUPB=IOCCUPB+1
           DOCCUPB = DOCCUPB + 
     & DABS(PRPAO(IBAS(I1)+I,IBAS(I2)+J,IZ))
          ELSE
           IEMPTYB=IEMPTYB+1
           DEMPTYB = DEMPTYB + 
     & DABS(PRPAO(IBAS(I1)+I,IBAS(I2)+J,IZ))
          ENDIF   
         ENDDO
         ENDDO
         ENDDO

      ENDDO

      DOCCUPB = DOCCUPB/DFLOAT(IOCCUPB)
      DEMPTYB = DEMPTYB/DFLOAT(IEMPTYB)

      IF (IPRINT.GE.3) THEN
       WRITE(LUPRI,'(/,1X,A,A16,/,A,2D12.7/)')
     &  'CHECK_LLOP:  For property ', PROPNAME,
     &  'we have (occup, empty): ',DOCCUPB,DEMPTYB   
      ENDIF

      IF (DABS(DEMPTYB).GE.DEB_THRS) THEN
       WRITE(LUPRI,'(/,1X,A,A16,/,A,2D12.7,A,D12.7)')
     &  'CHECK_LLOP:  For property ', PROPNAME,
     &  'we have (occup,empty): ',DOCCUPB,DEMPTYB,
     &  ' ... "zero" thershold is ',DEB_THRS   
       WRITE(LUPRI,'(2X,A)')
     & 'WARNING CHECK_LLOP: Non-zero sum of values in a'//
     & ' presumably zero matrix block!!!'
      ENDIF

      CALL QEXIT('CHECK_LLOP')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
C  /* Deck MAKE_FPM */      
       SUBROUTINE MAKE_FPM(TBUF,TMAT,VMAT,EIG,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C      get FP (kinetic) matrix in the ON "theta" basis into TBUF
C
C   On input:
C  -----------
C
C   TMAT - transformation matrix 4c SA-AO "xhi" -> 4c MO "theta"
C           ... always in non-linear (=ordinary) symmetry !
C
C  On output:   
C  -----------   
C
C   TMAT FP (kinetic) matrix in the SA-AO "xhi" basis 
C           (is also written into BSSMAT)
C
C
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C                                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0,
     &           THRNULL = 1.0D-10, THRNULL1 = 1.0D-4)
C  ...  THRNULL is round off error with 15 digits relative to 2c^2 = 3.7D4
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcborb.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbdhf.h"
      DIMENSION TBUF(*),VMAT(*),TMAT(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL DBGDIAG, DOJACO_SAVE, FNDLAB, SBSAVE

      CALL QENTER('MAKE_FPM')
#include "memint.h"


C ... get FP matrix in SA-AO into TBUF
      CALL FREEMT(TBUF,IPRHAM,WORK(KFREE),LFREE) 

C         Transform the free particle matrix
C         TBUF in AO basis
C         in modified Dirac form (orhonormal "theta" basis); 
C
C      VMAT =  TMAT+ * TBUF * TMAT 

      IVMAT = 1
      ITMAT = 1
      DO I = 1,NFSYM
      IF(NORB(I).GT.0) THEN
        IF (IPRHAM.GE.7) THEN
           CALL HEADER(
     &     'MAKE_FPM: Free particle matrix in SA-AO basis "xhi"'//
     &     ' (after call FREEMT)',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &                 NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
           CALL HEADER(
     &     'MAKE_FPM: ENTERING RKB transf. matrix "xhi"->"theta"'//
     &     ' (nonlinear symmetry!), TMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NORB(I),
     &                 NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                             NORB(I),NORB(I),
     &      TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &           VMAT(IVMAT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &           TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &           WORK(KFREE),LFREE,IPRHAM)
C
C       Zero the small elements in real part of WORK(KBUF)
C       which are non-zero because of round-off errors:
C       check that quaternion blocks are numerically zero
C
        DO J = 0,NORB(I)*NORB(I)-1
          IF (DABS(VMAT(IVMAT+J)) .LT. THRNULL) VMAT(IVMAT+J) = D0
        END DO

        IF (NZ .GT. 1) THEN
           IVMAT2 = IVMAT + NORB(I)*NORB(I)
           NVMAT2 = (NZ-1)*NORB(I)*NORB(I)
           TIMAGNORM =  DNRM2(NVMAT2,VMAT(IVMAT2),1)
        ELSE
           TIMAGNORM =  D0
        END IF

C  Write out TIMAGNORM

        IF (TIMAGNORM.GT.THRNULL1) THEN

          WRITE(LUPRI,'(2X,A,A,I1,A,I1,/,2X,A,D12.6)')
     & 'Free particle matrix in the ON basis "theta" ',
     & 'in fermion corep ',I,'/',NFSYM,
     & '*** Norm of the i,j,k imag. parts *** :',TIMAGNORM
          WRITE(LUPRI,'(/2X,A,D13.6)')
     & 'MAKE_FPM: WARNING, TIMAGNORM > THRNULL1 (=',THRNULL1

        ENDIF

        IF (IPRHAM.GE.6) THEN
          CALL HEADER(
     &'MAKE_FPM: Free particle matrix in the ON'//
     & ' RKB basis "theta",VMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &'*** Fermion corep ',I,'/',NFSYM
          WRITE(LUPRI,'(/A,1P,D12.6)')
     & ' Norm of the i,j,k imag. parts of this matrix:',TIMAGNORM
          CALL PRQMAT(VMAT(IVMAT),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IVMAT = IVMAT + NORB(I)*NORB(I)*NZ
        ITMAT = ITMAT + NFBAS(I,0)*NORB(I)*NZ

      ENDIF
      ENDDO


C  ... Td_kin is set to be pure real matrix !
      NZEKIN = 1 
      NZSAVE = NZ
      NZ = NZEKIN 

C     note: DFDIAG zeroes all of TBUF, even if NZ is reset to 1
C     (it uses N2BBASXQ); so we can continue with the quaternion TBUF afterwards.

      DOJACO_SAVE = DOJACO
      DOJACO = .TRUE.
C     ... use JACOBI diagonalization because it does not mix degenerate eigenvectors
C     Copy information from NBORB to the subblock array if we have not
C     partioned the orbital space otherwise.
C
      CALL INISUB

C     ... to tell DFDIAG to use boson symmetries
      SBSAVE = SUB_BL 
      SUB_BL = .TRUE.

      IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/,3X,A,A/)') 
     & 'MAKE_FPM: Diagonalizing (Jacobi) the real part',
     & ' of the full FP matrix in the ON "theta" basis !'
      ENDIF

C ... after diagonalization we get the IBEIG array which is used later !
      IF (IPRHAM.GE.5) THEN 
        DBGDIAG=.TRUE.
      ELSE
        DBGDIAG=.FALSE.
      ENDIF

C  DO copy FPM from VMAT to TBUF - do restore block !
      CALL DZERO(TBUF,N2BBASXQ)
     
C      ... do the diagonalization at once ! Block structure is needed !!!
C      ... needed NTMO, I2TMOT, IORB offsets ...
      CALL DFDIAG (VMAT,EIG,IBEIG,TBUF,DBGDIAG,
     &             WORK,KFREE,LFREE)

      IF (N2TMT.GT.N2BBASXQ) THEN
        CALL QUIT('DIAG-FPM: N2TMT.GT.N2BBASX!')
      ENDIF

      NZ = NZSAVE
      DOJACO = DOJACO_SAVE
      SUB_BL = SBSAVE

      CVAL2 = CVAL*CVAL

      IMAT = 1
      DO  I = 1,NFSYM
      IF(NORB(I).GT.0) THEN

C   ... do always check the RKB free-particle eigenvalues !

        IF (IPRHAM.GE.3) THEN
         CALL HEADER(
     & 'MAKE_FPM: Free particle (in ON basis "theta") eigenvalues :',-1)
          WRITE(LUPRI,'(3X,A,I2,A,I2)') '*** Fermion corep ',I,'/',NFSYM
          WRITE(LUPRI,'(A,3X,A,2X,A,2X,A,3X,A)')
     &    '   Shell no ','FP eigval',
     &    'FP eigval+2c^2','bos.symm.','diff.'
        ENDIF

        TOTDIFF  = D0
        ISHELL = -NPSH(I)
        DO J=1,NORB(I) 
          FPEV   =  EIG(IORB(I)+J)
          FPEV2C = EIG(IORB(I)+J)+CVAL2
C .... this is the boson symmetry of the free-particle eigenshell
          IBS    = IBEIG(IORB(I)+J)+1
C ...  positronic and electronic free particle eigenvalues differ by 2c^2
          DIFF = EIG(IORB(I)+J)+EIG(IORB(I)+NORB(I)-J+1)+CVAL2+CVAL2
          TOTDIFF = TOTDIFF + DABS(DIFF)
          IF (ISHELL.EQ.0) ISHELL = 1
          
C   Print free particle eigenvalues with differences e(+)-e(i)
          IF (IPRHAM.GE.3) THEN
           WRITE(LUPRI,'(I4,1X,I5,2D16.6,3X,I2,3X,D10.4)')
     &             J,ISHELL,FPEV,FPEV2C,IBS,DIFF
          ENDIF
          ISHELL = ISHELL + 1
        ENDDO
        TOTDIFF = TOTDIFF/DFLOAT(NORB(I)/2)
C ... do always print the average difference for RKB free-particle e-p pairing !
        IF (IPRHAM.GE.3) THEN
         WRITE(LUPRI,'(/2X,A,I1,A,I1,A,D12.6)')
     &   'MAKE_FPM: Average difference of the mirror e-p pairing of '//
     &   'free particle eigenvalues (fsym ',I,
     &   '/',NFSYM,') : ',TOTDIFF
        ENDIF

       IF (IPRHAM.GE.7) THEN
        CALL HEADER(
     & 'MAKE_FPM: |k> matrix (FPM eigenvectors)'//
     & ' in the ON basis "theta", TBUF ',-1)
        WRITE(LUPRI,'(/3X,A,I2/)') '*** Fermion corep ',I
        CALL PRQMAT(TBUF(IMAT),NORB(I),NORB(I),
     &                NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C-------------------------------------------------------------------------------
C  Backtransform ELECTRONIC and POSITRONIC solutions from the ON basis "theta"
C  to the AO SO-basis "xhi".
C (we later do the no-pair projection to the electronic solutions after
C  we have done the DKH folding of the nuclear attraction into the
C  electronic part).
C
C  VMAT = TBUF*TMAT 
C-------------------------------------------------------------------------------
        CALL BCKTRA(VMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),
     &              TBUF(IMAT),NORB(I),NORB(I),
     &              NORB(I),NZ,NORB(I),1,NFBAS(I,0),
     &              TMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),NZ,
     &              IPRHAM)

        IF(IPRHAM.GE.7) THEN
          CALL HEADER(
     &    'MAKE_FPM: FP |k> matrix, VMAT (output-TMAT)'//
     &    ' backtransformed to SA-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &    '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(VMAT(I2TMT(I)+1),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

         IMAT = IMAT + NORB(I)*NORB(I)*NZ

      ENDIF    
      ENDDO 

C     ... do export "k" in "xchi" in TMAT
      CALL DCOPY(N2TMT,VMAT,1,TMAT,1)

C  ... do save VMAT - "k" in "xhi" into the file BSSMAT
C ... save Lowdin matrix into BSSMAT
       CALL W2BSSMAT(2,.TRUE.,'FPMAT4C ','MAKE_FPM',
     &                TMAT,N2TMT,LUBSS,IPRHAM)

C ===============================================================
C        ... write help transformation matrix if needed ....
C ===============================================================
       ISIZE = NORB(1)*NORB(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ

       CALL W2BSSMAT(4,.TRUE.,'FPMAT_TH','MAKE_FPM',
     &       TBUF,ISIZE,LUBSS,IPRHAM)

       IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_FPM: 4c "k"-> 4c "theta" transf. matrix FPMAT_TH, TBUF'//
     & ' written to the BSSMAT file.'
       ENDIF

C    ... write also boson irreps...
       ISIZE = NORB(1)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)
       CALL W2BSSMAT(3,.TRUE.,'IBOS_ONB','MAKE_FPM',
     &                IBEIG,ISIZE,LUBSS,IPRHAM)
       IF (IPRHAM.GE.5) THEN
         CALL PRINT_IBEIG(IBEIG,'MAKE_FPM')
       ENDIF
       IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_FPM: IBOS_ONB (boson irreps) of the free-particle basis'//
     & ' were written to into the BSSMAT file.'
       ENDIF

       IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'MAKE_FPM: 4c "k"-xhi transf. matrix, TMAT, is '//
     & ' written to the BSSMAT file, EOFLABEL renewed.'
       ENDIF

      CALL QEXIT('MAKE_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
      SUBROUTINE GET_A_FAC(FACKIN,EIG,VMAT,TMAT,TBUF,
     &                      FPB,USEDF,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:  Construct (and check for PPB) the diagonal elements of the A factor 
C    ========
C                                                         
C    On input: 
C  ------------
C              EIG  - free particle eigenvalues  if FPB=.true.,otherwise zeroes 
C           
C              FPB
C
C   Output: 
C -------------
C            FACKIN - normalization factors   
C            TMAT - the transformation matrix...but it is possible to read it
C                any time from the open BSSMAT file...                               
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C                                                         
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, D2 = 2.0D0, THRNULL = 1.0D-10 ) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h" 
#include "dgroup.h"   
      DIMENSION EIG(*),VMAT(*),TMAT(*),TBUF(*),FACKIN(*),WORK(LWORK)    
      LOGICAL FPB,USEDF,FNDLAB
      CHARACTER*8 LBL 

      CALL QENTER('GET_A_FAC')
#include "memint.h"         


C      ... read TMAT - transformation matrix 4c AO-> 4c MO
      IF (FPB) THEN
        LBL = 'FPMAT4C '
      ELSE IF (.NOT.FPB) THEN
        LBL = 'SL_TM_4C'
      ELSE
       CALL QUIT('GET_A_FAC: wrong LBL !!!')
      ENDIF

      ISIZE = NFBAS(1,0)*NORB(1)*NZ 
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NORB(2)*NZ 
      CALL RFBSSMAT(LBL,.TRUE.,'GET_A_FAC',TMAT,ISIZE,LUBSS,IPRINT)

      CALL GTOVLT(TBUF,D0,IPRINT)

      CVAL2 = CVAL*CVAL

      IMAT = 1
      IFC = 0
      IBUF = 1
      IKUN1 = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

         IF(IPRINT.GE.7) THEN
          CALL HEADER(
     &    'GET_A_FAC: ON RKB/|k> matrix, TMAT (read from file),'//
     &    ' in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')
     &    '* Fermion irp no.',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &    'GET_A_FAC: Obtained (GTOVLT) LL Sovlp matrix'//
     &    ' in SO basis "xhi"',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

** DO C+(e) S(LL) C(e) into TMAT(KKUN) **

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NESH(I),NESH(I),
     &       TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       VMAT(IKUN1),NESH(I),NESH(I),1,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       TMAT(IMAT+NFBAS(I,0)*NPSH(I)),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRINT)

         IF(IPRINT.GE.3) THEN
          CALL HEADER(
     &    'GET_A_FAC: Diagonal (MUST BE!) RKB/<k_e|k_e>=A+.A'//
     &    ' overlap matrix, VMAT',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irep no.',I,'/',NFSYM
          CALL TEST_DIAG_MTX(VMAT(IKUN1),NESH(I),NESH(I),NESH(I),1)
          IF(IPRINT.GE.7) THEN
           CALL PRQMAT(VMAT(IKUN1),NESH(I),NESH(I),
     &          NESH(I),NESH(I),1,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF

*************************************************************
*        Extract the A kinematic factors into FACKIN(1,*)   *
*************************************************************
         IF (IPRINT.GE.3) THEN 
          CALL HEADER(
     &    'GET_A_FAC: "A" kinematic electronic factors',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irep no.',I,'/',NFSYM
          WRITE(LUPRI,'(/14X,A/)')
     &    'A^2_ovlp          A^2_calc          diff.'
         ENDIF

         SUMDIFF = D0
         DO J = 1, NESH(I)
          DAA = VMAT(IKUN1+(J-1)*NESH(I)+J-1)
          FACKIN(IFC+J) = DAA
C         --------------------------------------------
          X1 = (EIG(IORB(I)+NESH(I)+J)/CVAL2) + D2
          X2 = D2*(EIG(IORB(I)+NESH(I)+J)/CVAL2) + D2
          XD = (X1/X2)
C         --------------------------------------------
          DIFF = FACKIN(IFC+J)-XD
          SUMDIFF = SUMDIFF + DABS(DIFF)
          IF (IPRINT.GE.3) THEN
           WRITE(LUPRI,"(1X,I4,1X,2F18.10,4X,D11.5)")
     &     J,FACKIN(IFC+J),XD,DIFF
          ENDIF
         ENDDO 
         SUMDIFF = SUMDIFF/DFLOAT(NESH(I))         
         
         IF (IPRINT.GE.1) THEN
          WRITE(LUPRI,'(/2X,A,I1,A,I1,A,D14.8)')
     & 'GET_A_FAC: Average differences for A-factors (fsym ',
     & I,'/',NFSYM,') : ',SUMDIFF
          CALL FLSHFO(LUPRI)
         ENDIF

         IF (SUMDIFF.GE.THRNULL) THEN
          WRITE(LUPRI,*)
     &'  GET_A_FAC: Sum of differences, (A^2_ovlp-A^2_calc)/N, makes ',
     & SUMDIFF,' while "zero" limit is ',THRNULL
          CALL FLSHFO(LUPRI)
         ENDIF

C       ... off-sets
         IFC   = IFC + NESH(I)
         IKUN1 = IKUN1 + NESH(I)*NESH(I)
         IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
         IBUF = IBUF + NFBAS(I,0)*NFBAS(I,0)

      ENDIF
      ENDDO

      CALL QEXIT('GET_A_FAC')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck mksaomo_le  */
      SUBROUTINE MKSAOMO_Le(TBUF,TMAT,VMAT,FACKIN,FPB,
     &                    IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PURPOSE:
C    ==========
C        Construct half transformed overlap matrix, 
C         Saomo = Sao(LL) * C(L,el) * Ad(ee)^(-1) (FPB=.true.)
C         Saomo = Sao(LL) * C(L,el)   (FPB=.false.)
C
C    On input:
C   -----------
C       TMAT - transformation matrix 
C       FACKIN  Ad(ee) _ if FPB=.true.
C       TBUF -  array
C
C    On output:
C   ------------
C        none; data written to file
C
C      Written:  M.Ilias, Odense, 2003  (based on MKSAOMO by L. Visscher)
C                Last update by M.Ilias, Strasburg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TMAT(*),VMAT(*),TBUF(*),FACKIN(*),WORK(LWORK)
      LOGICAL  FNDLAB,FPB

      CALL QENTER('MKSAOMO_Le')
#include "memint.h"

C     Get AO_LL overlap matrix , S(LL), WITH SS part ZEROED, into VMAT
      CALL GTOVLT(VMAT,D0,IPRINT)

      IOFF = 1
      JOFFESH = 0
      ITMT = 1
      ISAO = 1
      DO I = 1,NFSYM
      IF ( NESH(I) .GT. 0 ) THEN

         IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'MKSAOMO_Le: ENTERING ON RKB/"k" matrix, TMAT,'//
     &    ' in 4c SA-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(ITMT+NFBAS(I,0)*NESH(I)),NFBAS(I,0),NESH(I),
     &               NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &   'MKSAOMO_Le: Read SOVLP matrix(LL) in SO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(ISAO),NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)
         ENDIF

C            SAOMO_Le = S_LL * TMAT
C             TBUF   =  VMAT * TMAT

         CALL QGEMM(NFBAS(I,0),NESH(I),NFBAS(I,0),D1,
     &         'N','N',IPQTOQ(1,0),VMAT(ISAO),
     &         NFBAS(I,0),NFBAS(I,0),1,
     &         'N','N',IPQTOQ(1,0),
     &         TMAT(ITMT+NFBAS(I,0)*NESH(I)),
     &         NFBAS(I,0),NESH(I),NZ,
     &         D0,IPQTOQ(1,0),TBUF(IOFF),
     &         NFBAS(I,0),NESH(I),NZ)
C

            IF ( IPRINT .GE. 6) THEN
             CALL HEADER('MKSAOMO_Le: Right index transf. overlap'//
     &         ' matrix Saomo_Le',-1)
            WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
             CALL PRQMAT(TBUF(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
             CALL FLSHFO(LUPRI)
            END IF

C            Now multiply by diagonal Ad^-1, FACKIN(j)  = Ad^2(j,j)
C               SAOMO_Le = SAOMO_Le * Ad^-1
          IF (FPB) THEN
           ICX = 0
           DO IZ = 1,NZ
           DO IESH = 1, NESH(I)
           DO IFBAS = 1, NFBAS(I,0) 
            TBUF(IOFF+ICX) = 
     &      TBUF(IOFF+ICX)/DSQRT(FACKIN(IESH+JOFFESH))
            ICX = ICX + 1
           ENDDO
           ENDDO 
           ENDDO 
          ENDIF

         IF ( IPRINT .GE. 5 ) THEN
         CALL HEADER('MKSAOMO_Le: FINAL transf. matrix '//
     &   'H2c SA_AO "xhi"  -> H2c ON RKB/FP "k", TBUF ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(TBUF(IOFF),NFBAS(I,0),NESH(I),
     &            NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
            CALL FLSHFO(LUPRI)
          END IF

        JOFFESH = JOFFESH + NESH(I)
        IOFF = IOFF + NFBAS(I,0)*NESH(I)*NZ
        ITMT = ITMT + NFBAS(I,0)*NORB(I)*NZ
        ISAO = ISAO + NFBAS(I,0)*NFBAS(I,0)

      ENDIF
      ENDDO

      ISIZE = NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NESH(2)*NZ
      CALL W2BSSMAT(1,.TRUE.,'TM2CAOFP','MKSAOMO_Le',
     &              TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     & 'MKSAOMO_Le: TM2CAOFP written into BSSMAT, EOFLABEL renewed.'
      ENDIF

      CALL QEXIT('MKSAOMO_Le')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
       SUBROUTINE GET_H1_ONB(TBUF,TMAT,VMAT,IBEIG,
     &          SF_BEG,FPB,LINAD,SAVESO1,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Get the H1 Hamiltonian in the FP ON "k" basis (FPB=.TRUE.)
C    ======== or in the ON RKB MO basis, 
C             if wished (SF_BEG=.true.),eliminate spin-orbit terms 
C             from the beginning using the IBEIG boson info
C
C  On input:
C -----------
C   IBEIG - boson info on MO shells (needed for the spin-free case)
C
C   TMAT  - "k" eigenvectors in "xhi" (when FPB=.true.), or RKB basis
C                                                         
C   FPB  - if .true., the "k" FP ON basis is used, otherwise SL-resorted RKB basis    
C
C   LINAD - flag annoucing additional adaptation due the linear symmetry 
C
C  On output:
C ------------
C     TBUF - H1/FD1 Hamiltonian in the FP ON "k" basis (spin-free,if SF_BEG .true.) 
C                             in RKB ON basis (when FPB=.false.) 
C 
C  Rewritten by Miro ILIAS, Strasburg, Oct. 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"  
C ... needed due to the BSS,USEDF 
#include "dcbham.h"  
      DIMENSION TBUF(*),TMAT(*),VMAT(*),IBEIG(*),WORK(LWORK)    
      LOGICAL SF_BEG, SAVESO1, FPB, LINAD, FNDLAB, EX
      CHARACTER*8 LBL

      CALL QENTER('GET_H1_ONB')
#include "memint.h"         

C=============================================================================================
C
C    In this case we do the 
C   transformation of the (converged) four-component Fock-Dirac operator in the RKB_MO basis. This is
C   the post-DC-SCF approach that gives several new two-component Hamiltonians.
C
C=============================================================================================
      IF (DO4C2C.AND..NOT.START2C) THEN

       IF (IPRHAM.GE.3) THEN
         CALL HEADER( 
     & 'GET_H1_ONB: Block diagonalizing the (converged) Fock-Dirac '//
     & 'operator in MO basis (file DFFOCK).',-1)
       ENDIF

C      ... entering transf. marix in TMAT
       CALL DF2ONB('DFFOCK',TBUF,TMAT,VMAT,LINAD,
     &             FPB,WORK(KFREE),LFREE,IPRHAM)

      ISIZE = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
      CALL  W2BSSMAT(1,.TRUE.,'FD4C_ONB','GET_H1_ONB',
     &                  TBUF,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     &  'GET_H1_ONB: 4c Fock-Dirac(FD4C_ONB) in ON RKB/"k" ready for'//
     &  ' BSS transformations was written (TBUF)'//
     &  ' into the BSSMAT file; EOFLABEL was renewed.'
      ENDIF

C ... now if special case - DFFOCK2 from outside
!       currently this cannot be activated any keyword (unless you know
!       where to "enter" the 6...) - stefan - april 2011
        IF (CONT2C.AND.I2CHAM.EQ.6) THEN
C      ... get transf. matrix into TMAT
          IF (FPB) THEN
            LBL = 'FPMAT4C '
          ELSE 
            LBL = 'SL_TM_4C'
          ENDIF

       ISIZE = NFBAS(1,0)*NORB(1)*NZ 
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NORB(2)*NZ 
       CALL RFBSSMAT(LBL,.TRUE.,'GET_A_FAC',TMAT,ISIZE,LUBSS,IPRHAM)

       CALL DF2ONB('DFFOCK2',TBUF,TMAT,VMAT,LINAD,
     &                FPB,WORK(KFREE),LFREE,IPRHAM)

C     ... save the total 4c FD in ON RKB/"k" into the BSSMAT file !!!
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL W2BSSMAT(1,.TRUE.,'FD4CXONB','GET_H1_ONB',
     &                TBUF,ISIZE,LUBSS,IPRHAM) 
        IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     &  'GET_H1_ONB:Special 4c Fock-Dirac(FD4CXONB)'//
     &  ' in ON RKB/"k" ready for'//
     &  ' BSS transformations was written (TBUF)'//
     &  ' into the BSSMAT file; EOFLABEL was renewed.'
        ENDIF

C        CALL QUIT('GET_H1_ONB: quit for I2CHAM=6!')
        ENDIF

      ENDIF

C================================================================================================
C
C  In this (second) case we do the block-diagonalization of the one-electron Dirac operator.
C
C================================================================================================
      IF (.NOT.(USEDF.AND..NOT.START2C)) THEN

       IF (IPRHAM.GE.3) THEN
         CALL HEADER( 
     & 'GET_H1_ONB: Block diagonalizing the one-electron Dirac '//
     & 'Hamiltonian.',-1)
       ENDIF

C      ... read TMAT - transformation matrix 4c AO-> 4c MO
      IF (FPB) THEN
        LBL = 'FPMAT4C '
      ELSE 
        LBL = 'SL_TM_4C'
      ENDIF

       ISIZE = NFBAS(1,0)*NORB(1)*NZ 
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NORB(2)*NZ 
       CALL RFBSSMAT(LBL,.TRUE.,'GET_H1_ONB',TMAT,ISIZE,LUBSS,IPRHAM)
       IF (IPRHAM.GE.3) THEN
         WRITE(LUPRI,'(/2X,A,A,A)')       
     & 'GET_H1_ONB: ',LBL,' 4c AO-> 4c MO transf. matrix'//
     & ' read from BSSMAT into TMAT.'         
       ENDIF

C First get the bare nucleus Dirac Hamiltonian (BNDH) in the SA-AO basis "xhi"
C into VMAT array

       !CALL MEMCHK('GET_H1_ONB:',WORK,1) ! hurts !

       BSS=.FALSE.
       CALL ONEFCK(VMAT,IPRHAM,WORK(KFREE),LFREE) 
       BSS=.TRUE.

C
C    transform the BNDH from its SO-basis "xhi" to the ON basis "|k>"
C    (reuse TBUF)
C  where H is the bare nucleus Hamiltonian in the SO basis "xhi"
C
C        TBUF = TMAT+ * VMAT  * TMAT; 
C  H1(in |k>) =   U1+ * H0 * U1 (notation in our paper)
C
C  if desired, eliminate spin-orbit components from H1 in C1 symmetry
C  - spin-free at the "beginning"
C

      IF (SF_BEG) THEN
         MZ = 1
      ELSE
         MZ = NZ
      END IF

#ifdef MOD_OOF
! ... perform the QTRANS inside the highly modularized routine
      call handle_matrixes(VMAT,TMAT,WORK(KFREE),LFREE)
#endif

      IOFT = 1
      IMAT = 1

      DO I = 1,NFSYM
      IF(NORB(I).GT.0) THEN

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'GET_H1_ONB: Entering the ON RKB/FP "k" matrix,'//
     &   ' TMAT, in 4c SO-AO basis "xhi" ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GET_H1_ONB: Read Dirac bare nucleus matrix,'//
     &   ' VMAT, in 4c SO-AO basis "xhi" ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       VMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &       TBUF(IOFT),NORB(I),NORB(I),MZ,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'GET_H1_ONB: Obtained Dirac H1 matrix in '//
     &  '4c "k"/ON basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C       ... update the offsets
        IOFT = IOFT + NORB(I)*NORB(I)*NZ
        IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ

      ENDIF
      ENDDO

      ENDIF ! end of branching point for the three different defining h1

C=====================================================================
C     ... do extract (ee) elements from TBUF into VMAT 
C         and eliminate SF terms to get pure SO1 matrix !
C=====================================================================
!     print *, 'bss: SAVESO1 ==> ',SAVESO1
      IF (SAVESO1) THEN
  
       IOFT = 1
       IEE  = 1
       DO I=1,NFSYM
       IF (NORB(I).GT.0) THEN

         N2ESHI = NESH(I)*NESH(I)
         N2T    = NORB(I)*NORB(I)
         DO IZ = 1, NZ
           IBUF = IEE + (IZ-1)*N2ESHI
           JH1EE = IOFT + (NPSH(I)*NORB(I))+NPSH(I)
     &          + (IZ-1)*N2T
           DO J=1,NESH(I)
           CALL DCOPY(NESH(I),TBUF(JH1EE),1,VMAT(IBUF),1)
           JH1EE = JH1EE + NORB(I)
           IBUF  = IBUF  + NESH(I)
           ENDDO
         ENDDO

         IF (IPRHAM.GE.4) THEN
           CALL HEADER(
     & 'GET_H1_ONB: (ee) part of the H1 matrix in ON RB/"k" basis '//
     & ' before eliminating spin-free terms (SAVESO1=.true.),VMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IEE),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL SPF_HK(VMAT(IEE),.FALSE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &            NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRHAM)

         IF (IPRHAM.GE.4) THEN
           CALL HEADER(
     & 'GET_H1_ONB: (ee) part of the H1 matrix in ON RB/"k" basis '//
     & ' after eliminating spin-free terms-pure SO1 Ham., VMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IEE),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
         ENDIF
         IEE  = IEE + NESH(I)*NESH(I)*NZ
         IOFT = IOFT + NORB(I)*NORB(I)*NZ

       ENDIF
       ENDDO

C  ...  save H_SO1(ee) terms into the BSSMAT file ...
      ISIZE = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NESH(2)*NZ)
      CALL W2BSSMAT(2,.TRUE.,'HSO1_ONB','GET_H1_ONB',
     &               VMAT,ISIZE,LUBSS,IPRHAM)      

      ISIZE=NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE=ISIZE+NORB(2)*NORB(2)*NZ
      CALL W2BSSMAT(3,.TRUE.,'H1_ONB  ','GET_H1_ONB',
     &               TBUF,ISIZE,LUBSS,IPRHAM)      
      IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'GET_H1_ONB: H1_SO1 & H1_ONB in ON RKB/"k" basis have been '//
     & 'written to the BSSMAT file, EOFLABEL is renewed.'
       CALL FLSHFO(LUPRI)
      ENDIF

      ENDIF

C .... if spin-free, eliminate spin-orbit dependent terms
      IF (SPINFR.AND.SF_BEG) THEN
        IOFT = 1
        DO I = 1, NFSYM
        IF (NORB(I).GT.0) THEN
         CALL SPF_HK( TBUF(IOFT),.TRUE.,IBEIG(IORB(I)+1),
     &                IBEIG(IORB(I)+1),
     &             NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRHAM )
         IF(IPRHAM.GE.5) THEN
          CALL HEADER(
     &  'GET_H1_ONB: H1 Dirac matrix in 4c ON RKB/"k" basis '//
     &  ' after elim. spin-orbit dependent terms "BEFORE", TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL FLSHFO(LUPRI)
         ENDIF
         IOFT = IOFT + NORB(I)*NORB(I)*NZ
        ENDIF
        ENDDO
      ENDIF

CMI ... control printout of IBEIG
      IF (IPRHAM.GE.5) THEN
        CALL PRINT_IBEIG(IBEIG,'GET_H1_ONB')
      ENDIF

C ... write the TBUF - H1 ON RKB/"k" to the file...
      ISIZE=NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE=ISIZE+NORB(2)*NORB(2)*NZ
      CALL W2BSSMAT(1,.TRUE.,'H1_ONBAS','GET_H1_ONB',
     &               TBUF,ISIZE,LUBSS,IPRHAM)      
      IF (IPRHAM.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'GET_H1_ONB:  H1_ONBAS (in ON RKB/"k") matrix is '//
     &' written to the BSSMAT file, EOFLABEL renewed.'
       CALL FLSHFO(LUPRI)
      ENDIF

      CALL QEXIT('GET_H1_ONB')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
      SUBROUTINE PRINT_IBEIG(IBEIG,CALLED_FROM)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Prints out IBEIG - boson info on MO's in this separate routine. It's important.
C 
C  Written by Miro Ilias, Sept.2010, Banska Bystrica
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
#include "dgroup.h"  
#include "dcborb.h"   
      DIMENSION IBEIG(*)
      CHARACTER*(*) CALLED_FROM

      WRITE(LUPRI,'(/,2X,A,A,A)')
     & '=====  PRINT_IBEIG called from ',CALLED_FROM,'  ====='
      DO I=1,NFSYM
         WRITE(LUPRI,'(2X,A,I2,A,I1,1X,A,I4)')
     &   '* boson symmetries for fermion corep ',
     &   I,' /',NFSYM,':  # of orbitals:',NORB(I)
         WRITE(LUPRI,*) ' ...negative ("positronic") states:'
         DO J=1,NPSH(I)
           WRITE(LUPRI,'(8X,I3,A,I2)')
     &     -(NPSH(I)-J+1),"--->",IBEIG(IORB(I)+J)
         ENDDO
         WRITE(LUPRI,*) ' ...positive (electronic) states:'
         DO J=1,NESH(I)
           WRITE(LUPRI,'(8X,I3,A,I2)')
     &     J,"--->",IBEIG(IORB(I)+NPSH(I)+J)
         ENDDO
      ENDDO
      CALL FLSHFO(LUPRI)
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
      SUBROUTINE DF2ONB(DFFNAME,TBUF,TMAT,VMAT,
     &                 LINAD,FPB,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Purpose: Read the converged Fock MO matrix from the DFFNAME file,
C           transform it to the 4c "k"/RKB ON basis to be ready for the block-diagonalization
C
C  On input: DFFNAME - name of the file
C            LINAD -  flag annoucing additional adaptation due the linear symmetry 
C            TMAT  - "k" eigenvectors in "xhi" (when FPB=.true.), or RKB basis (in one-step transformation)
C
C  On output: TBUF - resulting 4c matrix in "k"/RKB ON basis
C
C
C  Called from:  GET_H1_ONB
C
C  Written by Miro Ilias, dec.2005, Strasbourg
C  Last changes:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( D1 = 1.0D0, DM1= -1.0D0, D0 = 0.0D0,THRNULL = 1.0D-10) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
#include "dcbdhf.h"
      DIMENSION TBUF(*),TMAT(*),VMAT(*),WORK(LWORK)
      CHARACTER*(*) DFFNAME
      LOGICAL LINAD,FPB,EX,FNDLAB

      CALL QENTER('DF2ONB')
#include "memint.h"         

        IPRHAM = IPRINT ! no dcbham.h

C  ... check the presence of DFFOCK !
       INQUIRE(FILE=DFFNAME,EXIST=EX)
       IF (.NOT.EX) THEN
        write(lupri,*) 'DFFNAME=',DFFNAME
        WRITE(LUPRI,*)
     &  'DF2ONB: The required ',DFFNAME ,' file for the DC2IOTC'//
     &  ' transformation is not present !'
        CALL QUIT('DF2ONB: The file is not present !')
       ENDIF

C ...  read converged Fock-Dirac MO (may be LINSYM!) matrix from the DFFOCK file
       OPEN (LUFOCK, FILE=DFFNAME,FORM='UNFORMATTED',
     &       ACCESS='DIRECT',RECL=8*N2TMOTQ,STATUS='OLD')
C ... RESORT OF RKB_TMAT IS NEEDED FOR THE RKB based transformation !!!
       IREC = -1
       IF(MXDIIS.EQ.0) THEN
          IREC = 1
       ELSE IF (NITER.GT.0) THEN
          IREC = MOD(NITER-2,MXDIIS-1) + 1
       ELSE
          IREC = 1
       ENDIF
       IF (IREC.EQ.0) IREC = 1
       CALL READAC(LUFOCK,N2TMOTQ,VMAT,IREC)
C      ... we no longer need the DFFOCK file 
C       as the next BSS-SCF creates own new files
       CLOSE(LUFOCK,STATUS="DELETE")
       IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A,A,A)')
     & 'DF2ONB: Fock-Dirac 4c MO (RKB) matrix '//
     & 'was read from the ',DFFNAME,
     & ' file(deleted afterwards) file into VMAT'
       ENDIF
 
       IF (LINAD) THEN
C       ... do the transformation back to the nonlinsym RKB_MO basis...
C   =======================================================================
       ISIZE = NORB(1)*NORB(1)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)
       CALL RFBSSMAT('TMAT_LSY',.TRUE.,'DF2ONB',TMAT,ISIZE,LUBSS,IPRHAM)
       IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     & 'DF2ONB: 4c ON "theta"-> 4c ON LINSYM "theta" transf. '//
     & 'matrix TMAT_LSY was read from the BSSMAT file into TMAT.'
       ENDIF

C     ... do the reverse transformation (MOAO):  TMAT . VMAT . TMAT + = TBUF
       IOFT = 1
       DO I = 1, NFSYM
       IF (NORB(I).GT.0) THEN

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: Read (from DFFOCK) Fock-Dirac matrix in '//
     &  '4c ON MO basis "theta", VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IOFT),NORB(I),NORB(I),NORB(I),
     &              NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
   
C         ... do transform it to the "k" basis into TBUF !
        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: Transformation matrix '//
     &  ' 4c ON RKB "theta" -> LSYM "theta", TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('MOAO','S',D0,NORB(I),NORB(I),
     &       NORB(I),NORB(I),
     &       TBUF(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       VMAT(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: Dirac-Fock matrix in '//
     &  '4c ON RKB (non_linear symmetry!) basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &              NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        
         IOFT = IOFT + NORB(I)*NORB(I)*NZ
        ENDIF
        ENDDO

        ISIZE = NORB(1)*NORB(1)
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)
C         VMAT contains FD in ordinary RKB basis ... 
        CALL DCOPY(ISIZE,TBUF,1,VMAT,1)
        IF (IPRHAM.GE.2) THEN
           WRITE(LUPRI,'(2X,A)')
     &  'DF2ONB: The TBUF was copyied into VMAT'
        ENDIF
      ENDIF
C   ! end of linear, we have Fock-Dirac nonlinear_MO matrix in VMAT...

C   ... if we are in the "k" basis mode, we need transformation matrix
C       RKB "theta" -> FPM "k"
       IF (FPB) THEN

       ISIZE = NORB(1)*NORB(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
       CALL RFBSSMAT('FPMAT_TH',.TRUE.,'DF2ONB',TMAT,ISIZE,LUBSS,IPRHAM)
        IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     &    'DF2ONB: 4c "k"-> 4c "theta" transf. '//
     &    'matrix FPMAT_TH was read from the BSSMAT file into TMAT.'
        ENDIF

       ELSE 

C      ... when being in 1.step BSS, we need to transform to the AO basis first
C       (using "standard" RKB AO-MO transformation matrix) and 
C       then transform it to the SL-sorted MO basis...

       ISIZE = (NORB(1)*NFBAS(1,0)*NZ)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NFBAS(2,0)*NZ)
       CALL RFBSSMAT('TMAT4C  ',.TRUE.,'DF2ONB',TMAT,
     &                 ISIZE,LUBSS,IPRHAM)
        IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     &    'DF2ONB: 4c AO -> 4c "theta" transf. '//
     &    'matrix TMAT4C was read from the BSSMAT file into TMAT.'
        ENDIF
 
         CALL DZERO(TBUF,N2BBASXQ)
         CALL MKSAOMO(TBUF,TMAT,.FALSE.,IPRHAM,WORK(KFREE),LFREE)
         CALL DCOPY(N2BBASXQ,TBUF,1,TMAT,1)

C       ... do get AO Fock-Dirac matrix (MOAO): TMAT . VMAT . TMAT^+ = TBUF
          IBUF = 1
          IVMAT = 1
          ITMAT = 1
          DO I = 1,NFSYM
          IF (NORB(I).GT.0) THEN

         CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                           NORB(I),NORB(I),
     &      TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),
     &      VMAT(IVMAT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN

         CALL HEADER(
     &  'DF2ONB: Dirac-Fock matrix '//
     &  'in MO basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &  'DF2ONB:4c AO -> 4c SAOMO "theta" transf. matrix '//
     &  ', TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)


         CALL HEADER(
     &  'DF2ONB: Dirac-Fock matrix '//
     &  'retransformed into 4c SA-AO basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),
     &        NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF

            IBUF = IBUF + (NFBAS(I,0)*NFBAS(I,0)*NZ)
            IVMAT = IVMAT + (NORB(I)*NORB(I)*NZ)
            ITMAT = ITMAT + (NFBAS(I,0)*NORB(I)*NZ )
          ENDIF
          ENDDO

C       ... read SL_TM_4C into TMAT and do TMAT^+ . TBUF . TMAT = VMAT
        CALL DZERO(TMAT,N2BBASXQ)

        ISIZE = (NORB(1)*NFBAS(1,0)*NZ)
        IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NFBAS(2,0)*NZ)
        CALL RFBSSMAT('SL_TM_4C',.TRUE.,'DF2ONB',
     &                 TMAT,ISIZE,LUBSS,IPRHAM)
        IF (IPRHAM.GE.2) THEN
          WRITE(LUPRI,'(/2X,A)')
     &    'DF2ONB: 4c AO -> 4c "theta" SL_TM_4C transf. '//
     &    '  matrix was read from the BSSMAT file into TMAT.'
        ENDIF

          CALL DZERO(VMAT,N2BBASXQ)

          IBUF = 1
          IVMAT = 1
          ITMAT = 1
          DO I = 1,NFSYM
          IF (NORB(I).GT.0) THEN

         CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                           NORB(I),NORB(I),
     &      TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),
     &      VMAT(IVMAT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN

         CALL HEADER(
     &  'DF2ONB:4c AO -> 4c "SL theta" transf. matrix '//
     &  ', TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &  'DF2ONB: Dirac-Fock matrix '//
     &  'retransformed into 4c SA-AO basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,0),
     &        NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)


         CALL HEADER(
     &  'DF2ONB: Dirac-Fock matrix in '//
     &  '4c "SL-sorted" RKB_MO basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

          IBUF = IBUF + (NFBAS(I,0)*NFBAS(I,0)*NZ)
         IVMAT = IVMAT + (NORB(I)*NORB(I)*NZ)
         ITMAT = ITMAT + (NFBAS(I,0)*NORB(I)*NZ )
          ENDIF
          ENDDO

          CALL DCOPY(N2BBASXQ,VMAT,1,TBUF,1)

       ENDIF

       IOFT = 1
       DO I = 1, NFSYM
       IF (NORB(I).GT.0) THEN

         IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: Read (from DFFOCKx)/retransformed  '//
     &  'Fock-Dirac matrix in 4c ON MO basis "theta", VMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
   
C         ... do transform it to the "k" basis into TBUF !
        IF (FPB) THEN

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: Transformation matrix '//
     &  ' 4c ON  RKB "theta" -> FP "k" basis, TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('AOMO','S',D0,NORB(I),NORB(I),
     &       NORB(I),NORB(I),
     &       VMAT(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       TBUF(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       TMAT(IOFT),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: Dirac-Fock matrix in '//
     &  '4c ON FP "k" basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        ENDIF

        IF (.NOT.FPB.AND.IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'DF2ONB: 4c Dirac/Fock-Dirac matrix in '//
     &  '4c ON RKB/"k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        
        IOFT = IOFT + NORB(I)*NORB(I)*NZ

       ENDIF
       ENDDO

      CALL QEXIT('DF2ONB')
      RETURN
      END
#ifdef MOD_UNRELEASED
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hinf_bd */      
      SUBROUTINE GET_HINF_BD(TBUF,TMAT,VMAT,IBEIG,
     &                   GETPCTM,SF_BEG,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Get the H_inf 2c Hamiltonian in the ON RKB/"k" basis
C    ========  
C             - save H1 ON RKB/"k"_block diagonal (full, both pp and ee blocks)
C             - save 4c ON RKB/"k"->2c ON RKB/"k" (full) picture change transf.matrix
C             - save 4c "xhi" -> 2c "k" (only ee) picture change transf.matrix 
C             - save Hinf(2c) in "k" basis to the file
C
C  On input:
C --------------
C    TBUF - BNDH in ON RKB/"k" basis
C    VMAT ON RKB/"k" in SA-AO "xhi" basis
C
C  On output:
C --------------
C    TBUF - H2c (ee) in ON RKB/"k" basis set, needed in the subsequent H2CFINAL routine
C
C Written by:  M.Ilias, Strasbourg, February 2006   
C                 revisons: M.Ilias, Prievidza, November 2006
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( D1 = 1.0D0, DM1= -1.0D0, D0 = 0.0D0,THRNULL = 1.0D-10) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
      DIMENSION TBUF(*),TMAT(*),VMAT(*),IBEIG(*),
     &          WORK(LWORK)    
      DIMENSION IBL(1,4)
      LOGICAL   SPINFR_SAVE,DBGDIAG,GETPCTM,SF_BEG,FNDLAB

      CALL QENTER('GET_HINF_BD')
#include "memint.h"         

      NBL = 1

      IOFT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &  'GET_HINF_BD: Entering H1 Dirac matrix in 4c "k" basis '//
     &  'TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
       
C     .. starting point
        IBL(1,1)=1
        IBL(1,2)=NPSH(I)+1
C     ... ending point
        IBL(1,3)=NPSH(I)
        IBL(1,4)=NORB(I)

C     one disadvantage  ... this routine is taking too long time ...
        CALL QJACOBI(TBUF(IOFT),TMAT(IOFT),NORB(I),
     &               NZ,NBL,IBL,.TRUE.,IPRINT)

        IF(IPRINT.GE.5) THEN
         CALL HEADER(
     &  'GET_HINF_BD: Block-diagonalized full matrix H2c '//
     &  'in 4c ON RKB/"k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NORB(I),NORB(I),NORB(I),
     &               NORB(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'GET_HINF_BD: Full picture change tranf.matrix,  '//
     &  ' 4c ON RKB/"k" -> 2c ON RKB/"k", TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IOFT),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)

        ENDIF
        IOFT = IOFT + NORB(I)*NORB(I)*NZ
      ENDIF
      ENDDO

      ISIZE1=NFBAS(1,0)*NORB(1)*NZ
      ISIZE2=NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) THEN
        ISIZE1=ISIZE1+NFBAS(2,0)*NORB(2)*NZ
C       ... full matrix ....
        ISIZE2=ISIZE2+NORB(2)*NORB(2)*NZ
      ENDIF

C ... save he H1 block diagonalized Hamiltonian 
      CALL W2BSSMAT(1,.TRUE.,'H1ONB_BD','GET_HINF_FPM',
     &              TBUF,ISIZE2,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_BD: H1 ON RKB/"k" block-diagonalized matrix is '//
     & ' written to the BSSMAT file.'
      ENDIF

C     ... save the full picture change transformation matrix
      CALL W2BSSMAT(1,.TRUE.,'PCTM_ONB','GET_HINF_BD',
     &                 TMAT,ISIZE2,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'GET_HINF_BD: Full 4c ON RKB/"k" -> 2c ON/"k" picture change '//
     &'transformation matrix written to the BSSMAT file (PCTM_ONB).'
      ENDIF

C  ... get the picture change transformation matric 4c "xhi" -> 2c "k" K.U 
C  ... into the TBUF = VMAT*TMAT
      IMAT=1
      IPC =1 
      IOFT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GET_HINF_BD: Entering RKB/"k" matrix, VMAT,'//
     & ' in 4c SA-AO basis "xhi" ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IMAT),NFBAS(I,0),NORB(I),
     &                NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QGEMM(NFBAS(I,0),NESH(I),NORB(I),D1,
     &         'N','N',IPQTOQ(1,0),VMAT(IMAT),
     &         NFBAS(I,0),NORB(I),NZ,
     &         'N','N',IPQTOQ(1,0),
     &         TMAT(IPC+NORB(I)*NPSH(I)),
     &         NORB(I),NESH(I),NZ,
     &         D0,IPQTOQ(1,0),TBUF(IOFT),
     &         NFBAS(I,0),NESH(I),NZ)

        IF(IPRINT.GE.5) THEN
         CALL HEADER(
     & 'GET_HINF_BD: 4c "xhi" - 2c ON RKB/"k" picture change'//
     & ' transformation matrix ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IOFT),NFBAS(I,0),NESH(I),
     &          NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IMAT = IMAT + NFBAS(I,0)*NORB(I)*NZ
        IPC  = IPC  + NORB(I)*NORB(I)*NZ
        IOFT = IOFT + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO
 
C  ... save the picture change transformation matrix, TBUF...
      ISIZE1=NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE1=ISIZE1+NFBAS(2,0)*NESH(2)*NZ

      CALL W2BSSMAT(1,.TRUE.,'PCTM4C2K','GET_HINF_BD',
     &              TBUF,ISIZE1,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GET_HINF_BD: 4c SA-AO"xhi" -> 2c "k" picture change '//
     & 'transformation matrix (PCTM4C2K) is written to the BSSMAT file.'
        WRITE(LUPRI,'(2X,A)') 'EOFLABEL is added.'
      ENDIF

C     ... extract U1_ONBAS, TMAT -> TBUF
      IMAT = 1
      IBUF = 1
      DO I=1,NFSYM
       DO IZ=1,NZ
       DO J=1,NESH(I)
        CALL DCOPY(NORB(I),
     &    TMAT(IMAT+(NORB(I)*NPSH(I))+
     &        (IZ-1)*NORB(I)*NORB(I)+(J-1)*NORB(I)),1,
     &    TBUF(IBUF+(IZ-1)*NORB(I)*NESH(I)+(J-1)*NORB(I)),1)
       ENDDO
       ENDDO
       IF (IPRINT.GE.6) THEN
        CALL HEADER(
     &  'GET_HINF_BD: Extracted U1(+e) picture change '//
     &  'transformation matrix 4c->2c ON RKB/"k" basis, TBUF',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(IBUF),NORB(I),NESH(I),NORB(I),NESH(I),
     &                NZ,IPQTOQ(1,0),LUPRI)
       ENDIF
       IMAT = IMAT + NORB(I)*NORB(I)*NZ
       IBUF = IBUF + NORB(I)*NESH(I)*NZ
      ENDDO

      ISIZE = NESH(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NORB(2)*NZ
      CALL W2BSSMAT(1,.TRUE.,'U1_ONBAS','GET_HINF_BD',
     &              TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'GET_HINF_BD: 4c ON RKB/"k"-> 2c ON RKB/"k" (ee) picture change'//
     &' transf. matrix (U1_ONBAS) saved, EOFLABEL renewed '
      ENDIF

C  ... read the block diagonal matrix for extracting the (ee) block only !
       CALL RFBSSMAT('H1ONB_BD',.TRUE.,'GET_HINF_BD',TBUF,
     &               ISIZE2,LUBSS,IPRINT)
      
      IBUF = 1
      IMAT = 1
      DO I = 1, NFSYM

       IF(IPRINT.GE.5) THEN
         CALL HEADER(
     & 'GET_HINF_BD: Full H1 ON RKB/"k" BLOCK DIAGONALIZED '//
     & ' Hamiltonian, TBUF ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C  ... extract ee (bottom right block) ,  TBUF -> VMAT
       DO IZ=1,NZ
       DO J=1,NESH(I)
         CALL DCOPY(NESH(I),
     &   TBUF(IBUF+(NORB(I)*NPSH(I))+NESH(I)+
     &              (IZ-1)*NORB(I)*NORB(I)+(J-1)*NORB(I)),1,
     &   VMAT(IMAT+(IZ-1)*NESH(I)*NESH(I)+(J-1)*NESH(I)),1)
       ENDDO
       ENDDO

       IF(IPRINT.GE.5) THEN
        CALL HEADER(
     &  'GET_HINF_BD: Extracted H2c RKB/"k" (ee)'//
     &  ' Hamiltonian, VMAT ',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IMAT),NESH(I),NESH(I),
     &       NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUF = IBUF + NORB(I)*NORB(I)*NZ
       IMAT = IMAT + NESH(I)*NESH(I)*NZ
      ENDDO

      ISIZE3 = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE3 = ISIZE3 + NESH(2)*NESH(2)*NZ
      CALL W2BSSMAT(1,.TRUE.,'H2C_ONBA','GET_HINF_BD',
     &              VMAT,ISIZE3,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'GET_HINF_BD: H2c_infinite in ON RKB/"k" (ee) written to BSSMAT'//
     &' (H2C_ONBA), EOFLABEL renewed '
      ENDIF

C Copy VMAT into TBUF as TBUF is leaving this routine...
      CALL DCOPY(ISIZE3,VMAT,1,TBUF,1)

      CALL QEXIT('GET_HINF_BD')
      RETURN
      END
#endif
      SUBROUTINE DIAG_h_test(hmat,eigvl,eigvc,ndim,nz,ierr)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       write(6,*) ' DIAG_h_test: call to printaddr'
       call printaddr(hmat)
       call qdiag90(nz,ndim,hmat,ndim,ndim,eigvl,1,eigvc,ndim,ndim,ierr)
       END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hinf_r */      
      SUBROUTINE GET_HINF_R(TBUF,TMAT,VMAT,FACKIN,EIG,IBEIG,FPB,
     &                   YREQ1,GETPCTM,SF_BEG,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE: Get the H_inf 2c Hamiltonian in the  ON RKB/"k" basis
C  ========  
C  Use the method of "R" solving to get picture change (in ON basis)
C  transformation matrix U and after applying U^+ H1 U get the H_inf 2c
C  Hamiltonian in a given ON basis.
C
C    On input:
C  ------------
C    TBUF - BNDH in ON RKB / FP "k" basis
C    IBEIG - only helping array for QCHOLD/QCHOLS, its content is hurt !
C
C   On output:
C  ------------
C    TBUF - 2 comp. infinite order Hamiltonian in ON RKB / FP "k" basis (ee block)
C
C    Written by:  M.Ilias, HJ Aa Jensen, Odense, 2002,
C  last revison: M.Ilias, T.Saue, Strasbourg, february 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( D1 = 1.0D0, DM1= -1.0D0, D0 = 0.0D0,THRNULL = 1.0D-10) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   
      DIMENSION TBUF(*),TMAT(*),VMAT(*),FACKIN(*),EIG(*),IBEIG(*),
     &          WORK(LWORK)    
      LOGICAL  SPINFR_SAVE,DBGDIAG,GETPCTM,SF_BEG, YREQ1,
     &          FNDLAB,FPB,FP2TH,FP2TH1,FP2TH2,USEJACO

      CALL QENTER('GET_HINF_R')
#include "memint.h"         

      IBUF = 1
      IEIG = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN     

        IF(IPRINT.GE.5) THEN
         CALL HEADER(
     &  'GET_HINF_R (0): ENTERING 4c H1/FD1 bare nucl. matrix in '//
     &  ' ON RKB/FP "k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C  ... diagonalize the H1 in ON RKB/"k" basis
       CALL QDIAG(NZ,NORB(I),TBUF(IBUF),NORB(I),NORB(I),
     &            EIG(IEIG),1,TMAT(IBUF),NORB(I),NORB(I),
     &            WORK(KFREE),LFREE,IERR)

!       call print_x2cmat(TMAT(IBUF),NORB(i),NORB(i),nz,IPQTOQ(1,0),
!    &                    'bss - h1_4c-evc',LUPRI)
        IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &   'GET_HINF_R (0): Eigenvalues of H1 bare nucl. matrix in '//
     &   ' ON "k"/RKB basis, EIG',-1)
         IX=0
         DO J=1,NORB(I)
          WRITE(LUPRI,'(2X,I4,2X,F16.6)') J-1-NPSH(I)+IX,EIG(IEIG+J-1)
          IF (J.EQ.NESH(I)) IX=1
         ENDDO

         CALL HEADER(
     &   'GET_HINF_R (0): Eigenvectors of H1 bare nucl. matrix in '//
     &   'ON "k"/RKB basis, TMAT',-1)
         CALL PRQMAT(TMAT(IBUF),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUF = IBUF + NORB(I)*NORB(I)*NZ
       IEIG = IEIG + NORB(I)

      ENDIF
      ENDDO

C  .... next (control) cycle ... do transform HDBN "k" eigenvectors 
C       to the HDBN "theta" eigenvectors
      FP2TH=.FALSE.
      IF (FP2TH.AND.FPB) THEN
C ... read the FPMAT_TH 4c "k"-> 4c "theta" transf. matrix into TBUF
        REWIND(LUBSS)
        IF (FNDLAB('FPMAT_TH',LUBSS)) THEN
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL READT(LUBSS,ISIZE,TBUF)
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R (TEST): 4c "k"-> 4c "theta" transf. '//
     & '  matrix FPMAT_TH was read from the BSSMAT file into TBUF.'
        ELSE
         CALL QUIT(
     &  'GET_HINF_R: Error, EOFLABEL not found in BSSMAT file!')
        END IF

        IBUF = 1
        DO I=1,NFSYM
        IF (NORB(I).GT.0) THEN
          IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &   'GET_HINF_R (T1): Eigenvectors of H1 bare nucl. matrix in '//
     &   'ON "k"/RKB basis, TMAT',-1)
         CALL PRQMAT(TMAT(IBUF),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &  'GET_HINF_R (T1): 4c "k"-> 4c "theta" transf.'//
     & ' matrix FPMAT_TH  in TBUF',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C       ... VMAT: =  TBUF*TMAT
        CALL QGEMM(NORB(I),NORB(I),NORB(I),D1,
     &     'N','N',IPQTOQ(1,0),
     &     TBUF(IBUF),NORB(I),NORB(I),NZ,
     &     'N','N',IPQTOQ(1,0),
     &      TMAT(IBUF),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),VMAT(IBUF),NORB(I),NORB(I),NZ)

          IF(IPRINT.GE.5) THEN
            CALL HEADER(
     &  'GET_HINF_R (T1): Eigenvectors of H1("k") transformed'//
     &  ' to the "theta" basis, VMAT',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
            IBUF = IBUF + NORB(I)*NORB(I)*NZ
        ENDIF
        ENDDO
      ENDIF

      IBUF = 1
      IBUFE = 1
      IEIG = 1
      IBMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN     
       NREF   = NESH(I)
       N2REF  = NREF*NREF
       N2EQ   = NREF*NREF*NZ
       N2REFQ = N2REF*NZ

       IF (YREQ1) THEN
C ... make A and B matrixes as in the original (HJAaJ) approach...
       CALL DUNIT2(TBUF(IBUFE),NREF,NREF,NREF,NZ)
C ... do A := A - Yp*(Yp)+ into TBUF(IBUFE), A = 1 - Y-+.Y-+^+
       CALL QGEMM(NREF,NREF,NREF,DM1, 
     &    'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF+(NPSH(I)*NORB(I))),NORB(I),NORB(I),NZ,
     &    'H','N',IPQTOQ(1,0),
     &     TMAT(IBUF+(NPSH(I)*NORB(I))),NORB(I),NORB(I),NZ,
     &     D1,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C  -- make b matrix (store in VMAT); b := Yp * (Ye)+ = Y-+.Y++^+
       CALL QGEMM(NREF,NREF,NREF,D1,
     &     'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF+(NPSH(I)*NORB(I))),NORB(I),NORB(I),NZ,
     &     'H','N',IPQTOQ(1,0),
     &      TMAT(IBUF+(NPSH(I)*NORB(I))+NPSH(I)),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ)

      ELSE
C       .....  Trond's suggestion (YR1=.false.)  .....
C        A = Y--.Y--^+
       CALL QGEMM(NREF,NREF,NREF,D1, 
     &    'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NORB(I),NORB(I),NZ,
     &    'H','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NORB(I),NORB(I),NZ,
     &     D0,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C  -- make b matrix (store in VMAT); B := - Y--.Y+-^+
       CALL QGEMM(NREF,NREF,NREF,DM1,
     &     'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NORB(I),NORB(I),NZ,
     &     'H','N',IPQTOQ(1,0),
     &      TMAT(IBUF+NPSH(I)),NORB(I),NORB(I),NZ,
     &      D0,IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ)
       ENDIF

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GET_HINF_R (1): A matrix for getting R in ON basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GET_HINF_R (1): B matrix for getting R in ON basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUF = IBUF + NORB(I)*NORB(I)*NZ
       IBUFE = IBUFE + N2EQ
       IEIG = IEIG + NORB(I)
      ENDIF
      ENDDO

      IF (IPRINT.GE.5) THEN
        CALL PRINT_IBEIG(IBEIG,'GET_HINF_R (1)')
      ENDIF

C ... now we have eigenvectors: prepare for solving R equation !
C        (1-Y-+.Y-+^+) R = Y-+.Y++^+ or
C        (Y--.Y--^+) R = -Y--.Y+-^+

      IEIG = 1
      IBUFE =1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF   = NESH(I)
       N2REF = NREF*NREF
       N2EQ = NREF*NREF*NZ
       N2REFQ = N2REF*NZ

       CALL MEMGET('REAL',KW,NREF,WORK,KFREE,LFREE)  

C ... THRUNULL of 1.0D-20 OR 0.0D0 HAS NO INFLUENCE...
       JOB = 1
       CALL QCHOLD(TBUF(IBUFE),NREF,NZ,NREF,NREF,EIG(IEIG),
     &                  THRNULL,NEFF,JOB,IBEIG(IEIG))

       IF(NEFF.LT.NREF) THEN
         WRITE(LUPRI,'(/2A)') '* GET_HINF_R WARNING: ',
     &   'linear dependencies detected in Cholesky decomposition !'
         WRITE(LUPRI,'(11X,A,I5,A,I5)')
     &   'Reduced dimensionality: ',NREF,' --> ',NEFF
         CALL QUIT('GET_HINF_R: Cholesky decompos.problem !')
       ENDIF
       CALL QCHOLS(TBUF(IBUFE),NREF,NEFF,NREF,NZ,NREF,NREF,
     &           EIG(IEIG),VMAT(IBUFE),NREF,NREF,
     &                     TMAT(IBUFE),NREF,NREF,
     &           JOB,IBEIG(IEIG),WORK(KW))

C   ... print the R
       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'GET_HINF_R (1): Obtained R after CHoleski, TMAT',-1)
         IF (.NOT.YREQ1) THEN
           write(lupri,*) "Trond's form of the YR equation."
         ELSE
           write(lupri,*) "HJAaJ's form of the YR equation."
         ENDIF
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL MEMREL('W    ',WORK,KWORK,KW,KFREE,LFREE)

       IEIG = IEIG + NREF
       IBUFE = IBUFE + N2EQ

      ENDIF
      ENDDO

      IF (IPRINT.GE.5) THEN
        CALL PRINT_IBEIG(IBEIG,'GET_HINF_R (2)')
      ENDIF

C ... save the full matrix R, stored in TMAT, into the BSSMAT file
      ISIZE = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
      CALL W2BSSMAT(1,.TRUE.,'R_MTX   ','GET_HINF_R',
     &              TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &'R matrix (ON basis), TMAT, written to BSSMAT, EOFLABEL renewed'
       CALL FLSHFO(LUPRI)
      ENDIF

C =========== TESTING "k" AND RKB approaches ======================
C      ... does not work for R !!!
      FP2TH1=.FALSE.
      IF (FP2TH1.AND.FPB) THEN
C ... read the FPMAT_TH 4c "k"-> 4c "theta" transf. matrix into TBUF
        REWIND(LUBSS)
        IF (FNDLAB('FPMAT_TH',LUBSS)) THEN
        ISIZE = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
        CALL READT(LUBSS,ISIZE,TBUF)
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R (T2): 4c "k"-> 4c "theta" transf. matrix FPMAT_TH '//
     & ' was read from the BSSMAT file into TBUF.'
        ELSE
         CALL QUIT(
     &  'GET_HINF_R: Error, EOFLABEL not found in BSSMAT file!')
        END IF

        IMAT = 1
        IBUF = 1
        DO I=1,NFSYM
        IF (NORB(I).GT.0) THEN
C      ... do TBUF. TMAT. TBUF^+ = VMAT (MOAO transf.)
         IBUFA = IBUF + NPSH(I)*NORB(I) 

         IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'GET_H_INF_R (T2): R in ON "k" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TMAT(IMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &     'GET_H_INF_R (T2): Transf.mtx FPM'//
     &     ' "k"->RKB "theta" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &     'GET_H_INF_R (T2): Transf.mtx FPM'//
     &     ' "k"->RKB "theta" basis, ep-block ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TBUF(IBUFA),NPSH(I),NESH(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL QTRANS('MOAO','S',D0,NESH(I),NESH(I),
     &               NESH(I),NESH(I),
     &            VMAT(IMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TMAT(IMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

         IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'GET_H_INF_R (T2): R in ON "theta" (transformed'//
     &     ' from "k") basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(VMAT(IMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                   NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         IMAT = IMAT + NESH(I)*NESH(I)*NZ
         IBUF = IBUF + NORB(I)*NORB(I)*NZ

        ENDIF
        ENDDO
        CALL QUIT('R(SL): "k"->"theta" transform')
      ENDIF

C ...  generate the Omega+ matrix  
C ... first get TBUF = 1 + R^+.R matrix
      IBUFE = 1
      IEIG =  1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF   = NESH(I)
       N2REF = NREF*NREF
       N2EQ = N2REF*NZ

C ... first make the unit matrix TBUF(IBUFE)
       CALL DUNIT2(TBUF(IBUFE),NREF,NREF,NREF,NZ)

C ... then do TBUF = TBUF + R^+.R
       CALL QGEMM(NREF,NREF,NREF,D1, 
     &    'H','N',IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ,
     &    'N','N',IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ,
     &     D1,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C ... do print out 1 + R^+.R matrix
       IF (IPRINT.GE.7) THEN
         CALL HEADER( 'GET_HINF_R (2): 1 + R^+.R matrix, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... do diagonalize  1 + R^+.R, TBUF
       CALL QDIAG(NZ,NREF,TBUF(IBUFE),NREF,NREF,
     &                 EIG(IEIG),1,VMAT(IBUFE),NREF,
     &                 NREF,WORK(KFREE),LFREE,IERR)
       IF (IERR.NE.0) CALL QUIT('GET_HINF_R: QDIAG failed!')

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &  'GET_HINF_R (2): Eigenvalues of (1 + R+ * R)',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         DO J=1,NREF
           WRITE(LUPRI,'(I5,F25.10)') J,EIG(IEIG-1+J)
         ENDDO
       ENDIF

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'GET_HINF_R (2): Eigenvectors of (1 + R+ * R)',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... do copy eigenvectors into TBUF which is no longer needed
       CALL DCOPY(N2EQ,VMAT(IBUFE),1,TBUF(IBUFE),1)

C ... do rescale eigenvectors
       JEVEC = IBUFE
       DO IZ = 1, NZ
       DO J = 1,NREF
         FACTOR = D1 / DSQRT(EIG(IEIG-1+J))
C              ....  d^-1/2 * U(IZ)
         CALL DSCAL(NREF,FACTOR,TBUF(JEVEC),1)
         JEVEC = JEVEC + NREF 
       ENDDO
       ENDDO

C     now generate Omega+ = (X+ * X)^(-1/2) into TMAT (is reused)
C       .... U * d^-1/2 * U+ 
       CALL QGEMM(NREF,NREF,NREF,D1,  
     &   'N','N',IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ,
     &   'H','N',IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ,
     &    D0,IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ)

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'GET_HINF_R (2): Omega+ =  (1 + R+ * R)^(-1/2), TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUFE = IBUFE + N2EQ
       IEIG  = IEIG + NREF

      ENDIF
      ENDDO

C ... read R into VMAT
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
       CALL RFBSSMAT('R_MTX   ',.TRUE.,'GET_HINF_FPM',VMAT,
     &               ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'R matrix (in ON basis) was read from BSSMAT into VMAT.'
       ENDIF

C ...  make VMAT -> TBUF = (1 R)  !
      IBUFE  = 1
      IBUF2E = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN

       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GET_HINF_R (3): R matrix read from BSSMAT, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C       now insert unit matrix in bottom half of X:
C    ... do make (1 R) into TBUF (R is in VMAT)
       JOFF =  IBUFE + NREF
       DO IZ = 1,NZ
        DO J = 1,NREF
C      ... do copy NREF into TBUF
         IVMAT = IBUFE + (IZ-1)*N2REF + (J-1)*NREF
         IBUF  = IBUF2E + (IZ-1)*2*N2REF + (J-1)*NREF2
         CALL DCOPY(NREF,VMAT(IVMAT),1,TBUF(IBUF),1)
C      ... add 1 and O
         IBUF2 = IBUF + NREF
         CALL DZERO(TBUF(IBUF2),NREF)
         IF (IZ.EQ.1) TBUF(IBUF2+J-1) = D1
        END DO
       END DO

       IBUFE = IBUFE + N2EQ
       IBUF2E = IBUF2E + 2*N2EQ

      ENDIF
      ENDDO

      IBUFE  = 1
      IBUF2E = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &'GET_HINF_R (4): (1 + R) matrix, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C      and multiply it to form X Omega+ = (1  R)t * Omega+:
       CALL QGEMM(NREF2,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),TBUF(IBUF2E),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(IBUFE),NREF ,NREF,NZ,
     &          D0,IPQTOQ(1,0),VMAT(IBUF2E),NREF2,NREF,NZ)


       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GET_HINF_R (4): Final U1= (1 R)t*Omega+ picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUFE = IBUFE + N2EQ
       IBUF2E = IBUF2E + 2*N2EQ

      ENDIF
      ENDDO

C ... TEST: Do U+.U into TMAT
!     IBUF2E=1
!     IBUFE=1
!     IBUF4E=1
!     DO I=1,NFSYM
!     IF (NORB(I).GT.0) THEN
!       NREF = NESH(I)
!       NREF2 = 2*NESH(I)
!-----------------------------------------------------------------------
!     U^{+}(M,K)*U(K,N) =  U^+(NREF,2*NREF)*U(2*NREF,NREF)
!-----------------------------------------------------------------------
!                 M     N      K
!      CALL QGEMM(NREF,NREF,2*NREF,D1, 
!    &    'H','N',IPQTOQ(1,0),
!    &     VMAT(IBUF2E),2*NREF,NREF,NZ,
!    &    'N','N',IPQTOQ(1,0),
!    &     VMAT(IBUF2E),2*NREF,NREF,NZ,
!    &     D0,IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ)

!      IF (IPRINT.GE.7) THEN
!        CALL HEADER(
!    &   'GET_HINF_R : U1+ * U1 multiplied '//
!    &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
!        WRITE(LUPRI,'(3X,A,I1,A,I1)')
!    &  '*** Fermion corep ',I,'/',NFSYM
!        CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
!    &             NZ,IPQTOQ(1,0),LUPRI)
!      ENDIF

!         IBUFE =  IBUFE + NESH(I)*NESH(I)*NZ
!        IBUF2E = IBUF2E + 2*NESH(I)*NESH(I)*NZ
!        IBUF4E = IBUF4E + NREF2*NREF2*NZ

!     ENDIF
!     ENDDO

C =================================================================================
C ... save 4c ON-> 2c ON picture change transformation matrix into the file
C =================================================================================
       ISIZE = 2*NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + 2*NESH(2)*NESH(2)*NZ
       CALL W2BSSMAT(1,.TRUE.,'U1_ONBAS','GET_HINF_R',
     &               VMAT,ISIZE,LUBSS,IPRINT)

       IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R: U1_ONBAS (picture change transf.matrix '//
     & '4c ONB->2c ONB), VMAT, is written to the BSSMAT'//
     & ' file, EOFLABEL is renewed.'
       ENDIF

C =============================================================================
C ... do the  V U1("k") (Ad-1 V+) transformation to get the U1 in "theta" basis
C
C    V-left is full transf.matrix, U1 is ep.e matrix and V-right is only (ee) block !
C
C        We do not get U1' matrix identical with its RKB counterpart !
C  (checked with Octave for 1s(identical) and 2s basis (not identical) )
C =============================================================================
      FP2TH2=.FALSE.
      IF (FP2TH2.AND.FPB) THEN
       ISIZE = NORB(1)*NORB(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
       CALL RFBSSMAT('FPMAT_TH',.TRUE.,'GET_HINF_R',TBUF,
     &      ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R (T3): 4c "k"-> 4c "theta" transf. matrix FPMAT_TH '//
     & ' was read from the BSSMAT file into TBUF.'
        ENDIF

        IBUF = 1
        IBUF2E = 1
        JOFFESH = 0
        DO I=1,NFSYM
        IF (NORB(I).GT.0) THEN
            
         NREF = NESH(I)
         NREF2 = NORB(I)

         IF(IPRINT.GE.5) THEN
          CALL HEADER(
     &   'GET_HINF_R (T3): Final U1= (1 R)t*Omega+ picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &  'GET_HINF_R (T3): 4c "k"-> 4c "theta" transf.'//
     & ' matrix FPMAT_TH  in TBUF',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

             write(lupri,*) 'FACKIN:',
     &     ( FACKIN(JESH+JOFFESH),JESH=1,NESH(I))
C         .... renormalize only ee block ????
            DO IZ = 1,NZ
C           write(lupri,*) 'IZ=',IZ
            DO IESH = 1, NREF ! run over columns...
             IBUFX = IBUF+(NREF*NREF2)+
     &         ((IESH-1)*NREF2)+(IZ-1)*NORB(I)*NORB(I)
            DO JESH = 1, NREF2
              TBUFOLD = TBUF(IBUFX+JESH-1)
             write(lupri,*) 'elem:',IBUFX+JESH-1,
     &               TBUFOLD,TBUF(IBUFX+JESH-1),FACKIN(IESH+JOFFESH)
            ENDDO
            ENDDO 
            ENDDO 

        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &  'GET_HINF_R (T3): 4c "k"-> 4c "theta" transf.'//
     & ' matrix FPMAT_TH in TBUF after RENORMALIZATION',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C   DO the V U1 V+ into TMAT

        CALL QTRANS('MOAO','S',D0,NREF2,NREF,
     &                            NREF2,NREF,
     &            TMAT(IBUF2E),NREF2,NREF,NZ,IPQTOQ(1,0),
     &            VMAT(IBUF2E),NREF2,NREF,NZ,IPQTOQ(1,0),
     &            TBUF(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

        IF(IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GET_HINF_R (T3): Final U1 picture change transformation '//
     &   'matrix TRANSFORMED BACK from ON "k" to RKB basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

        ENDIF

         JOFFESH = JOFFESH + NESH(I)
         IBUF = IBUF + NORB(I)*NORB(I)*NZ
         IBUF2E = IBUF2E + NESH(I)*NORB(I)*NZ

        ENDIF
        ENDDO
        call quit('GET_HINF_R: Tra. U1 transformation to "k"')
      ENDIF
C     ... end of blind "if"

C ...  prepare the resulting Hamiltonian in ON RKB/"k" basis
C   after the picture change transformation
C ... read the Dirac Hamiltonian in the ON RKB/"k" basis into  TMAT
      ISIZE = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
      CALL RFBSSMAT('H1_ONBAS',.TRUE.,'GET_HINF_R',
     &               TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'GET_HINF_R: H1_ONBAS read into TMAT.'
      ENDIF

      IBUF = 1
      IE   = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GET_HINF_R (5): Final U1= (1 R)t*Omega+ picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(ITMAT),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GET_HINF_R (5): H1 Dirac bare nucleus matrix '//
     &   'in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUF),NREF2,NREF2,NREF2,NREF2,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QTRANS('AOMO','S',D0,NREF2,NREF2,NREF,NREF,
     &        TMAT(IBUF),NREF2,NREF2,NZ ,IPQTOQ(1,0),
     &        TBUF(IE),  NREF ,NREF ,NZ ,IPQTOQ(1,0),
     &        VMAT(ITMAT),NREF2,NREF ,NZ, IPQTOQ(1,0),
     &        VMAT(ITMAT),NREF2,NREF ,NZ, IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.4) THEN
           CALL HEADER(
     &   'GET_HINF_R (5): Obtained U1+ H1 U1 = H2c_infinite'//
     &   ' in ON RKB/"k" basis, TBUF (output)',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(IE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
        END IF

        IE  = IE + N2EQ
        ITMAT = ITMAT + 2*N2EQ 
        IBUF = IBUF + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO

C ... save H2c_ON basis into BSSMAT
      ISIZE = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
      CALL W2BSSMAT(1,.TRUE.,'H2C_ONBA','GET_HINF_R',
     &              TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GET_HINF_R:  H2c_infinite (in ON basis) written'//
     &  ' to BSSMAT, EOFLABEL renewed '
      ENDIF

      CALL QEXIT('GET_HINF_R')

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck h2cfinal */
      SUBROUTINE H2CFINAL(TBUF,TMAT,VMAT,FACKIN,IBEIG,WRITEH2C,
     &                    SAVESO1,SPFR,
     &                    SPFR_BEG,GETPCTM,FPB,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C     PURPOSE:
C    ==========
C        Do the transformation of the H2c(ee) from ON MO(ee) basis
C    to the final SA-AO(LL) basis.
C
C    On input:
C   -----------
C    TBUF -  H2c(ee) in ON MO basis RKB/"k"
C    WRITEH2C - variable that decides to assign the label of H_DK 1 el.
C              transformed elements.
C    IBEIG - array to be refilled with boson info (from BSSMAT)
C
C    On output: H2c in SA-AO is written to the file...
C   ------------
C
C    Written:  M.Ilias, Strasbourg 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TMAT(*),VMAT(*),TBUF(*),FACKIN(*),IBEIG(*),WORK(LWORK)
      CHARACTER*8 MYLABEL
      LOGICAL SAVESO1, SPFR,SPFR_BEG, 
     &          WRITEH2C, GETPCTM, TRAFPRKB, FPB

#include "memint.h"

C ... assumes that H2C_ONBA is in FP "k" basis !
      TRAFPRKB = .FALSE.
C .... NEED TO INTERCHANGE  VARIABLES !!!
      IF (TRAFPRKB.AND.FPB) THEN
C  ... read the transformation matrix "k" -> "theta"
      ISIZE = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
      CALL RFBSSMAT('FPMAT_TH',.TRUE.,'H2CFINAL',
     &               TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'H2CFINAL: FPMAT_TH read into TBUF.'
      ENDIF

      IVMAT = 1
      IBUF = 1
      ITMAT = 1
      JOFFESH = 0
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

         IBUFA = IBUF + (NORB(I)*NESH(I)) + NESH(I)

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in ON "k" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(VMAT(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

           CALL HEADER(
     &     'H2CFINAL: Transf.mtx FPM "k"->RKB "theta" basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

           write(lupri,*) 'FACKIN:',
     &     ( FACKIN(JESH+JOFFESH),JESH=1,NESH(I) )

C         .... renormalize only the "ee" block !
            DO IZ = 1,NZ
C           write(lupri,*) 'IZ=',IZ
            DO IESH = 1, NESH(I)
             IBUFX = IBUFA+((IESH-1)*NORB(I))+(IZ-1)*NORB(I)*NORB(I)
            DO JESH = 1, NESH(I) 
             TBUF(IBUFX+JESH-1) = 
     &       TBUF(IBUFX+JESH-1)/DSQRT(FACKIN(IESH+JOFFESH))
C            write(lupri,*) 'elem:',ICX,IBUFX+ICX,TBUF(IBUFX+ICX)
            ENDDO
            ENDDO 
            ENDDO 

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: RENORMALIZED FPM "k"->RKB "theta" basis'//
     &     ' ee block',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IBUFA),NESH(I),NESH(I),NORB(I),NORB(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

        CALL QTRANS('MOAO','S',D0,NESH(I),NESH(I),
     &               NESH(I),NESH(I),
     &            TMAT(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            VMAT(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            TBUF(IBUFA),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in ON "theta" (transformed'//
     &     ' from "k") basis ',-1)
           WRITE(LUPRI,'(/3X,A,I1,A,I1/)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TMAT(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                 NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

        JOFFESH = JOFFESH  + NESH(I)
        IVMAT = IVMAT + NESH(I)*NESH(I)*NZ
        IBUF = IBUF + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO 
      ENDIF ! end of test sequence

      IF (SPFR.OR.SAVESO1.OR.FPB.OR.SPFR_BEG) THEN
CMI ... take care of IBEIG, because it is used further 
      ISIZE =  NORB(1)            
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)
      CALL RFBSSMAT('IBOS_ONB',.TRUE.,'H2CFINAL',IBEIG,
     &                ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(2X,A)')                        
     &   'H2CFINAL: IBEIG was read from BSSMAT' 
         IF (IPRINT.GE.5) THEN
           CALL PRINT_IBEIG(IBEIG,'H2CFINAL')
         ENDIF
      ENDIF
      ENDIF

C    ... do add SO1  terms in MO into TBUF !!!
      IF (SAVESO1) THEN
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
       CALL RFBSSMAT('HSO1_ONB',.TRUE.,'H2CFINAL',TMAT,
     &               ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.2) THEN
          WRITE(LUPRI,'(/2X,A)')
     &    'H2CFINAL: HSO1_ONB read into TMAT.'
       ENDIF
      
       IF (IPRINT.GE.5) THEN
         CALL PRINT_IBEIG(IBEIG,'H2CFINAL')
       ENDIF

        IMAT = 1
        DO I=1,NFSYM
        IF (NESH(I).GT.0) THEN
         N2EQ = NESH(I)*NESH(I)*NZ
         IF (IPRINT.GE.5) THEN
           CALL HEADER(
     & 'H2CFINAL: H1_SO1 in ON RKB/"k" basis, TMAT',-1)
           WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TMAT(IMAT),NESH(I),NESH(I),
     &     NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

CMI      ... eliminate spin-orbit terms 
         IF (SPFR.AND..NOT.SPFR_BEG) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in RKB/|k> - before elimination '//
     &     ' of spin_orbit terms "AFTER"',-1)
           WRITE(LUPRI,'(3X,A,I2,A,I1)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
         ENDIF

          CALL SPF_HK(TBUF(IMAT),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &           NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

          IF (IPRINT.GE.5) THEN
          CALL HEADER(
     & 'H2CFINAL: H2c (SPIN-FREE!) in the ON RKB/"k"basis, TBUF',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TBUF(IMAT),NESH(I),NESH(I),
     &          NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C        ...  add H_SO1 into TBUF !!!
         CALL DAXPY(N2EQ,D1,TMAT(IMAT),1,TBUF(IMAT),1)

         IMAT = IMAT + N2EQ

        ENDIF
        ENDDO
      ENDIF

C  ... read transformation matrix 4c AO_LL -> 4c ON(ee)
      ISIZE = NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NESH(2)*NZ
      CALL RFBSSMAT('TM2CAOFP',.TRUE.,'H2CFINAL',TMAT,
     &               ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN 
         WRITE(LUPRI,'(/2X,A)')
     &   'H2CFINAL: TM2CAOFP read into TMAT.'
      ENDIF

C ....  do the transformation into VMAT
      IVMAT = 1
      IBUF = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NESH(I).GT.0) THEN
C  ********************************************************************
C    ... if wished, eliminate SO terms from H2FP(IH2FP) "at the end"
C  ********************************************************************
        IF (SPFR.AND..NOT.SPFR_BEG.AND..NOT.SAVESO1) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &     'H2CFINAL: H2c in RKB/|k> ON basis - before elimination '//
     &     ' of spin_orbit terms "AFTER"',-1)
           WRITE(LUPRI,'(3X,A,I2,A,I1)')
     &      '*** Fermion corep ',I,' /',NFSYM
            CALL PRQMAT(TBUF(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)

          ENDIF

          CALL SPF_HK(TBUF(IVMAT),.TRUE.,IBEIG(IORB(I)+1+NPSH(I)),
     &            IBEIG(IORB(I)+1+NPSH(I)),
     &           NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

          IF (IPRINT.GE.5) THEN
            CALL HEADER(
     &     'H2CFINAL: H2c in ON RKB/"k" - eliminated'//
     &     ' spin_orbit terms "AFTER"',-1)
            WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     &     '*** Fermion corep ',I,' /',NFSYM
           CALL PRQMAT(TBUF(IVMAT),NESH(I),NESH(I),NESH(I),NESH(I),
     &                  NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

        ENDIF
 
        IF (IPRINT.GE.4) THEN
          CALL HEADER(
     & 'H2CFINAL: H2c in ON basis',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TBUF(IVMAT),NESH(I),NESH(I),
     &          NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     & 'H2CFINAL: Transf.matrix 2c ON MO RKB/"k" -> 2c SA-AO, TMAT',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(TMAT(IBUF),NFBAS(I,0),NESH(I),
     &          NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C       Row dimension of H2AO: only L [NFBAS(I,1)] or L+S [NFBAS(I,0)]?
C       currently we use L+S, to use only L we should delete small Xhi
C       after MAKE_H2C /5-Aug-2002 hjaaj

        CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &               NESH(I),NESH(I),
     &            VMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &            TBUF(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &            TMAT(IBUF),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            TMAT(IBUF),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)

        IF (IPRINT.GE.4) THEN
          CALL HEADER(
     & 'H2CFINAL: FINAL H2c in SA-AO basis "xhi_LL", VMAT',-1)
          WRITE(LUPRI,'(/3X,A,I2,A,I1/)')
     & '*** Fermion corep ',I,' /',NFSYM
          CALL PRQMAT(VMAT(1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &          NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IVMAT = IVMAT + NESH(I)*NESH(I)*NZ
        IBUF = IBUF + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO

C ====================================================================
C  ... write the final H2c/FDAO_DK SA-AO(LL) atomic matrix elements
C     into the file BSSMAT...
C ====================================================================
      IF (WRITEH2C) THEN
        MYLABEL='H1AO_DK '
      ELSE
        MYLABEL='FDAO_DK '
      ENDIF
      CALL W2BSSMAT(1,.TRUE.,MYLABEL,'H2CFINAL',VMAT,
     &                N2BBASXQ,LUBSS,IPRINT)      
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A,A)')
     & 'H2CFINAL: ',MYLABEL,' AO_LL elements were written'//
     & ' into the BSSMAT file, EOFLABEL renewed.'
      ENDIF

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck getpctmat */      
      SUBROUTINE GETPCTMAT(TBUF,TMAT,VMAT,FPB,USEDF,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:  
C    ========  
C
C    On input/output - nothing, all quantities are read/written from/into BSSMAT.
C                                                         
C    FPB
C    USEDF                                                         
C                                                         
C         ...  construct the picture change transformation matrixes
C   -----------------------------------------------------------------------------
C
C       H4c AO ->  H2c AO
C
C    sfB: H4c AO ->  H4c "k"  &   H4c"k" -> H2c AO
C
C    sfE: H4c AO ->  H2c "k"  &   H2c"k" -> H2c AO
C
C    Written:  M.Ilias, Strasbourg 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!     use x2c_utils, only:
!    &    print_x2cmat 
#include "implicit.h"
#include "priunit.h"
      PARAMETER ( D1 = 1.00D00 , D0 = 0.00D00 )
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      DIMENSION TMAT(*),VMAT(*),TBUF(*),WORK(LWORK)
      LOGICAL  FNDLAB, FPB, USEDF
      CHARACTER*8 LABL

      CALL QENTER('GETPCTMAT')
#include "memint.h"

      ISIZE = NORB(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NESH(2)*NZ
      CALL RFBSSMAT('U1_ONBAS',.TRUE.,'GETPCTMAT',
     &               TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     & 'GETPCTMAT: U1_ONBAS (4c ON->2c ON) read into TBUF.'
      ENDIF

      IF (FPB) THEN
C      ... free-particle ("k") basis
        LABL = 'FPMAT4C '
      ELSE IF (.NOT.FPB) THEN
C        ... 
        LABL = 'SL_TM_4C'
      ELSE
        CALL QUIT('GETPCTMAT: blind branch for LABL !!!')
      ENDIF

      ISIZE2 = NFBAS(1,0)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + NFBAS(2,0)*NORB(2)*NZ
      IF (ISIZE2.NE.N2TMT)
     &     CALL QUIT('GETPCTMAT: ISIZE2.NE.N2TMT !!') 
      CALL RFBSSMAT(LABL,.TRUE.,'GETPCTMAT',
     &               VMAT,ISIZE2,LUBSS,IPRINT)

      IBUF = 1
      IVMAT = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & 'ON RKB/"k" H4c -> H2c, U1_ONBAS ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NORB(I),NESH(I),
     &                    NORB(I),NESH(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
          CALL HEADER(
     &'GETPCTMAT: 4c ON basis (RKB/"k") in SA-AO basis "xhi"',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NORB(I),
     &        NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C ...   construct picture change transf. matrix,  VMAT*TBUF = TMAT
      CALL QGEMM(NFBAS(I,0),NESH(I),NORB(I),D1, 
     &       'N','N',IPQTOQ(1,0),
     &        VMAT(IVMAT),NFBAS(I,0),NORB(I),NZ,
     &       'N','N',IPQTOQ(1,0),
     &        TBUF(IBUF),NORB(I),NESH(I),NZ,
     &        D0,IPQTOQ(1,0),
     &        TMAT(ITMAT),NFBAS(I,0),NESH(I),NZ)

        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & 'H4c SA-AO "xhi" -> H2c ON RKB/"k", TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),
     &    NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IBUF = IBUF + NORB(I)*NESH(I)*NZ
        IVMAT = IVMAT + NFBAS(I,0)*NORB(I)*NZ
        ITMAT = ITMAT + NFBAS(I,0)*NESH(I)*NZ
 
      ENDIF
      ENDDO

C   ==========================================================
C    ... write TMAT into the BSSMAT to be reread again later
C   ==========================================================
      ISIZE = NESH(1)*NFBAS(1,0)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NFBAS(2,0)*NZ)
      CALL W2BSSMAT(1,.TRUE.,'T4CAO2CK','GETPCTM',
     &               TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTM: PCTM 4c SA-AO "xhi" -> H2c ON RKB/"k"'//
     &  ' written to the BSSMAT file (T4CAO2CK), EOFLABEL renewed.'
      ENDIF

      ISIZE2 = NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + (NFBAS(2,0)*NESH(2)*NZ)
      CALL RFBSSMAT('TM2CAOFP',.TRUE.,'GETPCTMAT',
     &               VMAT,ISIZE2,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/2X,A)')
     & 'GETPCTMAT: TM2CAOFP (H2c "xhi_L" -> H2c(ee) ON) read into VMAT.'
      ENDIF

C ======================================================
C   ... construct TBUF.VMAT+ = TMAT (4c ON -> 2c AO)
C ======================================================
      IBUF = 1
      IVMAT = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'GETPCTMAT: Final transf.matrix H2c "xhi_L" - H2c(ee) ON  '//
     &   ' (TM2CAOFP), VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NESH(I),
     &         NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &   'GETPCTMAT: Picture change transf.matrix '//
     & ' H4c in ON -> H2c(ee) in ON, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NORB(I),NESH(I),
     &                    NORB(I),NESH(I),
     &                     NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C  ...  construct picture change transf matrix, TBUF*VMAT+ = TMAT
       CALL QGEMM(NORB(I),NFBAS(I,1),NESH(I),D1, 
     &           'N','N',IPQTOQ(1,0),
     &           TBUF(IBUF),NORB(I),NESH(I),NZ,
     &          'H','N',IPQTOQ(1,0),
     &           VMAT(IVMAT),NFBAS(I,0),NESH(I),NZ,
     &          D0,IPQTOQ(1,0),
     &          TMAT(ITMAT),NORB(I),NFBAS(I,1),NZ)

        IF (IPRINT.GE.5) THEN
         CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & 'H4c ON -> H2c SA-AO "xhi", TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NORB(I),NFBAS(I,1),
     &       NORB(I),NFBAS(I,1),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
        
        IBUF = IBUF + NORB(I)*NESH(I)*NZ
       ITMAT = ITMAT + NFBAS(I,1)*NORB(I)*NZ
       IVMAT = IVMAT + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO

C    ... write TMAT into the BSSMAT
      ISIZE = NORB(1)*NFBAS(1,1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NFBAS(2,1)*NZ
      CALL W2BSSMAT(1,.TRUE.,'T4CK2CAO','GETPCTM',
     &               TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTM: PCTM H4c ON RKB/"k" -> H2c SA-AO "xhi"'//
     &  ' written to the BSSMAT file(T4CK2CAO), EOFLABEL renewed.'
      ENDIF

C   ... reread T4CAO2CK back into TMAT
      ISIZE2 = NFBAS(1,0)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + NFBAS(2,0)*NESH(2)*NZ
      CALL RFBSSMAT('T4CAO2CK',.TRUE.,'GETPCTMAT',
     &               TMAT,ISIZE2,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'GETPCTMAT: T4CAO2CK (H2c "xhi_L" -> H2c(ee) ON) read into TMAT.'
      ENDIF

      IBUF = 1
      IVMAT = 1
      ITMAT = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF (IPRINT.GE.5) THEN

          CALL HEADER(
     & 'GETPCTMAT: Picture change transf.matrix '//
     & 'H4c SA-AO "xhi" -> H2c ON RKB/"k", TMAT',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),
     &    NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'GETPCTMAT: Transf.matrix H2c ON(e) -> H2c "xhi_L" '//
     &   ' (TM2CAOFP), VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NESH(I),
     &         NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C  ...  construct full picture change transfprmatin matrix, TMAT*VMAT+ = TBUF
      CALL QGEMM(NFBAS(I,0),NFBAS(I,1),NESH(I),D1, 
     &          'N','N',IPQTOQ(1,0),
     &           TMAT(ITMAT) ,NFBAS(I,0),NESH(I),NZ,
     &          'H','N',IPQTOQ(1,0),
     &           VMAT(IVMAT),NFBAS(I,0),NESH(I),NZ,
     &          D0,IPQTOQ(1,0),
     &          TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),NZ)

        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     & 'GETPCTMAT: Final picture change transf.matrix '//
     & '  H4c SA-AO "xhi" -> H2c SA-AO "xhi"',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),
     &              NFBAS(I,0),NFBAS(I,1),
     &             NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

!       debug print
!       call print_x2cmat(TBUF(IBUF),NFBAS(I,0),NFBAS(I,1),nz,
!    &                    IPQTOQ(1,0),'bss - pctmat AO l',6)

        IBUF =  IBUF + NFBAS(I,0)*NFBAS(I,1)*NZ
       ITMAT = ITMAT + NFBAS(I,0)*NESH(I)*NZ
       IVMAT = IVMAT + NFBAS(I,0)*NESH(I)*NZ

      ENDIF
      ENDDO

C ************************************************************************
C .. save the final picture change transformation matrix into the file
C ************************************************************************
      ISIZE = NFBAS(1,0)*NFBAS(1,1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NFBAS(2,0)*NFBAS(2,1)*NZ
      CALL W2BSSMAT(1,.TRUE.,'U_PICTRM','GETPCTM',
     &               TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GETPCTMAT: Final picture change transformation matrix '//
     &  ' - U_PICTRM - was written to the BSSMAT file, '//
     &  ' EOFLABEL was renewed.'
      ENDIF

      CALL QEXIT('GETPCTMAT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck cmp_eigval */      
       SUBROUTINE CMP_EIGVAL(VMAT,TMAT,TBUF,EIG,IBEIG,WRITEH2C,SPFR,
     &                       SPFR_BEG,ONLYSO1,SAVESO1,
     &                       USEDF,DO4C2C,START2C,
     &                       IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C   ==========
C    Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians;
C  transform each of them to the ON basis and do the diagonalization
C
C  On input: arrays to be utilized inside, and flags
C
C  On output: VMAT reorganized for H2c                                                    
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0 ) 
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"
#include "dummy.h"
   
      DIMENSION VMAT(*),TMAT(*),
     &          EIG(*),IBEIG(*),TBUF(*),WORK(LWORK)    
      DIMENSION NBORBX(4, 2, 0:2),NPSHX(2)
      LOGICAL SPFR,SPFR_BEG,ONLYSO1,SAVESO1,WRITEH2C,
     &        PURLO,USEDF,DO4C2C,START2C,WANT_LINEAR
      CHARACTER*8 LABELH1,LABELH2,LBL

      CALL QENTER('CMP_EIGVAL')
#include "memint.h"         

C ====================================================================
C   ... do the picture change of FD4C_ONB using 4c ONB -> 2c SA-AO
C         picture change transformation matrix
C ====================================================================
      IF (DO4C2C.AND..NOT.USEDF.AND.WRITEH2C.AND..NOT.START2C) THEN
C ...   read the pictre change matrix into TMAT
C    ... write TMAT into the BSSMAT  
       ISIZE = NORB(1)*NFBAS(1,1)*NZ  
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NFBAS(2,1)*NZ  
       CALL RFBSSMAT('T4CK2CAO',.TRUE.,'CMP_EIGVAL',TMAT,
     &               ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.3) THEN   
        WRITE(LUPRI,'(/2X,A)') 
     &  'CMP_EIGVAL: PCTM H4c ON RKB/"k" -> H2c SA-AO "xhi"'//  
     &  ' read from the BSSMAT file(T4CK2CAO) into TMAT.' 
       ENDIF

C     ... read the FD4C_ONB matrix which was prepared earlier
       ISIZE = NORB(1)*NORB(1)*NZ  
       IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ  
       CALL RFBSSMAT('FD4C_ONB',.TRUE.,'CMP_EIGVAL',VMAT,
     &               ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.3) THEN   
        WRITE(LUPRI,'(/2X,A)') 
     & 'CMP_EIGVAL: 4c Fock-Dirac(FD4C_ONB) in ON RKB/"k" was'//  
     & ' read from the BSSMAT file(FD4C_ONB) into VMAT.' 
       ENDIF

C      =========================================================================
C       ... do the picture change transformation  TBUF = TMAT^+ . VMAT . TMAT
C      =========================================================================
       IVMAT = 1
       ITMAT = 1
       IBUF = 1
       DO I=1,NFSYM
       IF (NORB(I).GT.0) THEN

        IBUF = 1 + I2BASX(I,I)

        CALL QTRANS('AOMO','S',D0,NORB(I),NORB(I),
     &      NFBAS(I,1),NFBAS(I,1),    
     &      VMAT(IVMAT),NORB(I),NORB(I),NZ,IPQTOQ(1,0), 
     &      TBUF(IBUF),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),   
     &      TMAT(ITMAT),NORB(I),NFBAS(I,1),NZ,IPQTOQ(1,0),    
     &      TMAT(ITMAT),NORB(I),NFBAS(I,1),NZ,IPQTOQ(1,0),    
     &      WORK(KFREE),LFREE,IPRINT)     

        IF (IPRINT.GE.5) THEN        

         CALL HEADER( 
     &  'CMP_EIGVAL: Read 4c Fock-Dirac matrix in '//     
     &  '4c ON RKB/FP "k" basis, VMAT',-1)    
         WRITE(LUPRI,'(3X,A,I1,A,I1)')    
     &   '*** Fermion corep ',I,'/',NFSYM   
         CALL PRQMAT(VMAT(IVMAT),NORB(I),NORB(I),    
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)   

         CALL HEADER( 
     &  'CMP_EIGVAL: Read H4c ON RKB/"k" -> H2c SA-AO "xhi" '//     
     &  'picture change transf. mtx, TMAT',-1)    
         WRITE(LUPRI,'(3X,A,I1,A,I1)')    
     &   '*** Fermion corep ',I,'/',NFSYM   
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,1),NORB(I),    
     &         NFBAS(I,1),NORB(I),NZ,IPQTOQ(1,0),LUPRI)   

         CALL HEADER( 
     &  'CMP_EIGVAL: Resulting pctra Fock-Dirac_2c matrix in'//     
     &  ' 2c SA-AO "xhi", TBUF',-1)    
         WRITE(LUPRI,'(3X,A,I1,A,I1)')    
     &   '*** Fermion corep ',I,'/',NFSYM   
         CALL PRQMAT(TBUF(IBUF),NFBAS(I,1),NFBAS(I,1),    
     &         NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)   
        ENDIF                                        

         IVMAT = IVMAT + NORB(I)*NORB(I)*NZ
         ITMAT = ITMAT + NFBAS(I,1)*NORB(I)*NZ
       ENDIF
       ENDDO

C     ... save Fock-Dirac matrix in 2c SA-AO "xhi" into BSSMAT
      CALL W2BSSMAT(1,.TRUE.,'FDAO_DK ','CMP_EIGVAL',
     &               TBUF,N2BBASXQ,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'CMP_EIGVAL: FD2c elements (FDAO_DK) were written'//
     &  ' into the BSSMAT file, EOFLABEL renewed.'
      ENDIF

C  ======== End of handling Fock-Dirac 2c =========
      ENDIF
C  ======== End of handling Fock-Dirac 2c =========

      WANT_LINEAR=.FALSE.
       IF (LINEAR) THEN
!     IF (WANT_LINEAR) THEN
          LABELH1='VMAT2C_L'
      ELSE
          LABELH1='VMAT2C  '
      ENDIF
      CALL RFBSSMAT(LABELH1,.TRUE.,'CMP_EIGVAL',
     &        VMAT,N2BBASX,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/2X,A,A)')                  
     &  'CMP_EIGVAL: Lowdin 2c transf. matrix was read'// 
     &  ' into VMAT from the BSSMAT file,label=',LABELH1           
       ENDIF

C ... read H_2c/FD_2c in SA-AO basis into TBUF
      IF (WRITEH2C) THEN
C      ... transformed  bare nucleus Dirac Hamiltonian
       LABELH2='H1AO_DK '
      ELSE
C      ... transformed converged Fock-Dirac Hamiltonian (both closed/open-shell)
       LABELH2='FDAO_DK '
      ENDIF
 
      CALL RFBSSMAT(LABELH2,.TRUE.,'CMP_EIGVAL',
     &              TBUF,N2BBASXQ,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A,A,A)')
     & 'CMP_EIGVAL: ',LABELH2,' was read from BSSMAT.'
      ENDIF

C   ... transform H2c AO through Lowdin matrix into ON MO ee basis

      IH2MO = 1
      IVMAT = 1
      IH1FP = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN
       IF (IPRINT.GE.4) THEN
C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'CMP_EIGVAL: Final H(2c)_LL in "xhi" SA-AO basis ',-1)
        IF (WRITEH2C) THEN
          WRITE(LUPRI,'(3X,A,/)')
     &    'This is 2c transformed Dirac bare nucleus,'//
     &    ' H1AO_DK.' 
        ELSE
          WRITE(LUPRI,'(3X,A,/)')
     &    'This is 2c transformed Fock-Dirac Hamiltonian,'//
     &    ' FDAO_DK.' 
        ENDIF

        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(1+I2BASX(I,I)),NFBAS(I,1),NFBAS(I,1),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

C    ... print out Lowdin matrix
        CALL HEADER(
     &  'CMP_EIGVAL: Adapted Lowd. transf. matrix (Le)'//
     &  ' for the 2c SCF procedure:',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NESH(I),NFBAS(I,0),
     &              NESH(I),1,IPQTOQ(1,0),LUPRI)

       ENDIF

C   ... transform H2c AO through Lowdin matrix into ON MO ee basis
        CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NESH(I),NESH(I),
     &    TBUF(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    TMAT(IH2MO),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &    VMAT(IVMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &    VMAT(IVMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &        WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &   'CMP_EIGVAL: H2c in (MO_ee) ON Lowdin basis "Theta" to'//
     &   ' be diagonalized, TMAT',-1)
          IF (.NOT.WRITEH2C) THEN
            WRITE(LUPRI,'(3X,A,/)')
     &      'This is FD2C_LMO, FD_2c in Lowdin matrix for subsequent'//
     &      ' postDHF->BSS SCF step' 
          ENDIF
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TMAT(IH2MO),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
C       ... offset for Lowdin matrix
       IVMAT = IVMAT + NFBAS(I,0)*NESH(I)
       IH1FP = IH1FP + NORB(I)*NORB(I)*NZ
      ENDIF
      ENDDO

      IF (.NOT.WRITEH2C) THEN
C ...  save FD_2c in Lowdin MO basis to the BSSMAT file
C ...  to be used at the post-DC-BSS-SCF level in DFSOLV !
       ISIZE = NESH(1)*NESH(1)*NZ
       IF (NFSYM.EQ.2) ISIZE=ISIZE+(NESH(2)*NESH(2)*NZ)
       CALL W2BSSMAT(1,.TRUE.,'FD2C_LMO','CMP_EIGVAL',
     &               TMAT,ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')                  
     &  'CMP_EIGVAL: FD_2c elements in Lowdin MO basis were written'// 
     &  ' into the BSSMAT file (FD2C_LMO), EOFLABEL renewed.'  
       ENDIF

      ENDIF 

C ... WE WANT THE EIGENVECTORS AS WELL .... ADD ONE MORE CYCLE ONLY 
C      FOR THE DIAGONALIZATION...REUSE TBUF for eigenvectors

      IH2MO = 1
      IEIG = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

C      ... diagonalize the H2c in MO(ee) basis !
!       DOJACO is defined in dcbgen.h, only for real matrixes
       IF (DOJACO) THEN
!        write(lupri,*) 'CMP_EIGVAL: RSJACO used for diagonalization' 
         CALL RSJACO(NESH(I),NESH(I),NESH(I),TMAT(IH2MO),
     &                  EIG(IEIG),1,1,0,TBUF(IH2MO))
       ELSE
        CALL QDIAG(NZ,NESH(I),TMAT(IH2MO),NESH(I),NESH(I), 
     &             EIG(IEIG),1,          
     &             TBUF(IH2MO),NESH(I),       
     &             NESH(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 1 !...')
       ENDIF

       IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &   'CMP_EIGVAL: Eigenvectors of the diagonalized H2c'//
     &   ' in ON MO_ee Lowdin basis "Theta": ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH2MO),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IEIG  = IEIG  + NESH(I)
       IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
      ENDIF
      ENDDO

C .... read and diagonalize H1 in ON RKB/FP "k" basis
C ... write the TBUF - H1 in RKB/"k" to the file...
      ISIZE=NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE=ISIZE+NORB(2)*NORB(2)*NZ
C ... get H1_ONB into TBUF
      LABELH1 = 'H1_ONBAS'
      IF (SAVESO1) LABELH1 = 'H1_ONB  '
      CALL RFBSSMAT(LABELH1,.TRUE.,'CMP_EIGVAL',TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
       WRITE(LUPRI,'(/2X,A,A,A)')
     & 'CMP_EIGVAL: ',LABELH1,' read into TBUF.'
      ENDIF

      IH1FP = 1
      IH2MO = 1
      IEIG4  = 1
      IEIG2  = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN
C        ... eliminate spin-free terms utilizing IBEIG(*)...
       IF (SPFR.AND..NOT.SPFR_BEG) THEN

         IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'CMP_EIGVAL: H4c in ON MO basis RKB/FP "k" before '//
     &    'elimination of SO terms (we have 2c spin-free "AFTER") ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH1FP),NORB(I),NORB(I),
     &             NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         CALL SPF_HK(TBUF(IH1FP),.TRUE.,IBEIG(IORB(I)+1),
     &    IBEIG(IORB(I)+1),
     &    NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRINT )

       ENDIF

        IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'CMP_EIGVAL: H4c in ON MO basis RKB/FP "k" before '//
     &    ' the diagonalization ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH1FP),NORB(I),NORB(I),
     &             NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QDIAG(NZ,NORB(I),TBUF(IH1FP),NORB(I),NORB(I), 
     &             TMAT(IEIG4),0,          
     &             DUMMY,NORB(I),       
     &             NORB(I),WORK(KFREE),LFREE,IERR)
        IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 2!')

        IF (IPRINT.GE.2) THEN
         CALL HEADER(
     &'CMP_EIGVAL: Eigenvalues of H4c and H2c '//
     &' Hamiltonians - comparison',-1)
         WRITE(LUPRI,'(8X,A,I1,A,I1/)')
     &  '*** Fermion corep ',I,'/',NFSYM
C ... first positronic states first
        WRITE(LUPRI,'(3X,I5,F25.10)')
     &  ( J-NPSH(I)-1,TMAT(IEIG4+J-1),J=1,NPSH(I))
C ... then continue with electronic states
        DO J = 1, NESH(I)
          VAL2C = EIG(IEIG2+J-1)
          VAL4C = TMAT(IEIG4+NPSH(I)+J-1)
        DIFF =  VAL4C - VAL2C
        WRITE(LUPRI,'(3X,I5,2F25.10,3X,D12.5)') 
     &    J, VAL4C,VAL2C,DIFF
        ENDDO
        write(lupri,*)
        ENDIF

        IEIG4  = IEIG4  + NORB(I)
        IEIG2  = IEIG2  + NESH(I)
        IH2MO = IH2MO + NESH(I)*NESH(I)*NZ 
        IH1FP = IH1FP + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO

C ****************************************************************************
C **
C **    Prepare the Fock-2c matrix in Lowdin MO to be used in restarting 
C **  the 2-component SCF. 
C **
C ****************************************************************************
      IF (DO4C2C.AND..NOT.USEDF.AND.WRITEH2C.AND..NOT.START2C) THEN

C      ... do prepare the FD4C_ONB 
C .... read Lowdin matrix into TMAT
        IF (LINEAR) THEN
          LABELH1='VMAT2C_L'
        ELSE
          LABELH1='VMAT2C  '
        ENDIF

        CALL RFBSSMAT(LABELH1,.TRUE.,'CMP_EIGVAL',
     &               TMAT,N2BBASX,LUBSS,IPRINT)

        IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/2X,A,A)')                  
     &  'CMP_EIGVAL: Lowdin 2c transf. matrix was read'// 
     &  ' into TMAT from the BSSMAT file,label=',LABELH1           
        ENDIF

        CALL RFBSSMAT('FDAO_DK ',.TRUE.,'CMP_EIGVAL',
     &                 VMAT,N2BBASXQ,LUBSS,IPRINT)

        IF (IPRINT.GE.3) THEN
          WRITE(LUPRI,'(/2X,A)')                  
     &  'CMP_EIGVAL: FD2c elements (2c SA-AO basis,FDAO_DK) were read'//
     &  ' into VMAT from the BSSMAT file.'           
        ENDIF
C       ... do the transformation ...
      IBUF  = 1
      ITMAT = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

       IF (IPRINT.GE.5) THEN
C    ... print out FH2c matrix in AO "xhi" basis
        CALL HEADER(
     &'CMP_EIGVAL: FD(2c)_LL in "xhi" SA-AO basis ' ,-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(1+I2BASX(I,I)),NFBAS(I,1),NFBAS(I,1),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)

C    ... print out Lowdin matrix
        CALL HEADER(
     &  'CMP_EIGVAL: Adapted Lowd. transf. matrix (Le)'//
     &  ' for the 2c SCF procedure:',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &        '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),NFBAS(I,0),
     &              NESH(I),1,IPQTOQ(1,0),LUPRI)

       ENDIF

C   ... transform H2c AO through Lowdin matrix into ON MO ee basis
       CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &               NESH(I),NESH(I),
     &    VMAT(1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    TBUF(IBUF),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &    TMAT(ITMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &    TMAT(ITMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &     WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &  'CMP_EIGVAL: FD2c in (MO_ee) ON Lowdin basis "Theta" to'//
     &  ' be diagonalized, TBUF',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(IBUF),NESH(I),NESH(I),
     &        NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUF = IBUF + NESH(I)*NESH(I)*NZ 
C       ... offset for Lowdin matrix
       ITMAT = ITMAT + NFBAS(I,0)*NESH(I)
      ENDIF
      ENDDO

       ISIZE = NESH(1)*NESH(1)*NZ 
       IF (NFSYM.EQ.2) ISIZE=ISIZE+(NESH(2)*NESH(2)*NZ)
       CALL W2BSSMAT(1,.TRUE.,'FD2C_LMO','CMP_EIGVAL',
     &                TBUF,ISIZE,LUBSS,IPRINT)
       IF (IPRINT.GE.3) THEN
         WRITE(LUPRI,'(/2X,A)')                  
     &  'CMP_EIGVAL: FD_2c elements in Lowdin MO basis were written'// 
     &  ' into the BSSMAT file, EOFLABEL renewed.'           
       ENDIF
C        ... do diagonalize the FD2C_LMO ...
      IEIG  = 1
      IH2MO = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

       CALL QDIAG(NZ,NESH(I),TBUF(IH2MO),NESH(I),NESH(I), 
     &             EIG(IEIG),1,          
     &             TBUF(IH2MO),NESH(I),       
     &             NESH(I),WORK(KFREE),LFREE,IERR)
       IF (IERR.NE.0) CALL QUIT('CMP_EIGVAL: FAILED QDIAG 1 !...')

       IF (IPRINT.GE.15) THEN
          CALL HEADER(
     &   'CMP_EIGVAL: Eigenvectors of the diagonalized FD2c'//
     &   ' in ON MO_ee Lowdin basis "Theta": ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(TBUF(IH2MO),NESH(I),NESH(I),
     &             NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IF (IPRINT.GE.2) THEN
          CALL HEADER(
     &   'CMP_EIGVAL: Eigenvalues of the diagonalized FD2c'//
     &   ' in ON MO_ee Lowdin basis "Theta": ',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         DO J=1,NESH(I)
          WRITE(LUPRI,'(3X,I5,F25.10)') 
     &    J,EIG(IEIG+J-1)
         ENDDO
       ENDIF
       IEIG  = IEIG  + NESH(I)
       IH2MO = IH2MO + (NESH(I)*NESH(I)*NZ) 
      ENDIF
      ENDDO
      ENDIF

      CALL QEXIT('CMP_EIGVAL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hdk2_fpm */      
      SUBROUTINE GET_HDK2_FPM(TBUF,TMAT,VMAT,EIG,IBEIG,
     &                     DK2SO1,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE:
C    ========
C
C  Do comparison of the eigenvalues of the Dirac and H2c Hamiltonians
C   trasnform each to the ON basis and do the diagonalization
C
C  On input:  TBUF H1 in FP "k" basis
C             EIG - free particle eigenvalues needed for denominators
C             IBEIG - boson irreps if spin-free, or DK2SO1
C                                                         
C  On output:  TBUF 
C             the desired h1_DK2 Hamiltonian in ON FP "k" basis is saved into BSSMAT
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJAaJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER(D1 = 1.0D0, D0 = 0.0D0, D2 = 2.0D0, DM1 = -1.0D0 ) 
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dcbgen.h"   
#include "dgroup.h"
   
      DIMENSION TBUF(*),TMAT(*),VMAT(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL SPFR_BEG, DK2SO1

      CALL QENTER('GET_HDK2_FPM')
#include "memint.h"         

C       ... second order Douglas-Kroll Hamiltonian, or DK2
C
C        TBUF contains the H_1 matrix in orthonormal 4-comp. |k> basis
C             (= free particle basis in the specifed basis set);
C             first NPSH positron solutions, then NESH electron solutions
C
C        H_DKH(2,ee) = H_1(ee)
C                    + 0.5 * ( R_DK2+(ep) H_1(pe) + H_1(ep) R_DK2(pe) )
C   where:
C
C        R_DK2(pe)(i,j) = H_1(pe)(i,j)/(eig_e(i) - eig_p(j))
C     
C     (taken from Barysz&Sadlej, JMS(Theochem) 573(2001) 181, eq. 46 )    
C

      IF (IPRINT.GE.2) THEN
        CALL HEADER(
     &'GET_HDK2_FPM: DKH2 Hamiltonian to be generated',-1)
      ENDIF

      CVAL2 = CVAL*CVAL
      JOFT1 = 1
      JBUF  = 1
      DO  I = 1,NFSYM  
      IF(NORB(I).GT.0) THEN

        N2T    = NORB(I)*NORB(I)
        N2ESHI = NESH(I)*NESH(I)
C
C         1) generate R_DK2(pe) in TBUF(pe) block:
C
         IOFT  = JOFT1 + NORB(I)*NPSH(I) - 1
         JR2PE = IOFT  + 1

         DO IE = NPSH(I) + 1,  NORB(I)
           DO IP = 1, NPSH(I)
C ... eig_p and eig_e are separated by 2c^2; use rather average
C ... eig(p) = -eig(e) - 2c^2
C         .... Calculate the average value
             EIG_P = 
     & ((-EIG(NORB(I)-IP+1+IORB(I))-(D2*CVAL2)) + 
     &   EIG(IP+IORB(I)))/D2

             EIG_E =
     & ( EIG(IE+IORB(I)) + 
     &   (-EIG(NORB(I)-IE+1+IORB(I))-(D2*CVAL2)) )/D2

             FAC = DM1/ (EIG_P - EIG_E)
             DO IZ = 0,NZ-1
               TBUF(IOFT+IP+IZ*N2T) =
     &         TBUF(IOFT+IP+IZ*N2T) * FAC
             END DO
           END DO
           IOFT = IOFT + NORB(I)
         END DO

          JH1EP = JOFT1 + NPSH(I)
 
C ----------------------------------------------------------------------------y
C  If desired, eliminate SO1 terms from both H1_ep and R_pe - DK2SO1 approach
C -----------------------------------------------------------------------------
          IF (DK2SO1) THEN
            IF (IPRINT.GE.2) THEN
              WRITE(LUPRI,'(3X,A)')
     & 'GET_HDK2_FPM: DKH2sfBSO1 option - elimination of spin-orbit'//
     & ' terms from H1(ep) and R(pe) blocks'
            ENDIF
           CALL SPF_HK( TBUF(JH1EP),.TRUE.,IBEIG(IORB(I)+NPSH(I)+1),
     &                 IBEIG(IORB(I)+1),
     &                 NESH(I),NPSH(I),NORB(I),NORB(I),NZ,IPRINT )
           CALL SPF_HK( TBUF(JR2PE),.TRUE.,IBEIG(IORB(I)+1),
     &                 IBEIG(IORB(I)+1+NPSH(I)),
     &                 NPSH(I),NESH(I),NORB(I),NORB(I),NZ,IPRINT )
          ENDIF

          IF (IPRINT.GE.6) THEN
           CALL HEADER(
     &    'GET_HDK2_FPM: 1. order wave function R_DK2(pe)'//
     &    ' in "|k>" basis',-1)
           WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM

           CALL PRQMAT(TBUF(JR2PE),NPSH(I),NESH(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
           CALL HEADER(
     &'GET_HDK2_FPM: H(ep) part of H1 in "|k>" basis',-1)
           WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(TBUF(JH1EP),NESH(I),NPSH(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C
C         2) Generate  "X" = H_1(ep) R_DK2(pe) into VMAT(JBUF)
C

          CALL QGEMM(NESH(I),NESH(I),NPSH(I),D1,
     &          'N','N',IPQTOQ(1,0),TBUF(JH1EP),NORB(I),NORB(I),NZ,
     &          'N','N',IPQTOQ(1,0),TBUF(JR2PE),NORB(I),NORB(I),NZ,
     &          D0,IPQTOQ(1,0),VMAT(JBUF), NESH(I),NESH(I),NZ)
C
C         3) Generate  X = 0.5 ( X + X(dagger)) 
C
          CALL QGETHE(VMAT(JBUF),NESH(I),NZ,NESH(I),NESH(I))

          IF (IPRINT.GE.7) THEN
            CALL HEADER(
     &'GET_HDK2_FPM: 2. order correction to H1 in "|k_e>" basis',-1)
           WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
            CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),NESH(I),
     &                  NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
*  Add to first order (ee-blocks!): H_1(ee) =  H_1(ee) + X

          DO IZ = 1,NZ
            IBUF  = JBUF + (IZ-1)*N2ESHI
            JH1EE = JOFT1 + NPSH(I)*NORB(I) + NPSH(I)
     &            + (IZ-1)*NORB(I)*NORB(I)
            DO J = 1, NESH(I)
              CALL DAXPY(NESH(I),D1,TBUF(JH1EE),1,VMAT(IBUF),1)
              JH1EE = JH1EE + NORB(I)
              IBUF  = IBUF  + NESH(I)
            END DO
          END DO

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &    'GET_HDK2_FPM: Full H2c_DKH2(or DKH2sfBSO1) Hamiltonian'//
     &    ' in ON FP "k" basis, VMAT',-1)
          WRITE(LUPRI,'(5X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),NESH(I),
     &                  NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF
C         ... update the offsets...
          JBUF  = JBUF  + NESH(I)*NESH(I)*NZ
          JOFT1 = JOFT1 + N2T*NZ
      ENDIF
      ENDDO

C ... do save Hamiltonian matrix (VMAT) into the file BSSMAT
      ISIZE = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
      CALL DCOPY(ISIZE,VMAT,1,TBUF,1)
      CALL W2BSSMAT(1,.TRUE.,'H2C_ONBA','GET_HDK2_FPM',
     &              TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'GET_HDK2_FPM: H2c_DKH2 Hamiltonian matrix (in "k" basis) is '//
     & ' written to the BSSMAT file, EOFLABEL renewed.'
      ENDIF

      CALL QEXIT('GET_HDK2_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*          
C  /* Deck get_hdk1_fpm */      
      SUBROUTINE GET_HDK1_FPM(TBUF,VMAT,IBEIG,ONLYSO1,
     &                         IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Do construct the simplest H1(ee) Hamiltonian in "k" basis
C    ========
C
C  On input:  TBUF -  H1 in "k" basis
C             IBEIG - boson irreps if pure spin-orbit integrals (ONLYSO1=.true.)
C             ONLYSO1 - flag for removing all scalar terms from Hamiltonian 
C
C    NB: ONLYSO1 is active only for two-step (free particle) type of transformation,
C                where we can get DKH Hamiltonians.
C                                                         
C  On output: TBUF - H2c in "k" basis
C           ...  the desired DK1 Hamiltonian in "k" basis is written to the file BSSMAT
C                                                         
C  Rewritten by Miro ILIAS, Strasburg, 2005 (after HJJ,MI, Odense,2002)
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"                                       
#include "priunit.h"   
      PARAMETER( DM1= -1.0D0 )
#include "dcbgen.h"  
#include "dcbbas.h"  
#include "dcborb.h"   
#include "dgroup.h"   

      DIMENSION TBUF(*),VMAT(*),IBEIG(*),WORK(LWORK)
      LOGICAL ONLYSO1

      CALL QENTER('GET_HDK1_FPM')
#include "memint.h"         
  
C  ... do copy the ee-part of the full H1 Hamiltonian into VMAT
      JBUF  = 1
      JOFT1 = 1
      DO I = 1, NFSYM
      IF (NESH(I).GT.0) THEN

        N2ESHI = NESH(I)*NESH(I)
        N2T    = NORB(I)*NORB(I)

        DO IZ = 1, NZ
         IBUF  = JBUF + (IZ-1)*N2ESHI
         JH1EE = JOFT1 + (NPSH(I)*NORB(I)) + NPSH(I)
     &          + (IZ-1)*N2T
         DO J = 1, NESH(I)
           CALL DCOPY(NESH(I),TBUF(JH1EE),1,VMAT(IBUF),1)
           JH1EE = JH1EE + NORB(I)
           IBUF  = IBUF  + NESH(I)
         ENDDO
        ENDDO

C      ... if desired, remove all scalar terms to have only SO1 integrals...
C      ======================================================================
        IF (ONLYSO1) THEN

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &  'GET_HDK1_FPM: Full H_DK1 matrix in ON FP "k" basis '//
     &  ' before eliminating spin-free terms (ONLYSO1=.true.)',-1)
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL SPF_HK(VMAT(JBUF),.FALSE.,IBEIG(IORB(I)+1+NPSH(I)),
     &               IBEIG(IORB(I)+1+NPSH(I)),
     &              NESH(I),NESH(I),NESH(I),NESH(I),NZ,IPRINT)

        ENDIF

          IF (IPRINT.GE.5) THEN
           CALL HEADER(
     &  'GET_HDK1_FPM: Final H_DK1-SO1 matrix in "k" basis ',-1)
           IF (ONLYSO1) WRITE(LUPRI,'(3X,A)')
     &  '...spin-free terms are deleted!'
           WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
           CALL PRQMAT(VMAT(JBUF),NESH(I),NESH(I),
     &         NESH(I),NESH(I), NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

C       ... update the offsets...
        JBUF  = JBUF  + N2ESHI*NZ
        JOFT1 = JOFT1 + N2T*NZ
      ENDIF 
      ENDDO

C ... do write the H1ee Hamiltonian into the file...
      ISIZE = NESH(1)*NESH(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NESH(2)*NESH(2)*NZ
      CALL DCOPY(ISIZE,VMAT,1,TBUF,1)
      
      CALL W2BSSMAT(1,.TRUE.,'H2C_ONBA','GET_HDK1_FPM',
     &              TBUF,ISIZE,LUBSS,IPRINT)

      IF (IPRINT.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'GET_HDK1_FPM: H2c_DKH1 Hamiltonian matrix (in "k" basis) is '//
     &  ' written to the BSSMAT file, EOFLABEL renewed.'
      ENDIF

      CALL QEXIT('GET_HDK1_FPM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck slsort */
      SUBROUTINE SLSORT(TMAT,TBUF,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE:
C ----------
C  Resort the entering RKB transformation matrix for the one-step (Trond's) BSS
C  transformation. Called from MODHAM, if needed (SLGROUP).
C  Routine is called when DC_RKB Hamiltonian is applied, for post-DCSCF transformations.
C
C  NOTE: Routine does the job only if TMAT contains
C        the "standard" symmetry (ie no linear symmetry) 
C        transformation matrix and when NSYM.GT.1.
C        If NSYM=1 (NZ=4) then only copy TMAT into TBUF, empty run.
C
C  Called from:  MODHAM/dirone.F
C
C  On input: 
C ===========
C  TMAT - the transformation matrix (is unchanged!)
C   which columns are to be resorted into order 
C   S1S2..S(nbsym) L1L2...L(nbsym) from the original order S1 L1 S2 L2 ...S(nbsym)L(nbsym)
C
C  On output:
C ============= 
C  TBUF - contains reordered S1S2..S(nbsym) L1L2...L(nbsym) RKB transformation matrix
C
C  Written by Miro Ilias (MI), July2005/Strasbourg
C  Last modifications:  MI, April2006/Strasbourg
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dgroup.h"    
#include "dcbbas.h"   
#include "dcborb.h"   

      DIMENSION TMAT(*),TBUF(*),WORK(LWORK)

      CALL QENTER('SLSORT')
#include "memint.h"

      NSYM = 4/NZ

      IF (IPRINT.GE.3) THEN
        CALL HEADER('***  Output from SLSORT  ***',-1)
        write(lupri,*)
     &  'Going to resort the RKB transformation matrix to'//
     &  ' the S-L columns order; NZ=',
     &   NZ,' NSYM=(4/NZ)=',NSYM,' and NFSYM=',NFSYM
        write(lupri,*) '============================================='//
     &  '======================================'
        DO I=1,NFSYM
          write(lupri,*) '   .....ifs=',I,'/',NFSYM
          write(lupri,*) '  S:',(NBORB(ISYM,I,1),ISYM=1,NSYM)
          write(lupri,*) '  L:',(NBORB(ISYM,I,2),ISYM=1,NSYM)
          write(lupri,*) 'S+L:',(NBORB(ISYM,I,0),ISYM=1,NSYM)
        ENDDO
      ENDIF

C     ... this check does not work here !!!
C     CALL MEMCHK('in SLSORT',WORK,1)

      IPOS_NFSYM = 1

      DO I = 1, NFSYM                      
      IF (NORB(I).GT.0) THEN                  
      DO IZ = 1, NZ                                      

       IPOS_IZ_NFSYM = IPOS_NFSYM + ((IZ-1)*NORB(I)*NFBAS(I,0))

       IPOS_LS = 0
       IPOS_S  = 0
       DO ISYM = 1, NSYM 
       IF (NBORB(ISYM,I,0).GT.0) THEN    
          NBS = NFBAS(I,0)*NBORB(ISYM,I,1)               
          IPOS_LS1 = IPOS_LS + IPOS_IZ_NFSYM   
          IPOS_S1  = IPOS_S  + IPOS_IZ_NFSYM 
          CALL DCOPY(NBS,TMAT(IPOS_LS1),1,TBUF(IPOS_S1),1)
C         ... skip the SL(isym) block
          IPOS_LS = IPOS_LS + (NFBAS(I,0)*NBORB(ISYM,I,0))
C         .... accumulate (and skip) S(isym) blocks
          IPOS_S  = IPOS_S  + (NFBAS(I,0)*NBORB(ISYM,I,1)) 
       ENDIF                                      
       ENDDO                                      
                                                    
C .. then do copy L columns            
       IPOS_LS = 0
       IPOS_L  = IPOS_S
       DO ISYM = 1, NSYM    
       IF (NBORB(ISYM,I,0).GT.0) THEN   
         NBS = NFBAS(I,0)*NBORB(ISYM,I,2)         
C        ... skip the previous S(isym) block
         IPOS_LS1 = IPOS_LS + NFBAS(I,0)*NBORB(ISYM,I,1) + IPOS_IZ_NFSYM
         IPOS_L1  = IPOS_L  + IPOS_IZ_NFSYM
         CALL DCOPY(NBS,TMAT(IPOS_LS1),1,TBUF(IPOS_L1),1)
C        ... skip the SL(isym) block
         IPOS_LS = IPOS_LS + (NFBAS(I,0)*NBORB(ISYM,I,0))
C         .... accumulate (and skip) L(isym) blocks
         IPOS_L  = IPOS_L + (NFBAS(I,0)*NBORB(ISYM,I,2)) 
       ENDIF                                                       
       ENDDO 

       ENDDO  ! of IZ

       IF (IPRINT.GE.17) THEN                                       
         CALL HEADER(
     &   'SLSORT: "SL" sorted final transformation matrix'//  
     &   ' 4c SA-AO -> RKB_MO, TBUF (output) ',-1)           
          WRITE(LUPRI,'(/A,I1,A,I1)')                         
     &   '* Fermion ircop no.',I,'/',NFSYM                     
          CALL PRQMAT(TBUF(IPOS_NFSYM),NFBAS(I,0),NORB(I),           
     &          NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)   
       ENDIF                                   

       IPOS_NFSYM = IPOS_NFSYM + (NFBAS(I,0)*NORB(I)*NZ)   

      ENDIF                                   
      ENDDO  

      WRITE(LUPRI,'(A)')
     &' * Applied SL-regrouping on AO2MO tranf.matrix in SLSORT.'

      CALL QEXIT('SLSORT')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      SUBROUTINE FDTR2C(VMAT,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE:
C ----------
C    Do the BSS transformation of 4comp Fock-Dirac AO matrixes to the
C  two-component form.
C
C   Fock-Dirac AO matrixes are read from the DFFCK1,
C  DFFCK2 files and in the latter case (2el.Fock matrixes) are written back.
C
C   All other necessary stuff (picture change transformation matrixes etc) 
C  is read from the BSSMAT file and was prepared before.
C
C  On input : only arrays (VMAT serves only as an array)
C
C  Called from MAKE_H2C
C
C  Written by Miro Ilias, July, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dcbgen.h"     
#include "dgroup.h"    
#include "dcbbas.h"   
#include "dcborb.h"   
#include "dcbdhf.h"   
#include "dcbham.h"    

      DIMENSION VMAT(*),IBEIG(*),WORK(LWORK)
      LOGICAL FPB, KEEPH1, BSS_SAVE
      CHARACTER*8 LABL

      CALL QENTER('FDTR2C')
#include "memint.h"

      IF (IPRHAM.GE.2) THEN
        CALL HEADER('***  Output from FDTR2C  ***',-1)
      ENDIF

      FPB = .NOT.NOPRTR

      CALL MEMGET('REAL',KF1,N2BBASXQ      ,WORK,KFREE,LFREE)

      ISIZE = NFBAS(1,0)*NFBAS(1,1)*NZ              
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NFBAS(2,0)*NFBAS(2,1)*NZ) 
      CALL RFBSSMAT('U_PICTRM',.TRUE.,'FDTR2C',VMAT,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.3) THEN
          WRITE(LUPRI,'(/2X,A)')                        
     &  'FDTR2C: Picture change transformation matrix (4c AO->2c AO)'//
     &  ' was read from the BSSMAT file.'     
      ENDIF
    
      IF (IPRHAM.GE.5) THEN
C ... print out symmetry blocked picture change transformation matrix                             
      IPCTM1 = 1
      WRITE(LUPRI,'(/A)')                                          
     &'FDTR2C: SYMMETRY BLOCKED picture change transformation'//
     &' matrix(4c SA AO-2c SA AO), VMAT:'
      DO I1 = 1, NFSYM                
        WRITE(LUPRI,'(/,3X,A,I1,A,I1)')                                 
     &  '*** Fermion corep ',I1,'/',NFSYM                         
        CALL PRQMAT(VMAT(IPCTM1),NFBAS(I1,0),NFBAS(I1,1),              
     &     NFBAS(I1,0),NFBAS(I1,1),NZ,IPQTOQ(1,0),LUPRI)              
        IPCTM1 = IPCTM1 + NFBAS(I1,0)*NFBAS(I1,1)*NZ                 
      ENDDO               
      ENDIF

C    .. allocate arrays for picture change transformation
      CALL MEMGET('REAL',KF2,N2BBASXQ*NFMAT,WORK,KFREE,LFREE)

      CALL OPNFIL(LUFCK1,'DFFCK1','OLD','FDTR2C')
      CALL REAFCK(LUFCK1,WORK(KF2),.TRUE.,1)
      CLOSE(LUFCK1,STATUS='KEEP')

      IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     &  'FDTR2C: READ (from DFFCK1) 1el converged Fock-Dirac matrix',-1)
        CALL PRQMAT(WORK(KF2),NTBAS(0),NTBAS(0),
     &         NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
      ENDIF

C .. do the picture change transformation first of Fock1=hD(1) (KF2) into (KF1)
C ===============================================================================
      IPCTM1 = 1
      DO I=1, NFSYM
         CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                NFBAS(I,1),NFBAS(I,1),
     &   WORK(KF2+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &   WORK(KF1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &   VMAT(IPCTM1),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &   VMAT(IPCTM1),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &   WORK(KFREE),LFREE,IPRHAM)

        IF (IPRHAM.GE.5) THEN
          CALL HEADER(
     &'FDTR2C: H2c_1el(=h_inf) constructed (pctm.)'//
     & ' again in SA-AO basis "xhi_L", symm.blocked ',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(KF1+I2BASX(I,I)),NFBAS(I,1),NFBAS(I,1),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF
C       ... offset for pict.change transf.matrix 4c AO "xhi" -> 2c  AO "xhi"
         IPCTM1 = IPCTM1 + NFBAS(I,0)*NFBAS(I,1)*NZ
      ENDDO

      CALL W2BSSMAT(1,.TRUE.,'FD1AO_DK','FDTR2C',
     &              WORK(KF1),N2BBASXQ,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN 
        WRITE(LUPRI,'(/2X,A)')                         
     &  'FDTR2C: FD1AO_DK AO_LL elements were written'//       
     &  ' into the BSSMAT file, EOFLABEL was renewed.'
      ENDIF

C ==================================================================
C    ... continue with 2el Fock-AO matrixes ... do the picture
C    change transformations of them if .NOT.CONT2C
C ==================================================================
      IF (.NOT.ONESYS.AND..NOT.CONT2C) THEN
        CALL OPNFIL(LUFCK2,'DFFCK2','OLD','PSISCF')
C          ... read NFMAT 2el Fock matrixes !
        CALL REAFCK(LUFCK2,WORK(KF2),.TRUE.,NFMAT)

C       ... run over all 2-electron Fock-Dirac AO matrixes ...
        DO IFMAT = 1, NFMAT

         IF2 = KF2 + (IFMAT-1)*N2BBASXQ

         IF (IPRHAM.GE.5) THEN
          CALL HEADER('FDTR2C: read full 2-el. Fock matrix',-1)
          WRITE(LUPRI,'(3X,A,I2,A,I2)')
     &    '>>> Fock matrix no. ',IFMAT,'/',NFMAT
          CALL PRQMAT(WORK(IF2),NTBAS(0),NTBAS(0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

         IPCTM1 = 1
         DO I=1, NFSYM
C         ... do the picture chnage transf; WORK(KFOCK2) -> WORK(KF2C2EL)
          CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                     NFBAS(I,1),NFBAS(I,1),
     &    WORK(IF2+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    WORK(KF1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    VMAT(IPCTM1),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &    VMAT(IPCTM1),NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &    WORK(KFREE),LFREE,IPRHAM)

        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'FDTR2C: 2l F2c=u2c constructed (pctm)'//
     &   ' again in SA-AO basis "xhi_L" ',-1)
         WRITE(LUPRI,'(3X,A,I2,A,I2)')
     &   '>>> Fock matrix no. ',IFMAT,'/',NFMAT
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KF1+I2BASX(I,I)),NFBAS(I,1),NFBAS(I,1),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
C       ... offset for pict.change transf.matrix 4c AO "xhi" -> 2c  AO "xhi"
         IPCTM1 = IPCTM1 + NFBAS(I,0)*NFBAS(I,1)*NZ
        ENDDO
C        ... copy the picture change transformed 2el Fock matrix back into array...
         CALL DCOPY(N2BBASXQ,WORK(KF1),1,WORK(IF2),1)
        ENDDO
C       ... after the picture change transformation write all 2-electron
C        Fock matrixes back into the file DFFCK2
         CALL WRIFCK(LUFCK2,WORK(KF2),NFMAT)
         CLOSE(LUFCK2,STATUS='KEEP')
         IF (IPRHAM.GE.2) THEN
         write(lupri,*)
     &'FDTR2C: Picture change transformed (using FD_4c) '//
     &' 2-electron Fock matrix(-es) were written to the DFFCK2 file!'
         ENDIF
      ENDIF
C     ... release the largest allocated array is it is no longer needed....
      CALL MEMREL('KF2.FDTR2C',WORK,KWORK,KF2,KFREE,LFREE)
C  ===============================================================
C
C       ...  prepare h_2c(+u2c_SO) Hamiltonian (when .NOT.SPINFR!)
C
C      u2c_SO comes from transformed converged Fock-Dirac Ham.
C
C  ===============================================================
      IF (.NOT.SPINFR) THEN
C     ... read bosonic symmetries on orbitals
      ISIZE =  NORB(1)            
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)
      CALL RFBSSMAT('IBOS_ONB',.TRUE.,'FDTR2C',IBEIG,
     &               ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.3) THEN
         WRITE(LUPRI,'(/2X,A)')                        
     &   'FDTR2C: IBEIG was read from BSSMAT'
         IF (IRPINT.GE.5) THEN
           CALL PRINT_IBEIG(IBEIG,'FDTR2C')
         ENDIF
       ENDIF

      ISIZE = NORB(1)*NORB(1)*NZ      
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ) 
      CALL MEMGET('REAL',KFDK, ISIZE,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KFDK1,ISIZE,WORK,KFREE,LFREE)

      ISIZE = NORB(1)*NFBAS(1,0)*NZ      
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NFBAS(2,0)*NZ 
      CALL MEMGET('REAL',KFPTM,ISIZE,WORK,KFREE,LFREE)

C       ... read the 1el Fock matrix 
C     ========================================================
      BSS_SAVE = BSS
      BSS = .FALSE.
      CALL ONEFCK(WORK(KF1),IPRHAM,WORK(KFREE),LFREE)
      BSS = BSS_SAVE

C     ... read the RKB/"k" transformation matrix from BSSMAT into VMAT ...
      IF (FPB) THEN
        LABL = 'FPMAT4C '
      ELSE IF (.NOT.FPB) THEN
        LABL = 'SL_TM_4C'
      ELSE
         CALL QUIT('FPMAT4C: wrong label !')
      ENDIF

      ISIZE = (NFBAS(1,0)*NORB(1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NFBAS(2,0)*NORB(2)*NZ)
      CALL RFBSSMAT(LABL,.TRUE.,'FDTR2C',
     &                VMAT,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/2X,A,A,A)')
     & 'FDTR2C:',LABL,' - AO2MO transformation matrix '//
     1 ' was read from BSSMAT.'
      ENDIF

C =======================================================================================
C    ... read converged Fock-Dirac matrix in the "k" basis from the BSSMAT file ....
C    ... read the SPECIAL converged Fock-Dirac matrix in the "k" basis from the BSSMAT file ....
C =======================================================================================
      IF (I2CHAM.EQ.6) THEN
        LABL = 'FD4CXONB'
      ELSE
        LABL = 'FD4C_ONB'
      ENDIF
      ISIZE = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ)
      CALL RFBSSMAT(LABL,.TRUE.,'FDTR2C',WORK(KFDK1),
     &               ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.3) THEN
        WRITE(LUPRI,'(/2X,A,A,A)')
     & 'FDTR2C: ',LABL,' was read from BSSMAT.'
      ENDIF

C ===============================================================
C   ... transform  the 1-el Dirac operator into the "k" basis
C ===============================================================
       IMAT = 1
      IFDK  = KFDK
      IFDK1 = KFDK1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'FDTR2C: FD_4c matrix in '//
     &   '4c ON RKB/"k" basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(IFDK1),NORB(I),NORB(I),
     &       NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'Read the ON RKB / FP "k" matrix,'//
     &   ' VMAT, in 4c SA-AO basis "xhi" ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IMAT),NFBAS(I,0),NORB(I),
     &       NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'Read Dirac bare nucleus matrix (from subr.ONEFCK)'//
     &   ' in 4c SO-AO basis "xhi" ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KF1+I2BASX(I,I)),NFBAS(I,0),NFBAS(I,0),
     &                NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NORB(I),NORB(I),
     &       WORK(KF1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &       WORK(IFDK),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &       VMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       VMAT(IMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'FDTR2C: Obtained Dirac H1 matrix in '//
     &  '4c ON "k"/RKB basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(IFDK),NORB(I),NORB(I),
     &     NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

         IMAT = IMAT + NORB(I)*NFBAS(I,0)*NZ
        IFDK  = IFDK + NORB(I)*NORB(I)*NZ
       IFDK1  = IFDK1 + NORB(I)*NORB(I)*NZ

      ENDIF
      ENDDO 

C    ----------------------------------------------------------
C       ... get U1=F1-H1 matrix in the ON-MO basis ...
C    ----------------------------------------------------------
       ISIZE = (NORB(1)*NORB(1)*NZ)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ)
       CALL DAXPY(ISIZE,DM1,WORK(KFDK),1,WORK(KFDK1),1)

       IFDK  = KFDK
       IFDK1 = KFDK1
       DO I=1,NFSYM
       IF (NORB(I).GT.0) THEN

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'FDTR2C: Obtained U1=F1-H1 matrix in '//
     &  '4c "k"/ON basis, TBUF',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &    '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(IFDK1),NORB(I),NORB(I),NORB(I),
     &                NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

C  =======================================================
C         ... remove SPIN-FREE TERMS FROM U1 ...
C  =======================================================
      CALL SPF_HK(WORK(IFDK1),.FALSE.,
     &      IBEIG(IORB(I)+1),IBEIG(IORB(I)+1),
     &      NORB(I),NORB(I),NORB(I),NORB(I),NZ,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'FDTR2C: Obtained U1_SO (after elim.of SF terms) matrix in '//
     &  '4c "k"/ON basis',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(IFDK1),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

          IFDK  = IFDK  + NORB(I)*NORB(I)*NZ
         IFDK1  = IFDK1 + NORB(I)*NORB(I)*NZ

       ENDIF
       ENDDO 

      KEEPH1=.FALSE.
      IF (.NOT.KEEPH1) THEN
C ... add U1_SO to the H1 
        ISIZE2 = NORB(1)*NORB(1)*NZ
        IF (NFSYM.EQ.2) ISIZE2 = ISIZE2 + NORB(2)*NORB(2)*NZ
        CALL DAXPY(ISIZE2,D1,WORK(KFDK1),1,WORK(KFDK),1)
CMI      .. DM1 give surprisingly results for B 3s1p atom ...
      ENDIF

C  =======================================================================
C  ... transform F1(=H1+U1_SO) from 4c "k" to 2c SA-AO basis "xhi" ...
C  =======================================================================
C   ... read the picture change transformation matrix 4c "k" -> 2c AO "xhi"
      ISIZE =  NORB(1)*NFBAS(1,1)*NZ                             
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NFBAS(2,1)*NZ       
      CALL RFBSSMAT('T4CK2CAO',.TRUE.,'FDTR2C',
     &              VMAT,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.3) THEN
      WRITE(LUPRI,'(/2X,A)')                      
     &'FDTR2C: T4CK2CAO (4c "k" -> 2c AO "xhi") was read into VMAT'
      ENDIF

      CALL DZERO(WORK(KF1),N2BBASXQ)

      IPCTM1 = 1
      IFDK = KFDK
      DO I=1,NFSYM

         IF (IPRHAM.GE.6) THEN
         CALL HEADER(
     &  'FDTR2C: The F1=H1+U1_SO Hamiltonain in '//
     &  '4c "k"/RKB basis',-1)
          WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &     '*** Fermion corep ',I,'/',NFSYM
          CALL PRQMAT(WORK(IFDK),NORB(I),NORB(I),NORB(I),
     &                 NORB(I),NZ,IPQTOQ(1,0),LUPRI)
C        ... must be the same as in GETPCTM routine ...
         CALL HEADER(
     & 'FDTR2C: Picture change transf.matrix '//
     & 'H4c ON -> H2c SA-AO "xhi", VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     & '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IPCTM1),NORB(I),NFBAS(I,1),
     &       NORB(I),NFBAS(I,1),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

C         ... do the picture change transf U^+.F1.U
          CALL QTRANS('AOMO','S',D0,
     &              NORB(I),NORB(I),
     &           NFBAS(I,1),NFBAS(I,1),
     &    WORK(IFDK),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &    WORK(KF1+I2BASX(I,I)),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &    VMAT(IPCTM1),NORB(I),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &    VMAT(IPCTM1),NORB(I),NFBAS(I,1),NZ,IPQTOQ(1,0),
     &    WORK(KFREE),LFREE,IPRHAM)

        IF(IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'FDTR2C: Obtained desired f=h+u_so matrix in '//
     &  '2c SA-AO basis "xhi".',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(WORK(KF1+I2BASX(I,I)),NFBAS(I,0),
     &          NFBAS(I,0),NTBAS(0),NTBAS(0),
     &          NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        IPCTM1 = IPCTM1 + NORB(I)*NFBAS(I,1)*NZ
         IFDK  = IFDK  + NORB(I)*NORB(I)*NZ

      ENDDO

C      ... write the F1=H1+U1_SO 2c into the BSSMAT file
C      to be used later in ONEFCK
      CALL W2BSSMAT(1,.TRUE.,'FD2C_USO','FDTR2C',
     &            WORK(KF1),N2BBASXQ,LUBSS,IPRHAM)
      IF (IPRHAM.GE.3) THEN
        WRITE(LUPRI,'(/2X,A)')                         
     &  'FDTR2C: FD2c(h1+u1_SO) SA-AO elements were written '//
     &  ' into the BSSMAT file (FD2C_USO), the EOFLABEL renewed.'
      ENDIF

C     =====   end of SPINFR =====
      ENDIF

C     ..  release all the allocated memmory in this subroutine at once ...!
      CALL MEMREL('FDTR2C',WORK,KWORK,KWORK,KFREE,LFREE)
      CALL QEXIT('FDTR2C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck fullpctm */
      SUBROUTINE FULLPCTM(VMAT,TMAT,TBUF,EIG,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C  PURPOSE: Generate the full picture change transformation matrix
C           using R (in BSSMAT), do compare it 
C 
C  On input - only arrays
C
C  On output - all is written into the BSSMAT file
C
C  Written by Miro Ilias, Strasbourg, 2005
C 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(THRNULL=1.0D-7,D0=0.00D00,D1=1.00D00,DM1=-1.00D00)
#include "dcbgen.h"     
#include "dgroup.h"    
#include "dcbbas.h"   
#include "dcborb.h"   

      DIMENSION VMAT(*),TMAT(*),TBUF(*),EIG(*),WORK(LWORK)

      CALL QENTER('FULLPCTM')
#include "memint.h"

C     ... read R into TMAT
      ISIZE = (NESH(1)*NESH(1)*NZ)                                  
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NESH(2)*NZ)        
      CALL RFBSSMAT('R_MTX   ',.TRUE.,'FULLPCTM',
     &               TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')                                      
     & 'FULLPCTM: R matrix (in ON basis) was read'//
     & ' from BSSMAT into TMAT.' 
      ENDIF

C   ... do generate the Omega- matrix
C ... do generate the Omega+ matrix  
C ... first get TBUF = 1 + R.R^+ matrix
      IBUFE = 1
      IEIG =  1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF   = NESH(I)
       N2REF = NREF*NREF
       N2EQ = N2REF*NZ

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &  'FULLPCTM: R matrix (read from BSSMAT), TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... first make the unit matrix TBUF(IBUFE)
       CALL DUNIT2(TBUF(IBUFE),NREF,NREF,NREF,NZ)

C ... then do TBUF = TBUF + R * R^+
       CALL QGEMM(NREF,NREF,NREF,D1, 
     &     'N','N',IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ,
     &     'H','N',IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ,
     &    D1,IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ)

C ... do print out 1 + R^+.R matrix
       IF (IPRINT.GE.7) THEN
         CALL HEADER( 'FULLPCTM: 1 + R.R^+ matrix, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... do diagonalize 1 + R.R^+, TBUF
       CALL QDIAG(NZ,NREF,TBUF(IBUFE),NREF,NREF,
     &                 EIG(IEIG),1,VMAT(IBUFE),NREF,
     &                 NREF,WORK(KFREE),LFREE,IERR)
       IF (IERR.NE.0) CALL QUIT('FULLPCTM: QDIAG failed!')

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &  'FULLPCTM: Eigenvalues of (1 + R+ * R)',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         DO J=1,NREF
           WRITE(LUPRI,'(I5,F25.10)') J,EIG(IEIG-1+J)
         ENDDO
       ENDIF

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'FULLPCTM: Eigenvectors of (1 + R+ * R)',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &              NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... do copy eigenvectors into TBUF which is no longer needed
       CALL DCOPY(N2EQ,VMAT(IBUFE),1,TBUF(IBUFE),1)

C ... do rescale eigenvectors
       JEVEC = IBUFE
       DO IZ = 1, NZ
       DO J = 1,NREF
         FACTOR = D1 / DSQRT(EIG(IEIG-1+J))
C              ....  d^-1/2 * U(IZ)
         CALL DSCAL(NREF,FACTOR,TBUF(JEVEC),1)
         JEVEC = JEVEC + NREF 
       ENDDO
       ENDDO

C     now generate Omega- = - (X+ * X)^(-1/2) into TMAT (is reused)
C       .... U * d^-1/2 * U+ 
       CALL QGEMM(NREF,NREF,NREF,DM1,  
     &   'N','N',IPQTOQ(1,0),TBUF(IBUFE),NREF,NREF,NZ,
     &   'H','N',IPQTOQ(1,0),VMAT(IBUFE),NREF,NREF,NZ,
     &    D0,IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ)

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'FULLPCTM: Omega- =  (1 + R * R^+)^(-1/2), TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUFE = IBUFE + N2EQ
       IEIG  = IEIG + NREF

      ENDIF
      ENDDO

C     ... read R into VMAT
      ISIZE = (NESH(1)*NESH(1)*NZ)                                  
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NESH(2)*NZ)        
      CALL RFBSSMAT('R_MTX   ',.TRUE.,'FULLPCTM',
     &                VMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.3) THEN   
       WRITE(LUPRI,'(/2X,A)')                                      
     & 'FULLPCTM: R matrix (in ON basis) was read'//
     & ' from BSSMAT into VMAT.' 
      ENDIF

C ... do make VMAT -> TBUF = [1 R]^t  !
      IBUFE  = 1
      IBUF2E = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN

       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.7) THEN
        CALL HEADER(
     &  'FULLPCTM: R matrix read from BSSMAT, VMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C      ... make -R^+ matrix into VMAT...
       CALL QAHM(VMAT(IBUFE),NREF,NZ,NREF,NREF)

       IF (IPRINT.GE.7) THEN
        CALL HEADER(
     &  'FULLPCTM: -R^+ matrix after QAHM, VMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IBUFE),NREF,NREF,NREF,NREF,
     &                 NZ,IPQTOQ(1,0),LUPRI)
       ENDIF


C       now insert unit matrix in UPPER (not bottom) half of X:
C    ... do make [1 -R^+]^t into TBUF (R is in VMAT)
       JOFF =  IBUFE + NREF
       DO IZ = 1, NZ
        DO J = 1, NREF
C      ... do copy NREF into TBUF
         IVMAT = IBUFE + (IZ-1)*N2REF + (J-1)*NREF
         IBUF  = IBUF2E + (IZ-1)*2*N2REF + (J-1)*NREF2
         CALL DCOPY(NREF,VMAT(IVMAT),1,TBUF(IBUF+NREF),1)
C      ... add 1 and O
         IBUF2 = IBUF 
         CALL DZERO(TBUF(IBUF2),NREF)
         IF (IZ.EQ.1) TBUF(IBUF2+J-1) = D1
        END DO
       END DO

       IBUFE = IBUFE + N2EQ
       IBUF2E = IBUF2E + 2*N2EQ

      ENDIF
      ENDDO

      IBUFE  = 1
      IBUF2E = 1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN

       NREF  = NESH(I)
       N2REF = NESH(I)*NESH(I)
       NREF2 = 2*NESH(I)
       N2EQ = NREF*NREF*NZ

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &'  FULLPCTM: [1 (-R^+)]^t matrix, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C      and multiply it to form X Omega- = (1  R)t * Omega-:
       CALL QGEMM(NREF2,NREF,NREF,D1,
     &          'N','N',IPQTOQ(1,0),TBUF(IBUF2E),NREF2,NREF,NZ,
     &          'N','N',IPQTOQ(1,0),TMAT(IBUFE),NREF ,NREF,NZ,
     &          D0,IPQTOQ(1,0),VMAT(IBUF2E),NREF2,NREF,NZ)


       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'FULLPCTM: Final U1-= (1 -R^+)t*Omega- picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUF2E),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUFE = IBUFE + N2EQ
       IBUF2E = IBUF2E + 2*N2EQ

      ENDIF
      ENDDO

C  ... do test UNIT MATRIX:  U1- . U1- ^ +
      IF (IPRINT.GE.7) THEN

      IBUF2E=1
      IBUFE=1
      IBUF4E=1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
        NREF = NESH(I)
        NREF2 = 2*NESH(I)

       CALL QGEMM(NREF,NREF,2*NREF,D1, 
     &    'H','N',IPQTOQ(1,0),
     &     VMAT(IBUF2E),2*NREF,NREF,NZ,
     &    'N','N',IPQTOQ(1,0),
     &     VMAT(IBUF2E),2*NREF,NREF,NZ,
     &     D0,IPQTOQ(1,0),TMAT(IBUFE),NREF,NREF,NZ)

       IF (IPRINT.GE.7) THEN
        CALL HEADER(
     &  'FULLPCTM : U1-^+ * U1- multiplied '//
     &  'transformation matrix in ON RKB/"k" basis, TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(IBUFE),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL TEST_UNIT_MTX(TMAT(IBUFE),NREF,NREF,NREF,NZ)

        IBUFE =  IBUFE + (NREF*NREF*NZ)
       IBUF2E = IBUF2E + (NREF2*NREF*NZ)
       IBUF4E = IBUF4E + (NREF2*NREF2*NZ)

      ENDIF
      ENDDO

      ENDIF

C   ... do construct the full picture change transformation matrix !!! ...
C   ... read U1+ pctm into TMAT
      ISIZE = (2*NESH(1)*NESH(1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (2*NESH(2)*NESH(2)*NZ)
      CALL RFBSSMAT('U1_ONBAS',.TRUE.,'FULLPCTM',
     &               TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'FULLPCTM: U1_ONBAS (+)  matrix (in ON basis)'//
     & ' was read from the BSSMAT file into TMAT.'
      ENDIF

C ===========================================================================
C ... reconstruct the FULL picture change transformation matrix into TMAT
C ===========================================================================
      CALL DCOPY(ISIZE,TMAT(1),1,TBUF(1),1)
      IBUF = 1
      ITMAT = 1
      DO I=1,NFSYM
         NREF = NESH(I)
         NREF2 = 2*NESH(I)
        IF (IPRINT.GE.7) THEN
          CALL HEADER(
     &   'FULLPCTM: U1(+) picture change (read from BSSMAT) '//
     &   'transformation matrix in ON RKB/"k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

          CALL HEADER(
     &   'FULLPCTM: U1(-) picture change (constructed) '//
     &   'transformation matrix in ON RKB/"k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IBUF),NREF2,NREF,NREF2,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)

        ENDIF

        ISIZE1 = (2*NESH(I)*NESH(I))
        ISIZE2 = (NORB(I)*NORB(I))
       DO IZ=1, NZ
         IBPOS = IBUF + ((IZ-1)*ISIZE1)
         ITPOS = ITMAT + ((IZ-1)*ISIZE2)
         ITPOS2 = ITMAT + ((IZ-1)*(ISIZE2))+ISIZE1
         CALL DCOPY(ISIZE1,TBUF(IBPOS),1,TMAT(ITPOS),1)
         CALL DCOPY(ISIZE1,VMAT(IBPOS),1,TMAT(ITPOS2),1)
       ENDDO

       IF (IPRINT.GE.7) THEN
        CALL HEADER(
     &  'FULLPCTM: Final constructed U1 picture change '//
     &  'transformation matrix in ON RKB/"k" basis, TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(ITMAT),NORB(I),NORB(I),NORB(I),NORB(I),
     &           NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       IBUF = IBUF + (2*NESH(I)*NESH(I)*NZ)
       ITMAT = ITMAT + (NORB(I)*NORB(I)*NZ)

      ENDDO

      IF (IPRINT.GE.6) THEN
        IBUF = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN
        NREF = NORB(I)
        CALL HEADER(
     &  'FULLPCTM: Final constructed U1 picture change '//
     &  'transformation matrix in ON RKB/"k" basis, TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(IBUF),NREF,NREF,NREF,NREF,
     &               NZ,IPQTOQ(1,0),LUPRI)
        IBUF = IBUF + (NREF*NREF*NZ)
      ENDIF
      ENDDO
      ENDIF

C =================================================================================
C   ... do save 4c ON-> 2c ON FULL picture change transformation matrix into the file
C =================================================================================
      ISIZE = (NORB(1)*NORB(1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ)
      CALL W2BSSMAT(1,.TRUE.,'U1_ONFUL','FULLPCTM',
     &                TMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'FULLPCTM: U1_ONFUL (FULL picture change transf.matrix '//
     & '4c ONB->2c ONB), TMAT, is written to the BSSMAT'//
     & ' file, EOFLABEL is renewed.'
      ENDIF

C    .. do the test U1+.U1 and U1.U1+ into TBUF...
      IF (IPRINT.GE.2) THEN

        IBUF=1
      DO I=1,NFSYM
      IF (NORB(I).GT.0) THEN
       NREF = NORB(I)
       CALL QGEMM(NREF,NREF,NREF,D1, 
     &    'H','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NREF,NREF,NZ,
     &    'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NREF,NREF,NZ,
     &     D0,IPQTOQ(1,0),TBUF(IBUF),NREF,NREF,NZ)

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'FULLPCTM: U1+ * U1 multiplied '//
     &   'transformation matrix in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QGEMM(NREF,NREF,NREF,D1, 
     &    'N','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NREF,NREF,NZ,
     &    'H','N',IPQTOQ(1,0),
     &     TMAT(IBUF),NREF,NREF,NZ,
     &     D0,IPQTOQ(1,0),TBUF(IBUF),NREF,NREF,NZ)

       CALL TEST_UNIT_MTX(TBUF(IBUF),NREF,NREF,NREF,NZ)

       IF (IPRINT.GE.5) THEN
         CALL HEADER(
     &   'FULLPCTM: U1 * U1+ multiplied '//
     &   'transformation matrix in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL TEST_UNIT_MTX(TBUF(IBUF),NREF,NREF,NREF,NZ)

       IBUF = IBUF + (NREF*NREF*NZ)

       ENDIF
       ENDDO

      ENDIF

C     ... Now do U1+ (TMAT) * H1(BUF) * U1 (TMAT) to get both H2c(ee), H2c(pp) into VMAT
C  ... read H1 into TBUF
C ... read the Dirac Hamiltonian in the ON RKB/"k" basis into  TMAT
      ISIZE = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
      CALL RFBSSMAT('H1_ONBAS',.TRUE.,'FULLPCTM',
     &               TBUF,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'FULLPCTM: H1_ONBAS read into TMAT.'
      ENDIF

      ITMAT = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN
       NREF  = NORB(I)
       N2REF = NORB(I)*NORB(I)
       N2REF = NREF*NREF*NZ

       IF (IPRINT.GE.7) THEN
         CALL HEADER(
     &   'FULLPCTM: TOTAL U1 picture change '//
     &   'transformation matrix in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NREF,NREF,NREF,NREF,
     &                  NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'FULLPCTM: H1 Dirac bare nucleus matrix '//
     &   'in ON RKB/"k" basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(ITMAT),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QTRANS('AOMO','S',D0,NREF,NREF,
     &           NREF,NREF,
     &      TBUF(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      VMAT(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,IPRINT)

       IF (IPRINT.GE.4) THEN
        CALL HEADER(
     & 'FULLPCTM: Obtained U1+ H1 U1 = H2c_inf e/p'//
     & ' in ON RKB/"k" basis, VMAT (output)',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(ITMAT),NREF,NREF,NREF,NREF,
     &           NZ,IPQTOQ(1,0),LUPRI)
       END IF

       ITMAT = ITMAT + N2REF 

      ENDIF
      ENDDO

C =================================================================================
C    .. do save 2c Hamiltonians (ee, pp) into BSSMAT to be later back-transformed
C =================================================================================
      ISIZE = (NORB(1)*NORB(1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ)
      CALL W2BSSMAT(1,.TRUE.,'H2CON_EP','FULLPCTM',
     &               VMAT,ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'FULLPCTM: H2c_infinite (in ON basis, ee and pp) was written'//
     &' into the BSSMAT file, EOFLABEL renewed '
      ENDIF

      CALL QEXIT('FULLPCTM')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck fdtr4c */
      SUBROUTINE FDTR4C(TBUF,TMAT,VMAT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE:   Do the reverse picture change transformation of the (converged)
C ----------   2c Fock matrix into 4c "Fock-Dirac" starting Fock matrix
C             (called from PSISCF/dirscf.F)
C
C  On input : only arrays
C
C  On output: FD4C_RKB written to the BSSMAT to be read later in the DC-SCF
C
C  Written by Miro Ilias, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dcbgen.h" 
#include "dgroup.h" 
#include "dcbdhf.h" 
#include "dcbham.h" 
#include "dcbbas.h" 
#include "dcborb.h" 

      DIMENSION TBUF(*),TMAT(*),VMAT(*),WORK(*)
      LOGICAL FPB, KEEPH1
      CHARACTER*8 LBL

      CALL QENTER('FDTR4C')
#include "memint.h"
        CALL MEMGET('REAL',K1,1,    WORK,KFREE,LFREE)

        CALL MEMCHK('FDTR4C 0.0.1',WORK,1)
C     ... do initiate the important FreeParticleBasis flag!
C   If .false., we are using the one-step RKB approach !
      FPB = .NOT.NOPRTR

C     ... read converged FD2C (from DFFOCK) in LMO into the VMAT !
C ...  read converged Fock-Dirac MO (may be LINSYM!) matrix from the DFFOCK file
       ISIZE = (NESH(1)*NESH(1)*NZ)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NESH(2)*NZ) 
       OPEN (LUFOCK, FILE='DFFOCK',FORM='UNFORMATTED',
     &       ACCESS='DIRECT',RECL=8*ISIZE,STATUS='OLD')
C    &        ACCESS='DIRECT',RECL=8*N2TMOTQ,STATUS='OLD')
C ... RESORT OF RKB_TMAT IS NEEDED FOR THE RKB based transformation !!!
       IREC = -1
       IF(MXDIIS.EQ.0) THEN
          IREC = 1
       ELSE
          IREC = MOD(NITER-2,MXDIIS-1) + 1
       ENDIF
C   ISIZE of 2c = N2TMOTQ
       CALL READAC(LUFOCK,ISIZE,VMAT,IREC)
C      ... we no longer need the DFFOCK file 
C       as the next BSS-SCF creates own new files
       CLOSE(LUFOCK,STATUS="DELETE")
       IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     & 'FDTR4C: Fock-Dirac 2c MO (LMO) matrix '//
     & 'was read from the DFFOCK(deleted afterwards) file into VMAT'
       ENDIF

C      ... read Lowdin 2c MO matrix ... to get 2c AO Fock matrix
C .... read Lowdin matrix into TMAT
        LBL='VMAT2C  '
        CALL RFBSSMAT(LBL,.TRUE.,'FDTR4C',TMAT,N2BBASX,LUBSS,IPRHAM)
        IF (IPRHAM.GE.2) THEN
          WRITE(LUPRI,'(/2X,A,A,A)')                  
     &  'FDTR4C: Lowdin 2c transf. matrix was read'// 
     &  ' into TMAT from the BSSMAT file, label ',LBL           
        ENDIF

C       ... control print out of the Lowdin 2c matrix ....
        IF (IPRHAM.GE.5) THEN
         ITMAT = 1
         DO I=1,NFSYM
         IF (NESH(I).GT.0) THEN
          CALL HEADER(
     &   'FDTR4C: Lowdin transformation matrix 2c_AO -> 2c_MO '//
     &   ' TMAT, read from BSSMAT ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),
     &       NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),LUPRI)
           ITMAT = ITMAT + (NFBAS(I,0)*NESH(I))
         ENDIF 
         ENDDO
        ENDIF 
        CALL MEMCHK('FDTR4C 0.1',WORK,1)

C   ... needs SaomoC matrix (this should work for linear symmetry as well) !
        CALL DZERO(TBUF,N2BBASXQ)
        NZT_SAVE = NZT
        NZT = 1
        CALL MKSAOMO(TBUF,TMAT,.TRUE.,IPRHAM,WORK(KFREE),LFREE)
        NZT = NZT_SAVE
        CALL DCOPY(N2BBASXQ,TBUF,1,TMAT,1)
      
C     =====================================================================================
C      ... do get FD2C in AO basis:  U(TMAT) . FD2CMO(VMAT) . U(TMAT)+  =  FD2CAO (TBUF)
C     =====================================================================================
       IVMAT = 1
       ITMAT = 1
       DO I=1,NFSYM
       IF (NESH(I).GT.0) THEN

        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'FDTR4C: Sovlp*Lowdin half-transformation'//
     &   ' matrix 2c_AO -> 2c_MO ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),
     &       NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &  'FDTR4C: Read FD_2C MO matrix (after BSS-SCF) "'//
     &  ', VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(IVMAT),NESH(I),NESH(I),
     &        NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

         CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &              NESH(I),NESH(I),
     &   TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &   VMAT(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &   TMAT(ITMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &   TMAT(ITMAT),NFBAS(I,0),NESH(I),1,IPQTOQ(1,0),
     &   WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     &  'FDTR4C: FD_2C (after BSS-SCF) matrix in 2c SA-AO basis'//
     &  ' ',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &       NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

        ITMAT = ITMAT + (NFBAS(I,0)*NESH(I))
        IVMAT = IVMAT + (NESH(I)*NESH(I)*NZ)

       ENDIF
       ENDDO
        CALL MEMCHK('FDTR4C 0.2',WORK,1)

C ======================================================================
C      ... do transform 2c FD AO marix into 2c ON RKB/FP "k" basis
C ======================================================================
        ISIZE = (NFBAS(1,0)*NESH(1)*NZ)
        IF (NFSYM.EQ.2) ISIZE = ISIZE + (NFBAS(2,0)*NESH(2)*NZ)
        CALL RFBSSMAT('TM2CAOFP',.TRUE.,'FDTR4C',
     &               TMAT,ISIZE,LUBSS,IPRHAM)
        IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')                  
     &  'FDTR4C: H2c SA_AO "xhi" -> H2c ON RKB/FP "k" '//
     &  ' transf.matrix (TM2CAOFP) read into TMAT.' 
        ENDIF
        CALL MEMCHK('FDTR4C 0.3',WORK,1)

C ==========================================================================
C      ... Do:     TMAT^+ . TBUF . TMAT = VMAT(FD_2c ON)
C ==========================================================================
C  ... do get FD2C in SA-AO basis:
C             U(TMAT) . FD2CMO(VMAT) . U(TMAT)+  =  FD2CAO (TBUF)
C ===========================================================================
       IVMAT = 1
       ITMAT = 1
       IBUF = 1
       DO I=1, NFSYM
       IF (NESH(I).GT.0) THEN

        IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &  'FDTR4C: H2c SA_AO "xhi" -> H2c ON RKB/FP "k" transformation'//
     &  ' matrix',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'FDTR4C: FD2C (after BSS-SCF) matrix in 2c SA-AO basis'//
     &   ' ',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &   '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(I2BASX(I,I)+1),NFBAS(I,0),NFBAS(I,0),
     &       NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                       NESH(I),NESH(I),
     &   TBUF(I2BASX(I,I)+1),NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),
     &   VMAT(IVMAT),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &   TMAT(ITMAT),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &   TMAT(ITMAT),NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &   WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     & 'FDTR4C: Transformed FD_2C matrix (after BSS-SCF) "'//
     & 'in ON RKB/"k" basis, VMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NESH(I),NESH(I),
     &      NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

         ITMAT = ITMAT + (NFBAS(I,0)*NESH(I)*NZ)
         IVMAT = IVMAT + (NESH(I)*NESH(I)*NZ)

       ENDIF
       ENDDO
        CALL MEMCHK('FDTR4C 0.4',WORK,1)

C ====================================================================
C   ... store FD_2C in ON basis into ee part of block diagonal H1
C ====================================================================
C      ... first read H_ee and H_pp into TMAT
      ISIZE = (NORB(1)*NORB(1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ)
      CALL RFBSSMAT('H2CON_EP',.TRUE.,'FDTR4C',TMAT,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')                  
     &  'FDTR4C: H2c (ee,pp) ON RKB/FP "k" '//
     &  ' transf.matrix read into TMAT.' 
      ENDIF

      ITMAT = 1
      IVMAT = 1
      DO I = 1, NFSYM
      IF (NORB(I).GT.0) THEN

        IF (IPRHAM.GE.5) THEN
          CALL HEADER(
     & 'FDTR4C: H_EE+H_PP block diagonal matrix"'//
     & ' in ON RKB/"k" basis (read from BSSMAT), TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(ITMAT),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'FDTR4C: Transformed FD_2C matrix (after BSS-SCF) "'//
     & 'in ON RKB/"k" basis, VMAT, to be inserted into EE',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NESH(I),NESH(I),
     &      NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF
 
C   ... do insert the FD_2C symm.block into 
C       the EE block of the EE+PP block-diagonal matrix...
        DO IZ=1,NZ
        DO J=1,NESH(I)
         IVMAT1 = IVMAT+((IZ-1)*NESH(I)*NESH(I))+((J-1)*NESH(I))
         ITMAT1 = ITMAT+((IZ-1)*NORB(I)*NORB(I)) + ((J-1)*NORB(I))
         CALL DCOPY(NESH(I),VMAT(IVMAT1),1,TMAT(ITMAT1),1)
        ENDDO
        ENDDO

        IF (INI2C.EQ.2) THEN
C       ... do zero the H1_PP block !
         DO IZ=1,NZ
         DO J=1,NESH(I)
!         ITMAT1 = ITMAT+((NORB(I)*NESH(I))+NESH(I)) + 
!    &     ((IZ-1)*NORB(I)*NORB(I)) + ((J-1)*NORB(I))
!         CALL DZERO(TMAT(ITMAT1),NESH(I))
         ENDDO
         ENDDO
!  !  CALL MEMREL('FDTR4C.X2C',WORK,KWORK,KWORK,KFREE,LFREE)
!     return

        IF (IPRHAM.GE.5) THEN
          CALL HEADER(
     & 'FDTR4C: FD_2C (after BSS-SCF)+"0" block diagonal matrix"'//
     & ' in ON RKB/"k" basis, TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(ITMAT),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        ELSE

        IF (IPRHAM.GE.5) THEN
          CALL HEADER(
     & 'FDTR4C: FD_2C (after BSS-SCF)+H_PP block diagonal matrix"'//
     & ' in ON RKB/"k" basis, TMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TMAT(ITMAT),NORB(I),NORB(I),
     &        NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
        ENDIF

        ENDIF

        ITMAT = ITMAT + (NORB(I)*NORB(I)*NZ)
        IVMAT = IVMAT + (NESH(I)*NESH(I)*NZ)

      ENDIF
      ENDDO
        CALL MEMCHK('FDTR4C 0.5',WORK,1)

C   ... do the reverse picture change transformation
C           U . H . U^+ 

C  ... read full picture change transf. matrix into VMAT
      ISIZE = (NORB(1)*NORB(1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NORB(2)*NORB(2)*NZ)
      CALL RFBSSMAT('U1_ONFUL',.TRUE.,'FDTR4C',VMAT,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'FDTR4C: U1_ONFUL (FULL picture change transf.matrix '//
     & '4c ONB->2c ONB), VMAT, was read from to the BSSMAT'
      ENDIF

      ITMAT = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN
       NREF  = NORB(I)
       N2REF = NORB(I)*NORB(I)
       N2REF = NREF*NREF*NZ

       IF (IPRHAM.GE.7) THEN
         CALL HEADER(
     &   'FDTR4C: TOTAL U1 4c->2c picture change '//
     &   'transformation matrix in ON RKB/"k" basis, VMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(VMAT(ITMAT),NREF,NREF,NREF,NREF,
     &               NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'FULLPCTM: FD_2C + H2c_PP '//
     &   'in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NREF,NREF,NREF,NREF,
     &             NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QTRANS('MOAO','S',D0,NREF,NREF,
     &                       NREF,NREF,
     &      TBUF(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      TMAT(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      VMAT(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      VMAT(ITMAT),NREF,NREF,NZ,IPQTOQ(1,0),
     &      WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.4) THEN
        CALL HEADER(
     & 'FULLPCTM: Obtained U1 (FD2C+HPP) U1+ = FD_4C+H1_PP'//
     & ' in ON RKB/"k" basis, TBUF(output)',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(TBUF(ITMAT),NREF,NREF,NREF,NREF,
     &           NZ,IPQTOQ(1,0),LUPRI)
       END IF

       ITMAT = ITMAT + N2REF 

      ENDIF
      ENDDO
        CALL MEMCHK('FDTR4C 0.6',WORK,1)

C  ... do transform FD_4C in ON RKB/"k" into SA-AO (VMAT) (still non_linear symmetry)
      IF (FPB) THEN
        LBL = 'FPMAT4C '
      ELSE
        LBL = 'SL_TM_4C'
      ENDIF 
C  ... read full picture change transf. matrix into VMAT
       ISIZE = (NFBAS(1,0)*NORB(1)*NZ)
       IF (NFSYM.EQ.2) ISIZE = ISIZE + (NFBAS(2,0)*NORB(2)*NZ)
       IF (ISIZE.NE.N2TMT) THEN
         CALL QUIT('FDTR4C>>> ISIZE.NE.N2TMT !')
       ENDIF
       CALL RFBSSMAT(LBL,.TRUE.,'FDTR4C',TMAT,ISIZE,LUBSS,IPRHAM)
       IF (IPRHAM.GE.2) THEN
       WRITE(LUPRI,'(/2X,A,A,A)')
     & 'FDTR4C:',LBL,' 4c "xhi" SA-AO -> 4c ON RKB/"k"'//
     & ' transformation matrix read into TMAT.'
       ENDIF

       CALL DZERO(VMAT,N2BBASXQ)
       CALL MKSAOMO(VMAT,TMAT,.FALSE.,IPRHAM,WORK(KFREE),LFREE)
       CALL DCOPY(N2BBASXQ,VMAT,1,TMAT,1)

      ITMAT = 1
      IBUF  = 1
      IVMAT = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN

       IF (IPRHAM.GE.7) THEN
         CALL HEADER(
     &   'FDTR4C: 4c "xhi" SA-AO -> 4c ON RKB/"k" '//
     &   'transformation matrix, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NORB(I),
     &       NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

         CALL HEADER(
     &   'FDTR4C: 4c (FD4C+H1_PP) '//
     &   'in ON RKB/"k" basis, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),
     &      NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

       CALL QTRANS('MOAO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                       NORB(I),NORB(I),
     &  VMAT(IVMAT),NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),
     &  TBUF(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &  TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &  TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &  WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.5) THEN
        CALL HEADER(
     & 'FDTR4C: 4c (FD4c+H1_4c_pp) in SA-AO matrix "xhi",'//
     & ' VMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NFBAS(I,0),
     &     NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)
       END IF

       ITMAT  = ITMAT + (NFBAS(I,0)*NORB(I)*NZ)
       IVMAT  = IVMAT + (NFBAS(I,0)*NFBAS(I,0)*NZ)
        IBUF  = IBUF + (NORB(I)*NORB(I)*NZ)

      ENDIF
      ENDDO
        CALL MEMCHK('FDTR4C 0.7',WORK,1)

C  ... FINALLY: do transform FD_4C in SA-AO into RKB MO ...

C ... read TMAT from AOMOMAT...
      CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','FDTR4C')
      CALL READT(LUTMAT,N2TMT,TMAT)
      CLOSE(LUTMAT,STATUS='KEEP')
C    
      ITMAT = 1
      IBUF  = 1
      IVMAT = 1
      DO I=1, NFSYM
      IF (NORB(I).GT.0) THEN

       IF (IPRHAM.GE.7) THEN
         CALL HEADER(
     &   'FDTR4C: 4c "xhi" SA-AO -> 4c ON RKB (for DC-SCF)'//
     &   'transformation matrix, TMAT',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NORB(I),
     &       NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     & 'FDTR4C: FD_2c+H2c_pp in 4c SA-AO matrix "xhi",'//
     & ' VMAT',-1)
        WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
        CALL PRQMAT(VMAT(IVMAT),NFBAS(I,0),NFBAS(I,0),
     &     NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &                         NORB(I),NORB(I),
     &  VMAT(IVMAT),NFBAS(I,0),NFBAS(I,0),NZ,IPQTOQ(1,0),
     &  TBUF(IBUF),NORB(I),NORB(I),NZ,IPQTOQ(1,0),
     &  TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &  TMAT(ITMAT),NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),
     &  WORK(KFREE),LFREE,IPRHAM)

       IF (IPRHAM.GE.5) THEN
         CALL HEADER(
     &   'FDTR4C: FD4C+H1PP '//
     &   'in ON_RKB basis, TBUF',-1)
         WRITE(LUPRI,'(3X,A,I1,A,I1)')
     &  '*** Fermion corep ',I,'/',NFSYM
         CALL PRQMAT(TBUF(IBUF),NORB(I),NORB(I),
     &      NORB(I),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       END IF

       ITMAT  = ITMAT + (NFBAS(I,0)*NORB(I)*NZ)
       IVMAT  = IVMAT + (NFBAS(I,0)*NFBAS(I,0)*NZ)
        IBUF  = IBUF + (NORB(I)*NORB(I)*NZ)

      ENDIF
      ENDDO
        CALL MEMCHK('FDTR4C 0.8',WORK,1)

C ... save (FD4C+H1PP) into the BSSMAT to be used as a starting Fock MO matrix
C ..  for DC-SCF ... FD4C_RKB
      ISIZE = NORB(1)*NORB(1)*NZ
      IF (NFSYM.EQ.2) ISIZE = ISIZE + NORB(2)*NORB(2)*NZ
      CALL W2BSSMAT(1,.FALSE.,'FD4C_RKB','FDTR4C',
     &              TBUF,ISIZE,LUBSS,IPRHAM)
      IF (IPRHAM.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     &   'FDTR4C: FD4C_RKB (TBUF) was written into BSSMAT,'//
     &   'EOFLABEL renewed.'
      ENDIF

C     ... this flag must be deactivated for BSS2DC calculations
      START2C = .FALSE.
        CALL MEMCHK('FDTR4C 0.9',WORK,1)

      CALL MEMREL('FDTR4C.X2C',WORK,KWORK,KWORK,KFREE,LFREE)

      CALL QEXIT('FDTR4C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMOTR2C(VMAT,TBUF,CMO,EIG,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE:
C  Generate 2c MO from previous converged 4c MO (read from DFCOEF) by using 
C  the 4c-2c picture change transformation matrix  :
C
C           CMO(2c) = U^{+}(4c->2c) * CMO(4c)
C
C  Resulting 2c MO's are ....
C
C  On entry: VMAT
C            CMO
C            TBUF
C            EIG
C            IBEIG
C
C  Called from:
C
C  Written by Miro ILIAS, Prievidza, 2009
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dcbgen.h" 
#include "dgroup.h" 
#include "dcbdhf.h" 
#include "dcbham.h" 
#include "dcbbas.h" 
#include "dcborb.h" 

      DIMENSION VMAT(*),CMO(*),EIG(*),TBUF(*),WORK(LWORK)
      DIMENSION IBEIG(*)
      LOGICAL FNDLAB,EX,SET_2C
      CHARACTER*8 LBL
#include "memint.h"

      IPRINT = MAX(IPRHAM,IPRSCF)

      IF (IPRINT.GE.5) THEN
        CALL HEADER('Output from CMOTR2C',-1)
      ENDIF

C ================================================================
C         ... get 4c MO from the "DFCOEF" file into CMO
C ================================================================
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CMOTR2C')
      IF(IPRINT.GE.1) THEN
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,15)
      ELSE
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')

C   ... first get the picture change transf.matrix into VMAT
      ISIZE = (NFBAS(1,0)*NFBAS(1,1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NFBAS(2,0)*NFBAS(2,1)*NZ)
      CALL RFBSSMAT('U_PICTRM',.TRUE.,'CMOTR4C',VMAT,
     &               ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'CMOTR2C: Picture change transformation'//
     &  ' matrix was read into VMAT.'
      ENDIF

      !CALL MEMCHK('CMOTR2C 0.1',WORK,1)

      CALL DZERO(TBUF,N2BBASXQ)
C .. do own picture change transformation TBUF=VMAT^{+} * CMO
      IPCTM = 1
      IMO  = 1
      IMOe = 1
      DO I=1,NFSYM
      IF (NESH(I).GT.0) THEN

       IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &  'CMOTR2C: U pict.change transf.matrix 4c AO-> 2c AO,VMAT',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(VMAT(IPCTM),NFBAS(I,0),NFBAS(I,1),              
     &    NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     &  'CMOTR2C: 4c MO coefficients, CMO(LS,e-p)',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(CMO(IMO),NFBAS(I,0),NORB(I),              
     &    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     &  'CMOTR2C: 4c MO coefficients, e-part only, CMO(LS,e)',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(CMO(IMO+NFBAS(I,0)*NESH(I)*NZ),NFBAS(I,0),NESH(I),
     &    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ===============================================================
C     ... do the matrix multiplication...VMAT^{+} * CMO = TBUF
C
C             [U^+](L,LS) * CMO_4c(LS,p-e) = CMO_2c(L,e)
C             [U^+](M,K)  * CMO_4c(K ,  N) = CMO_2c(M,N)
C ===============================================================
!                      M          N        K
        CALL QGEMM(NFBAS(I,1),NESH(I),NFBAS(I,0),D1,
     &  'H','N',IPQTOQ(1,0),VMAT(IPCTM),
     &   NFBAS(I,0),NFBAS(I,1),NZ,
     &  'N','N',IPQTOQ(1,0),
     &   CMO(IMO+(NFBAS(I,0)*NPSH(I)*NZ)),NFBAS(I,0),NORB(I),NZ,
     &   D0,IPQTOQ(1,0),TBUF(IMOe),
     &   NFBAS(I,0),NESH(I),NZ)

       IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &  'CMOTR2C: MO_2c(L(of LS),e) = [U^+](L,LS).MO_4c(LS,e), TBUF',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM  
        CALL PRQMAT(TBUF(IMOe),NFBAS(I,0),NESH(I), 
     &    NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

! zero CMO_2c, NZ=1 !!! (fix later)
!        DO II=1,NESH(I)
!         CALL DZERO(
!    &     TBUF(IMOe+((II-1)*NFBAS(I,0))+NFBAS(I,1)),NFBAS(I,2) )
!        ENDDO

!      IF (IPRINT.GE.5) THEN
!       CALL HEADER(
!    &  'CMOTR2C: MO_2c(LS,e) with zeroed S-fu, TBUF',-1)
!       WRITE(LUPRI,'(A,I1,A,I1)')
!    &  '* Fermion ircop no.',I,'/',NFSYM                          
!       CALL PRQMAT(TBUF(IMOe),NFBAS(I,0),NESH(I), 
!       CALL PRQMAT(TBUF(IMO),NFBAS(I,0),NORB(I), 
!    &    NFBAS(I,1),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
!    &    NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)
!      ENDIF


          IMOe =  IMOe + (NFBAS(I,0)*NESH(I)*NZ)
          IMO  =   IMO + (NFBAS(I,0)*NORB(I)*NZ)
       IPCTM  =  IPCTM + (NFBAS(I,0)*NFBAS(I,1)*NZ) 

      ENDIF
      ENDDO

      !CALL MEMCHK('CMOTR2C 0.11',WORK,1)

C     ... check the orthonormality of generated 2c MO's
      CALL NRMSTMO(TBUF,VMAT,CMO,WORK,LWORK,IPRINT,.FALSE.)


!-------------------------------------------------------------------------------------
!     We have some '2c' MO, construct the 2c Fock matrix to be diagonalized later
!--------------------------------------------------------------------------------------

       !SET_2C=.FALSE.
       SET_2C=.TRUE.

       IF (SET_2C) THEN
       NORB1_SAVE=NORB(1)
       NORB2_SAVE=NORB(2)
       NPSH1_SAVE = NPSH(1)
       NPSH2_SAVE = NPSH(2)

       NORB(1) = NESH(1)
       NORB(2) = NESH(2)
       NPSH(1) = 0
       NPSH(2) = 0

      NSYM = 4/NZ
      DO I=1,NFSYM
       NPSH(I) = 0
       NORB(I) = NESH(I) + NPSH(I)
C  ... new variable ....
       NTMO(I) = NORB(I)
       DO ISYM = 1, NSYM
         NBORB(ISYM,I,2) = 0
         NBORB(ISYM,I,0) = NBORB(ISYM,I,1) + NBORB(ISYM,I,2)
       ENDDO
      ENDDO
C     ... recalculate necessary off-sets
      CALL SETDC2(0)
      CALL SETDC1(0)
      ENDIF



!       we are now in 2c mode...
      !CALL MEMCHK('CMOTR2C 0.2',WORK,1)
!      ... generate 2c density matrix(-es) into CMO (inside)
      CALL DENMAT(CMO,TBUF,IPRINT)

!     IF (SET_2C) THEN
!     NORB(1) = NORB1_SAVE
!     NORB(2) = NORB2_SAVE
!     NPSH(1) = NPSH1_SAVE
!     NPSH(2) = NPSH2_SAVE
!     ENDIF

!     construct Fock matrix
!        CALL MKMOFK(FMO,FAO,DMAT,WORK(KTMAT),WORK(KCMOSAV),
!    &                WORK(KSAOMO),WORK(KDTEMP),WORK(KDEVEC),
!    &                WORK(KFREE),LFREE)




      !CALL MEMCHK('CMOTR2C',WORK,1)
!      call quit('.... stop in CMOTR2C ....')
      WRITE(LUPRI,*) 'Leaving CMOTR2C '    

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
      SUBROUTINE CMOTR4C(TBUF,TMAT,CMO,EIG,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE:   Do the reverse picture change transformation of 2c CMO into 4c CMO
C ----------   to have initial 4c MO's to restart continuing 4c DC-SCF
C             
C            This routine is called when beeing in the two-component mode,
C            because of off-set settings for REACMO
C
C  On input : only arrays
C
C  On output:  the initial 4c MO are written into BSSMAT
C
C  Called from PSISCF (dirscf.F)
C
C  Written by Miro Ilias, Sept 2005, Strasbourg
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dcbgen.h" 
#include "dgroup.h" 
#include "dcbdhf.h" 
#include "dcbham.h" 
#include "dcbbas.h" 
#include "dcborb.h" 

      DIMENSION TBUF(*),TMAT(*),CMO(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL FNDLAB,EX
      CHARACTER*8 LBL
#include "memint.h"

      IPRINT = MAX(IPRHAM,IPRSCF)

      IF (IPRINT.GE.5) THEN
        CALL HEADER('Output from CMOTR4C',-1)
      ENDIF

      IF (INI2C.EQ.3) THEN ! restart code INI2C described in dirrdn.F
C ================================================================
C         ... get 2c MO from the "DFCOEF" file into CMO
C ================================================================
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CMOTR4C')
      IF(IPRINT.GE.1) THEN
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,15)
      ELSE
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,14)
      ENDIF
      CLOSE(LUCOEF,STATUS='KEEP')

CMI .... do the picture change transformation ...
CMI   ... first get the picture change transf.matrix
      ISIZE = (NFBAS(1,0)*NFBAS(1,1)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NFBAS(2,0)*NFBAS(2,1)*NZ)
      CALL RFBSSMAT('U_PICTRM',.TRUE.,'CMOTR4C',TMAT,
     &               ISIZE,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
        WRITE(LUPRI,'(/2X,A)')
     &  'CMOTR4C: Picture change transformation'//
     &  ' matrix was read into TMAT.'
      ENDIF

C   then do TBUF := TMAT . CMO^ 
      IPCTM = 1
      IMO = 1
      DO I=1,NFSYM
      IF (NESH(I).GT.0) THEN

       IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &  'CMOTR4C: U matrix 4c AO-> 2c AO',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(TMAT(IPCTM),NFBAS(I,0),NFBAS(I,1),              
     &    NFBAS(I,0),NFBAS(I,1),NZ,IPQTOQ(1,0),LUPRI)

        CALL HEADER(
     &  'CMOTR4C: 2c MO coefficients, CMO',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(CMO(IMO),NFBAS(I,0),NESH(I),              
     &    NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C     ... do the matrix multiplication...TMAT * CMO = TBUF
C ===============================================================
        CALL QGEMM(NFBAS(I,0),NESH(I),NFBAS(I,1),D1,
     &  'N','N',IPQTOQ(1,0),TMAT(IPCTM),
     &   NFBAS(I,0),NFBAS(I,1),NZ,
     &  'N','N',IPQTOQ(1,0),
     &   CMO(IMO),NFBAS(I,0),NESH(I),NZ,
     &   D0,IPQTOQ(1,0),TBUF(IMO),
     &   NFBAS(I,0),NESH(I),NZ)

       IF (IPRINT.GE.5) THEN
        CALL HEADER(
     &  'CMOTR4C: U.CMO_2C starting MOs(LS,e), TBUF',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(TBUF(IMO),NFBAS(I,0),NESH(I),              
     &    NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

          IMO =   IMO + (NFBAS(I,0)*NESH(I)*NZ)
       IPCTM  = IPCTM + (NFBAS(I,0)*NFBAS(I,1)*NZ) 

      ENDIF
      ENDDO

C     ... check the orthonormality of generated starting MO's
      CALL NRMSTMO(TBUF,TMAT,CMO,WORK,LWORK,IPRINT,.FALSE.)

C ********************************************************
C *
C *       Restarting from unmodified 2c MO
C *
C ********************************************************
      ELSE IF (INI2C.EQ.4) THEN
C    ... read 2c CMO into TBUF ...
       CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CMOTR4C')
       IF(IPRINT.GE.1) THEN
        CALL REACMO(LUCOEF,'DFCOEF',TBUF,EIG,IBEIG,TOTERG,15)
       ELSE
        CALL REACMO(LUCOEF,'DFCOEF',TBUF,EIG,IBEIG,TOTERG,14)
       ENDIF
       CLOSE(LUCOEF,STATUS='KEEP')

      ELSE
        CALL QUIT('CMOTR4C: INI2C differ from 3 or 4!')
      ENDIF

C =========================================
C ... do save initial 4c MO into BSSMAT - TBUF
C =========================================
      ISIZE = (NESH(1)*NFBAS(1,0)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NFBAS(2,0)*NZ)
      CALL W2BSSMAT(1,.FALSE.,'CMO_LS_e','CMOTR4C',
     &               TBUF,ISIZE,LUBSS,IPRINT)

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck cmoini4c */
      SUBROUTINE NRMSTMO(TBUF,TMAT,CMO,WORK,LWORK,IPRINT,Le)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Purpose: CHECK the normalization of pseudo-4c(LS,e)/2c (L,e) CMO.
C
C  On input/output: TBUF - normalized MO (LS,e)/(L,e) pseudo-4comp coefficients
C
C  TMAT,CMO - temporary arrays (content modified)
C
C   Le 
C
C  Called from : CMOTR4C
C                CMOTR2C
C
C  Written by Miro ILIAS, February 2007, Tel Aviv
C                         Dec. 2009, Prievidza
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dcbgen.h" 
#include "dgroup.h" 
C#include "dcbdhf.h" 
C#include "dcbham.h" 
#include "dcbbas.h" 
#include "dcborb.h" 

      DIMENSION TBUF(*),CMO(*),TMAT(*),WORK(LWORK)
      LOGICAL Le

#include "memint.h"

C ... get full S matrix into TMAT
      IF (.NOT.Le) THEN
        CALL GTOVLT(TMAT,D1,IPRINT)
      ELSE
        CALL GTOVLT(TMAT,D0,IPRINT)
      ENDIF

C ...  transform S to the MO's          :   TBUF^+(e,LS) . TMAT(LS,LS) . TBUF(LS,e) = CMO(e,e) = 1
C ...  transform S to the MO's (Le=true):   TBUF^+(e,L) . TMAT(LS,LS) . TBUF(L,e) = CMO(e,e) = 1
      IMAT = 1
      IMAT_L = 1
      IBUF = 1
      IBUF_Le = 1
      IMO = 1
      DO I = 1, NFSYM
      IF (NESH(I).GT.0) THEN

       IF (Le) THEN

        CALL QTRANS('AOMO','S',D0,NFBAS(I,1),NFBAS(I,1),
     &       NESH(I),NESH(I),
     &       TMAT(IMAT),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       CMO(IMO),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &       TBUF(IBUF),
     &       NFBAS(I,1),NESH(I),NZ,IPQTOQ(1,0),
     &       TBUF(IBUF_Le),
     &       NFBAS(I,1),NESH(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRINT)

       ELSE

        CALL QTRANS('AOMO','S',D0,NFBAS(I,0),NFBAS(I,0),
     &       NESH(I),NESH(I),
     &       TMAT(IMAT),NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),
     &       CMO(IMO),NESH(I),NESH(I),NZ,IPQTOQ(1,0),
     &       TBUF(IBUF),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       TBUF(IBUF),
     &       NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),
     &       WORK(KFREE),LFREE,IPRINT)

       ENDIF

       IF (.NOT.Le) THEN
         IF (IPRINT.GE.5) THEN
          CALL HEADER(
     &    'NRMSTMO: Full Sovlp(LL+SS) matrix (from GTOVLT)'//
     &    ' in SA-AO basis "xhi"',-1)
       ELSE
          CALL HEADER(
     &    'NRMSTMO: Full Sovlp(LL,SS=0) matrix (from GTOVLT)'//
     &    ' in SA-AO basis "xhi"',-1)
       ENDIF
          WRITE(LUPRI,'(A,I1,A,I1)')'* Fermion irp no.',I,'/',NFSYM
          CALL PRQMAT(TMAT(IMAT),NFBAS(I,0),NFBAS(I,0),
     &        NFBAS(I,0),NFBAS(I,0),1,IPQTOQ(1,0),LUPRI)

          IF (.NOT.Le) THEN
            CALL HEADER(
     &     'NRMSTMO: Normalized U.CMO_2C starting MOs(LS,e), TBUF',-1)
            WRITE(LUPRI,'(A,I1,A,I1)')
     &      '* Fermion ircop no.',I,'/',NFSYM                          
            CALL PRQMAT(TBUF(IBUF),NFBAS(I,0),NESH(I),              
     &         NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ELSE
            CALL HEADER(
     &     'NRMSTMO: Normalized U+.CMO_4C starting MOs(L,e), TBUF',-1)
            WRITE(LUPRI,'(A,I1,A,I1)')
     &      '* Fermion ircop no.',I,'/',NFSYM                          
            CALL PRQMAT(TBUF(IBUF_Le),NFBAS(I,1),NESH(I),              
     &           NFBAS(I,1),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
          ENDIF

          CALL HEADER(
     &   'NRMSTMO: Resulting Sovlp(ee) in pseudo4c MO -'//
     &   ' must be unit matrix',-1)
          WRITE(LUPRI,'(A,I1,A,I1)')
     &   '* Fermion ircop no.',I,'/',NFSYM                          
          CALL PRQMAT(CMO(IMO),NESH(I),NESH(I),              
     &      NESH(I),NESH(I),NZ,IPQTOQ(1,0),LUPRI)
         ENDIF

CMI    ... check the unit matrix !!!
        WRITE(LUPRI,'(/2X,A)')
     &  'Check of the unit matrix, C+.S.C, C(LS,e) - new 2c4c MOs'
        CALL TEST_UNIT_MTX(CMO(IMO),NESH(I),NESH(I),NESH(I),NZ)

      ENDIF

       IMAT = IMAT + NFBAS(I,0)*NFBAS(I,0)
       IMAT_L = IMAT_L + NFBAS(I,1)*NFBAS(I,1)
       IBUF = IBUF + (NFBAS(I,0)*NESH(I)*NZ)
       IBUF_Le = IBUF_Le + (NFBAS(I,1)*NESH(I)*NZ)
       IMO = IMO + NESH(I)*NESH(I)*NZ

      ENDDO

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck cmoini4c */
      SUBROUTINE CMOINI4C(CMO,TMAT,EIG,IBEIG,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE:   After the reverse picture change transformation of 2c CMO into 4c CMO (in CMOTR4C)
C ----------  fill positronic MOs with zeroes for silly restarting the DC-SCF from BSS-SCF
C
C  On input: only arrays, but needed CMO 2c-4c from BSSMAT,
C             also DFCOEF file is needed
C
C  On output:  the DFCOEF file is renewed !
C
C  Written by Miro Ilias, August, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00,DM1=-1.00D00)
#include "dcbgen.h" 
#include "dgroup.h" 
#include "dcbdhf.h" 
#include "dcbham.h" 
#include "dcbbas.h" 
#include "dcborb.h" 

      DIMENSION CMO(*),TMAT(*),EIG(*),IBEIG(*),WORK(LWORK)
      LOGICAL FPB, EX, KEEPH1
      CHARACTER*8 LBL

      CALL QENTER('CMOINI4C')
#include "memint.h"

C       ... read LS_e initial MO's from BSSMAT file
      ISIZE = (NESH(1)*NFBAS(1,0)*NZ)
      IF (NFSYM.EQ.2) ISIZE = ISIZE + (NESH(2)*NFBAS(2,0)*NZ)

      CALL RFBSSMAT('CMO_LS_e',.FALSE.,'CMOINI4C',
     &        TMAT,ISIZE,LUBSS,IPRHAM)

      IF (IPRHAM.GE.2) THEN
         WRITE(LUPRI,'(/2X,A)')
     & 'CMOINI4C: CMO_LS_e read from BSSMAT'
      ENDIF

C       ... add zeroes to form LS_0e orbitals ...
      IMO = 1
      ITMAT = 1
      DO I = 1, NFSYM
      IF (NESH(I).GT.0) THEN

       IF (IPRSCF.GE.5) THEN
        CALL HEADER(
     & 'CMOINI4C: Molecular orbitals (LS_e) from'//
     & ' 2c->4c picture change transformation',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
        CALL PRQMAT(TMAT(ITMAT),NFBAS(I,0),NESH(I),              
     &  NFBAS(I,0),NESH(I),NZ,IPQTOQ(1,0),LUPRI)       
       ENDIF

        DO IZ = 1, NZ
         IMO2 = IMO + ((IZ-1)*NFBAS(I,0)*NORB(I))
         CALL DZERO(CMO(IMO2),(NPSH(I)*NFBAS(I,0)))
         ITMAT2 = ITMAT + ((IZ-1)*NFBAS(I,0)*NESH(I))
         ISIZE =  NFBAS(I,0)*NESH(I)
         IMO3   = IMO2 + ISIZE
         CALL DCOPY(ISIZE,TMAT(ITMAT2),1,CMO(IMO3),1)
        ENDDO 

        IF (IPRSCF.GE.5) THEN
         CALL HEADER(
     & 'CMOINI4C: Molecular orbitals (LS_0e) from'//
     & ' 2c->4c picture change transformation',-1)
        WRITE(LUPRI,'(A,I1,A,I1)')
     &  '* Fermion ircop no.',I,'/',NFSYM                          
         CALL PRQMAT(CMO(IMO),NFBAS(I,0),NORB(I),              
     &   NFBAS(I,0),NORB(I),NZ,IPQTOQ(1,0),LUPRI)       
        ENDIF

        ITMAT = ITMAT + (NESH(I)*NFBAS(I,0)*NZ)
        IMO = IMO + (NFBAS(I,0)*NORB(I)*NZ)

      ENDIF
      ENDDO

C    ... do renew the DFCOEF file !!!
      INQUIRE(FILE='DFCOEF',EXIST=EX)
      IF (.NOT.EX) THEN
       CALL QUIT(
     & 'CMOINI4C: DFCOEF from the the previous BSS-SCF not found!')
      ENDIF

C     ... do read the DFCOEF file from previous iterations ...
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','CMOINI4C')
      IF(IPRSCF.GE.1) THEN
         CALL REACMO(LUCOEF,'DFCOEF',TMAT,EIG,IBEIG,TOTERG,15)
      ELSE
          CALL REACMO(LUCOEF,'DFCOEF',TMAT,EIG,IBEIG,TOTERG,14)
      ENDIF

C ===================================================
C   Write 2c->4c transformed coefficients ()
C ===================================================
      CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG,.FALSE.)
      CLOSE(LUCOEF,STATUS='KEEP')
      WRITE (LUPRI,'(/2X,A/)')
     & '>>> The DFCOEF file was renewed for restarting 4c SCF'//
     & ' with pseudo-4c MO from previous 2c SCF cycle.'

      CALL QEXIT('CMOINI4C')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck h2ctoll */
      SUBROUTINE H2CTOLL(I2COFK,HMAT,WORK,LWORK,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: Read the H2c_AO integrals, which are in the BSSMAT file 
C             and rewrite them packed into the LL-matrix under a new label
C
C    Called from: MAKE_H2C
C
C    On input:  I2COFK - determines which label from BSSMAT to be read
C               HMAT - empty field (N2BBASXQ) to be utilized 
C
C    On output:  I2COFK set to 7 (see routine GETH2CAO)
C                BSS/DKn integrals in LL_only matrix are written to the BSSMAT
C            
C    Written by Miro Ilias, Strasbourg, 2006
C    Modifications:
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER(THRNULL=1.0D-7,D0=0.00D00)
#include "dcbgen.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "ccom.h"
C
      DIMENSION WORK(LWORK),HMAT(*)
      LOGICAL ISONBSSMAT
      CHARACTER*8 HAM2CL

      CALL QENTER('H2CTOLL')
#include "memint.h"

C.....calculate dimension of (LL block x NZ )
      N2BBASXQ_L = NTBAS(1)*NTBAS(1)*NZ
      CALL MEMGET('REAL',KBSSLL,N2BBASXQ_L,WORK,KFREE,LFREE)
      CALL DZERO(WORK(KBSSLL),N2BBASXQ_L)
C
C     Read original 2c hamiltonian matrix from BSSMAT
C     ===============================================
C
C ....define label for reading hamiltonian integrals from BSSMAT
      IF      (I2COFK.EQ.1) THEN
C ...... h2c from h_D(1)_4C ... DEFAULT ...
         HAM2CL = 'H1AO_DK '
      ELSE IF (I2COFK.EQ.2) THEN
C ...... h2c_CSC(+AMFI if present) from hD(1)_4C
         HAM2CL = 'H1AO_DK2'
      ELSE
         call quit('H2CTOLL: wrong I2COFK ! Should be 1 or 2 ! ')
      ENDIF
      CALL RFBSSMAT(HAM2CL,.TRUE.,'H2CTOLL',HMAT,N2BBASXQ,LUBSS,IPRINT)
      IF (IPRINT.GE.5) THEN
       CALL HEADER(
     & 'H2CTOLL: H2C_AO in full block, total matrix in SA-AO',-1)
       CALL PRQMAT(HMAT,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     Extract the LL block and write to file
C     ======================================
C
      CALL EXTRLL(0,HMAT,WORK(KBSSLL),IPRINT)

      IF(DOSPHE) THEN
C.......do spherical transformation
        CALL SPH_IOTC(WORK(KBSSLL),HMAT,NSPH,IPRINT,WORK(KFREE),LFREE)
        NDIM=NSPH*NSPH*NZ
        CALL W2BSSMAT(1,.TRUE.,'H2CAO_LL','H2CTOLL',HMAT,
     &                NDIM,LUBSS,IPRINT)
      ELSE
        NDIM = N2BBASXQ_L
        CALL W2BSSMAT(1,.TRUE.,'H2CAO_LL','H2CTOLL',WORK(KBSSLL),
     &                NDIM,LUBSS,IPRINT)
      ENDIF
      IF (IPRINT.GE.3) THEN
       IF (DOSPHE) THEN
        write(lupri,*)
     &'H2CTOLL: 2comp. Hamiltonian in LL block written to the BSSMAT'//
     &' in SPHERICAL representation.'
       ELSE
        write(lupri,*)
     &'H2CTOLL: 2comp. Hamiltonian in LL block written to the BSSMAT'//
     &' in CARTESIAN representation.'
       ENDIF
      ENDIF

      IF (IPRINT.GE.5) THEN
       CALL HEADER(
     & 'H2CAO_LL: H2C_AO in LL block, total matrix in AO',-1)
       CALL PRQMAT(WORK(KBSSLL),NTBAS(1),NTBAS(1),NTBAS(1),NTBAS(1),NZ,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF
C
C     In the case of bare nucleus screening correction, 
C     read matrix of screening integrals, extract LL block
C     and write to file
C     ====================================================
C
      IF (BNCRON.AND.ISONBSSMAT('BNCR2C_4','H2CTOLL',LUBSS,IPRINT)) THEN
       CALL RFBSSMAT('BNCR2C_4',.TRUE.,'H2CTOLL',
     &                HMAT,N2BBASXQ,LUBSS,IPRINT)
       CALL EXTRLL(0,HMAT,WORK(KBSSLL),IPRINT)
CMI    ... hoping it will work ...
       IF(DOSPHE) THEN
        CALL SPH_IOTC(WORK(KBSSLL),HMAT,NSPH,IPRINT,WORK(KFREE),LFREE)
        NDIM=NSPH*NSPH*NZ
        CALL W2BSSMAT(1,.TRUE.,'BNCR2C_2','H2CTOLL',HMAT,
     &                NDIM,LUBSS,IPRINT)
       ELSE
        NDIM = N2BBASXQ_L
        CALL W2BSSMAT(1,.TRUE.,'BNCR2C_2','H2CTOLL',WORK(KBSSLL),
     &                NDIM,LUBSS,IPRINT)
       ENDIF
      IF (IPRINT.GE.5) THEN
       CALL HEADER(
     & 'H2CAO_LL: h1bnc in LL block, total matrix in AO',-1)
       CALL PRQMAT(WORK(KBSSLL),NTBAS(1),NTBAS(1),NTBAS(1),NTBAS(1),NZ,
     &            IPQTOQ(1,0),LUPRI)
      ENDIF
      IF (IPRINT.GE.3) THEN
       IF (DOSPHE) THEN
        write(lupri,*)
     &'H2CTOLL: Screened 2comp. Ham in LL block written to the BSSMAT'//
     &' in SPHERICAL representation.'
       ELSE
        write(lupri,*)
     &'H2CTOLL: Screened 2comp. Ham in LL block written to the BSSMAT'//
     &' in CARTESIAN representation.'
       ENDIF
      ENDIF
      ENDIF

C     ... release the memmory..
      CALL MEMREL('H2CTOLL',WORK,KBSSLL,KBSSLL,KFREE,LFREE)

C    ==============================================
C         ... finally set up flag for ONEFCK...
C    ==============================================
      I2COFK = 7
C
      CALL QEXIT('H2CTOLL')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck extrll */
      SUBROUTINE EXTRLL(IREP,PMAT4C,PMAT2C,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  PURPOSE: Extract LL block of 4c PMAT4C into 2c PMAT2C;
C           we are in the SA-AO basis.
C
C  On input:   IREP - irreducible repres. of the operator
C              PMAT4C - four-component matrix
C
C  On output:  PMAT2C - two-comp. symm. blocked matrix
C
C  Written by Miro Ilias, Strasbourg, March 2006
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcbbas.h"
      DIMENSION PMAT4C(*),PMAT2C(*)

      CALL QENTER('EXTRLL')

      IOPSY = JBTOF(IREP,1)

      N2TOT4C = NTBAS(0)*NTBAS(0)
      N2TOT2C = NTBAS(1)*NTBAS(1)

      DO IZ = 1, NZ

         IZOFF4C = (IZ-1)*N2TOT4C + 1
         IZOFF2C = (IZ-1)*N2TOT2C + 1

      DO IFS1 = 1, NFSYM

         IFS2 = MOD(IFS1+IOPSY,2) + 1

         IF (IFS1.EQ.1.AND.IFS2.EQ.1) THEN
           IFS4C = IZOFF4C
           IFS2C = IZOFF2C
         ELSE IF (IFS1.EQ.1.AND.IFS2.EQ.2) THEN 
           IFS4C = IZOFF4C + (NTBAS(0)*NFBAS(1,0))
           IFS2C = IZOFF2C + (NTBAS(1)*NFBAS(1,1))
         ELSE IF (IFS1.EQ.2.AND.IFS2.EQ.1) THEN 
           IFS4C = IZOFF4C + NFBAS(1,0)
           IFS2C = IZOFF2C + NFBAS(1,1)
         ELSE IF (IFS1.EQ.2.AND.IFS2.EQ.2) THEN 
           IFS4C = IZOFF4C + (NTBAS(0)*NFBAS(1,0))+NFBAS(1,0)
           IFS2C = IZOFF2C + (NTBAS(1)*NFBAS(1,1))+NFBAS(1,1)
         ELSE
          write(lupri,*) 'EXTRLL: ifs1,ifs2=',IFS1,IFS2
          call QUIT('EXTRLL: wrong IFS1,IFS2 !')
         ENDIF

         DO I = 1, NFBAS(IFS2,1)
            IPOS4C = IFS4C + (I-1)*NTBAS(0)
            IPOS2C = IFS2C + (I-1)*NTBAS(1)
            CALL DCOPY(NFBAS(IFS1,1),PMAT4C(IPOS4C),1,PMAT2C(IPOS2C),1)
         ENDDO

      ENDDO
      ENDDO

      IF (IPRINT.GE.12) THEN
       CALL HEADER('EXTRLL: Total entering 4c SA-AO oper.',-1)

       CALL PRQMAT(PMAT4C,NTBAS(0),NTBAS(0),NTBAS(0),NTBAS(0),NZ,
     &            IPQTOQ(1,IREP),LUPRI)

       CALL HEADER('EXTRLL: Total picked LL part of SA-AO oper.',-1)
       CALL PRQMAT(PMAT2C,NTBAS(1),NTBAS(1),NTBAS(1),NTBAS(1),NZ,
     &            IPQTOQ(1,IREP),LUPRI)

       IF (IPRINT.GE.16) THEN
         DO IFS1 = 1, NFSYM
            IFS2 = MOD(IFS1+IOPSY,2) + 1
            CALL HEADER(
     & 'EXTRLL: Entering symm. blocked 2c SA-AO operator/4c metric',-1)
            WRITE(LUPRI,'(2X,A,I1,A,I1,A,I1)')
     & 'fermion symm1=',IFS1,' fermion symm2=',IFS2,'/',NFSYM
            CALL PRQMAT(PMAT4C(1+I2BASX(IFS1,IFS2)),
     &            NFBAS(IFS1,0),NFBAS(IFS2,0),
     &        NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,IREP),LUPRI)

         IFS2C = (IFS2-1)*NFBAS(IFS1,1)*NTBAS(1) +
     &            (IFS1-1)*NFBAS(IFS1,1)

            CALL HEADER(
     &  'EXTRLL: Picked symm. blocked 2c SA-AO oper.',-1)
            WRITE(LUPRI,'(2X,A,I1,A,I1,A,I1)')
     & 'fermion symm1=',IFS1,' fermion symm2=',IFS2,' of ',NFSYM
            CALL PRQMAT(PMAT2C(1+IFS2C),
     &        NFBAS(IFS1,1),NFBAS(IFS2,1),
     &        NTBAS(1),NTBAS(1),NZ,IPQTOQ(1,IREP),LUPRI)
         ENDDO
       ENDIF
      ENDIF

      CALL QEXIT('EXTRLL')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck amfiin */
      SUBROUTINE AMFIIN(H2SAO,NB,IPTYPAO,INDXLM,IPNUC,IZTNUCS,
     &                  IPTA,INDXLMA,IAMFU,IANUCS,ISTNUCTYP,CHPNUC,
     &                   NSPH,LABINT,ONLYSO1,ONLYTRA,
     &             IMFCHA,IC_AMFI,IQPOS,CVAL_AMFI,IPRINT,WORK,LWORK)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C 
C  Interface to B.Schimmelpfennig's atomic mean-field integral code, AMFI   
C
C  On input:  H2SAO - atomic integrals in spherical unsorted basis
C             NB - total number of basis functions(L+S)
C              
C          ... IZTNUCS - array for charges (int) of all nuclei of the system
C             
C          ... for each sperical AO in H2CAO matrix: 
C             IPTYPAO - typ of spherical AO
C             INDXLM - index of typ
C             IPNUC - index nuclei (1-NUCDEP)
C             ISTNUCTYP - integer array for storing ICENT, INDXLM 
C 
C          ... for each sperical AO in AMFI triangular matrixes
C             IPTA - typ of spherical AO for AMFI
C             INDXLMA - index of AMFI type
C             IAMFU - mapping function AMFI ints -> H2SAO positions
C
C             LABINT - spherical AO labels
C             ONLYSO1 - is .true. do only comparison of SO1 BSS/AMFI integrals 
C
C          ... for AMFI program (subroutine)
C              IMFCHA
C              IC_AMFI
C              CVAL (from dcbgen.h) - actual speed of light
C
C
C On output: H2SAO with mean-field contributions (only if wished)
C            - only comparison made
C
C     Written by Miro Ilias, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
      PARAMETER(THRNULL=1.0D-7,D0=0.00D00)
#include "maxorb.h"
#include "mxcent.h"
#include "dummy.h"
#include "nuclei.h"
#include "shells.h"
#include "chrnos.h"

C
      LOGICAL BREIT,PROPRI,ISONLIST,ONLYSO1,ONLYTRA,EXAMF,
     &        NOAMFC
      DIMENSION WORK(LWORK),H2SAO(NB,NB,4),IQPOS(1:8)  
      DIMENSION IPTYPAO(NSPH),INDXLM(NSPH),IPNUC(NSPH),
     &          IPTA(*),INDXLMA(*),IAMFU(NSPH),ISTNUCTYP(NSPH,2),
     &          IZTNUCS(NUCDEP),IANUCS(NUCDEP),CHPNUC(NUCDEP)
      DIMENSION ncontrac(0:4)

      CHARACTER*(*) LABINT(*)
      character*8 xa,ya,za
      DIMENSION xa(4),ya(4),za(5)
      CHARACTER*20 AFN
      CHARACTER*3 file_tag
#include "amfi_if.h"
c     logical dirbss_verbose

      CALL QENTER('AMFIIN')
#include "memint.h"

c     dirbss_verbose = .false.
      iprint_save = iprint
!     iprint      = 5

      CALL HEADER('Output from AMFIIN',-1)
      IF (IPRINT.GE.3) THEN
      write(lupri,*) '# individual centers   : ',NUCIND
      write(lupri,*) '# of large-comp. shells: ',NLRGSH
      write(lupri,*) 'charge of atoms        : ',(charge(i),i=1,nucind)
      write(lupri,*) 'Total number of atoms  : ',NUCDEP
C   ... fill IZTNUCS - LIST OF ALL NUCLEIS
      ENDIF

      ix      = 1
      icenter = 0
      do I=1,NONTYP
        ndeg = 0
        do J=1,NONT(I)
          icenter     = icenter + 1
          ndeg        = ndeg + nucdeg(icenter)
          do k = 1, nucdeg(icenter)
            IZTNUCS(IX) = NINT(CHARGE(icenter))
            IX          = IX + 1
          end do
        end do
        if(iprint.ge.3) then
          write(lupri,*) 'atomtype:',i,
     &    ' ==> # nuclei of charge ',IZTNUCS(IX-1), 
     &    ' ==',ndeg
     &     
        endif
      end do

C    .. make ready the IPTYPAO, INDXLM,IPNUC
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'AMFIIN: List of shells of the system:' 
      ENDIF
      IOFF = 1
      JOFF = 1
      IX = 1
      DO ISHELL = 1, NLRGSH
        L     = NHKT(ISHELL)
        ICENT = NCENT(ISHELL)
        NDEG  = NUCDEG(ICENT)
        INDX=NUMCF(ISHELL)
        DO IDEG = 1, NDEG
C        ... find out the nuclei
         INUCPOS = 0
         DO IC=1,ICENT
           IF (IC.LT.ICENT) THEN
             IDEG1=NUCDEG(IC)
           ELSE
             IDEG1=IDEG
           ENDIF
           DO IDG = 1, IDEG1
             INUCPOS = INUCPOS + 1
           ENDDO
         ENDDO

         IF (IPRINT.GE.3) THEN
          write(lupri,'(A,I2,A,I2,A,D10.4,A,I2,A,I2,A,I2,A,I3)')
     &  'Shell L=',NHKT(ISHELL),
     &  ' ICENT=',NCENT(ISHELL),' charge:',CHARGE(ICENT),
     &  ' IDEG/NDEG=',IDEG,'/',NUCDEG(ICENT),' centre=',INUCPOS,
     &  ' INDXLM=',INDX
         ENDIF
        
C      ... specify IPTYPAO,  INDX_L_ML and center
        IF (L.EQ.1) THEN
          IPTYPAO(IX)=1
C          write(lupri,*) 'IX=',IX,' >',WORK(KIPTYPAO+IX)
          INDXLM(IX)=INDX
          IPNUC(IX)=INUCPOS
          IX = IX + 1 
C          write(lupri,*) 'IX set to:',IX
        ELSE IF (L.EQ.2) THEN
          IPTYPAO(IX)=2 ! 2px
C          write(lupri,*) 'IX=',IX,' >',WORK(KIPTYPAO+IX)
          INDXLM(IX)=INDX
          IPNUC(IX)=INUCPOS
          IPTYPAO(IX+1)=3  ! 2py
C          write(lupri,*) 'IX=',IX+1,' >',WORK(KIPTYPAO+IX+1)
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS
          IPTYPAO(IX+2)=4  ! 2pz
C          write(lupri,*) 'IX=',IX+2,' >',WORK(KIPTYPAO+IX+2)
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS
          IX = IX + 3
C          write(lupri,*) 'IX set to:',IX
        ELSE IF (L.EQ.3) THEN
          IPTYPAO(IX  )=5 ! 3d2-
          INDXLM(IX+0)=INDX
          IPNUC(IX)=INUCPOS

          IPTYPAO(IX+1)=6 ! 3d1-
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS

          IPTYPAO(IX+2)=7 ! 3d0
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS

          IPTYPAO(IX+3)=8 ! 3d1+
          INDXLM(IX+3)=INDX
          IPNUC (IX+3)=INUCPOS

          IPTYPAO(IX+4)=9 ! 3d2+
          INDXLM(IX+4)=INDX
          IPNUC (IX+4)=INUCPOS

          IX = IX + 5
        ELSE IF (L.EQ.4) THEN

          IPTYPAO(IX  )=10 ! 4f3-
          INDXLM(IX+0)=INDX
          IPNUC (IX+0)=INUCPOS

          IPTYPAO(IX+1)=11 ! 4f2-
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS

          IPTYPAO(IX+2)=12 ! 4f1-
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS

          IPTYPAO(IX+3)=13 ! 4f0 
          INDXLM(IX+3)=INDX
          IPNUC (IX+3)=INUCPOS

          IPTYPAO(IX+4)=14 ! 4f1+
          INDXLM(IX+4)=INDX
          IPNUC (IX+4)=INUCPOS

          IPTYPAO(IX+5)=15 ! 4f2+
          INDXLM(IX+5)=INDX
          IPNUC (IX+5)=INUCPOS

          IPTYPAO(IX+6)=16 ! 4f3+
          INDXLM(IX+6)=INDX
          IPNUC (IX+6)=INUCPOS

          IX = IX + 7
        ELSE IF (L.EQ.5) THEN

          IPTYPAO(IX+0)=17 ! 5g4-
          INDXLM(IX+0)=INDX
          IPNUC (IX+0)=INUCPOS

          IPTYPAO(IX+1)=18 ! 5g3-
          INDXLM(IX+1)=INDX
          IPNUC (IX+1)=INUCPOS

          IPTYPAO(IX+2)=19 ! 5g2-
          INDXLM(IX+2)=INDX
          IPNUC (IX+2)=INUCPOS

          IPTYPAO(IX+3)=20 ! 5g1-
          INDXLM(IX+3)=INDX
          IPNUC (IX+3)=INUCPOS

          IPTYPAO(IX+4)=21 ! 5g0
          INDXLM(IX+4)=INDX
          IPNUC (IX+4)=INUCPOS

          IPTYPAO(IX+5)=22 ! 5g1+
          INDXLM(IX+5)=INDX
          IPNUC (IX+5)=INUCPOS

          IPTYPAO(IX+6)=23 ! 5g2+
          INDXLM(IX+6)=INDX
          IPNUC (IX+6)=INUCPOS

          IPTYPAO(IX+7)=24 ! 5g3+
          INDXLM(IX+7)=INDX
          IPNUC (IX+7)=INUCPOS

          IPTYPAO(IX+8)=25 ! 5g4+
          INDXLM(IX+8)=INDX
          IPNUC (IX+8)=INUCPOS

          IX            = IX + 9

        else if(l.eq.6.or.l.eq.7.or.l.eq.8)then

          select case(l)
            case(6) ! h
             is_degeneracy = 11
             ilast_IPTYPAO = 25
            case(7) ! i
             is_degeneracy = 13
             ilast_IPTYPAO = 35
            case(8) ! k
             is_degeneracy = 15
             ilast_IPTYPAO = 47
          end select
          do jjj = 1, is_degeneracy
            ilast_IPTYPAO     = ilast_IPTYPAO + 1
            IPTYPAO(IX+jjj-1) = ilast_IPTYPAO
            INDXLM (IX+jjj-1) = INDX
            IPNUC  (IX+jjj-1) = INUCPOS
          end do
          IX = IX + is_degeneracy
          write(lupri,'(A)') 'warning: old amfi routines in use:'//
     &    ' zero AMFI contributions from h, i, k, ... shells '
        else
          write(lupri,'(A)') 'error: old amfi routines in use:'//
     &    ' no AMFI contributions from l, m, ... shells implemented. '
          CALL QUIT('AMFI interface ONLY UP TO k FUNCTIONS!')
        ENDIF
      ENDDO
      ENDDO

      IF ((IX-1).NE.NSPH) THEN
        write(lupri,*) 'AMFIIN error: IX-1=',IX-1, 'NSPH=',NSPH
        CALL QUIT('AMFIIN: IX-1<>NSPH ! CHECK IT !')
      ENDIF

C ... rearange INDXLM for splitted basis sets...
      CALL UNITE_INDXS(NSPH,LABINT,IPNUC,IPTYPAO,
     &                   INDXLM,ISTNUCTYP,IPRINT)

C     Open Mean-field input file generated in READIN      
      LUTMP2 = -1
      LUMNFPRP = -1
      CALL GPOPEN(LUTMP2,'MNF.INP','OLD',' ','FORMATTED',
     &            IDUMMY,.FALSE.)

      BREIT = .FALSE.

C  ++++ Find out individual charges on atoms
      CALL DZERO(CHPNUC,NUCDEP)
      IF (IMFCHA.NE.0) THEN
C     ... find out degeneracy of that nuclei IATOM with charge ICHARGE
        IMINCHA=IZTNUCS(1)
        INDXMINCH=1
        ISUMCH=0
        DO I = 1, NUCDEP
C        .... NOTE: The minimal charge MUST NOT be zero !
         IF ( (IMINCHA.GT.IZTNUCS(I).AND.IZTNUCS(I).NE.0)
     &        .OR.(IMINCHA.EQ.0)) THEN
             IMINCHA=IZTNUCS(I)
             INDXMINCH = I
         ENDIF
         ISUMCH = ISUMCH + IZTNUCS(I)
        ENDDO

        IF (ISUMCH.LE.ABS(IMFCHA)) THEN
         CALL FLSHFO(LUPRI)
         WRITE(LUPRI,'(2X,A,I4,A,I4,A)') 
     &   'Total sum of nucl.charges:',ISUMCH,
     &   ' while given IMFCHA is',IMFCHA,
     &   ' => WARNING, charge too high!'
         CALL FLSHFO(LUPRI)
        ENDIF

C     .... find out the individual charges on atoms
          WRITE(LUPRI,*)
     &     'The total nonzero charge of the system:',IMFCHA
           CHSUM = D0
          DO I=1,NUCDEP
            CHPNUC(I)=DFLOAT(IZTNUCS(I))/DFLOAT(IZTNUCS(INDXMINCH))
            CHSUM=CHSUM+CHPNUC(I)
C           write(lupri,*) I,' -',CHPNUC(I)
          ENDDO
          FACCHA=DFLOAT(IMFCHA)/CHSUM
          write(lupri,*) ' factor is :',FACCHA
          CHSUM = D0
          DO I=1,NUCDEP
            CHPNUC(I)=FACCHA*CHPNUC(I)
            CHSUM = CHSUM + CHPNUC(I)
            write(lupri,*) I,'.atom-nucleus charge:',IZTNUCS(I),
     &      '  partial charge:',CHPNUC(I)
          ENDDO
          write(lupri,*) '     Sum of all charges (real):',CHSUM
          write(lupri,*) 'Total charge of the system is :',IMFCHA
          CALL FLSHFO(LUPRI)
          IF (NINT(CHSUM).NE.IMFCHA) CALL 
     &      QUIT('AMFIIN: error in partial charges !')
      ENDIF

C ===================================
C     We loop over unique centers
C ===================================
      REWIND (LUTMP2)
C     ... skip text
      READ(LUTMP2,*) 
C     ... read number of unique nuclei
      READ(LUTMP2,*) NUNUC
      write(lupri,'(/A,I2)')
     &'  *** number of unique nuclei (from file MNF.INP):',NUNUC

      TOTDIFX=D0
      TOTDIFY=D0
      TOTDIFZ=D0
      ITD=0
C *************************************************************************
C
C                 The cycle for AMFI over unique nuclei
C
C *************************************************************************
      icenter = 0 
      ndeg    = 0
      DO IATOM = 1, NUNUC ! note: NUNUC == NONTYP (atomic types)

       READ(LUTMP2,*) ICHARGE
       WRITE(LUPRI,'(/A,I2,A,I5)')
     &   '  *** calculate AMFI for atom type',
     & IATOM,' with atomic charge ',ICHARGE

C      ... determine total number of nuclei of type IATOM with charge ICHARGE
       CALL IZERO(IANUCS,NUCDEP)
       NICHARGE        = 0
       INDXMINCH       = 100000000 ! that is an unbelievable large # of atoms in a Dirac run -> sknecht feb 2013 fixme: challenge me!
       icenter_tmp     = 0
       DO I = 1, NONT(IATOM)
         icenter       =     icenter     + 1
         icenter_tmp   =     icenter_tmp + 1
         do j = 1, nucdeg(icenter)
!          print *, 'icenter, ndeg, j, iz...',
!    &               icenter, ndeg, j, IZTNUCS(ndeg+j)
           IF(IZTNUCS(ndeg+j) == ICHARGE)THEN
             NICHARGE         = NICHARGE + 1
             INDXMINCH        = min(ndeg+j,INDXMINCH)
             IANUCS(NICHARGE) = ndeg     + j
           ENDIF
         END DO
         ndeg          = ndeg        + nucdeg(icenter)
       END DO

       IF (NICHARGE.EQ.0) THEN
        CALL QUIT('AMFIIN: atom type with given charge not found!')
       ENDIF 
       IF (INDXMINCH > NUCDEP) THEN
        CALL QUIT('AMFIIN: atom type index out of bounds!')
       ENDIF 

 
       IACHNUC=NINT(CHPNUC(INDXMINCH))

       if (IPRINT.GE.2) then
          write(lupri,*) 'The AMFI unique nuclei Z=',ICHARGE,
     &                 ' index=',INDXMINCH,
     &                 ' is degenerate: ',NICHARGE
          write(lupri,*) ' partial charge is (real):',CHPNUC(INDXMINCH)
          write(lupri,*) ' partial charge is (round int):',IACHNUC
          CALL FLSHFO(LUPRI)
       endif

       WRITE(LUPRI,'(A,i4)')
     &   '  *** number of nuclei with identical atom type:',nicharge
       CALL FLSHFO(LUPRI)

C ... now find out if all these (AMFI) nuclei are to be skipped...
       NOAMFC = .FALSE.
       IXC = 0
       DO I = 1, NICHARGE
         DO K = 1, NAMFA
           IF (IAMFA(K).EQ.IANUCS(I)) IXC = IXC + 1
         ENDDO
       ENDDO
       IF (IXC.EQ.NICHARGE) NOAMFC = .TRUE.
       ! no 2e-SO corrections for 1e-atoms like H or hydrogen-like systems!
       if(icharge-iachnuc < 2)then
         write(lupri,*) ' no 2e-SO corrections for hydrogen'//
     &                  ' or hydrogen-like 1e-systems.'//
     &                  ' AMFI is skipped.'
         NOAMFC = .TRUE.
       end if

C     ... file name for AMFI integrals of center
!      FIXME: name tag of file cannot be charge otherwise we do not allow to have different basis sets 
!      on atoms with the same charge (which will be assigned different atom types)
!      sknecht - feb 2013
!      mi/aug2013: oups,right ! fix is below
       nuc_indx=icenter - icenter_tmp + 1
       write(LUPRI,'(3x,a,i2)') 'unique nuclei index:',nuc_indx
       if(nuc_indx < 10)then
         file_tag=CHRNOS(nuc_indx)
       else if(nuc_indx < 100)then
         file_tag=CHRNOS(nuc_indx/10)//CHRNOS(MOD(nuc_indx,10))
       else if(nuc_indx < 1000)then
         file_tag=CHRNOS(nuc_indx/100)//
     &            CHRNOS(MOD(nuc_indx,100)/10)//
     &            CHRNOS(MOD( MOD(nuc_indx,100), 10))
       else
        CALL QUIT('AMFIIN: nucl.indx > 1000! This is not supported!')
       end if

       IF (ICHARGE.LT.10) THEN
        AFN='AOPROPER_MNF.'//CHRNOS(ICHARGE)//'.'//file_tag
       ELSE IF (ICHARGE.LT.100) THEN
        AFN='AOPROPER_MNF.'//CHRNOS(ICHARGE/10)//CHRNOS(MOD(ICHARGE,10))
     &                     //'.'//file_tag
       ELSE IF (ICHARGE.LT.1000) THEN
        AFN='AOPROPER_MNF.'//CHRNOS(ICHARGE/100)//
     &                       CHRNOS(MOD(ICHARGE,100)/10)//
     &                       CHRNOS(MOD( MOD(ICHARGE,100), 10) )
     &                      //'.'//file_tag
       ELSE
        CALL QUIT('AMFIIN: nucl.charge > 1000! This is not supported!')
       ENDIF

       IF (NOAMFC) THEN
         WRITE(LUPRI,'(/3X,A/)')
     &   '    *** This (AMFI) unique nuclei is not to be calculated !'//
     &   ' Only pass (to read input basis) through the AMFI routine.'
       ELSE
         WRITE(LUPRI,'(2A)')
     &   '  *** file with AMFI integrals for this center: ',AFN
       ENDIF

C      ... decide whether to calculate or not the AOPROPER_MNF.XY file
        INQUIRE(FILE=AFN,EXIST=EXAMF)
        IF (.NOT.EXAMF) THEN
C         WRITE(LUPRI,'(4X,A)')
C    &  '...this file does not exist, AMFI is utilized ! '
          CALL GPOPEN(LUMNFPRP,AFN,'UNKNOWN',' ','UNFORMATTED',
     &             IDUMMY,.FALSE.)
        ELSE 
          CALL GPOPEN(LUMNFPRP,AFN,'OLD',' ','UNFORMATTED',
     &             IDUMMY,.FALSE.)
        ENDIF

C      ... find out the correct index of nuclear exponent
       EXPNUCL=D0
       IF (GAUNUC) THEN
        
        IF( (icenter - icenter_tmp + 1) == 0 )THEN
         CALL QUIT('AMFIIN: atom type minimum index == 0')
        ELSE
         if (IPRINT.GE.3) write(lupri,*) ' indx=',
     &   icenter - icenter_tmp + 1,' CHARGE=',
     &   CHARGE(icenter - icenter_tmp + 1),ICHARGE,
     &   ' nucl. exponent:',GNUEXP(icenter - icenter_tmp + 1)
         EXPNUCL = GNUEXP(icenter - icenter_tmp + 1)
        ENDIF

C     ... gauss.nucleus (for SO1) does not work in AMFI !!!!
C     ... Bernd Schim., help us !
        IF (IC_AMFI.EQ.1.OR.IC_AMFI.EQ.4.OR.IC_AMFI.EQ.5) THEN
          CALL
     &  QUIT('AMFIIN: GNUC is not working in AMFI!')
        ENDIF
       ENDIF

C     .... call the AMFI program of B.Schimmelpfennig 
C     ... if file with integrals is present, return back without calculations
       IF (.NOT.EXAMF.AND.NOAMFC) EXAMF=.TRUE. 
         CALL AMFI(LUTMP2,LUMNFPRP,BREIT,GAUNUC,EXPNUCL,
     &             CVAL_AMFI,EXAMF,IC_AMFI,IACHNUC,WORK(KFREE),LFREE)

C      ... read the AOPROPER.MNF file if AMFI has calculated it ..
      IF (.NOT.NOAMFC) THEN

      REWIND (LUMNFPRP)

      READ(LUMNFPRP) xa
C    numballcart... total no of functions per unique center...
      READ(LUMNFPRP) numballcart

C  ... allocate memmory for this AMFI center
      length3=((numballcart*numballcart)+numballcart)/2
      CALL MEMGET('REAL',KSOX,length3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSOY,length3,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KSOZ,length3,WORK,KFREE,LFREE)

!     call dzero(work(KSOX),length3)
!     call dzero(work(KSOY),length3)
!     call dzero(work(KSOZ),length3)

!     CALL MEMCHK('after AMFI 2',WORK,1)

      if (IPRINT.GE.3) WRITE(LUPRI,'(/,2X,A)')
     & 'AMFI order of (centre specific) basis function: '
      DO I = 1, numballcart
        READ(LUMNFPRP) IPTA(I),INDXLMA(I)
        if (IPRINT.GE.3) then
          write(lupri,'(a,i3,a,a4,a,i3,a)')
     &    ' amfi fu no:',I,' type(indx):',LABINT(IPTA(I)),
     &                '(',INDXLMA(I),')'
       endif
      ENDDO

      READ(LUMNFPRP)  Lhigh,(ncontrac(I),I=0,Lhigh) 

      if (IPRINT.GE.3) THEN
        write(lupri,'(1x,a,i2,a,6i4)') 'Lhigh=',Lhigh,' ncontrac:',
     &  (ncontrac(I),I=0,Lhigh)
      endif

      iscontr=0
      DO I=0,Lhigh
       iscontr=iscontr+(ncontrac(I)*(I+I+1))
      ENDDO

      IF (iscontr.NE.numballcart) THEN
        write(lupri,*) 'iscontr,numballcart:',iscontr,numballcart
        CALL QUIT('AMFIIN: iscontr <> numballcart!!')
      ENDIF

C       ... prepare the mapping array AMFI-H2C basis functions...
        CALL IZERO(IAMFU,NSPH)
        DO I=1, NICHARGE
        DO J=1, numballcart
           K=0
 23        K=K+1
            ISONLIST=.FALSE.
            ISONLIST=(IPTA(J).EQ.IPTYPAO(K)).AND.
     &             (INDXLMA(J).EQ.INDXLM(K)).AND.
     &              (IANUCS(I).EQ.IPNUC(K))

           ! write(lupri,*) IPTA(J),IPTYPAO(K)
           ! write(lupri,*) INDXLMA(J),INDXLM(K)
           ! write(lupri,*) IANUCS(I),IPNUC(K)
           ! write(lupri,*) 'ISONLIST=',ISONLIST

          IF (K.LT.NSPH.AND..NOT.ISONLIST) GOTO 23 

          IF (K.EQ.NSPH.AND..NOT.ISONLIST) THEN
            WRITE(LUPRI,*) 'I,J,K:',I,J,K
            WRITE(LUPRI,*) 'NICHARGE,numballcart:',NICHARGE,numballcart
            DO K1=1,NSPH
              WRITE(LUPRI,'(1x,a,i3)') 'K1=',K1
              WRITE(LUPRI,*) '   IPTA,IPTYPAO:',IPTA(K1),IPTYPAO(K1)
              WRITE(LUPRI,*) 'INDXLMA, INDXLM:',INDXLMA(K1),INDXLM(K1)
              WRITE(LUPRI,*) 'IANUCS(I),IPNUC:',IANUCS(I),IPNUC(K1)
            ENDDO
            CALL QUIT('AMFIIN: mapping failed ...')
          ELSE IF (ISONLIST.AND.K.LE.NSPH) THEN
            JI=J+(numballcart*(I-1))
            IAMFU(JI) = K 
          ELSE
           CALL QUIT('AMFIIN: mapping failed 2 ...')
          ENDIF
        ENDDO
        ENDDO

C   ... write out mapping function correspondance....
       IF (IPRINT.GE.2) THEN
        write(lupri,'(/A)')
     &  '***** Mapping spherical atomic basis'//
     &  ' functions "AMFI - H2c" *****'
       ENDIF
       IX = 0
       DO I = 1, NICHARGE
       DO J = 1, numballcart
        IX = IX + 1
        IF (IPRINT.GE.3) THEN
         write(lupri,'(I4,A,I3,A,I4,A,A6,A,I3,A,I4)')
     &   IX,'>  nucl:',IANUCS(I),' charge:',ICHARGE,
     &   ' basis fu: ',LABINT(IPTA(J)),'(',INDXLMA(J),') '//
     &   ' pos (H2cS):',IAMFU(IX)
        ENDIF

        IF ( ( ICHARGE.NE.IZTNUCS(IPNUC(IAMFU(IX))) )    .OR.
     & ( IANUCS(I).NE.IPNUC(IAMFU(IX)) )                 .OR.
     & ( LABINT(IPTA(J)).NE.LABINT(IPTYPAO(IAMFU(IX))) ) .OR.
     & ( INDXLMA(J).NE.INDXLM(IAMFU(IX)) )   )  THEN
         write(LUPRI,*) 'NICHARGE,numballcart=',NICHARGE,numballcart
         CALL QUIT(
     &   'AMFIIN: No compatibilities between H2cS_AO and AMFI !')
        ENDIF
       ENDDO
       ENDDO

C     if (dirbss_verbose) write(lupri,*)
      READ(LUMNFPRP) (WORK(KSOX+irun),irun=0,length3-1)
!     write(lupri,*) ' xpart: ',
!    &               (WORK(KSOX+irun),irun=0,length3-1)
      READ(LUMNFPRP) ya
      READ(LUMNFPRP) (WORK(KSOY+irun),irun=0,length3-1)
      READ(LUMNFPRP) (WORK(KSOZ+irun),irun=0,length3-1)

!     CALL MEMCHK('after AMFI2',WORK,1)
C ************************************************************************
C     ... do comparison of SO1 integrals if wished (ONLYSO1=.true.)
C       or add AMFI contributions into H2CAO
C ************************************************************************

      IF (ONLYSO1.AND.ICHARGE.GT.0) THEN
       WRITE(LUPRI,'(/2X,A,I5/)') 
     & 'AMFIIN:  Doing comparison of SO1 integrals for'//
     & '  nucleus of charge ',ICHARGE
      ENDIF

      IXPOS=IQPOS(4)
      IYPOS=IQPOS(3)
      IZPOS=IQPOS(2)

CMI   if (dirbss_verbose) write(lupri,*)
      IF (IPRINT.GE.3) THEN 
        WRITE(LUPRI,'(2X,A,3I2)')
     &  'AMFIIN: quaternions x,y,z - IQPOS(4-2):',
     &   IQPOS(4),IQPOS(3),IQPOS(2)
      ENDIF

      IF (IXPOS.EQ.IYPOS.OR.IXPOS.EQ.IZPOS.OR.IYPOS.EQ.IZPOS.OR.
     & (IXPOS+IYPOS+IZPOS).NE.9) THEN
        CALL QUIT('AMFIIN: Wrong x,y,z quaternions !')
      ENDIF   

      DO IANUC = 1, NICHARGE

       ISREMOVED = 0
       DO K=1,NAMFA
         IF (IANUCS(IANUC).EQ.IAMFA(K))then
          ISREMOVED = 1 
         end if
       ENDDO

       IF (ISREMOVED.NE.0) THEN
        WRITE(LUPRI,'(/2X,A,I2,A,I3,A/)')
     &  '      *** AMFIIN: nucleus ',IANUCS(IANUC),
     &  ' with charge ',ICHARGE,' NOT ADDED to the BSSn Hamiltonian.'
       ELSE
        WRITE(LUPRI,'(2X,A,I3,A,I3,A)')
     &  '    *** AMFIIN: ADDING nucleus  ',IANUCS(IANUC),
     &  ' with charge ',ICHARGE,' to the BSSn Hamiltonian.'

        IX = 0
        DO  I = 1, numballcart
        DO 111 J = 1, I
C       ... get SO integrals from AMFI  
         SOXINT=-WORK(KSOX+IX)
         SOYINT=-WORK(KSOY+IX)
         SOZINT=-WORK(KSOZ+IX)
!        if(IPTYPAO(IX+1) .gt. 25)then
!          SOXINT = 0.0d0
!          SOYINT = 0.0d0
!          SOZINT = 0.0d0
!        end if
C      ... indexes of H2CAO(I1,J1,2-4)
       I1 = IAMFU(I+(numballcart*(IANUC-1)))
       J1 = IAMFU(J+(numballcart*(IANUC-1)))

       if((DABS(SOXINT) < THRNULL).and.
     &    (DABS(SOYINT) < THRNULL).and.
     &    (DABS(SOZINT) < THRNULL))then
         ix = ix + 1
         goto 111
       end if

       if(iprint >= 5) write(lupri,*) ' hamiltonian indices 1 and 2',
     &                                  i1,j1

C      ... print out AMFI integrals
       IF (IPRINT.GE.3) THEN
        IF (DABS(SOXINT).GT.THRNULL) THEN
          write(lupri,*) 'X amfi unzero contrib: <',
     &   LABINT(IPTA(I)),'(',INDXLMA(I),')',
     &  '[',IANUCS(IANUC),']','|X|',LABINT(IPTA(J)),'(',INDXLMA(J),')',
     &   '[',IANUCS(IANUC),']> ',SOXINT
        ENDIF

        IF (DABS(SOYINT).GT.THRNULL) THEN
          write(lupri,*) 'Y amfi unzero contrib: <',
     &   LABINT(IPTA(I)),'(',INDXLMA(I),')',
     &  '[',IANUCS(IANUC),']','|Y|',LABINT(IPTA(J)),'(',INDXLMA(J),')',
     &   '[',IANUCS(IANUC),']> ',SOYINT
        ENDIF

        IF (DABS(SOZINT).GT.THRNULL) THEN
          write(lupri,*) 'Z amfi unzero contrib: <',
     &   LABINT(IPTA(I)),'(',INDXLMA(I),')',
     &  '[',IANUCS(IANUC),']','|Y|',LABINT(IPTA(J)),'(',INDXLMA(J),')',
     &   '[',IANUCS(IANUC),']> ',SOZINT
        ENDIF

       ENDIF

       IF (IPRINT.GE.5) THEN
       write(lupri,*) 'amfi: <',LABINT(IPTA(I)),'(',INDXLMA(I),')',
     &  '[',IANUCS(IANUC),']','|O|',LABINT(IPTA(J)),'(',INDXLMA(J),')',
     &   '[',IANUCS(IANUC),']> H2CAO: <',LABINT(IPTYPAO(I1)),
     &     '(',INDXLM(I1),
     &    ')','[',IPNUC(I1),']','|O|',LABINT(IPTYPAO(J1)),
     &   '(',INDXLM(J1),
     &   ')[',IPNUC(J1),']>'
       ENDIF

C    ... do only comparison of between AMFI SO1 - BSS SO1 integrals
C   ... only for nucleis with charge > 0 !!!
      IF (ONLYSO1.AND.ICHARGE.GT.0) THEN

       DIFSOX = SOXINT-H2SAO(I1,J1,IXPOS)
       DIFSOY = SOYINT-H2SAO(I1,J1,IYPOS)
       DIFSOZ = SOZINT-H2SAO(I1,J1,IZPOS)
     
       TOTDIFX=TOTDIFX + DABS(DIFSOX)
       TOTDIFY=TOTDIFY + DABS(DIFSOY)
       TOTDIFZ=TOTDIFZ + DABS(DIFSOZ)

       ITD = ITD + 1

       IF (DABS(DIFSOX).GE.THRNULL) THEN
           ratio=SOXINT/H2SAO(I1,J1,4)
           write(lupri,*) '<',LABINT(IPTYPAO(I1)),'(',
     &        INDXLM(I1),')','|SOX|',
     &         LABINT(IPTYPAO(J1)),'(',INDXLM(J1),')','>'
           write(lupri,*) SOXINT,H2SAO(I1,J1,4),DIFSOX,ratio
         ENDIF

         IF (DABS(DIFSOY).GE.THRNULL) THEN
           ratio=SOYINT/H2SAO(I1,J1,3)
           write(lupri,*) '<',LABINT(IPTYPAO(I1)),'(',
     &        INDXLM(I1),')','|SOY|',
     &         LABINT(IPTYPAO(J1)),'(',INDXLM(J1),')','>'
           write(lupri,*) SOYINT,H2SAO(I1,J1,3),DIFSOY,ratio
         ENDIF

         IF (DABS(DIFSOZ).GE.THRNULL) THEN
           ratio=SOZINT/H2SAO(I1,J1,2)
           write(lupri,*) '<',LABINT(IPTYPAO(I1)),'(',
     &        INDXLM(I1),')','|SOZ|',
     &         LABINT(IPTYPAO(J1)),'(',INDXLM(J1),')','>'
           write(lupri,*) SOZINT,H2SAO(I1,J1,2),DIFSOZ,ratio
         ENDIF

         IF (I1.NE.J1) THEN

            DIFSOX = SOXINT+H2SAO(J1,I1,IXPOS)
            DIFSOY = SOYINT+H2SAO(J1,I1,IYPOS)
            DIFSOZ = SOZINT+H2SAO(J1,I1,IZPOS)

            TOTDIFX=TOTDIFX + DABS(DIFSOX)
            TOTDIFY=TOTDIFY + DABS(DIFSOY)
            TOTDIFZ=TOTDIFZ + DABS(DIFSOZ)

            ITD = ITD + 1

          IF (DABS(DIFSOX).GE.THRNULL) THEN
           ratio=SOXINT/H2SAO(J1,I1,4)
           write(lupri,*) '<',LABINT(IPTYPAO(J1)),'(',
     &        INDXLM(J1),')','|SOX|',
     &         LABINT(IPTYPAO(I1)),'(',INDXLM(I1),')','>'
           write(lupri,*) SOXINT,-H2SAO(J1,I1,4),DIFSOX,ratio
          ENDIF

          IF (DABS(DIFSOY).GE.THRNULL) THEN
           ratio=SOYINT/H2SAO(J1,I1,3)
           write(lupri,*) '<',LABINT(IPTYPAO(J1)),'(',
     &        INDXLM(J1),')','|SOY|',
     &         LABINT(IPTYPAO(I1)),'(',INDXLM(I1),')','>'
           write(lupri,*) SOYINT,-H2SAO(J1,I1,3),DIFSOY,ratio
           ENDIF

           IF (DABS(DIFSOZ).GE.THRNULL) THEN
            ratio=SOZINT/H2SAO(J1,I1,2)
           write(lupri,*) '<',LABINT(IPTYPAO(J1)),'(',
     &        INDXLM(J1),')','|SOZ|',
     &         LABINT(IPTYPAO(I1)),'(',INDXLM(I1),')','>'
           write(lupri,*) SOZINT,-H2SAO(J1,I1,2),DIFSOZ,ratio
           ENDIF

         ENDIF
       ELSE IF (.NOT.ONLYTRA) THEN
C  ==================================================
C      ... add AMFI contribution to the H2CAO
C        when onlytra=.false.
C  ==================================================
!       print *, 'H2SAO(I1,J1,IZPOS),  H2SAO(J1,I1,IZPOS), SOZINT',
!    &            H2SAO(I1,J1,IZPOS),  H2SAO(J1,I1,IZPOS), SOZINT
!       print *, 'H2SAO(I1,J1,IYPOS),  H2SAO(J1,I1,IYPOS), SOYINT',
!    &            H2SAO(I1,J1,IYPOS),  H2SAO(J1,I1,IYPOS), SOYINT
!       print *, 'H2SAO(I1,J1,IXPOS),  H2SAO(J1,I1,IXPOS), SOXINT',
!    &            H2SAO(I1,J1,IXPOS),  H2SAO(J1,I1,IXPOS), SOXINT
        H2SAO(I1,J1,IZPOS) = H2SAO(I1,J1,IZPOS) + SOZINT 
        H2SAO(I1,J1,IYPOS) = H2SAO(I1,J1,IYPOS) + SOYINT 
        H2SAO(I1,J1,IXPOS) = H2SAO(I1,J1,IXPOS) + SOXINT 

        H2SAO(J1,I1,IZPOS) = H2SAO(J1,I1,IZPOS) - SOZINT 
        H2SAO(J1,I1,IYPOS) = H2SAO(J1,I1,IYPOS) - SOYINT 
        H2SAO(J1,I1,IXPOS) = H2SAO(J1,I1,IXPOS) - SOXINT 
       ENDIF
       IX = IX + 1
 111  continue
      ENDDO
      ENDIF
      if (IPRINT.GE.3) 
     &   write(lupri,*) '*******************************************'
      ENDDO ! of IANUC....

C     ... when you have integrals, compare them..
C        ... release memmory...
      CALL MEMREL('KSOX',WORK,KWORK,KSOX,KFREE,LFREE)

      ENDIF ! of NOAMFC
      CALL GPCLOSE(LUMNFPRP,'KEEP')

      ENDDO ! of IATOM

C    .. close the MNF.INP file
      CALL GPCLOSE(LUTMP2,'KEEP')

      IF (ONLYSO1) THEN
C       write(lupri,*) ITD,'>>>',TOTDIFX,TOTDIFY,TOTDIFZ
        TOTDIFX=TOTDIFX/DFLOAT(ITD)
        TOTDIFY=TOTDIFY/DFLOAT(ITD)
        TOTDIFZ=TOTDIFZ/DFLOAT(ITD)

        WRITE(LUPRI,'(/3X,A/)')
     &  '**** AMFIIN: average differences between AMFI and '//
     &  'BSS:SO1 spin-orbit integrals (point nucleus) ****'
        WRITE(LUPRI,'(2X,A,D20.10)') '....for SO1X integrals:',TOTDIFX
        WRITE(LUPRI,'(2X,A,D20.10)') '....for SO1Y integrals:',TOTDIFY
        WRITE(LUPRI,'(2X,A,D20.10)') '....for SO1Z integrals:',TOTDIFZ

        CALL QUIT('AMFIIN: Only BSS-AMFI SO1 integral comparisons !')
      ENDIF
      iprint = iprint_save

      CALL QEXIT('AMFIIN')
      RETURN
      END

      INTEGER FUNCTION INDX_LML(IL,IML,I)
CMI  ... called from AMFI/contandmult
#include "implicit.h"
#include "priunit.h"
      DIMENSION IL(*),IML(*)

      IX  =  I
       J  =  0
      LX  =  IL(I)
      MLX =  IML(I)

      IF (IX.GT.0) THEN
 20     CONTINUE
        IF (LX.EQ.IL(IX) .AND. MLX.EQ.IML(IX)) THEN
          IX = IX - 1
           J = J + 1
          IF (IX.GT.0) GOTO 20
        ENDIF
      ENDIF

      INDX_LML = J 

      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck UNITE_INDXS */
      SUBROUTINE UNITE_INDXS(NSPH,LABINT,IPNUC,
     &                   IPTYPAO,INDXLM,ISTNUCTYP,IPRINT)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
C  PURPOSE: If the basis set is SPLITTED on input (splitted blocks of exponents),
C    unite the INDXLM indexes so that their fit for placing correctly 
C    the AMFI integrals.
C
C   This is quite complicated algorithm.
C
C
C   Modified quantity on output: INDXLM
C
C
C  Written by Miro Ilias, Strasbourg, 2005
C             Last modifications: Miro Ilias. Tel Aviv, Jun 2007
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"

      DIMENSION IPNUC(NSPH),IPTYPAO(NSPH),INDXLM(NSPH),ISTNUCTYP(NSPH,2)
      CHARACTER LABINT(MXQN**2)*4
      LOGICAL ITWAS

      CALL QENTER('UNITE_INDXS')

      IF (IPRINT.GE.3) THEN
       write(lupri,'(/2X,A/)')
     & 'UNITE_INDXS: Entering order of basis'//
     & ' functions in spherical AO basis:'
       DO I=1, NSPH
        write(lupri,'(I3,A,I3,A,A4,A,I3,A)')
     &  I,'.fu - center:',IPNUC(I),' type:',
     &  LABINT(IPTYPAO(I)),'(',INDXLM(I),')'
       ENDDO
       write(lupri,*) ' '
      ENDIF

C    ... list of stored INUC&IPTYP
      IT = 0 
      IL = 0
      I = 0
 234  I = I + 1

C     .. skip out before getting out of bonds
      IF (I.GT.NSPH) GOTO 236

      INUC=IPNUC(I)
      IPTYP=IPTYPAO(I)

C     ...find out if INUC & IPTYP are not in the list
      ITWAS=.FALSE.
      DO K=1,IL
       IF (INUC.EQ.ISTNUCTYP(K,1).AND.IPTYP.EQ.ISTNUCTYP(K,2))
     &  ITWAS=.TRUE.
      ENDDO
      IF (ITWAS) GOTO 234

C     ... store INUC&IPTYP into the store_list_array
      IL=IL+1
      ISTNUCTYP(IL,1)=INUC 
      ISTNUCTYP(IL,2)=IPTYP 
C     ... check out indexes from J
      JJ = 0
      J = I-1
 235  J = J + 1

      IF (INUC.EQ.IPNUC(J).AND.IPTYP.EQ.IPTYPAO(J).AND.J.LE.NSPH) THEN
      JJ = JJ + 1
      IT = IT + 1
        IF (IPRINT.GE.4) THEN
          write(lupri,'(A,I3,A,I3,A,I3,A,I4,A,2I4)')
     &   'I=',I,' J=',J,
     &   ' for INUC=',INUC,' IPTYP=',IPTYP,
     &   ' indexes:', JJ,INDXLM(J)
        ENDIF

        IF (INDXLM(J).NE.JJ) THEN
         IF (IPRINT.GE.4) THEN
          write(lupri,'(4X,A,I3,A,I3)')
     &   'difference in indexes --> changing ',
     &   INDXLM(J),' to:',JJ
         ENDIF
C        ... CORE OF THE ROUTINE: change the index 
         INDXLM(J)=JJ
        ENDIF
      ENDIF
      IF (J.LT.NSPH) GOTO 235 
 236  CONTINUE
      IF (I.LT.NSPH) GOTO 234

      IF (IT.NE.NSPH) THEN
        write(lupri,*)
     & 'UNITE_INDXS error output:  IT=',IT,' NSPH=',NSPH
        CALL QUIT('UNITE_INDXS: IT<>NSPH !')
      ENDIF

C ==========================================================================
C     ... when indexes are changed, print out the AMFI function list....
C ==========================================================================
      IF (IPRINT.GE.3) THEN
        write(lupri,'(/2X,A)')
     &  'UNITE_INDXS: Order of (united) AMFI'//
     &  ' spherical atomic basis functions:'
       DO I=1, NSPH
         write(lupri,'(2X,I3,A,I3,A,A4,A,I3,A)')  I,
     &   '.>  center:',IPNUC(I),' atomic funct. type(indx) : ',
     &   LABINT(IPTYPAO(I)),'(',INDXLM(I),')'
       ENDDO
      ENDIF

      CALL QEXIT('UNITE_INDXS')
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck amfi_stuff */
      SUBROUTINE AMFI_STUFF(H2CAO,ONLYSO1,
     &   ONLYTRA,IC_AMFI,WORK,LWORK,IPRINT,IMFCHA)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C    PURPOSE: 
C   ==========
C   Add AMFI integrals to the selected two-component Hamiltonian
C   -  or make comparison of SO1 AMFI integrals with DKH1_SO1 counterparts
C   -  or do only  the CSC transformation of two-comp.Hamiltonian for test purposes
C
C On input: H2CAO - array to be utilized inside....
C                  
C           ONLYSO1 - flag for doing the AMFI-SO1 and BSS-SO1 integrals comparison
C           ONLYTRA - flag for doing  the CSC transformation (testing purposes)
C           IC_AMFI - flag for the AMFI stuff
C           IMFCHA  - total charge of the system - used to adapt the mean-field summation
C
C On output: H2CAO - CSC trasformed (+ AMFI, if wished) 2c SA-AO integrals.
C
C Written by Miro Ilias, Strasbourg, 2005
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
#include "maxaqn.h"
#include "maxorb.h"
#include "mxcent.h"
#include "nuclei.h"
#include "ccom.h"

      DIMENSION H2CAO(*),WORK(LWORK)
      CHARACTER LABINT(MXQN**2)*4
      LOGICAL ONLYSO1,ONLYTRA,TESTCSC

      CALL QENTER('AMFI_STUFF')
#include "memint.h"

C ========================================================
C       ... read H2C_AO into the H2CAO file
C ========================================================
      CALL RFBSSMAT('H1AO_DK ',.TRUE.,'AMFI_STUFF',H2CAO,
     &                  N2BBASXQ,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     & 'AMFI_STUFF: H1AO_DK integrals read into H2CAO.'
      ENDIF

      IF (IPRINT.GE.5) THEN
C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'AMFI_STUFF: READ TOTAL H(2c)_LL in "xhi" SA-AO basis ' ,-1)
        CALL PRQMAT(H2CAO(1),NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

C ... insert phases
      IREP=0
      IF (NZ.LT.4) THEN
      DO IZ = 1,NZ                                                 
         IQ = IPQTOQ(IZ,IREP)                                       
         CALL Q2BPHASE('F',IQ,1,H2CAO(1+N2BBASX*(IZ-1)))
      ENDDO            
      ENDIF           

C   ... transform to BU basis  
      CALL BSTOBU(H2CAO,NZ,WORK(KFREE),LFREE)

      CALL MEMGET('REAL',KBUF1,N2BBASX*4,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KBUF2,N2BBASX*4,WORK,KFREE,LFREE)

C  ... transform to the AO basis ... as in dirone.F/PR1CEX2
      IPAR = JBTOF(IREP,1)
      NB = NTBAS(0)
      CALL DZERO(WORK(KBUF1),4*N2BBASX)
      DO IZ =1,4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTSOAO(H2CAO(1+N2BBASX*(IPQ-1)),
     &             WORK(KBUF1+N2BBASX*(IZ-1)),
     &             NB,IREPD,IPRINT)
      ENDDO

      IF (IPRINT.GE.4) THEN
        CALL HEADER(
     &'AMFI_STUFF: TOTAL H(2c)_LL in "xhi" C AO basis ' ,-1)
C        ... use IQDEF instead of IPQTOQ only when being in AO basis!
        CALL PRQMAT(WORK(KBUF1),NB,NB,NB,NB,4,IQDEF,LUPRI)    
      ENDIF

      CALL SPHGEN(1,2,.FALSE.,WORK(KFREE),LFREE,IPRINT)

C      ... transform C AO, WORK(KBUF1) -> S AO, WORK(KBUF2) 
      CALL DZERO(WORK(KBUF2),4*N2BBASX)
      LCOMP=0 
      DO IZ = 1, 4
         CALL MTAOCS(WORK(KBUF1+N2BBASX*(IZ-1)),NB,
     &        WORK(KBUF2+N2BBASX*(IZ-1)),NB,LCOMP,
     &        WORK(KFREE),LFREE,IPRINT)
      ENDDO

      IF (IPRINT.GE.4) THEN
          CALL HEADER(
     &'AMFI_STUFF: TOTAL ENTERING H(2c)_LL in "xhi" S AO basis ' ,-1)
          CALL PRQMAT(WORK(KBUF2),NB,NB,NB,NB,4,IQDEF,LUPRI)    
        ENDIF

C ... find out total number of L-component spherical functions
      NSPH = NSPHCM(1)

      CALL MEMGET('INTE',KHINLM,NSPH,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPTYPAO,NSPH,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIPNUC,NSPH,WORK,KFREE,LFREE)
!     CALL MEMCHK('XXX',WORK,1)
      CALL SPHLAB(MXQN-1,LABINT)

C     ...allocations for AMFI integrals: typ and index
       CALL MEMGET('INTE',KAIPT,NSPH,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KAILM,NSPH,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KAMFU,NSPH,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KSTNUCTYP,2*NSPH,WORK,KFREE,LFREE)

       CALL MEMGET('INTE',KZTNUCS,NUCDEP,WORK,KFREE,LFREE)
       CALL MEMGET('INTE',KANUCS,NUCDEP,WORK,KFREE,LFREE)
C      .... charges per nucleis
       CALL MEMGET('REAL',KCHPNUC,NUCDEP,WORK,KFREE,LFREE)

!      CALL MEMCHK('before AMFI2',WORK,1)

       CVAL_AMFI = CVAL

       CALL AMFIIN(WORK(KBUF2),NB,WORK(KIPTYPAO),WORK(KHINLM),
     &  WORK(KIPNUC),WORK(KZTNUCS),WORK(KAIPT),WORK(KAILM),WORK(KAMFU),
     &  WORK(KANUCS),WORK(KSTNUCTYP),WORK(KCHPNUC),
     &  NSPH,LABINT,ONLYSO1,ONLYTRA,IMFCHA,IC_AMFI,JM4POS,
     &  CVAL_AMFI,IPRINT,WORK(KFREE),LFREE)
!     CALL MEMCHK('before S->C',WORK,1)

      IF (IPRINT.GE.4) THEN
        CALL HEADER(
     &  'AMFI_STUFF: TOTAL H(2c)_LL (after AMFIIN)'//
     &  ' in "xhi" S AO basis' ,-1)
        CALL PRQMAT(WORK(KBUF2),NB,NB,NB,NB,4,IQDEF,LUPRI)    
      ENDIF
C ===================================================================
C     now transform back: S AO, WORK(KBUF2) -> C AO, WORK(KBUF1) 
C ===================================================================
        CALL DZERO(WORK(KBUF1),N2BBASX*4)
C       Be very careful here as sphgen resets the common blocks for
C       spherical transformations. Need to be set back as soon as 
C       possible (or better revise the code to get rid of these
C       stupid common blocks, ASPG & LV after one day of debugging)
C     ... tranformation with appropriate factors for the backtransform !!
        CALL SPHGEN(1,2,.TRUE.,WORK(KFREE),LFREE,IPRINT)
        DO IZ = 1, 4
          CALL MTAOSC(WORK(KBUF1+N2BBASX*(IZ-1)),NB,
     &          WORK(KBUF2+N2BBASX*(IZ-1)),NB,LCOMP,
     &          WORK(KFREE),LFREE,IPRINT)
        ENDDO
C       reset the common blocks (see above)
        CALL SPHGEN(1,2,.FALSE.,WORK(KFREE),LFREE,IPRINT)

C   ... point to check H2c in C AO basis again...
       IF (IPRINT.GE.4) THEN
          CALL HEADER(
     &'AMFI_STUFF: TOTAL H(2c)_LL BACK in "xhi" C AO basis ' ,-1)
          CALL PRQMAT(WORK(KBUF1),NB,NB,NB,NB,4,IQDEF,LUPRI) 
       ENDIF

C ============================================================
C   ... FINALLY: transform BACK from C AO(BU) to the SA-AO (BS) basis
C            WORK(KBUF1) -> H2CAO(1)
C ============================================================
      CALL DZERO(H2CAO,N2BBASXQ)
      IREP=0
      IPAR = JBTOF(IREP,1)
      NB = NTBAS(0)
      DO IZ = 1, 4
         IREPD = IRQMAT(IZ,IREP)
         IQ    = IQMULT(1,JQBAS(IREPD,IPAR),IZ)
         IPQ   = IQTOPQ(IQ,IREP)
         CALL MTAOSO(WORK(KBUF1+N2BBASX*(IZ-1)),
     &             H2CAO(1+N2BBASX*(IPQ-1)),
     &             NB,IREPD,IPRINT)
      ENDDO

      CALL BUTOBS(H2CAO,NZ,WORK(KFREE),LFREE) 

      IF(NZ.LT.4) THEN                              
        DO IZ = 1,NZ                                
          IQ = IPQTOQ(IZ,IREP)                    
          CALL Q2BPHASE('F',IQ,1,H2CAO(1+N2BBASX*(IZ-1)))  
        ENDDO                                    
      ENDIF

      IF (IPRINT.GE.6) THEN
C    ... print out H2c matrix in AO "xhi" basis
        CALL HEADER(
     &'AMFI_STUFF: AGAIN TOTAL H(2c)_LL in "xhi" SA-AO basis ' ,-1)
        CALL PRQMAT(H2CAO,NTBAS(0),NTBAS(0),
     &              NTBAS(0),NTBAS(0),NZ,IPQTOQ(1,0),LUPRI)
       ENDIF

      CALL MEMREL('AMFI_STUFF',WORK,KWORK,KWORK,KFREE,LFREE)

C  SAVE THE CSC TRANSFORMED INTEGRALS 
      CALL W2BSSMAT(1,.TRUE.,'H1AO_DK2','AMFI_STUFF',
     &              H2CAO,N2BBASXQ,LUBSS,IPRINT)
      IF (IPRINT.GE.2) THEN
       WRITE(LUPRI,'(/2X,A)')
     &'AMFI_STUFF: C-S-C transformed relativistic atomic integrals'//
     &' were written into the BSSMAT file, EOFLABEL renewed.'
      ENDIF

      CALL QEXIT('AMFI_STUFF')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck sph_iotc */
      SUBROUTINE SPH_IOTC(HCAR,HSPH,NSPH,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Transform 2comp matrix from Cartesian to spherical Gaussians
C
C     Written by Trond Saue June 16 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER(D1=1.00D00,D0=0.00D00)
C
#include "ccom.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION HCAR(*),HSPH(*),WORK(LWORK)
#include "memint.h"
      IPRINT = 0
      CALL MEMGET('INTE',KIP  ,NBSYM*4,WORK,KFREE,LFREE)
      CALL SPH_IOT1(HCAR,HSPH,NSPH,WORK(KIP),IPRINT,
     &              WORK(KFREE),LFREE)
C
      CALL MEMREL('SPH_IOTC',WORK,KWORK,KWORK,KFREE,LFREE)      
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck sph_iot1 */
      SUBROUTINE SPH_IOT1(HCAR,HSPH,NSPHTOT,IP,IPRINT,
     &                    WORK,LWORK)
C***********************************************************************
C
C     Transform 2comp matrix from Cartesian to spherical Gaussians
C
C     On output NPSHTOT gives the total number of large components
C     in spherical Gaussian basis
C
C     Written by Trond Saue June 16 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
      PARAMETER(D1=1.00D00,D0=0.00D00)
C
#include "ccom.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION HCAR(*),HSPH(*),IP(4,NBSYM),WORK(LWORK)
#include "memint.h"
C     
C     Generate spherical transformation matrix for each irrep
C
      IC   = 1
      NBRP = 4/NZ
      ICAR = 1
      ISPH = 1
      NSPHTOT = 0
      DO IFRP = 1,NFSYM
        DO ISYM = 1,NBRP
          IBSYM = JFSYM(ISYM,IFRP)
          IBREP = IBSYM - 1
          NCAR  = NBBAS(IBREP,IC)
          NMAT = NBBAS(IBREP,IC)*NBBAS(IBREP,IC)
          CALL MEMGET('REAL',IP(1,IBSYM),NMAT,WORK,KFREE,LFREE)
          CALL SPHCAR(WORK(IP(1,IBSYM)),IBREP,IC,NCAR,NSPH,
     &                .FALSE.,IPRINT,WORK(KFREE),LFREE)
          NSPHTOT = NSPHTOT + NSPH
          IP(2,IBSYM)=NSPH
          IP(3,IBSYM)=ICAR
          IP(4,IBSYM)=ISPH
          ICAR = ICAR + NBBAS(IBREP,IC)
          ISPH = ISPH + NSPH
        ENDDO
      ENDDO
C
C     Transform matrix
C
      CALL SPH_IOT2(HCAR,NTBAS(1),HSPH,NSPHTOT,IP,
     &              IPRINT,WORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Sph_iotc2 */
      SUBROUTINE SPH_IOT2(HCAR,NCAR,HSPH,NSPH,IP,
     &                    IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Pointer IP:
C       IP(1,ISYM) - offset for spherical transformation matrix
C       IP(2,ISYM) - number of spherical components for ISYM
C       IP(3,ISYM) - pointer ICAR
C       IP(4,ISYM) - pointer ISPH
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D1=1.00D00,D0=0.00D00)
C
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
      DIMENSION HCAR(NCAR,NCAR,NZ),HSPH(NSPH,NSPH,NZ),
     &          IP(4,NBSYM),WORK(*)
C
      IC = 1
      NDIM = NSPH*NSPH*NZ
      CALL DZERO(HSPH,NDIM)
      DO 10 JSYM = 1,NBSYM
        JREP = JSYM-1
        NJC  = NBBAS(JREP,1)
        IF(NJC.EQ.0) GOTO 10
        NJS  = IP(2,JSYM)
        JT   = IP(1,JSYM)
        JC   = IP(3,JSYM)
        JS   = IP(4,JSYM)
        DO 20 ISYM = 1,NBSYM
          IREP = ISYM-1
          NIC  = NBBAS(IREP,1)
          IF(NIC.EQ.0) GOTO 20
          NIS  = IP(2,ISYM)
          IT   = IP(1,ISYM)
          IC   = IP(3,ISYM)
          IS   = IP(4,ISYM)
          CALL QTRANS('AOMO','S',D0,NIC,NJC,NIS,NJS,
     &            HCAR(IC,JC,1),NCAR,NCAR,NZ,IPQTOQ(1,0),
     &            HSPH(IS,JS,1),NSPH,NSPH,NZ,IPQTOQ(1,0),
     &            WORK(IT),NIC,NIS,1,IPQTOQ(1,0),
     &            WORK(JT),NJC,NJS,1,IPQTOQ(1,0),
     &            WORK(KFREE),LFREE,IPRINT)
 20     CONTINUE
 10   CONTINUE
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
C  /* Deck amfi_relscf_inp */
      SUBROUTINE AMFI_RELSCF_INP(WORD,RESET)
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C   Input reading for AMFI/RELSCF module.
C
C   Input data are under the *AMFI keyword within the **HAMILTONIAN section.
C
C   Called from HAMINP.
C
C   Written by Miro ILIAS, Prievidza, 2008
C   Last modifications: MI, Odense, June 2008
C                       S. Knecht, Duesseldorf, Jan 2009
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbham.h"
#include "amfi_if.h"
#include "nuclei.h"
      PARAMETER (NTABLE = 20)
      PARAMETER (D1 = 1.0D0,D0=0.0D0,D2=2.0D00,
     &           DTHRS=1.0D-13, DP5 = 0.50D00)
C
      LOGICAL SET, NEWDEF,RESET
      CHARACTER PROMPT*1, WORD*7, TABLE(NTABLE)*7, WORD1*7
C
      SAVE SET
      DATA TABLE /'.PRNT_A','.PRNT_S','.MXITER','.AMFICH',
     &            '.NOAMFC','.ONLSCF','.AMFLVL','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX',
     &            '.XXXXXX','.XXXXXX','.XXXXXX','.XXXXXX'/
      DATA SET/.FALSE./
C
      IF (SET) THEN
         IF (RESET) SET = .FALSE.
         RETURN
      END IF
      IF (RESET) THEN
         SET = .FALSE.
      ELSE
         SET = .TRUE.
      END IF

C     CALL TITLER('AMFI/RELSCF input reading','*',127)

C
C        *** Initialize the default values ***
C
        IPR_AMFI = 0
      IPR_RELSCF = 0

!    ... maximum number of iterations
      MXXTRP_INP     = 50

!    ... default overall mean-field charge
      IMFCH=0
!    ... all atomic centers are included for the AMFI contribution
      ALLAMFA  = .TRUE. ! ALLAMFA = [A]ll [AMF]i [A]toms
      ONRELSCF = .FALSE. ! ONRELSCF = [ON]ly RELSCF, i.e. stop after RELSCF run
      NAMFA    = 0
      call izero(iamfa,nucdep)
!     order of Amfi contributions 
      ISORDER_AMFI_in = -1
C
C     Process input for AMFI/RELSCF
C     ==================================
C
      NEWDEF = (WORD .EQ. '*AMFI  ')

      ICHANG = 0
      IF (NEWDEF) THEN

         CALL TITLER('AMFI/RELSCF input reading','*',127)

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
                     GOTO ( 1, 2, 3, 4,
     &                      5, 6, 7), I
                  END IF
  200          CONTINUE
               IF (WORD .EQ. '.OPTION') THEN
                 CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
                 GO TO 100
               END IF
               WRITE (LUPRI,'(/,3A,/)') ' Keyword "',WORD,
     *            '" not recognized in AMFI_RELSCF_INP.'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal keyword in AMFI_RELSCF_INP.')
 1             CONTINUE
C&&&  .PRNT_A:  print variable for AMFI
               READ(LUCMD,*) IPR_AMFI
               GO TO 100
 2             CONTINUE
C&&&  .PRNT_S:  print variable for RELSCF
                 READ(LUCMD,*) IPR_RELSCF
               GO TO 100
 3             CONTINUE
C&&&  .MXITER:  read maximum number of iterations for RELSCF
                 READ(LUCMD,*) MXXTRP_INP
                 GO TO 100
 4               CONTINUE
C&&&  .AMFICH:    Get the artificial charge of the molecule 
!        for mean-field summations of each atom in the system
!     ...  read integer number
                 READ(LUCMD,*) IMFCH
                 GO TO 100
 5               CONTINUE
C&&& .NOAMFC: specify centers for which AMFI is not to be added to the BSS/X2C Hamiltonian
                 ALLAMFA = .FALSE.
                 READ(LUCMD,*) NAMFA
                 READ(LUCMD,*) (IAMFA(I),I=1,NAMFA)
                 GO TO 100
 6               CONTINUE
C&&& .ONLSCF: stop after RELSCF run - orbital coefficients are saved to
C                file RELSCF_COEF
                 ONRELSCF = .TRUE.
                 GO TO 100
 7               CONTINUE
C&&& .AMFLVL: set up order of AMFI contributions
                 READ(LUCMD,*) ISORDER_AMFI_in
                 if(ISORDER_AMFI_in .lt. 2 .or. 
     &              ISORDER_AMFI_in .gt. 3)then
                   call quit('*** error in AMFI_RELSCF_INP: illegal 
     & value wrt order of AMFI contributions: range: 2 <= x <= 3. ***')
                 end if
                 GO TO 100
            ELSE IF (PROMPT .EQ. '*') THEN
               GO TO 300
            ELSE
               WRITE (LUPRI,'(/,3A,/)') ' Prompt "',WORD,
     &            '" not recognized .'
               CALL PRTAB(NTABLE,TABLE,WORD1//' input keywords',LUPRI)
               CALL QUIT('Illegal prompt in ....')
            END IF
      END IF

  300 CONTINUE
C
C     Print section
C     =============
C
      CALL PRSYMB(LUPRI,'=',75,0)
      WRITE(LUPRI,'(3X,A)')
     & 'Set-up for AMFI/RELSCF calculations'
      CALL PRSYMB(LUPRI,'=',75,0)

      IF (.NOT.NEWDEF) THEN
        WRITE(LUPRI,'(2x,a)')
     &  '...no reading under "*AMFI  ", thus default settings'
      ENDIF

C
      WRITE(LUPRI,'(A,I5)')
     &  ' * AMFI   code print level:',IPR_AMFI,
     &  ' * RELSCF code print level:',IPR_RELSCF,
     &  ' * RELSCF maximum number of iterations:',MXXTRP_INP

      IF (IMFCH.NE.0) THEN
        WRITE(LUPRI,'(1X,A,I4)')
     &  '* AMFI mean-field summations on individual atoms'//
     &  ' are modified due to the artificial charge of the system:',
     &  IMFCH
      ELSE
        WRITE(LUPRI,'(1X,A)')
     &  '* All AMFI mean-field summations are on neutral'//
     &  ' individual atoms.'
      ENDIF
      IF( ONRELSCF ) THEN
        WRITE(LUPRI,'(1X,A)')
     &  '* AMFI will stop after RELSCF run ...'
      END IF
      if(x2c)then
        if(isorder_amfi_in > 1) isorder_amfi_x2c = isorder_amfi_in
        write(lupri,'(1x,a,i3)')
     &  '* order of AMFI contributions to the X2C Hamiltonian:',
     &  isorder_amfi_x2c
        if(isorder_amfi_x2c .eq. 2)then
          write(lupri,'(1x,a)')
     &  ' --> adding spin-same orbit MFSSO2 terms.'
        else if(isorder_amfi_x2c .eq. 3)then
          write(lupri,'(1x,a)')
     &  ' --> adding spin-same + spin-other orbit MFSSO2 terms.'
        end if
      end if

C ... take care of AMFI centers if specified in the input
      IF (.NOT.ALLAMFA) THEN 
          WRITE(LUPRI,'(1X,A,I2,A,$)')
     &' * AMFI contributions NOT calculated for',NAMFA,' centers:'
          WRITE(LUPRI,*) ( IAMFA(I), I = 1, NAMFA )
          WRITE(LUPRI,'(4X,A)')
     &   'List of all nuclei (see "Cartesian Coordinates" above):'
          IX = 0
          DO I=1,NUCIND
          DO J=1,NUCDEG(I)
           IX = IX + 1
           write(lupri,'(5X,I2,A,I4,$)')
     &     IX,'. nuclei, charge:',NINT(CHARGE(I))
           IZX = 0
           DO K = 1, NAMFA
             IF (IAMFA(K).EQ.IX) THEN 
              write(lupri,'(A)')
     &        ' ... no AMFI contribution from this nuclei '
              IZX=1
             ENDIF
           ENDDO
          ENDDO
          ENDDO
C          ... do some check of entering values of NAMFA,IAMFA
          IF (NAMFA.GT.IX) THEN 
           write(lupri,*) ': NAMFA=',NAMFA,
     &     ' is greater than total number of nuclei=',IX
           CALL QUIT('AMFI_RELSCF: Wrong value of NAMFA !!!')
          ENDIF
          DO K=1,NAMFA
           IF (IAMFA(K).GT.IX) THEN
             write(lupri,*)
     &       K,'. AMFI-banned nucleus has number=',IAMFA(K),
     &       ' greater than total number of nuclei=',IX
           CALL QUIT('AMFI_RELSCF: Wrong value of IAMFA(K) !!!')
           ENDIF
          ENDDO
          IF (NAMFA.GE.2) THEN
C           ... check on identical nucleus numbers ....
            DO K=1,NAMFA-1
              DO L=K+1,NAMFA
               IF (IAMFA(K).EQ.IAMFA(L)) THEN
                write(lupri,*) K,'. and ',L,
     &          '. AMFI-banned nuclei have identical numbers ! =',
     &           IAMFA(K),IAMFA(L)
                 CALL QUIT(': Found 2 identical values'//
     &                     ' of IAMFA(K) !!!')
               ENDIF
              ENDDO
            ENDDO
          ENDIF
      ENDIF

C
  999 CONTINUE
C
      RETURN
      END
