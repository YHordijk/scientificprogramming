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
      SUBROUTINE PAMCIM ()
C

      use memory_allocator

      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ====================================================================
C     DIRECT RELATIVISTIC RAS CI PROGRAM
C
C     FIRST VERSION WRITTEN IN JULY 1991 BY LUUK VISSCHER
C     Large revisions in July 1995 by Bert de Jong
C
C     COMMON BLOCKS :
C
C     ABEL    :   MULTIPLICATION TABLE OF ABELIAN SUBGROUP
C     CHR     :   CHARACTER VARIABLES (REPR. NAMES)
C     CMTRIX  :   CI-MATRIX IN THE BASIS OF THE EXPANSION VECTORS
C     CONTROL :   NUMBER OF ITERATIONS, SELECTED VECTORS, ETC.
C     CRITERS :   CONVERGENCE CRITERIA
C     DELSP   :   REORDERING INFORMATION NECESSARY TO DELETE VIRTUALS
C     DGELM   :   DIAGONAL ELEMENTS OF THE LARGE CI-MATRIX
C     EIVALS  :   EIGEN VALUES, RESIDUES, DAVIDSON CORRECTION
C     FNAMES  :   FILE NAMES
C     FNUMBS  :   FILE NUMBERS
C     GMAT    :   G MATRIX ELEMENTS
C     RAS     :   NUMBER OF ORBITALS, ELECTRONS IN THE RAS SUBSPACES
C     RECIND  :   ADDRESSES OF STORED 1-ELECTRON EXCITATIONS
C     SHUF(I) :   WORKSPACE TO REORDER 2-EL. INTEGRALS (*)
C     TIMING  :   TIMING INFORMATION
C     TWOE    :   TWO ELECTRON INTEGRALS
C     VECTOR  :   EXPANSION VECTORS, SIGMA VECTORS (*)
C     VWEIGHT :   VERTEX WEIGHTS
C     IPOINT  :   POINTERS OF THE MAIN MEMORY BLOCK
C
C     BLOCK DATA :
C
C     DIRDAT  :   INITIALIZES FILES NAMES AND NUMBERS
C
C     SUBROUTINES AND FUNCTIONS :
C
C     DIRECT  :   MAIN PROGRAM
C     SETUP   :   GETS NECESSARY DATA TO SET UP CI SPACE
C     DRIVER  :   CONTROLS DIAGONALISATION PROCEDURE
C     ORBSIN  :   READS ORBITAL INFORMATION
C     USERIN  :   READS USER INPUT
C     SHUFFLE :   REORDERING 2-ELECTRON INTEGRALS
C     MAKEGM  :   MAKES GMATRIX
C     VERTEX  :   CALCULATES VERTEX WEIGHTS
C     DETIND  :   CALCULATES ADRESS OF A BIT REPRESENTATION
C     DIAGELM :   CALCULATES DIAGONAL ELEMENTS OF THE LARGE CI-MATRIX
C     TRIAL   :   READS TRIAL CI VECTORS
C     REAVEC  :   READS CI ROOTS FROM PREVIOUS RUN (RESTART)
C     EXPAND  :   CALCULATES SIGMA VECTOR
C     CALIDS  :   CALCULATE ADRESSES OF INTERACTING DETERMINANTS
C     MAKIDS  :   STORE 1-ELECTRON EXCITATIONS
C     GETIDS  :   GET ADRESSES OF INTERACTING DETERMINANTS
C     TWOLIN  :   READS (PART OF) 2-ELECTRON INTEGRALS
C     CINPRD  :   INNER PRODUCT OF TWO VECTORS
C     DIAG    :   DIAGONALISATION OF THE CI MATRIX
C     SELECV  :   SELECT THE APPROPRIATE VECTOR
C     WRIVEC  :   WRITE SELECTED ROOTS TO FILE
C     INITVEC :   INITIALISES TEMPORARY FILE FOR EXPANSION VECTORS
C     PUTVEC  :   GETS EXPANSION VECTOR FROM (TEMPORARY) FILE
C     GETVEC  :   GETS EXPANSION VECTOR FROM (TEMPORARY) FILE
C     RESIDUE :   CALCULATES RESIDUAL VECTOR
C     CRITER  :   CHECKS WHETHER CONVERGENCE OR MAXITER IS REACHED
C     DENSMT  :   CALCULATES DENSITY MATRIX
C     PRECOND :   PRECONDITIONS SELECTED RESIDUAL VECTOR
C     ORTHN   :   SCHMIDT ORTHONORMALIZATION OF THE SELECTED VECTORS
C     USEROUT :   FINAL OUTPUT
C     LEADING :   PRINT DETERMINANTS WITH COEFFICIENT HIGHER THEN COMIN
C     ZZERO   :   INITIALIZES A VECTOR TO ZERO
C     NOVERI  :   CALCULATES N OVER I
C
C     FILES :
C
C      8 (MRCONEE)  : 1-ELECTRON INTEGRALS
C      9 (MDCINT)   : 2-ELECTRON COULOMB INTEGRALS
C     10 (MDBINT)   : 2-ELECTRON BREIT INTEGRALS 
C     11 (MRCTRIV)  : TRIAL VECTORS
C     12 (MRCFINV)  : FINAL VECTORS
C     13 (MRCEXCS)  : COUPLING COEFFICIENTS (TEMPORARY FILE)
C     14 (MRCVECS)  : EXPANSION VECTORS (TEMPORARY FILE)
C     15 (MRCTWOE)  : SORTED TWO-ELECTRON INTEGRALS (TEMPORARY FILE)
C
C     =======================================================================
C
#include "files.inc"
#include "param.inc"
CMKN
#include "general.inc"
CMKN
      LOGICAL DONE
      NAMELIST /CIFOPR/ PROPER,NEOPER,NAMEE
CMKN
      real(8), allocatable :: CI(:)

      call legacy_lwork_get(MCORE)
      call alloc(CI,MCORE,id='CI work space in CIDIRR')
C
      MAXCORE = MCORE
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IERR = 0
      CALL ASNUNIT(5,'-Yon',IERR)
      IF (IERR.NE.0)  WRITE (IW,'(/A,I4/)') ' Assign error ',IERR
#endif
CMKN
C Modified by Malaya .K. Nayak on 4th Feb. 2011
C This part will calculate normal RASCI energy
CMKN
      PROPER = .FALSE.
      OPEN (5,FILE='DIRAC.INP')
      CALL SETUP   (CI,DONE)
      IF (DONE) RETURN
      CALL DRIVER  (CI)
      CALL USEROUT (CI)
C Modified by Malaya .K. Nayak on 4th Feb. 2011
C This part will calculate the RASCI property
      REWIND (5)
      READ (5,CIFOPR,END=700,ERR=700)
  700 CONTINUE
      IF (PROPER) THEN
        WRITE (6,'(//80A1)') ('#',I=1,80)
        WRITE (6,'(A30,//A30)')'Entering RASCI Property Module',
     &                       'Implemented by Malaya K. Nayak'
        DO IEOPER = 1,NEOPER
          NAMEA = NAMEE(IEOPER)
          CALL SETUP   (CI,DONE)
          RESTART = .TRUE.
          IF (DONE) RETURN
          MAXITER = 1
          CALL DRIVER  (CI)
        END DO
        WRITE(6,'(//10x,A41)')
     &    '-- Normal end of RASCI Property Module --'
      END IF
      call dealloc(CI)
C
C     Clean up : Throw away Scratch Files 
C
      CLOSE (MRCEXCS,STATUS='DELETE')
      CLOSE (MRCVECS,STATUS='DELETE')
      CLOSE (MRCTWOE,STATUS='DELETE')
C
      END
C
      SUBROUTINE SETUP (CI,DONE)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ============================================================
C     SETS UP CI SPACE: CALCULATES DIAGRAM, GENERATES DETERMINANTS
C     MAKES G-MATRIX, REORDERS 2-ELECTRON INTEGRALS
C     ============================================================
C
#include "param.inc"
#include "symm.inc"
#include "general.inc"
#include "mempoint.inc"
C
      DIMENSION MINE(16),MAXE(16),IGAS(16),NGASO(256)
      NAMELIST/GOSCIP/NELACT,NGAS,MINE,MAXE,NGASO,REFDET,IREFE,IPRNT,FCI
      REAL*8 CI(*)
      CHARACTER*8 TIMEX,DATEX*10
      LOGICAL DONE
C
      WRITE (6,'(//80A1)') ('#',I=1,80)
      CALL DAYTIME (DATEX,TIMEX)
      WRITE (6,1000) DATEX,TIMEX
      CALL CPUUSED(CPUSET)
C
C     ! In case we run GOSCIP: we're done after that calculation
      NGAS  = 1
      MINE  = 0
      MAXE  = 0
      IGAS  = 0
      NGASO = 0
      NELACT = 0
      REWIND (5)
      READ (5,GOSCIP,END=1,ERR=1)
 1    IF (NELACT.NE.0) THEN
        CALL PAMGOS(CI,MAXCORE,NELACT,NGAS,IGAS,MINE,MAXE,NGASO)
        DONE = .TRUE.
        RETURN
      ELSE
        DONE = .FALSE.
      ENDIF
C
      CALL SYMMIN
      CALL USERIN
      CALL ORBSIN
      CALL SORTMEM
      CALL SHUFFLE (CI(ITR),CI(ITI),CI(IIK),
     &              CI(IIL),CI(ITIR),CI(ITII))
      CALL VERTEX (CI(IVW1))
      CALL CIMEM
      CALL MAKEGM (CI(IGMR),CI(IGMI),CI(ITNR),
     &             CI(ITNI),CI(IIJE))
C
      WRITE (6,1010) REPNA(IREP),NREFDET,NDET,NROOTS,
     &               ISTART,SELECT,(NSEL(I),I=1,NROOTS)
      WRITE (6,1011) MAXITER,CONVERE,CONVERR,CPUMAX,RESTART
      IF (GETDET) THEN
         WRITE (6,1012) GETDET,COMIN,MAKENAT
      ELSE 
         WRITE (6,1013) GETDET,MAKENAT
      ENDIF
      WRITE (6,1020) NELEC,NORB,NORBR
      IF (NORBR(1).EQ.0) MAXH1=0
      IF (NORBR(3).EQ.0) MAXE3=0
      DO 20 J1H = 0, MAXH1
      J1 = NORBR(1) - J1H
      DO 10 J3 = 0, MAXE3
      J2 = NELEC - J1 - J3
      IF (J2.LT.0.OR.J2.GT.NORBR(2)) GO TO 10
      NCDET = BICO(NORBR(1),J1)
      NCDET = BICO(NORBR(2),J2) * NCDET
      NCDET = BICO(NORBR(3),J3) * NCDET

      WRITE (6,1030) J1,J2,J3,NCDET
   10 CONTINUE
   20 CONTINUE
      IF (IGENEX.EQ.1) THEN
         WRITE (6,1040)
      ELSE IF (IGENEX.EQ.2) THEN
         WRITE (6,1041)
      ELSE
         WRITE (6,1042)
         CALL QUIT('cidirr line 174')
      ENDIF
C
      CALL CPUUSED (SEC)
      CPUSET = SEC - CPUSET
      WRITE (6,1050) CPUSET
C
      CALL FLSHFO (6)
C
C     Version 1.0 : Abelian Symmetry
C                   Multiroot Optimization
C                   All CI and Sigma vectors in core memory
C                   Needs trial vectors from GOSCIP
C     Version 1.1 : Generates start vectors by taking lowest diagonal elements
C     Version 1.1.1 : Generates only limited set of 1-electron excitations
C                     in 2-electron part.
C     Version 1.1.2 : Option to write 1-electron excitation to file.
C     Version 1.1.3 : Keeps 1-el. excits in core if possible.
C                     Writes CI-Vectors to file to be able to restart.
C     Version 1.2.0 : Calculates Density Matrix
C                     Code adapted to Hewlett Packard Workstations
C     Version 1.2.1 : Special code for integrals with 3 or more external labels
C     Version 1.2.2 : Dynamic Memory Allocation checking.
C     Version 1.3   : Selecting leading determinants.  
C     Version 1.3.1 : Generates MFDVECA file with natural spinors
C     Version 1.3.2 : BugFix. Now generates totaly symmetric natural spinors 
C     Version 1.3.3 : Option to delete virtual spinors 
C     Version 1.4   : Reduction of memory required : write out CI vectors
C                     Deleted CNVE statements
C     Version 1.5   : Ported to IBM RS6000
C     Version 1.6   : Reads new integral format (reduced list)
C                     MRCTWOE is now a real scratch file
C                     Removed check on short vectors : Non-stride variant
C                     in Expand will always be used.
C     Version 1.6.1 : Real arithmetic versions of EXPAND and TWOLIN
C     Version 1.6.2 : Corrected restart bug introduced in version 1.4
C                     Improved calculation of Davidson correction
C                     Changed input defaults for LEADDET nad NATURAL
C     Version 1.6.3 : Corrected bug in EXPAND (1-electron contribution)
C                     Reinstalled ORBDEL option
C     Version 1.6.4 : Can now run with zero iterations to calculate natural
C                     spinors after a crash
C     Version 1.6.5 : Installed explicit real diagonalization routine
C     Version 2.0   : Large revisions ; Better memory management
C                     Storing as much coupling coeff. as possible
C                     More efficient expansion algorithm
C                     Faster calculation of coupling coeff.
C
 1000 FORMAT (' DIRect Relativistic CI program version 2.0'//
     &' Written by Luuk Visscher ',
     &//' This run started at ', A10,1X,A8)
 1010 FORMAT (/' Abelian Representation :',1X,A4/
     &' Number of reference determinants :',I10/
     &' Total number of determinants :',I10/
     &' Number of roots to be optimized :',I3/
     &' Start vectors method :',I3/
     &' Overlap Selection : ',L1/
     &' Position in COSCI space :',(10I4))
 1011 FORMAT (' Maximum number of iterations :',I5/
     &' Energy convergence required :',1X,G10.1/
     &' Residue convergence required :',G10.1/
     &' Maximum of CPU seconds to be used :',F16.1/
     &' Restart calculation : ',L1)
 1012 FORMAT (' Printing contributions of important spinors : ',L1/
     &' Contributions are printed for coefficients larger then :',G10.1/
     &' Generate natural spinors after CI calculation : ',L1)
 1013 FORMAT (' Printing contributions of important spinors : ',L1/
     &' Generate natural spinors after CI calculation : ',L1)
 1020 FORMAT (/' Total number of active electrons :',I5/
     &' Number of active orbitals :',I5//
     &' RAS1 :',I4/' RAS2 :',I4/' RAS3 :',I4//
     &' Allowed RAS configurations :'//
     &'  RAS1  RAS2  RAS3  Determinants (All Representations)')
 1030 FORMAT (3(2X,I4),I12)
 1040 FORMAT (/' 1 and 2-electron excitations will be evaluated ',
     &'directly')
 1041 FORMAT (/' 1-electron excitations will be written to file')
 1042 FORMAT (/' No valid value for IGENEX : choose 1 (direct) or',
     &' 2 (file)')
 1050 FORMAT (/' Setup Routine completed in ',F8.3,' SECONDS')
C
      RETURN
      END
C
      SUBROUTINE DRIVER (CI)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ======================================
C     DRIVES DAVIDSON DIAGONALISATION PROCES
C     ======================================
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
#include "mempoint.inc"
CMKN
#include "iterat.inc"
CMKN
C
      real*8  CI(*)
      LOGICAL ENDRUN
c
      MVEC = 0
      NVEC = 0
C
C     Fill Core, and if IGENEX=2 also MRCEXCS
C     with one-electron coupling coefficients
C
      CALL MAKIDS (CI(IPIL),CI(IPIR),CI(IIJR),
     &             CI(IPJL),CI(IPJR),CI(IPJV),
     &             CI(IVW1))
C
C     Construct diagonal elements CI matrix
C
      CALL DIAGELM (CI(IPDR),CI(IPDG),CI(ITNR),
     &              CI(ITNI),CI(IIJE),CI(IGMR),
     &              CI(IGMI),CI(IPIL),CI(IPIR),
     &              CI(IPJL),CI(IPJR),CI(IPJV),
     &              CI(IVW1),CI(IIJR))
C
      CALL INITVEC
C
      IF (RESTART) THEN
         CALL REAVEC (CI(IPVR),CI(IPVI),CI(IPDR),
     &                CI(IPDI),CI(IPIL))
      ELSE
         LCOREX = NCOREX / INTOWP(1) + MOD(NCOREX,INTOWP(1))
         CALL TRIAL (CI(IPVR),CI(IPVI),CI(IPIL),
     &               CI(IPIR),CI(IPJL+LCOREX),
     &               CI(IPJR+LCOREX),CI(IVW1),
     &               CI(IPDG))
      ENDIF
      IF (.NOT.PROPER) WRITE (6,1000)
C
      DO 200 ITER = 1, MAXITER
         CALL CPUUSED(CPITER)
         CALL EXPAND (CI(IPDR),CI(IPDI),CI(IPVR),
     &                CI(IPVI),CI(IPSR),CI(IPSI),
     &                CI(IPIL),CI(IPIR),CI(IPJL),
     &                CI(IPJR),CI(IPJV),CI(ITNR),
     &                CI(ITNI),CI(IIJE),CI(IGMR),
     &                CI(IGMI),CI(IVW1),
     &                CI(IVW1+IVWW),CI(IIJR))
CMKN
C Commenting this if statement will give all-states property
CMKN
      IF (PROPER) RETURN
CMKN
         CALL CINPRD (CI(IPVR),CI(IPVI),CI(IPSR),
     &                CI(IPSI))
CMKN
C This part is meant for all-states property
CMKN
      IF (PROPER) THEN
      WRITE(*,25)NAMEA
      DO I = 1,NVEC
      DO J = 1,NVEC
      WRITE(*,50)I,J,CMTR(I,J),CMTI(I,J)
      ENDDO
      ENDDO
      GO TO 300
      ENDIF
CMKN
C Rest of this subroutine will be called only for normal CI
CMKN
         CALL DIAG
         CALL SELECV
         CALL WRIVEC (CI(IPVR),CI(IPVI),CI(IPDR),
     &                CI(IPDI),CI(IPIL),CI(IPIR))
         CALL RESIDUE (CI(IPVR),CI(IPVI),CI(IPSR),
     &                 CI(IPSI))
         CALL CRITER  (ENDRUN)
         IF (ENDRUN) GOTO 300
         CALL PRECOND (CI(IPDG),CI(IPVR),CI(IPVI))
         CALL ORTHN (CI(IPVR),CI(IPVI))
  200 CONTINUE
  300 CONTINUE
      RETURN
C
CMKN
   25 FORMAT (/' ALL-STATES PROPERTY: ',A8,17X,'REAL',14X,'IMAGINARY')
   50 FORMAT (/' Property between states',2I4,' is:',2F20.12)
CMKN
 1000 FORMAT (//1X,'Iter',1X,'Root',1X,'Isel',13X,'Energy',
     & 7X,'Correlation',5X,'Convergence',3X,'Residue')
C
      END
C
      SUBROUTINE CALIDS_UP (NORB1,NORB2,NORB3,NORB4,ISYM,JDET,NVERT,
     &                      ILIND,IRIND,JLIND,JRIND,JVECL,IVWU,
     &                      IVDET,JVDET)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ======================================================================
C     GIVEN : - EXCITATION A+(NORB1)A(NORB2) = E(NORB1,NORB2)
C             - ABELIAN SYMMETRY CHARACTER OF BRA OR KET DETERMINANTS
C               ISYM>0 KET DETERMINANTS BELONG TO REPRESENTATION ISYM
C               ISYM<0 BRA DETERMINANTS BELONG TO REPRESENTATION ISYM
C     GENERATE : - ADRESSES OF ALL INTERACTING DETERMINANTS
C                  JLIND: TAIL PART OF ADRESSES OF BRA DETERMINANTS
C                  JRIND: TAIL PART OF ADRESSES OF KET DETERMINANTS
C                         PLUS SIGN OF INTERACTION MULTIPLIED BY HEAD PART
C     DIRECTION : WALK FROM BOTTOM TO TOP OF GRAPH    
C     ======================================================================
C
#include "param.inc"
#include "symm.inc"
C
      INTEGER JVECL(2,*)
      INTEGER ILIND(*),IRIND(*),JLIND(*),JRIND(*)
      INTEGER IVDET(2*MXREP,0:NORB,0:NELEC),
     &        JVDET(2*MXREP,0:NORB,0:NELEC)
      INTEGER IVWU(2*MXREP,0:NORB,0:NELEC)
      LOGICAL IJARR
C
      JDET = 0
      NVERT = 1
      JVECL(1,1) = 0
      JVECL(2,1) = 0
      JVECL(1,2) = 0
C     -------------------------------------------------------------
C     SET CONSTRAINTS IMPOSED BY THE PREVIOUS (IF ANY) EXCITATION :
C     A+(NORB3)A(NORB4)
C     -------------------------------------------------------------
      ICREA = NORB3 - 1
      IANNI = NORB4 - 1
      IF (IANNI.EQ.ICREA) IANNI = - 1
      NDORB = NORB1 - NORB2
      IF (NDORB.LT.0) THEN
         ITAIL = NORB2
         IHEAD = NORB1
      ELSE
         ITAIL = NORB1
         IHEAD = NORB2
      ENDIF
C
      ITSYM = IRPAMO(ITAIL)
      IHSYM = IRPAMO(IHEAD)
      IF (ISYM*NDORB.GT.0) THEN
         LSYM = MULTBI(IHSYM,ABS(ISYM))
         LSYM = MULTB(ITSYM,LSYM)
      ELSE
         LSYM = ABS(ISYM)
      ENDIF
C     -----------------------------
C     LAST SWITCH MUST END IN JXIND
C     -----------------------------
      IF (MOD(NORB-IHEAD,2).NE.0) THEN
         IJARR = .TRUE.
         JLIND(1) = 1
         JRIND(1) = 1
      ELSE
         IJARR = .FALSE.
         ILIND(1) = 1
         IRIND(1) = 1
      ENDIF
C
      DO 30 I = 0, NORB
        DO 20 J = MAX0(IORBR(I)-1,0), MIN0(NELEC,I+1)
          DO 10 IR = 1, NREP
            IVDET(IR,I,J) = 0
            JVDET(IR,I,J) = 0
   10     CONTINUE
   20   CONTINUE
   30 CONTINUE
C
      IVDET(LSYM,NORB,NELEC) = 0
      JVDET(LSYM,NORB,NELEC) = 1
C     --------------
C     MAKE ALL TAILS
C     --------------
      DO 300 I = NORB-1, ITAIL, -1
      IJARR = (.NOT.IJARR)
      JDET = 0
      JMIN = IORBR(I)
      JMIN = MAX0(NELEC+I-NORB,0)
      JMAX = MIN0(NELEC,I)
C     ----------------------------------------------------------------
C     LOOP OVER VERTICES : CHECK IF BOTH (PARALLEL) PATHS ARE POSSIBLE
C     ----------------------------------------------------------------
      DO 200 J = JMAX, JMIN, -1
      JLMIN = (NREP / 2 * MOD(J+1,2)) + 1
      JLMAX = JLMIN + (NREP / 2) - 1
      DO 100 JLR = JLMIN, JLMAX
      JRR = MULTBI(ITSYM,JLR)
      JRR = MULTB(IHSYM,JRR)
      ILR = MULTB(JLR,IRPAMO(I+1))
      IRR = MULTB(JRR,IRPAMO(I+1))
      IF (IVWU(JLR,I,J).NE.0.AND.IVWU(JRR,I,J).NE.0) THEN
C     -------------------------------------------
C     MAKE DETERMINANTS WHICH PASS THESE VERTICES
C     -------------------------------------------
      IVDET(JLR,I,J) = JDET
      IF (IJARR) THEN
         IF (J.NE.NELEC.AND.I.NE.IANNI) THEN
            DO 60 K = IVDET(ILR,I+1,J+1)+1,
     &                IVDET(ILR,I+1,J+1)+JVDET(ILR,I+1,J+1)
            JDET = JDET + 1
            JLIND(JDET) = ILIND(K) + IVWU(ILR,I,J+1)
            JRIND(JDET) = IRIND(K) + IVWU(IRR,I,J+1)
   60       CONTINUE
         ENDIF
         IF (I.NE.ICREA) THEN
            DO 70 K = IVDET(JLR,I+1,J)+1,
     &                IVDET(JLR,I+1,J)+JVDET(JLR,I+1,J)
            JDET = JDET + 1
            JLIND(JDET) = ILIND(K)
            JRIND(JDET) = IRIND(K)
   70       CONTINUE
         ENDIF
      ELSE
         IF (J.NE.NELEC.AND.I.NE.IANNI) THEN
            DO 61 K = IVDET(ILR,I+1,J+1)+1,
     &                IVDET(ILR,I+1,J+1)+JVDET(ILR,I+1,J+1)
            JDET = JDET + 1
            ILIND(JDET) = JLIND(K) + IVWU(ILR,I,J+1)
            IRIND(JDET) = JRIND(K) + IVWU(IRR,I,J+1)
   61       CONTINUE
         ENDIF
         IF (I.NE.ICREA) THEN
            DO 71 K = IVDET(JLR,I+1,J)+1,
     &                IVDET(JLR,I+1,J)+JVDET(JLR,I+1,J)
            JDET = JDET + 1
            ILIND(JDET) = JLIND(K)
            IRIND(JDET) = JRIND(K)
   71       CONTINUE
         ENDIF
      ENDIF
      JVDET(JLR,I,J) = JDET - IVDET(JLR,I,J)
      ENDIF
  100 CONTINUE
  200 CONTINUE
      IF (JDET.EQ.0) RETURN
  300 CONTINUE
C
      IF (NDORB.EQ.0) GOTO 5000
C
C     -----------------------------
C     SPLIT THE LEFT AND RIGHT PATH
C     -----------------------------
      I = ITAIL - 1
      IJARR = (.NOT.IJARR)   
      JDET = 0
      JMIN = IORBR(I)
      JMIN = MAX0(NELEC+I-NORB,0)
      JMAX = MIN0(NELEC-1,I-1)
C     -------------------------------------------------------
C     LOOP OVER VERTICES : CHECK IF THE SPLITTING IS POSSIBLE
C     -------------------------------------------------------
      DO 2200 J = JMAX, JMIN, -1
      JLMIN = (NREP / 2 * MOD(J+1,2)) + 1
      JLMAX = JLMIN + (NREP / 2) - 1
      DO 2100 JLR = JLMIN, JLMAX
      JRR = MULTB(IHSYM,JLR)
      ILR = MULTB(IRPAMO(I+1),JLR)
      IF (IVWU(JLR,I,J).NE.0.AND.IVWU(JRR,I,J+1).NE.0) THEN
C     -------------------------------------------
C     MAKE DETERMINANTS WHICH PASS THESE VERTICES
C     -------------------------------------------
      IVDET(JLR,I,J) = JDET
      IF (IJARR) THEN
         DO 2060 K = IVDET(ILR,I+1,J+1)+1,
     &               IVDET(ILR,I+1,J+1)+JVDET(ILR,I+1,J+1)
         JDET = JDET + 1
         JLIND(JDET) = ILIND(K) + IVWU(ILR,I,J+1)
         JRIND(JDET) = IABS(IRIND(K))
 2060    CONTINUE
      ELSE
         DO 2061 K = IVDET(ILR,I+1,J+1)+1,
     &               IVDET(ILR,I+1,J+1)+JVDET(ILR,I+1,J+1)
         JDET = JDET + 1
         ILIND(JDET) = JLIND(K) + IVWU(ILR,I,J+1)
         IRIND(JDET) = IABS(JRIND(K))
 2061    CONTINUE
      ENDIF
      JVDET(JLR,I,J) = JDET - IVDET(JLR,I,J)
      ENDIF
 2100 CONTINUE
 2200 CONTINUE
      IF (JDET.EQ.0) RETURN
C
C     -------------------
C     MAKE ALL LOOPBODIES
C     -------------------
      DO 3300 I = ITAIL - 2, IHEAD, - 1
      IJARR = (.NOT.IJARR)      
      JDET = 0
      JMIN = IORBR(I)
      JMIN = MAX0(NELEC+I-NORB,0)
      JMAX = MIN0(NELEC-1,I-1)
C     ----------------------------------------------------------------
C     LOOP OVER VERTICES : CHECK IF BOTH (PARALLEL) PATHS ARE POSSIBLE
C     ----------------------------------------------------------------
      DO 3200 J = JMAX, JMIN, -1
      JLMIN = (NREP / 2 * MOD(J+1,2)) + 1
      JLMAX = JLMIN + (NREP / 2) - 1
      DO 3100 JLR = JLMIN, JLMAX
      JRR = MULTB(IHSYM,JLR)
      ILR = MULTB(IRPAMO(I+1),JLR)
      IRR = MULTB(IRPAMO(I+1),JRR)
      IF (IVWU(JLR,I,J).NE.0.AND.IVWU(JRR,I,J+1).NE.0) THEN
C     -------------------------------------------
C     MAKE DETERMINANTS WHICH PASS THESE VERTICES
C     -------------------------------------------
      IVDET(JLR,I,J) = JDET
      IF (IJARR) THEN
         IF (I.NE.IANNI) THEN
            DO 3060 K = IVDET(ILR,I+1,J+1)+1,
     &                  IVDET(ILR,I+1,J+1)+JVDET(ILR,I+1,J+1)
               JDET = JDET + 1
               JLIND(JDET) = ILIND(K) + IVWU(ILR,I,J+1)
               JRIND(JDET) = -1 * ISIGN(1,IRIND(K)) * 
     &                       (IABS(IRIND(K)) + IVWU(IRR,I,J+2))
 3060       CONTINUE
         ENDIF
         IF (I.NE.ICREA) THEN
            DO 3070 K = IVDET(JLR,I+1,J)+1,
     &                  IVDET(JLR,I+1,J)+JVDET(JLR,I+1,J)
               JDET = JDET + 1
               JLIND(JDET) = ILIND(K)
               JRIND(JDET) = IRIND(K)
 3070       CONTINUE
         ENDIF
      ELSE
         IF (I.NE.IANNI) THEN
            DO 3061 K = IVDET(ILR,I+1,J+1)+1,
     &                  IVDET(ILR,I+1,J+1)+JVDET(ILR,I+1,J+1)
               JDET = JDET + 1
               ILIND(JDET) = JLIND(K) + IVWU(ILR,I,J+1)
               IRIND(JDET) = -1 * ISIGN(1,JRIND(K)) * 
     &                       (IABS(JRIND(K)) + IVWU(IRR,I,J+2))
 3061       CONTINUE
         ENDIF
         IF (I.NE.ICREA) THEN
            DO 3071 K = IVDET(JLR,I+1,J)+1,
     &                  IVDET(JLR,I+1,J)+JVDET(JLR,I+1,J)
               JDET = JDET + 1
               ILIND(JDET) = JLIND(K)
               IRIND(JDET) = JRIND(K)
 3071       CONTINUE
         ENDIF
      ENDIF
      JVDET(JLR,I,J) = JDET - IVDET(JLR,I,J)
      ENDIF
 3100 CONTINUE
 3200 CONTINUE
      IF (JDET.EQ.0) RETURN
 3300 CONTINUE
C     --------------------------------------------------
C     CLOSE THE LOOPS; PUT THE INDICES ON THE RIGHT SIDE
C     --------------------------------------------------
      I = IHEAD - 1
      JDET = 0
      IVERT = 0
      JMIN = IORBR(I)
      JMIN = MAX0(NELEC+I-NORB,0)
      JMAX = MIN0(NELEC,I)
C     -----------------------------------------------------
C     LOOP OVER VERTICES : CHECK IF THE CLOSURE IS POSSIBLE
C     -----------------------------------------------------
      DO 4200 J = JMAX, JMIN, -1
      JLMIN = (NREP / 2 * MOD(J+1,2)) + 1
      JLMAX = JLMIN + (NREP / 2) - 1
      DO 4100 JLR = JLMIN, JLMAX
      IRR = MULTB(IRPAMO(I+1),JLR)
      IF (IVWU(JLR,I,J).NE.0.AND.JVDET(JLR,I+1,J).GT.0) THEN
C     ----------------------------------------
C     MAKE DETERMINANTS WHICH PASS THIS VERTEX
C     ----------------------------------------
      IVERT = IVERT + 1
      JVECL(1,IVERT) = JDET
      JVECL(2,IVERT) = IVWU(JLR,I,J)
      IF (NDORB.GT.0) THEN
      DO 4060 K = IVDET(JLR,I+1,J)+1,
     &               IVDET(JLR,I+1,J)+JVDET(JLR,I+1,J)
         JDET = JDET + 1
         JLIND(JDET) = ILIND(K)
         JRIND(JDET) = ISIGN(1,IRIND(K)) *
     &                   (IABS(IRIND(K)) + IVWU(IRR,I,J+1))
 4060 CONTINUE
      ELSE
      DO 4080 K = IVDET(JLR,I+1,J)+1,
     &               IVDET(JLR,I+1,J)+JVDET(JLR,I+1,J)
         JDET = JDET + 1
         JLIND(JDET) = IABS(IRIND(K)) + IVWU(IRR,I,J+1)
         JRIND(JDET) = ISIGN(1,IRIND(K)) * ILIND(K)
 4080 CONTINUE
      ENDIF
      ENDIF
 4100 CONTINUE
 4200 CONTINUE
C
      IVERT = IVERT + 1
      JVECL(1,IVERT) = JDET
      JVECL(2,IVERT) = 0
      NVERT = IVERT
      RETURN
C
 5000 CONTINUE
C     -------------------------------
C     NO LOOP BODIES IF IHEAD = ITAIL
C     -------------------------------
      I = IHEAD - 1
      JDET = 0
      IVERT = 0
      IF(I.EQ.IANNI) GOTO 5300
      JMIN = IORBR(I)
      JMIN = MAX0(NELEC+I-NORB,0)
      JMAX = MIN0(NELEC-1,I)
C     -----------------------------------------------------
C     LOOP OVER VERTICES : CHECK IF THE CLOSURE IS POSSIBLE
C     -----------------------------------------------------
      DO 5200 J = JMAX, JMIN, -1
        JLMIN = (NREP / 2 * MOD(J+1,2)) + 1
        JLMAX = JLMIN + (NREP / 2) - 1
        DO 5100 JLR = JLMIN, JLMAX
          ILR = MULTB(IRPAMO(I+1),JLR)
          IF (IVWU(JLR,I,J).NE.0.AND.JVDET(ILR,I+1,J+1).NE.0) THEN
C           ----------------------------------------
C           MAKE DETERMINANTS WHICH PASS THIS VERTEX
C           ----------------------------------------
            IVERT = IVERT + 1
            JVECL(1,IVERT) = JDET
            JVECL(2,IVERT) = IVWU(JLR,I,J)
            KOFF = IVDET(ILR,I+1,J+1)
            DO K = 1,JVDET(ILR,I+1,J+1)
              JLIND(JDET+K) = ILIND(KOFF+K) + IVWU(ILR,I,J+1)
              JRIND(JDET+K) = JLIND(JDET+K)
            ENDDO
            JDET = JDET + JVDET(ILR,I+1,J+1)
          ENDIF
 5100   CONTINUE
 5200 CONTINUE
 5300 CONTINUE
      IVERT = IVERT + 1
      JVECL(1,IVERT) = JDET
      JVECL(2,IVERT) = 0
      NVERT = IVERT
C
      RETURN
      END
C
      SUBROUTINE CIMEM 
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ========================
C     GENERATE MEMORY POINTERS 
C     ========================
C
#include "param.inc"
#include "symm.inc"
#include "general.inc"
#include "mempoint.inc"
C
      MAXB = NROOTS + 1
      MAXC = NROOTS
      NORBC = NORB**3
C
C     store een deel van de common block data
C
      ITNR = IIJR + NORB**2 / INTOWP(1) + MOD(NORB**2,INTOWP(1))
      ITNR = ITNR + MOD(MDET*(NREP/2),INTOWP(1))
      IF (REALAR) THEN
         ITNI = ITNR 
         IIJE = ITNI + NORBC
         IGMR = IIJE + NORBC
         IGMI = IGMR
      ELSE
         ITNI = ITNR + NORBC
         IIJE = ITNI + NORBC
         IGMR = IIJE + NORBC
         IGMI = IGMR + NORBC
      ENDIF
C
      IPDG = IGMI + NORBC
C
      IPDR = IPDG + NDET
      IF (REALAR) THEN
         IPDI = IPDR
         IPVR = IPDR + MDET*MAXC
         IPVI = IPVR
         IPSR = IPVR + NDET*MAXB
         IPSI = IPSR
      ELSE
         IPDI = IPDR + MDET*MAXC
         IPVR = IPDI + MDET*MAXC
         IPVI = IPVR + NDET*MAXB
         IPSR = IPVI + NDET*MAXB
         IPSI = IPSR + NDET*MAXC
      ENDIF
C   
      IPIL = IPSI + NDET*MAXC
      NDIM = NDET / INTOWP(1) + MOD(NDET,INTOWP(1))
      IPIR = IPIL + NDIM
C     -----------------------------------------------------------
C     USE FREE MEMORY TO STORE ONE ELECTRON COUPLING COEFFICIENTS
C     -----------------------------------------------------------
      NDIM = 5*NDIM
      IF ((IPIR+5*NDIM).GT.MAXCORE)  THEN
        WRITE(6,*) '* NDET/MAXCORE :',NDET,MAXCORE
        WRITE(6,*) '* Program needs :',5*NDIM+IPIR
        CALL Q2JOB (6,'MEMORY PROBLEMS',' ',0)
      ENDIF
C
      LREST = MAXCORE - (IPIR+5*(NDET/INTOWP(1)+MOD(NDET,INTOWP(1))))
      NORBD = (LREST - MOD(LREST,4)) / 4
      NELEM = NORB * NORB
      LNORBC = (NORBD - MOD(NORBD,NELEM)) / NELEM
      LNORBC = LNORBC * INTOWP(1) 
      IF (LNORBC.LT.NORB) THEN
         WRITE(6,*) 'LNORBC/NORB :',LNORBC,NORB
         CALL Q2JOB (6,'MEMORY PROBLEMS',' ',0)
      ENDIF
C
      IPJL = IPIR + NDET / INTOWP(1) + MOD(NDET,INTOWP(1))
      IPJR = IPJL + NDET / INTOWP(1) + MOD(NDET,INTOWP(1)) + NORBD
      IPJV = IPJR + NDET / INTOWP(1) + MOD(NDET,INTOWP(1)) + NORBD
C
      RETURN
C
      ENTRY SORTMEM
C
C     MEMORY WILL FIRST BE USED AS SCRATCH SPACE FOR
C     SORTING OF INTEGRALS:
C
C     MEMORY IN 8 BYTE WORDS
C     M1 : LENGTH OF BLOCKS
C     M2 : NUMBER OF BLOCKS
C     M3 : MEMORY ACTUALLY USED
C
      M1 = NORB * NORB
      MORB = NORB
      M1A = MORB * MORB
C
      MCORE = (MAXCORE - 4*M1A)/2
      IF (BREIT) THEN
         M2 = MIN0(MCORE/(2*M1),M1)
         M3 = 2 * M1 * M2
      ELSE
         M2 = MIN0(MCORE/M1,M1)
         M3 = M1 * M2
      ENDIF
      ITR = 1
      ITI = ITR + M1A
      IIK = ITI + M1A
      IIL = IIK + M1A
      ITIR = IIL + M1A
      ITII = ITIR + M3
C
C     Also set pointers for vertex weights, etc.
C
      IVWS = 2*MXREP * (NORB + 1) * (NELEC + 1)
      IVWW = IVWS / INTOWP(1) + MOD(IVWS,INTOWP(1))
      IVWI = IVWW * INTOWP(1)
      IVW1 = 1               
      IIJR = IVW1 + 3*IVWW
C
      RETURN
C
      END
C
      SUBROUTINE CopyFF (IUN1,IUN2)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =====================
C     Copies File1 to File2
C     =====================
C
      Character*132 Card
C
   10 READ (IUN1,'(A)',END=20,ERR=10001) CARD
      WRITE (IUN2,'(A)') CARD
      IF (CARD(1:3).EQ.'END') GOTO 20
      GOTO 10
   20 CONTINUE
      RETURN
10001 CALL QUIT('Cannot read Unit')
      END
C
      SUBROUTINE CRITER (ENDRUN)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =================================================
C     CHECKS IF ONE OF ENDRUN CONDITIONS IS FULFILLED :
C     1) CONVERGENCE REACHED
C     2) MAXITER REACHED
C     3) MAXIMUM CPU TIME USED
C     =================================================
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
#include "iterat.inc"
C
      LOGICAL ENDRUN
      DIMENSION MSEL(N8)
C
      ENDRUN = .FALSE.
      CALL CPUUSED (CPUTOT)
      DO 10 IROOT = 1, NROOTS
      MSEL(IROOT) = NSEL(IROOT)
      CONV(IROOT) = EIGEN(NSEL(IROOT)) + ECORE - EIGEN2(IROOT)
      EIGEN2(IROOT) = EIGEN(NSEL(IROOT)) + ECORE
C
C  Determine whether the Davdidson correction can be calculated
C
      IF (ITER.EQ.1) THEN
         EIGEN1(IROOT) = EIGEN2(IROOT)
C
C  The start vector should be within the reference space
C
         IF (1.D0-REFVEC(IROOT).LT.1.D-12) THEN
             DCORR(IROOT) = .TRUE.
         ELSE
             DCORR(IROOT) = .FALSE.
         ENDIF
C
C  Do not calculate correction in non SDCI cases and silly cases
C
         IF (MAXE3.NE.2.OR.NELEC.EQ.2) DCORR(IROOT) = .FALSE.
C
C  The reference space must be pre-diagonalized
C
         IF (REFRDU(IROOT).GT.1.D-12) DCORR(IROOT) = .FALSE.
      ENDIF
C
      ECORR = EIGEN2(IROOT) - EIGEN1(IROOT)
      WRITE (6,1000) ITER,IROOT,NSEL(IROOT),
     &               EIGEN2(IROOT),ECORR,CONV(IROOT),RESIDU(IROOT)
C
      CALL FLSHFO (6)
C
      IF (ABS(CONV(IROOT)).LT.CONVERE.AND.ITER.NE.1) THEN
         MSEL(IROOT) = 0
      ENDIF

      IF (ABS(RESIDU(IROOT)).LT.CONVERR) THEN
         MSEL(IROOT) = 0
      ENDIF
   10 CONTINUE
C
      JROOT = 0
      DO 20 IROOT = 1, NROOTS
         IF (MSEL(IROOT).GT.0) THEN
            JROOT = JROOT + 1
          ENDIF
   20     CONTINUE
C
C     CHECK CRITERIA FOR ENDRUN :
C
      IF (JROOT.EQ.0) THEN
         WRITE (6,1010) 'Converged'
         ENDRUN = .TRUE.
      ENDIF
C
      IF (CPUTOT.GE.CPUMAX) THEN
         WRITE (6,1010) 'CPU time exceeded'
         ENDRUN = .TRUE.
      ENDIF
C
      IF (ITER.EQ.MAXITER) THEN
         WRITE (6,1010) 'Maxiter reached'
         ENDRUN = .TRUE.
      ENDIF
C
      IF (NVEC+NROOTS.GE.N2) THEN
         WRITE (6,1010) 'NVEC too large, increase N2'
         ENDRUN = .TRUE.
      ENDIF
C
      IF (ENDRUN) THEN
         DO IROOT = 1, NROOTS
	  IF (DCORR(IROOT)) THEN
            ECORR = EIGEN2(IROOT)-EIGEN1(IROOT)
	    DAVCOR = (1.D0-REFVEC(IROOT))*ECORR
            WRITE (6,1020) IROOT,EIGEN2(IROOT),ECORR,DAVCOR,
     &                     EIGEN2(IROOT)+DAVCOR,RESIDU(IROOT)
	  ELSE
	    WRITE (6,1021) IROOT,EIGEN2(IROOT),EIGEN2(IROOT)-
     &                     EIGEN1(IROOT),(1.D0-REFVEC(IROOT)),
     &                     RESIDU(IROOT)
	  ENDIF
         ENDDO
         WRITE (6,1030) CPUTOT
         RETURN
      ENDIF
C
C     CONTINUE OPTIMIZATION OF NON-CONVERGED ROOTS :
C
      JROOT = 0
      DO 30 IROOT = 1, NROOTS
         IF (MSEL(IROOT).GT.0) THEN
            JROOT = JROOT + 1
             NSEL(JROOT) = MSEL(IROOT)
          ENDIF
   30     CONTINUE
      MVEC = NVEC
      NVEC = NVEC + JROOT
C
      RETURN
C
 1000 FORMAT (3(1X,I4),1X,G22.12,1X,G18.12,2(1X,G10.4))
 1010 FORMAT (//' End of iterations : ',A25)
 1020 FORMAT (/' --- Root',I4,' ---'/
     &'@ (MR)CI-SD Energy :    ',F30.12/
     &'@ Correlation energy :  ',F30.12/
     &'@ Davidson correction : ',F30.12/
     &'@ (MR)CI-SD+Q energy:   ',F30.12/
     &'@ Residue :             ',F30.12)
 1021 FORMAT (/' --- Root',I4,' ---'/
     &'@ CI Total Energy :       ',F30.12/
     &'@ Relative to reference : ',F30.12/
     &'@ 1 - Reference weight :  ',F30.12/
     &'@ Residue :               ',F30.12)
 1030 FORMAT (/' Total CPU time :     ',F10.3)
C
      END
C
      SUBROUTINE DENSMT (DR,DI,VECTR,VECTI,ILIND,IRIND,
     &                   JLIND,JRIND,JVECL,IJREC,
     &                   DENSR,DENSI,IVWA,IVEC,DD)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =====================================
C     CALCULATES FIRST ORDER DENSITY MATRIX
C     =====================================
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
#include "files.inc"
#include "mempoint.inc"
C
      DIMENSION DENSR(NORB,NORB),DENSI(NORB,NORB),DD(NORB)
      DIMENSION ILIND(*),IRIND(*),JLIND(*),JRIND(*),IVWA(*)
      DIMENSION DR(MDET,*),DI(MDET,*),VECTR(NDET,*),VECTI(NDET,*)
      DIMENSION JVECL(2,*),IJREC(NORB,NORB)
C
C     -----------------------------------------
C     LOOP OVER SINGLE PARTICLE REPRESENTATIONS
C     -----------------------------------------
      DO 1200 I = 1, NORB
         DO 1100 J = 1, I   
            DENSR(I,J) = 0.D0
            IF (.NOT.REALAR) DENSI(I,J) = 0.D0
            IF (IRPAMO(J).NE.IRPAMO(I)) GOTO 1100
            IREC = 0
            IF (IGENEX.EQ.2.AND.IJREC(I,J).LT.0) IREC = -IJREC(I,J)
            CALL GETIDS (I,J,0,0,IREP,ND,NV,ILIND,IRIND,
     &                   JLIND(1),JRIND(1),JVECL(1,1),
     &                   IVWA(1),IVWA(1+IVWI),IREC)
            IVERT = 0
C     --------------------------------
C     START "WHILE LOOP" OVER VERTICES
C     --------------------------------
   30      IVERT = IVERT + 1
           JD1 = JVECL(1,IVERT) + 1 
           JD2 = JVECL(1,IVERT+1) 
           NDH = IABS(JVECL(2,IVERT))
           IF (NDH.EQ.0) GOTO 1100
C$DIR FORCE_PARALLEL
           DO 500 ID = JD1, JD2
              IS = ISIGN(1,JRIND(ID))
              ILD = JLIND(ID) - 1
              IRD = IABS(JRIND(ID)) - 1
              IF (REALAR) THEN
                 DO 398 IDH = 1, NDH
                    IL = ILD + IDH
                    IR = IRD + IDH
                    DENSR(I,J) = DENSR(I,J) +
     &              VECTR(IL,IVEC) * IS * VECTR(IR,IVEC) 
  398            CONTINUE
              ELSE
                 DO 400 IDH = 1, NDH
                    IL = ILD + IDH
                    IR = IRD + IDH
                    DENSR(I,J) = DENSR(I,J) +
     &              VECTR(IL,IVEC) * IS * VECTR(IR,IVEC) +
     &              VECTI(IL,IVEC) * IS * VECTI(IR,IVEC)
                    DENSI(I,J) = DENSI(I,J) +
     &              VECTR(IL,IVEC) * IS * VECTI(IR,IVEC) -
     &              VECTI(IL,IVEC) * IS * VECTR(IR,IVEC)
  400            CONTINUE
              ENDIF
  500         CONTINUE
C             ----------------
C             GOTO NEXT VERTEX
C             ----------------
           GOTO 30
 1100    CONTINUE
 1200 CONTINUE
      DO 1220 I = 1, NORB
         DO 1210 J = 1, I - 1
            DENSR(J,I) = DENSR(I,J)
 1210       CONTINUE
 1220    CONTINUE
      IF (.NOT.REALAR) THEN
         DO 1221 I = 1, NORB
            DO 1211 J = 1, I - 1
               DENSI(J,I) = - DENSI(I,J)
 1211       CONTINUE
 1221    CONTINUE
      ENDIF
      DO 1230 I = 1, NORB
 1230    DD(I) = DENSR(I,I)
      MVEC = NVEC
      RETURN
      END
C
      SUBROUTINE DETIND (ND,IDIND,JDET,JELEC,JSYM,IVWU)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =========================================
C     CALCULATE ADRESSES OF THE STRINGS IN JDET
C     =========================================
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
C
      DIMENSION IDIND(*),JDET(*),JELEC(*),JSYM(*)
      DIMENSION IVWU(2*MXREP,0:NORB,0:NELEC)
C
C
#if defined (SYS_CRAY) || defined (SYS_T3D)
      IAND (I, J) = AND (I, J)
      IOR (I, J) = OR (I, J)
      NOT(I) = COMPL (I)
      ISHFT (I, J) = SHIFT (I, J)
#endif
C
      DO 10 ID = 1, ND
      JELEC(ID) = NORBR(1)
      IDIND(ID) = 1
      JSYM(ID) = IRA1
   10 CONTINUE
C
      DO 30 I = 1, NORBR(2)
      MASK = ISHFT(1,I-1)
      II = NORBR(1) + I
      DO 20 ID = 1, ND
      IF (IAND(JDET(ID),MASK).NE.0) THEN
         JELEC(ID) = JELEC(ID) + 1
         JSYM(ID) = MULTB(IRPAMO(II),JSYM(ID))
         IDIND(ID) = IDIND(ID) + IVWU(JSYM(ID),II-1,JELEC(ID))
      ENDIF
   20 CONTINUE
   30 CONTINUE
C
      DO 40 ID = 1, ND
      IF (JSYM(ID).NE.IREP)
     & CALL Q2JOB (2,'TRIAL VECTORS HAVE WRONG SYMMETRY',' ',0)
   40 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE DIAG
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ----------------------
C     DIAGONALIZES CI MATRIX
C     ----------------------
C
#include "param.inc"
#include "general.inc"
#include "iterat.inc"
C
      DIMENSION WORKR(N2,N2),WORKI(N2,N2),WORK2(N2),WORK3(N2),WORK4(N2)
C
      DO 10 I = 1, NVEC
      DO 10 J = 1, NVEC
   10    EVECR (I,J) = CMTR (I,J)
C
      IF (REALAR) THEN
         CALL RS (N2,NVEC,EVECR,EIGEN,1,EVECI,WORKR,WORKI,IERR)
         DO I = 1, NVEC
         DO J = 1, NVEC
            EVECR(J,I) = EVECI(J,I)
         ENDDO
         ENDDO
C     ----
      ELSE
C     ----
         DO 11 I = 1, NVEC
         DO 11 J = 1, NVEC
   11       EVECI (I,J) = CMTI (I,J)
         CALL CH (N2,NVEC,EVECR,EVECI,EIGEN,1,WORKR,WORKI,WORK2,
     $            WORK3,WORK4,IFAIL)
         DO 20 I = 1, NVEC
         DO 20 J = 1, NVEC
         EVECR(I,J) = WORKR(I,J)
         EVECI(I,J) = WORKI(I,J)
   20    CONTINUE
         IF (IFAIL.NE.0) CALL QUIT('**ERROR-CANNOT DIAGONALISE CMTR**')
C     -----
      ENDIF
C     -----
      RETURN
      END
C
      SUBROUTINE DIAGELM (DR,DGEL,TINTR,TINTI,IJEX,GMATR,GMATI,
     &                    ILIND,IRIND,JLIND,JRIND,JVECL,IVWA,
     &                    IJREC)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =========================================
C     CALCULATES DIAGONAL ELEMENTS OF CI MATRIX
C     =========================================
C
#include "param.inc"
#include "general.inc"
#include "files.inc"
#include "symm.inc"
#include "mempoint.inc"
      DIMENSION DGEL(*),DR(*),TINTR(*),TINTI(*),GMATR(*),GMATI(*)
      DIMENSION IVWA(*),IJREC(NORB,NORB)
      DIMENSION IJEX(*),ILIND(*),IRIND(*),JLIND(*),JRIND(*),JVECL(2,*)
      LOGICAL ENDFIL
      real*8, allocatable :: DGER(:),DGEI(:),DI(:)
C
      REWIND (MRCTWOE)
      CALL ZZERO(NDET,DGEL)
      IF (PROPER) THEN
        allocate(DGER(NDET))
        allocate(DGEI(NDET))
        allocate(DI(NDET))
        CALL ZZERO(NDET,DGER)
        CALL ZZERO(NDET,DGEI)
      ENDIF
C     -------------------------------------------
C     THE COULOMB TYPE OF CONTRIBUTIONS : (II|JJ)
C     -------------------------------------------
 3000 CALL TWOLIN (ENDFIL,TINTR,TINTI,IJEX,IJ,NIJ)
      IF (ENDFIL) GOTO 3001
        I = (IJ-1)/NORB + 1
        J = MOD(IJ-1,NORB) + 1
        IF (I.NE.J) GOTO 3000
          CALL ZZERO (MDET,DR(1))
          IF(PROPER) THEN
            IF (.NOT.REALAR) CALL ZZERO (NDET,DI)
          ENDIF
          IF (IJREC(I,J).GT.0) THEN
             IVERT = IJREC(I,J) - 1
             ICOREX = 0
          ELSE
             IREC = -IJREC(I,J)
             CALL GETIDS (I,J,0,0,IREP,ND,NV,ILIND,IRIND,
     &                    JLIND(NCOREX),JRIND(NCOREX),
     &                    JVECL(1,NTVERT),IVWA(1),IVWA(1+IVWI),
     &                    IREC)
             IVERT = NTVERT - 1
             ICOREX = NCOREX - 1
          ENDIF
C         --------------------------------
C         START "WHILE LOOP" OVER VERTICES
C         --------------------------------
   30       CONTINUE
            IVERT = IVERT + 1
            JD1 = JVECL(1,IVERT) + 1 + ICOREX
            JD2 = JVECL(1,IVERT+1) + ICOREX
            NDH = IABS(JVECL(2,IVERT))
            IF (NDH.EQ.0) GOTO 1001
            DO 400 ID = JD1, JD2
               ILD = JLIND(ID) - 1
               DO 300 IDH = 1, NDH
                  IL = ILD + IDH
                  DGEL(IL) = DGEL(IL) + GMATR(IJ)
                  DR(IL) = 1.D0
                  IF (PROPER) THEN
                    DGER(IL) = DGER(IL) + GMATR(IJ)
                    IF (.NOT.REALAR) DGEI(IL) = DGEI(IL) + GMATI(IJ)
                    DI(IL) = 1.D0
                  END IF
  300             CONTINUE
  400          CONTINUE
C           --------------
C           DO NEXT VERTEX
C           --------------
            GOTO 30
C
 1001     DO 2000 IX = 1, NIJ
            KL = IJEX(IX)
            K = (KL-1)/NORB + 1
            L = MOD(KL-1,NORB) + 1
            IF (K.EQ.L) THEN
               IF (IJREC(K,L).GT.0) THEN
                  IVERT = IJREC(K,L) - 1
                  ICOREX = 0
               ELSE
                  IREC = -IJREC(K,L)
                  CALL GETIDS (K,L,I,J,IREP,ND,NV,ILIND,IRIND,
     &                         JLIND(NCOREX),JRIND(NCOREX),
     &                         JVECL(1,NTVERT),IVWA(1),IVWA(1+IVWI),
     &                         IREC)
                  IVERT = NTVERT - 1
                  ICOREX = NCOREX - 1
               ENDIF
 1030          IVERT = IVERT + 1
               JD1 = JVECL(1,IVERT) + 1 + ICOREX
               JD2 = JVECL(1,IVERT+1) + ICOREX
               NDH = IABS(JVECL(2,IVERT))
               IF (NDH.EQ.0) GOTO 2000
               DO 1400 ID = JD1, JD2
                 ILD = JLIND(ID) - 1
                 DO 1300 IDH = 1, NDH
                   IL = ILD + IDH
                   DGEL(IL) = DGEL(IL) + TINTR(IX) * DR(IL)
                   IF (PROPER) THEN
                     DGER(IL) = DGER(IL) + TINTR(IX)*DR(IL)
                     IF (.NOT.REALAR) 
     &               DGEI(IL) = DGEI(IL) + TINTI(IX)*DI(IL)
                   END IF
 1300            CONTINUE
 1400          CONTINUE
C               --------------
C               DO NEXT VERTEX
C               --------------
                GOTO 1030
C            -------------------
C            END IF BLOCK K.EQ.L
C            -------------------
            ENDIF
 2000       CONTINUE
C     -----------------------
C     GET NEXT INTEGRAL BLOCK
C     -----------------------
      GOTO 3000
C
 3001 REWIND (MRCTWOE)
C
C     -------------------------------------------
C     THE EXCHANGE TYPE OF CONTRIBUTIONS : (IJ|JI)
C     -------------------------------------------
C
 5000 CALL TWOLIN (ENDFIL,TINTR,TINTI,IJEX,IJ,NIJ)
      IF (ENDFIL) GOTO 5001
      I = (IJ-1)/NORB + 1
      J = MOD(IJ-1,NORB) + 1
      JI = (J-1) * NORB + I
      IF (I.EQ.J) GOTO 5000
      DO 4000 IX = 1, NIJ
         IF (IJEX(IX).EQ.JI) THEN
           IF (IJREC(I,J).GT.0) THEN
              IVERT = IJREC(I,J) - 1
              ICOREX = 0
           ELSE
              IREC = -IJREC(I,J)
              CALL GETIDS (I,J,0,0,IREP,ND,NV,ILIND,IRIND,
     &                     JLIND(NCOREX),JRIND(NCOREX),
     &                     JVECL(1,NTVERT),IVWA(1),IVWA(1+IVWI),
     &                     IREC)
              IVERT = NTVERT - 1
              ICOREX = NCOREX - 1
           ENDIF
C             --------------------------------
C             START "WHILE LOOP" OVER VERTICES
C             --------------------------------
 2030        IVERT = IVERT + 1
             JD1 = JVECL(1,IVERT) + 1 + ICOREX
             JD2 = JVECL(1,IVERT+1) + ICOREX
             NDH = IABS(JVECL(2,IVERT))
             IF (NDH.EQ.0) GOTO 5000
             DO 2400 ID = JD1, JD2
               IRD = IABS(JRIND(ID)) - 1
               DO 2300 IDH = 1, NDH
                 IR = IRD + IDH
                 DGEL(IR) = DGEL(IR) + TINTR(IX)
                 IF (PROPER) THEN
                   DGER(IR) = DGER(IR) + TINTR(IX)
                   IF (.NOT.REALAR) DGEI(IR) = DGEI(IR) + TINTI(IX)
                 END IF
 2300          CONTINUE
 2400        CONTINUE
C            --------------
C            DO NEXT VERTEX
C            --------------
             GOTO 2030
C        ---------------------
C        END IF BLOCK KL.EQ.JI
C        ---------------------
         ENDIF
 4000    CONTINUE
C     -----------------------
C     GET NEXT INTEGRAL BLOCK
C     -----------------------
      GOTO 5000
 5001 CONTINUE
      IF (PROPER) THEN
        WRITE(*,25)NAMEA
        DO IROOT = 1,NROOTS
          WRITE(*,50)IROOT,INDMIN(IROOT),DGER(INDMIN(IROOT)),
     &                                   DGEI(INDMIN(IROOT))
        END DO
        deallocate(DGER)
        deallocate(DGEI)
        deallocate(DI)
      END IF
      RETURN
   25 FORMAT (/' UN-CORRELATED PROPERTY: ',A8,14X,'REAL',14X,
     &          'IMAGINARY')
   50 FORMAT (/' Root',I4,' is determinant',I4,' Exp. Value:',
     &           F16.12,4X,F16.12)
      END
C
      BLOCK DATA DIRDAT
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "iterat.inc"
CMKN
      DATA EIGEN1/N8 * 0.D0/
      DATA MRCONEE,MDCINT,MDBINT,MRCTRIV,MRCFINV,MRCEXCS,MRCVECS,
     &     MRCTWOE,MDPROP /8,9,10,11,12,13,14,15,16/
      DATA FNAM /'MRCONEE','MDCINT','MDBINT','MDTRIV','MRCFINV',
     &     'MRCEXCS','MRCVECS','MRCTWOE','MDPROP'/
CMKN
      END
C
      SUBROUTINE EXPAND(DR,DI,VECTR,VECTI,SIGMAR,SIGMAI,
     &                  ILIND,IRIND,JLIND,JRIND,JVECL,
     &                  TINTR,TINTI,IJEX,GMATR,GMATI,
     &                  IVWU,IVWA,IJREC)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =======================
C     CALCULATES SIGMA VECTOR
C     =======================
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
#include "files.inc"
#include "mempoint.inc"
C
      DIMENSION ILIND(*),IRIND(*),JLIND(*),JRIND(*)
      DIMENSION JVECL(2,*),IVWA(*),IJREC(NORB,NORB),IJEX(*)
      DIMENSION IVWU(2*MXREP,0:NORB,0:NELEC)
      DIMENSION DR(MDET,*),DI(MDET,*),VECTR(NDET,*),VECTI(NDET,*)
      DIMENSION SIGMAR(NDET,*),SIGMAI(NDET,*),TINTR(*),TINTI(*)
      DIMENSION GMATR(*),GMATI(*)
      LOGICAL ENDFIL
C
      INTERN = NORBR(1) + NORBR(2)
      IF (MAXE3.GT.2) INTERN = NORB
      REWIND (MRCTWOE)
      DO 10 IVEC = 1, NVEC - MVEC
         CALL ZZERO (NDET,SIGMAR(1,IVEC))
         IF (.NOT.REALAR) CALL ZZERO (NDET,SIGMAI(1,IVEC))
   10    CONTINUE
C     -------------------------------------------------------------------
C     LOOP OVER THE ONE-ELECTRON INTEGRALS FOR IRPAMO(I).EQ.IRPAMO(J)
C     LOOP OVER THE TWO-ELECTRON INTEGRALS : BLOCKS (KL|IJ) WITH IJ FIXED
C     -------------------------------------------------------------------
 5000 CALL TWOLIN (ENDFIL,TINTR,TINTI,IJEX,IJ,NIJ)
      IF (ENDFIL) GOTO 5001
      I = (IJ-1)/NORB + 1
      J = MOD(IJ-1,NORB) + 1
      DO 20 IVEC = 1, NVEC - MVEC
         CALL ZZERO (MDET,DR(1,IVEC))
         IF (.NOT.REALAR) CALL ZZERO (MDET,DI(1,IVEC))
   20    CONTINUE
      IF (IJREC(I,J).GT.0) THEN
         IVERT = IJREC(I,J) - 1
         ICOREX = 0
      ELSE
         IREC = -IJREC(I,J)
         CALL GETIDS (I,J,0,0,IREP,ND,NV,ILIND,IRIND,
     &                JLIND(NCOREX),JRIND(NCOREX),
     &                JVECL(1,NTVERT),IVWU,
     &                IVWA(1),IREC)
         IVERT = NTVERT - 1
         ICOREX = NCOREX - 1
      ENDIF
C     --------------------------------
C     START "WHILE LOOP" OVER VERTICES
C     --------------------------------
   40   IVERT = IVERT + 1
        JD1 = JVECL(1,IVERT) + 1 + ICOREX
        JD2 = JVECL(1,IVERT+1) + ICOREX
        NDH = IABS(JVECL(2,IVERT))
        IF (NDH.EQ.0) GOTO 1001
C$DIR FORCE_PARALLEL
        DO 1600 IVEC = 1, NVEC - MVEC
           DO 1500 ID = JD1, JD2
              IS = ISIGN(1,JRIND(ID))
              ILD = JLIND(ID) - 1
              IRD = IABS(JRIND(ID)) - 1
              DO IDH = 1, NDH
                 IL = ILD + IDH
                 IR = IRD + IDH
                 DR(IL,IVEC) = VECTR(IR,IVEC) * IS
              ENDDO
              IF (.NOT.REALAR) THEN
                 DO IDH = 1, NDH
                    IL = ILD + IDH
                    IR = IRD + IDH
                    DI(IL,IVEC) = VECTI(IR,IVEC) * IS
                 ENDDO
              ENDIF
              IF (IRPAMO(I).EQ.IRPAMO(J)) THEN
                 DO IDH = 1 , NDH 
                    IL = ILD + IDH
                    SIGMAR(IL,IVEC) = SIGMAR(IL,IVEC) +
     &              DR(IL,IVEC) * GMATR(IJ) 
                    IF (.NOT.REALAR) THEN
                       SIGMAR(IL,IVEC) = SIGMAR(IL,IVEC) -
     &                 DI(IL,IVEC) * GMATI(IJ)
                       SIGMAI(IL,IVEC) = SIGMAI(IL,IVEC) +
     &                 DR(IL,IVEC) * GMATI(IJ) +
     &                 DI(IL,IVEC) * GMATR(IJ)
                    ENDIF
                 ENDDO
              ENDIF
 1500      CONTINUE
 1600   CONTINUE
C       ----------------
C       GOTO NEXT VERTEX
C       ----------------
        GOTO 40
C     -------------------------------------------------------------------
C     LOOP OVER THE SECOND EXCITATIONS : TWO-ELECTRON PART OF HAMILTONIAN
C     -------------------------------------------------------------------
 1001 DO 3000 IX = 1, NIJ
         KL = IJEX(IX)
         K = (KL-1)/NORB + 1
         L = MOD(KL-1,NORB) + 1
C        -------------------------------------------------------------
C        FIND OUT WHICH DETERMINANTS INTERACT :
C        1) SPECIAL CASE OF 3 OR MORE EXTERNAL LABELS : ONLY 1 VERTEX
C           K=>I>INTERN; L>INTERN, I<>L
C           J,K,L>INTERN L=<INTERN
C        2) ADDRESSES ARE IN CORE MEMORY
C        3) ADDRESSES ARE ON FILE OR NEED TO BE CALCULATED
C        -------------------------------------------------------------
         IF (I.GT.INTERN.AND.L.GT.INTERN.AND.I.NE.L) THEN
            ICOREX = NCOREX - 1
            IVERT = NTVERT - 1
            IOSYM = IRPAMO(I)
            KOSYM = IRPAMO(K)
            LOSYM = IRPAMO(L)
            IV1SYM = MULTBI(KOSYM,IREP)
            KREP = MULTB(LOSYM,IV1SYM)
            IVSYM = MULTBI(IOSYM,IV1SYM)
            IV2SYM = MULTB(LOSYM,IVSYM)
            JVECL(1,NTVERT) = 0
            JVECL(1,NTVERT+1) = 1
            JVECL(2,NTVERT) = IVWU(IVSYM,INTERN,NELEC-2)
            IF (K.EQ.I) JVECL(2,NTVERT) = 0
            JVECL(2,NTVERT+1) = 0
            JRIND(NCOREX) = ISIGN(1,L-I) * (IVWU(IREP,K-1,NELEC) 
     &                      + IVWU(IV1SYM,I-1,NELEC-1)+1)
            IF (I.LT.L) THEN
               JLIND(NCOREX) = IVWU(KREP,L-1,NELEC) +
     &                    IVWU(IV1SYM,I-1,NELEC-1) + 1
            ELSE
               JLIND(NCOREX) = IVWU(KREP,I-1,NELEC) +
     &                    IVWU(IV2SYM,L-1,NELEC-1) + 1
            ENDIF
         ELSE IF (I.LE.INTERN.AND.J.GT.INTERN.AND.
     &            K.GT.INTERN.AND.L.GT.INTERN) THEN
            ICOREX = NCOREX - 1
            IVERT = NTVERT - 1
            IOSYM = IRPAMO(I)
            KOSYM = IRPAMO(K)
            LOSYM = IRPAMO(L)
            IVSYM = MULTBI(KOSYM,IREP)
            KREP = MULTB(LOSYM,IVSYM)
            JVECL(1,NTVERT) = 0
            JVECL(1,NTVERT+1) = 1
            JVECL(2,NTVERT) = IVWU(IVSYM,INTERN,NELEC-1)
            JVECL(2,NTVERT+1) = 0
            JRIND(NCOREX) = IVWU(IREP,K-1,NELEC) + 1
            JLIND(NCOREX) = IVWU(KREP,L-1,NELEC) + 1
         ELSE IF (IJREC(L,K).GT.0) THEN
            IVERT = IJREC(L,K) - 1
            ICOREX = 0
         ELSE
            IREC = -IJREC(L,K)
            CALL GETIDS (L,K,I,J,IREP,ND,NV,ILIND,IRIND,
     &                   JLIND(NCOREX),JRIND(NCOREX),
     &                   JVECL(1,NTVERT),IVWU,
     &                   IVWA(1),IREC)
            IVERT = NTVERT - 1
            ICOREX = NCOREX - 1
         ENDIF
C        --------------------------------
C        START "WHILE LOOP" OVER VERTICES
C        --------------------------------
 1030      IVERT = IVERT + 1
           JD1 = JVECL(1,IVERT) + 1 + ICOREX
           JD2 = JVECL(1,IVERT+1) + ICOREX
           NDH = IABS(JVECL(2,IVERT))
           IF (NDH.EQ.0) GOTO 3000
C$DIR FORCE_PARALLEL
           DO 2600 IVEC = 1, NVEC - MVEC
              DO 2500 ID = JD1, JD2
                 IS = ISIGN(1,JRIND(ID))
                 ILD = IABS(JRIND(ID)) - 1
                 IRD = JLIND(ID) - 1
                 DO IDH = 1, NDH
                    IL = ILD + IDH
                    IR = IRD + IDH
                    SIGMAR(IL,IVEC) = SIGMAR(IL,IVEC) +
     &              DR(IR,IVEC) * IS * TINTR(IX) 
                 ENDDO
                 IF (.NOT.REALAR) THEN
                    DO IDH = 1, NDH
                       IL = ILD + IDH
                       IR = IRD + IDH
                       SIGMAR(IL,IVEC) = SIGMAR(IL,IVEC) -
     &                 DI(IR,IVEC) * IS * TINTI(IX)
                       SIGMAI(IL,IVEC) = SIGMAI(IL,IVEC) +
     &                 DR(IR,IVEC) * IS * TINTI(IX) +
     &                 DI(IR,IVEC) * IS * TINTR(IX)
                    ENDDO
                 ENDIF
 2500         CONTINUE
 2600      CONTINUE
C          ----------------
C          GOTO NEXT VERTEX
C          ----------------
           GOTO 1030
 3000 CONTINUE
C     -----------------------
C     GET NEXT INTEGRAL BLOCK
C     -----------------------
      GOTO 5000
 5001 CONTINUE
CMKN
C      ************************************************************
C      HERE IT CREATS THE CIMATRIX IN DETERMINANTAL BASIS
C      Written by (Malaya K. Nayak)
C      ************************************************************
C      WRITE(*,*)'IEOPER=',IEOPER,' MAXITER=',MAXITER,' ITER=',ITER
C      IF (ITER.EQ.1) THEN
C      IF (.NOT.PROPER) THEN
C      WRITE(*,*)'NVEC-MVEC=',NVEC-MVEC,' NDET=',NDET
C      IU=20+IEOPER
C      WRITE(*,*)'IU=',IU
C      OPEN (17, FILE='HCIMAT', STATUS='UNKNOWN')
C      WRITE(17,*)'DETIND VECIND   REAL-PART       IMAG-PART'
C      DO IVEC=1,NVEC-MVEC
C      DO I=1,NDET
C      IF (ABS(SIGMAR(I,IVEC)) .GT. 1.0D-3 .OR.
C     &    ABS(SIGMAI(I,IVEC)) .GT. 1.0D-3) THEN
C      IF (I .EQ. IVEC) THEN
C      WRITE(17,100)I,IVEC,SIGMAR(I,IVEC),SIGMAI(I,IVEC),'DIAGONAL'
C      ELSE
C      WRITE(17,100)I,IVEC,SIGMAR(I,IVEC),SIGMAI(I,IVEC),'   OFF  '
C      END IF
C      END IF
C      END DO
C      END DO
C      ************************************************************
C      ************************************************************
C      ELSE IF(PROPER)THEN
C      WRITE(*,*)'NVEC-MVEC=',NVEC-MVEC,' NDET=',NDET
C      IU=20+IEOPER
C      WRITE(*,*)'IU=',IU
C      OPEN (IU, FILE=NAMEE, STATUS='UNKNOWN')
C      WRITE(IU,*)'DETIND VECIND   REAL-PART       IMAG-PART'
C      DO IVEC=1,NVEC-MVEC
C      EDMR=0.0D0
C      DO I=1,NDET
C      EDMR=EDMR+SIGMAR(I,IVEC)*VECTR(I,IVEC)
C      IF (ABS(SIGMAR(I,IVEC)) .GT. 1.0D-3 .OR.
C     &    ABS(SIGMAI(I,IVEC)) .GT. 1.0D-3) THEN
C      IF (I .EQ. IVEC) THEN
C      WRITE(IU,100)I,IVEC,SIGMAR(I,IVEC),SIGMAI(I,IVEC),'DIAGONAL'
C      ELSE
C      WRITE(IU,100)I,IVEC,SIGMAR(I,IVEC),SIGMAI(I,IVEC),'   OFF  '
C      END IF
C      END IF
C      END DO
C      WRITE(*,200)EDMR
C      END DO
C      *************************************************************
C      ELSE
C      WRITE(*,*)' Error in reading property name ',NAMEE
C      CALL QUIT(' Property name did not matching ')
C      END IF
C      END IF
C      *************************************************************
C      END OF WRITING CI-MATRIX FOR DIFFERENT OPERATORS
C      *************************************************************
CMKN
      IF (PROPER) THEN
        WRITE(*,25)NAMEA
        DO IVEC=1,NVEC-MVEC
          PRPR=0.0D0
          PRPI=0.0D0
        DO I=1,NDET
        IF (REALAR) THEN
          PRPR=PRPR+SIGMAR(I,IVEC)*VECTR(I,IVEC)
        ELSE
          PRPR=PRPR+SIGMAR(I,IVEC)*VECTR(I,IVEC)
     &             +SIGMAI(I,IVEC)*VECTI(I,IVEC)
          PRPI=PRPI+SIGMAR(I,IVEC)*VECTI(I,IVEC)
     &             -SIGMAI(I,IVEC)*VECTR(I,IVEC)
        END IF
      END DO
      WRITE(*,50) IVEC,IVEC,PRPR,PRPI
      END DO
      END IF
C     **************************************************************
C
      DO IVEC = 1, NVEC - MVEC
      CALL PUTSIG (MVEC+IVEC,SIGMAR(1,IVEC),SIGMAI(1,IVEC))
      ENDDO
C
      RETURN
   25 FORMAT (/' CORRELATED PROPERTY: ',A8,17X,'REAL',14X,'IMAGINARY')
   50 FORMAT (/' State',I4,' is vector',I4,' Exp. Value:',2F20.12)
   51 FORMAT (/' State',I4,' is vector',I4,' Exp. Value:',2F14.8)
  100 FORMAT (2(2X,I4),2(3X,E13.6),2X,A8)
  200 FORMAT (/,20X,'PROPERTY =',F20.12,/)
      END
c
      SUBROUTINE GETIDS (I,J,K,L,ISYM,ND,NV,ILIND,IRIND,
     &                   JLIND,JRIND,JVECL,
     &                   IVWU,IVWR,IREC)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     GET 1-ELECTRON COEFFICIENTS :
C     CALCULATE IF IGENEX = 1
C     FROM FILE IF IGENEX = 2
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "mempoint.inc"
C
      DIMENSION ILIND(*),IRIND(*),JLIND(*),JRIND(*)
      DIMENSION JVECL(2,*)
      DIMENSION IVWU(2*MXREP,0:NORB,0:NELEC),IVWR(*)
C
      INTERN = NORBR(1) + NORBR(2)
      IF (MAXE3.GT.2) INTERN = NORB
      IF (IGENEX.EQ.1.OR.IREC.EQ.0) THEN
         CALL CALIDS_UP (I,J,K,L,ISYM,ND,NV,ILIND,IRIND,JLIND,
     &                   JRIND,JVECL,IVWU,IVWR(1),IVWR(1+IVWI))
         RETURN
      ELSEIF (IGENEX.EQ.2) THEN
C     ------------------------------------------
C     One cannot use the symmetry 
C
C     USE SYMMETRY : <I|E(i,j)|J> = <J|E(j,i)|I>
C
C     if pointgroup symmetry is used
C     ------------------------------------------
      READ (MRCEXCS,REC=IREC) INREC,NVERT,NDWR1,NDWR2,
     &     ((JVECL(IP,IV),IP=1,2),IV=1,NVERT),
     &     (JLIND(ID),JRIND(ID),ID=NDWR1,NDWR2)
  201 IF (INREC.EQ.0) RETURN
      IREC = INREC
      READ (MRCEXCS,REC=IREC) INREC,NVERT,NDWR1,NDWR2,
     &     (JLIND(ID),JRIND(ID),ID=NDWR1,NDWR2)
      GOTO 201
      ENDIF
      RETURN
      END
C
      SUBROUTINE GetVec (ivec,VecR,VecI)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
C
      Real*8 VecR(*),VecI(*),SigR(*),SigI(*)
C
      irec = ivec * 2 - 1
      IF (REALAR) THEN
         Read (MRCVECS,REC=IREC) (VECR(I),I=1,NDET)
      ELSE
         Read (MRCVECS,REC=IREC) (VECR(I),VECI(I),I=1,NDET)
      ENDIF
      Return
C
      ENTRY GetSig (ivec,Sigr,SigI)
      irec = ivec * 2 
      IF (REALAR) THEN
         Read (MRCVECS,REC=IREC) (SIGR(I),I=1,NDET)
      ELSE
         Read (MRCVECS,REC=IREC) (SIGR(I),SIGI(I),I=1,NDET)
      ENDIF
      Return
C
      END
C
      SUBROUTINE InitVec
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     Initialises temporary file for expansion vectors
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
C
      IF (REALAR) THEN
         IRECL = NDET * 8
      ELSE
         IRECL = NDET * 16
      ENDIF
      Open (MRCVECS,File=FNAM(7),ACCESS='DIRECT',RECL=IRECL,
     &      IOSTAT=IOS,ERR=100)
      Return
 100  CALL QUIT('Cannot open MRCVECS')
      End
C
      SUBROUTINE CINPRD (VECTR,VECTI,SIGMAR,SIGMAI)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =================================================================
C     COMPUTES INNER PRODUCT OF NEW SIGMA VECTORS WITH PREVIOUS VECTORS
C     =================================================================
C
#include "param.inc"
#include "general.inc"
#include "iterat.inc"
C
      DIMENSION VECTR(NDET,*),VECTI(NDET,*)
      DIMENSION SIGMAR(NDET,*),SIGMAI(NDET,*)
C
      KVEC = NROOTS + 1
C
C     Real or Complex arithmetic
C
      IF (REALAR) THEN
         DO 201 JVEC = 1, NVEC
         CALL GETVEC (JVEC,VECTR(1,KVEC),VECTI(1,KVEC))
         DO 101 IVEC = MAX0(1,JVEC-MVEC), NVEC - MVEC
         LVEC = IVEC + MVEC
C
         CMTR(LVEC,JVEC) = DDOT (NDET,SIGMAR(1,IVEC),1,VECTR(1,KVEC),1)
C64B     CMTR(LVEC,JVEC) = SDOT (NDET,SIGMAR(1,IVEC),1,VECTR(1,KVEC),1)
  101    CONTINUE
  201    CONTINUE
      ELSE
         DO 200 JVEC = 1, NVEC
         CALL GETVEC (JVEC,VECTR(1,KVEC),VECTI(1,KVEC))
         DO 100 IVEC = MAX0(1,JVEC-MVEC), NVEC - MVEC
         LVEC = IVEC + MVEC
C
         CMTR(LVEC,JVEC) = DDOT (NDET,SIGMAR(1,IVEC),1,VECTR(1,KVEC),1)
     &                   + DDOT (NDET,SIGMAI(1,IVEC),1,VECTI(1,KVEC),1)
         CMTI(LVEC,JVEC) = DDOT (NDET,SIGMAR(1,IVEC),1,VECTI(1,KVEC),1)
     &                   - DDOT (NDET,SIGMAI(1,IVEC),1,VECTR(1,KVEC),1)
C64B     CMTR(LVEC,JVEC) = SDOT (NDET,SIGMAR(1,IVEC),1,VECTR(1,KVEC),1)
C64B &                   + SDOT (NDET,SIGMAI(1,IVEC),1,VECTI(1,KVEC),1)
C64B     CMTI(LVEC,JVEC) = SDOT (NDET,SIGMAR(1,IVEC),1,VECTI(1,KVEC),1)
C64B &                   - SDOT (NDET,SIGMAI(1,IVEC),1,VECTR(1,KVEC),1)
  100    CONTINUE
  200    CONTINUE
      ENDIF
C
      RETURN
      END
C
      SUBROUTINE LEADING(IVWU)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
#include "files.inc"
#include "iterat.inc"
C
      PARAMETER (NDREC=1023)
C
      DIMENSION IVWU(2*MXREP,0:NORB,0:NELEC)
      DIMENSION IDREC(N8),IRAS1(N3),IRAS2E(N3),IRAS2H(N3),IRAS3(N3)
      DIMENSION COR(NDREC),COI(NDREC),IDETER(NDREC)
      LOGICAL GOVER
C
      WRITE(6,*)
      WRITE(6,*)
      WRITE(6,'(A35,F6.4)') 'Determinants with coefficient .GT. ',COMIN
#if defined (INT_STAR8)
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=8)
#else
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=4)
#endif
      READ (MRCFINV,REC=1,ERR=10000) LENREC
      CLOSE (MRCFINV)
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=LENREC)
      READ (MRCFINV,REC=1,ERR=10000) LENREC,NPDET,NPROOT,
     &                     (EIGEN1(IROOT),IROOT=1,NPROOT),
     &                     (IDREC(IROOT),IROOT=1,NPROOT)
      IF (NPDET.NE.NDET.OR.NPROOT.LT.NROOTS) GOTO 10000
      DO 7120 IROOT = 1,NROOTS
      WRITE(6,*)
      WRITE(6,'(A24,I4,A18,F16.8)') "Determinant(s) of root: ",
     $ IROOT," with eigenvalue: ",EIGEN1(IROOT)
      WRITE(6,*)
      IREC = IDREC(IROOT)
 7110 IF (REALAR) THEN
         READ (MRCFINV,REC=IREC,ERR=10000) INREC, NONZERO,
     &        (IDETER(IDNZ),IDNZ=1,NONZERO),
     &        (COR(IDNZ),IDNZ=1,NONZERO)
      ELSE
         READ (MRCFINV,REC=IREC,ERR=10000) INREC, NONZERO,
     &        (IDETER(IDNZ),IDNZ=1,NONZERO),
     &        (COR(IDNZ),COI(IDNZ),IDNZ=1,NONZERO)
      ENDIF
      DO 7111 IDNZ = 1, NONZERO
         COEFF = COR(IDNZ)**2
         IF (.NOT.REALAR) COEFF = COEFF + COI(IDNZ)**2
         IF (COEFF.GT.COMIN) THEN
            WRITE(6,'(A12,I8,A23,F16.8)') "Determinant ",IDETER(IDNZ),
     &       " with square of coeff. ",COEFF
            IA=0
            IB=0
            IC=0
            IG=0
            LORB=NORB
            IEL=NELEC
            IGAM=IREP
 7108       IF (LORB.EQ.0) GOTO 7109
            GOVER=IVWU(IGAM,LORB-1,IEL).GE.IDETER(IDNZ)
            IF (GOVER) THEN
               IF (LORB.LE.NORBR(1)) THEN
                  IA=IA+1
                  IRAS1(IA)=LORB
               ELSEIF (LORB.LE.(NORBR(1)+NORBR(2))) THEN
                  IB=IB+1
                  IRAS2H(IB)=LORB
               ENDIF
            ELSE
               IDETER(IDNZ)=IDETER(IDNZ)-IVWU(IGAM,LORB-1,IEL)
               IGAM=MULTBI(IRPAMO(LORB),IGAM)
               IEL=IEL-1
               IF (LORB.GT.(NORBR(1)+NORBR(2))) THEN
                  IG=IG+1
                  IRAS3(IG)=LORB
               ELSEIF (LORB.GT.NORBR(1)) THEN
                  IC=IC+1
                  IRAS2E(IC)=LORB
               ENDIF
            ENDIF
            LORB=LORB-1
            GOTO 7108
 7109       WRITE(6,12010) "Holes in RAS1: ",
     &           (IRAS1(I),REPNA(IRPAMO(IRAS1(I))),I=IA,1,-1)
            IF (IC.LE.IB) THEN
               WRITE(6,12010) "Elec. in RAS2: ",
     &              (IRAS2E(I),REPNA(IRPAMO(IRAS2E(I))),I=IC,1,-1)
            ELSE
               WRITE(6,12010) "Holes in RAS2: ",
     &              (IRAS2H(I),REPNA(IRPAMO(IRAS2H(I))),I=IB,1,-1)
            ENDIF
            WRITE(6,12010) "Elec. in RAS3: ",
     &           (IRAS3(I),REPNA(IRPAMO(IRAS3(I))),I=IG,1,-1)
         ENDIF
 7111 CONTINUE
      IF (INREC.NE.0) THEN
         IREC = INREC
         GOTO 7110
      ENDIF
 7120 CONTINUE
      CLOSE (MRCFINV,STATUS='KEEP')
C
12010 FORMAT(A15,6(I4,"(",A4,")"))
      RETURN
10000 CALL Q2JOB (2,'ERROR READING MRCFINV',' ',0)
      RETURN
      END
C
      SUBROUTINE MAKEGM(GMATR,GMATI,TINTR,TINTI,IJEX)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ================================================
C     READS 1- AND 2-ELECTRON INTEGRALS, MAKES GMATRIX
C     ================================================
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
#include "symm.inc"
C
      DIMENSION GMATR(*),GMATI(*),TINTR(*),TINTI(*),IJEX(*)
      real(8), allocatable :: PROP(:,:,:)
      real(8), allocatable :: AIJ(:)
      DIMENSION APHASE(2)
      LOGICAL ENDFIL
C
      OPEN (MRCONEE,FILE=FNAM(1),FORM='UNFORMATTED')
      REWIND (MRCONEE)
      READ (MRCONEE) morb
      READ (MRCONEE)
      READ (MRCONEE)
      READ (MRCONEE)
      READ (MRCONEE)
C     -------------------------------------------------------
C     Read one-electron integrals : use tint as scratch space
C     -------------------------------------------------------
      IF (REALAR) THEN
         READ (MRCONEE) (tintr(JI),tdum,JI=1,morb*morb)
      ELSE
         READ (MRCONEE) (tintr(JI),tinti(JI),JI=1,morb*morb)
      ENDIF
      CLOSE (MRCONEE)
C
C     -------------------------------------------------------------
C     Reorder one-electron integrals (some orbitals may be deleted)
C     -------------------------------------------------------------
      jid = 0
      do ii = 1, morb
	 do jj = 1, morb
	    jid = jid + 1
	    i = indmo(ii)
	    j = indmo(jj)
	    If (i.NE.0.AND.j.NE.0) Then
	       ji = norb*(i-1)+j
	       gmatr(ji) = tintr(jid)
!LV 28-9-2020: did not modify next line after MRCONEE change, code appeared to have already used Fortran ordering
	       IF (.NOT.REALAR) gmati(ji) = tinti(jid)
            EndIf
         enddo
      enddo
C
      REWIND (MRCTWOE)
    1 CALL TWOLIN (ENDFIL,TINTR,TINTI,IJEX,KL,NKL)
      IF (ENDFIL) GOTO 201
C
C     -------------------------------------
C     ADD CONTRIBUTION TO G-MATRIX ELEMENTS
C     -------------------------------------
C
      K = (KL-1)/NORB + 1
      L = MOD(KL-1,NORB) + 1
      DO 200 IX = 1, NKL
      IJ = IJEX(IX)
      I = (IJ-1)/NORB + 1
      J = MOD(IJ-1,NORB) + 1
      IF (K.EQ.J) THEN
         IL = (I-1) * NORB + L
         GMATR(IL) = GMATR(IL) - TINTR(IX)
         IF (.NOT.REALAR) GMATI(IL) = GMATI(IL) - TINTI(IX)
      ENDIF
  200 CONTINUE
      GOTO 1
  201 CONTINUE
CMKN
C     ******************************************
C     READ PROPERTY INTEGRALS FROM FILE 'MDPROP'
C     Proposed & Written by (Malaya K. Nayak)
C     ******************************************
      IF (PROPER) THEN
        allocate(PROP(2,MORB,MORB))
        allocate(AIJ(2*MORB*MORB))
        IF (REALAR) THEN
          IRW = 1
        ELSE
          IRW = 2
        END IF
        CALL PROP_INT(PROP,AIJ,MORB,IRW,APHASE)
C     *****************************************************
C     Reorder Property Integrals as per 'CI' Subroutine and
C     Replace the 1-electron integral by Property Integral
C     **********************************************
        IA = 0
        DO II=1,MORB
          DO JJ=1,MORB
            IA = IA+1
             I = INDMO(II)
             J = INDMO(JJ)
             IF (I.NE.0 .AND. J.NE.0) THEN
               IJ = (I-1)*NORB+J
               GMATR(IJ) = 0.0D0
               GMATI(IJ) = 0.0D0
               IF (REALAR) THEN
                 GMATR(IJ) = GMATR(IJ) + AIJ(IA)
               ELSE
                 GMATR(IJ) = GMATR(IJ) + AIJ(IRW*IA-1)
                 GMATI(IJ) = GMATI(IJ) + AIJ(IRW*IA)
              END IF
            END IF
          END DO
        END DO
        deallocate(PROP)
        deallocate(AIJ)
      END IF
C     *****************************************************
C     END OF READING & RE-ORDERING PROPERTY INTIGRALS
C     *****************************************************
CMKN
C
      RETURN
C
      END
C
      SUBROUTINE MAKIDS (ILIND,IRIND,IJREC,JLIND,JRIND,
     &                   JVECL,IVWA)
C
      IMPLICIT REAL*8 (A-H, O-Z)
      PARAMETER(DM=100.0D0)
C
C     GENERATE FILE WITH ALL 1-ELECTRON COUPLING COEFFICIENTS
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "mempoint.inc" 
C
      DIMENSION JVECL(2,*),JLIND(*),JRIND(*),ILIND(*),IRIND(*)
      DIMENSION IJREC(NORB,NORB),IVWA(*)
C
      NSTORE = 0
      NCOUPL = 0
      IREC = 0
      IVERT = 1
      ICOREX = 1
      LNCORE = LNORBC*NORB*NORB
      LEXREC = 512
C64B  LENREC = LEXREC * 8
      LENREC = LEXREC * 4
      IF (IGENEX.EQ.2) OPEN (MRCEXCS,FILE=FNAM(6),ACCESS='DIRECT',
     &                       RECL=LENREC)
      DO 200 I = 1, NORB
         DO 100 J = 1, NORB
            IJREC(I,J) = 0
            CALL GETIDS (I,J,0,0,IREP,ND,NVERT,ILIND,IRIND,
     &                   JLIND(ICOREX),JRIND(ICOREX),
     &                   JVECL(1,IVERT),IVWA(1),
     &                   IVWA(1+IVWI),0)
            NCOUPL = NCOUPL + ND
            IF (ND.LE.LNCORE) THEN
               NSTORE = NSTORE + ND
               JVERT = IVERT + NVERT
               IJREC(I,J) = IVERT
               DO 20 KVERT = IVERT, JVERT - 1
                  JVECL(1,KVERT) = JVECL(1,KVERT) + ICOREX - 1
   20             CONTINUE
               IVERT = JVERT
               ICOREX = ICOREX + ND
               LNCORE = LNCORE - ND
               GOTO 100
            ENDIF
            IF (IGENEX.EQ.1) GOTO 100
            IREC = IREC + 1
            IJREC(I,J) = - IREC
C           IJREC(J,I) = - IREC
            INREC = 0
            NSKIP = 2 * NVERT + 4
            NDWR1 = 1
            NDWR2 = MIN0(ND,NDWR1+(LEXREC-2-NSKIP)/3-1)
            IF (NDWR2.NE.ND) INREC = IREC + 1
            WRITE (MRCEXCS,REC=IREC) INREC,NVERT,NDWR1,NDWR2,
     &            ((JVECL(IP,IV),IP=1,2),IV=IVERT,IVERT+NVERT-1),
     &            (JLIND(ID),JRIND(ID),ID=NDWR1+ICOREX-1,
     &                                    NDWR2+ICOREX-1)
   99       IF (INREC.EQ.0) GOTO 100
            NSKIP = 4
            IREC = IREC + 1
            INREC = 0
            NDWR1 = NDWR2 + 1
            NDWR2 = MIN0(ND,NDWR1+(LEXREC-2-NSKIP)/3-1)
            IF (NDWR2.NE.ND) INREC = IREC + 1
            WRITE (MRCEXCS,REC=IREC) INREC,NVERT,NDWR1,NDWR2,
     &            (JLIND(ID),JRIND(ID),ID=NDWR1+ICOREX-1,
     &                                    NDWR2+ICOREX-1)
            GOTO 99
  100       CONTINUE
  200    CONTINUE
      IF (LNCORE.EQ.0) GOTO 500
      DO 300 I = 1, NORB
         DO 400 J = 1, NORB
            IF (IJREC(I,J).LE.0) THEN
               CALL GETIDS (I,J,0,0,IREP,ND,NVERT,ILIND,IRIND,
     &                      JLIND(ICOREX),JRIND(ICOREX),
     &                      JVECL(1,IVERT),IVWA(1),
     &                      IVWA(1+IVWI),0)
               IF (ND.LE.LNCORE) THEN
                  NSTORE = NSTORE + ND
                  NCOUPL = NCOUPL - ND
                  JVERT = IVERT + NVERT
                  IJREC(I,J) = IVERT
                  DO 30 KVERT = IVERT, JVERT - 1
                     JVECL(1,KVERT) = JVECL(1,KVERT) + ICOREX - 1
   30             CONTINUE
                  IVERT = JVERT
                  ICOREX = ICOREX + ND
                  LNCORE = LNCORE - ND
                  IF (LNCORE.EQ.0) GOTO 500
                  GOTO 400
               ENDIF
            ENDIF
  400    CONTINUE
  300 CONTINUE
  500 NCOREX = ICOREX
      NTVERT = IVERT
      RSTORE = REAL(NSTORE) / REAL(NCOUPL) * DM
      WRITE (6,1000) RSTORE
      RETURN
C
 1000 FORMAT (/,F6.2,' Percent of the Coupling ',
     & 'Coefficients are stored.')
C
      END
      SUBROUTINE NATORB(IROOT,DENSR,DENSI,WORK1,WORK2,WORK3)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     =====================================
C     CALCULATES NATURAL ORBITALS (SPINORS)
C     =====================================
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
C
      DIMENSION DENSR(NORB,NORB),DENSI(NORB,NORB)
      DIMENSION WORK1(*),WORK2(*),WORK3(*)
      Dimension NSinRep(16),Ndeg(16)
      Character*8 FilNNS,FilVCA,FilNSc
C
CLV   NATORB does not work with DIRAC, for use with MOLFDIR : remove
CLV   the RETURN statement of the next line.
#ifdef PRG_DIRAC
      WRITE(6,'(/A/A)')
     &  ' DIRRCI INFO: Generation of natural spinors requested, but',
     &  ' generation of natural spinors is not implemented in DIRAC.'
      RETURN
#else
      Write (6,1000)
      FilVCA = 'MFDVECA'
      FilNSc = 'SCRATCH'
      Open (20,File=FilVCA,Form='Formatted',STATUS='OLD',ERR=110)
      Open (22,File=FilNSc,Form='Formatted')
      MaxSinRep = 0
      Do I = 1, 16
	 NSinRep(I) = 0
	 Ndeg(I) = 0
      EndDo
      Do I = 1,16
        DO m=1,16
          IF (RepN(I)(1:4).EQ.RepN(m)(1:4)) Ndeg(I)=Ndeg(I)+1
	ENDDO
      ENDDO
      Do I = 1, Norb
	 NSinRep(IrpMO(I)) = NSinRep(IrpMO(I)) + 1
      EndDo
      Do I = 1, 16
	 If (MaxSinRep.LT.NSinRep(I)) MaxSinRep = NSinRep(I)
      EndDo
C
      Ioff1 = 1			     
      Ioff5 = Ioff1 + MaxSinRep
      Ioff2 = Ioff5 + MaxSinRep
      Ioff3 = Ioff2 + MaxSinRep
      Ioff4 = Ioff3 + MaxSinRep
      IF (REALAR) THEN
         ISize = Ioff4 + MaxSinRep**2
      ELSE
         ISize = Ioff4 + MaxSinRep
      ENDIF
      IF (ISize.GT.2*NORBD) THEN
         WRITE(6,10000) 
         Close (22,Status='Delete')
         RETURN
      ENDIF
      Ioff6 = Ioff5 + MaxSinRep   	
      Ioff7 = Ioff6 + MaxSinRep
      ISize = Ioff7 + 4*MaxSinRep**2
      IF (ISize.GT.2*NORBD) THEN
         WRITE(6,10000) 
         Close (22,Status='Delete')
         RETURN
      ENDIF
C
      Write (FilNNS,'(A6,I2)') 'MFDNAT',IROOT
      IF (Iroot.lt.10) FilNNS(7:7) = '_'
      Rewind (20)
      Open (21,File=FilNNS,Form='Formatted')
      Call CopyFF (20,21)
      Do 100 Irp = 1, NSymrp
         NS = NSinRep(Irp)
C        Skip degenerate representations.
	 If (Irp.GT.1) Then
	    If (RepN(Irp)(1:4).EQ.RepN(Irp-1)(1:4)) Goto 100
	 EndIf
C        Skip representations without active spinors
	 If (NSinRep(Irp).eq.0) GoTo 100
C        Average with degenerate representations 
	 Do nd=0,Ndeg(Irp)-1
	    k = 0
            iorb = 0
            Do i =1, Norb
               If (Irpmo(i).EQ.Irp+nd) Then
	         Work3(Ioff5+iorb) = Eps(i)
	         iorb = iorb + 1
	         jorb = 0 
	         Do j = 1, Norb
	           If (IrpMO(j).EQ.Irp+nd) Then
		     If (nd.EQ.0) Then
		       Work1(k+1) = 0.0D0
		       IF (.NOT.REALAR) Work2(k+1) = 0.0D0
                     Endif
	             jorb = jorb + 1
C----------------------------------------------------
C Use - (Density matrix) to get eigenvalues in right
C order from the diagonalisation routine
C----------------------------------------------------
	             Work1(k+1) = Work1(k+1) - Densr(j,i)*
     &                            (1.0D0/Ndeg(Irp))
                     IF (.NOT.REALAR) Work2(k+1) = Work2(k+1) -
     &                                Densi(j,i)*(1.0D0/Ndeg(Irp))
	             k = k + 1
                   EndIf
	         EndDo
	       EndIf
            EndDo
         Enddo
         IFAIL=0
         IF (REALAR) THEN
            CALL RS (NS,NS,Work1,Work3(Ioff1),1,Work3(Ioff4),
     &               Work3(Ioff2),Work3(Ioff3),IFAIL)
            DO I = 1, NS*NS
               Work1(I) = Work3(Ioff4+I-1)
            ENDDO
         ELSE
Clv         Use NAG if the library is available
Clv         Better : rewrite as call to CH
            IFAIL = 1
C           TOL = 1.0D-14
C           CALL F01BCF(NS,TOL,Work1,NS,Work2,NS, Work3(Ioff1),
C64B        CALL F01BCE(NS,TOL,Work1,NS,Work2,NS, Work3(Ioff1),
C    &                  Work3(Ioff2),Work3(Ioff3),Work3(Ioff4))
C           EPSI = X02AJF()
C           CALL F02AYF(NS,EPSI,Work3(Ioff1),Work3(Ioff2),
C64B        CALL F02AYE(NS,EPSI,Work3(Ioff1),Work3(Ioff2),
C    &                  Work1,NS,Work2,NS,IFail)
         ENDIF
         IF (IFAIL.NE.0) CALL QUIT('**ERROR-CANNOT DIAGONALISE CMTR**')
         Do IEIGEN = 1, NS
            Work3(Ioff1+IEigen-1) = - Work3(Ioff1+IEigen-1)
         EndDo
         Write (6,1010) Iroot, RepN(Irp)
         Write (6,1020) (Work3(Ioff1+IEigen-1),IEigen=1,NS)
         Call NSRepl (Irp,NS,Work3(Ioff1),Work1,Work2,RepN(Irp),
     &                Work3(Ioff5),Work3(Ioff6),Work3(Ioff7))
  100 Continue 
      Close (21)
      Write (6,1030) Iroot, FilNNS
      Close (22,Status='Delete')
  110 Continue
      RETURN
 1000 Format (//' Diagonalising (active part) of Density Matrix')
 1010 Format (/' Eigenvalues Density Matrix (Root',I3,', Repr. ',
     &         A14,')')
 1020 Format (6F12.8)
 1030 Format (/' Natural Spinors for Root',I3,' written to ',A8)
10000 Format (/' Cannot make the natural spinors due to ',
     &         'lack of memory')
#endif
      END
C
      INTEGER FUNCTION NOVERI (N,I)
C
      II = MIN0(I,N-I)
      NI = 1
      DO 1 J = 0, II - 1
    1 NI = (N-J) * NI
      DO 2 J = 1, II
    2 NI = NI / J
      NOVERI = NI
C
      RETURN
      END
C        
      SUBROUTINE NSRepl (Irp,NS,OccNr,VecR,VecI,
     &                   RepName,EigenV,IRecSP,WORK)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C    ===========================================
C     Replaces Active Spinors by Natural Spinors
C    ===========================================
C
#include "param.inc"
#include "symm.inc"
#include "general.inc"
C
      Dimension OccNr(NS),VecR(NS,NS),VecI(NS,NS),EigenV(NS)
      Dimension IRecSp(NS),WORK(*)
      Character Card*132,Format*20,RepName*14
C
      Do IS = 1, NS
	 IRecSp(IS) = 0
      EndDo
C     First Pass : Find the Active Spinors.
      Rewind (21)
      Read (21,'(A)') Card
      Read (21,'(A)') Card
      Read (21,'(A)') Card
      Read (21,'(A)') Card
      Read (21,'(A)') Card
      Read (Card,'(A)') Format
      Read (21,'(A)') Card
      Read (21,'(A)') Card
      Read (21,'(A)') Card
      Irec = 8
      IFound = - 1
   10 Read (21,'(A)',End=30) Card
      Irec = Irec + 1
      If (Card(6:10).EQ.'*****'.AND.Card(12:25).EQ.RepName) Then
         Read (Card(31:38),'(2I4)') NST,Nbas
	 IFound = 0
C
C        Allocate pointers to store two sets of NS vectors of Nbas length
C
         Ioff1 = 1
         Ioff2 = Ioff1 + Nbas * NS
	 Ioff3 = Ioff2 + Nbas * NS
	 Ioff4 = Ioff3 + Nbas * NS
      EndIf
      If (Card(2:3).EQ.'MO'.AND.IFound.GE.0) Then
         Read (Card,'(10X,F20.10)') EV
C        Check whether this is one of the active set.
C        If the representation conatins degenerate functions we make
C        sure that this one is not appointed to a previously found spinor.
         Do IS = 1, NS
            If (EV.EQ.EigenV(IS).AND.IrecSp(IS).EQ.0) Then
               IrecSp(IS) = Irec
	       IFound = Ifound + 1
	       GoTo 20
	     EndIf
          EndDo
   20     If (Ifound.EQ.NS) GoTo 30
       EndIf
       GoTo 10
C
   30 If (Ifound.NE.NS) Then
	 Print*, 'Cannot make natural spinors for repr. ',Repname
	 Print*, 'Make sure that original MFDVECA file is present'
	 Return
      Endif
C
C     Read in Canonical Spinors
C
      Do IS = 1, NS
         Rewind (21)
         Do IR = 1, IrecSp(IS)
            Read (21,'(A)') Card
         EndDo
	 Ioff1a = Ioff1 + (IS-1) * Nbas
	 Ioff2a = Ioff2 + (IS-1) * Nbas
         Read (21,Format) (Work(i+ioff1a-1),Work(i+ioff2a-1),i=1,nbas)
      EndDo
C
C     Make Natural Spinors
C
      Call ZZero (Nbas*Ns,Work(Ioff3))
      Call ZZero (Nbas*Ns,Work(Ioff4))
      Do IS = 1, NS
	 Ioff1a = Ioff1 + (IS-1) * Nbas
	 Ioff2a = Ioff2 + (IS-1) * Nbas
	 Do JS = 1, NS
	    Ioff3a = Ioff3 + (JS-1) * Nbas
	    Ioff4a = Ioff4 + (JS-1) * Nbas
	    Do I = 1, Nbas
	       Work(Ioff3a+i-1) = Work(Ioff3a+i-1)
     &                          + Work(Ioff1a+i-1) * VecR(IS,JS)
               IF (.NOT.REALAR) THEN
               Work(Ioff3a+i-1) = Work(Ioff3a+i-1) 
     &                          - Work(Ioff2a+i-1) * VecI(IS,JS)
	       Work(Ioff4a+i-1) = Work(Ioff4a+i-1)
     &                          + Work(Ioff1a+i-1) * VecI(IS,JS)
     &                          + Work(Ioff2a+i-1) * VecR(IS,JS)
               ENDIF
	     EndDo
          EndDo
      EndDo
C
C     Write the Natural Spinors.
C
      Do IS = 1, NS
         Rewind (21)
         Rewind (22)
         Do IR = 1, IrecSp(IS)-1
            Read (21,'(A)') Card
            Write (22,'(A)') Card
         EndDo
	 Read (21,'(A)') Card
	 Read (Card,'(4X,I2)') MO
	 Write(Card,1000) MO,IS,OccNr(IS),0.0
	 Write (22,'(A)') Card
	 Ioff3a = Ioff3 + (IS-1) * Nbas
	 Ioff4a = Ioff4 + (IS-1) * Nbas
         Write (22,Format) (Work(i+ioff3a-1),Work(i+ioff4a-1),i=1,nbas)
C
C        Skip just as many records on unit 21
C
         Read (21,Format) (Work(ioff1+i-1),i=1,2*nbas)
         Call CopyFF (21,22)
	 Rewind (21)
	 Rewind (22)
         Call CopyFF (22,21)
      EndDo
      Return
 1000 FORMAT(1X,'MO:',I2,' NS:',I2,18X,' OCC:',F15.10,' COPCOF:',
     +       G20.10)
      End
C
      SUBROUTINE SYMMIN
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "general.inc"
#include "iterat.inc"
C
    1 OPEN (MRCONEE,FILE=FNAM(1),FORM='UNFORMATTED')
      READ (MRCONEE) NORB,BREIT,ECORE,NFSYM,NZ,SPINFR
      REALAR = NZ.EQ.1.OR.SPINFR
      IF (NORB.GT.N3) CALL Q2JOB(3,'P.ORBSIN','    N3',NORB)
      READ (MRCONEE) NSYMRP,(REPN(IRP),IRP=1,NSYMRP)
      READ (MRCONEE) NREP,(REPNA(IRP),IRP=1,2*NREP)
      READ (MRCONEE) ((MULTB(I,J),I=1,2*NREP),J=1,2*NREP)
      READ (MRCONEE) (IRPMO(IMO),IRPAMO(IMO),EPS(IMO),IMO=1,NORB)
      CLOSE (MRCONEE)
      NREP = 2 * NREP
C
      WRITE(6,1200)
      WRITE(6,1220) BREIT,ECORE,NORB
C
C     ------------------------------------
C     MAKE INVERSE OF MULTIPLICATION TABLE
C     ------------------------------------
      DO 500 IR = 1, NREP
      DO 400 JR = 1, NREP
      KR = MULTB(IR,JR)
      MULTBI(IR,KR) = JR
  400 CONTINUE
  500 CONTINUE
C     ---------------------------------------
C     DETERMINE POSITION OF A1 REPRESENTATION
C     ---------------------------------------
      DO 40 IR = 1, NREP
      DO 39 JR = 1, NREP
      IF (MULTB(IR,JR).NE.JR) GOTO 40
   39 CONTINUE
      IRA1 = IR
      GOTO 41
   40 CONTINUE
      CALL Q2JOB (6,'A1 REPRESENTATION NOT FOUND',' ',0)
   41 CONTINUE
C
C
      RETURN
 1200 FORMAT(//1X,'Information read from MRCONEE ')
 1220 FORMAT(//1X,'Breit Interaction included :',T50,L1
     +        /1X,'Core Energy :',T30,G20.10
     +        /1X,'Number of orbitals :',T47,I4)
      END
C
      SUBROUTINE ORBSIN
C
C     Reorder the orbitals based on the RAS division
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "symm.inc"
#include "general.inc"
#include "iterat.inc"
C
      DO I = 1, NORB
         INDMO(I) = 0
      ENDDO
C
      WRITE (6,1200)
C
      MORB = 0
      DO IRAS = 1, 3
         DO IRP = 1, NSYMRP
            NIRP = 0
            DO I=1,NORB
               IF (IRPMO(I).EQ.IRP) THEN
                  NIRP = NIRP + 1
                  IF (NIRP.LE.NRAS(1,IRP).AND.IRAS.EQ.1) THEN
                     MORB = MORB + 1
                     INDMO(I) = MORB
                     WRITE(6,1230) I,MORB,REPN(IRPMO(I)),EPS(I)
                  ELSEIF (NIRP.GT.NRAS(1,IRP).AND.NIRP.LE.NRAS(2,IRP)
     &                                       .AND.IRAS.EQ.2) THEN
                     MORB = MORB + 1
                     INDMO(I) = MORB
                     WRITE(6,1230) I,MORB,REPN(IRPMO(I)),EPS(I)
                  ELSEIF (NIRP.GT.NRAS(2,IRP).AND.IRAS.EQ.3) THEN
                     MORB = MORB + 1
                     INDMO(I) = MORB
                     WRITE(6,1230) I,MORB,REPN(IRPMO(I)),EPS(I)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
      DO I = 1, NORB
         INDSRT(INDMO(I))= IRPMO(I)
      ENDDO
      DO I = 1, NORB
         IRPMO(I) = INDSRT(I)
      ENDDO
C
      DO I = 1, NORB
         INDSRT(INDMO(I))= IRPAMO(I)
      ENDDO
      DO I = 1, NORB
         IRPAMO(I) = INDSRT(I)
      ENDDO
C
      DO I = 1, NORB
         EPSSRT(INDMO(I))= EPS(I)
      ENDDO
      DO I = 1, NORB
         EPS(I) = EPSSRT(I)
      ENDDO
C
 300  CONTINUE
 310  CONTINUE
      RETURN
 1200 FORMAT (//1X,'Orbital',T10,'Index',
     +        T20,'Representation',T42,'Orbital Energy')
 1230 FORMAT(2X,I4,T10,1X,I4,T20,A14,T40,G20.10,5X,A)
      END
C
      SUBROUTINE ORTHN (VECTR,VECTI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     =========================================================
C     SCHMIDT ORTHONORMALIZATION OF THE NEW EXPANSION VECTORS
C     ON THE PREVIOUS ONES (2 TIMES TO AVOID NUMERICAL PROBLEMS)
C     ==========================================================
C
#include "param.inc"
#include "general.inc"
C
      DIMENSION VECTR(NDET,*),VECTI(NDET,*)
C
      KVEC = NROOTS + 1
      DO 1000 IVEC = 1, NVEC - MVEC
C
      RR = DDOT (NDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
C64B  RR = SDOT (NDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) RR = RR +
     &     DDOT (NDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
C64B &     SDOT (NDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
      RR = 1.D0 / DSQRT (RR)
C64B  RR = 1.D0 / SQRT (RR)
      CALL DSCAL (NDET,RR,VECTR(1,IVEC),1)
C64B  CALL SSCAL (NDET,RR,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) 
     &   CALL DSCAL (NDET,RR,VECTI(1,IVEC),1)
C64B &   CALL SSCAL (NDET,RR,VECTI(1,IVEC),1)
C
      DO 200 IONR = 1, 2
      DO 100 JVEC = 1, MVEC + IVEC - 1
      CALL GETVEC (JVEC,VECTR(1,KVEC),VECTI(1,KVEC))
      RR = DDOT (NDET,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
C64B  RR = SDOT (NDET,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) THEN
         RI = DDOT (NDET,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
         RI = RI - DDOT (NDET,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
         RR = RR + DDOT (NDET,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
C
C64B     RI = SDOT (NDET,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
C64B     RI = RI - SDOT (NDET,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
C64B     RR = RR + SDOT (NDET,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
      ENDIF
      CALL DAXPY (NDET,-RR,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
C64B  CALL SAXPY (NDET,-RR,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) THEN
         CALL DAXPY (NDET,-RI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
         CALL DAXPY (NDET,-RR,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
         CALL DAXPY (NDET, RI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
C
C64B     CALL SAXPY (NDET,-RI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,-RR,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET, RI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
      ENDIF
  100 CONTINUE
      RR = DDOT (NDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
C64B  RR = SDOT (NDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) RR = RR +
     &     DDOT (NDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
C64B &     SDOT (NDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
      RR = 1.D0 / DSQRT (RR)
C64B  RR = 1.D0 / SQRT (RR)
      CALL DSCAL (NDET,RR,VECTR(1,IVEC),1)
C64B  CALL SSCAL (NDET,RR,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) 
     &   CALL DSCAL (NDET,RR,VECTI(1,IVEC),1)
C64B &   CALL SSCAL (NDET,RR,VECTI(1,IVEC),1)
  200 CONTINUE
      IF (RR.GT.1.D6) THEN
         WRITE (*,*) ' NORM OF NEW VECTOR < 1.D-12 '
         WRITE (*,*)
         WRITE (*,*) '  -- SILLY END OF PROGRAM -- '
         CALL QUIT('dirrci line  2473')
      ENDIF
      CALL PUTVEC (MVEC+IVEC,VECTR(1,IVEC),VECTI(1,IVEC))
 1000 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE PRECOND (DGEL,VECTR,VECTI)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ==================================
C     PRECONDITIONING OF RESIDUAL VECTOR
C     ==================================
C
#include "param.inc"
#include "general.inc"
#include "iterat.inc"
C
      DIMENSION VECTR(NDET,*),VECTI(NDET,*)
      DIMENSION DGEL(NDET)
C
      DO 200 IVEC = 1, NVEC - MVEC
      JROOT = IVEC
      E = EIGEN(NSEL(JROOT))
C
      DO 100 I = 1, NDET
      DE = DGEL(I) - E
C64B  IF (ABS(DE).LT.1.E-8) DE = SIGN(1.E-8,DE)
      IF (ABS(DE).LT.1.D-8) DE = SIGN(1.D-8,DE)
      VECTR(I,IVEC) = VECTR(I,IVEC) / DE
      IF (.NOT.REALAR)
     &   VECTI(I,IVEC) = VECTI(I,IVEC) / DE
  100 CONTINUE
  200 CONTINUE
      RETURN
      END
C
      SUBROUTINE PutVec (ivec,VecR,VecI)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "general.inc"
#include "files.inc"
C
      Real*8 VecR(*),VecI(*),SigR(*),SigI(*)
C
      irec = ivec * 2 - 1
      IF (REALAR) THEN
         Write (MRCVECS,REC=IREC) (VECR(I),I=1,NDET)
      ELSE
         Write (MRCVECS,REC=IREC) (VECR(I),VECI(I),I=1,NDET)
      ENDIF
      Return
C
      ENTRY PutSig (ivec,Sigr,SigI)
      irec = ivec * 2 
      IF (REALAR) THEN
         Write (MRCVECS,REC=IREC) (SIGR(I),I=1,NDET)
      ELSE
         Write (MRCVECS,REC=IREC) (SIGR(I),SIGI(I),I=1,NDET)
      ENDIF
      Return
C
      END
C
      SUBROUTINE REAVEC (VECTR,VECTI,DR,DI,ILIND)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
#include "iterat.inc"
C
      DIMENSION VECTR(NDET,*),VECTI(NDET,*),DR(*),DI(*)
      DIMENSION ILIND(*)
      DIMENSION IDREC(N8)
C
#if defined (INT_STAR8)
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=8)
#else
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=4)
#endif
      READ (MRCFINV,REC=1,ERR=10000) LENREC
      CLOSE (MRCFINV)
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=LENREC)
      READ (MRCFINV,REC=1,ERR=10000) LENREC,NPDET,NPROOT,
     &                     (EIGEN1(IROOT),IROOT=1,NPROOT),
     &                     (IDREC(IROOT),IROOT=1,NPROOT)
      IF (NPDET.NE.NDET.OR.NPROOT.LT.NROOTS) GOTO 10000
      DO 20 IROOT = 1, NROOTS
      EIGEN2(IROOT) = EIGEN1(IROOT)
      IF (.NOT.PROPER) WRITE (6,1000) IROOT,EIGEN1(IROOT)
      IVEC = IROOT
      CALL ZZERO (NDET,VECTR(1,IVEC))
      IF (.NOT.REALAR) CALL ZZERO (NDET,VECTI(1,IVEC))
      IREC = IDREC(IROOT)
  10  IF (REALAR) THEN
        READ (MRCFINV,REC=IREC,ERR=10000) INREC, NONZERO,
     &       (ILIND(IDNZ),IDNZ=1,NONZERO),
     &       (DR(IDNZ),IDNZ=1,NONZERO)
      ELSE
        READ (MRCFINV,REC=IREC,ERR=10000) INREC, NONZERO,
     &       (ILIND(IDNZ),IDNZ=1,NONZERO),
     &       (DR(IDNZ),DI(IDNZ),IDNZ=1,NONZERO)
      ENDIF
      DO 11 IDNZ = 1, NONZERO
         VECTR(ILIND(IDNZ),IVEC) = DR(IDNZ)
  11     CONTINUE
      IF (.NOT.REALAR) THEN
         DO 12 IDNZ = 1, NONZERO
            VECTI(ILIND(IDNZ),IVEC) = DI(IDNZ)
  12     CONTINUE
      ENDIF
      IF (INREC.NE.0) THEN
         IREC = INREC
         GOTO 10
      ENDIF
      NVEC = IVEC
      CALL ORTHN (VECTR,VECTI)
  20  CONTINUE
      CLOSE (MRCFINV,STATUS='KEEP')
C
      RETURN
 1000 FORMAT (/' Vector #',I4,' read from MRCFINV; Eigenvalue ',F20.12)
10000 CALL Q2JOB (2,'ERROR READING MRCFINV',' ',0)
      RETURN
      END
C
      SUBROUTINE RESIDUE (VECTR,VECTI,SIGMAR,SIGMAI)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "general.inc"
#include "iterat.inc"
C
      DIMENSION VECTR(NDET,*),VECTI(NDET,*)
      DIMENSION SIGMAR(NDET,*),SIGMAI(NDET,*)
C
      DO IROOT = 1, NROOTS
      IVEC = IROOT
      CALL ZZERO (NDET,VECTR(1,IVEC))
      IF (.NOT.REALAR) CALL ZZERO (NDET,VECTI(1,IVEC))
      ENDDO
C
C     Calculate residue : R = HC - EC = Sig.B - E.D.B
C     Start with EC part.
      KVEC = NROOTS + 1
      DO 20 N = 1, NVEC
      CALL GETVEC (N,VECTR(1,KVEC),VECTI(1,KVEC))
      DO 10 IROOT = 1, NROOTS
      IVEC = IROOT
      ER = - EIGEN(NSEL(IROOT)) * EVECR(N,NSEL(IROOT))
      CALL DAXPY (NDET,ER,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
C64B  CALL SAXPY (NDET,ER,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) THEN
         EI = - EIGEN(NSEL(IROOT)) * EVECI(N,NSEL(IROOT))
         CALL DAXPY (NDET,EI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,EI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
         CALL DAXPY (NDET,ER,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,ER,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
         EI = - EI
         CALL DAXPY (NDET,EI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
C64B     CALL SAXPY (NDET,EI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
      ENDIF
   10 CONTINUE
   20 CONTINUE
C     Calculate HC part.
      DO 40 N = 1, NVEC
      CALL GETSIG (N,VECTR(1,KVEC),VECTI(1,KVEC))
      DO 30 IROOT = 1, NROOTS
      IVEC = IROOT
      ER = EVECR(N,NSEL(IROOT))
      CALL DAXPY (NDET,ER,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
C64B  CALL SAXPY (NDET,ER,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) THEN
         EI = EVECI(N,NSEL(IROOT))
         CALL DAXPY (NDET,ER,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,ER,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
         CALL DAXPY (NDET,EI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,EI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
         EI = - EI
         CALL DAXPY (NDET,EI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
C64B     CALL SAXPY (NDET,EI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
      ENDIF
   30 CONTINUE
   40 CONTINUE
C
C     Calculate norm of the residual vector
C
      DO 50 IROOT = 1, NROOTS
      IVEC = IROOT
      RES = DDOT (NDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
C64B  RES = SDOT (NDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) RES = RES +
     &      DDOT (NDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
C64B &      SDOT (NDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
      RESIDU(IROOT) = RES
   50 CONTINUE
C
C     Calculate norm of reference part of the residual vector
C
      DO 60 IROOT = 1, NROOTS
      IVEC = IROOT
      RES = DDOT (NREFDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
C64B  RES = SDOT (NREFDET,VECTR(1,IVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) RES = RES +
     &      DDOT (NREFDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
C64B &      SDOT (NREFDET,VECTI(1,IVEC),1,VECTI(1,IVEC),1)
      REFRDU(IROOT) = RES
   60 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SELECV
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ====================
C     SELECTS EIGENVECTORS
C     ====================
C
#include "param.inc"
#include "general.inc"
#include "iterat.inc"
C
      DO 20 IROOT = 1, NROOTS
         IF (SELECT) THEN
CTROND         IF (SELECT.OR.ITER.GT.2) THEN
         OVMAX = 0.D0
         DO 10 I = 1, NVEC
            OV = EVECR(IROOT,I)*EVECR(IROOT,I)
            IF (.NOT.REALAR)
     &      OV = OV + EVECI(IROOT,I)*EVECI(IROOT,I)
            IF (OV.GT.OVMAX) THEN
               NSEL(IROOT) = I
               OVMAX = OV
            ENDIF
   10       CONTINUE
         ELSE
            NSEL(IROOT) = IROOT
         ENDIF
   20 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE SHUFFLE (RKLR,RKLI,INDK,INDL,TIR,TII)
C
C     Modified to read in only real integrals when applicable.
C     Note that this is not done efficiently as the zero imaginary
C     parts are still processed. This will have to wait until a more 
C     general integral face is developed.
C     LV : 28-6-2002
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
#include "symm.inc"
#include "mempoint.inc"
C
C     ---------------------------------------------------------
C     CANONICAL INTEGRALS ARE GENERATED FROM THE ORIGINAL LIST.
C     DIAGONAL INTEGRALS (IJ=KL) NEED TO BE DIVIDED BY 2.
C     shift if necessary if spinors are deleted
C     ---------------------------------------------------------
C
      REAL*8 RKLR(*),RKLI(*),TIR(*),TII(*)
      INTEGER INDK(*),INDL(*)
      DIMENSION KR(-N3:N3)
C
      CHARACTER*10 DATEX,DATEXB,TIMEX*8,TIMEXB*8
      CHARACTER*10 NEXTFIL
      LOGICAL CHECK_INTS

      CHECK_INTS = .FALSE.
CMI   CHECK_INTS = .TRUE.
      TRESINT = 1.D-12
      NINT = 0
      KLS = 0
      KLF = M2
      MM = M1 * M2
      OPEN(MDCINT,FILE=FNAM(2),FORM='UNFORMATTED')
      Read (mdcint,err=10000,end=10000) datex,timex,nkr,
     & (kr(i),kr(-i),i=1,nkr)

      IF (CHECK_INTS) THEN
        write(6,*) 'MDCINT control printout; datex=',datex,
     &  ' timex=',timex,' nkr=',nkr
        do i = 1, nkr
         write(6,'(2x,a,i3,a,i3,a,i3)')
     &   'i=',i,' kr(i)=',kr(i),' kr(-i)=',kr(-i)
        enddo
      ENDIF

      IF (2*NKR.NE.NORB) 
     &   CALL QUIT('MRCONEE AND MDCINT ARE INCOMPATIBLE')
      IF (BREIT) THEN
      OPEN(MDBINT,FILE=FNAM(3),FORM='UNFORMATTED')
      Read (mdbint,err=10001,end=10001) datexb,timexb,nkrb
      IF (NKR.NE.NKRB) CALL QUIT('MDCINT AND MDBINT ARE INCOMPATIBLE')
      ENDIF
      IF (REALAR) CALL DCOPY(M1,0.0,0,RKLI,1)
      OPEN (MRCTWOE,FILE=FNAM(8),FORM='UNFORMATTED')
   10 REWIND (MDCINT)
      READ (MDCINT)
      IF (BREIT) REWIND (MDBINT)
      IF (BREIT) READ (MDBINT)
      MINT = 0
      MINTB = 0
      MAX = MM
      IF (BREIT) MAX = MAX * 2
      CALL DCOPY(MAX,0.D0,0,TIR,1)
      CALL DCOPY(MAX,0.D0,0,TII,1)
C
C     -------------------------
C     READ IN COULOMB INTEGRALS
C     -------------------------
C
  30  CONTINUE
      IF (REALAR) THEN
         read (mdcint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (rklr(inz),inz=1,nz)
      ELSE
         read (mdcint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (rklr(inz),rkli(inz),inz=1,nz)
      ENDIF

      IF (CHECK_INTS) THEN
       write(6,*) 'read ikr,jkr,nz:',ikr,jkr,nz
      ENDIF
C
      MINT = MINT + NZ
      IF (IKR.EQ.0) THEN
         IF (JKR.EQ.0) THEN
            NZ = -1
         ELSE
            READ (MDCINT) NEXTFIL
            CLOSE (MDCINT,STATUS='KEEP')
            OPEN(MDCINT,FILE=NEXTFIL,FORM='UNFORMATTED')
            GOTO 30
         ENDIF
      ENDIF

      IF (CHECK_INTS) THEN
        write(6,'(2x,a,3i4)') 'read IKR,JKR,NONZR:',IKR,JKR,NZ
        do i = 1, nz
          write(6,*) indk(i),indl(i),rklr(i)
        enddo
      ENDIF

      If (nz.EQ.-1) GoTo 40 ! End of File
      i = indmo(kr(ikr))
      j = indmo(kr(jkr))
      itr = indmo(kr(-ikr))
      jtr = indmo(kr(-jkr))
      SignIJ = SIGN(1,ikr) * SIGN(1,jkr)
      DO INZ = 1, NZ
         kkr = indk(inz)
         lkr = indl(inz)
         k = indmo(kr(kkr))
         l = indmo(kr(lkr))
         ktr = indmo(kr(-kkr))
         ltr = indmo(kr(-lkr))
         SignKL = SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral
         IF (I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
         TWOR = RKLR(INZ)
         TWOI = RKLI(INZ)
  
          IF (CHECK_INTS) THEN
           write(6,'(a,4i4,a,f12.6)') 'i,j,k,l=',i,j,k,l,
     &     ' 2-el. integral:',TWOR
          ENDIF

         IJ = (I-1) * NORB + J
         KL = (K-1) * NORB + L
         JI = (J-1) * NORB + I
         LK = (L-1) * NORB + K
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1) = TWOR
            TII(IJ+(KL-KLS-1)*M1) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1) = TWOR
            TII(JI+(LK-KLS-1)*M1) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1) = TWOR
            TII(KL+(IJ-KLS-1)*M1) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1) = TWOR
            TII(LK+(JI-KLS-1)*M1) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
C--> Time reversal left
         IF (ITR.NE.0.AND.JTR.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
         TWOR = RKLR(INZ) * SignIJ
         TWOI = RKLI(INZ) * SignIJ
         IJ = (JTR-1) * NORB + ITR
         KL = (K-1) * NORB + L
         JI = (ITR-1) * NORB + JTR
         LK = (L-1) * NORB + K
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1) = TWOR
            TII(IJ+(KL-KLS-1)*M1) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1) = TWOR
            TII(JI+(LK-KLS-1)*M1) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1) = TWOR
            TII(KL+(IJ-KLS-1)*M1) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1) = TWOR
            TII(LK+(JI-KLS-1)*M1) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
C--> Time reversal right
         IF (I.NE.0.AND.J.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
         TWOR = RKLR(INZ) * SignKL
         TWOI = RKLI(INZ) * SignKL
         IJ = (I-1) * NORB + J
         KL = (LTR-1) * NORB + KTR
         JI = (J-1) * NORB + I
         LK = (KTR-1) * NORB + LTR
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1) = TWOR
            TII(IJ+(KL-KLS-1)*M1) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1) = TWOR
            TII(JI+(LK-KLS-1)*M1) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1) = TWOR
            TII(KL+(IJ-KLS-1)*M1) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1) = TWOR
            TII(LK+(JI-KLS-1)*M1) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
C--> Time reversal both
         IF (ITR.NE.0.AND.JTR.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
         TWOR = RKLR(INZ) * SignIJ * SignKL
         TWOI = RKLI(INZ) * SignIJ * SignKL
         IJ = (JTR-1) * NORB + ITR
         KL = (LTR-1) * NORB + KTR
         JI = (ITR-1) * NORB + JTR
         LK = (KTR-1) * NORB + LTR
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1) = TWOR
            TII(IJ+(KL-KLS-1)*M1) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1) = TWOR
            TII(JI+(LK-KLS-1)*M1) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1) = TWOR
            TII(KL+(IJ-KLS-1)*M1) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1) = TWOR
            TII(LK+(JI-KLS-1)*M1) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
      ENDDO
C
      GOTO 30
C
   40 IF (MINT.EQ.0) WRITE (6,1010) 'MDCINT'
C
      IF (.NOT.BREIT) GOTO 61
C
C     -----------------------
C     READ IN BREIT INTEGRALS
C     -----------------------
C
  50  CONTINUE
      IF (REALAR) THEN
         read (mdbint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (rklr(inz),inz=1,nz)
      ELSE
         read (mdbint,ERR=10010,END=10010) ikr,jkr,nz,
     &                (indk(inz),indl(inz),inz=1,nz),
     &                (rklr(inz),rkli(inz),inz=1,nz)
      ENDIF
      MINTB = MINTB + NZ
      If (ikr.EQ.0) GoTo 60 ! End of File
      i = indmo(kr(ikr))
      j = indmo(kr(jkr))
      itr = indmo(kr(-ikr))
      jtr = indmo(kr(-jkr))
      SignIJ = - SIGN(1,ikr) * SIGN(1,jkr)
      DO INZ = 1, NZ
         kkr = indk(inz)
         lkr = indl(inz)
         k = indmo(kr(kkr))
         l = indmo(kr(lkr))
         ktr = indmo(kr(-kkr))
         ltr = indmo(kr(-lkr))
         SignKL = - SIGN(1,kkr) * SIGN(1,lkr)
C--> Original integral
         IF (I.NE.0.AND.J.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
         TWOR = RKLR(INZ)
         TWOI = RKLI(INZ)
         IJ = (I-1) * NORB + J
         KL = (K-1) * NORB + L
         JI = (J-1) * NORB + I
         LK = (L-1) * NORB + K
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
            TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
            TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
            TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
            TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
C--> Time reversal left
         IF (ITR.NE.0.AND.JTR.NE.0.AND.K.NE.0.AND.L.NE.0) THEN
         TWOR = RKLR(INZ) * SignIJ
         TWOI = RKLI(INZ) * SignIJ
         IJ = (JTR-1) * NORB + ITR
         KL = (K-1) * NORB + L
         JI = (ITR-1) * NORB + JTR
         LK = (L-1) * NORB + K
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
            TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
            TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
            TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
            TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
C--> Time reversal right
         IF (I.NE.0.AND.J.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
         TWOR = RKLR(INZ) * SignKL
         TWOI = RKLI(INZ) * SignKL
         IJ = (I-1) * NORB + J
         KL = (LTR-1) * NORB + KTR
         JI = (J-1) * NORB + I
         LK = (KTR-1) * NORB + LTR
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
            TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
            TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
            TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
            TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
C--> Time reversal both
         IF (ITR.NE.0.AND.JTR.NE.0.AND.KTR.NE.0.AND.LTR.NE.0) THEN
         TWOR = RKLR(INZ) * SignIJ * SignKL
         TWOI = RKLI(INZ) * SignIJ * SignKL
         IJ = (JTR-1) * NORB + ITR
         KL = (LTR-1) * NORB + KTR
         JI = (ITR-1) * NORB + JTR
         LK = (KTR-1) * NORB + LTR
         IF (KL.GT.KLS.AND.KL.LE.KLF.AND.IJ.GE.KL) THEN
            TIR(IJ+(KL-KLS-1)*M1+MM) = TWOR
            TII(IJ+(KL-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (LK.GT.KLS.AND.LK.LE.KLF.AND.JI.GE.LK) THEN
            TIR(JI+(LK-KLS-1)*M1+MM) = TWOR
            TII(JI+(LK-KLS-1)*M1+MM) = - TWOI
         ENDIF
         IF (IJ.GT.KLS.AND.IJ.LE.KLF.AND.KL.GE.IJ) THEN
            TIR(KL+(IJ-KLS-1)*M1+MM) = TWOR
            TII(KL+(IJ-KLS-1)*M1+MM) = TWOI
         ENDIF
         IF (JI.GT.KLS.AND.JI.LE.KLF.AND.LK.GE.JI) THEN
            TIR(LK+(JI-KLS-1)*M1+MM) = TWOR
            TII(LK+(JI-KLS-1)*M1+MM) = - TWOI
         ENDIF
         ENDIF ! check on deleted spinors
      ENDDO
C
      GOTO 50
C
   60 IF (MINTB.EQ.0) WRITE (6,1010) 'MDBINT'
C
C     -----------------------------------------
C     ADD BREIT INTEGRALS AND COULOMB INTEGRALS
C     -----------------------------------------
C
      DO IJKL = 1, (KLF-KLS) * M1
         TIR(IJKL) = TIR(IJKL) + TIR(IJKL+MM)
         TII(IJKL) = TII(IJKL) + TII(IJKL+MM)
      ENDDO
C
   61 CONTINUE
C
C     ----------------------------------
C     WRITE 2-ELECTRON INTEGRALS TO FILE
C     ----------------------------------
C
      IJKL = 0
      DO 200 KL = KLS + 1, KLF
      TIR (KL+(KL-KLS-1)*M1) = 0.5D0 * TIR(KL+(KL-KLS-1)*M1)
      TII (KL+(KL-KLS-1)*M1) = 0.5D0 * TII(KL+(KL-KLS-1)*M1)
      NKL = 0
      DO 100 IJ = 1, M1
         IJKL = IJKL + 1
         IF (ABS(TIR(IJKL))+ABS(TII(IJKL)).GT.TRESINT) THEN
            NKL = NKL + 1
            INDK(NKL) = IJ
            RKLR(NKL) = TIR(IJKL)
            RKLI(NKL) = TII(IJKL)
          ENDIF
  100     CONTINUE
      NINT = NINT + NKL
      IF (REALAR) THEN
         WRITE (MRCTWOE) KL,NKL,(INDK(IX),IX=1,NKL),(RKLR(IX),IX=1,NKL)
      ELSE
         WRITE (MRCTWOE) KL,NKL,(INDK(IX),IX=1,NKL),(RKLR(IX),IX=1,NKL),
     &                   (RKLI(IX),IX=1,NKL)
      ENDIF
  200 CONTINUE
      KLS = KLF
      KLF = MIN0 (M1, KLF+M2)
      IF (KLS.NE.M1) GOTO 10
      IF (.NOT.BREIT) THEN
         WRITE (6,1000) MINT,DATEX,TIMEX,NINT
      ELSE
         WRITE(6,1001) MINT,DATEX,TIMEX,MINTB,DATEXB,TIMEXB,NINT
      ENDIF
      IF (NINT.GT.0) THEN
         WRITE (MRCTWOE) -NINT,0
      ELSE
         WRITE (MRCTWOE) -1,0
      ENDIF
      IF (REALAR) THEN
         WRITE (6,1002)
      ELSE
         WRITE (6,1003)
      ENDIF
C
      CLOSE (MDCINT,STATUS='KEEP')
      IF (BREIT) CLOSE (MDBINT,STATUS='KEEP')
      RETURN
 1000 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Coulomb integrals written to MRCTWOE')
 1001 FORMAT (
     &/I8,' Unique Coulomb integrals read from MDCINT (',A10,1X,A8,')'
     &/I8,' Unique Breit integrals read from MDBINT (',A10,1X,A8,')'
     &/I8,' Combined integrals written to MRCTWOE')
 1002 FORMAT (/' Real arithmetic will be used in EXPAND and TWOLIN')
 1003 FORMAT (/' Complex arithmetic will be used')
 1010 FORMAT (//' CAUTION : No non-zero integrals found on ',A//)
10000 CALL QUIT('ERROR READING HEADER OF MDCINT')
10001 CALL QUIT('ERROR READING HEADER OF MDBINT')
10010 CALL QUIT('ERROR READING INTEGRALS FROM MDCINT')
10011 CALL QUIT('ERROR READING INTEGRALS FROM MDBINT')
      END
C
      SUBROUTINE TRIAL (VECTR,VECTI,IDIND,JDET,JELEC,JSIND,
     &                  IVWU,DGEL)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "general.inc"
#include "files.inc"
#include "symm.inc"
#include "iterat.inc"
C
      PARAMETER (NDREC=1023)
      PARAMETER (LREC=20*(NDREC+1))
C64B  PARAMETER (LREC=24*(NDREC+1))
C
      DIMENSION IDIND(*),JDET(*),JELEC(*)
      DIMENSION VECTR(NDET,*),VECTI(NDET,*),DGEL(NDET)
C
      CHARACTER*11 VECFIL
C
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=LREC)
C
      IF (ISTART.GT.1) GOTO 300
      VECFIL = FNAM(4)(1:6)//'_'//REPNA(IREP)
      IF (VECFIL(8:8).EQ.'''') VECFIL(8:8) = '1'
      IF (VECFIL(8:8).EQ.'"') VECFIL(8:8) = '2'
      IOS = 0
      OPEN (MRCTRIV,FILE=VECFIL,ACCESS='DIRECT',STATUS='OLD'
     &      ,RECL=4,ERR=300,IOSTAT=IOS)
C64B &      ,RECL=8,ERR=300,IOSTAT=IOS)
      IF (IOS.NE.0) GOTO 300
      READ (MRCTRIV,REC=1,ERR=300,IOSTAT=IOS) LENREC
      IF (IOS.NE.0) GOTO 300
      CLOSE (MRCTRIV)
      OPEN (MRCTRIV,FILE=VECFIL,ACCESS='DIRECT',RECL=LENREC,ERR=9999)
      READ (MRCTRIV,REC=1,ERR=9999) LENREC,NTRIV
      READ (MRCTRIV,REC=2,ERR=9999) NCOEF,(JDET(J),J=1,NCOEF)
C
C     ----------------------------------------------------------
C     CHECK IF THE TRIAL VECTORS CAN BE USED IN THIS CALCULATION
C     ----------------------------------------------------------
C
CTROND      NIDET = NOVERI(NORBR(2),NELEC-NORBR(1))
      NIDET = BICO(NORBR(2),NELEC-NORBR(1))
      IF (NIDET.LT.NCOEF.OR.NDET.LT.NCOEF)
     &CALL Q2JOB (2,'TRIAL VECTORS CAN NOT BE USED',' ',0)
      WRITE (6,1010)
      CALL DETIND (NCOEF,IDIND,JDET,JELEC,JSIND,IVWU)
C
C     ------------------------------------------------------------
C     WRITE REFERENCE DETERMINANTS TO FILE FOR DAVIDSON CORRECTION
C     ------------------------------------------------------------
C
      WRITE(MRCFINV,REC=2) NCOEF,(IDIND(J),J=1,NCOEF)
C
      DO 200 IROOT = 1, NROOTS
	 IVEC = IROOT
         CALL ZZERO (NDET,VECTR(1,IVEC))
         IF (.NOT.REALAR) CALL ZZERO (NDET,VECTI(1,IVEC))
         ISEL = NSEL(IROOT)
         IF (ISEL.GT.NTRIV) CALL Q2JOB
     &      (2,'TRIAL VECTOR NOT ON FILE',' ',0)
         IREC = ISEL + 2
         READ (MRCTRIV,REC=1,ERR=9999) LENREC,NTRIV,(DUM,I=1,ISEL-1),
     &                                 EIGEN1(IROOT)
         EIGEN2(IROOT) = EIGEN1(IROOT)
         IF (REALAR) THEN
            READ (MRCTRIV,REC=IREC,ERR=9999)
     &       (VECTR(IDIND(J),IVEC),VDUM,J=1,NCOEF)
         ELSE
            READ (MRCTRIV,REC=IREC,ERR=9999)
     &       (VECTR(IDIND(J),IVEC),VECTI(IDIND(J),IVEC),J=1,NCOEF)
         ENDIF
         NVEC = NVEC + 1
         CALL ORTHN (VECTR,VECTI)
         WRITE (6,1000) NSEL(IROOT),EIGEN1(IROOT)
  200    CONTINUE
      RETURN
 1000 FORMAT (/' Vector #',I4,' read from MRCTRIV; Eigenvalue ',F20.12)
 1010 FORMAT (/' WARNING !    Trial vector coefficients adapted'/
     &' Be sure that : NORBR(2,OLD)=NORBR(2,NEW) AND ',
     &'NORBR(1,OLD)=0')
C     ---------------------------------------------------------------
C     NO CORRECT FILE FOUND : MAKE START VECTOR BY LOOKING FOR LOWEST
C     DIAGONAL ELEMENTS.
C     ---------------------------------------------------------------
  300 CONTINUE
      WRITE (6,1020)
      IF (ISTART.NE.3) THEN
         DO 320 IROOT = 1, NROOTS
            IVEC = IROOT
            INDMIN(IROOT) = IDMIN(NDET,DGEL,1)
C64B        INDMIN(IROOT) = ISMIN(NDET,DGEL,1)
            EIGEN1(IROOT) = DGEL(INDMIN(IROOT)) + ECORE
            EIGEN2(IROOT) = EIGEN1(IROOT)
            CALL ZZERO (NDET,VECTR(1,IVEC))
            IF (.NOT.REALAR) CALL ZZERO (NDET,VECTI(1,IVEC))
            VECTR(INDMIN(IROOT),IVEC) = DGEL(INDMIN(IROOT))
            DGEL(INDMIN(IROOT)) = 0.D0
            WRITE (6,1030) IROOT,INDMIN(IROOT),EIGEN1(IROOT)
  320    CONTINUE
C
         DO 340 IROOT = 1, NROOTS
	    IVEC = IROOT
            DGEL(INDMIN(IROOT)) = VECTR(INDMIN(IROOT),IVEC)
            VECTR(INDMIN(IROOT),IVEC) = 1.D0
            CALL PUTVEC (IVEC,VECTR(1,IVEC),VECTI(1,IVEC))
  340    CONTINUE
         WRITE(MRCFINV,REC=2) NROOTS,(INDMIN(J),J=1,NROOTS)
         ISTART=2
      ELSE
         DO 350 IROOT = 1, NROOTS
            IVEC = IROOT
            INDMIN(IROOT) = IROOT
            CALL ZZERO (NDET,VECTR(1,IVEC))
            IF (.NOT.REALAR) CALL ZZERO (NDET,VECTI(1,IVEC))
            VECTR(IROOT,IVEC) = 1.D0
            CALL PUTVEC (IVEC,VECTR(1,IVEC),VECTI(1,IVEC))
            EIGEN1(IROOT) = DGEL(IROOT) + ECORE
            EIGEN2(IROOT) = EIGEN1(IROOT)
            WRITE (6,1030) IROOT,IROOT,EIGEN1(IROOT)
  350    CONTINUE
         ISTART=3
         WRITE(MRCFINV,REC=2) NROOTS,(INDMIN(J),J=1,NROOTS)
      ENDIF
      NVEC = NROOTS
C
      CLOSE(MRCFINV,STATUS='KEEP')
C
 1020 FORMAT
     & (/' Start vectors are single determinants')
 1030 FORMAT (/' Root',i4,' is determinant',i4,' Exp. Value',F20.12)
      RETURN
 9999 CALL Q2JOB (1,'FILE MRCTRIV',' ',0)
      RETURN
      END
C
      SUBROUTINE TWOLIN (ENDFIL,TINTR,TINTI,IJEX,KL,NKL)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ================================
C     READS BLOCK OF INTEGRALS (**|KL)
C     ** : ALL IJ >= KL
C     ================================
C
#include "param.inc"
#include "general.inc"
#include "files.inc"
C
      DIMENSION TINTR(*),TINTI(*),IJEX(*)
      LOGICAL ENDFIL
C
      ENDFIL = .FALSE.
      IF (REALAR) THEN
         READ (MRCTWOE) KL,NKL,(IJEX(IX),IX=1,NKL),(TINTR(IX),IX=1,NKL)
      ELSE
         READ (MRCTWOE) KL,NKL,(IJEX(IX),IX=1,NKL),(TINTR(IX),IX=1,NKL),
     &                  (TINTI(IX),IX=1,NKL)
      ENDIF
      IF (KL.LT.0) ENDFIL=.TRUE.
C     ******************************************************
C     THIS PART IS DONE TO MAKE ALL 2-ELECTRON INTEGRAL ZERO
C     Proposed & Implimented by (Malaya K. Nayak)
C     ******************************************************
      IF (PROPER) THEN
        CALL ZZERO(NKL,TINTR)
        CALL ZZERO(NKL,TINTI)
      END IF
C     ******************************************************
C     END OF MAKING ALL 2-ELECTRON INTEGRAL ZERO
C     ******************************************************
      RETURN
      END
C
      SUBROUTINE USERIN
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
#include "symm.inc"
C
      LOGICAL TOBE
      CHARACTER*4 IREPNA,VECFIL*11
      DIMENSION NRAS1(16),NRAS2(16)
C
      NAMELIST /RASORB/ NELEC,NRAS1,NRAS2,MINH1,MAXH1,MAXE3
      NAMELIST /DIRECT/ MAXITER,CPUMAX,CONVERR,CONVERE,RESTART
      NAMELIST /CIROOT/ NROOTS,NSEL,SELECT,IREPNA,ISTART
      NAMELIST /OPTIM/ IGENEX
      NAMELIST /LEADDET/ GETDET,COMIN
      NAMELIST /NATURAL/ MAKENAT
C
      NELEC = 0
      DO I = 1, NSYMRP
         NRAS1(I) = 0
         NRAS2(I) = 0
      ENDDO
      MAXH1 = 0
      MAXE3 = 0
      REWIND (5)
      READ (5,RASORB,END=100,ERR=100)
  100 MAXITER = MIN0(10,N2-1)
C
      NORBR(1) = 0
      NORBR(2) = 0
      DO I = 1, NSYMRP
         NORBR(1) = NORBR(1) + NRAS1(I)
         NORBR(2) = NORBR(2) + NRAS2(I)
         NRAS(1,I)= NRAS1(I)
         NRAS(2,I)= NRAS1(I) + NRAS2(I)
      ENDDO
      NORBR(3)= NORB - NORBR(1) - NORBR(2)
C
      IF (NORBR(1).EQ.0) MAXH1=0
      IF (NORBR(3).EQ.0) MAXE3=0
      RESTART = .FALSE.
      CPUMAX = 604800.D0
      CONVERR = 1.D-10
      CONVERE = 1.D-9
      REWIND (5)
      READ (5,DIRECT,END=200,ERR=200)
  200 SELECT = .FALSE.
      NROOTS = 1
      DO 210 ISEL = 1, N8
         NSEL(ISEL) = ISEL
  210    CONTINUE
C     --------------------------------------
C     FIND SUITABLE DEFAULT ABELIAN SYMMETRY
C     --------------------------------------
      IF (MOD(NELEC,2).EQ.0) THEN
         IREP1 = IRA1
         IREP2 = NREP
      ELSE
         IREP1 = 1
         IREP2 = IRA1 - 1
      ENDIF
      WRITE(6,'(A)') '* List of allowed representations:'
      WRITE(6,'(3X,A4)') (REPNA(IREP),IREP = IREP1,IREP2)
      DO 220 IREP = IREP1, IREP2
         VECFIL = FNAM(4)(1:6)//'_'//REPNA(IREP)
         IF (VECFIL(8:8).EQ.'''') VECFIL(8:8) = '1'
         IF (VECFIL(8:8).EQ.'"') VECFIL(8:8) = '2'
         OPEN (MRCTRIV,FILE=VECFIL,ACCESS='DIRECT',RECL=4,ERR=220)
C64B     OPEN (MRCTRIV,FILE=VECFIL,ACCESS='DIRECT',RECL=8,ERR=220)
         IOS = 0
         READ (MRCTRIV,REC=1,ERR=220,IOSTAT=IOS) LENREC
         IF (IOS.NE.0) GOTO 220
         TOBE = .TRUE.
         IF (TOBE) GOTO 221
  220    CONTINUE
C     -------------------------------------------------------------
C     NO TRIAL VECTORS FOUND : TAKE FIRST REPRESENTATION AS DEFAULT
C     -------------------------------------------------------------
      IREP = IREP1
  221 CONTINUE
      IREPNA = REPNA(IREP)
      ISTART=1
      REWIND (5)
      READ (5,CIROOT,END=300,ERR=300)
  300 CONTINUE
C     IF (IREPNA(1:1).EQ.'1') IREPNA(1:1) = ''''
C     IF (IREPNA(1:1).EQ.'2') IREPNA(1:1) = '"'
      DO 310 IREP = IREP1, IREP2
         IF (REPNA(IREP).EQ.IREPNA) GOTO 311
  310 CONTINUE
      CALL Q2JOB (2,'REPRESENTATION NOT SPANNED',' ',0)
  311 CONTINUE
      NORBRT = NORBR(1) + NORBR(2) + NORBR(3)
      IGENEX = 2
      REWIND (5)
      READ (5,OPTIM,END=400,ERR=400)
  400 GETDET=.TRUE.
      COMIN=0.1D0
      REWIND (5)
      READ (5,LEADDET,END=500,ERR=500)
  500 CONTINUE
#ifdef PRG_DIRAC
C     Make natural spinors does not work for DIRAC (only for MOLFDIR)
      MAKENAT=.FALSE.
#else
      MAKENAT=.TRUE.
#endif
      REWIND (5)
      READ (5,NATURAL,END=600,ERR=600)
  600 CONTINUE
      IF (NORBRT.NE.NORB) CALL Q2JOB (2,'SUM(NORBR).NE.NORB',' ',0)
      IF (NELEC.EQ.0)     CALL Q2JOB (2,'NELEC.EQ.0',' ',0)
      IF (NELEC.GT.NORB)  CALL Q2JOB (2,'NELEC.GT.NORB',' ',0)
      IF (NELEC.GT.N4)    CALL Q2JOB (3,'P.USERIN','    N4',NELEC)
      IF (NROOTS.GT.N8)   CALL Q2JOB (3,'P.USERIN','    N8',NROOTS)
      RETURN
      END
C
      SUBROUTINE USEROUT (CI)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "general.inc"
#include "symm.inc"
#include "files.inc"
#include "iterat.inc"
#include "mempoint.inc"
C
      REAL*8 CI(*)

      DIMENSION OV(N2)
C
      CALL REAVEC (CI(IPVR),CI(IPVI),CI(IPDR),CI(IPDI),CI(IPIL))
C
      DO 40 IROOT = 1, NROOTS
C        ------------------------------------
C        Calculate first order density matrix
C        ------------------------------------
         CALL DENSMT (CI(IPDR),CI(IPDI),CI(IPVR),
     &                CI(IPVI),CI(IPIL),CI(IPIR),
     &                CI(IPJL),CI(IPJR),
     &                CI(IPJV),CI(IIJR),CI(ITNR),
     &                CI(ITNI),CI(IVW1),IROOT,
     &                CI(IIJE+NORB*(IROOT-1)))
C        -------------------------
C        Calculate natural Spinors
C        -------------------------
         IF (MAKENAT) CALL NATORB(IROOT,CI(ITNR), 
     &                            CI(ITNI),CI(IGMR),
     &                            CI(IGMI),CI(IPJV))
   40 CONTINUE
C
      WRITE (6,1020) (IROOT,IROOT=1,NROOTS)
      DO 30 I = 1, NORB
         IRP = IRPAMO(I)
      WRITE (6,1030) I,REPNA(IRP),(CI(IIJE+NORB*(IROOT-1)+I-1),
     &                                      IROOT=1,NROOTS)
  30  CONTINUE
C
      DO 20 IROOT = 1, NROOTS
      ISEL = NSEL(IROOT)
      DO 10 I = 1, NVEC
      OV(I) = EVECR(I,ISEL)*EVECR(I,ISEL) 
      IF (.NOT.REALAR) 
     &   OV(I) = OV(I) + EVECI(I,ISEL)*EVECI(I,ISEL)
   10 CONTINUE
      WRITE (6,1010) IROOT,(OV(I),I=1,NROOTS)
   20 CONTINUE
C
C     Give an overview of the leading determinants
C
      IF (GETDET) CALL LEADING(CI(IVW1))
      WRITE (6,1000)
C
 1000 FORMAT(//10X,'-- Normal end of program --')
 1010 FORMAT(//' Contributions of start vectors to final CI vector',I4
     &/(5(2X,G14.8)))
 1020 FORMAT(/' Occupation numbers : '//' Orbital Representation',
     &5(3X,'Root',I3))
 1030 FORMAT(I6,8X,A4,5X,5(1X,F9.6))
C
      RETURN
      END
C
      SUBROUTINE VERTEX(IVWU)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
C     ========================================================
C     WALKS THROUGH DIAGRAM OF THIS RESTRICTED ACTIVE SPACE TO
C     CALCULATE VERTEX WEIGHTS
C     ========================================================
C
#include "param.inc"
#include "symm.inc"
#include "general.inc"
C
      DIMENSION IVWU(2*MXREP,0:NORB,0:NELEC)
C
      DO 30 I = 0, NORB
      DO 20 J = 0, NELEC
      DO 10 IR = 1, NREP
      IVWU (IR,I,J) = 0
   10 CONTINUE
   20 CONTINUE
   30 CONTINUE
C     ===============================
C     DO THE WALK FROM TOP TO BOTTOM
C     ===============================
      ICOL = MAXH1
      IF (NORBR(1).EQ.0) ICOL = NORBR(2) - NELEC + MAXE3
      DO 50 I = 0, ICOL 
      IVWU(IRA1,I,0) = 1
   50 CONTINUE
C
      IORB = 0
      IORBR(0)=0
C
      IF (NORBR(1).EQ.0) GOTO 1101
C
      DO 100 I = 1, NORBR(1)
      IORB = IORB + 1
      JMIN = MAX0( I - MAXH1, 1)
      IORBR(IORB) = MAX0( I - MAXH1, 0)
      JMAX = IORB
      DO 90  J = JMIN, JMAX
      DO 80 IR = 1, NREP
      JR = MULTB(IR,IRPAMO(IORB))
      IVWU(JR,IORB,J) = IVWU(JR,IORB-1,J) + IVWU(IR,IORB-1,J-1)
   80 CONTINUE
   90 CONTINUE
  100 CONTINUE
C
 1101 IF (NORBR(2).EQ.0) GOTO 2101
C
      DO 200 I = 1, NORBR(2)
      IORB = IORB + 1
      JMAX = MIN0 ( NELEC, IORB)
      JMIN = MAX0 ( NELEC - MAXE3 + I - NORBR(2), NORBR(1) - MAXH1)
      IORBR(IORB) = JMIN 
      DO 190  J = JMIN, JMAX
      IF (J.EQ.0) THEN
         IF (NORBR(1).GT.0) IVWU(IRA1,IORB,J) = 1
         GOTO 190
      ENDIF
      DO 180 IR = 1, NREP
      JR = MULTB(IR,IRPAMO(IORB))
      IVWU(JR,IORB,J) = IVWU(JR,IORB-1,J) + IVWU(IR,IORB-1,J-1)
  180 CONTINUE
  190 CONTINUE
  200 CONTINUE
C
 2101 IF (NORBR(3).EQ.0) GOTO 3001
C
      DO 300 I = 1, NORBR(3)
      IORB = IORB + 1
      JMIN = MAX0 ( NELEC + I - NORBR(3), NELEC - MAXE3)
      IORBR(IORB) = JMIN 
      JMAX = NELEC
      DO 290  J = JMIN, JMAX
      IF (J.EQ.0) THEN
         IVWU(IRA1,IORB,J) = 1
         GOTO 290
      ENDIF
      DO 280 IR = 1, NREP
      JR = MULTB(IR,IRPAMO(IORB))
      IVWU(JR,IORB,J) = IVWU(JR,IORB-1,J) + IVWU(IR,IORB-1,J-1)
  280 CONTINUE
  290 CONTINUE
  300 CONTINUE
      
C     ======================================================
C     DEFINE MAXIMUM NUMBER OF SYMMETRY ADAPTED DETERMINANTS 
C     ======================================================
 3001 NREFDET = IVWU(IREP,NORBR(1)+NORBR(2),NELEC)
      NDET = IVWU(IREP,IORB,NELEC)
      MDET = NDET
      DO 310 IR = 1, NREP
      IF (IVWU(IR,IORB,NELEC).GT.MDET) MDET = IVWU(IR,IORB,NELEC)
  310 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE WRIVEC (VECTR,VECTI,DR,DI,ILIND,IRIND)
C
      IMPLICIT REAL*8 (A-H, O-Z)
C
#include "param.inc"
#include "files.inc"
#include "general.inc"
#include "iterat.inc"
C
      DIMENSION VECTR(NDET,*),VECTI(NDET,*),DR(*),DI(*)
      DIMENSION ILIND(*),IRIND(*)
      DIMENSION IDREC(N8)
C
      PARAMETER (NDREC=1023)
      PARAMETER (LENREC=20*(NDREC+1))
C64B  PARAMETER (LENREC=24*(NDREC+1))
C
      OPEN (MRCFINV,FILE=FNAM(5),ACCESS='DIRECT',RECL=LENREC)
      IREC = 3
      IF (ISTART.EQ.1) THEN
         READ(MRCFINV,REC=2) NCOEF,(IRIND(J),J=1,NCOEF)
      ELSE
         READ(MRCFINV,REC=2) NCOEF,(INDMIN(J),J=1,NCOEF)
      ENDIF
      DO IROOT = 1, NROOTS
      IVEC = IROOT
      CALL ZZERO (NDET,VECTR(1,IVEC))
      IF (.NOT.REALAR) CALL ZZERO (NDET,VECTI(1,IVEC))
      ENDDO
C
      KVEC = NROOTS + 1
      DO 20 N = 1, NVEC
      CALL  GETVEC (N,VECTR(1,KVEC),VECTI(1,KVEC))
C
      DO 10 IROOT = 1, NROOTS
      IVEC = IROOT
      ER = EVECR(N,NSEL(IROOT))
      CALL DAXPY (NDET,ER,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
C64B  CALL SAXPY (NDET,ER,VECTR(1,KVEC),1,VECTR(1,IVEC),1)
      IF (.NOT.REALAR) THEN
         EI = EVECI(N,NSEL(IROOT))
         CALL DAXPY (NDET,EI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
         CALL DAXPY (NDET,ER,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
         CALL DAXPY (NDET,-EI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
C64B     CALL SAXPY (NDET,EI,VECTR(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,ER,VECTI(1,KVEC),1,VECTI(1,IVEC),1)
C64B     CALL SAXPY (NDET,-EI,VECTI(1,KVEC),1,VECTR(1,IVEC),1)
      ENDIF
   10 CONTINUE
   20 CONTINUE
      DO 30 IROOT = 1, NROOTS
      IVEC = IROOT
      NONZERO = 0
      IDREC(IROOT) = IREC
      DO 11 ID = 1, NDET
         IF (VECTR(ID,IVEC).NE.0.D0.OR.VECTI(ID,IVEC).NE.0.D0) THEN
            IF (NONZERO.EQ.NDREC) THEN
               IF (REALAR) THEN
                  WRITE (MRCFINV,REC=IREC) IREC + 1, NONZERO,
     &                  (ILIND(IDNZ),IDNZ=1,NONZERO),
     &                  (DR(IDNZ),IDNZ=1,NONZERO)
               ELSE
                  WRITE (MRCFINV,REC=IREC) IREC + 1, NONZERO,
     &                  (ILIND(IDNZ),IDNZ=1,NONZERO),
     &                  (DR(IDNZ),DI(IDNZ),IDNZ=1,NONZERO)
               ENDIF
               IREC = IREC + 1
               NONZERO = 0
            ENDIF
            NONZERO = NONZERO + 1
            ILIND(NONZERO) = ID
            DR(NONZERO) = VECTR(ID,IVEC)
            IF (.NOT.REALAR) DI(NONZERO) = VECTI(ID,IVEC)
         ENDIF
   11    CONTINUE
      IF (REALAR) THEN
         WRITE (MRCFINV,REC=IREC) 0, NONZERO,
     &         (ILIND(IDNZ),IDNZ=1,NONZERO),
     &         (DR(IDNZ),IDNZ=1,NONZERO)
      ELSE
         WRITE (MRCFINV,REC=IREC) 0, NONZERO,
     &         (ILIND(IDNZ),IDNZ=1,NONZERO),
     &         (DR(IDNZ),DI(IDNZ),IDNZ=1,NONZERO)
      ENDIF
      IREC = IREC + 1
C
C     Calculate norm reference part op the vector
C
      RES = 0.0D0
C64B  RES = 0.0E0
      IF (ISTART.EQ.1) THEN
        DO 12 IREF = 1, NCOEF
          RES = RES + VECTR(IRIND(IREF),IVEC) * VECTR(IRIND(IREF),IVEC)
          IF (.NOT.REALAR)
     &    RES = RES + VECTI(IRIND(IREF),IVEC) * VECTI(IRIND(IREF),IVEC)
  12    CONTINUE
      ELSE
        RES = RES + VECTR(INDMIN(IVEC),IVEC) * VECTR(INDMIN(IVEC),IVEC)
      ENDIF
      REFVEC(IROOT) = RES
C
   30 CONTINUE
C
      WRITE (MRCFINV,REC=1) LENREC,NDET,NROOTS,
     &                      ((EIGEN(NSEL(IROOT))+ECORE),IROOT=1,NROOTS),
     &                      (IDREC(IROOT),IROOT=1,NROOTS)
      CLOSE (MRCFINV,STATUS='KEEP')
C
      RETURN
      END
C     ***************************************************
C     OPEN AND READ PROPERTY INTEGRALS FROM FILE 'MDPROP'
C     Proposed & Written by (Malaya K. Nayak)
C     ***************************************************
C
      SUBROUTINE PROP_INT(PROP,AIJ,MORB,IRW,APHASE)
      IMPLICIT REAL*8 (A-H, O-Z)

#include "param.inc"
#include "files.inc"
#include "general.inc"
#include "symm.inc"

      PARAMETER (ACCUR=1.D-15)
      DIMENSION PROP(2,MORB,MORB),AIJ(2*MORB*MORB)
      DIMENSION APHASE(2)
      CHARACTER*32 ACHAR
      LOGICAL PHASE

      OPEN (MDPROP,FILE=FNAM(9),FORM='UNFORMATTED')
      INOP = 0
    1 READ (MDPROP,ERR=10,END=11) ACHAR
      IF (ACHAR(1:8).NE.'********'.OR.ACHAR(25:32).NE.NAMEA) GOTO 1
C      IF (IPRNT.GE.1) WRITE (IW,1000) NAMEA,ACHAR(9:16),ACHAR(17:24)
      WRITE (6,1000) NAMEA,ACHAR(9:16),ACHAR(17:24)
      READ (MDPROP) (((PROP(K,I,J),K=1,2),I=1,MORB),J=1,MORB)
      CLOSE (MDPROP,STATUS='KEEP')
      GOTO 12
   10 INOP = 1
      GOTO 12
   11 INOP = 2
   12 CONTINUE
C
      IF (INOP.EQ.1) GOTO 101
      IF (INOP.EQ.2) GOTO 102
CMKN
C     *************************************************************
C      IF ((ACHAR(25:32).EQ.'HYPER11') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER12') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER21') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER22')) THEN
C      WRITE(*,*)'ACHAR=',ACHAR
C      WRITE(*,*)'MORB=',MORB
C      MBY2 = MORB/2
C      DO I=1,MORB
C      DO J=1,MORB
C      TEMP(1,I,J) = 0.0D0
C      TEMP(2,I,J) = 0.0D0
C      END DO
C      END DO
C      DO I=1,MORB
C      DO J=1,MORB
C      IF (I .LE. MBY2 .AND. J .LE. MBY2) THEN
C      K = I
C      L = J+MBY2
C      TEMP(1,K,L) = PROP(1,I,J)
C      TEMP(2,K,L) = PROP(2,I,J)
C      IF ((ACHAR(25:32).EQ.'HYPER11') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER21')) THEN
C      TEMP(1,K,L) = -PROP(1,I,J)
C      TEMP(2,K,L) = -PROP(2,I,J)
C      END IF
C      ELSE IF (I .LE. MBY2 .AND. J .GT. MBY2) THEN
C      K = I
C      L = J-MBY2
C      TEMP(1,K,L) = PROP(1,I,J)
C      TEMP(2,K,L) = PROP(2,I,J)
C      ELSE IF (I .GT. MBY2 .AND. J .LE. MBY2) THEN
C      K = I
C      L = J+MBY2
C      TEMP(1,K,L) = PROP(1,I,J)
C      TEMP(2,K,L) = PROP(2,I,J)
C      IF ((ACHAR(25:32).EQ.'HYPER11') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER21')) THEN
C      TEMP(1,K,L) = -PROP(1,I,J)
C      TEMP(2,K,L) = -PROP(2,I,J)
C      END IF
C      ELSE IF (I .GT. MBY2 .AND. J .GT. MBY2) THEN
C      K = I
C      L = J-MBY2
C      TEMP(1,K,L) = PROP(1,I,J)
C      TEMP(2,K,L) = PROP(2,I,J)
C      END IF
C      END DO
C      END DO
C      DO I = 1,MORB
C      DO J = 1,MORB
C      PROP(1,I,J) = TEMP(1,I,J)
C      PROP(2,I,J) = TEMP(2,I,J)
C      WRITE(*,*)I,J,PROP(1,I,J),PROP(2,I,J)
C      END DO
C      END DO
C      END IF
C      ***********************************************************
C      IF ((ACHAR(25:32).EQ.'HYPER11') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER12') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER21') .OR.
C     &    (ACHAR(25:32).EQ.'HYPER22')) THEN
C
C      IF ((ACHAR(25:27).EQ.'HYP'.AND.ACHAR(29:30).NE.'13') .OR.
C     &    (ACHAR(25:32).EQ.'HYP'.AND.ACHAR(29:30).NE.'23')) THEN
C      ***********************************************************
CMKN
C      IF (ACHAR(25:27).EQ.'HYP'.AND.ACHAR(30:30).NE.'3') THEN
      IF ((ACHAR(28:30).EQ.'HYP'.OR.ACHAR(28:30).EQ.'ANA').AND.
     &    (ACHAR(25:25).NE.'Z')) THEN
      WRITE(6,*)'Diagonal and Off-diagonal blocks are exchanged'
      MBY2 = MORB/2
      DO I=1,MORB
      DO J=1,MORB
      TEMPR = 0.0D0
      TEMPI = 0.0D0
      IF (J .LE. MBY2) THEN
      K = I
      L = J+MBY2
      TEMPR = PROP(1,I,J)
      TEMPI = PROP(2,I,J)
      PROP(1,I,J) = PROP(1,K,L)
      PROP(2,I,J) = PROP(2,K,L)
C      IF (ACHAR(30:30).EQ.'2') THEN
      IF ((ACHAR(28:30).EQ.'HYP'.AND.ACHAR(25:25).EQ.'Y').OR.
     &    (ACHAR(28:30).EQ.'ANA'.AND.ACHAR(25:25).EQ.'X')) THEN
      PROP(1,K,L) = TEMPR
      PROP(2,K,L) = TEMPI
      ELSE
      PROP(1,K,L) = TEMPR
      PROP(2,K,L) =-TEMPI
      END IF
      END IF
      END DO
      END DO
C      IF (ACHAR(30:30).EQ.'1') WRITE(6,*)
      IF ((ACHAR(28:30).EQ.'HYP'.AND.ACHAR(25:25).EQ.'X').OR.
     &    (ACHAR(28:30).EQ.'ANA'.AND.ACHAR(25:25).EQ.'Y')) THEN
      WRITE(6,*)'Sign of the lower diagonal block is changed'
      END IF
      END IF
CMKN
CMKN
C      IF (ACHAR(25:27).EQ.'HYP'.AND.ACHAR(30:30).NE.'3') THEN
      IF ((ACHAR(28:31).EQ.'COMP').AND.(ACHAR(26:26).NE.'Z')) THEN
      WRITE(6,*)'Diagonal and Off-diagonal blocks are exchanged'
      MBY2 = MORB/2
      DO I=1,MORB
      DO J=1,MORB
      TEMPR = 0.0D0
      TEMPI = 0.0D0
      IF (J .LE. MBY2) THEN
      K = I
      L = J+MBY2
      TEMPR = PROP(1,I,J)
      TEMPI = PROP(2,I,J)
      PROP(1,I,J) = PROP(1,K,L)
      PROP(2,I,J) = PROP(2,K,L)
C      IF (ACHAR(30:30).EQ.'2') THEN
      IF ((ACHAR(28:31).EQ.'COMP'.AND.ACHAR(26:26).EQ.'Y')) THEN
      PROP(1,K,L) = TEMPR
      PROP(2,K,L) = TEMPI
      ELSE
      PROP(1,K,L) = TEMPR
      PROP(2,K,L) =-TEMPI
      END IF
      END IF
      END DO
      END DO
C      IF (ACHAR(30:30).EQ.'1') WRITE(6,*)
      IF ((ACHAR(28:31).EQ.'COMP'.AND.ACHAR(26:26).EQ.'X')) THEN
      WRITE(6,*)'Sign of the lower diagonal block is changed'
      END IF
      END IF
CMKN
C      ***********************************************************
C      WRITE(*,*)'ACHAR=',ACHAR
C      WRITE(*,*)'PROP=',(((PROP(K,I,J),K=1,2),I=1,MORB),J=1,MORB)
C      WRITE(*,*)'FIRST STRING=',((PROP(1,I,J),
C     &                  I=1,MORB),J=1,MORB)
C      WRITE(*,*)'SECOND STRING=',((PROP(2,I,J),
C     &                  I=1,MORB),J=1,MORB)
C      ***********************************************************
C      THIS PART TAKES CARE OF SIGMS_Z PART OF THE OPERATOR
C      IF THE OPERATOR IS DIAGONAL & IT HAS A SIGMA_Z PART
C      THIS IS VALID, ONLY WHEN WE USE THE 'C1' SYMMETRY
C      Proposed & Implimented by (Malaya K. Nayak)
C      ***********************************************************
C      IF ((ACHAR(25:32).EQ.'NEF 03') .OR.
C     &    (ACHAR(25:32).EQ.'NEF 06')) THEN
C      WRITE(*,*)'ACHAR=',ACHAR
C      MBY2 = MORB/2
C      DO I=1,MORB
C      DO J=1,MORB
C      IF (I .LE. MBY2 .AND. J .LE. MBY2) THEN
C      PROP(1,I,J) = PROP(1,I,J)*(1.0D0/2.0D0)
C      PROP(2,I,J) = PROP(2,I,J)*(1.0D0/2.0D0)
C      ELSE IF (I .GT. MBY2 .AND. J .GT. MBY2) THEN
C      PROP(1,I,J) = PROP(1,I,J)*(-1.0D0/2.0D0)
C      PROP(2,I,J) = PROP(2,I,J)*(-1.0D0/2.0D0)
C      END IF
C      END DO
C      END DO
C      END IF
C      ***********************************************************
C      WRITE(*,*)'ACHAR=',ACHAR
C      WRITE(*,*)'PROP=',(((PROP(K,I,J),K=1,2),I=1,MORB),J=1,MORB)
C      WRITE(*,*)'FIRST STRING=',((PROP(1,I,J),
C     &                  I=1,MORB),J=1,MORB)
C      WRITE(*,*)'SECOND STRING=',((PROP(2,I,J),
C     &                  I=1,MORB),J=1,MORB)
C      ***********************************************************
CMKN
      APHASE(1) = 1.0D0
      APHASE(2) = 0.0D0
      PHASE = .FALSE.
C
    2 AMAXI = 0.0D0
      IAMAX = 0
      JAMAX = 0
      IA = 0
      DO I = 1,MORB
      DO J = 1,MORB
      IA = IA+1
      IF (.NOT.REALAR) THEN
      AIJ(IRW*IA-1) = PROP(1,I,J)
      AIJ(IRW*IA)   = PROP(2,I,J)
      ELSE
CMKN Below two line with a (-ve) sign are original
C      AIJ(IA) = PROP(1,I,J)*APHASE(1) -
C     &          PROP(2,I,J)*APHASE(2)
      AIJ(IA) = PROP(1,I,J)*APHASE(1) +
     &          PROP(2,I,J)*APHASE(2)
CMKN Above two lines with a (+ve) sign are changed
      AIJI    = PROP(1,I,J)*APHASE(2) +
     &          PROP(2,I,J)*APHASE(1)
      IF (ABS(AIJI) .GT. ABS(AMAXI)) THEN
      AMAXR = AIJ(IA)
      AMAXI = AIJI
      IAMAX = I
      JAMAX = J
      END IF
      END IF
      END DO
      END DO
C
      IF (REALAR.AND.ABS(AMAXI).GT.ACCUR) THEN
      IF (PHASE) THEN
      WRITE (6,1010) IAMAX,JAMAX,AMAXI
      CALL QUIT('USE COMPLEX ARITHMETICS')
      ELSE
CMKN Below two lines without ABS are original
C      APHASE(1) = AMAXR
C      APHASE(2) = AMAXI
      APHASE(1) = ABS(AMAXR)
      APHASE(2) = ABS(AMAXI)
CMKN Above two lines with the ABS are changed
      ANORM = SQRT(AMAXR*AMAXR + AMAXI*AMAXI)
      APHASE(1) = APHASE(1)/ANORM
      APHASE(2) = APHASE(2)/ANORM
      PHASE = .TRUE.
C      IF (IPRNT .GT. 1) WRITE (IW,1020) NAMEA,APHASE
      WRITE (6,1020) NAMEA,APHASE
      GO TO 2
      END IF
      END IF
C
      RETURN
  101 WRITE (6,*) ' Error reading property ',NAMEA,' on file MDPROP'
      CALL QUIT(' Error reading property integrals')
  102 WRITE (6,*) ' Property ',NAMEA,' not found on file MDPROP'
      CALL QUIT(' Property integrals missing')
 1000 FORMAT (/' Read integral type ',A8,' created ',
     & A8,' storage info : ',A8)
 1010 FORMAT (/' Largest imaginary part of matrix element',2I5,F10.2)
 1020 FORMAT (/' Property ',A8,' scaled with phase factor ',2F10.4)
      END

