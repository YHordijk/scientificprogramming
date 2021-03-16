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

! FILE    : hercen.F
!
!------------------------------------------------------------------------------------
!  /* Deck ceninp */
      SUBROUTINE CENINP(LUINFO,WORK,LWORK,IQM,NBLCK,                       &
     &                  JCO,BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,   &
     &                  KATOM,KANG,KSETS,KBLOCK,KPRIM,ANG)
!*****************************************************************************
!
!       Loop over expansion centers and read coordinates and basis set information.
!
!       Revised version of BASINP, L. Visscher, december 2003.
!
!       Description of some of the arguments :
!
!       INTENT IN :
!       LUINFO : File from which the input is read (MOLECULE.INP, the .mol file in pam)
!       NONTYP : Number of different atom types that were defined by the user.
!       KATOM  : Dimension of the atom type arrays (maximum number of atom types)
!       KANG   : Dimension of the angular momentum arrays (maximum angular momentum)
!       KPRIM  : Dimension of the primitive function arrays (maximum number of primitives in a block)
!       KBLOCK : Dimension of the block arrays (maximum number of blocks that can be read)
!       KSETS  : Dimension of the basis set (maximum number of basis sets that we're supposed to read here)
!       ANG    : Convert coordinates from Angstrom to the internally used atomic units ?
!
!       INTENT INOUT (updated in routine):
!       NONT   : Number of symmetry unique atoms for each type.
!       IQM    : Highest anular momentum (+1 to start with 1 and not 0) for each type.
!       NBLCK  : Number of blocks read in total
!       JCO    : Number of blocks for each atom type and angular momentum
!       NUC    : Number of uncontracted functions "primitives" in each block
!       NRC    : Number of contracted functions "orbitals" in each block
!       ISGEN  : Defines generation of small component via kinetic balance in each (large component) block
!       ALPHA  : Value of the exponents of the primitive Gaussians in a block.
!       CPRIM  : Value of the contraction coefficients of the primitive Gaussians in a block.
!       CPRIMU : Value of the unnormalized contraction coefficients of the primitive Gaussians in a block.
!       SEG    : Is the segmented contraction scheme used in a block ?
!       BASREF : Lines of documentation for the basis set
!
!*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
#include "cbirea.h"
      DIMENSION WORK(LWORK)
      LOGICAL SEG,ANG,NOORBTS,CNTBAS
      DIMENSION IQM(KATOM,KSETS),JCO(KANG,KATOM,KSETS),                 &
     &          NBLCK(KATOM,KSETS),NUC(KBLOCK,KSETS),NRC(KBLOCK,KSETS), &
     &          SEG(KBLOCK,KSETS),                                      &
     &          ALPHA(KPRIM,KBLOCK,KSETS),                              &
     &          CPRIM(KPRIM,KPRIM,KBLOCK,KSETS),                        &
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,KSETS),ISGEN(KBLOCK,KSETS)
      CHARACTER*80 BASREF(10,KATOM,KSETS)
      LOGICAL GENERATE
      DIMENSION IBLCK(MXBSETS)
#include "molinp.h"
#include "nuclei.h"
#include "huckel.h"
#include "gencon.h"
#include "dcbgen.h"
!
!     Initialize variables that regard the whole basis set
!     e.g. whether it is contracted or not.
!
      NUCIND = 0
      NHTYP  = 0
      CALL IZERO (NBLCK,KATOM*KSETS)
      DO ISET = 1, KSETS
         DIRCON(ISET) = .FALSE.
         IBLCK(ISET) = 1
      END DO
!
!     Run over atomic centers
!     =======================
!
      DO I = 1, NONTYP
!
!*****************************************************************************
!       Card 6: ATOMIC DATA
!       Q     - nuclear charge
!       NONT  - number of symmetry independent atoms
!       IQM   - maximum angular quantum number (s=1,p=2,d=3 etc.)
!       NHTYP - maximum angular quantum number for ALL orbitals
!       JCO   - number of AO-blocks for each l-value
!*****************************************************************************
!
        CALL RDLINE(LUINFO)
!
!       Read/process data on centers/blocks
!       ===================================
!
        QEXP = 0.0D0
        READ (MLINE(NMLINE),'(BN,1X,F9.0,I5,F20.5)')  Q,NONT(I),QEXP
!
!       Read/process data on symmetry independent centers
!       =================================================
!
        IF (INPTST) THEN
           WRITE(LUPRI,'(/A,I5)') 'Atomic type: ',I
           WRITE(LUPRI,'(3X,A,F4.1)') 'Nuclear charge:  ',Q
           WRITE(LUPRI,'(3X,A,I5)') 'Number of symmetry '//             &
     &          'independent centres: ',NONT(I)
        END IF
!
!       Read coordinates of the centers of that belong to this atom type
!
        CALL CNTINP(LUINFO,NONT(I),Q,QEXP,ANG)
!
!       Loop over the basis sets that are to be associated with this atom type.
!
        DO ISET = 1, KSETS
           LWRK = LWORK
!          Tell the routine we're calling how many blocks remain.
           KAOVEC = KBLOCK - NBLCK(I,ISET)
           IBLOCK = IBLCK(ISET)
           CALL READ_BASIS (LUINFO,WORK,LWRK,Q,                         &
     &          IQM(I,ISET),NBLCK(I,ISET),CNTBAS,JCO(1,I,ISET),         &
     &          BASREF(1,I,ISET),                                       &
     &          NUC(IBLOCK,ISET),NRC(IBLOCK,ISET),SEG(IBLOCK,ISET),     &
     &          ALPHA(1,IBLOCK,ISET),CPRIM(1,1,IBLOCK,ISET),            &
     &          CPRIMU(1,1,IBLOCK,ISET),ISGEN(IBLOCK,ISET),             &
     &          KANG,ISET,KAOVEC,KPRIM,GENERATE)
           DIRCON(ISET) = DIRCON(ISET) .OR. CNTBAS
!
!          A set might be generated from a previous set for this atom type.
!
           IF (GENERATE) THEN
 
               IF (ISET.EQ.2) THEN
!                 Make (uncontracted) small component basis
!                 set by kinetic balance. We need to pass on the full
!                 array since the routine needs info from another set.
!                 IBLCK serves to find the starting points for this atom type
!                 in the block list
!
                  CALL KINBAL(I,IBLCK(1),IBLCK(2),IQM,JCO,NUC,NRC,SEG,  &
     &                 ALPHA,CPRIM,CPRIMU,NBLCK,ISGEN,KATOM,KANG,KBLOCK,&
     &                 KPRIM,INPTST,CNTBAS)
                  BASREF(1,I,ISET) = ' Derived from large component'
               ELSE
!                 TO DO : automatic generation of fit set ?
               END IF
           END IF
           print*,"Read basis set",iset,"for atom",i
           print*,"iqm",iqm(i,iset)
           print*,"iblock",iblock
           print*,"nblck",nblck(i,iset)
           print*,"generate",generate
           do ibx = 1, nblck(i,jset)
              print*,jset,ibx,"nuc",nuc(ibx,iset)
              print*,jset,ibx,"nrc",nrc(ibx,iset)
              print*,jset,ibx,"seg",seg(ibx,iset)
              print*,jset,ibx,"isgen",isgen(ibx,iset)
           enddo
        END DO
!
!       Update pointer to begin of set for the next atom
!
        IF (I.LT.NONTYP) THEN
           DO ISET = 1, KSETS
              IBLCK(ISET) = NBLCK(I,ISET) + 1
           END DO
        END IF
!
!     End of loop over unique atom types
      END DO
!
      RETURN
      END
!------------------------------------------------------------------------------------
!  /* Deck READ_BASIS */
      SUBROUTINE READ_BASIS (LUINFO,WORK,LWORK,Q,IQM,NBLCK,CNTBAS,      &
     &           JCO,BASREF,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,ISGEN,       &
     &           KANG,ISET,KAOVEC,KPRIM,GENERATE)
!*****************************************************************************
!
!       Read an expansion set for a given center
!       Note : NBLCK is the total number of blocks for this type of basis set that gets updated
!              NBLOCK is the number of blocks for this type and atom that is used locally
!
!       Most of the arguments form parts of the larger arrays described in CENINP,
!       the logical GENERATE tells the caller that a basis set is not defined explicity in input
!       and should be generated (e.g. the small component by the kinetic balance relation)
!
!*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      DIMENSION WORK(LWORK)
      LOGICAL   SEG,NOORBTS
      DIMENSION JCO(KANG),NUC(KAOVEC),NRC(KAOVEC),                      &
     &          SEG(KAOVEC),                                            &
     &          ALPHA(KPRIM,KAOVEC),                                    &
     &          CPRIM(KPRIM,KPRIM,KAOVEC),                              &
     &          CPRIMU(KPRIM,KPRIM,KAOVEC),ISGEN(KAOVEC)
      CHARACTER SPDCAR*1,BSET*5,BASFIL*20,BSKEYWORD*14
      CHARACTER*80 BASREF(10)
      PARAMETER (MAXBSFIL=5)
      CHARACTER*20 BSFIL(MAXBSFIL)
      LOGICAL   GENERATE, CNTBAS
      PARAMETER (DSM = 1.0D-30)

! cbirea.h : UNCONT
#include "molinp.h"
#include "nuclei.h"
#include "huckel.h"
#include "dcbgen.h"
#include "cbirea.h"

      GENERATE = .FALSE.
      CNTBAS = .FALSE.
      NHTYP  = 0
      IQM = 0
      NBLOCK = 0
      ISGEN(1:KAOVEC) = 0
!
!     Initialize reference strings for this atom type.
!
      BASREF(1) = ' Reference not found in input'
      DO IREF = 2, 10
         BASREF(IREF) = ' Not initialized'
      END DO
!
!     Read basis set info from file MOLECULE.INP
!     ==========================================
!
!     The input is handled by different routines, depending on the
!     keywords that are given. The first gives the type of the set
!     of functions that is read in. This shoold correspond to the set
!     parameter ISET that is given in the call to BASINP. Valid
!     names are :
!
!     ISET     BSKEYWORD   Type of functions
!       1      LARGE       Large component basis functions
!       2      SMALL       Small component basis functions
!       3      FTSET       Primary set for density fitting (total density matrix)
!       3      FTLRG       Primary set for density fitting (LL density matrix)
!       4      FTSML       Secundary set for density fitting (SS density matrix)
!     (note that FTLRG is an alias for FTSET at the moment)
!
!     The next keyword gives the type of basis set. Valid options are :
!
!     IBSFLAG  BSKEYWORD   Type of functions
!       1      (blank)     Read from input file in standard format
!       1      INTGRL      Read from input file in standard format
!       1      EXPLICIT    Read from input file in standard format
!       2      MOLFBAS     Read from specified basisfile in MOLFDIR format
!       3      BASIS       Read from specified basisfile in DALTON  format
!       4      GEOM        Read from input in even tempered format (alpha beta N IQM JCO())
!       4      EVENTEMP    Read from input in even tempered format (alpha beta N IQN JCO())
!       4      WELLTEMP    Read from input in well tempered format (alpha beta gamma delta N IQM JCO())
!       4      FAMILY      Read from input in family basis format (N IQM JCO())
!       4      DUALFAMILY  Read from input in dual family basis format (N IQM JCO())
!       5      MULTIBASIS  Read from multiple basis sets in DALTON format (used to augment sets)
!       6      KINBAL      Generate the small component basis by kinetic balance
!
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(A5,1X,A14)') BSET,BSKEYWORD
!
      SELECT CASE (BSET)
      CASE ("LARGE")
          JSET = 1
      CASE ("SMALL")
          JSET = 2
      CASE ("FTSET")
          JSET = 3
      CASE ("FTLRG")
          JSET = 3
      CASE ("FTSML")
          JSET = 4
      CASE DEFAULT
          JSET = 0
      END SELECT
!
!     Check the correct order of the input. We do supply defaults,
!     but want to detect trivial mistakes before doing so.
!
      IF (JSET.EQ.0.AND.ISET.EQ.1) THEN
         IF (.NOT.BASIS) THEN
!           With no general default one should always specify 
!           the large component basis
            WRITE(LUPRI,'(A,A,I5/A)')'Expected LARGE keyword, read',    &
     &      ' line no. :',NMLINE,MLINE(NMLINE)
           CALL QUIT('Error reading keyword in BASINP')
         END IF
      END IF
!
      IF (JSET.GT.0) THEN
!        Correct format, but it may not be the set we're looking for..
         IF (JSET.LT.ISET) THEN
            WRITE(LUPRI,'(A,A)')'Wrong order of keywords, specify',     &
     &      ' them in the order LARGE SMALL FTSET/FTLRG FTSML'
            CALL QUIT ('BASINP : Keywords in wrong order')
         ELSEIF (JSET.GT.ISET) THEN
!           We will handle this in a later call, make the select treat this as wrong line
            BSKEYWORD = "LATER"
         ENDIF
      END IF
!
!     Initialize basis set flag
!
      IBLANK = INDEX(BSKEYWORD," ")
      IF (IBLANK.NE.0) BSKEYWORD = BSKEYWORD(1:IBLANK-1)
!
      print*,bskeyword
      SELECT CASE (BSKEYWORD)
      CASE ("INTGRL")
           IBSFLAG = 1
           READ (MLINE(NMLINE),'(12X,12I5)') IQM,(JCO(K),K=1,IQM)
      CASE ("EXPLICIT")
           IBSFLAG = 1
           READ (MLINE(NMLINE),'(14X,12I5)') IQM,(JCO(K),K=1,IQM)
      CASE ("MOLFBAS")
           IBSFLAG = 2
           READ (MLINE(NMLINE),'(14X,A20)') BASFIL
      CASE ("BASIS")
           IBSFLAG = 3
           READ (MLINE(NMLINE),'(12X,A20)') BASFIL
      CASE ("WELLTEMP")
           IBSFLAG = 4
           READ (MLINE(NMLINE)(16:),*,ERR=9000)                         &
     &          (FAMPAR(K),K=1,4),NFAMEXP(1),IQM,(JCO(K),K=1,IQM)
           CALL GENFAMEXP(LUINFO,0)
      CASE ("EVENTEMP")
           IBSFLAG = 4
           READ (MLINE(NMLINE)(16:),*,ERR=9000)                         &
     &          (FAMPAR(K),K=1,2),NFAMEXP(1),IQM,(JCO(K),K=1,IQM)
           CALL GENFAMEXP(LUINFO,0)
      CASE ("GEOM")
           IBSFLAG = 4
           READ (MLINE(NMLINE)(11:),*,ERR=9000)                         &
     &          (FAMPAR(K),K=1,2),NFAMEXP(1),IQM,(JCO(K),K=1,IQM)
           FAMPAR(3) = 0.0D0
           FAMPAR(4) = 0.0D0
           CALL GENFAMEXP(LUINFO,0)
      CASE ("FAMILY")
           IBSFLAG = 4
           READ (MLINE(NMLINE),'(12X,12I5)',ERR=9000)                   &
     &          NFAMEXP(1),IQM,(JCO(K),K=1,IQM)
           CALL GENFAMEXP(LUINFO,1)
      CASE ("DUALFAMILY")
           IBSFLAG = 4
           READ (MLINE(NMLINE),'(16X,12I5)',ERR=9000)                   &
     &          NFAMEXP(1),NFAMEXP(2),IQM,(JCO(K),K=1,IQM)
           CALL GENFAMEXP(LUINFO,2)
      CASE ("MULTIBASIS")
           IBSFLAG = 5
           READ (MLINE(NMLINE)(18:),*,ERR=9000)                         &
     &          NFILES,(BSFIL(K),K=1,NFILES)
           IF ( NFILES .EQ. 0 ) THEN
              WRITE (LUPRI,*)'SPECIFY THE NUMBER OF BASIS SET FILES'
              GOTO 9000
           END IF
      CASE ("KINBAL")
           IF (ISET.EQ.2) THEN
              IBSFLAG = 6
           ELSE
              CALL QUIT ('KINBAL only defined for SC basis')
           END IF
      CASE DEFAULT
!         The default depends on the set that we are reading.
          SELECT CASE (ISET)
          CASE (1)
             IF (JSET.EQ.0) THEN
!               We have a default basis, the user omitted the line. Step back.
                IBSFLAG = 3
                NMLINE = NMLINE - 1
             ELSE
!               Assume it is like 'LARGE   2    1....'
                IBSFLAG = 1
                READ (MLINE(NMLINE),'(5X,12I5)',ERR=9000)               &
     &               IQM,(JCO(K),K=1,IQM)
             END IF
          CASE (2)
!            Assume kinetic balance, the user has omittted this line
             IBSFLAG = 6
             NMLINE = NMLINE - 1
          CASE (3,4)
!            Assume that no fit set was given, the user has omittted this line
             IBSFLAG = 0
             NMLINE = NMLINE - 1
          END SELECT
      END SELECT
!
!     We know what should be read, call the appropriate reader.
!
      SELECT CASE (IBSFLAG)
      CASE (0)
!
!        No basis set given, no need to read.
!        We have not set IQM yet in this case, do this here.
         IQM = 0
!
      CASE (1)
!
!         Read explicitly given basis set.
!
          CALL GTOINP(LUINFO,IQM,JCO,NUC,NRC,SEG,                       &
     &         ALPHA,CPRIM,CPRIMU,ISGEN,NBLOCK,KAOVEC,KPRIM)
          NBLCK = NBLCK + NBLOCK
          BASREF(1) = "Basis set typed explicitly in input file "
      CASE (2)
!
!         Molfdir basis set
!
          KBFPR = 1
          KBFCO = KBFPR + KPRIM
          KLAST = KBFCO + KPRIM*KPRIM
          IF (KLAST.GT.LWORK) CALL STOPIT('BASINP','BFGINP',KLAST,LWORK)
          CALL BFGINP(I,ISET,WORK(KBFPR),WORK(KBFCO),                   &
     &        IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,                   &
     &        NBLCK,ISGEN,KATOM,KANG,KAOVEC,KPRIM,BASFIL,               &
     &        BASREF)
      CASE (3)
!
!         Dalton basis set
!
          CALL BASLIB(IQM,JCO,NUC,NRC,SEG,                               &
     &         ALPHA,CPRIM,CPRIMU,                                       &
     &         NBLOCK,KAOVEC,KPRIM,Q,Q,DSM,                              &
     &         UNCONT,BASFIL,BASREF,IPREAD)
          NBLCK = NBLCK + NBLOCK
      CASE (4)
!
!         Generate well-temperered/family basis set
!
          CALL FAMBAS(LUINFO,IQM,JCO,NUC,NRC,SEG,                       &
     &         ALPHA,CPRIM,CPRIMU,ISGEN,                                &
     &         NBLOCK,KAOVEC,KPRIM)
          NBLCK = NBLCK + NBLOCK
          BASREF(1) =                                                   &
     &     "Well-tempered basis set typed explicitly in input file "
!
      CASE (5)
!
!         Read more than one Dalton type basis set
!
          CALL MULTI_BASIS(Q,IQM,JCO,NUC,NRC,SEG,                       &
     &         NFILES,BSFILES,ALPHA,CPRIM,CPRIMU,                       &
     &         NBLOCK,ISGEN,KAOVEC,KANG,KPRIM,BASREF)

      CASE (6)
!
!        Kinetic balance is special as it requires information
!        about the other basis that we do not have at this level.
!        Do this at a higher level.
!
         GENERATE = .TRUE.
      END SELECT
!
!     Check whether the basis set does not exceed the limit in angular momentum
!
      NOORBTS    = IQM.EQ.0
      IF (.NOT.NOORBTS) THEN
          NHTYP     = MAX(NHTYP,IQM)
          IF(NHTYP.GT.MXQN) GOTO 5000
      END IF
!
!     Check whether basis set is contracted
!     =====================================
!
      DO JBLCK = 1, NBLOCK
         IF ((NUC(JBLCK).NE.NRC(JBLCK)).AND.                            &
     &       (NRC(JBLCK).GT.0)) CNTBAS = .TRUE.
      END DO
!
      RETURN
!
!     Error messages:
!
 5000 CONTINUE
        WRITE (LUPRI,'(6X,A,I3,3A/9X,2A/9X,2(A,I3),A)')                 &
     &     '*  Input specifies highest orbital of atomic type ',        &
     &     I,' as "',SPDCAR(NHTYP-1),'".',                              &
     &     ' Highest allowed orbital in this version: ',                &
     &     SPDCAR(MXQN-1),                                              &
     &     ' Increase MXQN from',MXQN,' to',NHTYP/2,' and recompile.'
        CALL QUIT('Too high angular specified in input.')
 9000 CONTINUE
      WRITE(LUPRI,'(A,/,A,I3,/,A)')                                     &
     &     '*** ERROR in READ_BASIS ***',                               &
     &     'Error reading line ',NMLINE,                                &
     &     MLINE(NMLINE)
         CALL QUIT('*** ERROR in READIN ***')
      END
!------------------------------------------------------------------------------------
!  /* Deck MULTI_BASIS */
      SUBROUTINE MULTI_BASIS(Q,IQM,JCO,NUC,NRC,SEG,                     &
     &  NFILES,BSFIL,ALPHA,CPRIM,CPRIMU,                                &
     &  NBLOCK,ISGEN,KAOVEC,KANG,KPRIM,BASREF)
!*****************************************************************************
!
!       Driver to read a number of expansion sets. 
!       Original version by Andre Gomes, isolated and reorganized by L. Visscher.
!
!*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "maxaqn.h"
      PARAMETER (DSM = 1.0D-30)
#include "cbirea.h"
      LOGICAL SEG,NOORBTS,CNTBAS
      DIMENSION JCO(KANG),NUC(KAOVEC),NRC(KAOVEC),                      &
     &          SEG(KAOVEC),                                            &
     &          ALPHA(KPRIM,KAOVEC),                                    &
     &          CPRIM(KPRIM,KPRIM,KAOVEC),                              &
     &          CPRIMU(KPRIM,KPRIM,KAOVEC),ISGEN(KAOVEC)
      CHARACTER SPDCAR*1,BSET*5,BASFIL*20,BSKEYWORD*14
      CHARACTER*80 BASREF(10)
      CHARACTER*20 BSFIL(NFILES)
      LOGICAL GENERATE
#include "molinp.h"
#include "nuclei.h"
#include "huckel.h"
#include "dcbgen.h"
!     Scratch arrays needed for MULTIBASIS option
!     Should changed into allocatable arrays (either using fortran 90 allocate or
!     from work).
!
      PARAMETER (MAXBSFIL=5,KNG=MXAQN,KBLK=500,KPRM=30)
!
      INTEGER IQMT
      DIMENSION JCOT(KNG,MAXBSFIL)
      DIMENSION IBLKAGNMOM(KBLK)
      DIMENSION NUCT(KBLK)
      DIMENSION NRCT(KBLK)
      DIMENSION SEGT(KBLK)
      DIMENSION ALPHAT(KPRM,KBLK)
      DIMENSION CPRIMT(KPRM,KPRM,KBLK)
      DIMENSION CPRIMUT(KPRM,KPRM,KBLK)
      CHARACTER*80 BASREFT(10,MAXBSFIL)
      LOGICAL SEGT
!
! reading of basis sets from a number of files, including when the files
! do not necessarily have the same maximum and/or minimum orbital 
! angular momenta (last modif. 11/09/2003 by gomes) 
!
! we will move through the basis set files using the dummy variable iblock1
! 
       IBLOCK1 = 1
       MAXIQM  = 0
!
! we make sure we don't have any junk on the temporary variables 
!
       CALL IZERO(JCOT,KNG*NFILES)
       CALL IZERO(NUCT,KBLK)
       CALL IZERO(NRCT,KBLK)
!
       BASREF(1) = "Basis set read from multiple files"
!
       DO 463 M=1,NFILES
          KAOVEC1 = KAOVEC + 1 - IBLOCK1 
!
          BASREF(3*(M-1)+2) = " " 
          BASREF(3*(M-1)+3) = " ref. in file "//BSFIL(M) 
! 
          CALL BASLIB(IQMT,JCOT(1,M),                                   &
     &                 NUCT(IBLOCK1),                                   &
     &                 NRCT(IBLOCK1),SEGT(IBLOCK1),                     &
     &                 ALPHAT(1,IBLOCK1),                               &
     &                 CPRIMT(1,1,IBLOCK1),                             &
     &                 CPRIMUT(1,1,IBLOCK1),                            &
     &                 NBLOCK,KAOVEC1,KPRIM,Q,Q,1,DSM,                  &
     &                 BSFIL(M),BASREFT(1,M),IPREAD)
          BASREF(3*(M-1)+4) = BASREFT(1,M)
!
! to help in the reordering, we create a table mapping the orbital
! angular momentum of a block to the number of the block in the
! temp storage 
! 
          IDBLK = IBLOCK1
          IDANG = 1
          DO 464 IRR=1,IQMT
             DO 465 IRS=1,JCOT(IRR,M)
                IBLKAGNMOM(IDBLK) = IDANG
                IDBLK = IDBLK + 1
  465        CONTINUE 
             IDANG = IDANG + 1
  464     CONTINUE 
!
          IF ( IQMT .GT. MAXIQM ) THEN
             MAXIQM = IQMT
          END IF
          IBLOCK1 = IBLOCK1 + NBLOCK
  463  CONTINUE
!
! after reading all files, we reorder the data wrt orbital angular 
! momentum. in this reordering variable kt is the "final" index of a given block
!
       IQM = MAXIQM
!radovan: iblock is not set
       KT = IBLOCK
       DO 466 IQ=1,MAXIQM
       DO 467 IR=1,IBLOCK1 - 1
!
! in filling up the permanent storage with data, care should be taken to skip 
! zero-sized blocks that arise when we have only parts of a basis
! set in one file, as when the polarization functions are in one file and
! the scf set is in another...
!
            IF ( IBLKAGNMOM(IR) .EQ. IQ ) THEN
              IF ( NUCT(IR) .GT. 0 ) THEN

                   JCO(IQ)   = JCO(IQ)+ 1
                   NUC(KT)   = NUC(KT) + NUCT(IR)
                   NRC(KT)   = NRC(KT) + NRCT(IR)
                   SEG(KT)   = SEGT(IR)
                   ISGEN(KT) = 0

                   DO 468 IX=1, NUCT(IR)
                      ALPHA(IX,KT)      =  ALPHAT(IX,IR)
                      DO 468 JX=1, NRCT(IR) 
                       CPRIM(IX,JX,KT)  = CPRIMT(IX,JX,IR)
                       CPRIMU(IX,JX,KT) = CPRIMUT(IX,JX,IR)
  468              CONTINUE
                   KT = KT + 1
                END IF 
             END IF
  467     CONTINUE
  466  CONTINUE
!
! finally, we update nblck with the total number of non-zero sized
! blocks read, subtracting one only because on the loops above kt was 
! incremented after the last block was put in place ...
!
      NBLCK = KT - 1 
!
! ... and we are done
!
      RETURN
      END
