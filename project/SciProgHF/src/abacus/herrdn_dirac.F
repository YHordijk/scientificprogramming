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

C FILE : herrdn_dirac.F
C
C Purpose: special ".mol" input routines for DIRAC
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck bfginp */
      SUBROUTINE BFGINP(IAT,IBLCK,ICP,ALPIN,CPRIN,
     &              IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &              NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,FILNAM,BASREF)
C***********************************************************************
C
C      Read contracted basis set from a MOLFDIR type basis set file
C      (defined in FILNAM). IC tells us which component to read.
C      ICP = 1 : Read large component
C      ICP = 2 : Read small component
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "infpar.h"
#include "molinp.h"
C
      PARAMETER (DSM=1.0D-30)
C
      LOGICAL   SEG,TOBE,SEGIJ
      CHARACTER FILNAM*20,TEXT*80,DATE1*10,TIME1*8,FORMATE*10,FORMATC*10
      CHARACTER STRING*200, SZ(2)*5
      DIMENSION IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),NBLCK(KATOM,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2),ISGEN(KBLOCK),
     &          ALPIN(KPRIM),CPRIN(KPRIM*KPRIM)
      CHARACTER*80 BASREF(10,KATOM)
C     variable for slash in paths (/ or \\)
      CHARACTER SLASH
C     Check if we work on Linux or on Windows
      LOGICAL WIN
#include "dcbgen.h"
#include "dcbham.h"
#include "cbirea.h"
#include "ccom.h"
#include "gencon.h"
C
      DATA SZ / 'LARGE', 'SMALL' /
C
#include "ibtfun.h"
C
      LUBFIL = 55
      NTYP   = 0
      NBLOCK = IBLCK
#if defined (SYS_WINDOWS)
      WIN = .TRUE.
#else
      WIN = .FALSE.
#endif
C
C     We now deal with general contracted basis sets:
C     Use Unrestricted kinetic balance.
C     ===============================================
C
C     URKBAL = .TRUE.
C     This doesn't work because it is reset in HAMINP
C
C     Look if basis set file exists and, if so, read in data
C     ======================================================
C
C        ... this is the master, slaves will read from lines
C        read into MLINE(*) by master
C        
C      Linux:
C      The environment variable is a : separated string, e.g.:
C      /programs/Dirac/basis:/programs/Dalton/basis:.
C
C      Windows:
C      The environment variable is a ; separated string, e.g.:
C      C:\programs\Dirac\basis;C:\programs\Dalton\basis;.
C     
C      We loop over these directories in search for the basis set and
c      element given.
C
       ISTART = 1
C        
 10    CONTINUE
C        
C      Find next colon (: or ;) in BASDIR
       IF ( .NOT. WIN ) THEN
          ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ':' )
       ELSE
          ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ';' )
       END IF
C        
C      In there are no colons (: or ;), we are the last directory in BASDIR
       IF ( ICOLON .EQ. ISTART - 1 ) THEN
C           
C           Find end of string
C           
         ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:), ' ' )
C           
C           If end of string = start of string there are no more
C           directories, so we gracefully die...
C           
         IF ( ICOLON .EQ. ISTART )
     &        GOTO 5025
       END IF
       IEND = ICOLON - 1
       STRING = BASDIR(ISTART:IEND)
       ISTART = ICOLON + 1
C        
C     Add trailing slash if not present.
C     / for Linux and \\ for Windows
      IF ( .NOT. WIN ) THEN
          SLASH = '/'
      ELSE
          SLASH = '\\'
      END IF
C        
       IEND = INDEX ( STRING(1:), ' ') - 1
       IF ( STRING(IEND:IEND) .NE. SLASH ) THEN
          IEND = IEND + 1
          STRING(IEND:IEND) = SLASH
       END IF
       STRING = STRING(1:IEND)//FILNAM
       IEND = INDEX(STRING(1:),' ') - 1
CLV : Reactivate when debugging, otherwise see below.
CLV      WRITE(LUPRI,'(3A/,A1,A,A1)')
CLV  &        'BFGINP: Trying basis set file for ',
CLV  &        SZ(ICP),' component: ',
CLV  &        '"',STRING(1:IEND),'"'
C        
C        Inquire if "/path/basis-set" exists.
C        
       INQUIRE (FILE = STRING, EXIST = TOBE)
       IF (.NOT.TOBE) GOTO 10
       OPEN (LUBFIL,FILE=STRING,FORM='FORMATTED',STATUS='OLD')
       REWIND (LUBFIL)
      
C     Read info line of basis set file, copy this to BASREF array
C     ===========================================================
C
      CALL RDLINE(LUBFIL)
      READ (MLINE(NMLINE),'(A80)',END=5020) TEXT
C
      IF (ICP.EQ.1) THEN
C
C        We report the name of the file and the text line as both may give useful information
C
         BASREF(1,IAT) = "Basis set read from MOLFDIR file : "//FILNAM
         BASREF(2,IAT) = TEXT
      ELSE
C
C        Usually the small component comes from the same source as the LC, if not add additional info
C
         IF (BASREF(2,IAT).NE.TEXT) THEN
            BASREF(3,IAT) = 
     &   " Small component basis set read from MOLFDIR file : "//FILNAM
            BASREF(4,IAT) = TEXT
         ENDIF
      ENDIF
C
C     Read technical info line of basis set file
C     ==========================================
C
      CALL RDLINE(LUBFIL)
      READ (MLINE(NMLINE),'(A10,A8,2(A10,I4),F6.1,F20.5)',END=5020) 
     &      DATE1,TIME1,FORMATE,NE,FORMATC,NC,DCHRG,DEXPO
C
C     Read component, l-type and number of functions:
C     ICOMP  = 1 (large) or 2 (small)
C     ITYP   = l + 1
C     NPRIM  = # of primitives
C     NCONTR = # of contracted functions
C     NUNC   = # number of uncontracted functions
C     ===============================================
C
  301 CONTINUE
      CALL RDLINE(LUBFIL)
      READ (MLINE(NMLINE),'(6I4)',END=390,ERR=390)
     &      ICOMP,ITYP,NPRIM,NCONTR,IDUM,NUNC
      IF (ICOMP.EQ.0) GOTO 390 
C
C     Skip large component lines which we don't need in the small
C     component case
C     ===========================================================
C
      IF (ICOMP.NE.ICP) THEN
         IF (ICP.EQ.2) THEN
            NSTEP = NPRIM/NE + (NPRIM * NCONTR)/NC
            IF (MOD(NPRIM,NE).NE.0) NSTEP = NSTEP + 1
            IF (MOD(NCONTR*NPRIM,NC).NE.0) NSTEP = NSTEP + 1
            DO 310 INIT = 1, NSTEP
               READ (LUBFIL,*,END=5021)
  310       CONTINUE
            GOTO 301
         ELSE
            GOTO 390
         ENDIF
      ENDIF
      IQM(IAT,ICP) = ITYP
C
C     Check dimensions of arrays
C     ==========================
C
      IF (NPRIM.GT.KPRIM) GOTO 5040
      IF (NCONTR.GT.KPRIM) GOTO 5060
C
C     Read in primitives and contraction coefficients
C
C     Number of exponents and contraction coefficients on a line is 4
C     ===============================================================
C
      DO I = 1, NPRIM, 4
         CALL RDLINE(LUBFIL)
         READ (MLINE(NMLINE),FORMATE,ERR=5020,END=5020)
     &         (ALPIN(IND),IND=I,MIN(I+3,NPRIM))
      ENDDO
      DO I = 1, NPRIM*NCONTR, 4
         CALL RDLINE(LUBFIL)
         READ (MLINE(NMLINE),FORMATC,ERR=5020,END=5020)
     &         (CPRIN(IND),IND=I,MIN(I+3,NPRIM*NCONTR))
      ENDDO
C
C     We cannot split the primitives in blocks. However, we can
C     split up the contracted functions making that the leading
C     dimension. The number of contracted functions in a block
C     depends on the # of angular components and is set to a max.
C     of MXCON_PARAM (gencon.h)
C     ==========================================================
C
      MXDUM = MXCON_PARAM - MOD(MXCON_PARAM,(NPRIM*KHK(ITYP)))
      MXCONTR = MAX((MXDUM/(NPRIM*KHK(ITYP))),1)
      JCO(ITYP,IAT,ICP) = 1 + (NCONTR-1)/MXCONTR
C
C     Put the entire list of ITYP primitives and contracted functions
C     into KPRIM long blocks
C     ===============================================================
C
      ICONS = 1
      ICONE = 0
      DO 501 JBLCK = 1, JCO(ITYP,IAT,ICP)
         IF (JBLCK*MXCONTR.LE.NCONTR) THEN
            ICONE = ICONE + MXCONTR
         ELSE
            ICONE = NCONTR
         ENDIF
         IF(NBLOCK.GT.KBLOCK) GOTO 5000
         CALL DZERO(ALPHA (1,  NBLOCK,ICP),KPRIM)
         CALL DZERO(CPRIMU(1,1,NBLOCK,ICP),KPRIM*KPRIM)
         CALL DZERO(CPRIM (1,1,NBLOCK,ICP),KPRIM*KPRIM)
         IPRI = 0
         DO 502 JPRI = 1, NPRIM
            IPRI = IPRI + 1
            ALPHA(IPRI,NBLOCK,ICP) = ALPIN(JPRI)
            ICON = 0
            DO 503 JCON = ICONS, ICONE
               ICON = ICON + 1
               CPRIMU(IPRI,ICON,NBLOCK,ICP) = CPRIN(NPRIM*(JCON-1)+JPRI)
  503       CONTINUE
  502    CONTINUE
         NUCIJ = NPRIM
         NRCIJ = ICONE - ICONS + 1
         SEGIJ = NRCIJ .GE. 0
C
C        Remove exponents with all zero coefficients from block
C        This happens in SC in case of different primitives from
C        l+1 and l-1. It also happens with uncontracted functions
C
  600    CONTINUE
         DO 601 JPRI = 1, NUCIJ
            NZ = 0
            DO 602 JCON = 1, NRCIJ
               IF (ABS(CPRIMU(JPRI,JCON,NBLOCK,ICP)).GT.DSM) NZ = NZ + 1
  602       CONTINUE
            IF (NZ.EQ.0) THEN
               DO 603 IPRI = JPRI+1, NUCIJ
                  ALPHA(IPRI-1,NBLOCK,ICP) = ALPHA(IPRI,NBLOCK,ICP)
                  DO 604 JCON = 1, NRCIJ
                     CPRIMU(IPRI-1,JCON,NBLOCK,ICP) =
     &                      CPRIMU(IPRI,JCON,NBLOCK,ICP)
  604             CONTINUE
  603          CONTINUE
               ALPHA(NUCIJ,NBLOCK,ICP) = 0.0D0
               DO 605 JCON = 1, NRCIJ
                  CPRIMU(NUCIJ,JCON,NBLOCK,ICP) = 0.0D0
  605          CONTINUE
               NUCIJ = NUCIJ - 1
               GOTO 600
            ENDIF
  601    CONTINUE
C
C        Identify segmented contractions
C
         CALL SEGORB(SEGIJ,NUCIJ,NRCIJ,CPRIMU(1,1,NBLOCK,ICP),KPRIM,DSM)
         SEG(NBLOCK,ICP) = SEGIJ
C
C        Reorder primitive orbitals
C
         if (.not. ECPCALC) then
         CALL PRIORD(ALPHA(1,NBLOCK,ICP),CPRIMU(1,1,NBLOCK,ICP),NUCIJ,
     &               NRCIJ,SEG(NBLOCK,ICP),KPRIM,DSM)
         endif
C
C        Normalize orbitals
C
         CALL NRMORB(ITYP,NRCIJ,NUCIJ,ALPHA(1,NBLOCK,ICP),
     &               CPRIM(1,1,NBLOCK,ICP),
     &               CPRIMU(1,1,NBLOCK,ICP),KPRIM,NBLOCK)
C
         NUC(NBLOCK,ICP) = NUCIJ
         NRC(NBLOCK,ICP) = NRCIJ
C
         ISGEN (NBLOCK) = 0
C
         IF (INPTST)
     &      WRITE(LUPRI,'(2A,6I5)') SZ(ICP),' MOLFBAS Block:',
     &      NBLOCK,NUC(NBLOCK,ICP),NRC(NBLOCK,ICP)
         NBLOCK = NBLOCK + 1
         ICONS = ICONE + 1
  501 CONTINUE
      GOTO 301
  390 CONTINUE
      NBLCK(IAT,ICP) = NBLOCK - IBLCK
  391 CONTINUE
      CLOSE (LUBFIL)
      RETURN
 5000 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6/A)')
     *  ' >>> BFGINP error, number of AO-blocks            ',NBLOCK,
     *  '                   max. number (KBLOCK in READIN) ',KBLOCK,
     *  '                ** Increase WRKMEM and try again!'
        CALL QUIT('Too many AO-blocks')
 5020 CONTINUE
        WRITE (LUPRI,'(/A/A/2A)')
     &     'Error reading the following line :',MLINE(NMLINE),
     &     'in MOLFDIR type basis set file : ',FILNAM
        CALL QUIT('Error reading contents of MOLFDIR basis set file.')
 5021 CONTINUE
        WRITE (LUPRI,'(/A/A/2A)')
     &     'Error reading the line following this line :',MLINE(NMLINE),
     &     'in MOLFDIR type basis set file : ',FILNAM
        CALL QUIT('Error reading contents of MOLFDIR basis set file.')
 5025 CONTINUE
        WRITE (LUPRI,'(/2A)')
     &     'Cannot find MOLFDIR type basis set file : ',FILNAM
        CALL QUIT('Error finding MOLFDIR type basis set file.')
 5040 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BFGINP error, number of primitives           ',NPRIM,
     *  '                   current maximum number (MAXPRI)',KPRIM
        CALL QUIT('Too many primitives, increase MAXPRI in input.')
 5060 CONTINUE
        WRITE (LUPRI,'(/A,I6/A,I6)')
     *  ' >>> BFGINP error, number of contracted functions ',NCONTR,
     *  '                   current maximum number (MAXPRI)',KPRIM
        CALL QUIT('Too many contracted fun., increase MAXPRI in input.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck kinbal */
      SUBROUTINE KINBAL(IATOM,IBLL,IBLS,
     &     IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &     NBLCK,ISGEN,KATOM,KANG,KBLOCK,KPRIM,INPTST,
     &     CNTBAS)
C*****************************************************************************
C
C      Generate small component functions from large for
C      atom type IATOM using the kinetic balance relation
C
C           ISGEN - bit packed information on generation of small
C                   components using the kinetic balance relation
C             1 - downwards: small L+1 from large L
C             2 - upwards  : small L-1 from large L
C             0 - if not specified same as 3 (both down and up)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "ccom.h"
C
      LOGICAL SEG,INPTST,CNTBAS
      DIMENSION IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),NBLCK(KATOM,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2),ISGEN(KBLOCK)
#include "ibtfun.h"
      IDOWN(J) = IBTAND(1,J)
      IUP(J)   = IBTAND(1,IBTSHR(J,1))
C
      IBS = IBLS - 1
      IBU = IBLL + JCO(1,IATOM,1) - 1
      IBD = IBLL - 1
      JBS = IBS
      IQM(IATOM,2) = 0
C
C     Loop over possible range of (L+1) values for 
C     small components for atom type IATOM
C
      DO J = 1,(IQM(IATOM,1)+1)
        IBJ = IBS
C
C       Generate function from downwards
C
        IF(J.GT.1) THEN
          DO K = 1,JCO(J-1,IATOM,1)
            IBD = IBD + 1
            IF((IDOWN(ISGEN(IBD)).EQ.1).OR.(ISGEN(IBD).EQ.0))  THEN
              IF (CNTBAS) THEN
                JSBEG = 0
                JSEND = NUC(IBD,1)
                IBS = IBS + 1
                CALL GENSML(J,IBD,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                      SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
              ELSE
C               Collect the uncontracted small component
C               basis functions in chunks of about 50 basis
C               functions for higher efficiency in integral
C               evaluation.
                MXDUM = 50 - MOD(50, KHK(J))
                MAXS = MAX((MXDUM/KHK(J)),1)
                IBJCO=1+(NUC(IBD,1)-1)/MAXS
                JSEND = 0
                DO JBLCK=1, IBJCO
                  JSBEG = JSEND
                  IBS = IBS + 1
                  IF (JBLCK*MAXS.LE.NUC(IBD,1)) THEN
                    JSEND = JSEND + MAXS
                  ELSE
                    JSEND = NUC(IBD,1)
                  ENDIF
                  CALL GENSML(J,IBD,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                        SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
C
C       Generate function from upwards
C
        IF(J.LT.IQM(IATOM,1)) THEN
          DO K = 1,JCO(J+1,IATOM,1)
            IBU = IBU + 1
            IF((IUP(ISGEN(IBU)).EQ.1).OR.(ISGEN(IBU).EQ.0)) THEN
              IF (CNTBAS) THEN
                JSBEG = 0
                JSEND = NUC(IBU,1)
                IBS = IBS + 1
                CALL GENSML(J,IBU,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                      SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
              ELSE
C               See note above.
                MXDUM = 50 - MOD(50, KHK(J))
                MAXS = MAX((MXDUM/KHK(J)),1)
                IBJCO=1+(NUC(IBU,1)-1)/MAXS
                JSEND = 0
                DO JBLCK=1, IBJCO
                  JSBEG = JSEND
                  IBS = IBS + 1
                  IF (JBLCK*MAXS.LE.NUC(IBU,1)) THEN
                    JSEND = JSEND + MAXS
                  ELSE
                    JSEND = NUC(IBU,1)
                  ENDIF
                  CALL GENSML(J,IBU,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                        SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
                ENDDO
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        JCO(J,IATOM,2) = IBS - IBJ
        IF(JCO(J,IATOM,2).GT.0) IQM(IATOM,2) = J
      ENDDO
      NBLCK(IATOM,2) = IBS - JBS
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rkbbas */
      SUBROUTINE RKBBAS(IATOM,IBLL,IBLS,
     &     IQM,JCO,NUC,NRC,SEG,ALPHA,CPRIM,CPRIMU,
     &     NBLCK,KATOM,KANG,KBLOCK,KPRIM,INPTST,
     &     CNTBAS)
C*****************************************************************************
C
C      Generate small component functions atom type IATOM for use with
C      the modified Dirac Hamiltonian. 
C      This amounts to copying the large component basis information
C      into the small component one.
C
C      More elaborate schemes are necessary in the case of contracted basis sets
C
C       IBLL - index of first large component block of this atom type
C       IBLS - index of first small component block of this atom type
C       IQM             (KATOM)         - highest L-value
C       NBLCK           (KATOM)         - number of AO-blocks for a given
C                                         atomic type
C       JCO             (KANG,KATOM)    - number of AO-blocks for a given
C                                         atomic type and L-value
C       NUC             (KBLOCK)        - number of uncontracted functions
C                                         in a given AO-block
C       NRC             (KBLOCK)        - number of contracted functions
C                                         in a given AO-block
C       SEG             (KBLOCK)        - TRUE: segmented contraction
C       ALPHA           (KPRIM,KBLOCK)  - exponents
C       CPRIM           (KPRIM,KPRIM,KBLOCK)  - normalized contraction
C                                               coefficients
C       CPRIMU          (KPRIM,KPRIM,KBLOCK)  - contraction coefficients
C
C      Written by Trond Saue May 2 2012
C
C*******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "ccom.h"
C
      LOGICAL SEG,INPTST,CNTBAS
      DIMENSION IQM(KATOM,2),JCO(KANG,KATOM,2),
     &          NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),NBLCK(KATOM,2),
     &          CPRIMU(KPRIM,KPRIM,KBLOCK,2)
C
      IQM(IATOM,2)   = IQM(IATOM,1)
      NBLCK(IATOM,2) = NBLCK(IATOM,1)
      IBLOCK = IBLS
      DO I = 1, IQM(IATOM,2)
        JCO(I,IATOM,2) = JCO(I,IATOM,1)
        DO J = 1,JCO(I,IATOM,2)
          NUC(IBLOCK,2) = NUC(IBLOCK,1)
          NRC(IBLOCK,2) = NRC(IBLOCK,1)
          SEG(IBLOCK,2) = SEG(IBLOCK,1)
          ALPHA(1:NUC(IBLOCK,2),IBLOCK,2)
     &   =ALPHA(1:NUC(IBLOCK,1),IBLOCK,1)
          CPRIM(1:NUC(IBLOCK,2),1:NRC(IBLOCK,2),IBLOCK,2)
     &   =CPRIM(1:NUC(IBLOCK,1),1:NRC(IBLOCK,1),IBLOCK,1)
          CPRIMU(1:NUC(IBLOCK,2),1:NRC(IBLOCK,2),IBLOCK,2)
     &   =CPRIMU(1:NUC(IBLOCK,1),1:NRC(IBLOCK,1),IBLOCK,1)
          IBLOCK = IBLOCK + 1
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gensml */
      SUBROUTINE GENSML(IQM,IBL,IBS,IBJ,NUC,NRC,ALPHA,CPRIM,CPRIMU,
     &                  SEG,KBLOCK,KPRIM,INPTST,JSBEG,JSEND)
C*****************************************************************************
C
C     IBS - index of actual small component (SC) block
C     IBL - index of large component (LC) block from which it is generated
C     IBJ - offset for start of SC block with same IQM and atom type
C     IQM - (L+1)-value of actual SC block
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "infpar.h"
C
      PARAMETER (D1=1.0D0,SML = 1.0D-6)
C
      LOGICAL SEG,INPTST
      DIMENSION NUC(KBLOCK,2),NRC(KBLOCK,2),SEG(KBLOCK,2),
     &          ALPHA(KPRIM,KBLOCK,2),
     &          CPRIM(KPRIM,KPRIM,KBLOCK,2),CPRIMU(KPRIM,KPRIM,KBLOCK,2)
C
C
C     Small component basis set is always uncontracted,
C     hence NRC(IBS) = NUC(IBL) and SEG(IBS) = T
C
      SEG(IBS,2)   = .TRUE.
C
C     Remove small component basis functions already present.
C
      NUC_IBS = 0
      DO I = JSBEG + 1,JSEND
         ALPTEST = ALPHA(I,IBL,1)
         ALPINV  = D1/ALPTEST
C
C        Check previous shells:
C
         DO J = IBJ+1,IBS-1
            DO IJ = 1,NUC(J,2)
              TEST = ABS(ALPTEST-ALPHA(IJ,J,2))*ALPINV
              IF (TEST.LT.SML) THEN
C              WRITE(LUPRI,'(A,E24.16/8X,A,E24.16,A)') 
C     &           'WARNING from GENSML: Skipping ',ALPTEST,
C     &           'too close to ',ALPHA(IJ,J,2),' from previous shells'
                GO TO 100
              ENDIF
            END DO
         END DO
C
C        Check same shell (to eliminate duplication in
C        same shell as for example for Chlorine in
C        6-311G* basis set -- where some basis functions
C        are duplicated in order to have segmented contraction).
C        /HJAaJ D1c 2005
C     
         DO IJ = 1, NUC_IBS
            TEST = ABS(ALPTEST-ALPHA(IJ,IBS,2))*ALPINV
            IF (TEST.LT.SML) THEN
C              WRITE(LUPRI,'(A,E24.16/8X,A,E24.16,A)') 
C     &         'WARNING from GENSML: Skipping ',ALPTEST,
C     &         'too close to ',ALPHA(IJ,IBS,2),' from same shell'
              GO TO 100
            ENDIF
         END DO
C
         NUC_IBS = NUC_IBS + 1
         ALPHA(NUC_IBS,IBS,2) = ALPTEST
C
  100 CONTINUE
      END DO
      IF ((JSEND-JSBEG) .NE. NUC_IBS) THEN
         IF (.NOT.SLAVE) WRITE(LUPRI,'(A,I3,A,I3)')
     &   ' INFO: GENSML deleted',(JSEND-JSBEG) - NUC_IBS,
     &   ' small functions in a block out of',(JSEND-JSBEG)
      END IF
      IF (NUC_IBS .EQ. 0) THEN
C
C        We've deleted all functions from this block.
C        This might happen for (dual) family basis sets.
C        Decrement IBS and return.
C
         IBS = IBS - 1
         RETURN
      END IF
      NUC(IBS,2) = NUC_IBS
      NRC(IBS,2) = NUC_IBS
C
C     Contraction matrix for SMALL is just a unit matrix
C
      CALL DUNIT(CPRIMU(1,1,IBS,2),KPRIM)
      CALL NRMORB(IQM,NRC(IBS,2),NUC(IBS,2),ALPHA(1,IBS,2),
     &            CPRIM(1,1,IBS,2),CPRIMU(1,1,IBS,2),KPRIM,IBS)
C
C     Test print
C
      IF (INPTST) THEN
         WRITE(LUPRI,'(A,I4,A,I4,A)')
     &   'Contraction matrix (',NUC(IBS,2),'x',NUC(IBS,2),') :' 
         DO 400 L = 1, NUC(IBS,2)
            WRITE (LUPRI,'(E16.10)') ALPHA(L,IBS,2)
            WRITE (LUPRI,'(6(2X,F10.8))')
     $           (CPRIMU(L,M,IBS,2),M=1,NUC(IBS,2))
  400    CONTINUE
      END IF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck rkbrel */
      SUBROUTINE RKBREL(IRKB,IPRINT)
C******************************************************************************
C
C     Routine for Restricted Kinetic Balance (RKB) Relations.
C
C     Generates array IRKB running over shells
C       For large shells IRKB points to L+1 small partner
C       For small shells IRKB points to L-2 small partner
C       or to the L-1 large partner if the L-2 small partner does
C       not exist (pointer is negative in this case).
C
C     Two examples may make it easier to understand the code:
C     Consider the large component basis L1s1p1d:
C     a) the s and d exponents are different
C       The RKB array will be
C         IRKB(Ls   1) =   5 : pointing to Sp1
C         IRKB(Lp   2) =   7 : pointing to Sd
C         IRKB(Ld   3) =   8 : pointing to Sf
C         IRKB(Ss   4) =   0 : no (L-2) small partner, no (L-1) large partner 
C         IRKB(Sp1  5) =  -1 : no (L-2) small partner, pointing to Ls
C         IRKB(Sp2  6) =   0 : no (L-2) small partner, no (L-1) large partner 
C         IRKB(Sd   7) =   4 : pointing to (L-2) small partner Ss
C         IRKB(Sf   8) =   6 : pointing to (L-2) small partner Sp2
C
C     b) dual family basis; the s and exponents are identical
C
C         IRKB(Ls  1) =   5 : pointing to Sp
C         IRKB(Lp  2) =   6 : pointing to Sd
C         IRKB(Ld  3) =   7 : pointing to Sf
C         IRKB(Ss  4) =   0 : no (L-2) small partner, no (L-1) large partner 
C         IRKB(Sp  5) =  -1 : pointing to Ls
C         IRKB(Sd  6) =   4 : pointing to (L-2) small partner Ss
C         IRKB(Sf  7) =   5 : pointing to (L-2) small partner Sp
C
C     Consider next a 1s1p1d1f1g basis:
C     a) general case
C        IRKB(Ls   1) =   7 : pointing to Sp1, that is (L+1) small partner
C        IRKB(Lp   2) =   9 : pointing to Sd1, that is (L+1) small partner
C        IRKB(Ld   3) =  11 : pointing to Sf1, that is (L+1) small partner
C        IRKB(Lf   4) =  13 : pointing to Sg , that is (L+1) small partner
C        IRKB(Lg   5) =  14 : pointing to Sh , that is (L+1) small partner
C        IRKB(Ss   6) =   0 : no (L-2) small partner, no (L-1) large partner 
C        IRKB(Sp1  7) =  -1 : pointing to (L-1) large partner
C        IRKB(Sp2  8) =   0 : no (L-2) small partner, no (L-1) large partner 
C        IRKB(Sd1  9) =   6 : pointing to (L-2) small partner
C        IRKB(Sd2 10) =   0 : no (L-2) small partner, no (L-1) large partner 
C        IRKB(Sf1 11) =   8 : pointing to (L-2) small partner
C        IRKB(Sf2 12) =   0 : no (L-2) small partner, no (L-1) large partner 
C        IRKB(Sg  13) =  10 : pointing to (L-2) small partner
C        IRKB(Sh  14) =  12 : pointing to (L-2) small partner
C
C     b) dual family basis
C        IRKB(Ls  1) =   7  : pointing to Sp, that is (L+1) small partner
C        IRKB(Lp  2) =   8  : pointing to Sd, that is (L+1) small partner
C        IRKB(Ld  3) =   9  : pointing to Sf, that is (L+1) small partner
C        IRKB(Lf  4) =  10  : pointing to Sg, that is (L+1) small partner
C        IRKB(Lg  5) =  11  : pointing to Sh, that is (L+1) small partner
C        IRKB(Ss  6) =   0  : no (L-2) small partner, no (L-1) large partner 
C        IRKB(Sp  7) =  -1  : pointing to (L-1) large partner
C        IRKB(Sd  8) =   6  : pointing to (L-2) small partner
C        IRKB(Sf  9) =   7  : pointing to (L-2) small partner
C        IRKB(Sg 10) =   8  : pointing to (L-2) small partner
C        IRKB(Sh 11) =   9  : pointing to (L-2) small partner
C
C
C FIXME: In principle kinetic balance can be made for contracted
C FIXME: basis sets if the small component basis set exponents
C FIXME: are kinetically balanced. Currently, the code resets
C FIXME: ISPHTR if the small component basis set is explicitly given
C FIXME: (IBSFLAG .NE. 4, see INFO: above)./hjaaj+jth sep 2000
C
C     Written by Jon K. Laerdahl 05.12.1996
C     Major surgery by T.Saue May 5 1998
C     Changed definition of IRKB because of problems with contracted sets,
C     L. Visscher, Nov 16, 2000.
C
C******************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "shells.h"
#include "dcbbas.h"
#include "aovec.h"
      PARAMETER (SML = 1.0D-6)
#include "primit.h"
      DIMENSION IRKB(NPSHEL)
C
      IF (IPRINT.GE.10) WRITE(LUPRI,'(/1X,A)')
     &   'RKBREL: Constructing relations for scalar RKB transformation'
      CALL IZERO(IRKB,NPSHEL)
C
C     New code for generating RKBREL.
C     This should also work for contracted large component basis set.
C
C     Loop over large and small component shells.
C     For contracted basis sets, the same primitives are present in many
C     shells, and we skip shells with the same primitives.
C     A potential problem arises if you block a contracted basis set.
C
C
      iPrimTest = 0
C
      DO 100 ISHL = 1, KMAX
C
         NCENTI = NCENT(ISHL)
         NHKTI  = NHKT(ISHL)
C
         iPrimStart = JSTRT(ISHL) + 1
         iPrimEnd   = JSTRT(ISHL) + NUCO(ISHL)
C
C        Skip primitives we've already processed.
C
         IF (iPrimStart .EQ. iPrimTest) GO TO 100
         iPrimTest = iPrimStart
C
C
C        Loop over primitives in this shell.
C
         DO 200 IPRIM = iPrimStart, iPrimEnd
            PRIEXI = PRIEXP(IPRIM)
C
            JSH2M = 0
            jPrimTest = 0
C
C           Loop over small component shells
C           Find matching small component basis functions.
C
            DO 300 JSHL = NLRGSH+1, KMAX
               NCENTJ = NCENT(JSHL)
               NHKTJ  = NHKT(JSHL)
C
C              Must be on same atomic center.
C
               IF ( NCENTJ .NE. NCENTI ) GO TO 300
C
               jPrimStart = JSTRT(JSHL) + 1
               jPrimEnd   = JSTRT(JSHL) + NUCO(JSHL)
C
C              Skip primitives we've already processed.
C
               IF ( jPrimStart .EQ. jPrimTest) GO TO 300
               jPrimTest = jPrimStart
C
C              Loop over primitives.
C
               DO 400 JPRIM = jPrimStart, jPrimEnd
C
                  PRIEXJ = PRIEXP(JPRIM)
C
                  IF ( ABS ((PRIEXJ-PRIEXI)/PRIEXI) .LT. SML ) THEN
C
C
C                    We've got a match: same center, same exponent
C
C
C                    Match small component to L-2 small function....
C
                     IF(NHKTJ.EQ.NHKTI-2.AND.ISHL.GT.NLRGSH) THEN
C
C                       Check for two identical functions.
C                       This should never happen, since GENSML deletes
C                       functions it already knows.
C
                        IF(JSH2M.NE.0) THEN
                          WRITE(LUPRI,'(/A,2(/A,2I8,F20.10),2(/A,2I8))')
     &                    'RKBREL: Error, two identical S.C. functions',
     &                     ' ISHL,IPRIM,PRIEXI :',ISHL,IPRIM,PRIEXI,
     &                     ' JSHL,JPRIM,PRIEXJ :',JSHL,JPRIM,PRIEXJ,
     &                     ' NCENTI,NCENTJ     :',NCENTI,NCENTJ,
     &                     ' NHKTI,NHKTJ       :',NHKTI,NHKTJ
                          CALL QUIT
     &                   ('RKBREL: Error, two identical S.C. functions')
                        END IF
C
C                       Note that we use JSH2M and not JPRIM, since RKBLOW
C                       in dirone.F wants shell number, and not primitive
C                       number. This is not the same if L Component basis
C                       set is contracted.
C
                        JSH2M = JSHL + (JPRIM - jPrimStart)
                        IRKB(IPRIM) = JSH2M
                     ENDIF
C
C                    Match large component to L+1 small function and
C                    vice versa. The initialization of the SC is
C                    overwritten if an (L-1) SC does also exist. It
C                    serves to check whether this SC is used to balance
C                    by itself.
C
                     IF (NHKTJ.EQ.NHKTI+1.AND.ISHL.LE.NLRGSH) THEN
                        JSH1P = JSHL + (JPRIM - jPrimStart)
                        ISH1M = ISHL + (IPRIM - iPrimStart)
                        IRKB(IPRIM) = JSH1P
                        IRKB(JPRIM) = - ISH1M
                     ENDIF
                  ENDIF
 400           CONTINUE
 300        CONTINUE
 200     CONTINUE
 100  CONTINUE
C
      IF ( IPRINT .GE. 10 ) THEN
         CALL HEADER('RKBREL: RKB array',-1)
         DO I = 1,NPSHEL
            WRITE(6,'(A,I3,A,I3)') 'IRKB(',I,') = ',IRKB(I)
         END DO
      END IF
C
C
      RETURN
      END
C --- end of herrdn_dirac.F ---
