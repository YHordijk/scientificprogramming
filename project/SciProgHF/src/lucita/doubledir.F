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
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIRDIR(INSPC,IUTSPC,ISM,LUIN,LEVEL,E0,LUSCR,LUSCR2,
     &                  LU1)


* A vector is given in space INSPC . Calculate approximations
* to energies in larger space UTSPC
*
* Input vector is on file LUIN, but is assumed expanded to
* IUTSPC ( Suitable empty blocks have been inserted )
*
* LEVEL = 1 : Calculate norm of (H-E0)!0>
*                          !1> = (H0-E0)**(-1) (H-E0) !0>
*             Second order energy = <0!(H-E0) (H0-E0)**(-1) (H-E0)!0>
*             <1!1>,<0!1>
* The second order energy is only meaningful for spaces not included in !0>

*
* LEVEL = 2 : Calculate also <1!H!1>,
*
* Jeppe Olsen, August 1995 : Trying to get my miltibillion debt paid
*                            before going to Cambridge
*
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "strinp.inc"
#include "oper.inc"
#include "glbbas.inc"
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
C     COMMON/CECORE/ECORE,ECORE_ORIG
C     COMMON/CECORE/ECORE,ECORE_ORIG,ECORE_H
#include "cecore.inc"
*
      IDUM  = 0
      CALL MEMMAN(IDUM,IDUM,'MARK ',IDUM,'DIRDI')
*
      WRITE(6,'(20X,A)') '***********************************'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '*       Welcome to DIRDIR         *'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '*   vector free calculations      *'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '* Jeppe Olsen Version of oct.1995 *'
      WRITE(6,'(20X,A)') '*                                 *'
      WRITE(6,'(20X,A)') '***********************************'
      WRITE(6,'(20X,A)')
*
      WRITE(6,*) ' Input space ', INSPC
      WRITE(6,*) ' Output space ', IUTSPC
      WRITE(6,*) ' Level of calculation ', LEVEL
      WRITE(6,*)
*
      IATP = 1
      IBTP = 2
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
      NEL = NAEL + NBEL
*. Number of occupation classes in largest possible space ( compund space )
      IWAY = 1
      CALL OCCLS(1,NOCCLS,IOCCLS,NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),0,0)
C          OCCLS(IWAY,NOCCLS,IOCCLS,NEL,NGAS,IGSMIN,IGSMAX)
*. and then the occupation classes
      CALL MEMMAN(KLOCCLS,NGAS*NOCCLS,'ADDL  ',1,'KLOCCL')
      IWAY = 2
      CALL OCCLS(IWAY,NOCCLS,WORK(KLOCCLS),NEL,NGAS,
     &           IGSOCC(1,1),IGSOCC(1,2),0,0)
*. What classes of compound space are active in input and output space
      CALL MEMMAN(KLACCLI,NOCCLS,'ADDL  ',1,'ACCLI ')
      CALL MEMMAN(KLACCLO,NOCCLS,'ADDL  ',1,'ACCLO ')
*     OCCLS_IN_CI(NOCCLS,IOCCLS,ICISPC,NINCCLS,INCCLS)
      CALL OCCLS_IN_CI(NOCCLS,WORK(KLOCCLS),INSPC,NACTI,WORK(KLACCLI))
      CALL OCCLS_IN_CI(NOCCLS,WORK(KLOCCLS),IUTSPC,NACTO,WORK(KLACCLO))
*
*
*. Space for saving results : S11 : <1!1> in different classes
*.                            S01 : <0!1> in different classes
*.                            RES0 : norm of (H-E)!0> in different classes
*                             E2    : Second order energy contributions
*                                     from different classes
*
      CALL  MEMMAN(KLS11,NOCCLS,'ADDL  ',2,'KLS11 ')
      CALL  MEMMAN(KLS01,NOCCLS,'ADDL  ',2,'KLS01 ')
      CALL  MEMMAN(KLRES0,NOCCLS,'ADDL  ',2,'KLRES0')
      CALL  MEMMAN(KLE2 ,NOCCLS,'ADDL  ',2,'KLE2  ')
*. Calculate FI+FA if required
      IF(IPART.EQ.1) THEN
        WRITE(6,*) ' M-P Hamiltonian will be constructed '
        CALL COPVEC(WORK(KINT1O),WORK(KFI),NINT1)
        CALL FIFAM(WORK(KFI))
        CALL COPVEC(WORK(KFI),WORK(KFIO),NINT1)
        ECORE_H = 0.0D0
        IF(IUSE_PH.EQ.1) THEN
          CALL FI(WORK(KFI),ECORE_H,0)
C         ECORE = ECORE_ORIG + ECCP
        END IF
      END IF


        CALL DIRDIR1N(E0,ECORE,INSPC,IUTSPC,ISM,LUIN,
     &       WORK(KLS11),WORK(KLS01),WORK(KLRES0),WORK(KLE2),
     &       NOCCLS,WORK(KLOCCLS),NACTI,WORK(KLACCLI),
     &       NACTO,WORK(KLACCLO),LUSCR,LUSCR2,LEVEL,LU1)
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE DIRDIR1N(EREF,ECORE,INSPC,IIUTSPC,ISM,LUC,
     &             S11,S01,RES0,E2,
     &             NOCCLS,IOCCLS,NACTI,IACTI,
     &             NACTO,IACTO,LUSCR,LUSCR2,LEVEL,LU1)
*
* Part 1 of double direct CI approach
*
* First order correction to  the CI vector
* !1> = -(H0-E0) ** -1 ( H-E0 ) !
*
* is calculated in pieces and the information is processed to give
*
* S11 : Overlap matrix S11(IOCCLS) = <1(IOCCLS)!1(IOCCLS)>
* S01 : Overlap matrix S01(IOCCLS) = <0(IOCCLS)!1(IOCCLS)>
* RES0 : Norm of residual =  (H-E) !0> in different classes
* E2   : <0!(H-E0)(H0-E0)**-1(H-E0)!0> in different classes
*
*
* Jeppe Olsen, August 1995
*
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD, INPRDD
*
#include "mxpdim.inc"
#include "wrkspc.inc"
#include "cgas.inc"
#include "gasstr.inc"
#include "crun.inc"
#include "cicisp.inc"
#include "csm.inc"
#include "strbas.inc"
#include "stinf.inc"
#include "cstate.inc"
#include "orbinp.inc"
#include "lucinp.inc"
#include "glbbas.inc"
#include "strinp.inc"
#include "cprnt.inc"
*. Common block for communicating with sigma
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
#include "oper.inc"
C     COMMON/OPER/I12,IPERTOP,IAPR,MNRS1E,MXRS3E,IPART
      COMMON/CINTFO/I12S,I34S,I1234S,NINT1,NINT2,NBINT1,NBINT2
*.


*
      DIMENSION IOCCLS(NGAS,*)
      DIMENSION S11(NOCCLS),RES0(NOCCLS),E2(NOCCLS)
      DIMENSION S01(NOCCLS)
      IUTSPC = IIUTSPC

*
      NTEST =   20
*. LBLK should be defined outside for future use
      LBLK = -1


*
C     IF(NTEST.GE.3) THEN
        WRITE(6,*)
        WRITE(6,*) ' ================================= '
        WRITE(6,*) '        Welcome to DIRDIR1N '
        WRITE(6,*) ' ================================= '
        WRITE(6,*)
        WRITE(6,*) '     Input  CI space ', INSPC
        WRITE(6,*) '     Output CI space ', IUTSPC
        WRITE(6,*) '     Level parameter ', LEVEL
        WRITE(6,*)
        WRITE(6,*) ' Initial CI vector is in space ', INSPC
        WRITE(6,*)
        WRITE(6,*) ' Calculations to be performed : '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        IF(LEVEL.GE.1)
     &  WRITE(6,*) '     PT2 calculation in space .......', IUTSPC
        IF(LEVEL.GE.2)
     &  WRITE(6,*) '   + 1 CI iteration + PT3 in space ..', IUTSPC
        IF(LEVEL.EQ.3)
     &  WRITE(6,*) '   + PT4 in space ...................', IUTSPC
        IF(LEVEL.EQ.4)
     &  WRITE(6,*) '   + PT4 in space ...................', IUTSPC+1
        WRITE(6,*)
*
        IF(IC1DSC.EQ.0) THEN
          WRITE(6,*) ' First order correction is not stored '
        ELSE
          WRITE(6,*) ' First order correction stored on file ', LU1
        END IF
        WRITE(6,*) ' Scratch files ', LUSCR,LUSCR2
        WRITE(6,*) ' LU1 = ', LU1
*
C     END IF
*. Zero - before I forget
      ZERO = 0.0D0
      CALL SETVEC(S11,ZERO,NOCCLS)
      CALL SETVEC(S01,ZERO,NOCCLS)
      CALL SETVEC(RES0,ZERO,NOCCLS)
      CALL SETVEC(E2,ZERO,NOCCLS)
*
*. Transfer to /CANDS
      ICSPC = INSPC
      ISSPC = IUTSPC
*. zero and first order space
      ISPC0 = INSPC
      ISPC1 = IUTSPC
*
      ISSM = ISM
      ICSM = ISM
*
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'DIRDR1')
* The calculation of (H-E0) will be processed in blocks,
* allocate blocks for (H-E0)!0> and !0>
      WRITE(6,*) ' Size of BATCH ', LBLOCK
      CALL MEMMAN(KCB,LBLOCK,'ADDL  ',2,'KCB   ')
      CALL MEMMAN(KSB,LBLOCK,'ADDL  ',2,'KSB   ')
CM    CALL MEMMAN(KLB,LBLOCK,'ADDL  ',2,'KLB   ')
*
*. Structure and dimensions of output CI vector.
*
      IATP = 1
      IBTP = 2
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
C START
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
C END
*. a bit of local scratch
      WRITE(6,*) ' MXNTTS',MXNTTS
      NTTS = MXNTTS
      CALL MEMMAN(KLISCR,3*NTTS,'ADDL  ',1,'KLISCR')
      CALL MEMMAN(KLISC2,3*NTTS,'ADDL  ',1,'KLI2SC')
*
* Sigma information
*
      CALL MEMMAN(KLSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KLSBLTP,NSMST,'ADDL  ',2,'SBLTP ')
*
      CALL IAIBCM(IUTSPC,WORK(KLSIOIO))
      KSVST = 1
      WRITE(6,*) ' ISM = ',ISM
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLSBLTP),WORK(KSVST))
*
C     CALL MEMMAN(KSNOOS,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'SNOOS ')
C     CALL MEMMAN(KSIOOS,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'SIOOS ')
*
C     CALL MEMMAN(KSNOOSE,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'SNOOSE')
C     CALL MEMMAN(KSIOOSE,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'SIOOSE')
*. Combination form
C     CALL ZOOS(ISMOST(1,ISM),WORK(KLSBLTP),NSMST,WORK(KLSIOIO),
C    &          WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &          NOCTPA,NOCTPB,IDC,WORK(KSIOOS),WORK(KSNOOS),NSCMBC,0)
*. Sigma with expanded diagonal blocks
C     CALL ZOOS(ISMOST(1,ISM),WORK(KLSBLTP),NSMST,WORK(KLSIOIO),
C    &          WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &          NOCTPA,NOCTPB,IDC,WORK(KSIOOSE),WORK(KSNOOSE),NSCMBE,1)
*
* C information
*
      CALL MEMMAN(KLCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
      CALL MEMMAN(KLCBLTP,NSMST,'ADDL  ',2,'SBLTP ')
*
      CALL IAIBCM(INSPC,WORK(KLCIOIO))
      CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLCBLTP),WORK(KSVST))
*
C     CALL MEMMAN(KCNOOS,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'CNOOS ')
C     CALL MEMMAN(KCIOOS,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'CIOOS ')
*
C     CALL MEMMAN(KCNOOSE,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'CNOOSE')
C     CALL MEMMAN(KCIOOSE,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'CIOOSE')
*. Combination form
C     CALL ZOOS(ISMOST(1,ISM),WORK(KLCBLTP),NSMST,WORK(KLCIOIO),
C    &          WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &          NOCTPA,NOCTPB,IDC,WORK(KCIOOS),WORK(KCNOOS),NCCMBC,0)
*. Sigma with expanded diagonal blocks
C     CALL ZOOS(ISMOST(1,ISM),WORK(KLCBLTP),NSMST,WORK(KLCIOIO),
C    &          WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &          NOCTPA,NOCTPB,IDC,WORK(KCIOOSE),WORK(KCNOOSE),NSCMBE,1)
*. Additional block
      IF(LEVEL.LE.2) THEN
        CALL MXRESC(WORK(KLSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
        IF(IPRCIX.GE.2) THEN
          WRITE(6,*) 'DIRDIR : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
           WRITE(6,*) 'DIRDIR : MXADKBLK ', MXADKBLK
        END IF
        LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
        IF(IPRCIX.GE.2)
     &  WRITE(6,*) ' Space for resolution matrices ',LSCR2
        LSCR12 = MAX(LBLOCK,2*LSCR2)
        CALL MEMMAN(KVEC3,LSCR12,'ADDL  ',2,'KLB   ')
        KLB = KVEC3
      END IF
*
* Information on second order correction
*

      ISPC2 =  0
      IF(LEVEL.GE.3) THEN
*
        IF(LEVEL.EQ.3) THEN
         ISPC2 = IUTSPC
        ELSE
         ISPC2 = IUTSPC + 1
        END IF
*
        CALL MEMMAN(KLIOIO2,NOCTPA*NOCTPB,'ADDL  ',2,'IOIO2 ')
        CALL MEMMAN(KLBLTP2,NSMST,'ADDL  ',2,'BLTP2 ')
*
        CALL IAIBCM(ISPC2,WORK(KLIOIO2))
        KSVST = 1
        CALL ZBLTP(ISMOST(1,ISM),NSMST,IDC,WORK(KLBLTP2),WORK(KSVST))
*
C       CALL MEMMAN(KNOOS2,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'NOOS2 ')
C       CALL MEMMAN(KIOOS2,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'IOOS2 ')
*
C       CALL MEMMAN(KNOOS2E,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'NOOS2E')
C       CALL MEMMAN(KIOOS2E,NOCTPA*NOCTPB*NSMST,'ADDL  ',1,'IOOS2E')
*. Combination form
C       CALL ZOOS(ISMOST(1,ISM),WORK(KLBLTP2),NSMST,WORK(KLIOIO2),
C    &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &            NOCTPA,NOCTPB,IDC,WORK(KIOOS2),WORK(KNOOS2),NCMBC2,0)
*. Sigma with expanded diagonal blocks
C       CALL ZOOS(ISMOST(1,ISM),WORK(KLBLTP2),NSMST,WORK(KLIOIO2),
C    &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
C    &            NOCTPA,NOCTPB,IDC,WORK(KIOOS2E),WORK(KNOOS2E),
C    &            NCMB2E,1)
C?      write(6,*) ' Memcheck 4 '
C?      call memchk
C?      write(6,*) ' memchk passed '
*. Partitioning of second order correction
        NTTS = MXNTTS
        CALL MEMMAN(KLLBT2 ,NTTS  ,'ADDL  ',1,'LBT2  ')
        CALL MEMMAN(KLLEBT2 ,NTTS  ,'ADDL  ',1,'LEBT2 ')
        CALL MEMMAN(KLI1BT2,NTTS  ,'ADDL  ',1,'I1BT2 ')
        CALL MEMMAN(KLIBT2 ,8*NTTS,'ADDL  ',1,'IBT2  ')
*
C       CALL PART_CIV(IDC,WORK(KLBLTP2),WORK(KNOOS2E),WORK(KNOOS2),
C    &                NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLIOIO2),
C    &                ISMOST(1,ISM),
C    &                NBATCH2,WORK(KLLBT2),WORK(KLLEBT2),
C    &                WORK(KLI1BT2), WORK(KLIBT2),0)

C     PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &                  NSMST,MXLNG,IOCOC,ISMOST,
C    &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP)
        ITTSS_ORD = 2
        CALL PART_CIV2(IDC,WORK(KLBLTP2),WORK(KNSTSO(IATP)),
     &       WORK(KNSTSO(IBTP)),NOCTPA,NOCTPB,
     &       NSMST,LBLOCK,WORK(KLIOIO2),
     &       ISMOST(1,ISM),
     &       NBATCH2,WORK(KLLBT2),WORK(KLLEBT2),
     &       WORK(KLI1BT2), WORK(KLIBT2),0,ITTSS_ORD)
        WRITE(6,*) ' Number of batches for |2>   ', NBATCH2
        CALL MXRESC(WORK(KLIOIO2),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &              NSMST,NSTFSMSPGP,MXPNSMST,
     &              NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &              NELFSPGP,
     &              MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
        IF(IPRCIX.GE.2) THEN
          WRITE(6,*) 'DIRDIR : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                         MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
           WRITE(6,*) 'DIRDIR : MXADKBLK ', MXADKBLK
        END IF
        LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
        IF(IPRCIX.GE.2)
     &  WRITE(6,*) ' Space for resolution matrices ',LSCR2
        LSCR12 = MAX(LBLOCK,2*LSCR2)
        CALL MEMMAN(KVEC3,LSCR12,'ADDL  ',2,'KLB   ')
        KLB = KVEC3
       END IF
*


*. Obtain partitioning of (H-E)!0> vector
      NTTS = MXNTTS
      CALL MEMMAN(KLSLBT ,NTTS  ,'ADDL  ',1,'SLBT  ')
      CALL MEMMAN(KLSLEBT ,NTTS  ,'ADDL  ',1,'SLEBT ')
      CALL MEMMAN(KLSI1BT,NTTS  ,'ADDL  ',1,'SI1BT ')
      CALL MEMMAN(KLSIBT ,8*NTTS,'ADDL  ',1,'SIBT  ')
*
C     CALL PART_CIV(IDC,WORK(KLSBLTP),WORK(KSNOOSE),WORK(KSNOOS),
C    &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLSIOIO),
C    &              ISMOST(1,ISM),
C    &              NSBATCH,WORK(KLSLBT),WORK(KLSLEBT),
C    &              WORK(KLSI1BT), WORK(KLSIBT),0)
C     PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &                  NSMST,MXLNG,IOCOC,ISMOST,
C    &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP)
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KLSBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),
     &     NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLSIOIO),
     &     ISMOST(1,ISM),
     &     NSBATCH,WORK(KLSLBT),WORK(KLSLEBT),
     &     WORK(KLSI1BT), WORK(KLSIBT),0,ITTSS_ORD)
      WRITE(6,*) ' Number of batches for |1>   ', NSBATCH
*. Obtain partitionings of C vector
      CALL MEMMAN(KLCLBT ,NTTS  ,'ADDL  ',1,'CLBT  ')
      CALL MEMMAN(KLCLEBT ,NTTS  ,'ADDL  ',1,'CLEBT ')
      CALL MEMMAN(KLCI1BT,NTTS  ,'ADDL  ',1,'CI1BT ')
      CALL MEMMAN(KLCIBT ,8*NTTS,'ADDL  ',1,'CIBT  ')
*. for various partiotions
      CALL MEMMAN(KLLBT ,NTTS  ,'ADDL  ',1,'LBTC  ')
      CALL MEMMAN(KLLEBT,NTTS  ,'ADDL  ',1,'LECTC ')
      CALL MEMMAN(KLI1BT,NTTS  ,'ADDL  ',1,'I1BTC ')
      CALL MEMMAN(KLIBT ,8*NTTS,'ADDL  ',1,'IBTC  ')
      CALL MEMMAN(KLISCR,3*NTTS,'ADDL  ',1,'ISCR  ')
      CALL MEMMAN(KLISCR2,3*NTTS,'ADDL  ',1,'ISCR2 ')
*. The complete C vector (for disk reference )
C     CALL PART_CIV(IDC,WORK(KLCBLTP),WORK(KCNOOSE),WORK(KCNOOS),
C    &              NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
C    &              ISMOST(1,ISM),
C    &              NBATCHC,WORK(KLCLBT),WORK(KLCLEBT),
C    &              WORK(KLCI1BT),WORK(KLCIBT),1)
      ITTSS_ORD = 2
      CALL PART_CIV2(IDC,WORK(KLCBLTP),WORK(KNSTSO(IATP)),
     &     WORK(KNSTSO(IBTP)),
     &     NOCTPA,NOCTPB,NSMST,LBLOCK,WORK(KLCIOIO),
     &     ISMOST(1,ISM),
     &     NBATCHC,WORK(KLCLBT),WORK(KLCLEBT),
     &     WORK(KLCI1BT),WORK(KLCIBT),1,ITTSS_ORD)
C     PART_CIV2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
C    &                  NSMST,MXLNG,IOCOC,ISMOST,
C    &                  NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,ICOMP)
*. Total number of blocks in C
      NTBLOCKC = IFRMR(WORK(KLCI1BT),1,NBATCHC)
     &         + IFRMR(WORK(KLCLBT),1,NBATCHC) - 1
C     WRITE(6,*) ' DIRDIR1 : NTBLOCKC = ', NTBLOCKC
*. Total number of blocks in In Sigma ( First order correction )
      NTBLOCKS = IFRMR(WORK(KLSI1BT),1,NSBATCH)
     &         + IFRMR(WORK(KLSLBT),1,NSBATCH) - 1
C     WRITE(6,*) ' DIRDIR1 : NTBLOCKS = ', NTBLOCKS
*
* ======================================================
* Calculate first order correction to wave function and
* second order correction to energy
* ======================================================
*
       WRITE(6,*)
       WRITE(6,'(20X,A)') ' **********************************'
       WRITE(6,'(20X,A)') '  Second order energy calculation : '
       WRITE(6,'(20X,A)') ' **********************************'
       WRITE(6,*)
       WRITE(6,*)
*
       IF(IC1DSC.EQ.1) THEN
         CALL REWINE(LU1,-1)
       END IF
*
*. Start by studying whether |0> is an eigenfunction for H(apr).
*  This has some bearings upon how to evaluate (H0-E0)**-1
*. Calculate H(apr)|0> - <0!H(apr)!0> |0>
      CALL REWINE(LUC,-1)
      CALL REWINE(LUSCR,-1)
      ICSPC = ISPC0
      ISSPC = ISPC0
*. H(apr) |0>
      ZERO = 0.0D0
      CALL DIA0TRM_GAS(1,LUC,LUSCR,WORK(KCB),WORK(KSB),ZERO)
*. <0|H(apr)|0>
      ONE = 1.0D0
      HAPR00 = INPRDD(WORK(KCB),WORK(KSB),LUC,LUSCR,1,LBLK)
      CALL VECSMD(WORK(KCB),WORK(KSB),ONE,-HAPR00,
     &            LUSCR,LUC,LUSCR2,1,LBLK)
      XNORM = INPRDD(WORK(KCB),WORK(KSB),LUSCR2,LUSCR2,1,LBLK)
*
      IF(ABS(XNORM/HAPR00) .LE. 1.0D-12) THEN
       IHAPREIG = 1
      ELSE
       IHAPREIG = 0
      END IF
*. Choice of zero order energy
      IF(IE0AVEX.EQ.1) THEN
        E0 = HAPR00
      ELSE
        E0 = EREF
      END IF
*
      IF(IHAPREIG.EQ.0) THEN
*<0!(H(apr)-E0)**-1 |0>
        CALL DIA0TRM_GAS(2,LUC,LUSCR,WORK(KCB),WORK(KSB),-E0)
        HAPRM100  = INPRDD(WORK(KCB),WORK(KSB),LUC,LUSCR,1,LBLK)
        WRITE(6,*) ' HAPRM100', HAPRM100
      ENDIF
*. First order energy
      E1T = EREF - E0
*
      WRITE(6,*) ' Zero order energy  ', E0
      WRITE(6,*) ' First order energy ', E1T


      WRITE(6,*) '  HAPR00,  XNORM, IHAPREIG, HAPRM100 : ',
     &              HAPR00,  XNORM, IHAPREIG, HAPRM100
      CALL GFLUSH(6)
C     Call Abend2( ' Enforced stop in DIRDIR1N' )

*
*. Loop over partitionings of (H-E)!0> vector
      IFIRST = 1
      S01T = 0.0D0
      ICSPC = ISPC0
      ISSPC = ISPC1
      DO IBATCH = 1, NSBATCH
        IF(NTEST.GE.20) THEN
          WRITE(6,*) ' Information about sigma BATCH ', IBATCH
        END IF
*. Information about this batch
        NBLOCK = IFRMR(WORK(KLSLBT),1,IBATCH)
        NELMNT = IFRMR(WORK(KLSLEBT),1,IBATCH)
        IBLOCK_OFF = IFRMR(WORK(KLSI1BT),1,IBATCH)
*. Obtain blocks of first order corrections
        SIN = 0.0D0
        CALL GET_BATCH_OF_FIRST
     &       (NBLOCK,IBLOCK_OFF,NELMNT,WORK(KLSIBT),ISM,IUTSPC,
     &        E0,E1T,SIN,LUC,IDC,
     &        WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &        WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &        NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &        WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &        HAPR01,0,NTBLOCKS)
*. On output : LB : Batch of first order correction
*              SB : Batch of (H-E0)|0>
*
        IF(IC1DSC.EQ.1) THEN
*. Save batch of first order corrections on DISC '
C              EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
          CALL EXTRROW2(WORK(KLSIBT),8,IBLOCK_OFF,8,
     &                  NBLOCK,WORK(KLISCR))
          CALL TODSCN(WORK(KLB),NBLOCK,WORK(KLISCR),LBLK,LU1)
        END IF
*
        S01T = S01T + SOUT
*. Find contributions to
*. 1 : norm of nonorthogonalized first order correction vector
        CALL CLASS_PROD(WORK(KLB),WORK(KLB),NOCTPA,NOCTPB,
     &                  IBLOCK_OFF,NBLOCK,
     &                  WORK(KLSIBT),NOCCLS,IOCCLS,S11)
*. 2 : norm of residual vector
        CALL CLASS_PROD(WORK(KSB),WORK(KSB),NOCTPA,NOCTPB,
     &                  IBLOCK_OFF,NBLOCK,
     &                  WORK(KLSIBT),NOCCLS,IOCCLS,RES0)
*. 3 :  contributions to second order energy correction
        CALL CLASS_PROD(WORK(KSB),WORK(KLB),NOCTPA,NOCTPB,
     &                  IBLOCK_OFF,NBLOCK,
     &                  WORK(KLSIBT),NOCCLS,IOCCLS,E2)
      END DO
*. To put an end of to first order correction file
      IF(IC1DSC.EQ.1) THEN
        CALL ITODS(-1,1,LBLK,LU1)
      END IF
*. We wanted the norm of the residual so
      XNORMR = 0.0D0
      XNORM1 = 0.0D0
      DO JOCCLS = 1, NOCCLS
       XNORMR = XNORMR + RES0(JOCCLS)
       XNORM1 = XNORM1 + S11(JOCCLS)
       RES0(JOCCLS) = SQRT(RES0(JOCCLS))
      END DO
      XNORMR = SQRT(XNORMR)
      XNORM1 = SQRT(XNORM1)
*
      IF(IPART.EQ.1) IHYLLE2 = 1
*. (Not invoked for EN hamiltonian )
*
*. Calculate <1!H0-E0!1> to Obtain Hylleraas expression for second order
*  energy : <1!H0-E0!1> + 2 <0!H!1>
*
* |1> = |1'> - alpha1*(H apr -E )**-1 |0>
*
* <1!H0-E0!1> = <1!H apr -E0!1> = <1'!H apr -E0!1'>
*                               - <0|1'>**2 / <0!(H0-E0)-1!0>
      IF(IHYLLE2.EQ.1.AND.IPART.EQ.1) THEN
*
*
* Calculate <1'|H apr |'1> = sum(l,r) <0(l)|H apr |0(r)>
*                   = sum(l)    <0(l)|H apr |0(l)>
*                   + 2sum(r.lt.l) <0(l)|H apr |0(r)>
*
*. Loop over batches of <1'| and |1'>
        WRITE(6,*)
        WRITE(6,*) ' ==============================='
        WRITE(6,*) '  Hylleraas second order energy  '
        WRITE(6,*) ' ==============================='
        WRITE(6,*)
        CALL REWINE(LUSCR,-1)
        H011 = 0.0D0
*. Only one-electron part in sigma
        I12 = 1
*. Loop over left batches, i.e.  batches of H(apr)!1>
        DO ILBATCH = 1, NSBATCH
          WRITE(6,*) ' Starting outer loop for left  batch ', ILBATCH
          NLBLOCK = IFRMR(WORK(KLSLBT),1,ILBATCH)
          NLELMNT = IFRMR(WORK(KLSLEBT),1,ILBATCH)
          ILOFF = IFRMR(WORK(KLSI1BT),1,ILBATCH)
*. Loop over right batches, i.e.  batches of !1>
          DO IRBATCH  = 1, ILBATCH
            WRITE(6,*)
     &      '     Information about batches (left and right)',
     &      ILBATCH,IRBATCH
            NRBLOCK = IFRMR(WORK(KLSLBT),1,IRBATCH)
            NRELMNT = IFRMR(WORK(KLSLEBT),1,IRBATCH)
            IROFF = IFRMR(WORK(KLSI1BT),1,IRBATCH)
*
* Obtain first order correction |1'> in IRBATCH
*
            IF(IC1DSC.EQ.0) THEN
              SIN = 0.0D0
              CALL GET_BATCH_OF_FIRST
     &             (NRBLOCK,IROFF,NRELMNT,WORK(KLSIBT),ISM,IUTSPC,
     &              E0,E1T,SIN   ,LUC,IDC,
     &              WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &              WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &              NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &              WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &              HAPR01,0,NTBLOCKS)
            ELSE
*. Read in from DISC
              CALL GET_BLOCKS_FROM_DISC
     &        (LU1,NRBLOCK,IROFF,WORK(KLSIBT),NTBLOCKS,WORK(KLB),1)
            END IF
*
*. Transfer batch of first order correction to  LUSCR2
*  ( Sblock/Sblocks assumes that blocks are  on disc)
C              EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
          CALL EXTRROW2(WORK(KLSIBT),8,IROFF,8,
     &                  NRBLOCK,WORK(KLISCR))
          CALL REWINE(LUSCR2,-1)
          CALL TODSCNP(WORK(KLB),NRBLOCK,WORK(KLISCR),-1,LUSCR2)
          CALL ITODS(-1,1,-1,LUSCR2)
C         CALL ITODS(LREC(IREC),1,LBLK,LU)
*
*. Obtain (left) batch of H apr  |1(right_batch>
*
            I12 = 1
*. MP one-electron integrals
C           CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
            IRESTR = 1

            ICSPC = IUTSPC
            ISSPC = IUTSPC
            CALL SBLOCK(NLBLOCK,WORK(KLSIBT),ILOFF,WORK(KCB),WORK(KSB),
     &                  LUSCR2,IRESTR,0,1,IRBATCH,IRBATCH)
C                SBLOCK(NBLOCK,IBLOCK,IBOFF,CB,HCB,LUC,IRESTRICT,
C    &                  LUCBLK,ICBAT_RES,ICBAT_INI,ICBAT_END)
*. Add to previous obtained contributions to left batch of H apr !1>
            IF(IRBATCH.NE.1) THEN
              CALL REWINE(LUSCR,-1)
              CALL FRMDSC(WORK(KLB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
              ONE = 1.0D0
              CALL VECSUM(WORK(KSB),WORK(KSB),WORK(KLB),ONE,ONE,
     &                    NLELMNT)
            END IF
*. And transfer to disc to save for future generations
            CALL REWINE(LUSCR,-1)
            CALL TODSC(WORK(KSB),NLELMNT,-1,LUSCR)
C           CALL SBATCH_FROM_CBATCH(NLBLOCK,ILOFF,NRBLOCK,IROFF,
C    &                              WORK(KLSIBT),
C    &                              WORK(KLB),WORK(KSB),IRESTR)
*. Restore usual one-electron integrals
C           CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
          END DO
*         ^ End of loop over batches of !1'> (right batches)
*. Left Batch of H apr |1'> is now completed
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Batch of H0 |1> '
            CALL WRTTTS2(WORK(KSB),WORK(KLSIBT),NLBLOCK,ILOFF,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC)
          END IF
*
*. Obtain left batch of first order corrections
*
          ICSPC = INSPC
          ISSPC = IUTSPC
          IF(IC1DSC.EQ.0) THEN
            SIN = 0.0D0
            CALL GET_BATCH_OF_FIRST
     &           (NLBLOCK,ILOFF,NLELMNT,WORK(KLSIBT),ISM,IUTSPC,
     &            E0,E1T,SIN,LUC,IDC,
     &            WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &            WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &            HAPR01,0,NTBLOCKS)
          ELSE
*. Read in from DISC
              CALL GET_BLOCKS_FROM_DISC
     &        (LU1,NLBLOCK,ILOFF,WORK(KLSIBT),NTBLOCKS,WORK(KLB),1)
          END IF
*. Obtain saved sigma block
          CALL REWINE(LUSCR,-1)
          CALL FRMDSC(WORK(KSB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
          H011 = H011 + INPROD(WORK(KLB),WORK(KSB),NLELMNT)
        END DO
*       ^ End of loop over left batches
*
        IF(IRESTR.EQ.1) THEN
*. Only half of H011 was constructed,
          H011 = 2.0D0*H011
        END IF
      END IF
*     ^ End if for Hylleraas calculation
*
*
* Collect second order energy and print
*
      S11P = XNORM1 ** 2
      WRITE(6,*) ' H011 ', H011
      H0ME011P = H011 - E0*S11P
      WRITE(6,*) ' <1''!H0-E0!1''> ', H0ME011P
      IF(IHAPREIG.EQ.0) THEN
        H0ME011 = H0ME011P - S01T**2/HAPRM100
      ELSE
        H0ME011 = H0ME011P
      END IF
      WRITE(6,*)  ' <1!H0-E0!1> ',  H0ME011
*
      IF(IHAPREIG.EQ.0) THEN
        ALPHA1 = S01T/HAPRM100
      ELSE
        ALPHA1 = 0.0D0
      END IF
*
*. Orthogonalizetion term to second order energy
      E2ORT = ALPHA1*S01T
*. Total energy through second order
      E2T = 0.0D0
      DO JOCCLS = 1, NOCCLS
       E2T = E2T + E2(JOCCLS)
      END DO
      E2T = E2T + E2ORT
*. Hylleraas form of E2
      IF(IHYLLE2.EQ.1) E2TH = 2*E2T +  H0ME011
*. Total energy through second order
      E012T =  E0+ECORE+E1T+E2T
*. Total energy through second order , Hylleraas form
      IF(IHYLLE2.EQ.1) E012TH =  E0+ECORE+E1T+E2TH
*, Normalization factor
      FACN = 1.0D0/XNORM1
      WRITE(6,*) ' Occupation classes '
      CALL IWRTMA(IOCCLS,NGAS,NOCCLS,NGAS,NOCLLS)
*
C?    IF(NTEST.GE.1 ) THEN
        WRITE(6,*)
        WRITE(6,*) ' Contributions to second order energy :  '
        WRITE(6,*) ' ======================================= '
        WRITE(6,*)
        WRITE(6,'(6X,A)')
     &  'Class  Norm of !1(u)> Norm of (H-E)!0>   Second order energy'
        WRITE(6,'(6X,A)')
     &  ' ============================================================'
        DO JOCCLS = 1, NOCCLS
         WRITE(6,'(6x,I3,2X,E13.6,2X,E13.6,7X,E16.9)')
     &   JOCCLS,S11(JOCCLS),RES0(JOCCLS),E2(JOCCLS)
        END DO

        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Second order energy contribution ................... ',E2T
        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Second order energy Orthogonalization term ......... ',E2ORT
        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Energy through second order ........................ ',E012T
        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Norm of first order residual ....................... ',XNORMR
        WRITE(6,'(A,F18.9)')
     &  ' Norm of unorthogonalized first order vector ........ ',XNORM1
        WRITE(6,'(A,F18.9)')
     &  ' Overlap <0|(H0-E0)-1(H-E)|0> ....................... ',S01T
*
        IF(IHYLLE2.EQ.1) THEN
        WRITE(6,'(A,F18.9)')
     &  ' Hylleraas second order energy ...................... ',E2TH
        WRITE(6,'(A,F18.9)')
     &  ' Energy through second order ( Hylleraas )  ......... ',E012TH
        END IF
C     END IF
*
      IF(LEVEL.GE.2 .AND. IC1DSC.GT.0 .AND. ALPHA1 .NE. 0.0D0 ) THEN
*. Construct complete first order correction and store on disc
        WRITE(6,*)
        WRITE(6,*) ' ==============================================='
        WRITE(6,*) ' Construction of complete first order correction'
        WRITE(6,*) ' ==============================================='
        WRITE(6,*)
*
        HAPR01T = 0.0D0
        S11T = 0.0D0
        CALL REWINE(LU1,-1)
        DO IBATCH = 1, NSBATCH
          IF(NTEST.GE.20) THEN
            WRITE(6,*) ' Information about sigma BATCH ', IBATCH
          END IF
*. Information about this batch
          NBLOCK = IFRMR(WORK(KLSLBT),1,IBATCH)
          NELMNT = IFRMR(WORK(KLSLEBT),1,IBATCH)
          IBLOCK_OFF = IFRMR(WORK(KLSI1BT),1,IBATCH)
*. Obtain blocks of first order corrections
          CALL GET_BATCH_OF_FIRST
     &         (NBLOCK,IBLOCK_OFF,NELMNT,WORK(KLSIBT),ISM,IUTSPC,
     &          E0,E1T,ALPHA1,LUC,IDC,
     &          WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &          WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &          NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &          WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &          HAPR01,0,NTBLOCKS)
*. On output : LB : Batch of first order correction
          S11T = S11T + INPROD(WORK(KLB),WORK(KLB),NELMNT)
          HAPR01T = HAPR01T + HAPR01
*. Save batch of first order corrections on DISC '
          CALL EXTRROW2(WORK(KLSIBT),8,IBLOCK_OFF,8,
     &                  NBLOCK,WORK(KLISCR))
          CALL TODSCN(WORK(KLB),NBLOCK,WORK(KLISCR),LBLK,LU1)
        END DO
      END IF


      IF(LEVEL.GE.2) THEN
        WRITE(6,*)
        WRITE(6,'(20X,A)') ' *****************************'
        WRITE(6,'(20X,A)') '    One direct CI iteration : '
        WRITE(6,'(20X,A)') ' *****************************'
        WRITE(6,*)
        WRITE(6,*)
*
*
* Calculate <1|H|1> = sum(l,r) <0(l)|H|0(r)>
*                   = sum(l)    <0(l)|H|0(l)>
*                   + 2sum(r.lt.l) <0(l)|H|0(r)>
*
*. and <1|H0|1>
*. Loop over batches of <1| and |1>
        CALL REWINE(LUSCR,-1)
        H11 = 0.0D0
        H011 = 0.0D0
        IF(IC1DSC.EQ.0) THEN
          HAPR01T = 0.0D0
        END IF
        S11T = 0.0D0
        DO ILBATCH = 1, NSBATCH
          WRITE(6,*) ' Starting outer loop for left  batch ', ILBATCH
*. Construct this batch of first order corrections and store in LUSCR
          NLBLOCK = IFRMR(WORK(KLSLBT),1,ILBATCH)
          NLELMNT = IFRMR(WORK(KLSLEBT),1,ILBATCH)
          ILOFF = IFRMR(WORK(KLSI1BT),1,ILBATCH)
*
          IRESTRICT =  1
          IF(IRESTRICT.EQ.1) THEN
            IRMAX = ILBATCH
          ELSE
            IRMAX =  NSBATCH
          END IF
C         IRMAX = NSBATCH
          DO IRBATCH  = 1, IRMAX
*
            WRITE(6,*)
     &      '     Information about batches (left and right)',
     &      ILBATCH,IRBATCH
*
            NRBLOCK = IFRMR(WORK(KLSLBT),1,IRBATCH)
            NRELMNT = IFRMR(WORK(KLSLEBT),1,IRBATCH)
            IROFF = IFRMR(WORK(KLSI1BT),1,IRBATCH)
* calculate first order corrections in R batch
            ICSPC = INSPC
            ISSPC = IUTSPC
            CALL GET_BATCH_OF_FIRST
     &           (NRBLOCK,IROFF,NRELMNT,WORK(KLSIBT),ISM,IUTSPC,
     &            E0,E1T,ALPHA1,LUC,IDC,
     &            WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &            WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &            HAPR01,LU1,NTBLOCKS)
*. calculate those blocks of H |1(right_batch) that resides in
*  left batch
*
*. Transfer batch of first order correction to  LUSCR2
*  ( Sblock/Sblocks assumes that blocks are  on disc)
C                EXTRROW2(INMAT,IROW,ICOLOFF,NROW,NCOL,IOUTVEC)
            CALL EXTRROW2(WORK(KLSIBT),8,IROFF,8,
     &                    NRBLOCK,WORK(KLISCR))
            CALL REWINE(LUSCR2,-1)
            CALL TODSCNP(WORK(KLB),NRBLOCK,WORK(KLISCR),-1,LUSCR2)
            CALL ITODS(-1,1,-1,LUSCR2)
*
*. Obtain (left) batch of H apr  |1(right_batch>
*
            I12 = 2
            IPERTOP = 0
            IAPR = 0
            IRESTR = 1
*
            ICSPC = IUTSPC
            ISSPC = IUTSPC
            CALL SBLOCK(NLBLOCK,WORK(KLSIBT),ILOFF,WORK(KCB),WORK(KSB),
     &                  LUSCR2,IRESTR,0,1,IRBATCH,IRBATCH)
*. Add to previous obtained contributions to left batch of H !1>
            IF(IRBATCH.NE.1) THEN
              CALL REWINE(LUSCR,-1)
              CALL FRMDSC(WORK(KLB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
              ONE = 1.0D0
              CALL VECSUM(WORK(KSB),WORK(KSB),WORK(KLB),ONE,ONE,
     &                    NLELMNT)
            END IF
*. In last final batch, multiply with two to obtain full term
            IF(IRESTRICT.EQ.1.AND.IRBATCH.EQ.IRMAX) THEN
              TWO = 2.0D0
              CALL SCALVE(WORK(KSB),TWO,NLELMNT)
            END IF
*. And transfer to disc to save for future generations
            CALL REWINE(LUSCR,-1)
            CALL TODSCP(WORK(KSB),NLELMNT,-1,LUSCR)
          END DO
*         ^ End of loop over right batches
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Batch of H |1> '
            CALL WRTTTS2(WORK(KSB),WORK(KLSIBT),NLBLOCK,ILOFF,
     &                  NSMST,NOCTPA,NOCTPB,
     &                  WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC)
          END IF
*. construct left batch of first order corrections
          ICSPC = INSPC
          ISSPC = IUTSPC
          CALL GET_BATCH_OF_FIRST
     &         (NLBLOCK,ILOFF,NLELMNT,WORK(KLSIBT),ISM,IUTSPC,
     &          E0,E1T,ALPHA1,LUC,IDC,
     &          WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &          WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &          NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &          WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &          HAPR01,LU1,NTBLOCKS)
          IF(NTEST.GE.100) THEN
            WRITE(6,*) ' Left Batch of first order correction'
            CALL WRTTTS2(WORK(KLB),WORK(KLSIBT),NLBLOCK,ILOFF,
     &                   NSMST,NOCTPA,NOCTPB,
     &                   WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),IDC)
          END IF

          IF(IC1DSC.EQ.0)  HAPR01T = HAPR01T + HAPR01
          S11T = S11T + INPROD(WORK(KLB),WORK(KLB),NLELMNT)
* <1|H0!1>
          CALL COPVEC(WORK(KLB),WORK(KSB),NLELMNT)
          FACTOR = 0.0D0
          ITASK = 2
          ECOREX = 0.0D0
          JPERT = 1
          CALL DIATERM_GAS(FACTOR,ITASK,WORK(KSB),NLBLOCK,
     &         WORK(KLSIBT),ILOFF,JPERT,0,0)
          H011 = H011 + INPROD(WORK(KLB),WORK(KSB),NLELMNT)
*. <1|H|1>
          CALL REWINE(LUSCR,-1)
          CALL FRMDSC(WORK(KSB),NLELMNT,-1,LUSCR,IMZERO,IAMPACK)
*. Update H11
          H11 = H11 + INPROD(WORK(KLB),WORK(KSB),NLELMNT)
        END DO
*. Normalize to unit norm of first order correction
        XNORM = SQRT(S11T)
*
COLD    FACN = 1.0D0/XNORM1
        FACN = 1.0D0/XNORM
        H00 = E0+ E1T + ECORE
        H11 = H11/S11T + ECORE
        H01 = E2T/XNORM
C       H011 = H011 + ECORE*S11T
        WRITE(6,*) ' 2 by 2 matrix ( orthonormal vectors) '
        WRITE(6,*) '======================================'
        WRITE(6,*)
        WRITE(6,*) '          <0|H|0> '
        WRITE(6,*) '          <0|H|1>             <1|H|1> '
        WRITE(6,*)
        WRITE(6,'(4X,F22.12)') H00
        WRITE(6,'(4X,2F22.12)') H01,H11
*. In honour of B. Bechman :
        A =  1.0D0
        B = -(H00+H11)
        C = H00*H11-H01**2
*.
        EA = -B/(2*A) - SQRT(B**2 - 4*A*C)/(2*A)
        EB = -B/(2*A) + SQRT(B**2 - 4*A*C)/(2*A)
        ELOW = MIN(EA,EB)
        EHIGH = MAX(EA,EB)
        WRITE(6,'(A,2F22.12)') '   Eigenvalues ', ELOW,EHIGH
        WRITE(6,'(A)') '   Eigenvector for lowest eigenvalue '
        IF(H01.NE.0.0D0) THEN
          XLOW = (ELOW-H00)/H01
          XLNORM = SQRT(1.0D0+XLOW**2)
          FACTOR = 1.0D0/XLNORM
          WRITE(6,'(3X,2F22.12)')  FACTOR, XLOW*FACTOR
        ELSE
          IF(H00.LE.H11) THEN
            WRITE(6,'(A)') '1 ,  0 '
          ELSE
            WRITE(6,'(A)') '0 ,  1 '
          END IF
        END IF
*. Third order energy
        WRITE(6,*)
        WRITE(6,'(20X,A)')
     &  ' ****************************************'
        WRITE(6,'(20X,A)')
     &   '    Third order perturbation theory : '
        WRITE(6,'(20X,A)')
     &   ' ****************************************'
        WRITE(6,*)
        WRITE(6,*)
        H00 = E0
        H11 = H11 - ECORE
        WRITE(6,*) ' <1| H0 | 1> ', H011
        WRITE(6,*) '  <1| H0 -E0 | 1> ',H011-E0*S11T
        E3T = (H11-E1T)*S11T - H011
        WRITE(6,*) ' Third order energy correction ', E3T
        E0123T = E0 + ECORE + E1T +  E2T + E3T
        WRITE(6,*) ' Energy through third order ',
     &              E0+ECORE+E2T+E3T

        WRITE(6,*) '  <0| H apr | 1 > :' , HAPR01T
        WRITE(6,*) '   S11T = ', S11T
*
        XNORM1 = SQRT(S11T)
*
        WRITE(6,*) ' Contributions to third order energy :  '
        WRITE(6,*) ' ===================================== '
        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Third order energy contribution :................... ',E3T
        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Energy through third order :........................ ',E0123T
        WRITE(6,'(A,F18.9)')
     &  ' Norm of orthogonalized first order vector ..::...... ',XNORM1
      END IF
*
      IF(LEVEL.GE.3) THEN
        WRITE(6,*)
        WRITE(6,'(20X,A)')
     &  ' *****************************************'
        WRITE(6,'(20X,A)')
     &   '    Fourth order perturbation theory : '
        WRITE(6,'(20X,A)')
     &   ' *****************************************'
        WRITE(6,*)
        WRITE(6,*)
*. The fourth energy looks - in its most stable form as
*  E4 = <2|H0-E0|2> + 2<1|V-E1|2> - E2<1|1>
*. But the most easy - at hte moment is
*. E4 = <1|V-E1|2> - E2<1|1>
        IWAY = 3
        WRITE(6,*) ' expression used for fourth order energy '
        IF(IWAY.EQ.1) THEN
          WRITE(6,*) ' E4 = <2|H0-E0|2> + 2<1|V-E1|2> - E2<1|1> '
        ELSE IF(IWAY.EQ.2) THEN
          WRITE(6,*) ' E4 = <1|V-E1|2> - E2<1|1> '
        ELSE IF(IWAY.EQ.3) THEN
          WRITE(6,*) ' E4 =-<2|H0-E0|2> - E2<1|1> '
        END IF
*. Loop over batches of second order correction
        CALL REWINE(LUSCR,-1)
        H022 = 0.0D0
        V12 = 0.0D0
        S22U = 0.0D0
        S02T= 0.0D0
*
        DO IBATCH2 = 1, NBATCH2
          WRITE(6,*) ' Starting outer loop for batch ', IBATCH2
*. Construct this batch of first order corrections and store in LUSCR
          NBLOCK2 = IFRMR(WORK(KLLBT2) ,1,IBATCH2)
          NELMNT2 = IFRMR(WORK(KLLEBT2),1,IBATCH2)
          IOFF2 = IFRMR(WORK(KLI1BT2)  ,1,IBATCH2)
          write(6,*) ' NBLOCK2 NELMNT2 IOFF2 ',
     &                 NBLOCK2,NELMNT2,IOFF2
*
* calculate second order corrections in  batch
*
C         HAPR01TX = 0.0D0
          CALL GET_BATCH_OF_SECOND
     &         (NBLOCK2,IOFF2,NELMNT2,WORK(KLIBT2),ISM,ISPC2,
     &         NSBATCH,WORK(KLSLBT),WORK(KLSLEBT),WORK(KLSIBT),
     &         ISM,ISPC1,
     &         NBATCHC,WORK(KLCLBT),WORK(KLCIBT),ISPC0,
     &         E0,E1T,E2T,ALPHA1,LUC,LUSCR2,IDC,
     &         WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &         WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &         WORK(KLIBT),WORK(KLISCR),WORK(KLISCR2),
     &         WORK(KLB),WORK(KCB),WORK(KSB),SOUT,HAPR01T,LU1,
     &         NTBLOCKS)
          S22U = S22U + INPROD(WORK(KSB),WORK(KSB),NELMNT2)
          S02T = S02T + SOUT
*. Well now we have the second order correction in SB , save it !
          CALL REWINE(LUSCR,-1)
          CALL TODSC(WORK(KSB),NELMNT2,-1,LUSCR)
*. The term <2|H0-E0|2>
* (H0-E0)|2>
          CALL COPVEC(WORK(KSB),WORK(KLB),NELMNT2)
          FACTOR = - E0
          ITASK = 2
          ECOREX = 0.0D0
          IUTSPC = ISPC2
          JPERT = 1
          CALL DIATERM_GAS(FACTOR,ITASK,WORK(KLB),NBLOCK2,
     &         WORK(KLIBT2),IOFF2,JPERT,0,0)
C         DIATERM_GAS(FACTOR,ITASK,VEC,NBLOCK,IBLOCK,IOFF,
C    &                     ISM,ISPC,ECORE,I12)
          H022 = H022 + INPROD(WORK(KLB),WORK(KSB),NELMNT2)
*. The terms E2<1|2> and <1|V|2> : Loop over batches of first correction
          IF(IWAY.NE.3) THEN
          DO IBATCH1 = 1, NSBATCH
*
            WRITE(6,*)
     &      '     Information about batches (First - and second order)',
     &      IBATCH1,IBATCH2
*
            NBLOCK1 = IFRMR(WORK(KLSLBT),1,IBATCH1)
            NELMNT1 = IFRMR(WORK(KLSLEBT),1,IBATCH1)
            IOFF1 = IFRMR(WORK(KLSI1BT),1,IBATCH1)
* calculate first order corrections in batch
            ICSPC = ISPC0
            ISSPC = ISPC1
            CALL GET_BATCH_OF_FIRST
     &           (NBLOCK1,IOFF1,NELMNT1,WORK(KLSIBT),ISM,ISPC1,
     &            E0,E1T,ALPHA1,LUC,IDC,
     &            WORK(KLCIOIO),NOCTPA,NOCTPB,NSMST,
     &            WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
     &            NTBLOCKC,WORK(KLCIBT),WORK(KLIBT),WORK(KLI1BT),
     &            WORK(KLISCR),WORK(KSB),WORK(KCB),WORK(KLB),SOUT,
     &            HAPR01,LU1,NTBLOCKS)
*. First order correction is now in WORK(KLB)
*. save on LUSCR after second order batch
*. Position
            CALL REWINE(LUSCR,-1)
            CALL FRMDSC(WORK(KCB),NELMNT2,-1,LUSCR,IMZERO,IAMPACK)
*
            CALL TODSC(WORK(KLB),NELMNT1,-1,LUSCR)
*<2|H|1>
*
            ICSPC = ISPC1
            ISSPC = ISPC2
*
            ZERO = 0.0D0
            CALL SETVEC(WORK(KCB),ZERO,NELMNT2)
            I12 = 2
            IREST = 0
            CALL  SBATCH2_FROM_CBATCH(NBLOCK2,IOFF2,NBLOCK1,IOFF1,
     &                                 WORK(KLIBT2),WORK(KLSIBT),
     &                                 WORK(KLB),WORK(KCB),IREST)
*. Result is in in KCB
*. read second order batch into KSB
            CALL REWINE(LUSCR,-1)
            CALL FRMDSC(WORK(KSB),NELMNT2,-1,LUSCR,IMZERO,IAMPACK)
*
            V12 = V12 + INPROD(WORK(KSB),WORK(KCB),NELMNT2)
*. Obtain <1|H0-E1|2> and subtract
*. Second order blocks that are in first order  batch :
            CALL GET_TTS_BLK_IN_VECTOR(NBLOCK2,WORK(KLIBT2),IOFF2,
     &           WORK(KSB),NBLOCK1,WORK(KLSIBT),IOFF1,
     &           NBLOCK12,NELMNT12,WORK(KLIBT),WORK(KLB))
C           GET_TTS_BLK_IN_VECTOR
C    &             (NBLOCKI,IBLOCKI,IOFFI,VECI,
C    &              NBLOCKO,IBLOCKO,IOFFO,
C    &              NBLOCKIO,IBLCKIO,VECIO )
*. Extracted Second order blocks are now in KLB
*  H0+E1|2>
            FACTOR = E1T
            ITASK = 2
            ECOREX = 0.0D0
            JPERT = 1
            CALL DIATERM_GAS(FACTOR,ITASK,WORK(KLB),NBLOCK12,
     &           WORK(KLIBT),1,JPERT,0,0)
*. extract the corresponding first order blocks
*. Read in first order corrections in KCB
            CALL FRMDSC(WORK(KCB),NELMNT1,-1,LUSCR,IMZERO,IAMPACK)
*. And extract the common blocks
            CALL GET_TTS_BLK_IN_VECTOR(NBLOCK1,WORK(KLSIBT),IOFF1,
     &           WORK(KCB),NBLOCK12,WORK(KLIBT),1,
     &           NBLOCK12A,NELMNT12A,WORK(KLIBT),WORK(KSB))
*. Common   first order blocks are now in KSB
            V12 = V12 - INPROD(WORK(KLB),WORK(KSB),NELMNT12)
          END DO
          END IF
*.End of loop over batches of first order corrections
        END DO
*.End of loop over batches of second order corrections
*. The complete second order correction was not obtained,
*. It reads : |2> = |2'> - alpha2(H apr -E0 )-1 |0>
*. Only |2'> was obtained above
*. correct the V12 term
*. And <2|H0-E0|2>
        IF(IHAPREIG.NE.0) THEN
         V12 = V12
         ALPHA2 = 0.0D0
        ELSE
         ALPHA2 = S02T/HAPRM100
         V12 = V12 - ALPHA2*(-S02T+(E2T-HAPR01T)*HAPRM100)
         H022 = H022 -2*ALPHA2*S02T + ALPHA2**2 * HAPRM100
        END IF
*. and the norm
        S22 = S22U - S02T**2
        WRITE(6,'(A,F18.9)')
     &  ' Norm of orthogonalized second order correction vector ',
     &  SQRT( S22 )
*. The time has come to present the conclusions :
        WRITE(6,*)
        WRITE(6,*) ' Contributions to fourth order energy :  '
        WRITE(6,*) ' ======================================= '
        WRITE(6,*)
        IF(IWAY.EQ.1) THEN
          WRITE(6,'(A,F18.9)')
     &    '    <2|H0-E0|2> ..................... ', H022
          WRITE(6,'(A,F18.9)')
     &    '    2 * <1|V-E1|2> .................. ', 2.0D0*V12
          WRITE(6,'(A,F18.9)')
     &    '    -E2*<1|1> ....................... ', -E2T*S11T
          E4T = H022 + 2.0D0*V12 - E2T*S11T
        ELSE IF(IWAY.EQ.2) THEN
          WRITE(6,'(A,F18.9)')
     &    '    <1|V-E1|2> ...................... ', V12
          WRITE(6,'(A,F18.9)')
     &    '    -E2*<1|1> ....................... ', -E2T*S11T
          E4T = V12 - E2T*S11T
        ELSE IF (IWAY.EQ.3) THEN
          WRITE(6,'(A,F18.9)')
     &    '   -<2|H0-E0|2> ..................... ',-H022
          WRITE(6,'(A,F18.9)')
     &    '    -E2*<1|1> ....................... ', -E2T*S11T
          E4T = -H022 - E2T*S11T
        END IF
*
        WRITE(6,*)
        WRITE(6,'(A,F18.9)')
     &  ' Fourth order energy contribution ... ',E4T
        WRITE(6,*)
        E01234T =  E0+ECORE +E1T+E2T+E3T+E4T
        WRITE(6,'(A,F18.9)')
     &  ' Energy through fourth order ........ ',E01234T
        WRITE(6,*)

        WRITE(6,'(A,E18.9)')
     &  ' overlap between 2(unortho) and 0 ........ ',S02T

      END IF
*. End of fourth order part
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'DIRDR1')
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_BATCH_OF_FIRST
     &           (N1BLOCK,I1OFF,N1ELMNT,I1BLOCK,I1SM,I1SPC,
     &            E0,E1,ALPHA,LU0,IDC,
     &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &            N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,
     &            HAPR01,LU1,N1BLOCKT )
*
* Construct a part (a batch of blocks ) of the first order correction
*
*  |1> = -(H(apr)-E0) -1 (H-E0-E1)|0> - ALPHA * (H(apr)-E0)-1|0>
*
* Obtain  also  overlap <0|H(apr)|1> (Only  LU1 .EQ.0 )
*
* Jeppe Olsen, Written : October 1995
*              Modified : December 1995
*                         February 1996 : LU1 added
*
*
*
* =======
*  Input
* =======
*         N1BLOCK  : Number of blocks to be calculated
*         I1OFF    : Absolute number of first block to be calculated
*         N1ELMNT  : Number of ELEMENTS to be calculated
*         I1BLOCK  : List of all blocks in |1>
*         I1SM     : Symmetry of |1>
*         I1SPC    : Space of |1>
*         E0       : reference energy
*         SIN      : assumed overlap between |1> and |0>
*         LU0      : file containing 0
*         IDC      : packing in use
*         IOIO0    : allowed combination of alpha and beta supergroups for |0>
*         NOCTPA,NOCTPB : Number of alpha and beta types
*         NSMST   : Number of symmetries of strings
*         NSASO   : Number of alpha strings per sym and occupation type
*         NSBSO   : Number of beta  strings per sym and occupation type
*         N0BLOCK : Total number of blocks in |0>
*         I0BLOCK : Blocks in |0>
*         LU1     : file containing |1> ( if constructed previously)
*         N1BLOCKT: Total number of blocks in |1>
*
* ======
* Output
* ======
*
*
*     XB : Contains blocks of first order correction
*     SB : Contains blocks of (H-E)|0>                '
*     SOUT : Overlap <0| (H(apr)-E0) -1 (H-E0)|0>  (if LU1 = 0)
*     HAPR01 : Overlap <0|H(apr)|1>                (if LU1 = 0)
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IOIO0(NOCTPA,NOCTPB)
*. Block structure of |1>
      DIMENSION I1BLOCK(8,*)
*. Block structure of |0>
      DIMENSION I0BLOCK(8,*)
*
      DIMENSION SB(*),CB(*),XB(*)
*
#include "oper.inc"
C     COMMON/OPER/I12,IPERTOP,IAPR,MNRS1E,MXRS3E,IPART

      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_FI')
      NTEST = 3000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information from  GET_FIRST_ORDER_CORRECTION '
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks to be constructed ', N1BLOCK
        WRITE(6,*) ' TTSS for Blocks to be constructed '
        DO IBLOCK = I1OFF, I1OFF-1+N1BLOCK
          WRITE(6,'(10X,4I3,2I8)') (I1BLOCK(II,IBLOCK),II=1,4)
        END DO
        IF(NTEST.GE.2) WRITE(6,*) ' First order space : ', I1SPC
      END IF
      CALL GFLUSH(6)
*
      IF(LU1.GT.0) THEN
*
*. Easy Living, first order correction is stored on LU1
* ======================================================

C             GET_BLOCKS_FROM_DISC(LU,NBLOCK,IOFF,IBLOCK,NBLOCKT,C)
         CALL GET_BLOCKS_FROM_DISC
     &   (LU1,N1BLOCK,I1OFF,I1BLOCK,N1BLOCKT,XB,1)
*. Blocks should be scaled if IDC.NE.1 ???????
         HAPR01 = 0.0D0
*. We are also interested in
      ELSE
*
*. Well, not so easy living, I/WE ( maybe I am a parallel code by now)
*. have to construct the first order correction from SCRATCH
* ===================================================================
*.
*
*
* ============================
*. Obtain sigma blocks: H |0>
* ============================
*
      I12 = 2
      IPERTOP = 0
      IAPR = 0
C     write(6,*) ' call to SBLOCK : KVEC3 defined ?'
C     Call Abend2( 'before call to SBLOCK' )
      CALL SBLOCK(N1BLOCK,I1BLOCK,I1OFF,CB,SB,LU0,0,0,0,0,0)
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' Batch of blocks after SBLOCK'
        CALL WRTTTS(SB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
      IF(NTEST.GE.2) write(6,*) ' Returning from Sblock '
*. Sigma Blocks that are in actual C space
      CALL REDBLK(N1BLOCK,I1BLOCK,I1OFF,IOIO0,
     &            NOCTPA,NOCTPB,N10BLOCK,IXBLOCK,IX1)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from REDBLK '
*.Fetch the corresponding C blocks.
      IREW = 1
      ISCAL = 0
      CALL GET_TTS_BATCH(CB,N10BLOCK,IXBLOCK,N0BLOCK,I0BLOCK,
     &                   NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                   IDC,LU0,IX2,IREW,ISCAL)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from GET_TTS'
* ====================================
*. Subtract E0+E1 |0> to give (H-(E0+E1)) |0>
* ====================================
      FACTOR = - (E0 + E1)
      CALL ADDBLKV(SB,CB,FACTOR,
     &            N10BLOCK,IX1,I1BLOCK,I1OFF)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from ADDBLKV'
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' Batch of blocks after ADDBLKV'
        CALL WRTTTS(SB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
*
* ====================================
* obtain block of (H0-E0)-1**(H-E) !0>
* ====================================
*
      CALL COPVEC(SB,XB,N1ELMNT)
      FACTOR = - E0
      ECOREX = 0.0D0
      JPERT = 1
      CALL DIATERM_GAS(FACTOR,1,XB,N1BLOCK,
     &            I1BLOCK,I1OFF,JPERT,0,0)
      IF(NTEST.GE.2000) THEN
        WRITE(6,*) ' Batch of blocks after DIATERM_GAS'
        CALL WRTTTS(XB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from DIATERM'
      ONEM = -1.0D0
      CALL SCALVE(XB,ONEM,N1ELMNT)
* =====================================
*. overlap <0|(H(apr)-E0) -1 (E0-H)|0>
* =====================================
      I2EQ1 = 0
      SOUT = SCALAR_PRODUCT_OF_B
     &       (XB,CB,N1BLOCK,I1BLOCK,I2EQ1,N10BLOCK,IX1,I1OFF)
*
* =====================================
* Add -ALPHA * (H0-E0)**-1 |0>
* =====================================
*
      IF(ALPHA.NE.0.0D0) THEN
*. (H(apr)-E0)**-1|0>
        FACTOR = -E0
        ECOREX = 0
        JPERT = 1
        CALL DIATERM_GAS(FACTOR,1,CB,N10BLOCK,IXBLOCK,
     &                    1,JPERT,0,0)
        FACTOR = -ALPHA
        CALL ADDBLKV(XB,CB,FACTOR,N10BLOCK,IX1,I1BLOCK,I1OFF)
      END IF
*
*
* ================================
*. Contribution to <0| H(apr) |1>
* ================================
*
*. Fetch the C blocks again
      IREW = 1
      ISCAL = 0
      CALL GET_TTS_BATCH(CB,N10BLOCK,IXBLOCK,N0BLOCK,I0BLOCK,
     &                   NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                   IDC,LU0,IX2,IREW,ISCAL)
*. H(apr)|0>
      FACTOR = 0
      ECOREX = 0
      JPERT =  1
      CALL DIATERM_GAS(FACTOR,2,CB,N10BLOCK,IXBLOCK,
     &                  1,JPERT,0,0)
      HAPR01 = SCALAR_PRODUCT_OF_B
     &       (XB,CB,N1BLOCK,I1BLOCK,I2EQ1,N10BLOCK,IX1,I1OFF)
*.
*
      END IF
*. ( End of LU1=0/not 0 Branching)
*
*.
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' Batch of first order correction blocks '
        CALL WRTTTS(XB,I1BLOCK(1,I1OFF),N1BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_FI')
*
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE GET_BATCH_OF_SECOND
     &           (N2BLOCK,I2OFF,N2ELMNT,I2BLOCK,I2SM,I2SPC,
     &            N1BATCH,N1BLOCK,L1BLOCK,I1BLOCK,I1SM,I1SPC,
     &            N0BATCH,N0BLOCK,I0BLOCK,I0SPC,
     &            E0,E1,E2,SIN,LU0,LUSCR,IDC,
     &            IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &            IXBLOCK,IX1,IX2,SB,CB,XB,SOUT,HAPR01T,LU1,
     &            N1BLOCKT)
*
* Construct a part (a batch of blocks ) of the second order correction
*
*  |2> = (H0-E0) -1 {(E1-V)|1> + E2 * |0>}  - S02 |0>
*
* If LU1 .gt. 0, then it is assumed that |1> has been constructed and is
* stored on LU1
*
* Jeppe Olsen, October 1995
*              Feb 96 : LU1 added
*
*. Note : The terms S02 |0> is  not included here, since S02 is in genral not
*. Known
*. Usually our H0 has the form
*
* H0 = Q H apr Q + E0 P
* so the second order correction can be written
*
* |2> = Q (H apr -E0 ) -1 {(E1-H + Hapr)|1> + ( E2 - <0|H(apr)|1>) * |0>}
*
*
*
* =======
*  Input
* =======
*         N2BLOCK  : Number of blocks to be calculated
*         I2OFF    : Absolute number of first block to be calculated
*         N2ELMNT  : Number of ELEMENTS to be calculated
*         I2BLOCK  : List of all blocks in |2>
*         I2SM     : Symmetry of |2>
*         I2SPC    : Space of |2>
*         N1BATCH  : Number of batches in |1>
*         N1BLOCK  : Number of blocks in each batch of |1>
*         L1BLOCK  : Number of elements in each batch of |1>
*         I1BLOCK  : Blocks in |1>
*         E0       : reference energy
*         SIN      : assumed overlap between |1> and |0>
*         LU0      : file containing 0
*         IDC      : packing in use
*         IOIO0    : allowed combination of alpha and beta supergroups for |0>
*         NOCTPA,NOCTPB : Number of alpha and beta types
*         NSMST    : Number of symmetries of strings
*         NSASO    : Number of alpha strings per sym and occupation type
*         NSBSO    : Number of beta  strings per sym and occupation type
*         N0BLOCK  : Total number of blocks in |0>
*         I0BLOCK  : Blocks in |0>
*         E1,E2    : First and second order energy corrections
*         LUSCR    : Scratch file
*         LU1      : Nonvanishing indicates that |1> is stored on LU1
*         N1BLOCKT : Total number of blocks in |1>
*
* ======
* Output
* ======
*
*
*     XB : Contains blocks of first order correction
*     SOUT : Overlap between |2> as above  and |0>
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IOIO0(NOCTPA,NOCTPB)
*. BLock structure of |2>
      DIMENSION I2BLOCK(8,*)
*. Block structure of |1>
      DIMENSION I1BLOCK(8,*)
      DIMENSION N1BLOCK(*), L1BLOCK(*)
*. Block structure of |0>
      DIMENSION I0BLOCK(8,*)
*. Scratch
      DIMENSION IXBLOCK(8,*)
*
      DIMENSION SB(*),CB(*),XB(*)
* For communcating with MV7
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
#include "oper.inc"
C     COMMON/OPER/I12,IPERTOP,IAPR,MNRS1E,MXRS3E,IPART

*
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'GET_SE')
*
      WRITE(6,*) ' KVEC3 properly defined 3 block is used '
      Call Abend2( 'new implementation of kvec3' )
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' Information from  GET_SECOND_ORDER_CORRECTION '
        WRITE(6,*)
        WRITE(6,*) ' Number of blocks to be constructed ', N2BLOCK
        WRITE(6,*) ' TTSS for Blocks to be constructed '
        DO IBLOCK = I2OFF, I2OFF-1+N2BLOCK
          WRITE(6,'(10X,4I3,2I8)') (I2BLOCK(II,IBLOCK),II=1,4)
        END DO
        CALL GFLUSH(6)
      END IF
      IF(NTEST.GE.1)
     & WRITE(6,*) ' Entering GET_BATCH_OF_SECOND'
*
* ====================
* 1 : (V - E1 ) | 1 >
* ====================
*
*. Loop over partionings of |1>
      DO NN1BATCH = 1, N1BATCH
        IF(NN1BATCH.EQ.1) THEN
          I1OFF = 1
        ELSE
          I1OFF = I1OFF + NN1BLOCK
        END IF
        NN1BLOCK = N1BLOCK(NN1BATCH)
        LL1BLOCK = L1BLOCK(NN1BATCH)
*. Obtain this batch of first order corrections
C       INSPC = I0SPC
C       IUTSPC = I1SPC
*
        ICSPC = I0SPC
        ISSPC = I1SPC
*
        IF(NTEST.GE.200) write(6,*) ' NN1BATCH I1OFF NN1BLOCK',
     &               NN1BATCH,I1OFF,NN1BLOCK
        CALL GFLUSH(6)
        CALL GET_BATCH_OF_FIRST
     &       (NN1BLOCK,I1OFF,LL1BLOCK,I1BLOCK,I1SM,I1SPC,
     &        E0,E1,SIN,LU0,IDC,
     &        IOIO0,NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &        N0BLOCK,I0BLOCK,IXBLOCK,IX1,IX2,
     &        SB,XB,CB,SOUT2,
     &        HAPR01,LU1,N1BLOCKT)
*. First order correction returned in CB
        IF(NTEST.GE.5) WRITE(6,*)
     &   ' Returning from GET_BATCH_OF_FIRST'
* (E1-V)!|1> = (E1-H+H0+E1)|1>
*
*.1.1 : H |1(batch)>
*
*
         ICSPC = I1SPC
         ISSPC = I2SPC
*
         IF(NN1BATCH.EQ.1)  THEN
           ZERO = 0.0D0
           CALL SETVEC(SB,ZERO,N2ELMNT)
         ELSE
*. retrieve previous contributions to (V-E1) |1>
           CALL REWINE(LUSCR,-1)
           CALL FRMDSC(SB,N2ELMNT,-1,LUSCR,IMZERO,IAMPACK)
         END IF
         CALL COPVEC(CB,XB,LL1BLOCK)
C?       WRITE(6,*) ' GET_SEC : Just before SBATCH2 '
         CALL GFLUSH(6)
         I12 = 2
         IRESTR = 0
         CALL SBATCH2_FROM_CBATCH(N2BLOCK,I2OFF,NN1BLOCK,I1OFF,
     &                            I2BLOCK,I1BLOCK,XB,SB,IRESTR)
C?       WRITE(6,*) ' GET_SEC : Home form SBATCH2 '
         CALL GFLUSH(6)
*. SB contains now previous contributions to (V-E1) |1>
*. and H |1(batch)>, save it !
         CALL REWINE(LUSCR,-1)
         CALL TODSC(SB,N2ELMNT,-1,LUSCR)
*
*
*. 1.2 : H0+E1 |1(batch)>
*
* obtain blocks in current batch of |1> this is in the space of current
* batch of |2>.
C                GET_TTS_BLK_IN_VECTOR
C    &           (NBLOCKI,IBLOCKI,IOFFI,VECI,
C    &            NBLOCKO,IBLOCKO,IOFFO,
C    &            NBLOCKIO,IBLCOKIO,VECIO )
        CALL GET_TTS_BLK_IN_VECTOR
     &          (NN1BLOCK,I1BLOCK,I1OFF,CB,
     &          N2BLOCK,I2BLOCK,I2OFF,
     &          N12BLOCK,L12BLOCK,IXBLOCK,XB)
         IF(NTEST.GE.5)
     &   WRITE(6,*) ' GET___SEC : Home form GET_TTS '
*
        FACTOR =  E1
        ITASK = 2
        ECOREX = 0.0D0
        JPERT = 1
        CALL DIATERM_GAS(FACTOR,ITASK,XB,N12BLOCK,
     &       IXBLOCK,1,JPERT,0,0)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' H0+E1 |1 > '
          CALL WRTTTS(XB,I2BLOCK(1,I2OFF),N2BLOCK,
     &                NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
        END IF
*. XB contains now relevant block of (H0+E1), subtract from prev + H|1(batch)>
        ONEM = -1.0D0
        CALL ADD_TTS_BLK_TO_VECTOR
     &       (N12BLOCK,IXBLOCK,1,XB,
     &        N2BLOCK,I2BLOCK,I2OFF,SB,ONEM)
        IF(NTEST.GE.1000) THEN
          WRITE(6,*) ' (V - E1 ) | 1 > '
          CALL WRTTTS(SB,I2BLOCK(1,I2OFF),N2BLOCK,
     &                NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
        END IF
*. SB contains now sum(i1batch = 1, n1batch) (V-E1)|1(ibatch)>, save it
         CALL REWINE(LUSCR,-1)
         CALL TODSC(SB,N2ELMNT,-1,LUSCR)
      END DO
*. It was actually (E1-V)|1> we were interested in so
      CALL SCALVE(SB,ONEM,N2ELMNT)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' Final (E1 - V ) | 1 > '
        CALL WRTTTS(SB,I2BLOCK(1,I2OFF),N2BLOCK,
     &                NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
*
*==================================
* 2. Add (E2 - <0|H(apr)|1>) * |0>
*==================================
*
*. Blocks of |0> in current batch of |2>
      CALL REDBLK(N2BLOCK,I2BLOCK,I2OFF,IOIO0,
     &            NOCTPA,NOCTPB,N20BLOCK,IXBLOCK,IX1)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from REDBLK '
*.Fetch the corresponding C blocks.
      IREW = 1
      ISCAL=0
      CALL GET_TTS_BATCH(CB,N20BLOCK,IXBLOCK,N0BLOCK,I0BLOCK,
     &                   NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                   IDC,LU0,IX2,IREW,ISCAL)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from GET_TTS'
      FACTOR = E2 - HAPR01T
C?    WRITE(6,*) ' GET_SEC : FACTOR ', FACTOR
      CALL ADDBLKV(SB,CB,FACTOR,
     &            N20BLOCK,IX1,I2BLOCK,I2OFF)
      IF(NTEST.GE.1000) THEN
        WRITE(6,*) ' (E1 - V) | 1 > + E2 | 0 >'
        CALL WRTTTS(SB,I2BLOCK(1,I2OFF),N2BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)
      END IF
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from ADDBLKV'
*
* ==============================================
* 3 : Multiply with (H0-E0) to give
*     |2> = (H0-E0)-1 ((E1-V)|1> + E2 * |0>)
* ==============================================
      CALL COPVEC(SB,XB,N2ELMNT)
      FACTOR = - E0
      ECOREX = 0.0D0
      JPERT = 1
      CALL DIATERM_GAS(FACTOR,1,XB,N2BLOCK,
     &            I2BLOCK,I2OFF,JPERT,0,0)
      IF(NTEST.GE.10)
     &write(6,*) ' Returning from DIATERM'
*. overlap between |2> and |0>
      I2EQ1 = 0
      SOUT = SCALAR_PRODUCT_OF_B
     &       (XB,CB,N2BLOCK,I2BLOCK,I2EQ1,N20BLOCK,IX1,I2OFF)
*
      IF(NTEST.GE.200) THEN
        WRITE(6,*) ' Batch of second order correction blocks '
        CALL WRTTTS(XB,I2BLOCK(1,I2OFF),N2BLOCK,
     &              NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,IDC)

      END IF
*
      CALL MEMMAN(IDUM,IDUM,'FLUSM ',IDUM,'GET_SE')
      IF(NTEST.GE.1)
     & WRITE(6,*) ' Leaving GET_BATCH_OF_SECOND'
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE SBATCH2_FROM_CBATCH(NLBLOCK,ILOFF,NRBLOCK,IROFF,
     &                               ILBLOCK,IRBLOCK,
     &                               CB,HCB,IRESTRICT)
*
* Generate the parts of H times right batch that is in a given
* left batch
*
* The complete set of blocks of the vector is specified by IBLOCK
*
* The right batch is specified by IROFF,NRBLOCK
* The left  batch is specified by ILOFF,NLBLOCK
*
* the input right batch is assumed to reside in CB
*
* The output blocks are delivered in HCB
*
* The blocks are scaled and reformed to combination order
*
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "mxpdim.inc"
*
* =====
*.Input
* =====
*
*.Definition of c and sigma spaces
      COMMON/CANDS/ICSM,ISSM,ICSPC,ISSPC
*. Sigma blocks require
      INTEGER ILBLOCK(8,*)
*. C blocks included
      INTEGER IRBLOCK(8,*)
*
*./ORBINP/ : NACOB used
#include "orbinp.inc"
#include "cicisp.inc"
#include "strbas.inc"
#include "cstate.inc"
#include "strinp.inc"
#include "stinf.inc"
#include "csm.inc"
#include "wrkspc.inc"
#include "crun.inc"
#include "gasstr.inc"
#include "cgas.inc"
*. Used : NSMOB
#include "lucinp.inc"
#include "cprnt.inc"
#include "glbbas.inc"
#include "oper.inc"
*
      INTEGER ADASX,ASXAD,ADSXA,SXSXDX,SXDXSX
      COMMON/CSMPRD/ADASX(MXPOBS,MXPOBS),ASXAD(MXPOBS,2*MXPOBS),
     &              ADSXA(MXPOBS,2*MXPOBS),
     &              SXSXDX(2*MXPOBS,2*MXPOBS),SXDXSX(2*MXPOBS,4*MXPOBS)
      COMMON/HIDSCR/KLOCSTR(4),KLREO(4),KLZ(4),KLZSCR
*
      NTEST = 0
      IDUM = 0
      CALL MEMMAN(IDUM,IDUM,'MARK  ',IDUM,'SBFCB ')
C?    WRITE(6,*) ' LUC in SBLOCK ', LUC
*
* Info for this internal space
*. type of alpha and beta strings
      IATP = 1
      IBTP = 2
*. alpha and beta strings with an electron removed
      IATPM1 = 3
      IBTPM1 = 4
*. alpha and beta strings with two electrons removed
      IATPM2 = 5
      IBTPM2 = 6
*
*. Number of supergroups
      NOCTPA = NOCTYP(IATP)
      NOCTPB = NOCTYP(IBTP)
*. Offset for supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*
      NAEL = NELEC(IATP)
      NBEL = NELEC(IBTP)
* string sym, string sym => sx sym
* string sym, string sym => dx sym
      CALL MEMMAN(KSTSTS,NSMST ** 2,'ADDL  ',2,'KSTSTS')
      CALL MEMMAN(KSTSTD,NSMST ** 2,'ADDL  ',2,'KSTSTD')
      CALL STSTSM(WORK(KSTSTS),WORK(KSTSTD),NSMST)
*. Largest block of strings in zero order space
      MXSTBL0 = MXNSTR
*. Largest number of strings of given symmetry and type
      MAXA = 0
      IF(NAEL.GE.1) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM1)),NSMST*NOCTYP(IATPM1),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      IF(NAEL.GE.2) THEN
        MAXA1 = IMNMX(WORK(KNSTSO(IATPM2)),NSMST*NOCTYP(IATPM2),2)
        MAXA = MAX(MAXA,MAXA1)
      END IF
      MAXB = 0
      IF(NBEL.GE.1) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM1)),NSMST*NOCTYP(IBTPM1),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      IF(NBEL.GE.2) THEN
        MAXB1 = IMNMX(WORK(KNSTSO(IBTPM2)),NSMST*NOCTYP(IBTPM2),2)
        MAXB = MAX(MAXB,MAXB1)
      END IF
      MXSTBL = MAX(MAXA,MAXB,MXSTBL0)
      IF(IPRCIX.GE.2 ) WRITE(6,*)
     &' Largest block of strings with given symmetry and type',MXSTBL
*. Largest number of resolution strings and spectator strings
*  that can be treated simultaneously
      MAXI = MIN( MXINKA,MXSTBL)
      MAXK = MIN( MXINKA,MXSTBL)
*Largest active orbital block belonging to given type and symmetry
      MXTSOB_L = 0
      DO IOBTP = 1, NGAS
        DO IOBSM = 1, NSMOB
          MXTSOB_L = MAX(MXTSOB_L,NOBPTS(IOBTP,IOBSM))
        END DO
      END DO
C?    WRITE(6,*) ' MXTSOB_L = ', MXTSOB_L
      MAXIJ = MXTSOB_L ** 2
*.Local scratch arrays for blocks of C and sigma
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' ICISTR,LBLOCK ',ICISTR,LBLOCK

*.vectors able to hold strings of given sym and type
C     MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
C     LSCR3 = MAX(MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
C     CALL MEMMAN(KI1,LSCR3       ,'ADDL  ',1,'I1    ')
C     CALL MEMMAN(KXI1S,LSCR3       ,'ADDL  ',2,'XI1S  ')
*
C     CALL MEMMAN(KI2,LSCR3       ,'ADDL  ',1,'I2    ')
C     CALL MEMMAN(KXI2S,LSCR3       ,'ADDL  ',2,'XI2S  ')
*
C     CALL MEMMAN(KI3,MAXIK*MXTSOB_L,'ADDL  ',1,'I3    ')
C     CALL MEMMAN(KXI3S,MAXIK*MXTSOB_L,'ADDL  ',2,'XI3S  ')
*
C     CALL MEMMAN(KI4,MAXIK*MXTSOB_L,'ADDL  ',1,'I4    ')
C     CALL MEMMAN(KXI4S,MAXIK*MXTSOB_L,'ADDL  ',2,'XI4S  ')
*.SCRATCH space for integrals
* A 4 index integral block with four indeces belonging OS class
      INTSCR = MXTSOB_L ** 4
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Integral scratch space ',INTSCR
      CALL MEMMAN(KINSCR,INTSCR,'ADDL  ',2,'INSCR ')
*. Arrays giving allowed type combinations '
      CALL MEMMAN(KCIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'CIOIO ')
      CALL MEMMAN(KSIOIO,NOCTPA*NOCTPB,'ADDL  ',2,'SIOIO ')
*. Offsets for alpha and beta supergroups
      IOCTPA = IBSPGPFTP(IATP)
      IOCTPB = IBSPGPFTP(IBTP)
*. sigma needed for MXRESC
C          IAIBCM(ICISPC,IAIB)
      CALL IAIBCM(ISSPC,WORK(KSIOIO))
      CALL IAIBCM(ICSPC,WORK(KCIOIO))
*. Arrays giving block type
      CALL MEMMAN(KCBLTP,NSMST,'ADDL  ',2,'CBLTP ')
*. Arrays for additional symmetry operation
      IF(IDC.EQ.3.OR.IDC.EQ.4) THEN
        CALL MEMMAN(KSVST,NSMST,'ADDL  ',2,'SVST  ')
        CALL SIGVST(WORK(KSVST),NSMST)
      ELSE
         KSVST = 1
      END IF
*
*.scratch space for projected matrices and a CI block
*
*. Scratch space for CJKAIB resolution matrices
*. Size of C(Ka,Jb,j),C(Ka,KB,ij)  resolution matrices
      CALL MXRESC(WORK(KSIOIO),IOCTPA,IOCTPB,NOCTPA,NOCTPB,
     &            NSMST,NSTFSMSPGP,MXPNSMST,
     &            NSMOB,MXPNGAS,NGAS,NOBPTS,IPRCIX,MAXK,
     &            NELFSPGP,
     &            MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL,MXADKBLK)
      IF(IPRCIX.GE.2) THEN
        WRITE(6,*) 'SBLOCK : MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL',
     &                       MXCJ,MXCIJA,MXCIJB,MXCIJAB,MXSXBL
        WRITE(6,*) 'SBLOCK : MXADKBLK ', MXADKBLK
      END IF
      LSCR2 = MAX(MXCJ,MXCIJA,MXCIJB,MXCIJAB)
      IF(IPRCIX.GE.2)
     &WRITE(6,*) ' Space for resolution matrices ',LSCR2
*
      IF(IPRCIX.GE.2)  WRITE(6,*) ' LSCR2 = ', LSCR2
*
      LSCR12 = MAX(LBLOCK,2*LSCR2)
      CALL MEMMAN(KC2,LSCR12,'ADDL  ',2,'KC2   ')
      KCJRES = KC2
      KSIRES = KC2 + LSCR2
*
      KSSCR = KSIRES
      KCSCR = KCJRES
*
*.vectors able to hold strings of given sym and type
      MAXIK = MAX(MAXI,MAXK)
*. I1 and Xi1s must also be able to hold largest st block
      LSCR3 = MAX(MXADKBLK,MAXIK*MXTSOB_L*MXTSOB_L,MXSTBL0)
      CALL MEMMAN(KI1  ,LSCR3,'ADDL  ',1,'I1    ')
      CALL MEMMAN(KXI1S,LSCR3,'ADDL  ',2,'XI1S  ')
*
      CALL MEMMAN(KI2  ,LSCR3,'ADDL  ',1,'I2    ')
      CALL MEMMAN(KXI2S,LSCR3,'ADDL  ',2,'XI2S  ')
*
      CALL MEMMAN(KI3  ,LSCR3,'ADDL  ',1,'I3    ')
      CALL MEMMAN(KXI3S,LSCR3,'ADDL  ',2,'XI3S  ')
*
      CALL MEMMAN(KI4  ,LSCR3,'ADDL  ',1,'I4    ')
      CALL MEMMAN(KXI4S,LSCR3,'ADDL  ',2,'XI4S  ')
*
      CALL ZBLTP(ISMOST(1,ICSM),NSMST,IDC,WORK(KCBLTP),WORK(KSVST))
*.Some TTS arrays
      NOOS = NSMCI*NOCTPA*NOCTPB
      NTTS = MXNTTS
*
      CALL MEMMAN(KOOS1,NOOS,'ADDL  ',2,'OOS1  ')
      CALL MEMMAN(KOOS2,NOOS,'ADDL  ',2,'OOS2  ')
      CALL MEMMAN(KOOS3,NOOS,'ADDL  ',2,'OOS3  ')
      CALL MEMMAN(KOOS4,NOOS,'ADDL  ',2,'OOS4  ')
      CALL MEMMAN(KOOS5,NOOS,'ADDL  ',2,'OOS5  ')
*. Space for four blocks of string occupations and arrays of
*. reordering arrays
      LZSCR = (MAX(NAEL,NBEL)+1)*(NOCOB+1) + 2 * NOCOB
      LZ    = (MAX(NAEL,NBEL)) * NOCOB
      DO I1234 = 1, 4
        CALL MEMMAN(KLOCSTR(I1234),MAX_STR_OC_BLK,'ADDL  ',1,'KLOCS ')
        CALL MEMMAN(KLREO(I1234),MAX_STR_SPGP,'ADDL  ',1,'KLREO ')
        CALL MEMMAN(KLZ(I1234),LZ,'ADDL  ',1,'KLZ   ')
      END DO
      CALL MEMMAN(KLZSCR,LZSCR,'ADDL  ',1,'KLZSCR')
*
      IF(I12.EQ.2) THEN
        IDOH2 = 1
      ELSE
        IDOH2 = 0
      END IF
*. Place perturbation integrals over one body integrals
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
*. Scale and reformat
      IF(IDC.EQ.2) THEN
*. scale
        CALL SCDTTS(CB,IRBLOCK(1,IROFF),NRBLOCK,NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
        CALL SCDTTS(HCB,ILBLOCK(1,ILOFF),NLBLOCK,NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &              IDC,2,NTEST)
*. reform
        CALL RFTTS(CB,WORK(KC2),IRBLOCK(1,IROFF),NRBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &             IDC,PSSIGN,2,NTEST)
        CALL RFTTS(HCB,WORK(KC2),ILBLOCK(1,ILOFF),NLBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &             IDC,PSSIGN,2,NTEST)
      END IF
*
**      CALL SBLOCKS2(NLBLOCK,ILBLOCK(1,ILOFF),NRBLOCK,IRBLOCK(1,IROFF),
**     &            CB,HCB,WORK(KC2),
**     &            WORK(KCIOIO),ISMOST(1,ICSM),WORK(KCBLTP),
**     &            NACOB,WORK(KNSTSO(IATP)),WORK(KNSTSO(IBTP)),
**     &            NAEL,IATP,NBEL,IBTP,
**     &            IOCTPA,IOCTPB,NOCTPA,NOCTPB,
**     &            NSMST,NSMOB,NSMSX,NSMDX,NOBPTS,IOBPTS,MXPNGAS,
**     &            ITSOB,MAXIJ,MAXK,MAXI,INSCR,LBLOCK,
**     &            LBLOCK,WORK(KINSCR),WORK(KCSCR),WORK(KSSCR),
**     &            SXSTSM,WORK(KSTSTS),WORK(KSTSTD),SXDXSX,
**     &            ADSXA,ASXAD,NGAS,NELFSPGP,IDC,
**     &            WORK(KOOS1),WORK(KOOS2),WORK(KOOS3),WORK(KOOS4),
**     &            WORK(KOOS5),WORK(KI1),WORK(KXI1S),
**     &            WORK(KI2),WORK(KXI2S),IDOH2,MXPOBS,WORK(KSVST),
**     &            PSSIGN,IPRDIA,LUC,ICJKAIB,WORK(KCJRES),
**     &            WORK(KSIRES),WORK(KI3),WORK(KXI3S),
**     &            WORK(KI4),WORK(KXI4S),MXSXST,MXSXBL,
**     &            MOCAA,MOCAB,IAPR,IRESTRICT)
*
      IF(IDC.EQ.2) THEN
*. reform
        CALL RFTTS(HCB,CB,ILBLOCK(1,ILOFF),NLBLOCK,
     &             1,NSMST,NOCTPA,NOCTPB,
     &             WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &             IDC,PSSIGN,1,NTEST)
*. scale
        CALL SCDTTS(HCB,ILBLOCK(1,ILOFF),NLBLOCK,NSMST,NOCTPA,NOCTPB,
     &              WORK(KNSTSO(IATP)), WORK(KNSTSO(IBTP)),
     &              IDC,1,NTEST)
      END IF

*. restore order
      IF(IPERTOP.NE.0) CALL SWAPVE(WORK(KINT1),WORK(KFI),NINT1)
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(IDUM ,IDUM,'FLUSM ',2,'SBFCB ')
      RETURN
      END
