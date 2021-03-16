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

      SUBROUTINE RECP_LNK_READCP(LUINFO,NUC1,NONT0,IPREAD)
!     Read ECP parameters
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
#include "mxcent.h"
#include "molinp.h"
#include "nuclei.h"
#include "recpval.h"
#include "inc_mxvalue.h"
#include "argoscom.h"
!     global variables
      INTEGER LUINFO,NUC1,NONT0,IPREAD
      REAL*8  Q,Q1
!     local variables
      INTEGER I,J,K,IOS,RECPIN1,LUBAS,INTQ0,ISTART,ECPINDEX
      LOGICAL ECPLIB
      CHARACTER*3 ECPKEYWORD,SOREP_CHR
      CHARACTER*30 ECPNAM,ECPTMP
      CHARACTER*80 ELINE

      RECPIN1 = 51
      ECPLIB = .FALSE.

!     Save NONT0 
      ECPNONT(NUC1) = NONT0
      ECPINDEX = 0

!     print *,'Save ECPNONT(',NUC1,')',NONT0

!     Check 'ECP' keyword
!     ===================
      CALL RDLINE(LUINFO)
      READ (MLINE(NMLINE),'(BN,A3)',IOSTAT=IOS) ECPKEYWORD
      IF (IOS.NE.0) GOTO 8000

      IF (ECPKEYWORD.NE.'ECP') THEN
!        Discard the line we just read.
         NMLINE = NMLINE - 1
         GOTO 9000
      ENDIF
      
!     Set RECP_SET/ Open 'RECP_INPUT2' file
!     =====================================
      IF (.NOT.ALLOCATED(RECP_SET)) THEN
         ALLOCATE(RECP_SET(4))
         RECP_SET(1) = 0   ! Number of RECP 
         RECP_SET(3) = 0   ! lines of RECP read  
         RECP_SET(4) = 0   ! last number of NUC1 for  
!        Allocate variables for RECP parameters
         CALL RECP_NTRA_READCP(MXCENT)
         RECP_CORE = 0
         AREP_ANG  = 0
         SOREP_ANG = 0 
      ENDIF

!     Check for the last RECP Nuc. center 
      RECP_SET(4) = MAX(RECP_SET(4),NUC1)

!     Check types of RECP read
!     ========================
      IF (MLINE(NMLINE)(4:6).EQ.'LIB') THEN

!        Read From library
!        -----------------
!        CALL GETENV ('BASDIR', BASDIR)
         ECPLIB = .TRUE.
         READ(MLINE(NMLINE),'(7X,A30)') ECPNAM
         ISTART = INDEX(ECPNAM,' ')
         ECPTMP = ECPNAM(1:ISTART - 1)   ! ecp name for search
         INTQ0  = NINT(CHARGE(ECPNUC+1)) ! charge index for search
         LUBAS  = 20                     ! file unit for open
         CALL RECP_LNK_FINDEL(ECPTMP,INTQ0,LUBAS)
         CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
         READ (ELINE,'(3X,I3,I2,I2)',IOSTAT=IOS) 
     &        RECP_CORE(NUC1), AREP_ANG(NUC1), SOREP_ANG(NUC1)
         IF (IOS.NE.0) GOTO 8000
!        CALL RECP_LNK_FINPOS(NUC1,CORE0,AREP0,SOREP0,LUBAS)
      ELSE 
!        Explicit RECP parameter
!        -----------------------
!        Read #CORE, #AREP
         AREP_ANG(NUC1) = 0
         READ (MLINE(NMLINE),'(3X,I3,I2)',IOSTAT=IOS)
     &        RECP_CORE(NUC1), AREP_ANG(NUC1)
         IF (IOS.NE.0) GOTO 8000
!        Read #SOREP
         READ (MLINE(NMLINE),'(8X,A3)',IOSTAT=IOS) SOREP_CHR
         IF (IOS.NE.0) GOTO 8000
         IF (SOREP_CHR.EQ.'   ') THEN  ! Skip SOREP (blank)
            SOREP_ANG(NUC1) = 0
         ELSE
            READ (MLINE(NMLINE),'(8X,I2)',IOSTAT=IOS)SOREP_ANG(NUC1)
            IF (IOS.NE.0) GOTO 8000
         ENDIF
      ENDIF

!     Count the number of RECP   
      RECP_SET(1) = RECP_SET(1) + 1

!     Reset charge value
!     ==================
      Q1 = CHARGE(ECPNUC+1)
!     print *,'ECPNUC',ECPNUC,'Q',Q1,'CHARGE',CHARGE(ECPNUC+1)
!     print *,'CORE',RECP_CORE(NUC1)
      DO I = 1, NONT0
         ECPNUC = ECPNUC + 1
         IF (ECPNUC.GT.MXCENT) GOTO 8300
         CHARGE(ECPNUC) = Q1 - RECP_CORE(NUC1)*1.0D0
         IF (CHARGE(ECPNUC).LT.0) GOTO 8310
!        GNUEXP(ECPNUC) = 6.66D+37
      ENDDO
      IF (IPREAD .GT. 0) THEN
         WRITE(LUPRI,'(2X,A,I3,A)') 
     &     '* This atomic center has RECP with ',RECP_CORE(NUC1),
     &     ' core electrons.'
         WRITE(LUPRI,'(2X,A,F8.3)') 
     &     '  The charge value is changed as',CHARGE(ECPNUC)
      END IF

!     Write RECP information into file
!     ================================
!     WRITE (RECPIN1,'(4(I3))') 
!    &      RECP_CORE(NUC1),(AREP_ANG(NUC1)-1),SOREP_ANG(NUC1),NUC1

      IF (.NOT.ALLOCATED(RECP_MCRS)) THEN 
         ALLOCATE(RECP_MCRS(msup)) ! maximum value of msup is set in "inc_mxvalue.h"
         DO I = 1, msup
            RECP_MCRS(I) = 0
         ENDDO
      ENDIF

!     RECP_MCRS(NUC1) = RECP_SET(1)
      DO I = 1, (NUC1-1)
         ECPINDEX = ECPINDEX + ECPNONT(I)  
      ENDDO
      DO I = (ECPINDEX+1), (ECPINDEX+ECPNONT(NUC1))
         RECP_MCRS(I) = RECP_SET(1) 
      ENDDO
!     ---------
!     Read AREP
!     ---------
      IF (AREP_ANG(NUC1).EQ.0) GOTO 9000
      DO I = 1, AREP_ANG(NUC1)
!        Read the number of element in each block 
         IF (ECPLIB) THEN
            CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
            READ (ELINE,*,IOSTAT=IOS) AREP_BLK(NUC1,I)
         ELSE
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS) AREP_BLK(NUC1,I)
         ENDIF
         IF (IOS.NE.0) GOTO 8000
!        WRITE(RECPIN1,*) AREP_BLK(NUC1,I)

!        Read in each block 
         DO J = 1, AREP_BLK(NUC1,I)
            RECP_SET(3) = RECP_SET(3) + 1
            IF (ECPLIB) THEN
               CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
               READ (ELINE,*,IOSTAT=IOS)
!    &              AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
     &              RECPIN_R(RECP_SET(3)),RECPIN_E(RECP_SET(3)),
     &              RECPIN_C(RECP_SET(3))
            ELSE
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),*,IOSTAT=IOS)
!    &              AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
     &              RECPIN_R(RECP_SET(3)),RECPIN_E(RECP_SET(3)),
     &              RECPIN_C(RECP_SET(3))
            ENDIF
!           WRITE(RECPIN1,*)
!    &           AREP_R(NUC1,I,J),AREP_E(NUC1,I,J),AREP_C(NUC1,I,J)
            IF (IOS.NE.0) GOTO 8000
         ENDDO
      ENDDO

!     ---------------------
!     Read SOREP (if exist)
!     ---------------------
      IF (SOREP_ANG(NUC1).EQ.0) GOTO 9000 
      DO I = 1, SOREP_ANG(NUC1)
!        Read in each block         
         SOREP_BLK(NUC1,I) = -1     ! initial value
         IF (ECPLIB) THEN
            CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
            READ (ELINE,*,IOSTAT=IOS) SOREP_BLK(NUC1,I)
         ELSE
            CALL RDLINE(LUINFO)
            READ (MLINE(NMLINE),*,IOSTAT=IOS) SOREP_BLK(NUC1,I)
         ENDIF
         IF (IOS.NE.0) GOTO 8000 
!        WRITE(RECPIN1,*) SOREP_BLK(NUC1,I)

!        Read in each block         
         DO J = 1, SOREP_BLK(NUC1,I)
            RECP_SET(3) = RECP_SET(3) + 1
            IF (ECPLIB) THEN
               CALL RECP_LNK_RDECPLINE(ELINE,LUBAS)
               READ (ELINE,*,IOSTAT=IOS)
!    &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
     &              RECPIN_R(RECP_SET(3)),RECPIN_E(RECP_SET(3)),
     &              RECPIN_C(RECP_SET(3))
            ELSE
               CALL RDLINE(LUINFO)
               READ (MLINE(NMLINE),*,IOSTAT=IOS)
!    &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
     &              RECPIN_R(RECP_SET(3)),RECPIN_E(RECP_SET(3)),
     &              RECPIN_C(RECP_SET(3))
            ENDIF
!           WRITE(RECPIN1,*)
!    &           SOREP_R(NUC1,I,J),SOREP_E(NUC1,I,J),SOREP_C(NUC1,I,J)
            IF (IOS.NE.0) GOTO 8000 
         ENDDO
      ENDDO
      GOTO 9000

! Error message
 8000 CONTINUE
      WRITE(LUPRI,*) 'Error in reading line'
      WRITE(LUPRI,*) 'Line',NMLINE,': ',MLINE(NMLINE)
      CALL QUIT('Error in reading line')
 8300 CONTINUE     
      WRITE(LUPRI,*) 'ECPNUC exceeded MXCENT'
      CALL QUIT('Error in reading line')
 8310 CONTINUE     
      WRITE(LUPRI,*) 'The nuclear charge value is less than 0'
      WRITE(LUPRI,*) 'Total nuclear charge   : ',Q1
      WRITE(LUPRI,*) 'Number of core elecron : ',RECP_CORE(NUC1)
      CALL QUIT('Wrong charge or core-electron number')
 9000 CONTINUE
      IF(ECPLIB) CLOSE(UNIT=LUBAS,STATUS='KEEP')
      END SUBROUTINE RECP_LNK_READCP

      SUBROUTINE RECP_LNK_CHECKMOL(DOCART)
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
      LOGICAL DOCART

      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

!     Check if Cartesian gaussians is used 
      IF (.NOT.DOCART) THEN
         WRITE(LUPRI,'(A)') ' >>> ERROR in RECP_LNK_CHECKMOL'
         WRITE(LUPRI,'(A)') ' ECP calculation is only allowed for the'
         WRITE(LUPRI,'(A)') ' Cartesian Gaussians in this version'
         CALL QUIT('Error in setting Cartesian Gaussians') 
      ENDIF

 9000 CONTINUE
      END

      SUBROUTINE RECP_LNK_IRREP(MAXREP)
      USE RECP_NTR
      IMPLICIT NONE
#include "priunit.h"
#include "argoscom.h"
#include "argoscomch.h"
#include "pgroup.h"
#include "inc_geninp.h"
      INTEGER I,J,RECPIN1,MAXREP,ND(8)
      RECPIN1 = 51

      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

      CALL RECP_NTRA_IRREP(MAXREP)

!     Point group symmetry
!     ====================
      RECP_MAXREP(1) = (MAXREP+1)
      DO I = 1, RECP_MAXREP(1)
         RECP_ITYP(I) = itypag(I-1)
      ENDDO

!     P1 X P2 = P3
!     ============
      RECP_NDPT(1) = ndptag 
      DO I = 0, MAXREP
         DO J = 0, MAXREP
            IF (itypag(I).EQ.IRTAB(IGRPAG*10 + J)) 
     &         RECP_IRREP(J+1) = (I+1)
         ENDDO
      ENDDO
 
 9000 CONTINUE
      END SUBROUTINE RECP_LNK_IRREP

      SUBROUTINE RECP_LNK_RDGEO(SYMTXT,NSYMOP,KASYM,blk,
     &           IQM,JCO,KATOM,KANG)
!     Write geometry data for RECP integral
      USE RECP_NTR
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "nuclei.h"
      INTEGER      KATOM,KANG
      DIMENSION    IQM(KATOM),JCO(KANG,KATOM)
!     INTEGER      IQM(KATOM),JCO(KANG,KATOM)
      CHARACTER*1  KASYM(3,3)
      CHARACTER*11 SYMTXT
      REAL*8       XYZ3(3,8,MAX_NONT,MXATOM)
      INTEGER      MJ2(MXATOM,MAX_NONT),ISEND(8,MXATOM,8,MXATOM),NFUNC
      INTEGER      I,J,K,N,ICENT,NSYMOP,IGENQ(3),NGENP,blk,JBLOCK,KBLOCK

      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

!     IF (SYMTXT(1:1).EQ.'D' ) GOTO 8000
!     IF (SYMTXT(1:1).EQ.'C' ) GOTO 8000

      CALL RECP_LNK_PATM(ISEND,XYZ3,MJ2,NSYMOP,KASYM)
      CALL RECP_LNK_PGEN(ISEND,NSYMOP,NGENP,MJ2)


!     Set geometry data
!     =================
      CALL RECP_NTRA_GEO(MXATOM,KATOM,KANG,NONTYP,NONT,
     &                   MJ2,IQM,JCO,NGENP)
      ICENT  = 0
      JBLOCK = 0
      KBLOCK = 0
      DO I = 1,NONTYP ! atomic type
         DO N = 1,NONT(I) ! symmetry-independent center
            NFUNC = 0
            ICENT = ICENT + 1
            DO J = 1, IQM(I)
               NFUNC = NFUNC + JCO(J,I)
            ENDDO

            RECPIN_NAMN(ICENT)     = NAMN(ICENT)      ! Name of Nuc.center
            RECPIN_BASISBLK(ICENT) = NFUNC
            RECPIN_NUC2(ICENT)     = (MJ2(I,N)-1)
            RECPIN_CHARGE(ICENT)   = CHARGE(ICENT)

            DO J = 1, MJ2(I,N)-1
!              RECPIN_NAMN(J,ICENT)= NAMN(ICENT)
               RECPIN_GEO(1,J,ICENT) = XYZ3(1,J,N,I)  ! X
               RECPIN_GEO(2,J,ICENT) = XYZ3(2,J,N,I)  ! Y
               RECPIN_GEO(3,J,ICENT) = XYZ3(3,J,N,I)  ! Z
            ENDDO

!           Write interchange
            IGENQ(1)=2
            IGENQ(2)=3
            IGENQ(3)=5
            IF (MJ2(I,N).GT.2) THEN
               DO J=1,NGENP
                  DO K = 1,(MJ2(I,N)-1)
                     RECPIN_ICA(K,ICENT,J) = ISEND(K,N,IGENQ(J),I)
                  ENDDO
               ENDDO
            ENDIF

!           Match basis-block and AO2SO-matrix
            DO J = 1, IQM(I)
               JBLOCK = JBLOCK + 1
               DO K = 1, JCO(J,I)
                  KBLOCK = KBLOCK + 1
                  RECPIN_BLKMATCH(1,KBLOCK)= RECPIN_NCONS0(K,J,I) 
                  RECPIN_BLKMATCH(2,KBLOCK)= JBLOCK 
!                 print *, 'RECPIN_BLKMATCH',KBLOCK,
!    &              RECPIN_BLKMATCH(1,KBLOCK),RECPIN_BLKMATCH(2,KBLOCK)
               ENDDO
            ENDDO

!           Index for mcrs
!           WRITE(RECPIN1,'(I3)') RECP_MCRS(I)
         ENDDO
      ENDDO
!     DEALLOCATE ( RECPIN_NCONS0 ) 

!     General options
!     ===============

!     RECP_SET(1),RECP_SET(2)
!     Number of RECP parameters : nmcrs - RECP_SET(1)
!     Sum of all JCOs           : ncons - RECP_SET(2)

      RECPIN_GENVAL(1) = NGENP 
      RECPIN_GENVAL(2) = NUCIND
      RECPIN_GENVAL(3) = blk 
      GOTO 9000

!     Error message
!     =============
!8000 CONTINUE
!     WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
!     WRITE(LUPRI,'(X,A,A)') 'RECP integral is not allowed to ',
!    &                       'calculate linear symmetries.'
!     CALL QUIT('RECP_LNK: No linear symmetries are allowed')
 9000 CONTINUE
      END


      SUBROUTINE RECP_LNK_PATM(ISEND,XYZ3,MJ2,NSYMOP,KASYM)
      IMPLICIT NONE
#include "mxcent.h"
#include "nuclei.h"
      CHARACTER*1 KASYM(3,3)
      INTEGER   LGEN(3,3),NIREP,LTRANS(3,8),ICLAT,IREP1,ICENT
      REAL*8    DISTSQ,NCLOSE,sqdif,tolsq
      INTEGER   I,J,K,J2,J3,NSYMOP,  MJ2(MXATOM,MAX_NONT),
     &          ISEND(8,MXATOM,8,MXATOM)
      REAL*8    XYZ1(3,8,MXATOM), XYZ2(3,8,MXATOM,MAX_NONT),
     &          XYZ3(3,8,MAX_NONT,MXATOM),XYZ4(3)

      tolsq=1.0d-12

!     read rotation 
!     -------------      
      DO I = 1, NSYMOP
         DO J = 1, 3
            LGEN(J,I) = 1
         ENDDO
         DO J = 1, 3
            IF(KASYM(J,I).EQ.'X') LGEN(1,I)=-1
            IF(KASYM(J,I).EQ.'Y') LGEN(2,I)=-1
            IF(KASYM(J,I).EQ.'Z') LGEN(3,I)=-1
         ENDDO
!        print *,'LGEN',LGEN(1,i),LGEN(2,i),LGEN(3,i)
      ENDDO
 
!     generate symetry point group operations
!     ---------------------------------------
      DO j=1,3
         LTRANS(j,1)=1
         IF (NSYMOP.GE.1) THEN
            LTRANS(j,2)=LGEN(j,1)
            IF (NSYMOP.GE.2) THEN
               LTRANS(j,3)=LGEN(j,2)
               LTRANS(j,4)=LGEN(j,2)*LGEN(j,1)
               IF (NSYMOP.GE.3) THEN
                  LTRANS(j,5)=LGEN(j,3)
                  LTRANS(j,6)=LGEN(j,3)*LGEN(j,1)
                  LTRANS(j,7)=LGEN(j,3)*LGEN(j,2)
                  LTRANS(j,8)=LGEN(j,3)*LGEN(j,2)*LGEN(j,1)
               ENDIF
            ENDIF
         ENDIF
      ENDDO
 
!     set coordinate
!     --------------
      ICENT=1
      DO I=1,NONTYP
         DO J=1,NONT(I)
            DO K=1,3
               XYZ1(K,J,I)=CORD(k,ICENT) ! CORD from DIRAC
               IF(NSYMOP.eq.0)THEN
                 XYZ3(k,1,j,i)=XYZ1(K,J,I)
                 MJ2(i,j)=2
               ENDIF
            ENDDO
            ICENT=ICENT+1
         ENDDO
      ENDDO
      IF(NSYMOP.eq.0) GOTO 999
 
!     operation on to origin coordinate
!     ---------------------------------
      NIREP = 2 ** NSYMOP
 
      do i = 1, NONTYP 
         do j = 1, NONT(I)   
            MJ2(i,j)=1       
            do j2=1,NIREP
               do k=1,3
                  XYZ2(k,j2,j,i)=LTRANS(k,j2)*XYZ1(k,j,i)
                  XYZ4(k)=XYZ2(k,j2,j,i)
               enddo
!              check 
               do J3=1,j2-1
                  DISTSQ=sqdif(XYZ4(1),XYZ2(1,J3,j,i),3)
                  if (DISTSQ.lt.tolsq) then
                     goto 50
                  endif
               enddo
               do k=1,3
                  XYZ3(k,MJ2(i,j),j,i)=XYZ2(k,j2,j,i)
               enddo
               MJ2(i,j)=MJ2(i,j)+1
 50            continue
            enddo
         enddo
      enddo
 
!     operation  
!     ---------
      do i=1,NONTYP
         do IREP1=1,NIREP
            do j=1,NONT(i)
               do j2=1,MJ2(i,j)-1   ! natom(nc)
                  do k=1,3
                     XYZ4(k)=LTRANS(k,IREP1)*XYZ3(k,j2,j,i)
                  enddo
                  
                  NCLOSE=10000      ! closest atom
                  do j3=1,MJ2(i,j)-1
                     DISTSQ=sqdif(XYZ4(1),XYZ3(1,j3,j,i),3)
                     if (DISTSQ.lt.NCLOSE) then
                        NCLOSE=DISTSQ
                        ICLAT=j3    ! closest atom so far
                     endif
                  enddo
                  ISEND(j2,j,IREP1,i)=ICLAT
               enddo
            enddo
         enddo
      enddo
 999  continue
      END

      FUNCTION sqdif( a, b, n )
!     returns the sum of the squares of the differences between
!     the elements of the vectors a(:) and b(:)
      IMPLICIT NONE
      INTEGER n,it
      REAL*8 a(n), b(n), sqdif
      sqdif = 0d0
      DO it = 1, n
         sqdif = sqdif + (a(it) - b(it))**2
      ENDDO
      RETURN
      END

      SUBROUTINE RECP_LNK_PGEN(ISEND,NSYMOP,NGENP,MJ2)
      IMPLICIT NONE
#include "mxcent.h"
#include "nuclei.h"
      INTEGER   NGENP,IGENP(8),IGENQ(3),NSYMOP
      INTEGER   MJ2(MXATOM,MAX_NONT), ISEND(8,MXATOM,8,MXATOM)
      INTEGER   I,J,IGEN,JGEN,KGEN,IATOM,IATOM2

      NGENP=0
      IGENQ(1)=2
      IGENQ(2)=3
      IGENQ(3)=5
      DO JGEN=1,NSYMOP
         IGEN=IGENQ(JGEN)
         DO I=1,NUCIND         ! NUCIND=nsag
            DO J=1,NONT(I)
               DO IATOM=1,MJ2(I,J)-1          !MJ2(I,J)-1 = nc(I)
                  IF (ISEND(IATOM,J,IGEN,I).NE.IATOM) GOTO 550
               ENDDO
            ENDDO
         ENDDO
         GOTO 559

 550     CONTINUE               ! this is not the identity
         DO KGEN=1,NGENP
            DO I=1,NUCIND
               DO J=1,NONT(I)
                  DO IATOM=1,MJ2(I,J)-1
                     IF (ISEND(IATOM,J,IGEN       ,I) .NE.
     &                   ISEND(IATOM,J,IGENP(KGEN),I))
     &               GOTO 551
                  ENDDO
               ENDDO
            ENDDO
            GOTO 559
 551        CONTINUE            ! this is not generator KGEN
         ENDDO

         IF (NGENP.GE.2) THEN
            DO I=1,NUCIND
               DO J=1,NONT(I)
                  DO IATOM=1,MJ2(I,J)-1
                     IATOM2 = ISEND(IATOM,J,IGENP(2),I)
                     IF (ISEND(IATOM ,J,IGEN    ,I) .NE.
     &                   ISEND(IATOM2,J,IGENP(1),I)) 
     &                  GOTO 552
                  ENDDO
               ENDDO
            ENDDO
            GOTO 559
         ENDIF

 552     CONTINUE     ! this is not a product of previous generators
         NGENP=NGENP+1
         IGENP(NGENP)=IGENQ(JGEN)
 559     CONTINUE
      ENDDO
      END

      SUBROUTINE RECP_LNK_RDORB(IQM,JCO,NUC,NRC,
     &           ALPHA,CPRIM,CPRIMU,KATOM,KANG,KBLOCK,KPRIM)
!     Write Basis-set information for RECP integral
!     From SUBROUTINE ORBPRO 
      USE RECP_NTR
#include "implicit.h"
#include "mxcent.h"
#include "nuclei.h"
#include "dcbgen.h"
      INTEGER   KATOM,KANG,KBLOCK,KPRIM
      DIMENSION IQM(KATOM),
     &          JCO(KANG,KATOM),NUC(KBLOCK),NRC(KBLOCK),
     &          ALPHA(KPRIM,KBLOCK),
     &          CPRIM(KPRIM,KPRIM,KBLOCK),CPRIMU(KPRIM,KPRIM,KBLOCK)
      INTEGER I,J,K,L,M,N,RECPIN1,KBCH

!     Skip if 'RECP_SET' is not setted
      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000

!     Allocate variables
!     ==================
      CALL RECP_NTRA_BASIS(KATOM,KANG,KBLOCK,NONTYP,IQM,JCO,NUC,NRC)

!     Set BASIS variables
!     ===================
      KBCH = 0
      RECP_SET(2) = 0
      
      RECPIN_NONTYP(1) = NONTYP
      DO I = 1,NONTYP          ! atomic type
         RECPIN_IQM(I) = IQM(I)
         DO J = 1,IQM(I)       ! angular momentum value
            RECPIN_JCO(J,I) = JCO(J,I)
!           print *,'RECPIN_JCO(',J,',',I,')',RECPIN_JCO(J,I)

            DO K = 1, JCO(J,I) ! blocks of given atom and L-value
               KBCH = KBCH + 1
               RECPIN_NCONS0(K,J,I) = KBCH
!              ----------------------------------------------------
!              Set numbers to read in each block
!              =================================
               RECPIN_NUC(KBCH) = NUC(KBCH)
               RECPIN_NRC(KBCH) = NRC(KBCH)
               RECPIN_LMNP(KBCH)= J
!              print *,NUC(KBCH), K, NRC(KBCH)

               DO L = 1, NUC(KBCH)
                  RECPIN_ALPHA(L,KBCH) = ALPHA(L,KBCH)
                  DO M = 1, NRC(KBCH)
                     RECPIN_CPRIMU(L,M,KBCH) = CPRIMU(L,M,KBCH)
                  ENDDO
!                 print*,ALPHA(L,KBCH),
!    &                  (CPRIMU(L,M,KBCH),M=1,NRC(KBCH))
               ENDDO
!              ----------------------------------------------------
            ENDDO
         ENDDO
      ENDDO
      RECP_SET(2) = KBCH

 9000 CONTINUE
      END


      SUBROUTINE RECP_LNK_CTRAN(blk)
      USE RECP_NTR
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "symmet.h"
#include "pincom.h"
#include "aosotr.h"
#include "sphtrm.h"
#include "huckel.h"
#include "cbirea.h"
#include "argoscom.h"
#include "argoscomch.h"
#include "ibtfun.h"
      INTEGER blk
!     Skip if RECP is not set
      IF (.NOT.ALLOCATED(RECP_SET)) GOTO 9000
!     Allocate CTRANARG
      CALL RECP_NTRA_AOTOSO1

      blk = 0
      DO ICENT = 1, NUCIND   ! run over nuclear center
         DO IA = 1, KMAX     ! run over shells
!           Set irowmax
            irowmax(ICENT,NHKT(IA)) = nucdeg(ICENT) * KHKT(IA)
            IF (NCENT(IA).EQ.ICENT) THEN
               IF ((IA-1).EQ.0) THEN  !(IA-1) is bigger than 0
                  CALL RECP_LNK_CTRAN1(ICENT,IA,blk)
               ELSEIF (NCENT(IA).NE.NCENT(IA-1)) THEN
                  CALL RECP_LNK_CTRAN1(ICENT,IA,blk)
               ELSEIF (NHKT(IA).NE.NHKT(IA-1)) THEN
                  CALL RECP_LNK_CTRAN1(ICENT,IA,blk)
               ENDIF
            ENDIF
         ENDDO
      ENDDO   
!     Set AO to SO transformation matrices
      CALL RECP_LNK_CTRAN2
      DEALLOCATE (CTRANARG)
 9000 CONTINUE
      END 


      SUBROUTINE RECP_LNK_CTRAN1(ICENT,IA,blk)
      USE RECP_NTR
#include "implicit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "shells.h"
#include "symmet.h"
#include "argoscom.h"
#include "ibtfun.h"
      INTEGER ICENT,IA,blk,ila,NA,LA,IVARB,JKB,KB,irow,icol
      blk = blk +1
      DO NA = 1, KHKT(IA)   ! run over components
         ila = 0
         DO LA = 0, MAXREP  ! run over symmetry operator
            IVARB = IBTXOR(LA,ISYMAO(NHKT(IA),NA))
!           If orbital contributes to this representation:
!           See IVARB in SYMPRO (herrdn.F)
            IF (IBTAND(ISTBAO(IA),IVARB) .EQ. 0) THEN
               ila = ila + 1
               JKB = 0
!              loop over symmetry dependent centers
               DO KB = 0, MAXREP
                  IF (IBTAND(KB,ISTBAO(IA)) .EQ. 0) THEN
                     JKB  = JKB+ 1
                     irow = NA + KHKT(IA)*(ila-1)
                     icol = NA + KHKT(IA)*(JKB-1)
                     CTRANARG(ICENT,NHKT(IA),irow,icol)
     &                = PT(IBTAND(KB,IVARB))
                     ireparg(ICENT,NHKT(IA),irow) = LA + 1
                     MAXANG(ICENT) = NHKT(IA)
                  END IF
               ENDDO
            END IF
         ENDDO
      ENDDO
      END SUBROUTINE RECP_LNK_CTRAN1

      SUBROUTINE RECP_LNK_CTRAN2
!     Set AO to SO transformation matrices
      USE RECP_NTR
      IMPLICIT NONE
#include "mxcent.h"
#include "nuclei.h"
#include "argoscom.h"
      INTEGER I,J,K,L,AORDS1,IRREP1
      CALL RECP_NTRA_AOTOSO(NUCIND,MAXANG,irowmax)
      AORDS1 = 0
      DO I = 1, NUCIND
         DO J = 1, MAXANG(I)
            AORDS1 = AORDS1 + 1
!           Set list of irrep. in i-th AO reduction set
            RECPIN_MXROW(AORDS1) = irowmax(I,J)
            DO IRREP1 = 1,irowmax(I,J)
               RECPIN_LA1(IRREP1,AORDS1) = ireparg(I,J,IRREP1)
            ENDDO
!           Set AO-SO matrices
            DO K = 1, irowmax(I,J)
               DO L = 1, irowmax(I,J)
                  RECPIN_AO2SO(L,K,AORDS1)=nint(CTRANARG(I,J,K,L)) 
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      END SUBROUTINE RECP_LNK_CTRAN2


      SUBROUTINE RECP_LNK_FINDEL(BASNAM,INTQ0,LUBAS)
!     From FINDEL
#include "implicit.h"
#include "priunit.h"
#include "dummy.h"
C*********************************************************
C* Variable declarations:                                *
C*     BASNAM - The name of the basis file               *
C*     STRING - A character variable that helps bullet-  *
C*              proofing the subroutine.                 *
C*     EXST - Logical variable that helps inquire if     *
C*            there exists a file with that name.        *
C*     SLASH - variable for slash in paths (/ or \\)     *
C*     WIN - Check if we work on Linux or on Windows     *
C*********************************************************
      CHARACTER*(*) BASNAM
C     long because contains the full dir-path
      CHARACTER*200 STRING
      CHARACTER SIGN
      LOGICAL EXST, EXST2, USE_BASNAM
      INTEGER INTQ0
      CHARACTER SLASH
      LOGICAL WIN
#include "inftap.h"
#if defined (PRG_DIRAC)
#include "dcbgen.h"
#else
#include "gnrinf.h"
#endif
#if defined (SYS_WINDOWS)
      WIN = .TRUE.
#else
      WIN = .FALSE.
#endif
C
C     Linux:
C     The environment variable is a : separated string, e.g.:
C     /programs/Dirac/basis:/programs/Dalton/basis:.
C
C     Windows:
C     The environment variable is a ; separated string, e.g.:
C     C:\programs\Dirac\basis;C:\programs\Dalton\basis;.
C     
C     We loop over these directories in search for the basis set and
c     element given.
C
C     Depending on the machine, the string has trailing blanks
C     (ASCII char=32) or null elements (ASCII=0).
      ISTART = 1
      IBLANK=INDEX(BASDIR,CHAR(32))
      INULL =INDEX(BASDIR,CHAR(0))
      ISTOP=MAX(IBLANK,INULL)
      IF(ISTOP.EQ.0) THEN
C       If string is empty add current directory 
        BASDIR(1:1) = '.'
        ISTOP = 1
      ELSE
        ISTOP=ISTOP-1
      ENDIF
C     
 10   CONTINUE
C     
C     Find next colon (: or ;) in BASDIR
      IF ( .NOT. WIN ) THEN
          ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:ISTOP), ':' )
      ELSE
          ICOLON = ISTART - 1 + INDEX (BASDIR(ISTART:ISTOP), ';' )
      END IF
C     
C     In there are no colons (: or ;), we are the last directory in BASDIR
      IF ( ICOLON .EQ. ISTART - 1 ) ICOLON = ISTOP + 1
      IEND = ICOLON - 1
      STRING = BASDIR(ISTART:IEND)//' '
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
!     stefan: make check with idb indicates a possible access penalty if
!     BASDIR is not set...
      if(iend.lt.1)then
        write(lupri,'(/a)') '  *** error in FINDEL: environment '//
     &        'variable BASDIR is not set'
        call quit('*** error in FINDEL: BASDIR is not set')
      end if
      IF ( STRING(IEND:IEND) .NE. SLASH ) THEN
         IEND = IEND + 1
         STRING(IEND:IEND) = SLASH
      END IF
      STRING = STRING(1:IEND)//BASNAM
      IEND = INDEX(STRING(1:),' ') - 1
Cluuk Activate when debugging basis set search (should not be on by
Cdefault)
CLuuk WRITE(LUPRI,'(A/,A1,A,A1)')
CLuuk&     ' Searching basis : trying file: ','"',STRING(1:IEND),'"'
C
C     Inquire if "/path/basis-set" exists.
C
      INQUIRE (FILE = STRING, EXIST = EXST)
CMI  ... on certain systems (like PC-linux,OpenMPI-gfortran) string
C     beginning with '/' is NOT the proper file name ! therefore
C     continue in inquiring:
      EXST2=EXST
      USE_BASNAM=.FALSE.
      IF (.NOT.EXST2) THEN
        INQUIRE (FILE = BASNAM, EXIST = EXST)
        IF (EXST) THEN
          WRITE(LUPRI,"(/2A/2A)")
     &     "Info: used proper ECP name : ",BASNAM,
     &     "Info: instead of           : ",STRING
          USE_BASNAM=.TRUE.
        ENDIF
      ENDIF
      IF (EXST) THEN
         IF (.NOT.USE_BASNAM) THEN
            OPEN(LUBAS, FILE=STRING)
         ELSE
            OPEN(LUBAS, FILE=BASNAM)
         ENDIF
C     
C
C        Searching the file for the element.
C
         REWIND(LUBAS)  
 20      CONTINUE
         READ(LUBAS,'(A88)', IOSTAT = IOERR, ERR = 2000, END = 200)
     +       STRING
         READ (STRING, '(A1)',IOSTAT=IOS) SIGN
         IF ((SIGN .EQ. 'a') .OR. (SIGN .EQ. 'A')) THEN
            READ (STRING, 110,IOSTAT=IOS) SIGN, NUCEL
            IF (IOS.NE.0) THEN
               CALL QUIT('FINDEL: Error in reading NUCEL')
            ENDIF
 110        FORMAT (BN, A1, I4)
            IF (INTQ0.EQ. NUCEL) THEN
               RETURN
            ELSE
               GOTO 20
            END IF
         ELSE
            GOTO 20
         END IF
      END IF
      ISTART = ICOLON + 1
      IF(ISTART.GT.ISTOP) GOTO 300
      GOTO 10
 300  CONTINUE
      LSTR = 1
      DO J = LEN(BASNAM), 1, -1
         IF ( BASNAM(J:J) .NE. ' ' ) THEN
            LSTR = J
            GOTO 301
         END IF
      END DO
 301  CONTINUE
      WRITE (LUPRI,'(/,3A)')
     +     'Basis "',BASNAM(1:LSTR),'" doesn''t exist'
      WRITE(LUPRI,'(2X,A,A)') 'ECP name, STRING=',STRING
      WRITE(LUPRI,'(2X,A,A)') 'ECP name, BASNAM=',BASNAM
      IF (BASNAM(1:6) .EQ. 'Huckel') THEN
         WRITE (LUPRI,'(A)')
     +        'NB! Basis Huckel must be available to use Huckel guess!'
      END IF
      CALL QUIT('Non-existing ECP')
C
C Error messages 
C
 2000 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(2A/A,I5)')
     &   'FINDEL: Error when reading from ECP file ',STRING(1:IEND),
     &   ' IOSTAT =',IOERR
      CALL QUIT('I/O error in FINDEL in HERBAS')
 200  CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(I3,2A)') INTQ0,
     &     ' is an unsupported element in',BASNAM
      CALL QUIT('Unsupported ECP element')
      END

      SUBROUTINE RECP_LNK_FINPOS(NUC1,CORE0,AREP0,SOREP0,LUBAS)
!     From FINPOS
      IMPLICIT NONE
#include "priunit.h"
      INTEGER NUC1,CORE0,AREP0,SOREP0,LUBAS
      INTEGER NUC0,IOS
      CHARACTER SIGN
      CHARACTER*88 STRING

!     Read / Find nuclear charge index
 1000 CONTINUE
      READ (LUBAS,'(A88)',ERR=8100,END=8110) STRING
      READ (STRING,'(A1)') SIGN
      IF ((SIGN.EQ.'a') .OR. (SIGN.EQ.'A')) THEN
         READ (STRING,'(2X,I3)') NUC0
         IF (NUC0.NE.NUC1) GOTO 1000
!        Read / Skip if #/!/$ 
 1100    CONTINUE
         READ (LUBAS,'(A88)',ERR=8100,END=8110) STRING
         READ (STRING,'(A1)') SIGN
         IF ((SIGN.EQ.'#').OR.(SIGN.EQ.'!').OR.(SIGN.EQ.'$')) GOTO 1100
!        Read / Find CORE0,AREP0,SOREP0
         READ (STRING,'(3X,I3,I2,I2)',IOSTAT=IOS) CORE0,AREP0,SOREP0 
         IF (IOS.NE.0) GOTO 8100
         GOTO 9000
      ELSE
         GOTO 1000
      ENDIF
 
 8100 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(A)') 'Error in reading following line'
      WRITE(LUPRI,'(A)') STRING
      CALL QUIT('Error in reading ECP file')
 8110 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE (LUPRI,'(/A)') 'Error in reading file, your ECP '//
     &     'is not complete.'
      CALL QUIT('Incomplete ECP file')

 9000 CONTINUE
      END

      SUBROUTINE RECP_LNK_RDECPLINE(ELINE,LUBAS)
      IMPLICIT NONE
#include "priunit.h"
      INTEGER LUBAS
      CHARACTER SIGN
      CHARACTER*80 ELINE

 1000 CONTINUE
      READ (LUBAS,'(A80)',ERR=8100,END=8110) ELINE
      READ (ELINE,'(A1)') SIGN 
      IF ((SIGN.EQ.'#').OR.(SIGN.EQ.'!').OR.(SIGN.EQ.'$')) GOTO 1000
      GOTO 9000
!     Error message
 8100 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(A)') 'Error in reading following line'
      WRITE(LUPRI,'(A)') ELINE
      CALL QUIT('Error in reading ECP file')
 8110 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(A)') 'Error in reading file, your ECP '//
     &     'is not complete.'
      CALL QUIT('Incomplete ECP file')
 9000 CONTINUE
      END


      SUBROUTINE RECP_LNK_RDMENU(WORD0)
! Read menu file
      USE RECP_OUTPUT
      IMPLICIT NONE
#include "priunit.h"
#include "argoscom.h"
      INTEGER I,IOS,NTABLE
      PARAMETER (NTABLE=3)
      CHARACTER WORD0*7,WORD*7,PROMPT*1,TABLE(NTABLE)*7
      DATA TABLE /'.PRINT ','.DEBUG ','XXXXXXX'/

! Original parameter
      ECPINT_PRINT = 1
      ECPINT_DEBUG = 0

! Check options in *ECP
! ---------------------
 2000 CONTINUE
      READ (LUCMD,'(A7)',END=9900,ERR=9000) WORD
      CALL UPCASE(WORD)
      IF (WORD(1:1) .EQ. '!') THEN
         GOTO 2000
      ELSEIF (WORD(1:1) .EQ. '#') THEN
         GOTO 2000
      ELSEIF (WORD(1:1) .EQ. '*') THEN
         WORD0 = WORD
         GOTO 9900
      ELSEIF (WORD(1:1) .EQ. '.') THEN
         DO I = 1, NTABLE
            IF (WORD.EQ.TABLE(I)) THEN
               GOTO(3010,3020,9100),I
            ENDIF
         ENDDO
      ELSE
         GOTO 9100
      ENDIF

! Reading vaues : .PRINT  
! ----------------------
 3010 CONTINUE
      READ (LUCMD,'(A7)',ERR=9000) WORD
      CALL UPCASE(WORD)
      IF (WORD(1:1).EQ.'!') GOTO 3010
      IF (WORD(1:1).EQ.'#') GOTO 3010
      READ (WORD,*,ERR=3011) ECPINT_PRINT
      WRITE (LUPRI,'(/,A,I4)') ' * Print level in ECP : ',ECPINT_PRINT
      GOTO 2000
!     Error message
 3011 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(A,I7)') ' .PRINT value error : ',ECPINT_PRINT
      GOTO 9000

! Reading vaues : .DEBUG 
! ----------------------
 3020 CONTINUE
      READ (LUCMD,'(A7)',ERR=9000) WORD
      CALL UPCASE(WORD)
      IF (WORD(1:1).EQ.'!') GOTO 3020
      IF (WORD(1:1).EQ.'#') GOTO 3020
      READ (WORD,*,ERR=3021) ECPINT_DEBUG
      WRITE (LUPRI,'(/,A)') ' * Debug level in ECP : ',ECPINT_PRINT
      GOTO 2000
!     Error message
 3021 CONTINUE
      WRITE(LUPRI,'(/,A,I7)') ' .DEBUG value error : ',ECPINT_DEBUG
      GOTO 9000
     a
! Error message
! -------------
 9000 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(2A)') ' Error in reading menu file : ',WORD
      CALL QUIT('Error in reading line')
 9100 CONTINUE
      WRITE(LUPRI,'(/,X,A)') '=== RECP Error ==='
      WRITE(LUPRI,'(3A)') ' Keyword ',WORD,
     &      ' is not recognized in ECP'
      CALL QUIT('Wrong input')
 9900 CONTINUE
      RETURN
      END

