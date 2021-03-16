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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
      subroutine grtest()
C                                                                      C
C     T. Helgaker sep 99                                               C
C     DIRAC adaption by T.Saue                                         C
C     Modified for grid testing by Olav Fossgaard July 2002.           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      use dft_cfg
      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
 
      PARAMETER (DP5 = 0.5D0)

#include "dcbgen.h"
#include "dcbham.h"
#include "inforb.h"
#include "nuclei.h"

      real(8), allocatable :: smat(:)
      real(8), allocatable :: gao(:)
      real(8), allocatable :: gao1(:)
      real(8), allocatable :: gao2(:)
      real(8), allocatable :: gab1(:)
!     real(8), allocatable :: ncnt(:) stefan: FIXME - as far as i can see this array should be integer
      integer, allocatable :: ncnt(:)
      real(8), allocatable :: buf(:)
      integer, allocatable :: iclab(:)

      lunit_grid = 44

      call alloc(smat,  nbast*nbast*2)
      call alloc(gao,   nbast)
      call alloc(gao1,  nbast*3)
      call alloc(gao2,  nbast*6)
      call alloc(gab1,  nbast*3)
      call alloc(ncnt,  nbast)
      call alloc(buf,   nbast)
      call alloc(iclab, nbast)

      call grtes1(smat,
     &            gao,
     &            gao1,
     &            gao2,
     &            gab1,
     &            ncnt,
     &            buf,
     &            iclab,
     &            iprham)

      call dealloc(smat)
      call dealloc(gao)
      call dealloc(gao1)
      call dealloc(gao2)
      call dealloc(gab1)
      call dealloc(ncnt)
      call dealloc(buf)
      call dealloc(iclab)

      end subroutine

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C /* Deck grtes1 */
      SUBROUTINE GRTES1(SMAT,GAO,GAO1,GAO2,GAB1,
     &                  NCNT,BUF,ICLAB,IPRINT)
C                                                                      C
C     T. Helgaker sep 99                                               C
C     Modified for grid testing by Olav Fossgaard July 2002.           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

         use dft_cfg
      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)

      LOGICAL SSOVER,PPOVER,DDOVER,FFOVER

      CHARACTER*12 MXLB1,MXLB1S,MXLB1P,MXLB1D,MXLB1F,MXLB2,MXLB2S,
     &     MXLB2P,MXLB2D,MXLB2F,RMXLB1,RMXLB2
      CHARACTER*1 ORB1,ORB2
      INTEGER CENT1,CENT2,MAXI,MAXIS,MAXIP,MAXID,MAXIF,MAXJ,MAXJS,
     &     MAXJP,MAXJD,MAXJF,MAXIR,MAXJR,NLSOS,NSSOS        
      DOUBLE PRECISION ERRAV,ERRAB,RMS,ERRMX,MAXS,MAXP,MAXD,MAXF,
     &     RERRAV,RERRAB,RERRMX,SMATDF, SMATAN, SMATFR

#include "dcbgen.h"
#include "dcbham.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "maxorb.h" 
#include "dcblab.h"
#include "maxaqn.h" 
#include "symmet.h"
#include "nuclei.h"
#include "orgcom.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))
      DIMENSION SMAT(NBAST,NBAST,2),ICLAB(NBAST)
#else
      DIMENSION SMAT(NTBAS(0),NTBAS(0),2),ICLAB(NTBAS(0))
#endif

      DIMENSION GAO(*), GAO1(*), 
     &          GAO2(*), GAB1(*),
     &          NCNT(*), BUF(*)
      real(8), allocatable :: kx(:)
      real(8), allocatable :: ky(:)
      real(8), allocatable :: kz(:)
      real(8), allocatable :: kw(:)

#include "dcbibt.h"



C
C     Map basis functions onto atomic centers in array ICLAB
C

C     Large components:
      I = 0
      DO IFRP = 1,NFSYM
         IOFF = IBAS(IFRP)
         DO ILBAS  = IOFF+1,IOFF+NFBAS(IFRP,1)
            I = I + 1
            IOLAV  = IATTR(IPLAB(ILBAS,2),2)
            ICENT  = JGET(IOLAV)
            ICLAB(I) = ICENT
         ENDDO
      ENDDO

C     Small components
      J = I
      DO IFRP = 1,NFSYM
         IOFF = IBAS(IFRP) + NFBAS(IFRP,1)
         DO ISBAS  = IOFF+1,IOFF+NFBAS(IFRP,2)
            J = J + 1
            IOLAV  = IATTR(IPLAB(ISBAS,2),2)
            ICENT  = JGET(IOLAV)
            ICLAB(J) = ICENT
         ENDDO   
      ENDDO   

C      DO I = 1,NTBAS(0)
C         WRITE(6,*) 'ICLAB(',I,') = ',ICLAB(I)
C      ENDDO

C
C     Get analytical SMAT(I,J,2)
C
      NELM = NTBAS(0)*NTBAS(0)
      smat = 0.0d0
      CALL GTOVLX(SMAT(1,1,2),SSMTRC)
C
C     integrate
C
      lunit_grid = 44
      REWIND lunit_grid
      NPNTS = 0
  100 CONTINUE
      READ(lunit_grid) NPOINT
      IF (NPOINT.GT.0) THEN
         NPNTS = NPNTS + NPOINT
         call alloc(kx, npoint)
         call alloc(ky, npoint)
         call alloc(kz, npoint)
         call alloc(kw, npoint)
         CALL GRTES3(kx, ky, kz, kw,
     &               lunit_grid,NPOINT)

      CALL GRTES2(kx, ky, kz, kw,
     &               NPOINT,
     &               SMAT,
     &               GAO,GAO1,GAO2,GAB1,
     &               NCNT,BUF)

         call dealloc(kx)
         call dealloc(ky)
         call dealloc(kz)
         call dealloc(kw)
         GO TO 100
      ELSE IF (NPOINT .EQ.0 ) THEN
         GO TO 100
      ELSE
         GO TO 200
      END IF
  200 CONTINUE

C
C     Subtract the analytical SMAT(I,J,2) from the numerical SMAT(I,J,1)
C     The numerical SMAT is overwritten.
C
      CALL DAXPY(NELM,-1.0D0,SMAT(1,1,2),1,SMAT(1,1,1),1)
C
C     Analyze the difference matrix:
C
C     ERRAV = average error
C     ERRAB = absolute average error
C     RMS   = root mean square deviation
C     ERRMX = largest error
C     (MAXI, MAXJ) = overlap integral with error ERRMAX

      ERRAV = 0.0D0
      ERRAB = 0.0D0
      RMS   = 0.0D0
      ERRMX = 0.0D0
      MAXI  = 0
      MAXJ  = 0
C
C     Analysis of interatomic (s|s), (p|p), (d|d) and (f|f) overlap errors
C
      SSOVER = .FALSE.
      PPOVER = .FALSE.
      DDOVER = .FALSE.
      FFOVER = .FALSE.

      MAXS  = 0.0D0
      MAXIS = 0
      MAXJS = 0
      MAXP  = 0.0D0
      MAXIP = 0
      MAXJP = 0
      MAXD  = 0.0D0
      MAXID = 0
      MAXJD = 0
      MAXF  = 0.0D0
      MAXIF = 0
      MAXJF = 0
C
C     Relative values (only calculated for diagonal elements)
C
      RERRAV = 0.0D0
      RERRAB = 0.0D0
      RERRMX = 0.0D0
      MAXIR  = 0
      MAXJR  = 0
C
      SMATDF = 0.0D0
      SMATAN = 0.0D0
      SMATFR = 0.0D0
C
C     Loop over SO's
C
      DO IC = 1,MC
         DO IRP = 0,MAXREP
            DO I = 1,NBBAS(IRP,IC)
               IS = IBBAS(IRP,IC) + I
               CENT1 = ICLAB(IS)
               ORB1 = PLABEL(IPLAB(IS,2),2)(9:9)
               DO J = 1,NBBAS(IRP,IC)
                  JS = IBBAS(IRP,IC) + J
                  CENT2 = ICLAB(JS)
                  ORB2 = PLABEL(IPLAB(JS,2),2)(9:9)
C                  WRITE(6,*)'CENTER(',IS,',',JS,') = ',
C     &                 '(',CENT1,',',CENT2,')'
                  SMATDF = SMAT(IS,JS,1)
                  SMATAN = SMAT(IS,JS,2)
C
                  ERRAV = ERRAV + SMATDF
                  ERRAB = ERRAB + ABS(SMATDF)
                  RMS = RMS + SMATDF*SMATDF
                  
                  IF(CENT1.NE.CENT2) THEN 
C     Find largest internuclear (s|s) error
                     IF(ORB1.EQ.'s'.AND.ORB2.EQ.'s') THEN
                        SSOVER = .TRUE.
                        IF(ABS(SMATDF).GT.ABS(MAXS)) THEN
                           MAXS = SMATDF
                           MAXIS = IS
                           MAXJS = JS
                        ENDIF
C     Find largest internuclear (p|p) error
                     ELSEIF(ORB1.EQ.'p'.AND.ORB2.EQ.'p') THEN
                        PPOVER = .TRUE.
                        IF(ABS(SMATDF).GT.ABS(MAXP)) THEN
                           MAXP = SMATDF
                           MAXIP = IS
                           MAXJP = JS
                        ENDIF
C     Find largest internuclear (d|d) error
                     ELSEIF(ORB1.EQ.'d'.AND.ORB2.EQ.'d') THEN
                        DDOVER = .TRUE.
                        IF(ABS(SMATDF).GT.ABS(MAXD)) THEN
                           MAXD = SMATDF
                           MAXID = IS
                           MAXJD = JS
                        ENDIF
C     Find largest internuclear (f|f) error
                     ELSEIF(ORB1.EQ.'f'.AND.ORB2.EQ.'f') THEN
                        FFOVER = .TRUE.
                        IF(ABS(SMATDF).GT.ABS(MAXF)) THEN
                           MAXF = SMATDF
                           MAXIF = IS
                           MAXJF = JS
                        ENDIF
                     ENDIF 
                  ENDIF
C     Find largest error
                  IF(ABS(SMATDF) .GT. ABS(ERRMX)) THEN
                     ERRMX = SMATDF
                     MAXI = IS
                     MAXJ = JS
                  END IF
C   
                  IF(JS .EQ. IS) THEN
                     SMATFR = SMATDF/SMATAN
                     RERRAV = RERRAV + SMATFR
                     RERRAB = RERRAB + ABS(SMATFR)
                     IF (ABS(SMATFR) .GT. ABS(RERRMX)) THEN
                        RERRMX = SMATFR
                        MAXIR = IS
                        MAXJR = JS
                     END IF   
                  END IF
C   
               END DO
            END DO
         END DO
      END DO
C
C     Get labels of maximum error overlap SO's
C
      MXLB1 = PLABEL(IPLAB(MAXI,2),2)
      MXLB2 = PLABEL(IPLAB(MAXJ,2),2) 
      MXLB1S= PLABEL(IPLAB(MAXIS,2),2)
      MXLB2S= PLABEL(IPLAB(MAXJS,2),2)
      MXLB1P= PLABEL(IPLAB(MAXIP,2),2)
      MXLB2P= PLABEL(IPLAB(MAXJP,2),2)
      MXLB1D= PLABEL(IPLAB(MAXID,2),2)
      MXLB2D= PLABEL(IPLAB(MAXJD,2),2)
      MXLB1F= PLABEL(IPLAB(MAXIF,2),2)
      MXLB2F= PLABEL(IPLAB(MAXJF,2),2) 
      RMXLB1 = PLABEL(IPLAB(MAXIR,2),2)
      RMXLB2 = PLABEL(IPLAB(MAXJR,2),2)
C
C     Calculate average errors:
C
      NLSOS = 0
      NSSOS = 0

      DO I = 1, MAXREP+1
         NLSOS = NLSOS + NCOS(I,1)*NCOS(I,1)
         NSSOS = NSSOS + NCOS(I,2)*NCOS(I,2) 
      END DO
      
      IF(LEVYLE) THEN
         ERRAV = ERRAV/DBLE(NLSOS)
         ERRAB = ERRAB/DBLE(NLSOS)
         RMS = DSQRT(RMS/DBLE(NLSOS))
         RERRAV = RERRAV/DBLE(NTBAS(1))
         RERRAB = RERRAB/DBLE(NTBAS(1))
      ELSE
         ERRAV = ERRAV/(DBLE(NLSOS) + DBLE(NSSOS))
         ERRAB = ERRAB/(DBLE(NLSOS) + DBLE(NSSOS))
         RMS = DSQRT(RMS/(DBLE(NLSOS) + DBLE(NSSOS)))
         RERRAV = RERRAV/DBLE(NTBAS(0))
         RERRAB = RERRAB/DBLE(NTBAS(0))
      END IF
C
C     Print section
C
      CALL HEADER('Output from GRTEST', -1)
C
      WRITE(LUPRI,'(/A)') '* Errors: '
C
      WRITE (LUPRI,'( 2X,A,E16.10)')
     &      ' Average SMAT error               : ', ERRAV
      WRITE (LUPRI,'( 2X,A,E16.10)')
     &      ' Absolute average SMAT error      : ', ERRAB      
      WRITE (LUPRI,'( 2X,A,E16.10)')
     &      ' RMS deviation                    : ', RMS
      WRITE (LUPRI,'(/A,E16.10)')
     &    '   Largest error                    : ', ERRMX
      WRITE(LUPRI,'(2X,A,I4,A,I4,A,A12,A,A12,A)')
     &      ' Position and label               : (',
     &      MAXI,',',MAXJ,') = (',MXLB1,',',MXLB2,')'
      IF(SSOVER) THEN
         WRITE (LUPRI,'(/A,E16.10)')
     &      '   Largest internuclear (s|s) error : ', MAXS
         WRITE(LUPRI,'(2X,A,I4,A,I4,A,A12,A,A12,A)')
     &        ' Position and label               : (',
     &        MAXIS,',',MAXJS,') = (',MXLB1S,',',MXLB2S,')'
      ENDIF
      IF(PPOVER) THEN
         WRITE (LUPRI,'(/A,E16.10)')
     &        '   Largest internuclear (p|p) error : ', MAXP
         WRITE(LUPRI,'(2X,A,I4,A,I4,A,A12,A,A12,A)')
     &        ' Position and label               : (',
     &        MAXIP,',',MAXJP,') = (',MXLB1P,',',MXLB2P,')'
      ENDIF
      IF(DDOVER) THEN
         WRITE (LUPRI,'(/A,E16.10)')
     &        '   Largest internuclear (d|d) error : ', MAXD
         WRITE(LUPRI,'(2X,A,I4,A,I4,A,A12,A,A12,A)')
     &        ' Position and label               : (',
     &        MAXID,',',MAXJD,') = (',MXLB1D,',',MXLB2D,')'
      ENDIF
      IF(FFOVER) THEN
         WRITE (LUPRI,'(/A,E16.10)')
     &        '   Largest internuclear (f|f) error : ', MAXF
         WRITE(LUPRI,'(2X,A,I4,A,I4,A,A12,A,A12,A)')
     &        ' Position and label               : (',
     &        MAXIF,',',MAXJF,') = (',MXLB1F,',',MXLB2F,')'
      ENDIF

C
      WRITE(LUPRI,'(/A)')'* Relative errors (diagonal elements only) : '
C
      WRITE (LUPRI,'( 2X,A,E16.10)')
     &      ' Average SMAT error               : ', RERRAV
      WRITE (LUPRI,'( 2X,A,E16.10)')
     &      ' Absolute average SMAT error      : ', RERRAB      
      WRITE (LUPRI,'(/A,E16.10)')
     &      '   Largest error                    : ', ERRMX
      WRITE(LUPRI,'(2X,A,I4,A,I4,A,A12,A,A12,A)')
     &      ' Position and label               : (',
     &      MAXIR,',',MAXJR,') = (',RMXLB1,',',RMXLB2,')'

      IF (dft_cfg_integration_check_level .GE. 2) THEN
         CALL HEADER('Errors in numerical SMAT ', -1)
         CALL OUTPUT(SMAT(1,1,1),1,NTBAS(0),1,NTBAS(0),NTBAS(0),
     &        NTBAS(0),1,LUPRI)
      END IF

      RETURN 
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck grtes2 */                                                C
      SUBROUTINE GRTES2(CORX,CORY,CORZ,WEIGHT,NPOINT,
     &               SMAT,
     &               GAO,GAO1,GAO2,GAB1,
     &               NCNT,BUF)
C                                                                      C
C     Written by Trond Saue and T.Helgaker                             C
C     Modified for grid testing by Olav Fossgaard July 2002.           C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      use dft_cfg
      use dirac_ao_eval

#include "implicit.h"
#include "priunit.h"

      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, D2 = 2.0D0, D3 = 3.0D0,
     &           D4 = 4.0D0, DP3 = D1/D3, DP5 = 0.5D0)

#include "mxcent.h"
#include "nuclei.h"
#include "inforb.h"
#include "orgcom.h"
#include "functionals.h"
      
      DIMENSION SMAT(*)
      DIMENSION CORX(NPOINT),CORY(NPOINT),CORZ(NPOINT),
     &            WEIGHT(NPOINT),
     &            GAO(NBAST), GAO1(NBAST,3), GAO2(NBAST,6),
     &            GAB1(NBAST,3), 
     &            NCNT(NBAST), BUF(NBAST,4)
      DIMENSION RHG(3)
      nder = 0
      if (need_ao_order(1, 0)) nder = 1

      DO 300 IPNT = 1, NPOINT
         WGHT = WEIGHT(IPNT)
C
C       AOs
C       ===
C
        call getsos(gao,gao1,dummy,dummy,
     &               ncnt,
     &              CORX(IPNT),CORY(IPNT),CORZ(IPNT),
     &               BUF,
     &               nbast,nder,.false.,0)
C
C       Test accuracy of numerical integration
C       ======================================
C
        CALL GRTES4(WGHT,0,SMAT,GAO)

  300 CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck grtes3 */
      SUBROUTINE GRTES3(X,Y,Z,WT,LUQUAD,NPOINT)
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 
#include "priunit.h"
      DIMENSION X(NPOINT),Y(NPOINT),Z(NPOINT),WT(NPOINT)
      READ(LUQUAD) X,Y,Z,WT
      
      RETURN
      END


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C     /* Deck grtes4 */
      SUBROUTINE GRTES4(VXC,IREPOP,EXCMAT,GAO)      
C
C     Construct Kohn-Sham matrix                                       C
C                                                                      C
C     Written by T.Saue and T.Helgaker                                 C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h"
#include "priunit.h"

      PARAMETER(D0=0.0D0,D2=2.0D0)

#include "dcbham.h"
#include "dgroup.h"
#include "dcbbas.h"

#if defined (VAR_G77) || defined (SYS_SX)
      EQUIVALENCE (NBAST,NTBAS(0))

      DIMENSION GAO(NBAST),EXCMAT(NBAST,NBAST)
#else
      DIMENSION GAO(NTBAS(0)),EXCMAT(NTBAS(0),NTBAS(0))
#endif
      
#include "ibtfun.h"

      DO 10 IRP1 = 0,NBSYM-1
        DO 20 IRP2 = 0,NBSYM-1
          IF(IBTXOR(IRP1,IRP2).NE.IREPOP) GOTO 20
          DO 30 IC = 1,MC
            NB1 = NBBAS(IRP1,IC)
            IF(NB1.EQ.0) GOTO 30
             NB2 = NBBAS(IRP2,IC)
            IF(NB2.EQ.0) GOTO 30
            JB1 = IBBAS(IRP1,IC)
            JB2 = IBBAS(IRP2,IC)
            DO J = 1, NBBAS(IRP2,IC)
              GVXC = VXC*GAO(JB2+J)
              CALL DAXPY(NBBAS(IRP1,IC),GVXC,
     &                   GAO(JB1+1),1,EXCMAT(JB1+1,JB2+J),1)            
            ENDDO
   30     CONTINUE
   20   CONTINUE
   10 CONTINUE

      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck getwt1 */
      SUBROUTINE GETWT1(X,Y,Z,WT,IPRINT)
C                                                                      C
C     Interface routine for grid points                                C
C     T. Helgaker                                                      C
C                                                                      C
C     Modified for GETW1 by O. Fossgaard - May 2003                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"

      PARAMETER (MAXCEN=100)
      PARAMETER (D0=0.0D0)
      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)

#include "dftwrk.h"

      DIMENSION WT(NATOMS)

#include "ibtfun.h"

      real(8), allocatable :: psmu(:)
      real(8), allocatable :: rj(:)
      real(8), allocatable :: rij(:)
      real(8), allocatable :: aij(:)

C
C     ********************************
C     ***** Set up common blocks *****
C     ********************************
C
C     INFOA
C     =====
CTROND: Modified to handle symmetry dependent atoms
C       as well
      NAT = 0
      NUM = NBASIS
      DO ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         ZCNT   = CHARGE(ICENT)
         IF(NAMN(ICENT)(1:2).NE.'Gh') THEN
           DO ISYMOP = 0, MAXOPR         
           IF (IBTAND(ISYMOP,MULCNT) .EQ. 0) THEN
             NAT = NAT + 1
             ZAN(NAT) = ZCNT
             C(1,NAT) = PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
             C(2,NAT) = PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
             C(3,NAT) = PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
           ENDIF
           ENDDO
         ENDIF
      ENDDO
C
C     *******************************************
C     ***** Prepare arguments for CONSTRUCT *****
C     *******************************************
C
      NNAT  = NAT*(NAT-1)/2
      N2AT = NAT*NAT 

      call alloc(psmu, nat)
      call alloc(rj,   nat)
      call alloc(rij,  nnat)
      call alloc(aij,  nnat)
C
C     Calculate weights in point (X,Y,Z)
C     ==================================
C
      CALL WT1CON(X,Y,Z,WT,
     &            psmu,rj,rij,aij,iprint)

      call dealloc(psmu)
      call dealloc(rj)
      call dealloc(rij)
      call dealloc(aij)
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck wt1con */
      SUBROUTINE WT1CON(X,Y,Z,WT,PSMU,RJ,
     &                    RIJ,AIJ,IPRINT)
C                                                                      C
C      Creates a quadrature grid and writes it to disk                 C
C      Original version by C.W.Murray. Rewritten by                    C
C      A.M.Lee and D.J.Tozer.                                          C
C                                                                      C
C      Adapted for Dalton by T. Helgaker                               C
C      Adapted for DIRAC by T.Saue                                     C
C                                                                      C
C      Modified for WT1CON by C. Fossgaard - May 2003                  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 
#include "priunit.h" 

      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0,D3=3.0D0)

      PARAMETER (MAXCEN=100)

#include "dftwrk.h"

      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)

#include "dummy.h"
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "nuclei.h"
#include "symmet.h"
#include "ccom.h"

      DIMENSION PSMU(NAT),RJ(NAT),WT(NATOMS)
C
C     CARE: Some equivalencing through argument list!
C
      DIMENSION RIJ(NAT*(NAT-1)/2),AIJ(NAT*(NAT-1)/2)
C
C     Bragg contains Bragg-Slater radii for atoms. 
C     Following Becke, the hydrogen radius is 0.35. 
C     For the noble gases, the values have been guessed.
C
      DIMENSION BRAGG(0:103)
      DIMENSION XPASC(20)

C Trond Saue:
C     The below data gives atomic radii in Angstroms and stems from table I of 
C     J.C.Slater: "Atomic Radii in Crystals"
C     J.Chem.Phys. 41(1964) 3199-3204
C     Values for elements marked with an asterisk has been
C     guessed/interpolated
C
      DATA BRAGG/ 0.75D0,
C       H      He*   
     &  0.35D0,0.35D0,
C       Li     Be     B      C      N      O      F      Ne*
     &  1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.50D0,0.45D0,
C     &  1.45D0,1.05D0,0.85D0,0.70D0,0.65D0,0.60D0,0.886403D0,0.45D0,
C       Na     Mg     Al     Si     P      S      Cl     Ar*
     &  1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,1.00D0,1.00D0,
C     &  1.80D0,1.50D0,1.25D0,1.10D0,1.00D0,1.00D0,0.786807D0,1.00D0,
C       K      Ca     Sc     Ti     V      Cr     Mn     Fe     Co 
     &  2.20D0,1.80D0,1.60D0,1.40D0,1.35D0,1.40D0,1.40D0,1.40D0,1.35D0,
C       Ni     Cu     Zn     Ga     Ge     As     Se     Br     Kr*
     &  1.35D0,1.35D0,1.35D0,1.30D0,1.25D0,1.15D0,1.15D0,1.15D0,1.10D0,
C       Rb     Sr     Y      Zr     Nb     Mo     Tc     Ru     Rh
     &  2.35D0,2.00D0,1.80D0,1.55D0,1.45D0,1.45D0,1.35D0,1.30D0,1.35D0,
C       Pd     Ag     Cd     In     Sn     Sb     Te     I      Xe*
     &  1.40D0,1.60D0,1.55D0,1.55D0,1.45D0,1.45D0,1.40D0,1.40D0,1.40D0,
C       Cs     Ba     La    
     &  2.60D0,2.15D0,1.95D0,
C       Ce     Pr     Nd     Pm     Sm     Eu     Gd
     &  1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.85D0,1.80D0,
C       Tb     Dy     Ho     Er     Tm     Yb     Lu
     &  1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
C       Hf     Ta     W      Re     Os     Ir     Pt     Au     Hg
     &  1.55D0,1.45D0,1.35D0,1.30D0,1.30D0,1.35D0,1.35D0,1.35D0,1.50D0,
C       Tl     Pb*    Bi     Po     At*    Rn*
     &  1.90D0,1.75D0,1.60D0,1.90D0,1.50D0,1.50D0,
C       Fr*    Ra     Ac     
     &  2.15D0,2.15D0,1.95D0,
CTROND rad(U): 1.75 --> 1.37D0
C       Th     Pa     U      Np     Pu     Am     Cm*     
     &  1.80D0,1.80D0,1.37D0,1.75D0,1.75D0,1.75D0,1.75D0,
CTROND       Th     Pa     U      Np     Pu     Am     Cm*     
CTROND     &  1.80D0,1.80D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,
C       Bk*    Cf*    Es*    Fm*    Md*    No*    Lw*
     &  1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0,1.75D0/ 

#include "ibtfun.h"

C
C     USE INPUTTED INFORMATION IN COMMON BLOCK DFTWRK TO DEFINE THE TYPE
C     OF QUADRATURE
C
C     Experimental value of NTRANS
C
      NTRANS=10
C
C     Form Pascals triangle in XPASC for fuzzy Voronoi polyhedra code
C
      ISIGN=-1
      DO 5 I=NTRANS,1,-1
         ISIGN=-ISIGN
         XPASC(I+1) = ISIGN*FACULT(NTRANS)/(FACULT(I)*FACULT(NTRANS-I))
         XPASC(1)   = 1.0D0
5     CONTINUE

      APASC=0.0D0
      DO 6 I=1,NTRANS+1
         XPASC(I)=XPASC(I)/DBLE(2*I-1)
         APASC=APASC+XPASC(I)
6     CONTINUE
      APASC=0.5D0/APASC

      ITEMP=0
      DO 7 INA=1,NAT
      DO 7 JNA=1,INA-1
          ITEMP=ITEMP+1
          RIJ(ITEMP) = 1.0D0/(DSQRT((C(1,INA)-C(1,JNA))**2
     &                             +(C(2,INA)-C(2,JNA))**2
     &                             +(C(3,INA)-C(3,JNA))**2))
          CHI=BRAGG(NINT(ZAN(INA)))/BRAGG(NINT(ZAN(JNA)))
          TEMP=(CHI-1.0D0)/(CHI+1.0D0)
          AIJ(ITEMP)=TEMP/(TEMP*TEMP-1.0D0)
          IF (AIJ(ITEMP).GT.0.5D0) AIJ(ITEMP)=0.5D0
          IF (AIJ(ITEMP).LT.-0.5D0) AIJ(ITEMP)=-0.5D0
7     CONTINUE
C
C     Loop over atoms
C     ===============
C 
      MXOPR = 0 
      NATOM = 1
      WEIGHT = 0.0D0
CTROND: define scaling factor for Bragg: set to one !!!!
      TFAC = 1.0D0

C      DO 10 IATOM = 1, NUCIND
      DO 10 IATOM = 1,NATOMS
C     DO 10 ISYMOP = 0, MXOPR         
C     IF (IBTAND(ISYMOP,ISTBNU(IATOM)).EQ.0) THEN 
C        NDEG  = MULT(ISTBNU(IATOM))
         IF(NAMN(IATOM)(1:2).EQ.'Gh') GOTO 10
C         MULA  = ISTBNU(IATOM)
C         NDEG  = MULT(MULA)
C
C        Becke partitioning
C        ==================
C
         CALL WGHT1(RJ,PSMU,RIJ,AIJ,APASC,XPASC,
     &              X,Y,Z,WEIGHT,
     &              NTRANS,NATOM)
C         IF(NDEG.GT.0) THEN
C          Scale weights with multiplicity of center
C           FAC = dble(NDEG) 
           WT(NATOM) = WEIGHT
C           WRITE(6,*) 'WT(',NATOM,') = ', WT(NATOM)
C         ENDIF
      
C      NATOM = NATOM + NDEG
      NATOM = NATOM + 1

C      END IF

10    CONTINUE

      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  /* Deck wght1 */                                                   C
      SUBROUTINE WGHT1(RJ,PSMU,RIJ,AIJ,APASC,XPASC,X,Y,Z,WEIGHT,
     &                 NTRANS,NATOM)
C                                                                      C
C      Written by C. W. Murray                                         C
C      WGHT1 calculates the weights associated with                    C
C      the Becke partitioning amongst the atoms.                       C
C                                                                      C         
C      Modified for WGHT1 by O. Fossgaard - May 2003                   C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#include "implicit.h" 

      PARAMETER (MAXCEN=100)
      COMMON /INFOA/ NAT,NUM,ZAN(MAXCEN),C(3,MAXCEN)
      DIMENSION RJ(NAT),
     &          PSMU(NAT),
     &          RIJ(NAT*(NAT-1)/2),
     &          AIJ(NAT*(NAT-1)/2),
     &          XPASC(20)

C
C     WEIGHT FUNCTION BY BECKE PARTITIONING
C
      DO 100 INA=1,NAT
         RJ(INA)=DSQRT((C(1,INA)-X)**2
     &                  +(C(2,INA)-Y)**2
     &                  +(C(3,INA)-Z)**2)
         PSMU(INA)=1.0D0
  100 CONTINUE

      ITEMP=0
      DO 200 INA=1,NAT
      DO 200 JNA=1,INA-1
         ITEMP=ITEMP+1
            XMU=(RJ(INA)-RJ(JNA))*RIJ(ITEMP)
            XMUIJ=XMU+AIJ(ITEMP)*(1-XMU*XMU)
            XMUIJ2=XMUIJ*XMUIJ
            XMUIJN=XMUIJ
            ACCUM=0.0D0
  210    CONTINUE
         DO 220 I=1,NTRANS+1
            ACCUM=ACCUM+XPASC(I)*XMUIJN
            XMUIJN=XMUIJN*XMUIJ2
  220    CONTINUE
            PSMU(INA)=PSMU(INA)*(0.5D0-APASC*ACCUM)
            PSMU(JNA)=PSMU(JNA)*(0.5D0+APASC*ACCUM)
  230    CONTINUE
  200 CONTINUE

      ACCUM = 0.0D0
      DO 300 INA=1,NAT
         ACCUM=ACCUM+PSMU(INA)
  300 CONTINUE
C
C     Contract ACCUM into PSMU
C
      DO 400 INA=1,NAT
         PSMU(INA)=PSMU(INA)/ACCUM
C         WRITE(6,*) 'RJ(', INA, ') = ', RJ(INA)
C         WRITE(6,*) 'PSMU(',INA,') = ', PSMU(INA)
  400 CONTINUE

      WEIGHT=PSMU(NATOM)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1wt  */
      SUBROUTINE GT1WT()
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C     Modified for GT1WT by O. Fossgaard - May 2003
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
C
#include "dcbgen.h"
#include "dcbwt1.h"
#include "nuclei.h"
      real(8), allocatable :: kwt(:)
      real(8), allocatable :: kcoor(:)
C
      CALL QENTER('GT1WT')
C
      CALL TITLER('Integration weights along bonds','*',122)
C
      call alloc(kwt, natoms) 
      call alloc(kcoor, 6) 

      CALL GT1BW1(KWT, KCOOR, 0)
      call dealloc(kwt)
      call dealloc(kcoor)
C
 9999 CALL QEXIT('GT1WT')
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck gt1bw1 */
      SUBROUTINE GT1BW1(WT,COORD,IPRINT)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
C
#include "symmet.h"
#include "nuclei.h"
      DIMENSION COORD(3,2),WT(NATOMS)
      CHARACTER FILNAM*10
#include "ibtfun.h"
C
      LUWT = 1
C
C     Loop over pairs of nuclear centers
C
      DO ICENTA = 1, NUCIND
         MULA  = ISTBNU(ICENTA)
         IDEGA = 1
         COORD(1,1) = CORD(1,ICENTA)
         COORD(2,1) = CORD(2,ICENTA)
         COORD(3,1) = CORD(3,ICENTA)
         DO ICENTB = ICENTA+1, NUCIND
           MULB = ISTBNU(ICENTB)
           IDEGB = 0 
           DO ISYMOP = 0, MAXOPR
             IF (IBTAND(ISYMOP,MULB) .EQ. 0) THEN
               IDEGB = IDEGB + 1
               COORD(1,2)=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENTB)
               COORD(2,2)=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENTB)
               COORD(3,2)=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENTB)
               WRITE(FILNAM,'(A2,2(A3,I1))') 
     &           'WT',NAMN(ICENTA),IDEGA,NAMN(ICENTB),IDEGB
               DO I = 1,10
                 IF(FILNAM(I:I).EQ.' ') FILNAM(I:I)='_'
               ENDDO
               OPEN(LUWT,FILE =FILNAM,STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
               CALL GT1W2(ICENTA,ICENTB,WT,COORD,LUWT,IPRINT)
               CLOSE(LUWT,STATUS='KEEP')
             ENDIF
           ENDDO
         ENDDO
         IF (MULT(MULA).GT.1) THEN
           IDEGB = 1
           DO ISYMOP = 1, MAXOPR
             IF (IBTAND(ISYMOP,MULA) .EQ. 0) THEN
               IDEGB = IDEGB + 1
               COORD(1,2)=PT(IBTAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENTA)
               COORD(2,2)=PT(IBTAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENTA)
               COORD(3,2)=PT(IBTAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENTA)
               WRITE(FILNAM,'(A2,2(A3,I1))') 
     &           'WT',NAMN(ICENTA),IDEGA,NAMN(ICENTA),IDEGB
               DO I = 1,10
                 IF(FILNAM(I:I).EQ.' ') FILNAM(I:I)='_'
               ENDDO
               OPEN(LUWT,FILE =FILNAM,STATUS='NEW',
     &              ACCESS='SEQUENTIAL',FORM = 'FORMATTED')
               CALL GT1W2(ICENTA,ICENTB,WT,COORD,LUWT,IPRINT)
               CLOSE(LUWT,STATUS='KEEP')
              ENDIF
           ENDDO
         ENDIF 
      ENDDO
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck Gt1w2 */
      SUBROUTINE GT1W2(ICENTA,ICENTB,WT,COORD,LUWT,IPRINT)
C***********************************************************************
C
C     Get density in one dimension, that is along bonds
C     (ghost centers may be defined to get other lines)
C
C     Written by Trond Saue May 9 2000
C
C***********************************************************************
#include "implicit.h"
#include "codata.h"
#include "priunit.h"
#include "consts.h"
#include "mxcent.h"
#include "nuclei.h"
C
      DIMENSION COORD(3,2),WT(NATOMS)
#include "dcbwt1.h"
C
C     Calculate distance between centers
C
      DX = (COORD(1,2)-COORD(1,1))
      DY = (COORD(2,2)-COORD(2,1))
      DZ = (COORD(3,2)-COORD(3,1))
      DISTAN = SQRT(DX*DX+DY*DY+DZ*DZ)
      DX     = DX/DISTAN
      DY     = DY/DISTAN
      DZ     = DZ/DISTAN
      ASTEP  = DSTEP/XTANG
      NPOINT = INT(DISTAN/ASTEP)
      IF(MOD(NPOINT,2).EQ.1) NPOINT = NPOINT+1
      ASTEP  = DISTAN/dble(NPOINT)
C
C     Calculate density along between centers:
C     If bond distance is 2R, then points are calculated between (-R,3R)
C 
      IP1 = -NPOINT/2
      IP2 = NPOINT - IP1 
      DO IP = IP1,IP2
        PP  = dble(IP)*ASTEP
        PX  = COORD(1,1) + PP*DX
        PY  = COORD(2,1) + PP*DY
        PZ  = COORD(3,1) + PP*DZ

C     Get weights for each atom at (PX,PY,PZ)

        wt = 0.0d0

        CALL GETWT1(PX,PY,PZ,WT,IPRINT)
        PP = XTANG*PP
C        WRITE(LUWT,'(F10.6,3E20.10)') PP,(WT(I),I=1,NATOMS)
        WRITE(LUWT,'(F10.6,3E20.10)') PP,WT(ICENTA),WT(ICENTB)
      ENDDO
C
      END
