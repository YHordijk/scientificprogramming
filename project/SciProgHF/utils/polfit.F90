!      Copyright (c) 2018 by the authors of DIRAC.
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

      PROGRAM POLFIT
#include "implicit.h"
#include "priunit.h"
      PARAMETER(MORDER=16,MP=50)
      PARAMETER(MDIM=MORDER+1)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0)
      PARAMETER(NITER=10)
      DIMENSION A(MP,MDIM),B(MP),C(MDIM),D(MDIM)
      DIMENSION X(MP),Y(MP)
      CHARACTER RESFIL*40,REPLY*1,POTFIL*44
      IUNIT = 1
      ILOGG = 2
      WRITE(LUPRI,'(A)') 'Name input file with potential curve(A16)'
      READ(LUSTDIN,'(A)') POTFIL
      i=LNBLNK(POTFIL)
      RESFIL = POTFIL(1:i)//'.polfit'
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',             &
      ACCESS='SEQUENTIAL',ERR=5)
      OPEN(ILOGG,FILE=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',         &
      ACCESS='SEQUENTIAL')
      GOTO 6
 5    CONTINUE
      STOP 'Error opening POTFIL'
 6    CONTINUE
      NP = 0
   10 CONTINUE  
! Read points
      NP = NP + 1
      IF(NP.GT.MP) STOP 'Too many points !!!'
      READ(IUNIT,*,END=20) X(NP),Y(NP)
      GOTO 10
 20   CONTINUE
      NP = NP-1
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(2(4X,E20.10))') (X(I),Y(I),I=1,NP)
 30   CONTINUE
      WRITE(*,'(A)') 'Give order of polynomial:'
      READ(LUSTDIN,*) NORDER
      WRITE(*,'(I3)') NORDER
      IF(NP.LE.NORDER) THEN
        WRITE(*,'(A,I3)')                                               &
     'Too few points for a polynomial fit of order',NORDER
        WRITE(*,'(A)') '  Try again.'
        GOTO 30
      ENDIF
      NDIM = NORDER + 1
      CALL POLSVD(NDIM,NP,A,B,X,Y,C,D,CHISQ,ISKIP)
!      CALL CURVE(NORDER,NP,A,B,X,Y,VAR,IBUF)
!* Estimate fit:
      WRITE(LUPRI,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(3X,A)') '* Coefficients:'
      WRITE(LUPRI,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',           &
       C(I),I=1,NDIM)  
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(LUPRI,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
!
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(3X,A)') '* Coefficients:'
      WRITE(ILOGG,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',           &
       C(I),I=1,NDIM)  
      WRITE(LUPRI,'(8X,A,4X,2A)')                                       &
    'X',' Predicted Y',' Relative error'
      DO I = 1,NP
        YP = POLVAL(NORDER,C,X(I))
        DEV = (YP-Y(I))/Y(I)
        WRITE(LUPRI,'(3X,F6.3,4X,2E12.4)') X(I),YP,DEV
      ENDDO
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
      WRITE(LUPRI,'(A)') '* Do you want to find min. of curve (y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 664
!* Find smallest y-value and take corresponding x-value as start value
!* for Newton-Raphsons method      
!      XVAL = X(MIND(NP,Y))
      IYMIN = IDMIN(NP,Y,1)
      XMIN = X(IYMIN)
      CALL NEWRAP(NORDER,DTOL,C,XMIN,NITER,IERR)
!* Give coordinates of local minimum
      IF(IERR.EQ.1) THEN
        WRITE(*,'(A)') 'Newtons method failed'
      ELSE
        YMIN = POLVAL(NORDER,C,XMIN)
        WRITE(LUPRI,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',     &
         XMIN,',',YMIN,')'
        WRITE(ILOGG,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',     &
         XMIN,',',YMIN,')'
      ENDIF
 664  CONTINUE
      WRITE(LUPRI,'(A)')                                                &
   '* Do you want a function value at some point(y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 665
      WRITE(LUPRI,'(A)')                                                &
   '* Give point:'
      READ(LUSTDIN,*) XMIN   
      YMIN = POLVAL(NORDER,C,XMIN)      
      WRITE(LUPRI,'(3X,F6.3,4X,E22.16)') XMIN,YMIN
      WRITE(ILOGG,'(3X,A,F6.3,4X,E22.16)')                              &
   'Selected value: ',XMIN,YMIN
      GOTO 664
!* Calculate derivatives:
 665  CONTINUE
      WRITE(LUPRI,'(A)')                                                &
   '* Do you want derivatives of curve at minimum(y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
      TMP = POLNDER(NORDER,C,XMIN,2)
      WRITE(LUPRI,'(A,1P,E18.10)') '2nd derivative at minimum: ',TMP
      WRITE(ILOGG,'(A,1P,E18.10)') '2nd derivative at minimum: ',TMP
      TMP = POLNDER(NORDER,C,XMIN,3)
      WRITE(LUPRI,'(A,1P,E18.10)') '3rd derivative at minimum: ',TMP
      WRITE(ILOGG,'(A,1P,E18.10)') '3rd derivative at minimum: ',TMP
      TMP = POLNDER(NORDER,C,XMIN,4)
      WRITE(LUPRI,'(A,1P,E18.10)') '4th derivative at minimum: ',TMP
      WRITE(ILOGG,'(A,1P,E18.10)') '4th derivative at minimum: ',TMP
 666  CONTINUE
      WRITE(LUPRI,'(A)')                                                &
   '* Do you want derivatives of curve at some point(y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 667
      WRITE(LUPRI,'(A)') 'Give maximal order of derivative..'
      READ(LUSTDIN,*) MAXDER
      MAXDER = MIN(MAXDER,NORDER)
      WRITE(LUPRI,'(A)')                                                &
   '* Give point:'
      READ(LUSTDIN,*) XMIN   
      WRITE(LUPRI,'(A,1P,E18.10)') 'X = ',XMIN
      WRITE(ILOGG,'(A,1P,E18.10)') 'X = ',XMIN
      DO I = 1,MAXDER
        TMP = POLNDER(NORDER,C,XMIN,I)
        WRITE(LUPRI,'(I2,A,1P,E18.10)') I,'. derivative : ',TMP
        WRITE(ILOGG,'(I2,A,1P,E18.10)') I,'. derivative : ',TMP
      ENDDO
      GOTO 666
 667  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want a new polynomial fit ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.EQ.'y') GOTO 30
      END
!
!
