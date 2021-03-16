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

      PROGRAM MX2FIT
!***********************************************************************
!
!     Utility program for triatomics MX2
!     This program will do a quartic fit of potential curve,
!     then calculate spectroscopic constants:
!       - equilibrium distance
!       - harmonic frequency
!       - anharmonic constant
!       - mean displacement in harmonic ground state
!
!     Written by T. Saue 2005
!
!***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(MDIM=20,MP=50,MNDIM=(MDIM*(MDIM+1)/2))
      PARAMETER(MITER=20,MORDER=MDIM+1)
      PARAMETER(M2DIM=MDIM*MDIM)
      DIMENSION A(MP,MDIM),B(MP),C(MDIM),D(MDIM)
      DIMENSION X(MP),Y(MP),YB(MP),XB(MP)
      CHARACTER POTFIL*6,RESFIL*6,REPLY*1
      COMMON/INFO/IUNIT,ILOGG
!
!     Heading
!     =======
!
      CALL TITLER                                                       &
   ('MX2FIT for triatomics MX2: Written by T. Saue ','*',110)
      IUNIT = 1
      ILOGG = 2
      WRITE(LUPRI,'(A)') '* Name output file (A6)'
      READ(LUSTDIN,'(A6)') RESFIL
      OPEN(ILOGG,FILE=RESFIL,STATUS='UNKNOWN',FORM='FORMATTED',         &
      ACCESS='SEQUENTIAL')
!
!
!     Get masses
!     ==========
!
 10   CONTINUE
      WRITE(LUPRI,'(A)')                                                &
   'Select one of the following:',                                      &
   '  1. Select masses of the most abundant isotopes.',                 &
   '  2. Employ user-defined atomic masses.'
      READ(LUSTDIN,*) I
      IF(I.EQ.1) THEN
        WRITE(LUPRI,'(A)') '* Give charge of central atom M:'
        READ(LUSTDIN,*) IZ
        CALL GETMSS(IZ,AM,ABUND,INFO)
        IF(INFO.EQ.0) THEN
          WRITE(LUPRI,'(A,F12.4)')                                      &
     '* Mass     :', AM,                                                &
     '* Abundance:', ABUND
        ELSE
          WRITE(LUPRI,'(A)') 'Routine failed. Give mass in Daltons:'
          READ(LUSTDIN,*) AM
        ENDIF
        WRITE(LUPRI,'(A)') '* Give charge of ligand atom X :'
        READ(LUSTDIN,*) IZ
        CALL GETMSS(IZ,BM,ABUND,INFO)
        IF(INFO.EQ.0) THEN
          WRITE(LUPRI,'(A,F12.4)')                                      &
     '* Mass     :', BM,                                                &
     '* Abundance:', ABUND
        ELSE
          WRITE(LUPRI,'(A)') 'Routine failed. Give mass in Daltons:'
          READ(LUSTDIN,*) BM
        ENDIF
      ELSEIF(I.EQ.2) THEN
        WRITE(LUPRI,'(A)') '* Give mass of central atom M in Daltons:'
        READ(LUSTDIN,*) AM
        WRITE(LUPRI,'(A)') '* Give mass of ligand atom X in Daltons:'
        READ(LUSTDIN,*) BM
      ELSE
        WRITE(LUPRI,'(A)') ' You stupid fart ! Try again !'
        GOTO 10
      ENDIF
      WRITE(ILOGG,'(A/)') '* MASSES:'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom M      :',AM
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom O      :',BM
!**********************************************************************
!     Select vibrational mode
!**********************************************************************
 20   CONTINUE
      WRITE(LUPRI,'(A)') '* Select vibrational mode to analyze:'
      WRITE(LUPRI,'(A)')                                                &
   '1. Symmetric stretch    <-- -->',                                   &
   '2. Asymmetric stretch   <--  <-',                                   &
   '3. Angular deformation         ',                                   &
   '4. Quit....'
      READ(LUSTDIN,*) I
      IF    (I.EQ.1) THEN
        CALL SYMSTR(A,B,C,D,X,Y,MORDER,MP,MITER,AM,BM,IBUF)
      ELSEIF(I.EQ.2) THEN
        CALL ASMSTR(A,B,C,D,X,Y,YB,XB,MORDER,MP,MITER,AM,BM,IBUF)
      ELSEIF(I.EQ.3) THEN
        CALL ANGDEF(A,B,C,D,X,Y,MORDER,MP,MITER,AM,BM,IBUF)
      ELSEIF(I.EQ.4) THEN      
        GOTO 30
      ELSE
        WRITE(LUPRI,'(A)') ' You stupid fart ! Try again !'
      ENDIF
      GOTO 20
 30   CONTINUE
!
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SYMSTR(A,B,C,D,X,Y,MORDER,MP,MITER,AM,BM,IBUF)
      use codata
#include "implicit.h"
#include "priunit.h"
!
      PARAMETER(ANGM = 1.0D-10,XTCM = 2.19474625D+05)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      PARAMETER(D0 = 0.0D0,DMC = 0.01D0)
      PARAMETER(D2PI = D2*PI)
      PARAMETER(NITER = 50)
      DIMENSION A(*),B(*),C(*),D(*),IBUF(*)
      DIMENSION X(MP),Y(MP)
      CHARACTER POTFIL*6,REPLY*1
      DOUBLE PRECISION XTM,AMUKG
!
      call set_codata_values(CODSET)
      XTM   = XTANG*ANGM
      AMUKG = UMASS
!
      IUNIT = 1
      ILOGG = 2
!     Reduced mass
      UM = BM
!
!     Read points
!     ===========
!
      WRITE(LUPRI,'(A,F12.4)')                                          &
    '* Symmetric stretch  <-- M -->. Reduced mass:',UM
      WRITE(LUPRI,'(A)') 'Symmetry coordinate: [1/SQRT(2)](r12+r32)'
 10   CONTINUE
      WRITE(LUPRI,'(A)') 'Name input file with potential curve(A6)'
      READ(LUSTDIN,'(A)') POTFIL
      WRITE(LUPRI,'(A)') 'Format: ',                                    &
         'lines of  r12(Angstrom)   energy(Hartree)'
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',             &
      ACCESS='SEQUENTIAL',ERR=20)
      GOTO 30
 20   CONTINUE
      STOP 'Error opening POTFIL'
 30   CONTINUE
      NP = 0
      DO I = 1,MP
        READ(IUNIT,*,END=40) X(I),Y(I)
        NP = NP + 1
      ENDDO
 40   CONTINUE
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
!
!     Polynomial fit
!     ==============
!
      WRITE(ILOGG,'(/A/)')                                              &
    '* SYMMETRIC STRETCH  <-- M -->'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Reduced mass:',UM
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(F7.4,4X,F20.10)') (X(I),Y(I),I=1,NP)
!
!     Generate symmetry coordinate and convert to atomic units
!
      UM = XFAMU*UM
      FAC = SQRT(D2)/XTANG
      CALL DSCAL(NP,FAC,X,1)
 60   CONTINUE
      WRITE(LUPRI,'(A)') 'Polynomial fit: Give order of polynomial'
      READ(LUSTDIN,*) NORDER
      IF(NP.LE.NORDER) THEN
        WRITE(LUPRI,'(A,I3)')                                           &
     'Too few points for a polynomial fit of order',NORDER
        GOTO 666
      ENDIF
      IF(NORDER.GT.MORDER) THEN
        WRITE(LUPRI,'(A)') 'Order to large ....'
        GOTO 666
      ENDIF
      NDIM = NORDER + 1
      CALL POLSVD(NDIM,NP,A,B,X,Y,C,D,CHISQ,ISKIP)
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
!TROND      DO I = 1,NP
!TROND        YP = (POLVAL(NORDER,C,X(I))-Y(I))/Y(I)
!TROND        WRITE(LUPRI,'(3X,F6.3,4X,E12.4)') XTANG*X(I),YP
!TROND      ENDDO
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
!
!     Find minimum
!     ============
!
      WRITE(LUPRI,'(A)') '* Do you want to find a minimum ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 100
!* Find smallest y-value and take corresponding x-value as start value
!* for Newton-Raphsons method   
      IYOPT = IDMIN(NP,Y,1)
      XOPT = X(IYOPT)
      CALL NEWRAP(NORDER,DTOL,C,XOPT,NITER,IERR)
      IF(IERR.EQ.1) THEN
        WRITE(ILOGG,'(A)') 'Newtons method failed'
        STOP 'Did not find minimum. Newtons method failed...'
      ELSE
        YOPT = POLVAL(NORDER,C,XOPT)
        REQ = XOPT/FAC
        WRITE(LUPRI,'(A,F18.5,A)')                                      &
     '* Local minimum   :',REQ,' Angstroms'
        WRITE(ILOGG,'(A,F18.5,A)')                                      &
     '* Local minimum   :',REQ,' Angstroms'
        WRITE(LUPRI,'(A,1P,E18.10,A)')                                  &
     '* Expected energy :',YOPT,' Hartrees'
        WRITE(ILOGG,'(A,1P,E18.10,A)')                                  &
     '* Expected energy :',YOPT,' Hartrees'
      ENDIF
!
!     Harmonic constant
!     ==================
!
      WRITE(LUPRI,'(A)') '* Do you want harmonic constant ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 100
      V2  = POL2DER(NORDER,C,XOPT)
      FORCE = V2*XTJ/XTM/XTM
      FREQ  = SQRT(V2/UM)
      WRITE(LUPRI,*) V2,FORCE,FREQ
      WRITE(LUPRI,'(A,2(F14.4,A))') '* Force constant:',FORCE,' N/m = ',&
     FORCE*DMC, ' mDyne/A'
      WRITE(LUPRI,'(A,E10.3,A,3X,F14.4,A)') '* Frequency:',FREQ*XTHZ,   &
     'Hz --> ', FREQ*XTCM,'cm-1'
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,F7.4)') '* Req = ',REQ
      WRITE(ILOGG,'(A,2(F14.4,A))') '* Force constant:',FORCE,' N/m = ',&
     FORCE*DMC, ' mDyne/A'
      WRITE(ILOGG,'(A,E10.3,A,3X,F14.4,A)') '* Frequency:',FREQ*XTHZ,   &
     'Hz --> ', FREQ*XTCM,'cm-1'
!
!     Derivatives of curve
!     ====================
!
 100  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want derivatives of curve (y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
      WRITE(LUPRI,'(A)') 'Derivatives to what order(max.6) ?'
      READ(LUSTDIN,*) NDER
      DO J = 1,NP
        WRITE(ILOGG,'(F6.3,3X,1P,6E14.6)')                              &
    X(J),(POLNDER(NORDER,B,X(J),I),I = 1,NDER)
      ENDDO
      
 666  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want a new polynomial fit ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.EQ.'y') GOTO 60
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ASMSTR(A,B,C,D,X,Y,YB,XB,MORDER,MP,MITER,              &
                   AM,BM,IBUF)
      use codata
#include "implicit.h"
#include "priunit.h"
!
      PARAMETER(ANGM = 1.0D-10,XTCM=2.19474625D+05)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      PARAMETER(D0 = 0.0D0,DMC = 0.01D0)
      PARAMETER(D2PI = D2*PI)
      PARAMETER(NITER = 50)
      DIMENSION A(*),B(*),C(*),D(*),IBUF(*)
      DIMENSION X(MP),Y(MP),YB(MP),XB(MP)
      CHARACTER POTFIL*6,REPLY*1
      DOUBLE PRECISION XTM,AMUKG
!
      call set_codata_values(CODSET)
      XTM   = XTANG*ANGM
      AMUKG = UMASS
!
      IUNIT = 1
      ILOGG = 2
!     Reduced mass
      UM = (AM*BM)/(D2*BM+AM)
!
!     Read points
!     ===========
!
      WRITE(LUPRI,'(/A,F12.4)')                                         &
    '* Asymmetric stretch  <-- M <--. Reduced mass:',UM
      WRITE(LUPRI,'(A)') 'Symmetry coordinate: [1/SQRT(2)](r12-r32)'
 10   CONTINUE
      WRITE(LUPRI,'(A)') 'Name input file with potential curve(A6)'
      READ(LUSTDIN,'(A)') POTFIL
      WRITE(LUPRI,'(A)') 'Format: ',                                    &
         'lines of  r12(Angstrom)  r32(Angstrom)  energy(Hartree)'
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',             &
      ACCESS='SEQUENTIAL',ERR=20)
      GOTO 30
 20   CONTINUE
      STOP 'Error opening POTFIL'
 30   CONTINUE
      NP = 0
      DO I = 1,MP
        READ(IUNIT,*,END=40) X(I),XB(I),Y(I)
        NP = NP + 1
      ENDDO
 40   CONTINUE
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
      NP2 = NP+NP
      IF(NP2.GT.MP) STOP 'Too many points..'
      WRITE(ILOGG,'(/A/)')                                              &
    '* ASYMMETRIC STRETCH  <-- M <--'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Reduced mass:',UM
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(2F7.4,4X,F20.10)')                                  &
      (X(I),XB(I),Y(I),I=1,NP)
!
!     Generate symmetry coordinate and convert to atomic units
!
      UM = XFAMU*UM
      FCC = D1/SQRT(D2)/XTANG
      J = 0
      DO I = 1,NP
        X(I) = FCC*(X(I)-XB(I))
        IF(ABS(X(I)).GT.DTOL) THEN
          J = J + 1
          X(J+NP) = -X(I)
          Y(J+NP) =  Y(I)
        ENDIF
      ENDDO
      NP2 = NP + J
!
!     Polynomial fit
!     ==============
!
 60   CONTINUE
      WRITE(LUPRI,'(A,I5)') 'Total number of points :',NP2
      WRITE(LUPRI,'(A)') 'Polynomial fit: Give order of polynomial'
      READ(LUSTDIN,*) NORDER
      IF(NP2.LE.NORDER) THEN
        WRITE(LUPRI,'(A,I3)')                                           &
     'Too few points for a polynomial fit of order',NORDER
        GOTO 666
      ENDIF
      IF(NORDER.GT.MORDER) THEN
        WRITE(LUPRI,'(A)') 'Order to large ....'
        GOTO 666
      ENDIF
      NDIM = NORDER + 1
      CALL POLSVD(NDIM,NP2,A,B,X,Y,C,D,CHISQ,ISKIP)
!* Estimate fit:
      WRITE(LUPRI,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(3X,A)') '* Coefficients:'
      WRITE(LUPRI,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',           &
       C(I),I=1,NDIM)  
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP2
      WRITE(LUPRI,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
!
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(3X,A)') '* Coefficients:'
      WRITE(ILOGG,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',           &
       C(I),I=1,NDIM)  
!TROND      DO I = 1,NP
!TROND        YP = (POLVAL(NORDER,C,X(I))-Y(I))/Y(I)
!TROND        WRITE(LUPRI,'(3X,F6.3,4X,E12.4)') XTANG*X(I),YP
!TROND      ENDDO
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP2
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
!
!     Harmonic constant
!     ==================
!
      WRITE(LUPRI,'(A)') '* Do you want harmonic constant ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 100
      XOPT = D0
      V2  = POL2DER(NORDER,C,XOPT)
      FORCE = V2*XTJ/XTM/XTM
      FREQ  = SQRT(V2/UM)
      WRITE(LUPRI,'(A,2(F14.4,A))') '* Force constant:',FORCE,' N/m = ',&
     FORCE*DMC, ' mDyne/A'
      WRITE(LUPRI,'(A,E10.3,A,3X,F14.4,A)') '* Frequency:',FREQ*XTHZ,   &
     'Hz --> ', FREQ*XTCM,'cm-1'
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,F7.4)') '* Req = ',REQ
      WRITE(ILOGG,'(A,2(F14.4,A))') '* Force constant:',FORCE,'N/m = ', &
     FORCE*DMC, ' mDyne/A'
      WRITE(ILOGG,'(A,E10.3,A,3X,F14.4,A)') '* Frequency:',FREQ*XTHZ,   &
     'Hz --> ', FREQ*XTCM,'cm-1'
!
!     Derivatives of curve
!     ====================
!
 100  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want derivatives of curve (y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
      WRITE(LUPRI,'(A)') 'Derivatives to what order(max.6) ?'
      READ(LUSTDIN,*) NDER
      DO J = 1,NP2
        WRITE(ILOGG,'(F6.3,3X,1P,6E14.6)')                              &
    X(J),(POLNDER(NORDER,B,X(J),I),I = 1,NDER)
      ENDDO
      
 666  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want a new polynomial fit ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.EQ.'y') GOTO 60
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ANGDEF(A,B,C,D,X,Y,MORDER,MP,MITER,AM,BM,IBUF)
      use codata
#include "implicit.h"
#include "priunit.h"
!
      PARAMETER(ANGM = 1.0D-10,XTCM=2.19474625D+05)
      PARAMETER(DTOL = 1.0D-4,DM1 = -1.0D0,D1 = 1.0D0,D2 = 2.0D0)
      PARAMETER(D0 = 0.0D0,DMC = 0.01D0)
      PARAMETER(D2PI = D2*PI,CNV = PI/360.0D0)
      PARAMETER(NITER = 50)
      DIMENSION A(*),B(*),C(*),D(*),IBUF(*)
      DIMENSION X(MP),Y(MP)
      CHARACTER POTFIL*6,REPLY*1
      DOUBLE PRECISION XTM,AMUKG
!
      call set_codata_values(CODSET)
      XTM   = XTANG*ANGM
      AMUKG = UMASS
!
      IUNIT = 1
      ILOGG = 2
!     Reduced mass
      UM = (AM*BM)/(D2*BM+AM)
!
!     Read points
!     ===========
!
      WRITE(LUPRI,'(A,F12.4)')                                          &
    '* Bending mode. Reduced mass:',UM
      WRITE(LUPRI,'(A)') '  Symmetry coordinate: SQRT(2)*req*cos[phi/2]'
 10   CONTINUE
      WRITE(LUPRI,'(A)') 'Name input file with potential curve(A6)'
      READ(LUSTDIN,'(A)') POTFIL
      WRITE(LUPRI,'(A)') 'Format: ',                                    &
         'one line with r_eq (Angstrom)',                               &
         'lines of  phi(degree)  energy(Hartree)'
      WRITE(LUPRI,'(A)') 'Note that bending angle phi should be',       &
           'in degrees and energies in Hartrees !'
      OPEN(IUNIT,FILE=POTFIL,STATUS='OLD',FORM='FORMATTED',             &
      ACCESS='SEQUENTIAL',ERR=20)
      GOTO 30
 20   CONTINUE
      STOP 'Error opening POTFIL'
 30   CONTINUE
      READ(IUNIT,*) REQ
      NP = 0
      DO I = 1,MP
        READ(IUNIT,*,END=40) X(I),Y(I)
        NP = NP + 1
      ENDDO
 40   CONTINUE
      WRITE(LUPRI,'(A,I5)') '* Number of points read:',NP
      NP2 = NP+NP
      IF(NP2.GT.MP) STOP 'Too many points..'
      WRITE(ILOGG,'(/A/)') '* BENDING:'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Reduced mass:',UM
      WRITE(ILOGG,'(A,F6.3)') '*Using req: ',REQ
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,A6)') '* Points read from ',POTFIL
      WRITE(ILOGG,'(F8.3,4X,F20.10)') (X(I),Y(I),I=1,NP)
!
!     Generate symmetry coordinate and convert to atomic units
!
      UM = XFAMU*UM
      FAC = SQRT(D2)/XTANG
      J  = 0
      DO I = 1,NP
        X(I) = FAC*REQ*COS(CNV*X(I))
        IF(ABS(X(I)).GT.DTOL) THEN
          J = J + 1
          X(J+NP) = -X(I)
          Y(J+NP) =  Y(I)
        ENDIF
      ENDDO
      NP2 = NP + J
!
!     Polynomial fit
!     ==============
!
 60   CONTINUE
      WRITE(LUPRI,'(A,I5)') 'Total number of points :',NP2
      WRITE(LUPRI,'(A)') 'Polynomial fit: Give order of polynomial'
      READ(LUSTDIN,*) NORDER
      IF(NP2.LE.NORDER) THEN
        WRITE(LUPRI,'(A,I3)')                                           &
     'Too few points for a polynomial fit of order',NORDER
        GOTO 666
      ENDIF
      IF(NORDER.GT.MORDER) THEN
        WRITE(LUPRI,'(A)') 'Order to large ....'
        GOTO 666
      ENDIF
      NDIM = NORDER + 1
      CALL POLSVD(NDIM,NP2,A,B,X,Y,C,D,CHISQ,ISKIP)
!* Estimate fit:
      WRITE(LUPRI,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(LUPRI,'(3X,A)') '* Coefficients:'
      WRITE(LUPRI,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',           &
       C(I),I=1,NDIM)  
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(LUPRI,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP2
      WRITE(LUPRI,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
!
      WRITE(ILOGG,'(72A1)') ('=',I=1,72)
      WRITE(ILOGG,'(3X,A,I3)') '* Polynomial fit of order:',NORDER
      WRITE(ILOGG,'(3X,A)') '* Coefficients:'
      WRITE(ILOGG,'(3X,A,I3,A,1P,E14.6)') ('c(',(I-1),'):  ',           &
       C(I),I=1,NDIM)  
!TROND      DO I = 1,NP
!TROND        YP = (POLVAL(NORDER,C,X(I))-Y(I))/Y(I)
!TROND        WRITE(LUPRI,'(3X,F6.3,4X,E12.4)') XTANG*X(I),YP
!TROND      ENDDO
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square :  ',CHISQ
      WRITE(ILOGG,'(3X,A,E9.4)') '* Chi square per point:  ',CHISQ/NP2
      WRITE(ILOGG,'(3X,A,I3)') '* Number of singularities(SVD): ',ISKIP
!
!     Find minimum
!     ============
!
      WRITE(LUPRI,'(A)') '* Do you want to find a minimum ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 100
      IYOPT = IDMIN(NP,Y,1)
      XOPT = X(IYOPT)
      CALL NEWRAP(NORDER,DTOL,C,XOPT,NITER,IERR)
      IF(IERR.EQ.1) THEN
        WRITE(ILOGG,'(A)') 'Newtons method failed'
        STOP 'Did not find minimum. Newtons method failed...'
      ELSE
        YOPT = POLVAL(NORDER,C,XOPT)
        RPHI = ACOS(XOPT/(FAC*REQ))/CNV
        WRITE(LUPRI,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',     &
         RPHI,',',YOPT,')'
        WRITE(ILOGG,'(72A1)') ('-',I=1,72)
        WRITE(ILOGG,'(A,F17.10,A,F17.10,A)') '* Local minimum:  (',     &
         RPHI,',',YOPT,')'
      ENDIF
!
!     Harmonic constant
!     ==================
!
      WRITE(LUPRI,'(A)') '* Do you want harmonic constant ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 100
      V2  = POL2DER(NORDER,C,XOPT)
!TROND
      WRITE(6,*) 'test..',XOPT,V2      
      FORCE = V2*XTJ/XTM/XTM
      FREQ  = SQRT(V2/UM)
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,F7.2)') '* phieq = ',RPHI
      WRITE(LUPRI,'(A,2(F14.4,A))') '* Force constant:',FORCE, 'N/m = ',&
     FORCE*DMC, ' mDyne/A'
      WRITE(LUPRI,'(A,E10.3,A,3X,F14.4,A)') '* Frequency:',FREQ*XTHZ,   &
     'Hz --> ', FREQ*XTCM,'cm-1'
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      WRITE(ILOGG,'(A,2(F14.4,A))') '* Force constant:',FORCE,'N/m = ', &
     FORCE*DMC, ' mDyne/A'
      WRITE(ILOGG,'(A,E10.3,A,3X,F14.4,A)') '* Frequency:',FREQ*XTHZ,   &
     'Hz --> ', FREQ*XTCM,'cm-1'
!
!     Derivatives of curve
!     ====================
!
 100  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want derivatives of curve (y/n) ?'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.NE.'y') GOTO 666
      WRITE(LUPRI,'(A)') 'Derivatives to what order(max.6) ?'
      READ(LUSTDIN,*) NDER
      DO J = 1,NP
        WRITE(ILOGG,'(F6.3,3X,1P,6E14.6)')                              &
    X(J),(POLNDER(NORDER,B,X(J),I),I = 1,NDER)
      ENDDO
      
 666  CONTINUE
      WRITE(LUPRI,'(A)') '* Do you want a new polynomial fit ?(y/n)'
      READ(LUSTDIN,'(A1)') REPLY
      IF(REPLY.EQ.'y') GOTO 60
      END
