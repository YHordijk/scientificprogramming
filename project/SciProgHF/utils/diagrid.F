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

      PROGRAM DIAGRID
!***********************************************************************
!
!     Utility program for diatomics:
!       generates even-spaced grid within classical turning points
!
!***********************************************************************
      use codata
#include "implicit.h"
#include "priunit.h"
      INTEGER ILOGG
      CHARACTER PROJECT*40,REPLY*1
      ILOGG = 1
!
!     Heading
!
      CALL TITLER('DIAGRID for diatomics : Written by T. Saue ','*',110)
      WRITE(LUPRI,'(A)') 'Give project name (A40)'
      READ(LUSTDIN,*) PROJECT
      i=LNBLNK(PROJECT)
      OPEN(ILOGG,FILE=PROJECT(1:i)//'.diagrid',STATUS='UNKNOWN',          &
      FORM='FORMATTED',ACCESS='SEQUENTIAL')
      WRITE(ILOGG,*) '* Input data:'
!
!     Equilibrium bond length
!     =======================
!
      WRITE(LUPRI,'(A)') 'Give equilibrium bond length in Angstroms: '
      READ(LUSTDIN,*) REQ
      WRITE(LUPRI,'(A)') 'Give harmonic frequency in cm-1'
      READ(LUSTDIN,*) WEQ
      WRITE(LUPRI,'(A)') 'Select vibrational quantum number'
      READ(LUSTDIN,*) NU
      WRITE(ILOGG,*) ' - Equilibrium bond length (Angstrom)   ',REQ
      WRITE(ILOGG,*) ' - Harmonic frequency (cm-1)            ',WEQ
      WRITE(ILOGG,*) ' - Selected vibrational quantum number  ',NU
!
!     Get masses
!     ==========
!
 10   CONTINUE
      WRITE(LUPRI,'(A)')                                                &
   'Select one of the following:',                                      &
   '  1. Select masses of the most abundant isotopes.',                 &
   '  2. Employ user-defined atomic masses.',                           &
   '  3. Normal modes: Give reduced mass.'
      READ(LUSTDIN,*) I
      IF(I.EQ.1) THEN
        WRITE(LUPRI,'(A)') '* Give charge of atom A:'
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
        WRITE(LUPRI,'(A)') '* Give charge of atom B :'
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
        WRITE(LUPRI,'(A)') '* Give mass of atom A in Daltons:'
        READ(LUSTDIN,*) AM
        WRITE(LUPRI,'(A)') '* Give mass of atom B in Daltons:'
        READ(LUSTDIN,*) BM
      ELSEIF(I.EQ.3) THEN
        WRITE(LUPRI,'(A)') '* Give reduced mass in Daltons:'
        READ(LUSTDIN,*) UM
        GOTO 50
      ELSE
        WRITE(LUPRI,'(A)') ' You stupid fart ! Try again !'
        GOTO 10
      ENDIF
      UM = D1/(D1/AM + D1/BM)
      WRITE(ILOGG,'(A)') '* MASSES:'
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom A      :',AM
      WRITE(ILOGG,'(3X,A,F8.4)') 'Atom B      :',BM
 50   CONTINUE
      WRITE(ILOGG,*) ' - Reduced mass:',UM
      WRITE(ILOGG,'(72A1)') ('-',I=1,72)
      UM = XFAMU*UM
      CLOSE(ILOGG,STATUS='KEEP')
      END
