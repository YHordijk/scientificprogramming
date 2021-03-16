      PROGRAM ORBPRI
C*****************************************************************************
C
C     Utility program for reading orbital-file in GRASP and print radial
C     functions P and Q to formatted file for visualization.
C
C     Written by T.Saue July 6 1995
C
C*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(NW=40, NPT=4000)
      PARAMETER(LUIN=5,LUPRI=6,LUORB=10)
C
      CHARACTER REPLY*1,NH*2,FILNAM*8,ATOM*3
      DIMENSION NH(NW),NP(NW),NAK(NW),E(NW),PZ(NW),QZ(NW),XCAMAX(NW),
     &          MCOW(NW),COWMAN(NW),PF(NPT,NW),QF(NPT,NW),R(NPT)
C*****************************************************************************
C
C     NH   - Orbital label
C     NP   - Principal quantum number
C    NAK   - Angular(kappa) quantum number
C     N    - number of grid points
C     Z    - atomic number
C     H    - step size in logarithmic variable t
C     R    - values of radial variable at grid points
C    RNT   - first grid point in Bohrs
C     E    - Orbital eigenvalue
C     PZ   - leading term in small r expansion of PF
C     QZ   - leading term in small r expansion of QF
C   XCAMAX - maximu exchange value
C    MCOW  - flag for Cowan and Mann modification
C   COWMAN - Mann and Cowna factor
C     PF   - large component of wavefunction
C     QF   - small component of wavefunction
C
C*****************************************************************************
C
C     Read orbital file
C
      WRITE(LUPRI,'(A)') '*** OUTPUT from ORBPRI ***'
      WRITE(LUPRI,'(A)') 'Give unit number of orbital file to read'
      READ(LUIN,*) IORB
      OPEN(IORB,STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL')
      IW = 1
 10   CONTINUE
      READ (IORB,END=20) NH(IW),NP(IW),NAK(IW),N,Z,H,RNT,E(IW),PZ(IW),
     &             QZ(IW),XCAMAX(IW),MCOW(IW),COWMAN(IW)
      READ (IORB) (PF(I,IW),I = 1,N),(QF(I,IW),I = 1,N)
      IW = IW + 1
      GOTO 10
 20   CONTINUE
      IW = IW - 1
      CLOSE(LUORB,STATUS='KEEP')
C
C     Give status of file
C
      WRITE(LUPRI,'(A)') '*** OUTPUT from ORBPRI ***'
      WRITE(LUPRI,'(A,F9.1)') '* Atomic number:          ',Z
      WRITE(LUPRI,'(A,I9)')   '* Number of grid points:  ',N
      WRITE(LUPRI,'(A,E9.3)') '* First grid point:       ',RNT
      WRITE(LUPRI,'(A,E9.3)') '* Grid step size :        ',H
      WRITE(LUPRI,'(A)') '   '
      WRITE(LUPRI,'(A,I9)')   '* Number of orbitals read:',IW
      WRITE(LUPRI,'(A)') '   '
      WRITE(LUPRI,'(A)') 'Do you want info on orbitals (y/n) ?'
      READ(LUIN,'(A1)') REPLY
      IF(REPLY.EQ.'N'.OR.REPLY.EQ.'n') GOTO 999
      WRITE(LUPRI,'(A)')
     &'Orb.        E                   P0                     Q0'
      DO I = 1,IW
        WRITE(LUPRI,'(I2,1X,I2,A2,3X,3F20.12)') I,NP(I),NH(I),-E(I),
     &        PZ(I),QZ(I)
      ENDDO
C
C     Ask for printing of orbitals
C
      WRITE(LUPRI,'(A)')
     &  'Select one of the following:',
     &  '  1. Print radial functions P and Q.',
     &  '  2. Print radial functions R(large) and R(small).',
     &  '  3. Print radial distribution of individual orbitals.',
     &  '  4. Exit.'
      READ(LUIN,*) IOPT
      IF(IOPT.EQ.1) THEN
        FILNAM(7:8)='PQ'
      ELSEIF(IOPT.EQ.2) THEN
        FILNAM(7:8)='RR'
      ELSEIF(IOPT.EQ.3) THEN
        FILNAM(7:8)='rh'
      ELSE
        GOTO 999
      ENDIF
C
C     Atom name
C
      WRITE(LUPRI,'(A)') 'Name atom (A3)'
      READ(LUIN,'(A3)') ATOM
      IF(ATOM(2:2).EQ.' ') ATOM(2:2) = '_'
      IF(ATOM(3:3).EQ.' ') ATOM(3:3) = '_'
C
C
C     Set up grid
C
      EPH = EXP(H)
      R(1) = RNT
      DO I = 2,N
         R(I) = R(I-1)*EPH
      ENDDO
      WRITE(LUPRI,'(A,F9.3)') 'Outer radial point = :',R(N)
C
C     Request what orbital to print
C
 30   CONTINUE
      WRITE(LUPRI,'(A)') 'Give index of orbital to print(-1 to quit):'
      READ(LUIN,*) JW
      IF(JW.LT.1.OR.JW.GT.IW) GOTO 999
      IDIG = ICHAR('0') + NP(JW)
      FILNAM(1:6) = ATOM//CHAR(IDIG)//NH(JW)
      IF(FILNAM(6:6).EQ.' ') FILNAM(6:6)='_'
      OPEN(LUORB,FILE=FILNAM,STATUS='UNKNOWN',ACCESS='SEQUENTIAL',
     &     FORM='FORMATTED')
C.....Print radial functions P and Q
      IF(IOPT.EQ.1) THEN
        PMAX = 0.0D0
        PMIN = 0.0D0
        QMAX = 0.0D0
        QMIN = 0.0D0
        DO I = 1,N
          PMAX = MAX(PMAX,PF(I,JW))
          PMIN = MIN(PMIN,PF(I,JW))
          QMAX = MAX(QMAX,QF(I,JW))
          QMIN = MIN(QMIN,QF(I,JW))
          WRITE(LUORB,'(F24.16,2E24.16)') R(I),PF(I,JW),QF(I,JW)
        ENDDO
        WRITE(LUPRI,'(2(3X,A,E10.4))') 'PMIN :',PMIN,'PMAX :',PMAX
        WRITE(LUPRI,'(2(3X,A,E10.4))') 'QMIN :',QMIN,'QMAX :',QMAX
C.....Print radial functions R(large) and R(small)
      ELSEIF(IOPT.EQ.2) THEN
        PMAX = 0.0D0
        PMIN = 0.0D0
        QMAX = 0.0D0
        QMIN = 0.0D0
        DO I = 1,N
          RL = PF(I,JW)/R(I)
          RS = QF(I,JW)/R(I)
          PMAX = MAX(PMAX,RL)
          PMIN = MIN(PMIN,RL)
          QMAX = MAX(QMAX,RS)
          QMIN = MIN(QMIN,RS)
          WRITE(LUORB,'(F24.16,2E24.16)') R(I),RL,RS
        ENDDO
        WRITE(LUPRI,'(2(3X,A,E10.4))') 'RLMIN :',PMIN,'RLMAX :',PMAX
        WRITE(LUPRI,'(2(3X,A,E10.4))') 'RSMIN :',QMIN,'RSMAX :',QMAX
C.....Print density of individual orbitals
      ELSEIF(IOPT.EQ.3) THEN
        DO I = 1,N
          RC  = R(I)
          RHO = 2.0D0*(PF(I,JW)*PF(I,JW)+QF(I,JW)*QF(I,JW))
          WRITE(LUORB,'(F24.16,E24.16)') RC,RHO
        ENDDO
      ENDIF
      CLOSE(LUORB,STATUS='KEEP')
      GOTO 30
C
C     Exit program
C
 999  CONTINUE
      WRITE(LUPRI,'(A)') '*** EXITING ORBPRI ***'
      END
