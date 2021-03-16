      PROGRAM CF_ADDLABELS
!***********************************************************************
!
!     Add labels to DFCOEF
!
!     Written by T. Saue Jan 27 2021
!
!***********************************************************************
   implicit none        

   integer, parameter   :: LU1   = 1
   integer, parameter   :: LU2   = 2
   integer, parameter   :: LUPRI = 6
   character(len=80)    :: ARGC
   character(len=74)    :: TEXT
   real(8)              :: TOTERG
   integer              :: IDIM(3,2),NFSYM,NC,NO,IFRP,NORB,NZ,I,J
   integer              :: NARG,JARG,POS,FTELL
   real(8), allocatable :: CMO(:),EIG(:)
   integer, allocatable :: IBEIG(:)
   narg = command_argument_count()
   IF (NARG .NE. 1) THEN
      CALL get_command_argument(0,ARGC)
      JARG = LEN_TRIM(ARGC)
      WRITE(6,'(/3A/)') 'Usage: ',ARGC(1:JARG),' <coefficientfile>'
      STOP 'Wrong number of arguments'
   END IF
   CALL get_command_argument(1,ARGC)
   JARG = LEN_TRIM(ARGC)
   OPEN(LU1,FILE=ARGC(1:JARG),STATUS='UNKNOWN',FORM='UNFORMATTED')
   OPEN(LU2,FILE=ARGC(1:JARG)//'.new',STATUS='NEW',FORM='UNFORMATTED')   
   READ(LU1,END=10,ERR=20) TEXT,NFSYM,((IDIM(I,J),I = 1,3),J=1,NFSYM),TOTERG
! Calculate dimensions of coeffcients and eigenvalues
   NC = 0
   NO = 0
   DO IFRP = 1,NFSYM
     NORB = IDIM(1,IFRP)+IDIM(2,IFRP)
     NO   = NO + NORB
     NC   = NC + NORB*IDIM(3,IFRP)
   ENDDO
!  We do not know NZ; use ftell to find out
   NZ = ftell(LU1)   
   READ(LU1)
   NZ = (ftell(LU1)-NZ-8)/8/NC
   BACKSPACE LU1
   CALL NEWLAB('INFO    ',LU2,LUPRI)      
   WRITE(LU2) TEXT,NFSYM,NZ,((IDIM(I,J),I = 1,3),J=1,NFSYM),TOTERG
   NC = NC*NZ
!  Coefficients  
   allocate(CMO(NC))
   CALL READT(LU1,NC,CMO)
   CALL NEWLAB('COEFS   ',LU2,LUPRI)      
   WRITE(LU2) CMO
   deallocate(CMO)
!  Eigenvalues
   allocate(EIG(NO))
   CALL READT(LU1,NO,EIG)
   CALL NEWLAB('EVALS   ',LU2,LUPRI)            
   WRITE(LU2) EIG
   deallocate(EIG)
!  Supersymmetry
   allocate(IBEIG(NO))   
   CALL READI(LU1,NO,IBEIG)
   CALL NEWLAB('SUPERSYM',LU2,LUPRI)                  
   WRITE(LU2) IBEIG
   deallocate(IBEIG)
   CALL NEWLAB('EOFLABEL',LU2,LUPRI)                              
   CLOSE(LU1)
   CLOSE(LU2)
   GOTO 30
10 CONTINUE
   CALL QUIT('CF_ADDLABELS: END OF FILE reading TEXT')
20 CONTINUE
   CALL QUIT('CF_ADDLABELS: ERROR reading TEXT')
30 CONTINUE
   
 END PROGRAM CF_ADDLABELS
