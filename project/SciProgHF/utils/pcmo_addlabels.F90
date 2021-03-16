      PROGRAM PCMO_ADDLABELS
!***********************************************************************
!
!     Add labels to DFPCMO
!
!     Written by T. Saue March 3 2021
!
!***********************************************************************
   implicit none        

   integer, parameter   :: LU1   = 1
   integer, parameter   :: LU2   = 2
   integer, parameter   :: LUPRI = 6
   character(len=80)    :: ARGC
   character(len=74)    :: TEXT
   real(8)              :: TOTERG
   integer              :: IDIM(3,2),NFSYM,NC,NO,IFRP,NORB,NZ,IERR
   integer              :: NARG,JARG,I,J
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
   OPEN(LU1,FILE=ARGC(1:JARG),STATUS='UNKNOWN',FORM='FORMATTED')
   OPEN(LU2,FILE=ARGC(1:JARG)//'.new',STATUS='NEW',FORM='FORMATTED')   
   READ(LU1,'(A74)',IOSTAT=IERR) TEXT
   IF (IERR /= 0) CALL QUIT('error reading DFPCMO: title line')   
   READ(LU1,*,IOSTAT=IERR) NFSYM,((IDIM(I,J),I = 1,3),J=1,NFSYM)
   IF (IERR /= 0) CALL QUIT('error reading DFPCMO: symmetry information')
   READ(LU1,'(E24.16)',IOSTAT=IERR) TOTERG
   IF (IERR /= 0) CALL QUIT('error reading DFPCMO: energy')
! Calculate dimensions of coefficients and eigenvalues
   NC = 0
   NO = 0
   DO IFRP = 1,NFSYM
     NORB = IDIM(1,IFRP)+IDIM(2,IFRP)
     NO   = NO + NORB
     NC   = NC + NORB*IDIM(3,IFRP)
   ENDDO
!  We do not know NZ; ask user
   WRITE(6,*) 'Give parameter NZ'
   READ(5,*) NZ
!
   WRITE(LU2,'(A8)') 'INFO    '   
   WRITE(LU2,'(A74)') TEXT
   WRITE(LU2,'(8(X,I0))') NFSYM,NZ,((IDIM(I,J),I = 1,3),J=1,NFSYM)
   WRITE(LU2,'(E24.16)') TOTERG
   NC = NC*NZ
!  Coefficients  
   allocate(CMO(NC))
   READ(LU1,'(6F22.16)',iostat=ierr) CMO
   IF (IERR /= 0) CALL QUIT('error reading DFPCMO: coefficients')
   WRITE(LU2,'(A8)') 'COEFS   '      
   WRITE(LU2,'(6F22.16)') CMO
   deallocate(CMO)
!  Eigenvalues
   allocate(EIG(NO))
   READ(LU1,'(6E22.12)',iostat=ierr) EIG
   IF (IERR /= 0) CALL QUIT('error reading DFPCMO: eigenvalues')
   WRITE(LU2,'(A8)') 'EVALS   '            
   WRITE(LU2,'(6E22.12)') EIG
   deallocate(EIG)
!  Supersymmetry
   allocate(IBEIG(NO))   
   READ(LU1,*,iostat=ierr) IBEIG
   IF (IERR /= 0) CALL QUIT('error reading DFPCMO: ibeig')
   WRITE(LU2,'(A8)') 'SUPERSYM'
   WRITE(LU2,'(66(X,I0))') IBEIG
   CLOSE(LU1)
   CLOSE(LU2)
   WRITE(6,*) 'In the atomic case KAPPA-information must be added.'
   WRITE(6,*) 'The easiest solution in this case is to regenerate DFPCMO from scratch.'
 END PROGRAM PCMO_ADDLABELS
