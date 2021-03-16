      PROGRAM BASIS_GBS
C***********************************************************************
C
C     Purpose:
C     Read an uncontracted basis set on DALTON form from formatted file
C     and output it to GBS format
C
C***********************************************************************
      CHARACTER*20 BASFIL
      INTEGER IUNIT
      WRITE(6,*) 'Give name of basis file to read (A20)'
      READ(5,'(A20)') BASFIL
      IUNIT=1
      OPEN(IUNIT,FILE=BASFIL,STATUS='OLD',FORM='FORMATTED',             &
      ACCESS='SEQUENTIAL',ERR=20)
      READ (IUNIT)
      READ (IUNIT)
      
 20   CONTINUE
      WRITE(6,'(A,A20)') 'ERROR reading ',BASFIL
      STOP 'Error reading basis file'
      END
      
