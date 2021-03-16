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
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!
!      ****  Diagonalization Routines for the RELADC package  ****
!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
! this is the real version of a full diagonaliation routine
!
      SUBROUTINE FULLDIAR(IOUNIT, FILEN4, FILEN5, INTBUF, LADC, KREP)

      implicit none

      INTEGER                       :: IOUNIT,INTBUF,LADC,KREP
      CHARACTER(len=6)              :: FILEN4, FILEN5
!
!---------------Description--------------------------------------------
!
!  For not too large symmetric matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      REAL*8                          :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      REAL*8, allocatable             :: A(:,:),WORK(:),EVL(:),BUF(:)
      INTEGER, allocatable            :: IOI(:),IOJ(:)
      INTEGER                         :: NBUFS,LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=9)                :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering FULL diagonalizer (debug)+')

      IF(LADC.gt.5000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Skipping this step...'
        return
      ENDIF

      LUN_DIA=IOUNIT
      LUN_OFF=IOUNIT+1
      LUN_OUT=IOUNIT+2




      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A7,I2)') 'FULLEV.',KREP
      ELSE
        WRITE(longfilename,'(A6,A2,I1)') 'FULLEV','.0',KREP
      ENDIF

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(work(lwork))
      
      A=0.0d0
!_________________________ reading section
!|
!|
!|
!|
      OPEN(LUN_DIA,FILE=FILEN4,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_DIA); REWIND(LUN_OFF)
      READ(LUN_DIA,ERR=88,END=99) (A(IXX,IXX),IXX=1,LADC),NBUFS
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GE.IROW) CALL QUIT('ICOL.GE.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_DIA)
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling DSYEV with lwork=',lwork
      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )
!     write(*,*) 'AFTER  DSYEV with optimal lwork=',work(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAR (lapack dsyev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential', 
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)

!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(work)

      RETURN

 88   STOP 'Error reading the diagonal (FULLDIAR)'
 99   STOP 'Unexpected end of diagonal file (FULLDIAR)'
 199  STOP 'Error reading off-diagonal elements (FULLDIAR)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAC(IOUNIT, FILEN4, FILEN5, INTBUF, LADC, KREP)
      implicit none

      INTEGER                       :: IOUNIT,INTBUF,LADC,KREP
      CHARACTER(len=6)              :: FILEN4, FILEN5
!
!---------------Description--------------------------------------------
!
!  For not too large Hermitian matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      COMPLEX*16                      :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      COMPLEX*16,   allocatable       :: A(:,:),CWORK(:),BUF(:)
      REAL*8,       allocatable       :: EVL(:),RWORK(:)
      INTEGER,      allocatable       :: IOI(:),IOJ(:)
      INTEGER                         :: NBUFS,LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=9)                :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering FULL diagonalizer (debug)+')

      IF(LADC.gt.5000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Skipping this step...'
        return
      ENDIF

      LUN_DIA=IOUNIT
      LUN_OFF=IOUNIT+1
      LUN_OUT=IOUNIT+2



      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A7,I2)') 'FULLEV.',KREP
      ELSE
        WRITE(longfilename,'(A6,A2,I1)') 'FULLEV','.0',KREP
      ENDIF

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(rwork(lwork))
      allocate(cwork(lwork))
      
      A=(0.0d0,0.0d0)
!_________________________ reading section
!|
!|
!|
!|
      OPEN(LUN_DIA,FILE=FILEN4,FORM='UNFORMATTED',STATUS='UNKNOWN')
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_DIA); REWIND(LUN_OFF)
      READ(LUN_DIA,ERR=88,END=99) (A(IXX,IXX),IXX=1,LADC),NBUFS
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GE.IROW) CALL QUIT('ICOL.GE.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_DIA)
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling ZHEEV with lwork=',lwork
      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK,
     &            LWORK, RWORK, INFO )
!     write(*,*) 'AFTER  ZHEEV with optimal lwork=',cwork(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential',
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,2f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)
!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(rwork)
      deallocate(cwork)

      RETURN

 88   STOP 'Error reading the diagonal (FULLDIAC)'
 99   STOP 'Unexpected end of diagonal file (FULLDIAC)'
 199  STOP 'Error reading off-diagonal elements (FULLDIAC)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAR2(IOUNIT, FILEN5, INTBUF, NBUFS, LADC, KREP)

      use memory_allocator
      use adc_fano_exchange
      use adc_mat

      implicit none

      INTEGER                       :: IOUNIT,INTBUF,NBUFS,LADC,KREP
      CHARACTER(len=6)              :: FILEN5
!
!---------------Description--------------------------------------------
!
!  For not too large symmetric matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      REAL*8                          :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      REAL*8, allocatable             :: A(:,:),WORK(:),EVL(:),BUF(:)
      INTEGER, allocatable            :: IOI(:),IOJ(:)
      INTEGER                         :: LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=10)               :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering (new format) REAL FULL diagonalizer+')
      write(*,*) 'Number of buffers to read:',nbufs

      IF (reladc_md_isfano) THEN
        IF(LADC.gt.8000) THEN
        write(*,*) 'Matrix dimension too large for direct'
        write(*,*) 'diagonalization! Exiting from this...'
        CALL QUIT('No full diagonalization of final states.')
        END IF
      ELSE
        IF(LADC.gt.5000) THEN
          write(*,*) 'Matrix dimension too large for direct'
          write(*,*) 'diagonalization! Exiting from this...'
          return
        END IF
      ENDIF

      LUN_OFF=IOUNIT
      LUN_OUT=IOUNIT+1



      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A8,I2)') 'FULLEVD.',KREP
      ELSE
        WRITE(longfilename,'(A7,A2,I1)') 'FULLEVD','.0',KREP
      ENDIF

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(work(lwork))
      
      A=0.0d0
!_________________________ reading section
!|
!|
!|
!|
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_OFF)
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)
!         write(*,'(A,2I4,F20.12)') 'Q<--',irow,icol,S

          IF(ICOL.GT.IROW) CALL QUIT('ICOL.GT.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section


!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling DSYEV with lwork=',lwork
      CALL DSYEV( JOBZ, UPLO, LADC, A, LADC, EVL, WORK, LWORK, INFO )
!     write(*,*) 'AFTER  DSYEV with optimal lwork=',work(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in DSYEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAR (lapack dsyev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 20):'
      DO I=1,min(20,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential',
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)
!
!  in case of a Fano run, write eigenvectors to array
!
      IF(reladc_md_isfano) THEN
        CALL alloc(fin_evecs,ladc,ladc,id='array of Fano final states')
        CALL alloc(fin_energies,ladc,id='2h1p energies of final states')
        WRITE(*,*) 'Saving final state vectors and energies to arrays'
        
        fin_energies(1:ladc)     = evl(1:ladc)
        fin_evecs(1:ladc,1:ladc) = a(1:ladc,1:ladc)

      END IF
!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(work)

      RETURN

 199  STOP 'Error reading off-diagonal elements (FULLDIAR)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FULLDIAC2(IOUNIT, FILEN5, INTBUF, NBUFS, LADC, KREP)

      use memory_allocator
      use adc_fano_exchange
      use adc_mat

      implicit none

      INTEGER                       :: IOUNIT,INTBUF,NBUFS,LADC,KREP
      CHARACTER(len=6)              :: FILEN5
!
!---------------Description--------------------------------------------
!
!  For not too large Hermitian matrices the full diagonalization is
!  performed for debugging purposes. The required arrays are allocated
!  dynamically in order to achieve maximal independence.
!
!  The only specific entities are two available file handles,
!  the length of the write chunks and the length of the ADC matrix.  
!
!---------------Local variables--------------------------------------
!
      COMPLEX*16                      :: S
      CHARACTER(len=1)                :: JOBZ,UPLO
      REAL*8                          :: AUTOEV = 27.2113957D0
      COMPLEX*16,   allocatable       :: A(:,:),CWORK(:),BUF(:)
      REAL*8,       allocatable       :: EVL(:),RWORK(:)
      INTEGER,      allocatable       :: IOI(:),IOJ(:)
      INTEGER                         :: LWORK,TOTELE
      INTEGER                         :: I,J,K,IXX,IREC,NACT,JDUMMY
      INTEGER                         :: IROW,ICOL,INFO
      INTEGER                         :: LUN_DIA,LUN_OFF,LUN_OUT
      CHARACTER(len=10)                :: longfilename
!
!---------------Executable code--------------------------------------
!
!  do the allocations and
!  construct ADC matrix by reading nonzero elements from file
!

      CALL PST('Entering (new format) FULL complex diagonalizer+')
      write(*,*) 'Number of buffers to read:',nbufs

      IF (reladc_md_isfano) THEN                             
        IF(LADC.gt.8000) THEN                                
        write(*,*) 'Matrix dimension too large for direct'   
        write(*,*) 'diagonalization! Exiting from this...'   
        CALL QUIT('No full diagonalization of final states.')
        END IF                                               
      ELSE
        IF(LADC.gt.5000) THEN                                
          write(*,*) 'Matrix dimension too large for direct' 
          write(*,*) 'diagonalization! Exiting from this...' 
          return                                             
        END IF                                               
      ENDIF

      LUN_OFF=IOUNIT
      LUN_OUT=IOUNIT+1



      IF(KREP.GT.9) THEN
        WRITE(longfilename,'(A8,I2)') 'FULLEVD.',KREP
      ELSE
        WRITE(longfilename,'(A7,A2,I1)') 'FULLEVD','.0',KREP
      ENDIF

      allocate(a(ladc,ladc))
      allocate(evl(ladc))
      allocate(buf(intbuf))
      allocate(ioi(intbuf))
      allocate(ioj(intbuf))
      lwork=10*ladc
      allocate(rwork(lwork))
      allocate(cwork(lwork))
      
      A=(0.0d0,0.0d0)
!_________________________ reading section
!|
!|
!|
!|
      OPEN(LUN_OFF,FILE=FILEN5,FORM='UNFORMATTED',STATUS='UNKNOWN')

      REWIND(LUN_OFF)
      TOTELE = 0
      DO IREC = 1,NBUFS
        READ(LUN_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                        (IOI(IXX),IXX=1,INTBUF),
     &                        (IOJ(IXX),IXX=1,INTBUF),
     &                         NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          S    = BUF(K)

          IF(ICOL.GT.IROW) CALL QUIT('ICOL.GT.IROW should not happen')
 
          A(IROW,ICOL)=S
          TOTELE=TOTELE + 1
        ENDDO
      ENDDO
      WRITE(*,555) TOTELE
      CLOSE(LUN_OFF)
!|
!|
!|
!|________________________ end reading section

!
!  call symmetric diagonalizer
!  calculate eigenvectors being placed in A
!
      JOBZ='V'
      UPLO='L'

!     write(*,*) 'Calling ZHEEV with lwork=',lwork
      CALL ZHEEV( JOBZ, UPLO, LADC, A, LADC, EVL, CWORK,
     &            LWORK, RWORK, INFO )
!     write(*,*) 'AFTER  ZHEEV with optimal lwork=',cwork(1)

      IF(INFO.NE.0) THEN
           write(*,*) 'Error in ZHEEV. Info:',INFO
           CALL QUIT ('*** Error in FULLDIAC (lapack zheev) ***')
      ENDIF

      WRITE(*,*) 'Eigenvalues (lowest 10):'
      DO I=1,min(10,LADC)
        WRITE(*,*) I,EVL(I),EVL(I)*AUTOEV
      ENDDO

      WRITE(*,*) 'Writing to output file  ',longfilename
      open(lun_out,file=longfilename,access='sequential',
     &     form='formatted')
      do i=1,ladc
        write(lun_out,'(A,i6,3x,f17.12,A)')
     &         '----- eigenvalue: ',i,evl(i)*autoev,'  eV'
        do j=1,ladc
          write(lun_out,'(i6,3x,2f17.12)') j,a(j,i)
        enddo
      enddo
      write(lun_out,'(A)') '     end of output'
      close(lun_out)
!
!  in case of a Fano run, write eigenvectors to array
!
      IF(reladc_md_isfano) THEN
        CALL alloc(cfin_evecs,ladc,ladc,id='array of Fano final states')
        CALL alloc(fin_energies,ladc,id='2h1p energies of final states')
        WRITE(*,*) 'Saving final state vectors and energies to arrays'
        
        fin_energies(1:ladc)     = evl(1:ladc)
        cfin_evecs(1:ladc,1:ladc) = a(1:ladc,1:ladc)

      END IF
!
!  release all memory
!
      deallocate(a)
      deallocate(evl)
      deallocate(buf)
      deallocate(ioi)
      deallocate(ioj)
      deallocate(rwork)
      deallocate(cwork)

      RETURN

 199  STOP 'Error reading off-diagonal elements (FULLDIAC)'

 555  FORMAT (1X,'Total number of elements in OD-file:',T40,I8)

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE DIAG_LANC(IW,DOINCORE)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!! this is the main driver routine for the Lanczos routines  !!!!!
!!!!!! now also working with dynamic memory allocation  !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!

      use adc_mat      !obtain matrix parameters via descriptor file
      use adc_cfg

      IMPLICIT INTEGER (A-Z)

      INTEGER                            ::  IW
      LOGICAL                            ::  DOINCORE
C
C---------------Description--------------------------------------------
C
C  This routine executes diagonalizer on any matrix (SIP/DIP) coming
C  
C  Infohandle is the base file handle.
C  the subsequent handles INFOHANDLE+1, +2 ... are available and
C  can be used. The allocation is as follows:
c
C  INFOHANDLE+1:   packed ADC matrix
C  INFOHANDLE+2:   old format: ADC matrix diagonal elements.
C  INFOHANDLE+3:   a free handle for the XLANCZOS routines
C  INFOHANDLE+4:   a free handle

C
C---------------Local variables--------------------------------------
C
c *** for infofile
      INTEGER INFOHANDLE
      INTEGER IONIZLEVEL
      INTEGER INTBUF
      INTEGER KREP
      INTEGER RCW
      INTEGER IOLDNEW
      INTEGER LADC
      CHARACTER*6 FILEADC
      CHARACTER*6 FILEDIA
      CHARACTER*6 FILECNF 
      CHARACTER*5 NMSPEC 
      INTEGER IRECL
      INTEGER NITER
      INTEGER NOCC
      INTEGER NEIGENV
      REAL*8 EVLLIM
      INTEGER IVALID   ! number of valid elements for Lanczos incore diagonalization.
c
c *** others
      LOGICAL LFLAG,CARITH
      COMPLEX*16 CDUMMY
      REAL*8 RDUMMY
      INTEGER NINCWORDS
C
C--------------- common variables--------------------------------------
C
      COMMON/MATCHUNKS/NBUFS
C
C--------------- arrays for dynamical memory allocation, replacing memopts
C
      integer, allocatable, dimension(:)    :: WRK_IDI, WRK_IDJ,
     &          WRK_IND, WRK_IVX, WRK_LSP, WRK_MLT, WRK_OI, WRK_OJ

      real*8, allocatable, dimension(:)     :: WRK_AEV, WRK_BUF,
     &     WRK_DIAG, WRK_EVL, WRK_IDR,  WRK_P, WRK_R,
     &     WRK_RWORK1, WRK_T, WRK_WORK1, WRK_ZEV, WRK_Z
C
C---------------Executable code--------------------------------------
c
c default to a predefined buffer length. Will be changed
c in the future because we need a flexible buffer size,
c
      INTBUF = 5*1024*1024  !output buffer is now 5 MWORDS (40 MByte)
      IF(DOINCORE) THEN
        write(iw,*) 'Incore Lanczos currently not supported.'
        DOINCORE = .false.
      ENDIF
c
c  transfer required parameters from matrix descriptor module
c  to actual modules.
c
      IONIZLEVEL   =  reladc_md_ionizl
      INTBUF       =  reladc_md_intbuf
      INFOHANDLE   =  reladc_md_iobase
      KREP         =  reladc_md_desrep
      RCW          =  reladc_md_rcw
      IOLDNEW      =  reladc_md_ioldnew
      LADC         =  reladc_md_matdim
      FILEADC      =  reladc_md_fileadc
      FILEDIA      =  reladc_md_filediag
      FILECNF      =  reladc_md_filecnf
      NMSPEC       =  reladc_md_nmspec
      IRECL        =  reladc_md_irecl
      NITER        =  reladc_md_lnzitr
      NOCC         =  reladc_md_nmain
      NBUFS        =  reladc_md_nbufs    
      NEIGENV      =  reladc_md_neigenv
      IF (IONIZLEVEL.eq.1) THEN
        EVLLIM     =  reladc_sipprnt
      ELSE IF (IONIZLEVEL.eq.2) THEN
        EVLLIM     =  reladc_dipprnt
      ELSE IF (IONIZLEVEL.eq.3) THEN
        EVLLIM     =  reladc_excprnt
      END IF
c
c  open the ADC Matrix and configuration files and do diagonalization
c
      OPEN(INFOHANDLE+1,
     %     FILE=FILEADC,FORM='UNFORMATTED',STATUS='UNKNOWN')
      IF(IOLDNEW.EQ.1) OPEN(INFOHANDLE+2,
     %     FILE=FILEDIA,FORM='UNFORMATTED',STATUS='UNKNOWN')

c
c this is the absolute lower bound to which we will return as soon as the
c diagonalization has been finished. We must not go further down because
c some buffers are needed in the calling code.
c
      NP = NOCC
      CARITH=.FALSE.
      IF(RCW.EQ.2) CARITH=.TRUE.

      CALL PST('    Initializing iterative diagonalizer+')
!
!  if repacking was successful we work from here with smaller ADC matrix
!  (in the double ionization case).
!
!     IF(DOINCORE) THEN
!       IF(IOLDNEW.eq.1) THEN
!         WRITE(IW,*) 'Incore calculation for the old-format'
!         WRITE(IW,*) 'ADC matrix not possible. Switching back to'
!         WRITE(IW,*) 'Out of core diagonalization'
!         DOINCORE = .FALSE.
!       ELSE
!         WRITE(IW,*) 'Trying (partial) incore diagonalization.'
!       ENDIF
!     ENDIF


      WRITE(IW,*)
      WRITE(IW,*) 'ADC matrix is (1=SIP,2=DIP,3=EXC):',IONIZLEVEL
      WRITE(IW,*) 'Write buffer length is:           ',INTBUF
      WRITE(IW,*) 'Final states have symmetry:       ',KREP
      WRITE(IW,*) 'ADC matrix is (1=real,2=complex): ',RCW
      WRITE(IW,*) 'ADC matrix format (1=old/2=new):  ',IOLDNEW
      WRITE(IW,*) 'ADC matrix dimension (# of rows): ',LADC
      WRITE(IW,*) 'Packed matrix stored in:          ',FILEADC
      WRITE(IW,*) 'Diagonal elem. (old) stored in:   ',FILEDIA
      WRITE(IW,*) 'Configuration data stored in:     ',FILECNF
      WRITE(IW,*) 'Eigenvalues written to:           ',NMSPEC,'.#IREP'
      WRITE(IW,*) 'Record length of Config file:     ',IRECL
      WRITE(IW,*) 'Number of requ. Lanc iterations:  ',NITER
      WRITE(IW,*) 'Half band width (1h/2h/hp states):',NOCC
      WRITE(IW,*) 'Number of chunks in ADC file:     ',NBUFS
      WRITE(IW,*) 'Number of requested eigenvectors: ',NEIGENV
      WRITE(IW,*) 'Printing limit for eigenvalues:   ',EVLLIM
C
C  allocate required work space for band Lanczos iterations
C  real/complex case are accounted for by the RCW variable
C  wrk_diag is always real since we have a hermitian matrix
C  disk buffers are already allocated.
C
C  the OLD/NEW matrix formats are accounted for by the IOLDNEW
C  variable. In this case the file handle for the diagfile carries
C  the number -255
C

      ITAP1=INFOHANDLE+1    ! handle for the packed ADC matrix
      ITAP2=-255
      IF(IOLDNEW.EQ.1) ITAP2=INFOHANDLE+2  ! old matrix format !!


      ITAP3=INFOHANDLE+3    ! free handle
      ITAP5=INFOHANDLE+5    ! free handle
      ITAP6=INFOHANDLE+6    ! free handle

      N=2*NP+1
      allocate(WRK_IVX(N));WRK_IVX = 0

      N=(2*NP+1)*LADC*RCW
      allocate(WRK_P(N));WRK_P = 0.0d0

      N=(NP+1)*NITER*RCW
      allocate(WRK_T(N));WRK_T = 0.0d0

      N=LADC*RCW
      allocate(WRK_R(N));WRK_R = 0.0d0
      allocate(WRK_DIAG(N)); WRK_DIAG = 0.0d0

      N=INTBUF*RCW
      allocate(WRK_BUF(N)); WRK_BUF = 0.0d0

      N=INTBUF
      allocate(WRK_OI(N)); WRK_OI = 0
      allocate(WRK_OJ(N)); WRK_OJ = 0
c
c  remember memory pointers before we do the huge incore parts.
c
!     IF(DOINCORE) THEN
!       KINCSAV = KFREE
!       WRITE(IW,*) 'Memory pointers before incore diag:',
!    &               KFREE,LFREE
!       WRITE(IW,*) 'Saving free memory:',KINCSAV
!       WRITE(IW,*) 'Available buffer for incore diag:'
!       WRITE(IW,*) (LFREE-10),' words'
!       if(RCW.eq.1) then
!          N=(LFREE - 10)/3
!          WRITE(IW,*) 'We can store',N,' real numbers'
!       else
!          N=(LFREE - 10)/4
!          WRITE(IW,*) 'We can store',N,' complex numbers'
!       endif
!       IF(N.lt.INTBUF) THEN
!         WRITE(IW,*) 'Buffer for inc. diag too small!'
!         WRITE(IW,*) 'Switching back to out-of-core diag.'
!         DOINCORE = .FALSE.
!       ELSE
!         CALL MEMOPT('REAL',WRK_IDR,N*RCW,CC,KFREE,LFREE,ALLOC_ME)
!         CALL MEMOPT('INTE',WRK_IDI,N,CC,KFREE,LFREE,ALLOC_ME)
!         CALL MEMOPT('INTE',WRK_IDJ,N,CC,KFREE,LFREE,ALLOC_ME)
!         CALL MEMCHK('ICD',CC,KFRSAV)
!         NINCWORDS = N
!         IF(CARITH) THEN
!           CALL INITINCORE_C(ITAP1,ITAP3,NBUFS,NINCWORDS,INTBUF,
!    &                        CC(WRK_BUF),CC(WRK_OI),CC(WRK_OJ),
!    &                        CC(WRK_IDR),CC(WRK_IDI),CC(WRK_IDJ),
!    &                        IVALID)
!         ELSE
!           CALL INITINCORE_R(ITAP1,ITAP3,NBUFS,NINCWORDS,INTBUF,
!    &                        CC(WRK_BUF),CC(WRK_OI),CC(WRK_OJ),
!    &                        CC(WRK_IDR),CC(WRK_IDI),CC(WRK_IDJ),
!    &                        IVALID)
!         ENDIF
!       ENDIF
!     ENDIF
!
C  enter real resp. complex Band Lanczos iterations
C  for symmetric resp. Hermitian matrices
C  In the first run: no eigenvectors because no eigenvalues are
C  available
C
      LFLAG = .FALSE.

!     IF(DOINCORE.EQV..true.) THEN
!       IF(CARITH.EQV..TRUE.) THEN
!       CALL ZLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
!    &                   CC(WRK_P), CC(WRK_T), CC(WRK_R),
!    &                   CC(WRK_DIAG), CDUMMY, CDUMMY, CC(WRK_BUF),
!    &                   CC(WRK_OI), CC(WRK_OJ), CDUMMY,LFLAG,
!    &                   INTBUF,ITAP1, ITAP2, ITAP3, NINCWORDS,
!    &                   cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),ivalid)
!       CALL MEMCHK('AFTER COMPLEX INCORE-LANCZOS',CC,KFRSAV)
!       ELSE
!       CALL RLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
!    &          CC(WRK_P), CC(WRK_T), CC(WRK_R),
!    &          CC(WRK_DIAG), RDUMMY, RDUMMY, CC(WRK_BUF),
!    &          CC(WRK_OI), CC(WRK_OJ), RDUMMY,LFLAG,
!    &          INTBUF,ITAP1, ITAP2, ITAP3, NINCWORDS,
!    &          cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),ivalid)
!       CALL MEMCHK('AFTER REAL INCORE-LANCZOS',CC,KFRSAV)
!       ENDIF
!       CALL MEMREL('After incore Lanczos',CC,1,KINCSAV,KFREE,LFREE)
!       WRITE(IW,*)  'Memory after incore Lanczos released.'
!     ELSE 
        IF(CARITH.EQV..TRUE.) THEN
        CALL ZLANCZOS(KREP,LADC, NP, NITER, NEIGENV, WRK_IVX,
     &          WRK_P, WRK_T, WRK_R,
     &          WRK_DIAG, CDUMMY, CDUMMY, WRK_BUF,
     &          WRK_OI, WRK_OJ, CDUMMY,LFLAG,
     &          INTBUF,ITAP1, ITAP2, ITAP3)
        ELSE
        CALL RLANCZOS(KREP,LADC, NP, NITER, NEIGENV, WRK_IVX,
     &          WRK_P, WRK_T, WRK_R,
     &          WRK_DIAG, RDUMMY, RDUMMY, WRK_BUF,
     &          WRK_OI, WRK_OJ, RDUMMY,LFLAG,
     &          INTBUF,ITAP1, ITAP2, ITAP3)
        ENDIF
!     ENDIF
C
C  from here the Lanczos matrix T is available. we acquire work
C  space for the LAPACK Band-diagonalizer.
C  ITAP5 will be provided as a free file handle. 
C  It is only used locally and reusable after finishing RBANDDIA/ZBANDDIA
C  This file handle is used for storing the T eigenvector matrix
C
      N=NITER
      allocate(WRK_EVL(N))
      N=NITER*NITER*RCW
      allocate(WRK_ZEV(N))
      N=4*NITER*RCW
      allocate(WRK_WORK1(N))
      N=6*NITER
      allocate(WRK_RWORK1(N))

      IF(CARITH) THEN
        CALL ZBANDDIA(ITAP5,NITER,NP,NP+1,NEIGENV,WRK_T,WRK_EVL,
     &           WRK_ZEV,WRK_WORK1,WRK_RWORK1)
      ELSE
        CALL RBANDDIA(ITAP5,NITER,NP,NP+1,NEIGENV,WRK_T,WRK_EVL,
     &           WRK_ZEV,WRK_RWORK1)
      ENDIF
c
c remove spurious/multiple eigenvalues and
c analyze eigenvectors of band matrix
c open corresponding configuration file
c
      OPEN(ITAP6,FILE=FILECNF,
     &     ACCESS='DIRECT',RECL=IRECL,STATUS='UNKNOWN')

      N=2*NP*NITER*RCW
      allocate(WRK_Z(N))
      N=NITER
      allocate(WRK_IND(N))
      allocate(WRK_MLT(N))
      allocate(WRK_LSP(N))

      IF(CARITH) THEN
         CALL ZBLANA(IW,KREP,niter,np,nocc,ionizlevel,
     &        wrk_evl,wrk_zev,wrk_z,
     &        wrk_t,wrk_ind,wrk_mlt,wrk_lsp,
     &        ITAP5,ITAP6,EVLLIM,NMSPEC)
      ELSE
         CALL RBLANA(IW,KREP,niter,np,nocc,ionizlevel,
     &        wrk_evl,wrk_zev,wrk_z,
     &        wrk_t,wrk_ind,wrk_mlt,wrk_lsp,
     &        ITAP5,ITAP6,EVLLIM,NMSPEC)
      ENDIF

      deallocate(WRK_IND)
      deallocate(WRK_MLT)
      deallocate(WRK_LSP)
      deallocate(WRK_Z)
c
c  close configuration file
c
      CLOSE(ITAP6)
c
c check if user wants to calculate eigenvectors of full ADC matrix
c NGT is the actual number of nonspurious/unique eigenvalues.
c *** important: the corresponding intensities and corrections for
c the multiple eigenvalues are accounted for in the short eigenvector
c array z. This array has to be transferred to the eigenvector
c calculation.
C The corresponding indices sit in the array IND2 which has to be
C transferred to the eigenvector routine.
c
c
      IF(NEIGENV.GT.0) THEN

        LFLAG = .TRUE.

        N=LADC*NEIGENV*RCW
        allocate(WRK_AEV(N))
c
c  next redo Lanczos now with the FULL eigenvectors of T and create
c  the LONG eigenvectors for the TOTAL ADC matrix
c  this is executed if asked for long eigenvectors.
c  file 1 and 2 (the ADC matrix/diag files) are still open.
c

!       IF(DOINCORE.EQV..true.) THEN
!         IF(CARITH.EQV..TRUE.) THEN
!           CALL ZLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
!    &                       CC(WRK_P), CC(WRK_T), CC(WRK_R),
!    &                       CC(WRK_DIAG), CC(WRK_ZEV), CC(WRK_AEV),
!    &                       CC(WRK_BUF),CC(WRK_OI), CC(WRK_OJ),
!    &                       CC(WRK_EVL),LFLAG, INTBUF,
!    &                       ITAP1, ITAP2, ITAP3, NINCWORDS,
!    &                       cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),
!    &                       ivalid)
!         ELSE
!           CALL RLANCZOS_IC(KREP,LADC, NP, NITER, NEIGENV, CC(WRK_IVX),
!    &                       CC(WRK_P), CC(WRK_T), CC(WRK_R),
!    &                       CC(WRK_DIAG), CC(WRK_ZEV), CC(WRK_AEV),
!    &                       CC(WRK_BUF),CC(WRK_OI), CC(WRK_OJ),
!    &                       CC(WRK_EVL),LFLAG, INTBUF,
!    &                       ITAP1, ITAP2, ITAP3, NINCWORDS,
!    &                       cc(wrk_idr), cc(wrk_idi), cc(wrk_idj),
!    &                       ivalid)
!         ENDIF
!         CALL MEMREL('After incore Lanczos',CC,1,KINCSAV,KFREE,LFREE)
!         WRITE(IW,*)  'Memory after incore Lanczos released.'
!       ELSE
          IF(CARITH.EQV..TRUE.) THEN
            CALL ZLANCZOS(KREP,LADC, NP, NITER, NEIGENV, WRK_IVX,
     &                    WRK_P, WRK_T, WRK_R,
     &                    WRK_DIAG, WRK_ZEV, WRK_AEV,
     &                    WRK_BUF,WRK_OI, WRK_OJ,
     &                    WRK_EVL,LFLAG, INTBUF,
     &                    ITAP1, ITAP2, ITAP3)
          ELSE
            CALL RLANCZOS(KREP,LADC, NP, NITER, NEIGENV, WRK_IVX,
     &                    WRK_P, WRK_T, WRK_R,
     &                    WRK_DIAG, WRK_ZEV, WRK_AEV,
     &                    WRK_BUF,WRK_OI, WRK_OJ,
     &                    WRK_EVL,LFLAG, INTBUF,
     &                    ITAP1, ITAP2, ITAP3)
          ENDIF

!       ENDIF    ! doincore

        deallocate(WRK_AEV)

      ENDIF    !neigenv.gt.0
c
c  close and delete the (diag) and the packed ADC matrix file
c
      CLOSE(ITAP1,STATUS='DELETE')
      IF(IOLDNEW.EQ.1) CLOSE(ITAP2,STATUS='DELETE')

c
c  from here the file handles itap1 - itap3 are available again and we
c  can perform the detailed analysis of the long ADC eigenvectors together
c  with their configurations
c
c  All the memory is released and we do the analysis with a F90 routine.
c  irecl still contains the record length of the configuration file.
c
      itap1=infohandle+1
      itap2=infohandle+2
      itap3=infohandle+3
c
c  call long vector analyzer
c  handles 1-3 are exclusively used there, all files are closed
c  on exit.
c
      IF(NEIGENV.GT.0) THEN
        IF(CARITH) THEN
          CALL LVECANALY_C(itap1,itap2,itap3,krep,ladc,
     &                     nocc,irecl,filecnf,ionizlevel)
        ELSE
          CALL LVECANALY_R(itap1,itap2,itap3,krep,ladc,
     &                     nocc,irecl,filecnf,ionizlevel)
        ENDIF
      ENDIF

      deallocate(WRK_EVL)
      deallocate(WRK_ZEV)
      deallocate(WRK_WORK1)
      deallocate(WRK_RWORK1)

      deallocate(WRK_IVX)
      deallocate(WRK_P)
      deallocate(WRK_T)
      deallocate(WRK_R)
      deallocate(WRK_DIAG)
      deallocate(WRK_BUF)
      deallocate(WRK_OI)
      deallocate(WRK_OJ)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RLANCZOS(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                    P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                    CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN)
C
      IMPLICIT INTEGER (A-Z)
 
      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      REAL*8 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      REAL*8 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
      REAL*8 CPUTIM1,CPUTIM2,timer1,timer2,timertot
C
C---------------Description--------------------------------------------
C
C  Perform the real Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is asked for we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      REAL*8 Y,DDOT,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
c
c  create symmetry-specific name for the long eigenvector file.
c  The eigenvector file of the Lanczos matrix needed for the
c  creation of the long eigenvectors is generically stored.
c
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering real Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering real Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A30,I7)') ' Symmetry:',IRREP
      WRITE(IW,'(A30,I7)') ' Matrix dimension:',LADC
      WRITE(IW,'(A30,I7)') ' Band width:',NP
      WRITE(IW,'(A30,I7)') ' Number of iterations:',NITER
      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
c       WRITE(IW,'(A30,A)') ' Long eigenvectors in:',LEVECFN
      ENDIF
      WRITE(IW,*)
c
c  ... clear large eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=0.0D0
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TMATEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF

c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=0.0D0
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize cartesian start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=0.0D0
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=0.0D0
      ENDDO
      DO I=1,NP
        P(I,I+NP)=1.0D0
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO
c
c  ... and time count
c
      timertot=0.0d0
      CALL CPUUSED(CPUTIM1)


c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
        CALL CPUUSED(CPUTIM2)
        timer1 = cputim2 - cputim1
        write(*,*) 'Time per 100 iterations:',timer1
        timertot = timertot + timer1
        cputim1 = cputim2
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
c
c  since for the double ionization the diagonal elements are
c  not stored separately, we define an alternative matrix-vector
c  multiplication routine which directly multiplies the matrix.
c  If IO_DIA has value of -255 (which is surely no valid file handle)
c  we know that the alternative routine should be accessed.
c
c  att: we need the buffer number in RMATBLK_ND! In order to avoid a lot
c  of parameter rewriting we communicate this via a common block.
c
        IF(IO_DIA.eq.-255) THEN
          CALL RMATBLK_ND(LADC, P(1,IVX(NP+1)), R, BUF,
     &                    OI, OJ, INTBUF, IO_OFF)
        ELSE
          CALL RMATBLK(LADC, P(1,IVX(NP+1)), R, DIAG, BUF,
     &                 OI, OJ, INTBUF, IO_DIA, IO_OFF)
        ENDIF

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-T(ITROW,ITCOL)*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = DDOT(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = DDOT(LADC,R,1,R,1)
        Y = DSQRT(Y)
        T(NP+1,IL) = Y

        Y = 1.0d0/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED', 
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      write(*,*) 'Total time spent in Matrix * Vector:',timertot
      write(*,*) 'Average time for Matrix * Vector:',
     &            timertot/dble(niter)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZLANCZOS(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                    P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                    CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      COMPLEX*16 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      COMPLEX*16 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
C
C---------------Description--------------------------------------------
C
C  Perform the complex Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is demanded we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      COMPLEX*16 Y,ZDOTC,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
C
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering complex Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering complex Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A30,I7)') ' Symmetry:',IRREP
      WRITE(IW,'(A30,I7)') ' Matrix dimension:',LADC
      WRITE(IW,'(A30,I7)') ' Band width:',NP
      WRITE(IW,'(A30,I7)') ' Number of iterations:',NITER
      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
c       WRITE(IW,'(A30,A)') ' Long eigenvectors in:',LEVECFN
      ENDIF
      WRITE(IW,*)
c
c  ... clear approx eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILe='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF
c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=(0.0D0,0.0D0)
      ENDDO
      DO I=1,NP
        P(I,I+NP)=(1.0D0,0.0D0)
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
        IF(IO_DIA.eq.-255) THEN
          CALL CMATBLK_ND(LADC, P(1,IVX(NP+1)), R, BUF,
     &                 OI, OJ, INTBUF, IO_OFF)
        ELSE
          CALL CMATBLK(LADC, P(1,IVX(NP+1)), R, DIAG, BUF,
     &                 OI, OJ, INTBUF, IO_DIA, IO_OFF)
        ENDIF


c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-DCONJG(T(ITROW,ITCOL))*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = ZDOTC(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = ZDOTC(LADC,R,1,R,1)
c
        Y = SQRT(Y)
        T(NP+1,IL) = Y

        Y = (1.0d0,0.0d0)/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED',
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RLANCZOS_IC(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                       P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                       CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN,
     &                       NINCW,BUFC,OIC,OJC,IVALID)
C
      IMPLICIT INTEGER (A-Z)
 
      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      REAL*8 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      REAL*8 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
      INTEGER NINCW
      REAL*8 BUFC(NINCW)
      INTEGER OIC(NINCW),OJC(NINCW),IVALID
      REAL*8 CPUTIM1,CPUTIM2,timer1,timer2,timertot
C
C---------------Description--------------------------------------------
C
C  Perform the real Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is asked for we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      REAL*8 Y,DDOT,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
c
c  create symmetry-specific name for the long eigenvector file.
c  The eigenvector file of the Lanczos matrix needed for the
c  creation of the long eigenvectors is generically stored.
c
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering real IC-Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering real IC-Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A40,I12)') ' Symmetry:',IRREP
      WRITE(IW,'(A40,I12)') ' Matrix dimension:',LADC
      WRITE(IW,'(A40,I12)') ' Band width:',NP
      WRITE(IW,'(A40,I12)') ' Number of iterations:',NITER
      WRITE(IW,'(A40,I12)') ' Number of words for incore:',NINCW
      WRITE(IW,'(A40,I12)') ' Number of valid elements incore:',IVALID
      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
c       WRITE(IW,'(A30,A)') ' Long eigenvectors in:',LEVECFN
      ENDIF
      WRITE(IW,*)

      IF(IO_DIA.ne.-255)
     &  CALL QUIT('No incore version for old-format matrices!')
c
c  ... clear large eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=0.0D0
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TMATEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF

c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=0.0D0
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize cartesian start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=0.0D0
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=0.0D0
      ENDDO
      DO I=1,NP
        P(I,I+NP)=1.0D0
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO
c
c  ... and time count
c
      timertot=0.0d0
      CALL CPUUSED(CPUTIM1)

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
        CALL CPUUSED(CPUTIM2)
        timer1 = cputim2 - cputim1
        write(*,*) 'Time per 100 iterations:',timer1
        timertot = timertot + timer1
        cputim1 = cputim2
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
c
c  since for the double ionization the diagonal elements are
c  not stored separately, we define an alternative matrix-vector
c  multiplication routine which directly multiplies the matrix.
c  If IO_DIA has value of -255 (which is surely no valid file handle)
c  we know that the alternative routine should be accessed.
c
c  att: we need the buffer number in RMATBLK_ND! In order to avoid a lot
c  of parameter rewriting we communicate this via a common block.
c  We transfer the huge memory chunk to the incore Mat/Vec routine.
c
        CALL RMATBLK_IC(LADC, P(1,IVX(NP+1)), R, BUF,
     &                  OI, OJ, INTBUF, IO_OFF,
     &                  NINCW,BUFC,OIC,OJC,IVALID)

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-T(ITROW,ITCOL)*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = DDOT(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = DDOT(LADC,R,1,R,1)
        Y = DSQRT(Y)
        T(NP+1,IL) = Y

        Y = 1.0d0/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED', 
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      write(*,*) 'Total time spent in Matrix * Vector:',timertot
      write(*,*) 'Average time for Matrix * Vector:',
     &            timertot/dble(niter)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZLANCZOS_IC(IRREP, LADC, NP, NITER, NEIGENV, IVX,
     &                       P, T, R, DIAG, Z, AEV, BUF, OI, OJ, CEVL, 
     &                       CFLAG, INTBUF, IO_OFF, IO_DIA, IO_LAN,
     &                       NINCW,BUFC,OIC,OJC,IVALID)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER IRREP,LADC,NP,NITER,NEIGENV,IVX(2*NP+1)
      INTEGER INTBUF,IO_OFF,IO_DIA,IO_LAN
      COMPLEX*16 P(LADC,2*NP+1), T(NP+1,NITER), R(LADC), DIAG(LADC)
      COMPLEX*16 Z(NITER,NITER), AEV(LADC,NEIGENV), BUF(INTBUF)
      REAL*8 CEVL(NITER)
      INTEGER OI(INTBUF),OJ(INTBUF)
      LOGICAL CFLAG
      INTEGER NINCW
      COMPLEX*16 BUFC(NINCW)
      INTEGER OIC(NINCW),OJC(NINCW),IVALID
      REAL*8 CPUTIM1,CPUTIM2,timer1,timer2,timertot

C
C---------------Description--------------------------------------------
C
C  Perform the complex Lanczos iteration, generate matrix T_ij
C  the arrays BUF,OI,OJ are buffers for the adc matrix reading process
C  the file handles are IO_DIA and IO_OFF
C  If only eigenvalues are to be calculated, resp. in the first call of
C  this routine the Z-array is empty and will not be referenced !
C  In an eigenvector run Z contains the eigenvectors of the band matrix
C  diagonalization and is multiplied with the Lanczos vectors in order
C  to obtain the approximative eigenvectors of the original problem.
C  the variable .flag. determines if the additional eigenvector run
C  should be performed which is equivalent to the presence of the Z
C  vectors.
C  We write out the Lanczos matrix to file because it will be needed
C  for possible eigenvector calculations.
C  If an eigenvalue calculation is demanded we have to read in the
C  matrix of purged eigenvectors from ZMAT.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables and function definitions ------------
C
      COMPLEX*16 Y,ZDOTC,Z1
      CHARACTER*10 LEVECFN
      real*8 rdummy
C
C---------------Executable code--------------------------------------
C
      IF(IRREP.GT.9) THEN
        WRITE(LEVECFN,'(A8,I2)') 'LONGEVC.',IRREP
      ELSE
        WRITE(LEVECFN,'(A9,I1)') 'LONGEVC.0',IRREP
      ENDIF

      IF(CFLAG) THEN
        CALL PST('Entering complex IC Lanczos for Eigenvectors+')
      ELSE
        CALL PST('Entering complex IC Lanczos for Eigenvalues+')
      ENDIF

      WRITE(IW,'(A30,I7)') ' Symmetry:',IRREP
      WRITE(IW,'(A30,I7)') ' Matrix dimension:',LADC
      WRITE(IW,'(A30,I7)') ' Band width:',NP
      WRITE(IW,'(A30,I7)') ' Number of iterations:',NITER
      WRITE(IW,'(A40,I12)') ' Number of words for incore:',NINCW

      IF(CFLAG) THEN
        WRITE(IW,'(A30,I7)') ' Requested eigenvectors:',NEIGENV
      ENDIF

      IF(IO_DIA.ne.-255)
     &  CALL QUIT('No incore version for old-format matrices!')
c
c  ... clear approx eigenvector array if requested
c
      IF(CFLAG) THEN
        DO I=1,NEIGENV
          DO J=1,LADC
            AEV(J,I)=(0.0D0,0.0D0)
          ENDDO
        ENDDO
      ENDIF
c
c  ... read in purged eigenvector matrix TMATEVC for calculation of
c      the long eigenvectors. These are folded into the Lanczos
c      iterations. The third available file handle is used for that.
c      After the read process the handle can be reused and file is
c      deleted because this is the last position in the code where
c      it is needed.
c
      IF(CFLAG) THEN
        OPEN(IO_LAN,FILe='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        REWIND(IO_LAN)
        READ(IO_LAN) N1
        READ(IO_LAN) N2
        WRITE(IW,*) 'Found',N1,' Lanczos eigenvectors of length',
     &              N2,' on TMATEVC'
        IF(NEIGENV.gt.N1) THEN
          WRITE(*,*) 'More EVCs requested than available.'
          WRITE(*,*) 'Adjusting to the value given above.'
          NEIGENV = N1
        ENDIF
        IF(N2.NE.NITER)
     &    CALL QUIT('Incompatible T matrix eigenvector lengths!')
        DO I=1,NEIGENV
          READ(IO_LAN) CEVL(I)
          READ(IO_LAN) (Z(IX,I),IX=1,N2)
        ENDDO
        CLOSE(IO_LAN,STATUS='DELETE')
      ENDIF
c
c  ... clear Lanczos matrix
c
      DO I=1,NITER
        DO J=1,NP+1
          T(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
c
c  ... clear r/p-vector array and initialize start vectors
c      start vectors are orthonormal and sit in the fields
c      np+1,np+2,...,2*np
c
      DO I=1,2*NP+1
        DO J=1,LADC
          P(J,I)=(0.0D0,0.0D0)
        ENDDO
      ENDDO
      DO J=1,LADC
        R(J)=(0.0D0,0.0D0)
      ENDDO
      DO I=1,NP
        P(I,I+NP)=(1.0D0,0.0D0)
      ENDDO

c
c  ... initialize vector mapping array
c
      DO I=1,2*NP+1
        IVX(I)=I
      ENDDO
c
c  ... and time count
c
      timertot=0.0d0
      CALL CPUUSED(CPUTIM1)

c----------------------------- lanczos loop start ----------------------

      DO 99 IL = 1,NITER
c
c
      IF(mod(il,100).eq.0) THEN
        WRITE(IW,*) '---  Band Lanczos iteration ',IL,' of',niter
        CALL CPUUSED(CPUTIM2)
        timer1 = cputim2 - cputim1
        write(*,*) 'Time per 100 iterations:',timer1
        timertot = timertot + timer1
        cputim1 = cputim2
      ENDIF
c
c  form r = Hq_i - sum t_ji*q_j
c  the j index in the sum runs from i-p up to i-1
c  and the loop is only executed if previous vectors exist.
c  the actual vector to be multiplied by A always sits in
c  P(*,np+1)
c
      CALL CMATBLK_IC(LADC, P(1,IVX(NP+1)), R, BUF,
     &                OI, OJ, INTBUF, IO_OFF,
     &                NINCW,BUFC,OIC,OJC,IVALID)

c  subtract all the previous components

        ICC = 1
        DO IP = MAX0(IL-NP,1), IL-1
          ITROW = ICC+1
          ITCOL = IL - ICC
          DO J=1,LADC
            R(J) = R(J)-DCONJG(T(ITROW,ITCOL))*P(J,IVX(NP-ICC+1))
          ENDDO
          ICC = ICC + 1
        ENDDO
c
c  form the t_ji = q_j^+ * r and r=r-t_ji*q_j for all
c  j=i...i+p-1
c
        DO IP = 1,NP
          Y = ZDOTC(LADC,P(1,IVX(NP+IP)),1,R,1)
          T(IP,IL) = Y
          DO J=1,LADC
            R(J) = R(J) - Y*P(J,IVX(NP+IP))
          ENDDO
        ENDDO
c
c  now r contains the whole RHS of the band-lanczos equation
c  and will be normalized
c  form t_i+p,i = |r|, q_i+p,i = r/|r|
c
        Y = ZDOTC(LADC,R,1,R,1)
c
        Y = SQRT(Y)
        T(NP+1,IL) = Y

        Y = (1.0d0,0.0d0)/Y
        DO J=1,LADC
          P(J,IVX(2*NP+1)) = R(J) * Y
        ENDDO
c
c  if applicable form the approximate eigenvector
c  with the first Lanczos vector Q_1
c  Important: the ith ready Lanczos vector is available at
c  P(*,IVX(NP+1)) exactly at this position of the
c  iteration. In case of an eigenvector run its contributions
c  to ALL eigenvectors have to be taken into account. After
c  cycling the contributions of the i+1. lanczos vector will
C  be folded into the final eigenvectors
c
        IF(CFLAG) THEN
          DO IOM = 1,NEIGENV
            Z1 = Z(IL,IOM)
            DO ICOL = 1,LADC
              AEV(ICOL,IOM) = AEV(ICOL,IOM) + Z1*P(ICOL,IVX(NP+1))
            ENDDO
          ENDDO
        ENDIF
c
c  now cycle the mapping array IVX
c  and start new main iteration
c
        CALL IVXROT(IVX,2*NP+1)

 99   CONTINUE
c
c  write out long ADC eigenvectors when eigenvector run finished
c
      IF(CFLAG) THEN
        OPEN(IO_LAN, FILE=LEVECFN, FORM='UNFORMATTED',
     &       ACCESS='SEQUENTIAL', STATUS='UNKNOWN')
        WRITE(IO_LAN) LADC,NEIGENV
        DO J=1,NEIGENV
          WRITE(IO_LAN) CEVL(J)
          WRITE(IO_LAN) (AEV(I,J),I=1,LADC)
        ENDDO
        CLOSE(IO_LAN)
        WRITE(IW,*) ' ADC LONG EIGENVECTORS WRITTEN.'
      ENDIF

      write(*,*) 'Total time spent in Matrix * Vector:',timertot
      write(*,*) 'Average time for Matrix * Vector:',
     &            timertot/dble(niter)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBLANA(IW,IRREP,N,NP,NHOLES,IONIZLEVEL,
     &                  D,ZLONG,Z,T,IND,MLT,LSP,
     &                  IOUNIT,IOUNITC,ELIM,NMSPEC)

      IMPLICIT INTEGER (A-Z)

      INTEGER IW,IRREP,N,NP,NHOLES,IONIZLEVEL
      REAL*8 D(N)
      REAL*8 ZLONG(N,N),Z(2*NP,N),T(NP+1,N)
      INTEGER IND(N),MLT(N),LSP(N)
      INTEGER IOUNIT,IOUNITC
      REAL*8 ELIM
      CHARACTER*5 NMSPEC
C
C---------------Description--------------------------------------------
C
C  Analyses the eigenvalues/vectors of the Band matrix T and
C  treats multiple/spurious eigenvalues. All this for the
C  current main program symmetry stored in IRREP
c  (H.-D. Meyer's and Thomas Sommerfeld's formulation)
c
c  IOUNIT is a free handle for creating the XSPEC files
c  IOUNITC is an occupied handle referring to the config file.
c
c  eigenvalues come in D(*), complete eigenvectors in ZLONG(N,N) !
c  the first and last NP components will be needed for the analysis 
c  only. in order to use the unmodified code we copy those components
c  to the short vector array z(2*np,n)
c
c  the number N corresponds to the number of L-iterations == NITER
c
c  *** attention *** in the analysis some components of z are added
c  together or multiplied by some scaling factor stemming from the
c  multiplicity of some eigenvales. These modifications should also be applied
c  to the LONG eigenvector as well in order to get the right results.
c  If approx. eigenvectors are calculated the LONG eigenvector array
c  will be used because we need EVERY component of Z not only the
c  first and last NP ones. This causes trouble because according to
c  Dieter's analysis only the hole parts (the first NP ones) are weighted
c  and for the other components no strict modification scheme
c  is given. This means that we should determine the corrected long
c  Z vector by INVERSE INTERATION. There no artefacts will occur because
c  multiple and spurious eigenvalues are artefacts from LANCZOS and not
c  from the physical problem.
c
c
c  the band-matrix T is needed for the error analysis
c  according to Parlett: "The symmetric eigenvalue problem",
c  Chapter 13, sect 13-10-2 where a formula for the accuracy
c  of a Rayleigh-Ritz pair is given for the case of band-Lanczos.
C  From that one sees that the unused lower right triangle of T
C  together with the last NP components of the calculated eigenvectors
c  are needed in order to estimate the quality and convergence of 
c  the corresponding eigenvector.
c
c  the variable cntev actually counts the printed eigenvalues.
c  This is of use because each eigenvalue sequence then gets a contiguous
c  numbering. Numbering according to threshold will have gaps.
C
      REAL*8 AA,AK,A1K,SUM,SUP
      REAL*8 AUEV,THRS,OVL
      REAL*8 X1,Y1

      REAL*8 XM,DDOT
      DATA AUEV/27.2113957D0/

      CHARACTER*8 SPECFILE
      REAL*8 THTHR,WTHR
      CHARACTER*44 FIELD
C
C---------------Executable code--------------------------------------
C


      CALL PST('Processing the spectrum (real)+')
      WTHR=0.01D0
      WRITE(IW,*) 'Upper limit for printing eigenvalues on screen:',ELIM
      WRITE(IW,*) 'Threshold for considering eigenvectors:',WTHR
      WRITE(IW,*) 'Number of iterations:',N
      WRITE(IW,*) 'Band width:',NP

c
c  copy the first and last NP components of the long eigenvectors
c  to Z in order to apply dieters code easily.
c
      do i=1,n
        IND(i)=0
        MLT(i)=0
        LSP(i)=0
      enddo

c  prepare short array z
         
      do i=1,n
      do j=1,np
        z(j,i)=zlong(j,i)
        z(j+np,i)=zlong(n-np+j,i)
      enddo
      enddo

C     ----------------------------------------------------------
C     DETERMINE THE INDEX-VECTOR 'IND' IN ORDER TO REMOVE
C     MULTIPLE AND SPURIOUS EIGENVALUES.
C     ----------------------------------------------------------

      NP1 = NP + 1
      NP2 = 2 * NP

      NSP = 0
      NMU = 0
      ISP = 0
      JND = 0
      J   = 0
C
  50  J   = J + 1
      IF( J .GT. n )                               GOTO  40

      AA=DABS(DDOT(NP,Z(1,J),1,Z(1,J),1))

      NSP = NSP + 1
      ISP = ISP + 1
C.....TEST OF SPURIOUS EIGENVALUE.
      IF( AA .LT. 1.D-14 )                     GOTO  50
      NSP = NSP - 1
      ISP = ISP - 1
      JND = JND + 1
      IND(JND) = J
      MLT(JND) = 0
      LSP(JND) = ISP
      ISP      = 0
C
      IF( J .EQ. n )                               GOTO  40
      THRS = 5.D-7*( 1.D0+DABS(d(J)) )
C.....TEST OF MULTIPLE EIGENVALUE.
      IF( DABS( d(J) - d(J+1) ) .GT. THRS )           GOTO  50
      SUM = AA
      SUP = AA
      KK  = 0
      IA  = 0
C
      DO  41  K = 1, n - J
      OVL = 1.D-5
      IF( DABS( d(J+IA) - d(J+K) ) .GT. THRS )        GOTO  42

      AK=DABS(DDOT(NP,Z(1,J+K),1,Z(1,J+K),1))

      IF( AK .LE. 1.D-14 )  THEN
         KK     = K
         GOTO  41
      ELSE IF( AK/SUP + SUP/AK  .LT.  2.D+6 )    THEN
           OVL = 1.D-7
           IF( DABS(d(J+IA)-d(J+K)) .GT. 1.D-2*THRS )   GOTO  42
      END IF

      A1K=dabs(ddot(NP,Z(1,J+K),1,Z(1,J+IA),1))

      A1K = A1K*A1K/(SUP*AK)
      IF( DABS( A1K-1.D0 ) .LE. OVL )            THEN
         SUM = SUM + AK
         IF( AK .GT. SUP )  THEN
           SUP = AK
           IA  =  K
         END IF
         KK = K
      ELSE
        GOTO  42
      END IF
  41  CONTINUE
C
  42  d(J) = d(J+IA)
      DO  43  I = 1, NP
         Z(   I,J) = Z(   I,J+IA)*DSQRT(SUM/SUP)
         Z(NP+I,J) = Z(NP+I,J+IA)
  43  CONTINUE
C
      NMU = NMU + KK
      J   =  J  + KK
      MLT(JND) =  KK
      GOTO  50
C
  40  CONTINUE
C
C.....REMOVE SPURIOUS AND MULTIPLE EIGENVALUES.
C  ATT!
C  The complete eigenvectors also have to be purged!
C  In the Z-array only the upper and lower NP components
C  are available!
C
      DO  48  J = 1, JND
         I = IND(J)
         d(J) = d(I)
         DO  49  K = 1, NP2
            Z(K,J) = Z(K,I)
 49      continue
         DO  67  K = 1, N
            ZLONG(K,J) = ZLONG(K,I)
 67      continue
 48   CONTINUE

      write (IW,'(/,A)') '  Information on purging:'
      PRINT 905, JND, N, int(100.0d0*dble(jnd)/dble(n))
 905  FORMAT(4x,'There are',I5,' relevant residues out of',I5, 
     $     '   (',I3,'%)')
      write(6,'(4x,A,I4,6x,A,I4)') 'N_spurios = ', NSP,
     $     'N_multiple = ', NMU
c     write(6,*) 'Threshold for vector printing:',WTHR
C
C.....TRANSFORM THE EIGENVECTORS(FIRST NP ROWS) TO OBTAIN THE RESIDUES.
C.....STORE THE RESULT IN THE FIRST NP ROWS OF THE EIGENVECTOR MATRIX Z.
C.....MULTIPLY :  EIGENVECTORS(LAST ROWS) * 'ERROR-TRIANGLE'
C.....STORE THIS RESULT IN THE LAST ROWS OF THE EIGENVECTOR MATRIX Z.
C
      DO  30  J  = 1, JND
       DO  35  I  = 1, NP
         XM = 0.0D0
         DO  36  K = I, NP
           KK = K + NP
           IC = K + n - NP
           JC = I - K + NP1
           XM = XM + Z(KK,J)*T(JC,IC)
  36     CONTINUE
  35   Z(I+NP,J) = XM
  30  CONTINUE     
c
c    * * * * * *     S P E C T R U M     * * * * * * * * * * * *
c  spectrum is purged and available: all eigenvalues with corresponding
c  intensities obtained as squares of the hole/hole components
c  (pole strengths). IMPORTANT: the band width of the Lanczos matrix
c  can be larger than the number of actual hole states due to stability
c  reasons. If we sum the intensities we just look at the hole states !
c
c  print eigenvalues and accuracy in a more convenient format
c  threshold applies only to the screen output. pole strength = square of
c  eigenvector component. Write all eigenvalues and pole strength
c  to file named in SPECFILE
c
      if(irrep.gt.9) then
        write(SPECFILE,'(a5,a1,i2)') NMSPEC,'.',irrep
      else
        write(SPECFILE,'(a5,a2,i1)') NMSPEC,'.0',irrep
      endif
      OPEN(IOUNIT,FILE=SPECFILE,FORM='FORMATTED',STATUS='UNKNOWN')

      cntev=0
      DO J = 1, JND 
c
c  pole strength calculation
c
        x1 = 0.0d0
        do i = 1, nholes
           x1 = x1 + z(i,j)*z(i,j)
        enddo
c
c  error analysis
c
        y1 = 0.0d0
        do i = 1, np
          y1 = y1 + z(np+i,j)*z(np+i,j)
        enddo
        y1 = dsqrt(y1)
c
c  everything goes to file, no threshold for pole strength and
c  no upper limit for the eigenvalues. That means: all eigenvalues
c  together with their pole strengths and error estimates are written
c  to the XSPEC.XX files.
c
        IF(IONIZLEVEL.EQ.1) THEN
c
c  take care of single ionizations
c
           isipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) isipcount = isipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)')
     &        d(j)*AUEV,x1,y1,'@',isipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,F12.5,1X,"(",F7.5,")",A8,A22)')
     &         iths,':  ',z(iths,j),z(iths,j)*z(iths,j),'   conf:',field
             endif
           ENDDO

        ELSE
c
c  take care of double ionizations
c
           idipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) idipcount = idipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)') 
     &        d(j)*AUEV,x1,y1,'@',idipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=z(iths,j)*z(iths,j)
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,F12.5,1X,"(",F7.5,")",A8,A44)') 
     &         iths,':  ',z(iths,j),z(iths,j)*z(iths,j),'  conf: ',field
             endif
           ENDDO          
        ENDIF
c
c  print on screen only those eigenvalues that are below the given limit.
c  use cntev for contiguous numbering !
c
c  ATT !!! For the two-hole final states an output of the squared
c  coefficient is hiding information! We therefore write out the
c  coefficients themselves!
c
        if(d(j)*AUEV.le.ELIM) THEN
          cntev = cntev + 1
          write(iw,388)
 388      format(//2x,60('-')/'  no.    eigenvalue (eV)   mean-res.',
     &    3x,'error estimate'/2x,60('-')/)
          write(iw,389) cntev,'@', d(j)*AUEV, x1, y1, ind(j)
 389      format(I5,A1,2x,F14.8,4x,F8.4,4x,1PE11.2,'  (',I5,')')
          IF(IONIZLEVEL.EQ.2) THEN
            write(iw,*) '       -----Contributing 2h configurations:'
            DO iths=1,nholes
              ththr=z(iths,j)*z(iths,j)
              if(ththr.gt.0.0001) then
                read(iounitc,rec=iths) field
                write(iw,'(12X,I5,A3,F12.5,A22)') 
     &             iths,':  ',z(iths,j),field
              endif
            ENDDO          
          ELSE
            write(iw,'(100I8)') (ix,ix=1,nholes)
            write(iw,'(100F8.4)') 
     &           (z(ix,j)*z(ix,j),ix=1,nholes)
          ENDIF
        endif
      ENDDO

      CLOSE(IOUNIT)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBLANA(IW,IRREP,N,NP,NHOLES,IONIZLEVEL,
     &                  D,ZLONG,Z,T,IND,MLT,LSP,
     &                  IOUNIT,IOUNITC,ELIM,NMSPEC)

      IMPLICIT INTEGER (A-Z)

      INTEGER IW,IRREP,N,NP,NHOLES,IONIZLEVEL
      REAL*8 D(N)
      COMPLEX*16 ZLONG(N,N),Z(2*NP,N),T(NP+1,N)
      INTEGER IND(N),MLT(N),LSP(N)
      INTEGER IOUNIT,IOUNITC
      REAL*8 ELIM
      CHARACTER*5 NMSPEC
C
C---------------Description--------------------------------------------
C
C  Analyses the eigenvalues/vectors of the Band matrix T and
C  treats multiple/spurious eigenvalues. All this for the
C  current main program symmetry stored in IRREP
c  (H.-D. Meyer's and Thomas Sommerfeld's formulation)
c
c  IOUNIT is a free handle for creating the XSPEC files
c  IOUNITC is an occupied handle referring to the config file.
c
c  eigenvalues come in D(*), complete eigenvectors in ZLONG(N,N) !
c  the first and last NP components will be needed for the analysis 
c  only. in order to use the unmodified code we copy those components
c  to the short vector array z(2*np,n)
c
c  the number N corresponds to the number of L-iterations == NITER
c
c  *** attention *** in the analysis some components of z are added
c  together or multiplied by some scaling factor stemming from the
c  multiplicity of some eigenvales. These modifications should also be applied
c  to the LONG eigenvector as well in order to get the right results.
c  If approx. eigenvectors are calculated the LONG eigenvector array
c  will be used because we need EVERY component of Z not only the
c  first and last NP ones. This causes trouble because according to
c  Dieter's analysis only the hole parts (the first NP ones) are weighted
c  and for the other components no strict modification scheme
c  is given. This means that we should determine the corrected long
c  Z vector by INVERSE INTERATION. There no artefacts will occur because
c  multiple and spurious eigenvalues are artefacts from LANCZOS and not
c  from the physical problem.
c
c
c  the band-matrix T is needed for the error analysis
c  according to Parlett: "The symmetric eigenvalue problem",
c  Chapter 13, sect 13-10-2 where a formula for the accuracy
c  of a Rayleigh-Ritz pair is given for the case of band-Lanczos.
C  From that one sees that the unused lower right triangle of T
C  together with the last NP components of the calculated eigenvectors
c  are needed in order to estimate the quality and convergence of 
c  the corresponding eigenvector.
c
c  the variable cntev actually counts the printed eigenvalues.
c  This is of use because each eigenvalue sequence then gets a contiguous
c  numbering. Numbering according to threshold will have gaps.
C
C---------------Local variables--------------------------------------
C
      REAL*8 AA,AK,A1K,SUM,SUP
      REAL*8 AUEV,THRS,OVL
      REAL*8 X1,Y1

      COMPLEX*16 XM,ZDOTC
      DATA AUEV/27.2113957D0/

      CHARACTER*8 SPECFILE
c variables for hole-configuration analysis
      REAL*8 THTHR,WTHR
      CHARACTER*44 FIELD
C
C---------------Executable code--------------------------------------
C


      CALL PST('Processing the spectrum (complex)+')
      WTHR=0.01D0
      WRITE(IW,*) 'Upper limit for printing eigenvalues on screen',ELIM
      WRITE(IW,*) 'Threshold pole strength for eigenvectors:',WTHR
      WRITE(IW,*) 'Number of iterations:',N
      WRITE(IW,*) 'Band width:',NP
c
c  copy the first and last NP components of the long eigenvectors
c  to Z in order to apply dieters code easily.
c
      do i=1,n
        IND(i)=0
        MLT(i)=0
        LSP(i)=0
      enddo

c  prepare short array z
         
      do i=1,n
      do j=1,np
        z(j,i)=zlong(j,i)
        z(j+np,i)=zlong(n-np+j,i)
!       write(*,*) 'prepare:',i,z(j,i),z(j+np,i)
      enddo
      enddo

C     ----------------------------------------------------------
C     DETERMINE THE INDEX-VECTOR 'IND' IN ORDER TO REMOVE
C     MULTIPLE AND SPURIOUS EIGENVALUES.
C     ----------------------------------------------------------

      NP1 = NP + 1
      NP2 = 2 * NP

      NSP = 0
      NMU = 0
      ISP = 0
      JND = 0
      J   = 0
C
  50  J   = J + 1
      IF( J .GT. n )                               GOTO  40

      AA=ABS(ZDOTC(NP,Z(1,J),1,Z(1,J),1))

      NSP = NSP + 1
      ISP = ISP + 1
C.....TEST OF SPURIOUS EIGENVALUE.
      IF( AA .LT. 1.D-14 )                     GOTO  50
      NSP = NSP - 1
      ISP = ISP - 1
      JND = JND + 1
      IND(JND) = J
      MLT(JND) = 0
      LSP(JND) = ISP
      ISP      = 0
C
      IF( J .EQ. n )                               GOTO  40
      THRS = 5.D-7*( 1.D0+DABS(d(J)) )
C.....TEST OF MULTIPLE EIGENVALUE.
      IF( DABS( d(J) - d(J+1) ) .GT. THRS )           GOTO  50
      SUM = AA
      SUP = AA
      KK  = 0
      IA  = 0
C
      DO  41  K = 1, n - J
      OVL = 1.D-5
      IF( DABS( d(J+IA) - d(J+K) ) .GT. THRS )        GOTO  42

      AK=ABS(ZDOTC(NP,Z(1,J+K),1,Z(1,J+K),1))

      IF( AK .LE. 1.D-14 )  THEN
         KK     = K
         GOTO  41
      ELSE IF( AK/SUP + SUP/AK  .LT.  2.D+6 )    THEN
           OVL = 1.D-7
           IF( DABS(d(J+IA)-d(J+K)) .GT. 1.D-2*THRS )   GOTO  42
      END IF

      A1K=abs(zdotc(NP,Z(1,J+K),1,Z(1,J+IA),1))

      A1K = A1K*A1K/(SUP*AK)
      IF( DABS( A1K-1.D0 ) .LE. OVL )            THEN
         SUM = SUM + AK
         IF( AK .GT. SUP )  THEN
           SUP = AK
           IA  =  K
         END IF
         KK = K
      ELSE
        GOTO  42
      END IF
  41  CONTINUE
C
  42  d(J) = d(J+IA)
      DO  43  I = 1, NP
         Z(   I,J) = Z(   I,J+IA)*DCMPLX(DSQRT(SUM/SUP),0.0d0)
         Z(NP+I,J) = Z(NP+I,J+IA)
  43  CONTINUE
C
      NMU = NMU + KK
      J   =  J  + KK
      MLT(JND) =  KK
      GOTO  50
C
  40  CONTINUE
C
C.....REMOVE SPURIOUS AND MULTIPLE EIGENVALUES.
C  ATT!
C  The complete eigenvectors also have to be purged!
C  In the Z-array only the upper and lower NP components
C  are available!
C
      DO  48  J = 1, JND
         I = IND(J)
         d(J) = d(I)
         DO  49  K = 1, NP2
            Z(K,J) = Z(K,I)
 49      CONTINUE
         DO  67  K = 1, N
            ZLONG(K,J) = ZLONG(K,I)
 67      CONTINUE
 48   CONTINUE

      write (IW,'(/,A)') '  Information on purging:'
      PRINT 905, JND, N, int(100.0d0*dble(jnd)/dble(n))
 905  FORMAT(4x,'There are',I5,' relevant residues out of',I5, 
     $     '   (',I3,'%)')
      write(IW,'(4x,A,I4,6x,A,I4)') 'N_spurios = ', NSP,
     $     'N_multiple = ', NMU
c     write(IW,*) 'Threshold for vector printing:',WTHR
C
C
C.....TRANSFORM THE EIGENVECTORS(FIRST NP ROWS) TO OBTAIN THE RESIDUES.
C.....STORE THE RESULT IN THE FIRST NP ROWS OF THE EIGENVECTOR MATRIX Z.
C.....MULTIPLY :  EIGENVECTORS(LAST ROWS) * 'ERROR-TRIANGLE'
C.....STORE THIS RESULT IN THE LAST ROWS OF THE EIGENVECTOR MATRIX Z.
C
      DO  30  J  = 1, JND
       DO  35  I  = 1, NP
         XM = (0.0D0,0.0D0)
         DO  36  K = I, NP
           KK = K + NP
           IC = K + n - NP
           JC = I - K + NP1
           XM = XM + Z(KK,J)*T(JC,IC)
  36     CONTINUE
       Z(I+NP,J) = XM
  35   CONTINUE
  30  CONTINUE     

c
c    * * * * * *     S P E C T R U M     * * * * * * * * * * * *
c  spectrum is purged and available: all eigenvalues with corresponding
c  intensities obtained as squares of the hole/hole components
c  (pole strengths). IMPORTANT: the band width of the Lanczos matrix
c  can be larger than the number of actual hole states due to stability
c  reasons. If we sum the intensities we just look at the hole states !
c
c
c  print eigenvalues and accuracy in a more convenient format
c  threshold applies only to the screen output. pole strength = square of
c  eigenvector component. Write all eigenvalues and pole strength
c  to file named in SPECFILE
c
      if(irrep.gt.9) then
        write(SPECFILE,'(a5,a1,i2)') NMSPEC,'.',irrep
      else
        write(SPECFILE,'(a5,a2,i1)') NMSPEC,'.0',irrep
      endif
      open(iounit,file=SPECFILE,form='FORMATTED',status='UNKNOWN')
      
      cntev=0
      DO J = 1, JND 
c
c  pole strength calculation
c
        x1 = 0.0d0
        do i = 1, nholes
           x1 = x1 + DBLE( z(i,j)*dconjg(z(i,j)) )
        enddo
c
c  error analysis
c
        y1 = 0.0d0
        do i = 1, np
          y1 = y1 + DBLE( z(np+i,j)*dconjg(z(np+i,j)) )
        enddo
        y1 = dsqrt(y1)
c
c  everything goes to file, no threshold for pole strength and
c  no upper limit for the eigenvalues. That means: all eigenvalues
c  together with their pole strengths and error estimates are written
c  to the XSPEC.XX files.
c
        IF(IONIZLEVEL.EQ.1) THEN
c
c  take care of single ionizations
c
           isipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*dconjg(z(iths,j))
             if(ththr.gt.0.001) isipcount = isipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)')
     &        d(j)*AUEV,x1,y1,'@',isipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=dble(z(iths,j)*dconjg(z(iths,j)))
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,2F12.5,1X,"(",F7.5,")",A8,A22)')
     &         iths,':  ',z(iths,j),dble(z(iths,j)*dconjg(z(iths,j))),
     &         '   conf:',field
             endif
           ENDDO

        ELSE
c
c  take care of double ionizations
c
           idipcount=0
           DO iths=1,nholes
             ththr=z(iths,j)*dconjg(z(iths,j))
             if(ththr.gt.0.001) idipcount = idipcount + 1
           ENDDO 
           write(iounit,'(2x,F14.8,4x,2F12.6,4x,A1,I4,A,I3)')
     &        d(j)*AUEV,x1,y1,'@',idipcount,'  File: ',irrep
           DO iths=1,nholes
             ththr=dble(z(iths,j)*dconjg(z(iths,j)))
             if(ththr.gt.0.001) then
               read(iounitc,rec=iths) field
               write(iounit,'(10X,I5,A3,2F12.5,1X,"(",F7.5,")",A8,A44)')
     &            iths,':  ',z(iths,j),ththr,'  conf: ',field
             endif
           ENDDO
        ENDIF
c
c  print on screen only those eigenvalues that are below the given limit.
c  use cntev for contiguous numbering ! Again: we now write out the
c  coefficients themselves for the two-particle case!
c
        if(d(j)*AUEV.le.ELIM) then
          cntev = cntev + 1
          write(iw,388)
 388      format(//2x,60('-')/'  no.    eigenvalue (eV)   mean-res.',
     &    3x,'error estimate'/2x,60('-')/)
          write(iw,389) cntev,'@', d(j)*AUEV, x1, y1, ind(j)
 389      format(I5,A1,2x,F14.8,4x,F8.4,4x,1PE11.2,'  (',I5,')')

          IF(IONIZLEVEL.EQ.2) THEN
            write(iw,*) '       -----Contributing 2h configurations:'
            DO iths=1,nholes
              ththr=z(iths,j)*dconjg(z(iths,j))
              if(ththr.gt.0.0001) then
                read(iounitc,rec=iths) field
                write(iw,'(12X,I5,A3,2F12.5,A22)') 
     &             iths,':  ',z(iths,j),field
              endif
            ENDDO    
          ELSE
            write(iw,'(100I8)') (ix,ix=1,nholes)
            write(iw,'(100F8.4)') 
     &           (dble(z(ix,j)*dconjg(z(ix,j))),ix=1,nholes)
          ENDIF
        endif

      ENDDO

      CLOSE(IOUNIT)

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RBANDDIA(IOUNIT,N,NP,LDAB,NEIGENV,AB,EVL,Z,RWORK)
C
      use adc_mat
      use memory_allocator
      use adc_fano_diag

      IMPLICIT INTEGER (A-Z)

      INTEGER IOUNIT,N,NP,LDAB
      INTEGER NEIGENV
      REAL*8 AB(LDAB,*),Z(N,*)
      REAL*8 RWORK(6*N),EVL(N)
      
C
C---------------Description--------------------------------------------
C
C  Perform the diagonalization of the real symmetric band matrix
C  by calling the appropriate LAPACK routine. This routine
C  serves as an interface in order to avoid too many command line
C  arguments.
C  The Lanczos matrix T from the caller comes in the AB array.
C     ****  ATTENTION ****  due to the usual programming conventions
C  of the LAPACK routines the original matrix T is DESTROYED !!!
C
C  important: we used to save the eigenvectors of the Lanc matrix
C  *after* purging. But when one analyzes the rebuilding of the long ADC
C  vectors inconsistencies may arise. We therefore store *all* evecs of
C  L and deal with multitudes separately.

C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER INFO,LDZ
      CHARACTER*1 JOBZ,UPLO
      REAL*8 Y
      REAL*8 AUEV
      INTEGER                            :: counter

      INTEGER, ALLOCATABLE, DIMENSION(:) :: positions
      INTEGER                            :: nr_vecs
      INTEGER                            :: j, k
      REAL                               :: upper, lower

C
C---------------Executable code--------------------------------------
C
C we always compute the eigenvectors of the band matrix, because
C they are needed in the further analysis
C These eigenvectors have the FULL length of the Lanczos matrix. 
C but not the length of the original ADC matrix!
C In the original formulation only
C the first and last NP components were available and the whole
C eigenvector was calculated by inverse iteration. We do not pursue this
C direction due to serious problems with the inverse interation
C algorithm
C
C N is the dimension of the band matrix (not the number of sidebands!)

      JOBZ='V'
      UPLO='L'
      LDZ=N

      CALL DSBEV(JOBZ, UPLO, N, NP, AB, LDAB, EVL, Z, LDZ,
     &           RWORK, INFO)

      WRITE(iw,*) 'info = ', info

      IF(INFO.NE.0) THEN
        CALL QUIT('Error in DSBEV')
      ENDIF

      WRITE(IW,*)
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*) 'First 100 RAW Lanczos EVAL/EVECS + intensity.'
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*)
      WRITE(IW,*) ' No.         (a.u.)      (eV)'


      AUEV=27.2113957D0

      DO I=1,MIN(100,N)
        Y = 0.0D0
        DO J=1,NP
          Y = Y + Z(J,I)*Z(J,I)
        ENDDO
        WRITE(IW,111) I,EVL(I),EVL(I)*AUEV,DABS(Y)
      ENDDO

 111  FORMAT(I5,2x,2F12.6,3x,'mean-res',2x,F10.4)
c
c eigenvectors of the Lanczos matrix are ready. Save them *all* to disk
c Since each final state symmetry is treated immediately after
c calculation, TMATVEC does not get an individual symmetry label.
c

      IF (reladc_md_isfano) THEN
        CALL alloc(positions,N,id="array of position numbers fano Lanc")
        CALL fano_find_in_vec(N,EVL,Z,positions,nr_vecs,iw)

        NEIGENV = nr_vecs
        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO counter=1, nr_vecs
          j = positions(counter)
          WRITE(IOUNIT) EVL(J) ! eigenvalue
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        ENDDO
        CLOSE(IOUNIT)
      
        WRITE(*,*) 'Fano set of Lanczos eigenvectors saved to TMATEVC'

        CALL dealloc(positions)

      ELSE
        CALL alloc(positions,N,id='chosen eigenvectors for analysis')
        nr_vecs = 0

        upper = reladc_md_eeigv_upper
        lower = reladc_md_eeigv_lower
        WRITE(IW,*) 'lower = ', lower, ' upper = ', upper

        DO j = 1, N
          IF ((EVL(j).GE.lower).AND.(EVL(j).LE.upper)) THEN
            nr_vecs = nr_vecs + 1
            positions(nr_vecs) = j
          END IF
        END DO

        NEIGENV = nr_vecs

        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors (all!)
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO k = 1, nr_vecs
          j = positions(k)
          WRITE(IOUNIT) EVL(J) ! eigenvalue
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        END DO
        CLOSE(IOUNIT)
    
        WRITE(iw,*) NEIGENV, 'Lanczos eigenvectors saved to TMATEVC'
        WRITE(iw,*) 'Energy range in a.u.: ', lower, '--', upper

        CALL dealloc(positions)
      END IF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE ZBANDDIA(IOUNIT,N,NP,LDAB,NEIGENV,AB,EVL,Z,WORK,RWORK)
C
      use adc_mat
      use memory_allocator
      use adc_fano_diag
!
      IMPLICIT INTEGER (A-Z)

      INTEGER IOUNIT, N, NP, LDAB
      INTEGER NEIGENV
      COMPLEX*16 AB(LDAB,*),Z(N,*),WORK(4*N)
      DOUBLE PRECISION RWORK(6*N),EVL(N)
      
C
C---------------Description--------------------------------------------
C
C  Perform the diagonalization of the Hermitian band matrix
C  by calling the appropriate LAPACK routine. This routine
C  is interfacing in order to avoid too many command line
C  arguments.
C  The Lanczos matrix T from the caller comes in the AB array.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER INFO,LDZ
      CHARACTER*1 JOBZ,UPLO
      COMPLEX*16 Y
      REAL*8 AUEV

      INTEGER                            :: nr_vecs            
      INTEGER                            :: counter
      INTEGER, ALLOCATABLE, DIMENSION(:) :: positions
      INTEGER                            :: j, k
      REAL                               :: upper, lower

C
C---------------Executable code--------------------------------------
C
C we always compute the eigenvectors of the band matrix, because
C they are needed in the further analysis
C These eigenvectors have FULL length. In the original formulation only
C the first and last NP components were available and the whole
C eigenvector was calculated by inverse iteration.
C

      JOBZ='V'
      UPLO='L'
      LDZ=N

      CALL ZHBEV(JOBZ, UPLO, N, NP, AB, LDAB, EVL, Z, LDZ,
     &           WORK, RWORK, INFO)

      IF(INFO.NE.0) THEN
        CALL QUIT('Error in ZHBEV')
      ENDIF

      WRITE(IW,*)
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*) 'First 100 RAW Lanczos EVAL/EVECS + intensity.'
      WRITE(IW,*) '---------------------------------------------------'
      WRITE(IW,*)
      WRITE(IW,*) ' No.         (a.u.)      (eV)'


      AUEV=27.2113957D0

      DO I=1,MIN(100,N)
        Y = (0.0D0,0.0D0)
        DO J=1,NP
          Y = Y + Z(J,I)*DCONJG(Z(J,I))
        ENDDO
        WRITE(IW,111) I,EVL(I),EVL(I)*AUEV,ABS(Y)
      ENDDO

 111  FORMAT(I5,2x,2F12.6,3x,'mean-res',2x,F10.4)
c
c eigenvectors of the Lanczos matrix are ready. Save the chosen to disk
c Since each final state symmetry is treated immediately after
c calculation, TMATVEC does not get an individual symmetry label.
c
                                                                  
      IF (reladc_md_isfano) THEN                                  
        CALL alloc(positions,N,id="array of position numbers fano Lanc")
        CALL cfano_find_in_vec(N,EVL,Z,positions,nr_vecs,iw)    
                                                                  
        NEIGENV = nr_vecs                                      
        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO counter=1, nr_vecs                                  
          j = positions(counter)                                  
          WRITE(IOUNIT) EVL(J) ! eigenvalue                       
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        ENDDO                                                     
        CLOSE(IOUNIT)                                             
                                                                  
        WRITE(*,*) 'Fano set of Lanczos eigenvectors saved to TMATEVC'
                                                                  
        CALL dealloc(positions)                                   
                                                                  
      ELSE                                                        
        CALL alloc(positions,N,id='chosen eigenvectors for analysis')
        nr_vecs = 0

        upper = reladc_md_eeigv_upper
        lower = reladc_md_eeigv_lower

        DO j = 1, N
          IF ((EVL(j).GE.lower).AND.(EVL(j).LE.upper)) THEN
            nr_vecs = nr_vecs + 1
            positions(nr_vecs) = j
          END IF
        END DO

        NEIGENV = nr_vecs

        OPEN(IOUNIT,FILE='TMATEVC',FORM='UNFORMATTED',STATUS='UNKNOWN')
        WRITE(IOUNIT) nr_vecs   ! number of computed eigenvectors (all!)
        WRITE(IOUNIT) N   ! length of eigenvectors (= matrix dimension)
        DO k = 1, nr_vecs
          j = positions(k)
          WRITE(IOUNIT) EVL(J) ! eigenvalue
          WRITE(IOUNIT) (Z(IX,J),IX=1,N)  !corresponding eigenvector
        END DO
        CLOSE(IOUNIT)
    
        WRITE(iw,*) NEIGENV, 'Lanczos eigenvectors saved to TMATEVC'
        WRITE(iw,*) 'Energy range in a.u.: ', lower, '--', upper

        CALL dealloc(positions)

      END IF

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE IVXROT(IVX,N)
      IMPLICIT INTEGER (A-Z)
C
C  This routine cycles the mapping function according to a new
C  macro iteration IL in the calling program
C
      INTEGER IVX(N)
      INTEGER ITEMP

      ITEMP = IVX(1)
      DO I=1,N-1
        IVX(I)=IVX(I+1)
      ENDDO
      IVX(N)=ITEMP

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMATBLK(LADC,PO,PN,DIAG,BUF,IOI,IOJ,INTBUF,
     &                   IO_DIA,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_DIA,IO_OFF
      REAL*8 PO(LADC), PN(LADC), DIAG(LADC)
      REAL*8 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)
C
C---------------Description--------------------------------------------
C
C  Perform the real symmetric multiplication of the ADC matrix
C  with the Krylov vectors. The ADC matrix is real symmetric
C  This has to be taken into account.
C
C---------------Local variables--------------------------------------
C
      INTEGER NBUFS,IREC,NACT,JDUMMY,IROW,ICOL
      REAL*8 A
C
C---------------Executable code--------------------------------------

C
C  the diagonal elements
C
      REWIND(IO_DIA)
      READ(IO_DIA,ERR=88,END=99) (DIAG(IXX),IXX=1,LADC),NBUFS

      DO J = 1, LADC
        PN(J) = PO(J) * DIAG(J)
      ENDDO

C  for the off-diagonal elements we note that only the lower triangular
C  part of the ADC-matrix is stored. the missing contribution is from
C  the symmetric upper triangular part !

c  .. read off-diagonal elements and multiply on the fly

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)
            PN(IROW) = PN(IROW) + PO(ICOL) * A
            PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      RETURN

 88   CALL QUIT('Diagonal read error in RMATBLK')
 99   CALL QUIT('Unexpected end in RMATBLK')
 199  CALL QUIT('Off-diagonal read error in RMATBLK')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMATBLK_ND(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_OFF
      REAL*8 PO(LADC), PN(LADC)
      REAL*8 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)

C
C---------------Description--------------------------------------------
C
C  Perform the real symmetric multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix already and
C  a single multiplication suffices.
c
c  att! the required number of buffers is communicated through this
c  common block.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      REAL*8 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c  since the diagonal is now *contained* in the matrix
c
c  new vector has to be zeroed!
      DO J = 1, LADC
        PN(J) =  0.0d0
      ENDDO

c  .. read off-diagonal elements and multiply on the fly

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in RMATBLK_ND')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE RMATBLK_IC(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF,
     &                      nincw,bufc,ioic,iojc,ivalid)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_OFF
      REAL*8  PO(LADC), PN(LADC)
      REAL*8  BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)
      integer nincw
      real*8  bufc(nincw)
      integer ioic(nincw),iojc(nincw),ivalid

C
C---------------Description--------------------------------------------
C
C  Perform the real symmetric multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix already and
C  a single multiplication suffices.
c
c  att! we have the icore version here.  Matmul splits in an incore
c  part and an external part for those buffers that can not be held
c  in memory.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      REAL*8 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c  since the diagonal is now *contained* in the matrix

c  new vector has to be zeroed!

      DO J = 1, LADC
        PN(J) =  0.0d0
      ENDDO
c
c  do incore part. In this part of the program there *IS* an
c  incore contribution. So we do not need to check this condition
c  any further.
c
      DO K = 1, IVALID
        IROW = ioic(K)
        ICOL = iojc(K)
        A    = bufc(K)
        PN(IROW) = PN(IROW) + PO(ICOL) * A
        IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
      ENDDO
c
c  next do out of core part (if there is any)
c  remark: if matrix is completely incore the following loop
c  will never be executed. (nbufs = 0).
c  Under the IO_OFF file handle we find the *truncated* ADC
c  matrix only containing the ooc buffers. (there is no lseek
c  for this type of sequential files)
c
      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * A
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in RMATBLK_IC')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATBLK(LADC,PO,PN,DIAG,BUF,IOI,IOJ,INTBUF,
     &                   IO_DIA,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_DIA,IO_OFF
      COMPLEX*16 PO(LADC), PN(LADC), DIAG(LADC)
      COMPLEX*16 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)

C
C---------------Description--------------------------------------------
C
C  Perform the complex multiplication of the ADC matrix
C  with the Krylov vectors. The ADC matrix is Hermitian !
C  This has to be taken into account. The diagonal elements
C  are necessarily real.
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/complex.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER NBUFS,IREC,NACT,JDUMMY,IROW,ICOL
      COMPLEX*16 A
C
C---------------Executable code--------------------------------------

C
C  the diagonal elements
C
      REWIND(IO_DIA)
      READ(IO_DIA,ERR=88,END=99) (DIAG(IXX),IXX=1,LADC),NBUFS

      DO J = 1, LADC
        PN(J) = PO(J) * DIAG(J)
      ENDDO

C  for the off-diagonal elements we note that only the lower triangular
C  part of the ADC-matrix is stored. the missing contribution is from
C  the complex conjugate upper triangular part !

c  .. read off-diagonal elements

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=199) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)
            PN(IROW) = PN(IROW) + PO(ICOL) * A
            PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
        ENDDO
      ENDDO

      RETURN

 88   CALL QUIT('Diagonal read error in CMATBLK')
 99   CALL QUIT('Unexpected end in CMATBLK')
 199  CALL QUIT('Off-diagonal read error in CMATBLK')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATBLK_ND(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER LADC,INTBUF,IO_OFF
      COMPLEX*16 PO(LADC), PN(LADC)
      COMPLEX*16 BUF(INTBUF)
      INTEGER IOI(INTBUF),IOJ(INTBUF)

C
C---------------Description--------------------------------------------
C
C  Perform the complex hermitian multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix and
C  a single multiplication suffices.
c
c  att! the required number of buffers is communicated through this
c  common block.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      COMPLEX*16 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c
c  new vector has to be zeroed!
      DO J = 1, LADC
        PN(J) =  (0.0d0,0.0d0)
      ENDDO


c  .. read off-diagonal elements and multiply on the fly

      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in CMATBLK_ND')

      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMATBLK_IC(LADC,PO,PN,BUF,IOI,IOJ,INTBUF,IO_OFF,
     &                      nincw,bufc,ioic,iojc,ivalid)
C
      IMPLICIT INTEGER (A-Z)

      INTEGER     LADC,INTBUF,IO_OFF
      complex*16  PO(LADC), PN(LADC)
      complex*16  BUF(INTBUF)
      INTEGER     IOI(INTBUF),IOJ(INTBUF)
      integer     nincw
      complex*16  bufc(nincw)
      integer     ioic(nincw),iojc(nincw),ivalid

C
C---------------Description--------------------------------------------
C
C  Perform the complex hermitian multiplication of the ADC matrix
C  with the Krylov vectors. Here *NO* separate diagonal entries are
C  processed ! They are contained in the matrix already.
c
c  att! we have the icore version here.  Matmul splits in an incore
c  part and an external part for those buffers that can not be held
c  in memory.
C
C---------------Common Blocks--------------------------------------
C
      INTEGER NBUFS
      COMMON/MATCHUNKS/NBUFS
C
C---------------Local variables--------------------------------------
C
      INTEGER IREC,NACT,JDUMMY,IROW,ICOL
      complex*16 A
C
C---------------Executable code--------------------------------------

C  for the matrix elements we note that only the lower triangular
C  part *AND* the diagonal of the ADC-matrix is stored. 
C  the missing contribution is from the symmetric upper triangular part !

c  the if takes care of not double counting the diagonal.
c  since the diagonal is now *contained* in the matrix

c  new vector has to be zeroed!

      DO J = 1, LADC
        PN(J) = (0.0d0,0.0d0)
      ENDDO
c
c  do incore part. In this part of the program there *IS* an
c  incore contribution. So we do not need to check this condition
c  any further.
c
      DO K = 1, IVALID
        IROW = ioic(K)
        ICOL = iojc(K)
        A    = bufc(K)
        PN(IROW) = PN(IROW) + PO(ICOL) * A
        IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
      ENDDO
c
c  next do out of core part (if there is any)
c  remark: if matrix is completely incore the following loop
c  will never be executed. (nbufs = 0).
c  Under the IO_OFF file handle we find the *truncated* ADC
c  matrix only containing the ooc buffers. (there is no lseek
c  for this type of sequential files)
c
      REWIND(IO_OFF)
      DO IREC = 1,NBUFS
        READ(IO_OFF,ERR=889) (BUF(IXX),IXX=1,INTBUF),
     &                       (IOI(IXX),IXX=1,INTBUF),
     &                       (IOJ(IXX),IXX=1,INTBUF),
     &                        NACT,JDUMMY
        DO K = 1, NACT
          IROW = IOI(K)
          ICOL = IOJ(K)
          A    = BUF(K)

          PN(IROW) = PN(IROW) + PO(ICOL) * A
          IF(IROW.NE.ICOL) PN(ICOL) = PN(ICOL) + PO(IROW) * DCONJG(A)
        ENDDO
      ENDDO

      RETURN

 889  CALL QUIT('Matrix read error in CMATBLK_IC')

      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WCONDAT_S(LUN,FILENAME,IRECL,DESREP,LENADCM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes explicit irep names and spinor numbers for the various
C     final state configurations belonging to symmetry DESREP.
C     Avoiding too complicated things we fix IRECL in this routine
c     no matter which value is transferred to this routine. Therefore
c     IRECL is an OUTPUT parameter.
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*6 FILENAME
      INTEGER LUN,IRECL,DESREP,LENADCM
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER lcount
      CHARACTER*44 FIELD44
      INTEGER bl_conf
      CHARACTER*4 bl_irna,bl_abs
      
C
C---------------Executable code--------------------------------------
C
C  Code for ** single ionization **:  1h/2h1p
C
c fix record length for the direct access file and set some
c blank variables. Then open file and initialize configuration
c counter.
c *ATT* Is has turned out advantageous that storing the *absolute*
c spinor numbers is superior to storing the relative ones for the
c subsequent postprocessing!
c
c this will be followed by the single ionization case for the moment
c (1h and 2h1p) configurations
c
      IRECL = 44
      BL_CONF = 0
      BL_IRNA = '    '
      BL_ABS  = 'abs.'
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
      LCOUNT = 1

C_________________________ SIPS ______________________________
C|
C|
C|
        WRITE(IW,*) 'Writing SIP configs to ',FILENAME
        WRITE(IW,*) 'Final state symmetry: ',DESREP

c ** start with 1h configurations

        DO k=1,no(DESREP)
          koff=io(DESREP)+k
          WRITE(FIELD44,381) 
     &            k,repna(DESREP),
     &            koff,repna(DESREP),
     &            bl_conf,bl_irna,
     &            bl_conf,bl_irna 
          WRITE(LUN,REC=lcount) FIELD44
          lcount = lcount + 1
        ENDDO

c ** continue with 2h1p configurations
c ** write absolute spinor numbers

        DO 60 KLREP=1,NREP
          AREP=MULTB(DESREP,KLREP+NREP,2)
          DO 45 LREP=1,NREP
           KREP=MULTB(LREP,KLREP+NREP,2)
           IF(KREP.LT.LREP) GOTO 45
           DO L=1,NO(LREP)
             LOFF=IO(LREP)+L
             KMIN=1
             IF(KREP.EQ.LREP) KMIN = L + 1
             DO K=KMIN,NO(KREP)
               KOFF=IO(KREP)+K
               DO A=1,NV(AREP)
                 AOFF=IO(NREP+1)+IV(AREP)+A
c                WRITE(FIELD44,381) 
c    &             K,REPNA(KREP),
c    &             L,REPNA(LREP),
c    &             A,REPNA(AREP),
c    &             BL_CONF,BL_IRNA
                 WRITE(FIELD44,381) 
     &             KOFF,REPNA(KREP),
     &             LOFF,REPNA(LREP),
     &             AOFF,REPNA(AREP),
     &             BL_CONF,BL_IRNA
                 WRITE(LUN,REC=lcount) FIELD44
                 LCOUNT=LCOUNT+1
               ENDDO
             ENDDO
           ENDDO
 45       CONTINUE
 60     CONTINUE

C ** file is open only within this subroutine

        CLOSE(LUN)

C ** correct configuration counter

        LCOUNT = LCOUNT - 1
        IF(LCOUNT.NE.LENADCM) THEN
          WRITE(IW,*) 'LCOUNT,LENADCM:',LCOUNT,LENADCM
          CALL QUIT('Error in WCONDAT_S for single ionizations')
        ENDIF
C|
C|
C|____________________________________________________________

 381  FORMAT(4(I3,'(',A4,')  '))

      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE WCONDAT_D(LUN,FILENAME,IRECL,DESREP,LENADCM)
C
      IMPLICIT INTEGER (A-Z)
C
C---------------Description--------------------------------------------
C
C     Writes explicit irep names and spinor numbers for the various
C     final state configurations belonging to  symmetry DESREP.
C     Avoiding too complicated things we fix IRECL in this routine
c     no matter which value is transferred to this routine. Therefore
c     IRECL is an OUTPUT parameter.
C
C---------------Calling variables--------------------------------------
C
      CHARACTER*6 FILENAME
      INTEGER LUN,IRECL,DESREP,LENADCM
C
C---------------Common Blocks--------------------------------------
C
#include "../relccsd/files.inc"
#include "../relccsd/symm.inc"
C
C---------------Local variables--------------------------------------
C
      INTEGER lcount
      CHARACTER*44 FIELD44
      INTEGER bl_conf
      CHARACTER*4 bl_irna,bl_abs

      integer                                 :: MAXOCCSP
      integer, allocatable, dimension(:,:)    :: IROOO
      
C
C---------------Executable code--------------------------------------
C
C  Do the hole space configuration analysis for the
C  DIP case exclusively in this routine.
c  Initialize arrays:
c
      MAXOCCSP = IO(NREP+1)    ! IO array comes from symm.inc
      allocate(IROOO(MAXOCCSP,2))
      IROOO = 0

      IX=1
      Do KRP=1,NREP
        Do K=1,NO(KRP)
          IROOO(IX,1)=K
          IROOO(IX,2)=KRP
          IX = IX + 1
        Enddo
      Enddo
C
c fix record length for the direct access file and set some
c blank variables. Then open file and initialize configuration
c counter.
 
      IRECL = 44
      BL_CONF = 0
      BL_IRNA = '    '
      BL_ABS  = 'abs.'
      OPEN(LUN,FILE=FILENAME,ACCESS='DIRECT',RECL=IRECL,
     &     STATUS='UNKNOWN')
      LCOUNT = 1

C____________ START WITH DIPS ______________________________
C|
C|
C|
      WRITE(IW,*) 'Writing DIP configs to ',FILENAME
      WRITE(IW,*) 'Final state symmetry: ',DESREP

      DO 14 JRP = 1, NREP
        IRP = MULTB(JRP,DESREP+NREP,2)
        IF (IRP.LT.JRP) GOTO 14
        DO J = 1, NO(JRP)
          IMIN = 1
          IF (IRP.EQ.JRP) IMIN = J + 1
          DO I = IMIN, NO(IRP)
c
c  -- additionally determine absolute spinor numbers
c  -- in the active space for mulpop
c
            ISPINOR=IO(IRP)+I
            JSPINOR=IO(JRP)+J
            WRITE(FIELD44,381)
     &         I,REPNA(IRP),
     &         J,REPNA(JRP),
     &         ISPINOR,BL_ABS,
     &         JSPINOR,BL_ABS 
            WRITE(LUN,REC=LCOUNT) FIELD44
            LCOUNT = LCOUNT + 1
c
          ENDDO
        ENDDO
  14  CONTINUE
      IF( (LCOUNT-1).NE.NOOT(DESREP)) THEN
        WRITE(IW,*) 'LCOUNT,LENADCM:',LCOUNT,NOOT(DESREP)
        CALL QUIT('Error in WCONDAT_D!')
      ENDIF
c
c  determine satellite configurations
c  Since we work in a Kramers-paired formalism (closed shell)
c  then number of virtuals in RRP is equal to the number in its
c  Kramers irrep (RRP*). For the output RRP has to be the Kramers
c  partner since internally the cc of this orbital is used !
c
      DO K=1,MAXOCCSP
        KSP = IROOO(K,1)
        KRP = IROOO(K,2)
        DO J=K+1,MAXOCCSP
          JSP = IROOO(J,1)
          JRP = IROOO(J,2)
          JKRP = MULTB(JRP,KRP,1)
          DO I=J+1,MAXOCCSP
            ISP = IROOO(I,1)
            IRP = IROOO(I,2)
            IJKRP = MULTB(IRP,JKRP+NREP,1)
            RRP = MULTB(DESREP+NREP,IJKRP,2)
c determine Kramers-partner of RRP (is cc of RRP)
            RRPKR = MULTB(RRP,1+NREP,2)
            DO RSP = 1,NV(RRP)
            WRITE(FIELD44,381) ISP,REPNA(IRP),
     &                         JSP,REPNA(JRP),
     &                         KSP,REPNA(KRP),
     &                         RSP,REPNA(RRPKR)
            WRITE(LUN,REC=LCOUNT) FIELD44
            LCOUNT = LCOUNT + 1

            ENDDO
          ENDDO
        ENDDO
      ENDDO
      CLOSE(LUN)
C|
C|
C|____________________________________________________________

 381  FORMAT(4(I3,'(',A4,')  '))

      deallocate(IROOO)

      RETURN
      END
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
