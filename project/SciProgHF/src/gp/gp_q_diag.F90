!
! F90 file with general purpose routines for quaternion and other diagonalizations
!


#ifdef MOD_UNRELEASED

!> Program for testing (and improving) the quaternion matrix Jacobi (block) diagonalization procedure
!> written by Miro ILIAS, Strasbourg,2005, improved in Tel Aviv, 2012, 
!> Banska Bystrica, August 2014

module qjacobi_mod
      !> process also eigenvectors
      LOGICAL :: DOEVE 
      !> check the presence of the file with input data
      LOGICAL::  ISFILE 
      LOGICAL :: FROMFILE
      real(kind=8), allocatable :: A(:,:,:),EV(:,:,:)
      integer, allocatable :: IBL(:,:)
      integer :: N,NZ,NBL,IPRINT
      integer, parameter :: LUF=10
      real*8, PARAMETER :: D0=0.0D0
contains

subroutine GET_N_NZ
!> get in the basic variables:  N,NZ,NBL,IPRINT
implicit none
#include"priunit.h"
      INQUIRE(FILE="matrix_input",EXIST=ISFILE)
      IF (ISFILE) THEN
         write(lupri,*) '...going to read the "matrix_input" file' 
         open(luf,file="matrix_input",status="old")
         read(luf,*)
         read(luf,*) N,NZ
         read(luf,*)
         read(luf,*) IPRINT
         read(luf,*)
         read(luf,*) NBL
      ELSE
         write(lupri,*) 'No "matrix_input" file found, reading in from terminal.'
         write(lupri,*) 'Enter the N>0,NZ=1,2 or 4 :'
         read(*,*) N,NZ ! TODO: verify entering data
         NBL = 0
         write(lupri,*) 'Enter the NBL (zero for no blocks) :'
         read(*,*) NBL
         write(lupri,*) 'Enter IPRINT, print level :'
         read(*,*) IPRINT
      ENDIF
end subroutine

subroutine PREPMATRIX
!> ... prepare the testing matrix - allocate
implicit none
#include"priunit.h"
 integer :: I,J,IBL1,IBL2,IBL3,IBL4,I1,J1,IZ,istat
 real(kind=8) :: AR, AI,AJ,AK
      FROMFILE=.FALSE.
      DOEVE = .TRUE.
      IF (ISFILE) THEN
         write(lupri,*) 'PRINT level=',IPRINT
         write(lupri,*) 'Number of blocks=',NBL
         IF (NBL.GT.0) THEN
!         ... read the indexes of blocks
          DO I = 1, NBL
           read(luf,*) IBL1,IBL2,IBL3,IBL4
              IBL(I,1) = IBL1
              IBL(I,2) = IBL2 
              IBL(I,3) = IBL3 
              IBL(I,4) = IBL4
           IF (IPRINT.GE.2) THEN
            write(lupri,*) 'Read block coordinates:', & 
            IBL(I,1),IBL(I,2),IBL(I,3),IBL(I,4)
           ENDIF
          ENDDO
         ENDIF
         read(luf,*)
         read(luf,*) FROMFILE
         read(luf,*)
      ELSE
         IF (NBL.GT.0) THEN
          DO I=1,NBL
            write(*,*) 'Enter ',I,'. block coordinates:'
            read(*,*)  IBL(I,1),IBL(I,2),IBL(I,3),IBL(I,4)
          ENDDO
         ELSE
           write(*,*) 'No blocks....full diagonalization'
         ENDIF
      ENDIF

      DO I = 1, N
      DO J = I, N
          IF (FROMFILE) THEN
            AI = D0
            AJ = D0
            AK = D0
            IF (NZ.EQ.1) THEN
             read(luf,*) I1,J1,AR
            ELSE IF (NZ.EQ.2) THEN
             read(luf,*) I1,J1,AR,AI
            ELSE IF (NZ.EQ.4) THEN
             read(luf,*) I1,J1,AR,AI,AJ,AK
             If (IPRINT.GE.2) THEN
              write(lupri,*) 'read i1,j1,ar,ai,aj,ak:',I1,J1,AR,AI,AJ,AK
             EndIf
            ENDIF
            IF (I1.NE.I.OR.J1.NE.J) THEN
             write(lupri,*) '   read I,J:',I1,J1
             write(lupri,*) 'desired I,J:',I,J
             CALL QUIT('...wrong indexes !')
            ENDIF
            IF (I1.EQ.J1) THEN
             IF (AI.NE.D0.OR.AJ.NE.D0.OR.AK.NE.D0) THEN
              write(lupri,*) 'I1,J1 --> AI,AJ,AK:',I1,J1,AI,AJ,AK
              call QUIT('diagonal i,j,k elements must be zero !!!')
             ENDIF
            ENDIF
          ELSE
!          ... fill with random numbers
           IF (I.NE.J) THEN
            AR = RANDOMNUMBER(I,J)
            AI = RANDOMNUMBER(I,J)
            AJ = RANDOMNUMBER(I,J)
            AK = RANDOMNUMBER(I,J)
           ELSE
            AR = RANDOMNUMBER(I,J)
            AI = D0
            AJ = D0
            AK = D0
           ENDIF
          ENDIF

!       ....now fill the whole matrix 
          IF (NZ.EQ.1) THEN
            A(I,J,1) = AR
            A(J,I,1) = AR
          ELSE IF (NZ.EQ.2) THEN
            A(I,J,1) = AR
            A(J,I,1) = AR
            A(I,J,2) = AI
            A(J,I,2) = -AI
          ELSE IF (NZ.EQ.4) THEN
            A(I,J,1) =  AR
            A(J,I,1) =  AR
            A(I,J,2) =  AI
            A(J,I,2) = -AI
            A(I,J,3) =  AJ
            A(J,I,3) = -AJ
            A(I,J,4) =  AK
            A(J,I,4) = -AK
            if (IPRINT.GE.10) then
             write(lupri,*) 'A filled:', & 
             I,J,':',A(I,J,1),A(I,J,2),A(I,J,3),A(I,J,4)
            endif
          ENDIF
          IF (IPRINT.GE.3) THEN
            write(lupri,'(a,i3,a,i3,a,4d12.5)')  & 
           'random number in matrix element A(',I,';',J,')=',  & 
           (A(I,J,IZ),IZ=1,NZ)
          ENDIF
      ENDDO 
      ENDDO 
      IF (ISFILE) CLOSE(LUF,STATUS='KEEP')
end subroutine

REAL*8 FUNCTION RANDOMNUMBER(I,J)
!> Function used to generate pseudorandom numbers, also according to 
!> entering indexes I,J
! See also http://gcc.gnu.org/onlinedocs/gfortran/RAND.html
implicit none
#include"priunit.h"
integer, intent(in) :: I,J
integer (kind=4) :: rand_param=1
real(kind=8), external :: RAND
  !print *,'rand=',RAND(rand_param)
  !RANDOMNUMBER = DFLOAT(3*I+4*J)*RAND(rand_param) 
 ! miro: deactivated RAND, as this is not universally working
  RANDOMNUMBER = DFLOAT(3*I+4*J)
END FUNCTION
end module qjacobi_mod


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!  /* Deck qjacobi */
        SUBROUTINE QJACOBI(A,EV,N,NZ,NBL,IBL,DOEVE,IPRINT) 
!*************************************************************************
!
!  PURPOSE:  Do the selective (or full) diagonalization of 
!            the entering hermitian (real,complex,quaternion) matrix
!            using the Jacobi diagonalization method, which is  extended 
!            for quaternion matrixes.
!
!  On input:  A(N,N,NZ) ... quaternion matrix (is destroyed after!)
!             DOEVE ... if .true., calculate eigenvectors,
!                   and for IPRINT.GE.2 also do the consistency check 
!
!             NBL=0 ... do the traditional diagonalization (parameter IBL substitute with dummy variable,IDUMMY)
!             NBL>0 ... do the block diagonalization, where blocks to be zeored are specified in IBL
!
!             IPRINT - print level 
!  1 ... print timing and basic info
!  2 ... print out entering matrix into text file, reconstruct it
!        back at using eigenvalues and eigenvectors (if DOEVE=.true.)
!  3 ... + print out entering and final matrix, print out result fo Ph*P,P*Ph,
!          print out iteration process, eigenvalues (if DOEVE=.true.)
!  4 ... +  print eigenvectors (if DOEVE=.true.)
!  5 ... +  print whole matrix in each iteration
!  6 ... +  print whole P after iter., print A before and after iter.
!  7 ... +  print Ppq
!                       
!  On output: A - diagonal (NBL=0) or block diagonal (NBL>1) matrix
!             so eigenvalues are on the diagonal of the (destroyed) A matrix
!
!             EV(N,N,NZ) - matrix of eigenvectors if DOEVE=.TRUE.
!
!  Written: M Ilias, Strasbourg, May 2005;
!  MI, August 2014 - rearragned for F90 
!      
!**************************************************************************
use include_dgroup_h, only: get_IPQTOQ
implicit none                   
#include"priunit.h"
  real*8,PARAMETER :: D1=1.0D00, DM1=-1.00D00, D0=0.0D00,D2=2.0D0, DP5=0.50D0
  integer, intent(in) :: N, NZ, NBL
  real*8, intent(inout) :: A(N,N,NZ)
  real*8, intent(inout) :: EV(N,N,NZ)
  integer, intent(in) :: IBL(NBL,1:4), IPRINT
  real*8 :: SUMSQA(2)
  LOGICAL::  DOEVE, FINISHED, WASELIM, ELIM(2)
  CHARACTER*9 :: FILENAME
  CHARACTER  :: SECTID*12,WALLTID*12, CPUTID*12
  real*8, parameter :: ZEROTHRS= 1.0D-12, ZEROTHRS1 = 1.0D-9  
  !      .... file number and file name of the check file
  integer, parameter :: LUF = 33 
  DATA FILENAME / 'QJACO.MTX' /

  real*8 :: CPU1,WALL1,SUMSQ,APQR,APPR,BDA,PHI,CU,SU,APQJ,AQQR,CU2,SU2
  real*8 :: TOAIPR,TOAIQR,TOPIPR,TOAIPI,TOAIQI,APQI,XAI,TOAIPJ,CUSU
  real*8 :: TOPIQR,TOPIPI,TOPIQI,TOAIPK,TOPIQJ,TOPIPK,TOPIQK,XAJ,APQK
  real*8 :: CPU2,WALL2,DPHKJR,TOAIQJ,TOAIQK,TOPIPJ,XAK,CPU,WALL,AIJ_R
  real*8 :: APPI,TOQIPR,AIJFR,AIJ_I,TOQIPI,AIJ_J,AIJ_K
  real*8 :: DPHKJI,DPHKJJ,AIJFI,TOQIPK,TOQIPJ,DPHKJK,AIJFJ,AIJFK
  integer :: IP,IQ,I,J,K,L,IZ,ITER,I1,J1,IZ1

      IF (IPRINT.GE.2) THEN
       CALL HEADER(                                                     & 
 'QJACOBI, selective diagonalization of the quaternion matrix.',-1)
       write(lupri,'(2X,A,I5)')                                         &
   'size of the square hermitian/real symmetric matrix, N=',N
       write(lupri,'(2X,A,I2)')                                         &
   'type of matrix {real(NZ=1), complex(NZ=2) or'//                     &
   ' quaternion(NZ=4)}, NZ=',NZ
       write(lupri,'(2X,A,I5)')                                         &
   'PRINT level:',IPRINT
       CALL FLSHFO(LUPRI)
       CALL GETTIM(CPU1,WALL1)
      ENDIF

!     ... do test the blocks
      IF (NBL.GT.0) THEN
       CALL TESTBLOCKS(N,NZ,NBL,IBL,IPRINT)
      ENDIF

      CALL FINDSUMBL(A,N,NZ,NBL,IBL,IP,IQ,SUMSQ,IPRINT)
      SUMSQA(1)=SUMSQ
      SUMSQA(2)=SUMSQ

!      ... do make unit matrixes for accumulation of future eigenvector
      IF (DOEVE) THEN
        CALL DUNIT(EV(1,1,1),N)
        IF (NZ.GE.2) THEN 
         DO IZ=2,NZ
          CALL DZERO(EV(1,1,IZ),N*N)
         ENDDO
        ENDIF
      ENDIF

      IF (IPRINT.GE.2) THEN

        IF (DOEVE) THEN
            OPEN(LUF,FILE=FILENAME,STATUS="UNKNOWN",FORM="FORMATTED")
            write(lupri,'(2X,A,I2,A,A,A)')                              &
       'The file #',LUF,' with name:',FILENAME,                         &
       ' has been open.'

            IF (IPRINT.GE.3) THEN
             WRITE(LUPRI,'(/,2X,A)')                                    &
  'QJACOBI: ENTERING hermitian(real,complex,quaternion) matrix A'
             CALL PRQMAT(A,N,N,N,N,NZ,get_IPQTOQ(1,0),LUPRI) ! IPQTOQ in dgroup.h, this contain NZ
            ENDIF

            DO I=1,N
            DO J=1,N
            IF (IPRINT.GE.3) THEN
             write(lupri,'(A,I3,A,I3,A,4D15.8)')   &
          '>>>>>  A(',I,',',J,') =',(A(I,J,IZ),IZ=1,NZ)
            ENDIF
!            ... write matrix elements to file for later checking...
             DO IZ=1,NZ
              write(LUF,1000) I,J,IZ,A(I,J,IZ)
             ENDDO
!            write(luf,*) I,J,(A(I,J,IZ),IZ=1,NZ)
            ENDDO
            ENDDO
        ENDIF
         write(lupri,'(2X,A,A)')                                        &
  '...entering matrix has been written to the file ',FILENAME
      ENDIF

!  ... format of elements in the file...
 1000 FORMAT(2I5,2X,I1,2X,D30.17) 
  
         FINISHED = .false.
         ITER = 0
         WASELIM = .FALSE.
         ELIM(1) = .FALSE.
         ELIM(2) = .FALSE.
! *****************************************************************
! ***
! ***
! ***           Start the iteration cycle
! ***
! ***
! *****************************************************************
 10      CONTINUE

! ============================================================================
!
!                   Real (symmetric) matrix
! 
! ============================================================================
         IF (NZ.EQ.1) THEN
!          ... eliminate real part of P,Q
25         CONTINUE

!           ... write out the A_pq  elements !
           IF (IPRINT.GE.6) THEN
             write(lupri,'(/A,2I5)')                                    &
    "======  BEFORE elimination App,Apq,Aqp,Aqq r:",IP,IQ
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1)
           ENDIF

           APQR =  A(IP,IQ,1)

           IF (DABS(APQR).GT.ZEROTHRS) THEN

            ITER = ITER + 1
            WASELIM = .TRUE.

            APPR =  A(IP,IP,1)
            AQQR =  A(IQ,IQ,1)

!          ... an alternative to the calculation using function DATAN ...
!           APQSQR = APQR*APQR
!           DIFF = APPR - AQQR
!           SIGROT = D1
!           IF (DIFF.LT.D0) THEN
!            SIGROT = -SIGROT
!            DIFF = -DIFF
!           ENDIF
!           TEMP = DIFF + DSQRT(DIFF*DIFF + 4.0*APQSQR)
!           TANA = D2*SIGROT*APQR/TEMP
!           CU =  D1/DSQRT(D1+TANA*TANA)
!           SU = -CU*TANA

!         ... do calculate Cu,Su...
            BDA =  D2*APQR/(AQQR-APPR)
            PHI = DP5*DATAN(BDA)
            CU = DCOS(PHI)
            SU = DSIN(PHI)

            CU2 = CU*CU
            SU2 = SU*SU
            CUSU = CU*SU

!           ... Update REAL diagonal  MTX elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQR 
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQR

!...        do eliminate REAL off-diagonal pq, qp elements
!               .... (1,i).r = (1,i) = (0,i)
            A(IP,IQ,1)=D0
!           A(IQ,IP,1)=D0

!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After r-elim App,Apq,Aqp,Aqq r:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1)
            ENDIF

!   ...   do modify the MTX matrix alement A_ip, A_iq (A_pi,A_qi)
            DO I=1,IP-1
!               I < IP, I < IQ
               TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               A(I,IP,1) =  TOAIPR 
               A(I,IQ,1) =  TOAIQR
!              A(IP,I,1) =  TOAIPR
!              A(IQ,I,1) =  TOAIQR
            ENDDO
            DO I=IP+1,IQ-1
!               I > IP, I < IQ
!C             TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIPR = A(IP,I,1)*CU - A(I,IQ,1)*SU
!C             TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIQR = A(IP,I,1)*SU + A(I,IQ,1)*CU
!              A(I,IP,1) =  TOAIPR 
               A(I,IQ,1) =  TOAIQR
               A(IP,I,1) =  TOAIPR
!              A(IQ,I,1) =  TOAIQR
            ENDDO
            DO I=IQ+1,N
!               I > IP, I > IQ
!C             TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIPR = A(IP,I,1)*CU - A(IQ,I,1)*SU
!C             TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIQR = A(IP,I,1)*SU + A(IQ,I,1)*CU
!              A(I,IP,1) =  TOAIPR 
!              A(I,IQ,1) =  TOAIQR
               A(IP,I,1) =  TOAIPR
               A(IQ,I,1) =  TOAIQR
            ENDDO

!          ... do update the eigenvectors...
            IF (DOEVE) THEN
               DO I=1,N
                TOPIPR = EV(I,IP,1)*CU - EV(I,IQ,1)*SU
                TOPIQR = EV(I,IP,1)*SU + EV(I,IQ,1)*CU
                EV(I,IP,1) = TOPIPR
                EV(I,IQ,1) = TOPIQR
               ENDDO

              IF (IPRINT.GE.7) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvector after Apq_"r" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r:',                &
          EV(I,J,1)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

           ENDIF

           ENDIF

! ============================================================================
!
!                   Complex (hermitian) matrix
! 
! ============================================================================
         ELSE IF (NZ.EQ.2) THEN
!          ... eliminate real part of P,Q
15         CONTINUE

!           ... write out the A_pq  elements !
           IF (IPRINT.GE.6) THEN
             write(lupri,'(/A,2I5)')                                    &
    "======  BEFORE elimination App,Apq,Aqp,Aqq r,i:",IP,IQ
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2)
           ENDIF

           APQR =  A(IP,IQ,1)

           IF (DABS(APQR).GT.ZEROTHRS) THEN

            ITER = ITER + 1
            WASELIM = .TRUE.

            APPR =  A(IP,IP,1)
            AQQR =  A(IQ,IQ,1)

!         ... do calculate Cu,Su...
            BDA =  D2*APQR/(AQQR-APPR)

            PHI = DP5*DATAN(BDA)
            CU = DCOS(PHI)
            SU = DSIN(PHI)
            CU2 = CU*CU
            SU2 = SU*SU
            CUSU = CU*SU

!           ... Update REAL diagonal  MTX elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQR 
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQR

!...        do eliminate REAL off-diagonal pq, qp elements
!               .... (1,i).r = (1,i) = (0,i)
            A(IP,IQ,1)=D0
!           A(IQ,IP,1)=D0

!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After r-elim App,Apq,Aqp,Aqq r,i:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2)
            ENDIF

!   ...   do modify the MTX matrix alement A_ip, A_iq (A_pi,A_qi)
            DO I=1,IP-1
!               I < IP, I < IQ
               TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIPI = A(I,IP,2)*CU - A(I,IQ,2)*SU
               TOAIQI = A(I,IP,2)*SU + A(I,IQ,2)*CU
                A(I,IP,1) =  TOAIPR 
                A(I,IQ,1) =  TOAIQR
                A(I,IP,2) =  TOAIPI
                A(I,IQ,2) =  TOAIQI
!               A(IP,I,1) =  TOAIPR
!               A(IQ,I,1) =  TOAIQR
!               A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
            ENDDO
            DO I=IP+1,IQ-1
!              I > IP, I < IQ
!              TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIPR = A(IP,I,1)*CU - A(I,IQ,1)*SU
!              TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIQR = A(IP,I,1)*SU + A(I,IQ,1)*CU
!              TOAIPI =  A(I,IP,2)*CU - A(I,IQ,2)*SU
               TOAIPI = -A(IP,I,2)*CU - A(I,IQ,2)*SU
!              TOAIQI =  A(I,IP,2)*SU + A(I,IQ,2)*CU
               TOAIQI = -A(IP,I,2)*SU + A(I,IQ,2)*CU
!               A(I,IP,1) =  TOAIPR 
                A(I,IQ,1) =  TOAIQR
!               A(I,IP,2) =  TOAIPI
                A(I,IQ,2) =  TOAIQI
                A(IP,I,1) =  TOAIPR
!               A(IQ,I,1) =  TOAIQR
                A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
            ENDDO
            DO I=IQ+1,N
!               I > IP, I > IQ
!              TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIPR = A(IP,I,1)*CU - A(IQ,I,1)*SU
!              TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIQR = A(IP,I,1)*SU + A(IQ,I,1)*CU
!              TOAIPI =  A(I,IP,2)*CU - A(I,IQ,2)*SU
               TOAIPI = -A(IP,I,2)*CU + A(IQ,I,2)*SU
!              TOAIQI =  A(I,IP,2)*SU + A(I,IQ,2)*CU
               TOAIQI = -A(IP,I,2)*SU - A(IQ,I,2)*CU
!               A(I,IP,1) =  TOAIPR 
!               A(I,IQ,1) =  TOAIQR
!               A(I,IP,2) =  TOAIPI
!               A(I,IQ,2) =  TOAIQI
                A(IP,I,1) =  TOAIPR
                A(IQ,I,1) =  TOAIQR
                A(IP,I,2) = -TOAIPI
                A(IQ,I,2) = -TOAIQI
            ENDDO

!          ... do update the eigenvectors...
            IF (DOEVE) THEN
               DO I=1,N
                 TOPIPR = EV(I,IP,1)*CU - EV(I,IQ,1)*SU
                 TOPIQR = EV(I,IP,1)*SU + EV(I,IQ,1)*CU

                 TOPIPI = EV(I,IP,2)*CU - EV(I,IQ,2)*SU
                 TOPIQI = EV(I,IP,2)*SU + EV(I,IQ,2)*CU

                 EV(I,IP,1) = TOPIPR
                 EV(I,IQ,1) = TOPIQR

                 EV(I,IP,2) = TOPIPI
                 EV(I,IQ,2) = TOPIQI

               ENDDO

              IF (IPRINT.GE.7) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvector after Apq_"r" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r,i,j,k:',          &
          EV(I,J,1),EV(I,J,2)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

           ENDIF

           ENDIF

!     =========================================
!          ... eliminate "i" part of Apq, Aqp
!     =========================================

16         APQI =  A(IP,IQ,2)

           IF (DABS(APQI).GT.ZEROTHRS) THEN

            ITER = ITER + 1
            WASELIM = .TRUE.

            APPR =  A(IP,IP,1)
            AQQR =  A(IQ,IQ,1)

!         ... do calculate Cu,Su...
           BDA =  D2*APQI/(AQQR-APPR)

           PHI = DP5*DATAN(BDA)
           CU = DCOS(PHI)
           SU = DSIN(PHI)
           CU2 = CU*CU
           SU2 = SU*SU
           CUSU = CU*SU


!           ... Update (only REAL) diagonal elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQI 
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQI

!...        do eliminate "i" off-diagonal IMAX,JMAX elements
!           A(IP,IQ,2)=D0
!           A(IQ,IP,2)=D0

!      ...  update-rearrange off diag.element  Apq 
! ...          (1,i).i = (i,-1) = (i,0)
            XAI = A(IP,IQ,1)

            A(IP,IQ,2)= XAI
!           A(IQ,IP,2)=-XAI
               
!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After i-elim App,Apq,Aqp,Aqq r,i:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2)
            ENDIF

!      ...  update Aip, Aiq (i<>p,q) after elimination!
! .....   A_ip(n+1) = A_ip(n).i.c - A_iq(n).s    
!         A_iq(n+1) = A_ip(n).i.s + A_iq(n).c
! 
!                  (1,i).i = (i,-1)
            DO I=1,IP-1
!               I < IP, I < IQ
               TOAIPR = -A(I,IP,2)*CU - A(I,IQ,1)*SU
               TOAIQR = -A(I,IP,2)*SU + A(I,IQ,1)*CU
               TOAIPI =  A(I,IP,1)*CU - A(I,IQ,2)*SU
               TOAIQI =  A(I,IP,1)*SU + A(I,IQ,2)*CU
               A(I,IP,1) = TOAIPR 
               A(I,IQ,1) = TOAIQR
               A(I,IP,2) = TOAIPI
               A(I,IQ,2) = TOAIQI
!              A(IP,I,1) = TOAIPR
!              A(IQ,I,1) = TOAIQR
!              A(IP,I,2) = -TOAIPI
!              A(IQ,I,2) = -TOAIQI
            ENDDO
            DO I=IP+1,IQ-1
!              I > IP, I < IQ
!              TOAIPR = -A(I,IP,2)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,2)*CU - A(I,IQ,1)*SU
!              TOAIQR = -A(I,IP,2)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,2)*SU + A(I,IQ,1)*CU
!              TOAIPI = A(I,IP,1)*CU - A(I,IQ,2)*SU
               TOAIPI = A(IP,I,1)*CU - A(I,IQ,2)*SU
!              TOAIQI = A(I,IP,1)*SU + A(I,IQ,2)*CU
               TOAIQI = A(IP,I,1)*SU + A(I,IQ,2)*CU
!              A(I,IP,1) = TOAIPR 
               A(I,IQ,1) = TOAIQR
!              A(I,IP,2) = TOAIPI
               A(I,IQ,2) = TOAIQI
               A(IP,I,1) = TOAIPR
!              A(IQ,I,1) = TOAIQR
               A(IP,I,2) = -TOAIPI
!              A(IQ,I,2) = -TOAIQI
            ENDDO
            DO I=IQ+1,N
!              I > IP, I > IQ
!              TOAIPR = -A(I,IP,2)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,2)*CU - A(IQ,I,1)*SU
!              TOAIQR = -A(I,IP,2)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,2)*SU + A(IQ,I,1)*CU
!              TOAIPI = A(I,IP,1)*CU - A(I,IQ,2)*SU
               TOAIPI = A(IP,I,1)*CU + A(IQ,I,2)*SU
!              TOAIQI = A(I,IP,1)*SU + A(I,IQ,2)*CU
               TOAIQI = A(IP,I,1)*SU - A(IQ,I,2)*CU
!              A(I,IP,1) = TOAIPR 
!              A(I,IQ,1) = TOAIQR
!              A(I,IP,2) = TOAIPI
!              A(I,IQ,2) = TOAIQI
               A(IP,I,1) = TOAIPR
               A(IQ,I,1) = TOAIQR
               A(IP,I,2) = -TOAIPI
               A(IQ,I,2) = -TOAIQI
            ENDDO

!      ...  update P after "i" elimination!
            IF (DOEVE) THEN
              
               DO I=1,N

!          ... P_ip(n+1) = P_ip(n).i.c - P_iq(n).s
!              P_iq(n+1) = P_ip(n).i.s + P_iq(n).c   
!
!               ... (1,i).i  =  (i,-1)

                 TOPIPR = -EV(I,IP,2)*CU - EV(I,IQ,1)*SU
                 TOPIQR = -EV(I,IP,2)*SU + EV(I,IQ,1)*CU

                 TOPIPI =  EV(I,IP,1)*CU - EV(I,IQ,2)*SU
                 TOPIQI =  EV(I,IP,1)*SU + EV(I,IQ,2)*CU

                 EV(I,IP,1) = TOPIPR
                 EV(I,IQ,1) = TOPIQR

                 EV(I,IP,2) = TOPIPI
                 EV(I,IQ,2) = TOPIQI

             ENDDO

             IF (IPRINT.GE.7) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvector after Apq_"i" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r,i:',              &
          EV(I,J,1),EV(I,J,2)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

         ENDIF

         ENDIF

         ELSE IF (NZ.EQ.4) THEN
! =======================================================================
!
!                   Quaternion matrix
!
! =======================================================================

!          ... eliminate real part of P,Q
 5         CONTINUE

!           ... write out the A_pq  elements !
           IF (IPRINT.GE.6) THEN
             write(lupri,'(/A,2I5)')                                    &
    "======  BEFORE elimination App,Apq,Aqp,Aqq r,i,j,k:",IP,IQ
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2),A(IP,IP,3),A(IP,IP,4)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2),A(IP,IQ,3),A(IP,IQ,4)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2),A(IQ,IP,3),A(IQ,IP,4)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2),A(IQ,IQ,3),A(IQ,IQ,4)
           ENDIF

           APQR =  A(IP,IQ,1)

           IF (DABS(APQR).GT.ZEROTHRS) THEN

            ITER = ITER + 1
            WASELIM = .TRUE.

            APPR =  A(IP,IP,1)
            AQQR =  A(IQ,IQ,1)

!         ... do calculate Cu,Su...
            BDA =  D2*APQR/(AQQR-APPR)

            PHI = DP5*DATAN(BDA)
            CU = DCOS(PHI)
            SU = DSIN(PHI)
            CU2 = CU*CU
            SU2 = SU*SU
            CUSU = CU*SU

!           ... Update REAL diagonal  MTX elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQR 
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQR

!...        do eliminate REAL off-diagonal pq, qp elements
!               .... (1,i,j,k).r = (1,i,j,k) = (0,i,j,k)
            A(IP,IQ,1)=D0
!           A(IQ,IP,1)=D0

!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After r-elim App,Apq,Aqp,Aqq r,i,j,k:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2),A(IP,IP,3),A(IP,IP,4)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2),A(IP,IQ,3),A(IP,IQ,4)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2),A(IQ,IP,3),A(IQ,IP,4)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2),A(IQ,IQ,3),A(IQ,IQ,4)
            ENDIF

!   ...   do modify the MTX matrix alement A_ip, A_iq (A_pi,A_qi)
            DO I=1, IP-1
!            I < IP, I < IQ
               TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIPI = A(I,IP,2)*CU - A(I,IQ,2)*SU
               TOAIQI = A(I,IP,2)*SU + A(I,IQ,2)*CU
               TOAIPJ = A(I,IP,3)*CU - A(I,IQ,3)*SU
               TOAIQJ = A(I,IP,3)*SU + A(I,IQ,3)*CU
               TOAIPK = A(I,IP,4)*CU - A(I,IQ,4)*SU
               TOAIQK = A(I,IP,4)*SU + A(I,IQ,4)*CU
                A(I,IP,1) =  TOAIPR 
                A(I,IQ,1) =  TOAIQR
                A(I,IP,2) =  TOAIPI
                A(I,IQ,2) =  TOAIQI
                A(I,IP,3) =  TOAIPJ
                A(I,IQ,3) =  TOAIQJ
                A(I,IP,4) =  TOAIPK
                A(I,IQ,4) =  TOAIQK
!               A(IP,I,1) =  TOAIPR
!               A(IQ,I,1) =  TOAIQR
!               A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
!               A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
!               A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IP+1,IQ-1
!               I > IP,  I < IQ
!              TOAIPR = A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIPR = A(IP,I,1)*CU - A(I,IQ,1)*SU
!              TOAIQR = A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIQR = A(IP,I,1)*SU + A(I,IQ,1)*CU
!              TOAIPI =  A(I,IP,2)*CU - A(I,IQ,2)*SU
               TOAIPI = -A(IP,I,2)*CU - A(I,IQ,2)*SU
!              TOAIQI =  A(I,IP,2)*SU + A(I,IQ,2)*CU
               TOAIQI = -A(IP,I,2)*SU + A(I,IQ,2)*CU
!              TOAIPJ =  A(I,IP,3)*CU - A(I,IQ,3)*SU
               TOAIPJ = -A(IP,I,3)*CU - A(I,IQ,3)*SU
!              TOAIQJ =  A(I,IP,3)*SU + A(I,IQ,3)*CU
               TOAIQJ = -A(IP,I,3)*SU + A(I,IQ,3)*CU
!              TOAIPK =  A(I,IP,4)*CU - A(I,IQ,4)*SU
               TOAIPK = -A(IP,I,4)*CU - A(I,IQ,4)*SU
!              TOAIQK =  A(I,IP,4)*SU + A(I,IQ,4)*CU
               TOAIQK = -A(IP,I,4)*SU + A(I,IQ,4)*CU
!               A(I,IP,1) =  TOAIPR 
                A(I,IQ,1) =  TOAIQR
!               A(I,IP,2) =  TOAIPI
                A(I,IQ,2) =  TOAIQI
!               A(I,IP,3) =  TOAIPJ
                A(I,IQ,3) =  TOAIQJ
!               A(I,IP,4) =  TOAIPK
                A(I,IQ,4) =  TOAIQK
                A(IP,I,1) =  TOAIPR
!               A(IQ,I,1) =  TOAIQR
                A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IQ+1,N
!               I > IP, I > IQ
!              TOAIPR =  A(I,IP,1)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,1)*CU - A(IQ,I,1)*SU
!              TOAIQR =  A(I,IP,1)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,1)*SU + A(IQ,I,1)*CU
!              TOAIPI =  A(I,IP,2)*CU - A(I,IQ,2)*SU
               TOAIPI = -A(IP,I,2)*CU + A(IQ,I,2)*SU
!              TOAIQI =  A(I,IP,2)*SU + A(I,IQ,2)*CU
               TOAIQI = -A(IP,I,2)*SU - A(IQ,I,2)*CU
!              TOAIPJ =  A(I,IP,3)*CU - A(I,IQ,3)*SU
               TOAIPJ = -A(IP,I,3)*CU + A(IQ,I,3)*SU
!              TOAIQJ =  A(I,IP,3)*SU + A(I,IQ,3)*CU
               TOAIQJ = -A(IP,I,3)*SU - A(IQ,I,3)*CU
!              TOAIPK =  A(I,IP,4)*CU - A(I,IQ,4)*SU
               TOAIPK = -A(IP,I,4)*CU + A(IQ,I,4)*SU
!              TOAIQK =  A(I,IP,4)*SU + A(I,IQ,4)*CU
               TOAIQK = -A(IP,I,4)*SU - A(IQ,I,4)*CU
!               A(I,IP,1) =  TOAIPR 
!               A(I,IQ,1) =  TOAIQR
!               A(I,IP,2) =  TOAIPI
!               A(I,IQ,2) =  TOAIQI
!               A(I,IP,3) =  TOAIPJ
!               A(I,IQ,3) =  TOAIQJ
!               A(I,IP,4) =  TOAIPK
!               A(I,IQ,4) =  TOAIQK
                A(IP,I,1) =  TOAIPR
                A(IQ,I,1) =  TOAIQR
                A(IP,I,2) = -TOAIPI
                A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
                A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
                A(IQ,I,4) = -TOAIQK
            ENDDO

!          ... do update the eigenvectors...
            IF (DOEVE) THEN
               DO I=1,N

                 TOPIPR = EV(I,IP,1)*CU - EV(I,IQ,1)*SU
                 TOPIQR = EV(I,IP,1)*SU + EV(I,IQ,1)*CU

                 TOPIPI = EV(I,IP,2)*CU - EV(I,IQ,2)*SU
                 TOPIQI = EV(I,IP,2)*SU + EV(I,IQ,2)*CU

                 TOPIPJ = EV(I,IP,3)*CU - EV(I,IQ,3)*SU
                 TOPIQJ = EV(I,IP,3)*SU + EV(I,IQ,3)*CU

                 TOPIPK = EV(I,IP,4)*CU - EV(I,IQ,4)*SU
                 TOPIQK = EV(I,IP,4)*SU + EV(I,IQ,4)*CU

                 EV(I,IP,1) = TOPIPR
                 EV(I,IQ,1) = TOPIQR

                 EV(I,IP,2) = TOPIPI
                 EV(I,IQ,2) = TOPIQI

                 EV(I,IP,3) = TOPIPJ
                 EV(I,IQ,3) = TOPIQJ

                 EV(I,IP,4) = TOPIPK
                 EV(I,IQ,4) = TOPIQK

               ENDDO

              IF (IPRINT.GE.7) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvector after Apq_"r" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r,i,j,k:',          &
          EV(I,J,1),EV(I,J,2),EV(I,J,3),EV(I,J,4)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

           ENDIF

           ENDIF

!     =========================================
!          ... eliminate "i" part of Apq, Aqp
!     =========================================

 6         APQI =  A(IP,IQ,2)

           IF (DABS(APQI).GT.ZEROTHRS) THEN

            ITER = ITER + 1
            WASELIM = .TRUE.

            APPR =  A(IP,IP,1)
            AQQR =  A(IQ,IQ,1)

!         ... do calculate Cu,Su...
           BDA =  D2*APQI/(AQQR-APPR)

           PHI = DP5*DATAN(BDA)
           CU = DCOS(PHI)
           SU = DSIN(PHI)
           CU2 = CU*CU
           SU2 = SU*SU
           CUSU = CU*SU


!           ... Update (only REAL) diagonal elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQI 
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQI

!...        do eliminate "i" off-diagonal IMAX,JMAX elements
!           A(IP,IQ,2)=D0
!           A(IQ,IP,2)=D0

!      ...  update-rearrange off diag.element  Apq 
! ...          (1,i,j,k).i = (i,-1,-k,+j) = (i,0,-k,j)
            XAI = A(IP,IQ,1)
            XAJ = A(IP,IQ,4)
            XAK = A(IP,IQ,3)

            A(IP,IQ,2)= XAI
            A(IP,IQ,3)= XAJ
            A(IP,IQ,4)=-XAK

!           ... the same for qp elements...
!           A(IQ,IP,2)=-XAI
!           A(IQ,IP,3)=-XAJ
!           A(IQ,IP,4)= XAK
               
!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After i-elim App,Apq,Aqp,Aqq r,i,j,k:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2),A(IP,IP,3),A(IP,IP,4)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2),A(IP,IQ,3),A(IP,IQ,4)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2),A(IQ,IP,3),A(IQ,IP,4)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2),A(IQ,IQ,3),A(IQ,IQ,4)
            ENDIF

!      ...  update Aip, Aiq (i<>p,q) after elimination!
! .....   A_ip(n+1) = A_ip(n).i.c - A_iq(n).s    
!         A_iq(n+1) = A_ip(n).i.s + A_iq(n).c
! 
!                  (1,i,j,k).i = (i,-1,-k,j)
            DO I=1,IP-1
!             I  <  IP,  I < IQ
               TOAIPR = -A(I,IP,2)*CU - A(I,IQ,1)*SU
               TOAIQR = -A(I,IP,2)*SU + A(I,IQ,1)*CU
               TOAIPI = A(I,IP,1)*CU - A(I,IQ,2)*SU
               TOAIQI = A(I,IP,1)*SU + A(I,IQ,2)*CU
               TOAIPJ = A(I,IP,4)*CU - A(I,IQ,3)*SU
               TOAIQJ = A(I,IP,4)*SU + A(I,IQ,3)*CU
               TOAIPK = -A(I,IP,3)*CU - A(I,IQ,4)*SU
               TOAIQK = -A(I,IP,3)*SU + A(I,IQ,4)*CU
                A(I,IP,1) = TOAIPR 
                A(I,IQ,1) = TOAIQR
                A(I,IP,2) = TOAIPI
                A(I,IQ,2) = TOAIQI
                A(I,IP,3) = TOAIPJ
                A(I,IQ,3) = TOAIQJ
                A(I,IP,4) = TOAIPK
                A(I,IQ,4) = TOAIQK
!               A(IP,I,1) = TOAIPR
!               A(IQ,I,1) = TOAIQR
!               A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
!               A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
!               A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IP+1,IQ-1
!                I > IP, I < IQ
!              TOAIPR = -A(I,IP,2)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,2)*CU - A(I,IQ,1)*SU
!              TOAIQR = -A(I,IP,2)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,2)*SU + A(I,IQ,1)*CU
!              TOAIPI = A(I,IP,1)*CU - A(I,IQ,2)*SU
               TOAIPI = A(IP,I,1)*CU - A(I,IQ,2)*SU
!              TOAIQI = A(I,IP,1)*SU + A(I,IQ,2)*CU
               TOAIQI = A(IP,I,1)*SU + A(I,IQ,2)*CU
!              TOAIPJ =  A(I,IP,4)*CU - A(I,IQ,3)*SU
               TOAIPJ = -A(IP,I,4)*CU - A(I,IQ,3)*SU
!              TOAIQJ =  A(I,IP,4)*SU + A(I,IQ,3)*CU
               TOAIQJ = -A(IP,I,4)*SU + A(I,IQ,3)*CU
!              TOAIPK = -A(I,IP,3)*CU - A(I,IQ,4)*SU
               TOAIPK =  A(IP,I,3)*CU - A(I,IQ,4)*SU
!              TOAIQK = -A(I,IP,3)*SU + A(I,IQ,4)*CU
               TOAIQK =  A(IP,I,3)*SU + A(I,IQ,4)*CU
!               A(I,IP,1) = TOAIPR 
                A(I,IQ,1) = TOAIQR
!               A(I,IP,2) = TOAIPI
                A(I,IQ,2) = TOAIQI
!               A(I,IP,3) = TOAIPJ
                A(I,IQ,3) = TOAIQJ
!               A(I,IP,4) = TOAIPK
                A(I,IQ,4) = TOAIQK
                A(IP,I,1) = TOAIPR
!               A(IQ,I,1) = TOAIQR
                A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IQ+1,N
!             I > IP,  I > IQ
!              TOAIPR = -A(I,IP,2)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,2)*CU - A(IQ,I,1)*SU
!              TOAIQR = -A(I,IP,2)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,2)*SU + A(IQ,I,1)*CU
!              TOAIPI = A(I,IP,1)*CU - A(I,IQ,2)*SU
               TOAIPI = A(IP,I,1)*CU + A(IQ,I,2)*SU
!              TOAIQI = A(I,IP,1)*SU + A(I,IQ,2)*CU
               TOAIQI = A(IP,I,1)*SU - A(IQ,I,2)*CU
!              TOAIPJ =  A(I,IP,4)*CU - A(I,IQ,3)*SU
               TOAIPJ = -A(IP,I,4)*CU + A(IQ,I,3)*SU
!              TOAIQJ =  A(I,IP,4)*SU + A(I,IQ,3)*CU
               TOAIQJ = -A(IP,I,4)*SU - A(IQ,I,3)*CU
!              TOAIPK = -A(I,IP,3)*CU - A(I,IQ,4)*SU
               TOAIPK =  A(IP,I,3)*CU + A(IQ,I,4)*SU
!              TOAIQK = -A(I,IP,3)*SU + A(I,IQ,4)*CU
               TOAIQK =  A(IP,I,3)*SU - A(IQ,I,4)*CU
!               A(I,IP,1) = TOAIPR 
!               A(I,IQ,1) = TOAIQR
!               A(I,IP,2) = TOAIPI
!               A(I,IQ,2) = TOAIQI
!               A(I,IP,3) = TOAIPJ
!               A(I,IQ,3) = TOAIQJ
!               A(I,IP,4) = TOAIPK
!               A(I,IQ,4) = TOAIQK
                A(IP,I,1) = TOAIPR
                A(IQ,I,1) = TOAIQR
                A(IP,I,2) = -TOAIPI
                A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
                A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
                A(IQ,I,4) = -TOAIQK
            ENDDO


!      ...  update P after "i" elimination!
            IF (DOEVE) THEN
              
               DO I=1,N

!          ... P_ip(n+1) = P_ip(n).i.c - P_iq(n).s
!              P_iq(n+1) = P_ip(n).i.s + P_iq(n).c   
!
!               ... (1,i,j,k).i  =  (i,-1,-k,j)

                 TOPIPR = -EV(I,IP,2)*CU - EV(I,IQ,1)*SU
                 TOPIQR = -EV(I,IP,2)*SU + EV(I,IQ,1)*CU

                 TOPIPI =  EV(I,IP,1)*CU - EV(I,IQ,2)*SU
                 TOPIQI =  EV(I,IP,1)*SU + EV(I,IQ,2)*CU

                 TOPIPJ =  EV(I,IP,4)*CU - EV(I,IQ,3)*SU
                 TOPIQJ =  EV(I,IP,4)*SU + EV(I,IQ,3)*CU

                 TOPIPK = -EV(I,IP,3)*CU - EV(I,IQ,4)*SU
                 TOPIQK = -EV(I,IP,3)*SU + EV(I,IQ,4)*CU

                 EV(I,IP,1) = TOPIPR
                 EV(I,IQ,1) = TOPIQR

                 EV(I,IP,2) = TOPIPI
                 EV(I,IQ,2) = TOPIQI

                 EV(I,IP,3) = TOPIPJ
                 EV(I,IQ,3) = TOPIQJ

                 EV(I,IP,4) = TOPIPK
                 EV(I,IQ,4) = TOPIQK

             ENDDO

             IF (IPRINT.GE.6) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvector after Apq_"i" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r,i,j,k:',          &
          EV(I,J,1),EV(I,J,2),EV(I,J,3),EV(I,J,4)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

         ENDIF

         ENDIF

!     =========================================
!          ... eliminate "j" part of P,Q
!     =========================================

 7        APQJ =  A(IP,IQ,3)

          IF (DABS(APQJ).GT.ZEROTHRS) THEN

            ITER = ITER + 1
            WASELIM = .TRUE.

            APPR =  A(IP,IP,1)
            AQQR =  A(IQ,IQ,1)

!         ... do calculate Cu,Su...
           BDA =  D2*APQJ/(AQQR-APPR)

           PHI = DP5*DATAN(BDA)
           CU = DCOS(PHI)
           SU = DSIN(PHI)
           CU2 = CU*CU
           SU2 = SU*SU
           CUSU = CU*SU

!           ... Update (only REAL) diagonal elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQJ 
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQJ

!...        do eliminate "j" off-diagonal IMAX,JMAX elements
!           A(IP,IQ,3)=D0
!           A(IQ,IP,3)=D0

!      ...  update off diag. Apq ...(1,i,j,k).j = (j,k,0,-i)
            XAI = A(IP,IQ,4)
            XAJ = A(IP,IQ,1)
            XAK = A(IP,IQ,2)

            A(IP,IQ,2)=-XAI
            A(IP,IQ,3)= XAJ
            A(IP,IQ,4)= XAK

!           ... the same for Aqp elements...
!           A(IQ,IP,2)= XAI
!           A(IQ,IP,3)=-XAJ
!           A(IQ,IP,4)=-XAK

!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After j-elim App,Apq,Aqp,Aqq r,i,j,k:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2),A(IP,IP,3),A(IP,IP,4)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2),A(IP,IQ,3),A(IP,IQ,4)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2),A(IQ,IP,3),A(IQ,IP,4)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2),A(IQ,IQ,3),A(IQ,IQ,4)
            ENDIF

!      ...  update Aip, Aiq (i<>p,q) after elimination!
!
! .....   A_ip(n+1) = A_ip(n).j.c - A_iq(n).s  
!         A_iq(n+1) = A_ip(n).j.s + A_iq(n).c
!             
!                 (1,i,j,k).j = (j,k,-1,-i)

            DO I=1,IP-1
!            I < IP,  I < IQ
               TOAIPR = -A(I,IP,3)*CU - A(I,IQ,1)*SU
               TOAIQR = -A(I,IP,3)*SU + A(I,IQ,1)*CU
               TOAIPI = -A(I,IP,4)*CU - A(I,IQ,2)*SU
               TOAIQI = -A(I,IP,4)*SU + A(I,IQ,2)*CU
               TOAIPJ = A(I,IP,1)*CU - A(I,IQ,3)*SU
               TOAIQJ = A(I,IP,1)*SU + A(I,IQ,3)*CU
               TOAIPK = A(I,IP,2)*CU - A(I,IQ,4)*SU
               TOAIQK = A(I,IP,2)*SU + A(I,IQ,4)*CU
                A(I,IP,1) = TOAIPR 
                A(I,IQ,1) = TOAIQR
                A(I,IP,2) = TOAIPI
                A(I,IQ,2) = TOAIQI
                A(I,IP,3) = TOAIPJ
                A(I,IQ,3) = TOAIQJ
                A(I,IP,4) = TOAIPK
                A(I,IQ,4) = TOAIQK
!               A(IP,I,1) = TOAIPR
!               A(IQ,I,1) = TOAIQR
!               A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
!               A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
!               A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IP+1,IQ-1
!               I > IP, I < IQ
!              TOAIPR = -A(I,IP,3)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,3)*CU - A(I,IQ,1)*SU
!              TOAIQR = -A(I,IP,3)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,3)*SU + A(I,IQ,1)*CU
!              TOAIPI = -A(I,IP,4)*CU - A(I,IQ,2)*SU
               TOAIPI =  A(IP,I,4)*CU - A(I,IQ,2)*SU
!              TOAIQI = -A(I,IP,4)*SU + A(I,IQ,2)*CU
               TOAIQI =  A(IP,I,4)*SU + A(I,IQ,2)*CU
!              TOAIPJ = A(I,IP,1)*CU - A(I,IQ,3)*SU
               TOAIPJ = A(IP,I,1)*CU - A(I,IQ,3)*SU
!              TOAIQJ = A(I,IP,1)*SU + A(I,IQ,3)*CU
               TOAIQJ = A(IP,I,1)*SU + A(I,IQ,3)*CU
!              TOAIPK =  A(I,IP,2)*CU - A(I,IQ,4)*SU
               TOAIPK = -A(IP,I,2)*CU - A(I,IQ,4)*SU
!              TOAIQK =  A(I,IP,2)*SU + A(I,IQ,4)*CU
               TOAIQK = -A(IP,I,2)*SU + A(I,IQ,4)*CU
!               A(I,IP,1) = TOAIPR 
                A(I,IQ,1) = TOAIQR
!               A(I,IP,2) = TOAIPI
                A(I,IQ,2) = TOAIQI
!               A(I,IP,3) = TOAIPJ
                A(I,IQ,3) = TOAIQJ
!               A(I,IP,4) = TOAIPK
                A(I,IQ,4) = TOAIQK
                A(IP,I,1) = TOAIPR
!               A(IQ,I,1) = TOAIQR
                A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IQ+1,N
!               I > IP, I > IQ
!              TOAIPR = -A(I,IP,3)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,3)*CU - A(IQ,I,1)*SU
!              TOAIQR = -A(I,IP,3)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,3)*SU + A(IQ,I,1)*CU
!              TOAIPI = -A(I,IP,4)*CU - A(I,IQ,2)*SU
               TOAIPI =  A(IP,I,4)*CU + A(IQ,I,2)*SU
!              TOAIQI = -A(I,IP,4)*SU + A(I,IQ,2)*CU
               TOAIQI =  A(IP,I,4)*SU - A(IQ,I,2)*CU
!              TOAIPJ =  A(I,IP,1)*CU - A(I,IQ,3)*SU
               TOAIPJ =  A(IP,I,1)*CU + A(IQ,I,3)*SU
!              TOAIQJ =  A(I,IP,1)*SU + A(I,IQ,3)*CU
               TOAIQJ =  A(IP,I,1)*SU - A(IQ,I,3)*CU
!              TOAIPK =  A(I,IP,2)*CU - A(I,IQ,4)*SU
               TOAIPK = -A(IP,I,2)*CU + A(IQ,I,4)*SU
!              TOAIQK =  A(I,IP,2)*SU + A(I,IQ,4)*CU
               TOAIQK = -A(IP,I,2)*SU - A(IQ,I,4)*CU
!               A(I,IP,1) = TOAIPR 
!               A(I,IQ,1) = TOAIQR
!               A(I,IP,2) = TOAIPI
!               A(I,IQ,2) = TOAIQI
!               A(I,IP,3) = TOAIPJ
!               A(I,IQ,3) = TOAIQJ
!               A(I,IP,4) = TOAIPK
!               A(I,IQ,4) = TOAIQK
                A(IP,I,1) = TOAIPR
                A(IQ,I,1) = TOAIQR
                A(IP,I,2) = -TOAIPI
                A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
                A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
                A(IQ,I,4) = -TOAIQK
            ENDDO

!      ...  update P after the A_pq"j" elimination!
            IF (DOEVE) THEN

               DO I=1,N
!          ... P_ip(n+1) = P_ip(n).j.c - P_iq(n).s
!              P_iq(n+1) = P_ip(n).j.s + P_iq(n).c 
!
!               (1,i,j,k).j =  (j, k, -1, -i)

                 TOPIPR = -EV(I,IP,3)*CU - EV(I,IQ,1)*SU
                 TOPIQR = -EV(I,IP,3)*SU + EV(I,IQ,1)*CU

                 TOPIPI = -EV(I,IP,4)*CU - EV(I,IQ,2)*SU
                 TOPIQI = -EV(I,IP,4)*SU + EV(I,IQ,2)*CU

                 TOPIPJ = EV(I,IP,1)*CU - EV(I,IQ,3)*SU
                 TOPIQJ = EV(I,IP,1)*SU + EV(I,IQ,3)*CU

                 TOPIPK = EV(I,IP,2)*CU - EV(I,IQ,4)*SU
                 TOPIQK = EV(I,IP,2)*SU + EV(I,IQ,4)*CU

                 EV(I,IP,1) = TOPIPR
                 EV(I,IQ,1) = TOPIQR

                 EV(I,IP,2) = TOPIPI
                 EV(I,IQ,2) = TOPIQI

                 EV(I,IP,3) = TOPIPJ
                 EV(I,IQ,3) = TOPIQJ

                 EV(I,IP,4) = TOPIPK
                 EV(I,IQ,4) = TOPIQK

             ENDDO

             IF (IPRINT.GE.7) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvector after Apq_"j" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r,i,j,k:',          &
          EV(I,J,1),EV(I,J,2),EV(I,J,3),EV(I,J,4)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

            ENDIF

         ENDIF

!        ... after "j" elim.check the "i" element !
         IF (DABS(A(IP,IQ,2)).GE.ZEROTHRS) GOTO 6

!     =========================================
!          ... eliminate "k" part of P,Q
!     =========================================

 8        APQK =  A(IP,IQ,4)

          IF (DABS(APQK).GE.ZEROTHRS) THEN

           ITER = ITER + 1
            WASELIM = .TRUE.

           APPR =  A(IP,IP,1)
           AQQR =  A(IQ,IQ,1)

!         ... do calculate Cu,Su...
           BDA =  D2*APQK/(AQQR-APPR)

           PHI = DP5*DATAN(BDA)
           CU = DCOS(PHI)
           SU = DSIN(PHI)
           CU2 = CU*CU
           SU2 = SU*SU
           CUSU = CU*SU

!           ... Update (only REAL) diagonal elements
            A(IP,IP,1)= APPR*CU2 + AQQR*SU2 - D2*CUSU*APQK
            A(IQ,IQ,1)= APPR*SU2 + AQQR*CU2 + D2*CUSU*APQK

!...        do eliminate the last  A_pq "k" off-diagonal  elements
!           A(IP,IQ,4)=D0
!           A(IQ,IP,4)=D0

!      ...  update off diag. Apq ...(1,i,j,k).k = (k,-j,i,-1)=(k,-j,i,0)
!               (1,i,j,k).k = (k, -j, i, -1)
            XAI = A(IP,IQ,3)
            XAJ = A(IP,IQ,2)
            XAK = A(IP,IQ,1)

            A(IP,IQ,2)=  XAI
            A(IP,IQ,3)= -XAJ
            A(IP,IQ,4)=  XAK

!           ... the same for Aqp elements...
!           A(IQ,IP,2)= -XAI
!           A(IQ,IP,3)=  XAJ
!           A(IQ,IP,4)= -XAK

!           ... check the A_pq  elements !
            IF (IPRINT.GE.6) THEN
             write(lupri,'(/A)')                                        &
    "======   After k-elim App,Apq,Aqp,Aqq r,i,j,k:"
             write(lupri,*) '  A(',IP,',',IP,'):',                      &
         A(IP,IP,1),A(IP,IP,2),A(IP,IP,3),A(IP,IP,4)
             write(lupri,*) '  A(',IP,',',IQ,'):',                      &
         A(IP,IQ,1),A(IP,IQ,2),A(IP,IQ,3),A(IP,IQ,4)
             write(lupri,*) '  A(',IQ,',',IP,'):',                      &
         A(IQ,IP,1),A(IQ,IP,2),A(IQ,IP,3),A(IQ,IP,4)
             write(lupri,*) '  A(',IQ,',',IQ,'):',                      &
         A(IQ,IQ,1),A(IQ,IQ,2),A(IQ,IQ,3),A(IQ,IQ,4)
            ENDIF

!      ...  update Aip, Aiq (i<>p,q) after elimination!
! .....   A_ip(n+1) = A_ip(n).k.c - A_iq(n).s  
!         A_iq(n+1) = A_ip(n).k.s + A_iq(n).c
!             
!                 (1,i,j,k).k = (k,-j,i,-1)

            DO I=1,IP-1
!             I < IP, I < IQ
               TOAIPR = -A(I,IP,4)*CU - A(I,IQ,1)*SU
               TOAIQR = -A(I,IP,4)*SU + A(I,IQ,1)*CU
               TOAIPI =  A(I,IP,3)*CU - A(I,IQ,2)*SU
               TOAIQI =  A(I,IP,3)*SU + A(I,IQ,2)*CU
               TOAIPJ = -A(I,IP,2)*CU - A(I,IQ,3)*SU
               TOAIQJ = -A(I,IP,2)*SU + A(I,IQ,3)*CU
               TOAIPK =  A(I,IP,1)*CU - A(I,IQ,4)*SU
               TOAIQK =  A(I,IP,1)*SU + A(I,IQ,4)*CU
                A(I,IP,1) = TOAIPR 
                A(I,IQ,1) = TOAIQR
                A(I,IP,2) = TOAIPI
                A(I,IQ,2) = TOAIQI
                A(I,IP,3) = TOAIPJ
                A(I,IQ,3) = TOAIQJ
                A(I,IP,4) = TOAIPK
                A(I,IQ,4) = TOAIQK
!               A(IP,I,1) = TOAIPR
!               A(IQ,I,1) = TOAIQR
!               A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
!               A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
!               A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IP+1,IQ-1
!               I > IP, I < IQ
!              TOAIPR = -A(I,IP,4)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,4)*CU - A(I,IQ,1)*SU
!              TOAIQR = -A(I,IP,4)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,4)*SU + A(I,IQ,1)*CU
!              TOAIPI =  A(I,IP,3)*CU - A(I,IQ,2)*SU
               TOAIPI = -A(IP,I,3)*CU - A(I,IQ,2)*SU
!              TOAIQI =  A(I,IP,3)*SU + A(I,IQ,2)*CU
               TOAIQI = -A(IP,I,3)*SU + A(I,IQ,2)*CU
!              TOAIPJ = -A(I,IP,2)*CU - A(I,IQ,3)*SU
               TOAIPJ =  A(IP,I,2)*CU - A(I,IQ,3)*SU
!              TOAIQJ = -A(I,IP,2)*SU + A(I,IQ,3)*CU
               TOAIQJ =  A(IP,I,2)*SU + A(I,IQ,3)*CU
!              TOAIPK =  A(I,IP,1)*CU - A(I,IQ,4)*SU
               TOAIPK =  A(IP,I,1)*CU - A(I,IQ,4)*SU
!              TOAIQK =  A(I,IP,1)*SU + A(I,IQ,4)*CU
               TOAIQK =  A(IP,I,1)*SU + A(I,IQ,4)*CU
!               A(I,IP,1) = TOAIPR 
                A(I,IQ,1) = TOAIQR
!               A(I,IP,2) = TOAIPI
                A(I,IQ,2) = TOAIQI
!               A(I,IP,3) = TOAIPJ
                A(I,IQ,3) = TOAIQJ
!               A(I,IP,4) = TOAIPK
                A(I,IQ,4) = TOAIQK
                A(IP,I,1) = TOAIPR
!               A(IQ,I,1) = TOAIQR
                A(IP,I,2) = -TOAIPI
!               A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
!               A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
!               A(IQ,I,4) = -TOAIQK
            ENDDO
            DO I=IQ+1,N
!              I > IP,  I > IQ
!              TOAIPR = -A(I,IP,4)*CU - A(I,IQ,1)*SU
               TOAIPR =  A(IP,I,4)*CU - A(IQ,I,1)*SU
!              TOAIQR = -A(I,IP,4)*SU + A(I,IQ,1)*CU
               TOAIQR =  A(IP,I,4)*SU + A(IQ,I,1)*CU
!              TOAIPI =  A(I,IP,3)*CU - A(I,IQ,2)*SU
               TOAIPI = -A(IP,I,3)*CU + A(IQ,I,2)*SU
!              TOAIQI =  A(I,IP,3)*SU + A(I,IQ,2)*CU
               TOAIQI = -A(IP,I,3)*SU - A(IQ,I,2)*CU
!              TOAIPJ = -A(I,IP,2)*CU - A(I,IQ,3)*SU
               TOAIPJ =  A(IP,I,2)*CU + A(IQ,I,3)*SU
!              TOAIQJ = -A(I,IP,2)*SU + A(I,IQ,3)*CU
               TOAIQJ =  A(IP,I,2)*SU - A(IQ,I,3)*CU
!              TOAIPK =  A(I,IP,1)*CU - A(I,IQ,4)*SU
               TOAIPK =  A(IP,I,1)*CU + A(IQ,I,4)*SU
!              TOAIQK =  A(I,IP,1)*SU + A(I,IQ,4)*CU
               TOAIQK =  A(IP,I,1)*SU - A(IQ,I,4)*CU
!               A(I,IP,1) = TOAIPR 
!               A(I,IQ,1) = TOAIQR
!               A(I,IP,2) = TOAIPI
!               A(I,IQ,2) = TOAIQI
!               A(I,IP,3) = TOAIPJ
!               A(I,IQ,3) = TOAIQJ
!               A(I,IP,4) = TOAIPK
!               A(I,IQ,4) = TOAIQK
                A(IP,I,1) = TOAIPR
                A(IQ,I,1) = TOAIQR
                A(IP,I,2) = -TOAIPI
                A(IQ,I,2) = -TOAIQI
                A(IP,I,3) = -TOAIPJ
                A(IQ,I,3) = -TOAIQJ
                A(IP,I,4) = -TOAIPK
                A(IQ,I,4) = -TOAIQK
            ENDDO

!      ...  update P after the A_pq"k" elimination!
            IF (DOEVE) THEN
               DO I=1,N

!          ... P_ip(n+1) = P_ip(n).k.c - P_iq(n).s
!              P_iq(n+1) = P_ip(n).k.s + P_iq(n).c 
!
!               (1,i,j,k).k = (k, -j, i, -1)

                 TOPIPR = -EV(I,IP,4)*CU - EV(I,IQ,1)*SU
                 TOPIQR = -EV(I,IP,4)*SU + EV(I,IQ,1)*CU

                 TOPIPI =  EV(I,IP,3)*CU - EV(I,IQ,2)*SU
                 TOPIQI =  EV(I,IP,3)*SU + EV(I,IQ,2)*CU

                 TOPIPJ = -EV(I,IP,2)*CU - EV(I,IQ,3)*SU
                 TOPIQJ = -EV(I,IP,2)*SU + EV(I,IQ,3)*CU

                 TOPIPK =  EV(I,IP,1)*CU - EV(I,IQ,4)*SU
                 TOPIQK =  EV(I,IP,1)*SU + EV(I,IQ,4)*CU

                 EV(I,IP,1) = TOPIPR
                 EV(I,IQ,1) = TOPIQR

                 EV(I,IP,2) = TOPIPI
                 EV(I,IQ,2) = TOPIQI

                 EV(I,IP,3) = TOPIPJ
                 EV(I,IQ,3) = TOPIQJ

                 EV(I,IP,4) = TOPIPK
                 EV(I,IQ,4) = TOPIQK

             ENDDO

             IF (IPRINT.GE.7) THEN
!              ... print out actual P eigenvector matrix
              write(lupri,'(/2X,A)')                                    &
          'Actual eigenvectors after Apq_"k" elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')_r,i,j,k:',          &
          EV(I,J,1),EV(I,J,2),EV(I,J,3),EV(I,J,4)
              ENDDO
              ENDDO
              write(lupri,*)
             ENDIF

            ENDIF

           ENDIF

          ELSE
            write(lupri,*) 'Wrong NZ!'
            CALL QUIT('QJACOBI: Wrong NZ !')
          ENDIF

! .... If wished, print out whole matrix 
           IF (IPRINT.GE.5.AND.WASELIM) THEN
             write(lupri,'(/A,2I5)')                                    &
  "Whole matrix A AFTER elimination of Apq,Aqp elements, p,q=",         &
   IP,IQ
            DO I=1,NZ
            DO J=1,NZ
             write(lupri,*)                                             &
        ('A(',I,',',J,' [',IZ,')=',A(I,J,IZ),IZ=1,NZ)
            ENDDO
            ENDDO
            write(lupri,*)
           ENDIF

           IF (IPRINT.GE.6.AND.WASELIM) THEN
!              ... print out actual P eigenvector matrix
             write(lupri,'(/2X,A)')                                     &
        'Actual eigenvectors after the Apq elimination:'
              DO I=1,N
              DO J=1,N
               write(lupri,*) '      P(',I,',',J,')=',                  &
          ( EV(I,J,IZ),IZ=1,NZ )
              ENDDO
              ENDDO
              write(lupri,*)
          ENDIF

!        ... after elimination of Apq find indexes of the next largest element
          IF (WASELIM) THEN
           CALL FINDSUMBL(A,N,NZ,NBL,IBL,IP,IQ,SUMSQ,IPRINT)
           SUMSQA(1)=SUMSQA(2)
           SUMSQA(2)=SUMSQ
          ENDIF

          ELIM(1) = ELIM(2)
          ELIM(2) = WASELIM 
          WASELIM = .FALSE.

          IF (IPRINT.GE.3) THEN
            write(lupri,'(3X,A,I6,A,D12.7)')                            &
       '...ITER=',ITER,' SUMSQ=',SUMSQ
          ENDIF
          IF (.NOT.ELIM(1).AND..NOT.ELIM(2)) FINISHED=.TRUE.

      IF (.NOT.FINISHED) GOTO 10

      IF (IPRINT.GE.2) THEN
       CALL GETTIM(CPU2,WALL2)
       CPU    = CPU2 - CPU1
       WALL   = WALL2 - WALL1
       CPUTID = SECTID(CPU)
       WALLTID = SECTID(WALL)
       write(lupri,'(2X,A,I7)')                                         &
   '>>>> QJACOBI finished, # of iterations: ',ITER
       write(lupri,'(2X,A,2D12.5,/,2X,A,D12.5)')                        &
       'last 2 sums of largets elements: ',SUMSQA(2),SUMSQA(1),         &
       '"zero" off-diag elem. threshold: ',ZEROTHRS
       WRITE(LUPRI,'(2X,A,A12,A1,A12)')                                 &
         '>>>> CPU/wall  time used in QJACOBI:',                        &
              CPUTID,'/',WALLTID
       CALL FLSHFO(LUPRI)
      ENDIF
!C *********************************************************************************

!C DO the symmetrization of the (diagonal/block-diagonal) A matrix
      DO I = 1,N-1
      DO J = I+1, N
       A(J,I,1)=A(I,J,1)
       DO IZ=2,NZ
        A(J,I,IZ)=-A(I,J,IZ)
       ENDDO
      ENDDO
      ENDDO

! ==== Finished === 
! ===  if full diagonalization, do sort eigenvalues on the main diagonal...
       IF (NBL.EQ.0) THEN
          DO I=1,N-1
            APPR = A(I,I,1)
            DO J=I+1, N
             IF (APPR.GT.A(J,J,1)) THEN
!              ... swich I-J columns and I-J rows as well
               DO K=1,N
               DO IZ=1,NZ
                 APPI=A(I,K,IZ)
                 A(I,K,IZ)=A(J,K,IZ)
                 A(J,K,IZ)=APPI
               ENDDO
               ENDDO

               DO K=1,N
               DO IZ=1,NZ
                 APPI=A(K,I,IZ)
                 A(K,I,IZ)=A(K,J,IZ)
                 A(K,J,IZ)=APPI
               ENDDO
               ENDDO
               APPR = A(I,I,1)
!            ... switch I,J columns of P(eigenvectors) !
              IF (DOEVE) THEN
               DO K=1,N
                DO IZ=1,NZ
                 APQR = EV(K,I,IZ)
                 EV(K,I,IZ) = EV(K,J,IZ)
                 EV(K,J,IZ) = APQR
                ENDDO
               ENDDO
              ENDIF

             ENDIF
            ENDDO
          ENDDO

          IF (IPRINT.GE.3) THEN
           write(lupri,'(/3X,A)')                                       &
        '>>>> Eigenvalues (after sorting) <<<<<'
           DO I=1,N
            write(lupri,*)                                              &
       '  ',I,'. >',(A(I,I,IZ),IZ=1,NZ)
           ENDDO
          ENDIF
          write(lupri,*)
      ENDIF

      IF (IPRINT.GE.3) THEN
          write(lupri,'(/3X,A)')                                        &
 '>>>> QJACOBI: RESULTING MATRIX after unit.transformations <<<<'
         DO I=1,N
         DO J=1,N
           IF (I.NE.J) THEN
          write(lupri,'(A,I4,A,I4,A,4F16.8)')                           &
      'offdiag elem, A(',I,',',J,')=',(A(I,J,IZ),IZ=1,NZ)
           ELSE
          write(lupri,'(A,I4,A,I4,A,4F16.8)')                           &
      'diagonal elem,A(',I,',',J,')=',(A(I,J,IZ),IZ=1,NZ)
           ENDIF
         ENDDO
         ENDDO
         write(lupri,*)
      ENDIF

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    .... if wished, do check: P_H*P = P*P_H = 1
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++
         IF (IPRINT.GE.2.AND.DOEVE) THEN
!         ... rewind file
          REWIND(LUF)

         IF (IPRINT.GE.4) THEN
           write(lupri,'(/3X,A)') '**** Eigenvectors:'
         ENDIF

! ***********************************************
! ***************  Real matrixes ****************
! ***********************************************
         IF (NZ.EQ.1) THEN

          DO I=1,N
          DO J=1,N
!          ....  P_H*P 
           TOPIPR = D0
!          ....  P*P_H 
           TOQIPR = D0

           DO K=1, N

!          ....  P_H*P 
            TOPIPR = TOPIPR +                                           &
                EV(K,I,1)*EV(K,J,1) 

!          ....  P*P_H 
            TOQIPR =   TOQIPR +                                         &
                EV(I,K,1)*EV(J,K,1) 

           ENDDO

!          ... do print eignevector...
          IF (IPRINT.GE.4) THEN
          write(lupri,'(A,I3,A,I3,A,F18.12)') '** P(',I,',',J,')_r:',   &
      EV(I,J,1)
          ENDIF

!         ... print of unit matrixes
          IF (IPRINT.GE.3) THEN
           write(lupri,'(A,I3,A,I3,A,F18.12,A,F18.12)')                 &
       '>> P^H*P(',I,',',J,')_r:',                                      &
      TOPIPR,' >> P*P^H:',                                              &
      TOQIPR
          ENDIF
!         ... check the unit matrix...
          IF (I.EQ.J) THEN
             IF (DABS(TOPIPR-D1).GT.ZEROTHRS.OR.                        &
            DABS(TOQIPR-D1).GT.ZEROTHRS                                 &
                   ) THEN
         write(lupri,*) 'WARNING!;  P^H*P(',I,',',J,')_r,i:',           &
      TOPIPR,' >> P*P^H:',                                              &
      TOQIPR
             ENDIF
          ELSE
             IF (   DABS(TOPIPR).GT.ZEROTHRS.OR.                        &
               DABS(TOQIPR).GT.ZEROTHRS                                 &
         ) THEN

         write(lupri,*) 'WARNING!;  P^H*P(',I,',',J,')_r:',             &
      TOPIPR,' >> P*P^H:',TOQIPR

             ENDIF
          ENDIF

          ENDDO
          ENDDO

!         ...  Do the check P.D.PH = P.(D.PH) = Sum(k) P_ik. (D.PH)_kj
!          (D.PH)_kj = Sum(l) D_kl.PH_lj = Sum(l) D_kl.(P_jl)*

          IF (IPRINT.GE.3) THEN
           WRITE(LUPRI,'(/2X,A)')                                       &
  'The P*A(diagonalized)*P^H = A(reconstructed) check:'
          ENDIF

          DO I=1,N
          DO J=1,N
            AIJ_R = D0
            DO K=1,N
             DPHKJR=D0
             DO L=1,N
!     (D.PH)_kj = Sum(l) D_kl.PH_lj = Sum(l) D_kl.(P_jl)*
               DPHKJR = DPHKJR                                          &
                  + A(K,L,1)*EV(J,L,1)

             ENDDO
               AIJ_R = AIJ_R                                            &
              + EV(I,K,1)*DPHKJR
            ENDDO

!         ... read elements from the file...
           READ(LUF,1000) I1,J1,IZ1,AIJFR
!          READ(LUF,*) I1,J1,AIJFR

           IF (IPRINT.GE.3) THEN
           write(lupri,'(A,I3,A,I3,A,D12.5)')                           &
      '>>>>>  reconstructed A(',I,',',J,')_r =',                        &
      AIJ_R
           write(lupri,'(A,I3,A,I3,A,D12.5,A,D12.5)')                   &
      '>>>>> entering elem. A(',I,',',J,')_r =',                        &
      AIJFR,' difference: ',AIJ_R-AIJFR
           ENDIF

!           ... check the differences...
           IF (DABS(AIJ_R-AIJFR).GT.ZEROTHRS1                           &
           ) THEN
           write(lupri,*)                                               &
      'WARNING ! reconstr.  A(',I,',',J,')_r =',                        &
      AIJ_R
           write(lupri,*)                                               &
      '>>>>> entering elem. A(',I,',',J,')_r =',                        &
      AIJFR
           write(lupri,*) ' difference is (in abs):',                   &
       DABS(AIJ_R-AIJFR)

           ENDIF

            ENDDO
            ENDDO

            CLOSE(LUF,STATUS="KEEP")
            write(lupri,'(2X,A,I2,A)')                                  &
       'The auxiliary file (lu=#,',LUF,') was closed, checking'//       &
       ' of A=P.Ad.P^H done'

! ****************************************************
! *************** Complex  matrixes ******************
! ****************************************************
         ELSE IF (NZ.EQ.2) THEN

          DO I=1,N
          DO J=1,N
!          ....  P_H*P 
           TOPIPR = D0
           TOPIPI = D0
!          ....  P*P_H 
           TOQIPR = D0
           TOQIPI = D0

           DO K=1, N

!          ....  P_H*P 
            TOPIPR = TOPIPR +                                           &
                EV(K,I,1)*EV(K,J,1) +                                   &
                EV(K,I,2)*EV(K,J,2) 

            TOPIPI =   TOPIPI                                           &
                - EV(K,I,2)*EV(K,J,1)                                   &
                + EV(K,I,1)*EV(K,J,2) 

!          ....  P*P_H 
            TOQIPR =   TOQIPR +                                         &
                EV(I,K,1)*EV(J,K,1) +                                   &
                EV(I,K,2)*EV(J,K,2) 

            TOQIPI =   TOQIPI                                           &
                - EV(I,K,1)*EV(J,K,2)                                   &
                + EV(I,K,2)*EV(J,K,1) 

           ENDDO

!          ... do print eignevector...
          IF (IPRINT.GE.4) THEN
          write(lupri,*) '** P(',I,',',J,')_r,i:',                      &
      EV(I,J,1),EV(I,J,2)
          ENDIF

!         ... print of unit matrixes
          IF (IPRINT.GE.4) THEN
         write(lupri,*) '>> P^H*P(',I,',',J,')_r,i:',                   &
      TOPIPR,TOPIPI,' >> P*P^H:',                                       &
      TOQIPR,TOQIPI
          ENDIF
!         ... check the unit matrix...
          IF (I.EQ.J) THEN
             IF (DABS(TOPIPR-D1).GT.ZEROTHRS.OR.                        &
               DABS(TOPIPI).GT.ZEROTHRS.OR.                             &
            DABS(TOQIPR-D1).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPI).GT.ZEROTHRS                                 &
                   ) THEN

         write(lupri,*) 'WARNING!;  P^H*P(',I,',',J,')_r,i:',           &
      TOPIPR,TOPIPI,' >> P*P^H:',                                       &
      TOQIPR,TOQIPI
             ENDIF
          ELSE
             IF (   DABS(TOPIPR).GT.ZEROTHRS.OR.                        &
               DABS(TOPIPI).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPR).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPI).GT.ZEROTHRS                                 &
         ) THEN

         write(lupri,*) 'WARNING!;  P^H*P(',I,',',J,')_r,i:',           &
      TOPIPR,TOPIPI,' >> P*P^H:',                                       &
      TOQIPR,TOQIPI

             ENDIF
          ENDIF

          ENDDO
          ENDDO

!         ...  Do the check P.D.PH = P.(D.PH) = Sum(k) P_ik. (D.PH)_kj
!          (D.PH)_kj = Sum(l) D_kl.PH_lj = Sum(l) D_kl.(P_jl)*

          IF (IPRINT.GE.3) THEN
           WRITE(LUPRI,'(/2X,A)')                                       &
  'The P*A(diagonalized)*P^H = A(reconstructed) check:'
          ENDIF

          DO I=1,N
          DO J=1,N
            AIJ_R = D0
            AIJ_I = D0
            DO K=1,N
             DPHKJR=D0
             DPHKJI=D0
             DO L=1,N

!     (D.PH)_kj = Sum(l) D_kl.PH_lj = Sum(l) D_kl.(P_jl)*

               DPHKJR = DPHKJR                                          &
                  + A(K,L,1)*EV(J,L,1)                                  &
                  + A(K,L,2)*EV(J,L,2)

               DPHKJI = DPHKJI                                          &
                  - A(K,L,1)*EV(J,L,2)                                  &
                  + A(K,L,2)*EV(J,L,1)

             ENDDO
!         ...  Do the check P.D.PH = A_ij =  P.(D.PH)_ij = Sum(k) P_ik. (D.PH)_kj
! q = (a, b, c, d), p = (x, y, z, w); 
!       = (a*x - b*y - c*z - d*w,
!           a*y + b*x + c*w - d*z,
!           a*z - b*w + c*x + d*y,
!           a*w + b*z - c*y + d*x)

               AIJ_R = AIJ_R                                            &
              + EV(I,K,1)*DPHKJR                                        &
              - EV(I,K,2)*DPHKJI

               AIJ_I = AIJ_I                                            &
              + EV(I,K,1)*DPHKJI                                        &
              + EV(I,K,2)*DPHKJR

             ENDDO

!         ... read elements from the file...
           READ(LUF,1000) I1,J1,IZ1,AIJFR
           READ(LUF,1000) I1,J1,IZ1,AIJFI
!          READ(LUF,*) I1,J1,AIJFR,AIJFI

           IF (IPRINT.GE.3) THEN
           write(lupri,*)                                               &
      '>>>>>  reconstructed A(',I,',',J,')_r,i =',                      &
      AIJ_R,AIJ_I
           write(lupri,*)                                               &
      '>>>>> entering elem. A(',I,',',J,')_r,i =',                      &
      AIJFR,AIJFI,' difference:',                                       &
        AIJ_R-AIJFR, AIJ_I-AIJFI
           ENDIF

!           ... check the differences...
           IF (DABS(AIJ_R-AIJFR).GT.ZEROTHRS1.OR.                       &
          DABS(AIJ_I-AIJFI).GT.ZEROTHRS1                                &
           ) THEN
           write(lupri,*)                                               &
      'WARNING ! reconstr.  A(',I,',',J,')_r,i =',                      &
      AIJ_R,AIJ_I
           write(lupri,*)                                               &
      '>>>>> entering elem. A(',I,',',J,')_r,i =',                      &
      AIJFR,AIJFI
           write(lupri,*) ' difference is (in abs):',                   &
       DABS(AIJ_R-AIJFR),DABS(AIJ_I-AIJFI)

           ENDIF

            ENDDO
            ENDDO
            CLOSE(LUF,STATUS="KEEP")
            write(lupri,'(/,A)')                                        &
       'The file was closed, comparison A=P.Ad.P^H was done'

! ****************************************************
! ************* Quaternion matrixes ******************
! ****************************************************
         ELSE IF (NZ.EQ.4) THEN

          DO I=1,N
          DO J=1,N
!          ....  P_H*P 
           TOPIPR = D0
           TOPIPI = D0
           TOPIPJ = D0
           TOPIPK = D0
!          ....  P*P_H 
           TOQIPR = D0
           TOQIPI = D0
           TOQIPJ = D0
           TOQIPK = D0

           DO K=1, N

!          ....  P_H*P 
            TOPIPR = TOPIPR +                                           &
                EV(K,I,1)*EV(K,J,1) +                                   &
                EV(K,I,2)*EV(K,J,2) +                                   &
                EV(K,I,3)*EV(K,J,3) +                                   &
                EV(K,I,4)*EV(K,J,4) 

            TOPIPI =   TOPIPI                                           &
                - EV(K,I,2)*EV(K,J,1)                                   &
                - EV(K,I,3)*EV(K,J,4)                                   &
                + EV(K,I,1)*EV(K,J,2)                                   &
                + EV(K,I,4)*EV(K,J,3) 

            TOPIPJ =  TOPIPJ                                            &
                + EV(K,I,1)*EV(K,J,3)                                   &
                + EV(K,I,2)*EV(K,J,4)                                   &
                - EV(K,I,3)*EV(K,J,1)                                   &
                - EV(K,I,4)*EV(K,J,2) 

            TOPIPK =   TOPIPK                                           &
                + EV(K,I,1)*EV(K,J,4)                                   &
                - EV(K,I,2)*EV(K,J,3)                                   &
                + EV(K,I,3)*EV(K,J,2)                                   &
                - EV(K,I,4)*EV(K,J,1) 

!          ....  P*P_H 
            TOQIPR =   TOQIPR +                                         &
                EV(I,K,1)*EV(J,K,1) +                                   &
                EV(I,K,2)*EV(J,K,2) +                                   &
                EV(I,K,3)*EV(J,K,3) +                                   &
                EV(I,K,4)*EV(J,K,4) 

            TOQIPI =   TOQIPI                                           &
                - EV(I,K,1)*EV(J,K,2)                                   &
                + EV(I,K,2)*EV(J,K,1)                                   &
                - EV(I,K,3)*EV(J,K,4)                                   &
                + EV(I,K,4)*EV(J,K,3) 

            TOQIPJ =   TOQIPJ                                           &
                - EV(I,K,1)*EV(J,K,3)                                   &
                + EV(I,K,2)*EV(J,K,4)                                   &
                + EV(I,K,3)*EV(J,K,1)                                   &
                - EV(I,K,4)*EV(J,K,2) 

            TOQIPK =  TOQIPK                                            &
                - EV(I,K,1)*EV(J,K,4)                                   &
                - EV(I,K,2)*EV(J,K,3)                                   &
                + EV(I,K,3)*EV(J,K,2)                                   &
                + EV(I,K,4)*EV(J,K,1) 
           ENDDO

!          ... do print eignevector...
          IF (IPRINT.GE.4) THEN
          write(lupri,*) '** P(',I,',',J,')_r,i,j,k:',                  &
      EV(I,J,1),EV(I,J,2),EV(I,J,3),EV(I,J,4)
          ENDIF

!         ... print of unit matrixes
          IF (IPRINT.GE.3) THEN
         write(lupri,*) '>> P^H*P(',I,',',J,')_r,i,j,k:',               &
      TOPIPR,TOPIPI,TOPIPJ,TOPIPK,' >> P*P^H:',                         &
      TOQIPR,TOQIPI,TOQIPJ,TOQIPK
          ENDIF
!         ... check the unit matrix...
          IF (I.EQ.J) THEN
             IF (DABS(TOPIPR-D1).GT.ZEROTHRS.OR.                        &
               DABS(TOPIPI).GT.ZEROTHRS.OR.                             &
               DABS(TOPIPJ).GT.ZEROTHRS.OR.                             &
               DABS(TOPIPK).GT.ZEROTHRS.OR.                             &
            DABS(TOQIPR-D1).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPI).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPJ).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPK).GT.ZEROTHRS) THEN

         write(lupri,*) 'WARNING!;  P^H*P(',I,',',J,')_r,i,j,k:',       &
      TOPIPR,TOPIPI,TOPIPJ,TOPIPK,' >> P*P^H:',                         &
      TOQIPR,TOQIPI,TOQIPJ,TOQIPK
             ENDIF
          ELSE
             IF (   DABS(TOPIPR).GT.ZEROTHRS.OR.                        &
               DABS(TOPIPI).GT.ZEROTHRS.OR.                             &
               DABS(TOPIPJ).GT.ZEROTHRS.OR.                             &
               DABS(TOPIPK).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPR).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPI).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPJ).GT.ZEROTHRS.OR.                             &
               DABS(TOQIPK).GT.ZEROTHRS) THEN

         write(lupri,*) 'WARNING!;  P^H*P(',I,',',J,')_r,i,j,k:',       &
      TOPIPR,TOPIPI,TOPIPJ,TOPIPK,' >> P*P^H:',                         &
      TOQIPR,TOQIPI,TOQIPJ,TOQIPK

             ENDIF
          ENDIF

          ENDDO
          ENDDO

!         ...  Do the check P.D.PH = P.(D.PH) = Sum(k) P_ik. (D.PH)_kj
!          (D.PH)_kj = Sum(l) D_kl.PH_lj = Sum(l) D_kl.(P_jl)*

          IF (IPRINT.GE.3) THEN
           WRITE(LUPRI,'(/2X,A)')                                       &
  'The P*A(diagonalized)*P^H = A(reconstructed) check:'
          ENDIF

          DO I=1,N
          DO J=1,N

            AIJ_R = D0
            AIJ_I = D0
            AIJ_J = D0
            AIJ_K = D0

            DO K=1,N

             DPHKJR=D0
             DPHKJI=D0
             DPHKJJ=D0
             DPHKJK=D0

             DO L=1,N

!     (D.PH)_kj = Sum(l) D_kl.PH_lj = Sum(l) D_kl.(P_jl)*

               DPHKJR = DPHKJR                                          &
                  + A(K,L,1)*EV(J,L,1)                                  &
                  + A(K,L,2)*EV(J,L,2)                                  &
                  + A(K,L,3)*EV(J,L,3)                                  &
                  + A(K,L,4)*EV(J,L,4)

               DPHKJI = DPHKJI                                          &
                  - A(K,L,1)*EV(J,L,2)                                  &
                  + A(K,L,2)*EV(J,L,1)                                  &
                  - A(K,L,3)*EV(J,L,4)                                  &
                  + A(K,L,4)*EV(J,L,3)

               DPHKJJ = DPHKJJ                                          &
                  - A(K,L,1)*EV(J,L,3)                                  &
                  + A(K,L,2)*EV(J,L,4)                                  &
                  + A(K,L,3)*EV(J,L,1)                                  &
                  - A(K,L,4)*EV(J,L,2)

               DPHKJK = DPHKJK                                          &
                  - A(K,L,1)*EV(J,L,4)                                  &
                  - A(K,L,2)*EV(J,L,3)                                  &
                  + A(K,L,3)*EV(J,L,2)                                  &
                  + A(K,L,4)*EV(J,L,1)
             ENDDO
!         ...  Do the check P.D.PH = A_ij =  P.(D.PH)_ij = Sum(k) P_ik. (D.PH)_kj

               AIJ_R = AIJ_R                                            &
              + EV(I,K,1)*DPHKJR                                        &
              - EV(I,K,2)*DPHKJI                                        &
              - EV(I,K,3)*DPHKJJ                                        &
              - EV(I,K,4)*DPHKJK

               AIJ_I = AIJ_I                                            &
              + EV(I,K,1)*DPHKJI                                        &
              + EV(I,K,2)*DPHKJR                                        &
              + EV(I,K,3)*DPHKJK                                        &
              - EV(I,K,4)*DPHKJJ

               AIJ_J = AIJ_J                                            &
              + EV(I,K,1)*DPHKJJ                                        &
              - EV(I,K,2)*DPHKJK                                        &
              + EV(I,K,3)*DPHKJR                                        &
              + EV(I,K,4)*DPHKJI

               AIJ_K = AIJ_K                                            &
              + EV(I,K,1)*DPHKJK                                        &
              + EV(I,K,2)*DPHKJJ                                        &
              - EV(I,K,3)*DPHKJI                                        &
              + EV(I,K,4)*DPHKJR

             ENDDO

!         ... read elements from the file...
           READ(LUF,1000) I1,J1,IZ1,AIJFR
           READ(LUF,1000) I1,J1,IZ1,AIJFI
           READ(LUF,1000) I1,J1,IZ1,AIJFJ
           READ(LUF,1000) I1,J1,IZ1,AIJFK
!          READ(LUF,*) I1,J1,AIJFR,AIJFI,AIJFJ,AIJFK

           IF (IPRINT.GE.3) THEN
           write(lupri,*)                                               &
      '>>>>>  reconstructed A(',I,',',J,')_r,i,j,k =',                  &
      AIJ_R,AIJ_I,AIJ_J,AIJ_K
           write(lupri,*)                                               &
      '>>>>> entering elem. A(',I,',',J,')_r,i,j,k =',                  &
      AIJFR,AIJFI,AIJFJ,AIJFK,' difference:',                           &
      AIJ_R-AIJFR,AIJ_I-AIJFI,                                          &
      AIJ_J-AIJFJ,AIJ_K-AIJFK
           ENDIF

!           ... check the differences...
           IF (DABS(AIJ_R-AIJFR).GT.ZEROTHRS1.OR.                       &
          DABS(AIJ_I-AIJFI).GT.ZEROTHRS1.OR.                            &
          DABS(AIJ_J-AIJFJ).GT.ZEROTHRS1.OR.                            &
          DABS(AIJ_K-AIJFK).GT.ZEROTHRS1) THEN
           write(lupri,*)                                               &
      'WARNING ! reconstr.  A(',I,',',J,')_r,i,j,k =',                  &
      AIJ_R,AIJ_I,AIJ_J,AIJ_K
           write(lupri,*)                                               &
      '>>>>> entering elem. A(',I,',',J,')_r,i,j,k =',                  &
      AIJFR,AIJFI,AIJFJ,AIJFK
           write(lupri,*) ' difference is (in abs):',                   &
       DABS(AIJ_R-AIJFR),DABS(AIJ_I-AIJFI),                             &
       DABS(AIJ_J-AIJFJ),DABS(AIJ_K-AIJFK)

           ENDIF

            ENDDO
            ENDDO
            CLOSE(LUF,STATUS="KEEP")
            write(lupri,'(/,A/)')                                       &
       'The file was closed, comparison A=P.Ad.P^H was done'
         ENDIF

         ENDIF

         RETURN
         END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
!  /* Deck FINDSUMBL  */
         SUBROUTINE FINDSUMBL(A,N,NZ,NBL,IBL,IP,IQ,SUMSQ,IPRINT)
!****************************************************************************
!
!  ... find the largest off diagonal element (in block) in the upper triangular
!       and the sum of absolute values of its elements (r,i,j,k)
!
!     On output:  IP < IQ
!                 SUMSQ
!
!     Called from QJACOBI
!
!     Written by Miro Ilias, Strasbourg, 2006
!
!****************************************************************************
#include"implicit.h"
#include"priunit.h"
         DIMENSION A(N,N,NZ)
         DIMENSION IBL(NBL,1:4)

         SUMSQ=0.0D0
         EMAX = 0.0D0
         IF (NBL.EQ.0) THEN
!         ... no blocks... run over all off-diagonal elements
!         ... for finding largest elem, do make the initialization of 1,2 elem...
           DO IZ=1,NZ
            EMAX = EMAX + DABS(A(1,2,IZ))
           ENDDO
           IMAX = 1
           JMAX = 2

          DO I=1,N-1
          DO J=I+1,N

          EMAX1 = 0.D0
          DO IZ=1,NZ
           EMAX1 = EMAX1 + DABS(A(I,J,IZ))
          ENDDO

          IF (EMAX1.GT.EMAX) THEN
            IMAX = I
            JMAX = J
            EMAX = EMAX1
          ENDIF

          ENDDO
          ENDDO
         ELSE IF (NBL.GE.1) THEN
!         ... for finding largest elem, do make the initialization elem...
           DO IZ=1,NZ
            EMAX = EMAX + DABS(A(IBL(1,1),IBL(1,2),IZ))
           ENDDO
           IMAX = IBL(1,1)
           JMAX = IBL(1,2)
!          ... run only over specified blocks !!
          DO IB=1,NBL
           DO I=IBL(IB,1),IBL(IB,3)
           DO J=IBL(IB,2),IBL(IB,4)
            EMAX1 = 0.D0
            DO IZ=1, NZ
             EMAX1 = EMAX1 + DABS(A(I,J,IZ))
            ENDDO
            IF (EMAX1.GT.EMAX) THEN
              IMAX = I
              JMAX = J
              EMAX = EMAX1
            ENDIF
           ENDDO
           ENDDO
          ENDDO
         ELSE
          CALL QUIT('FINDSUMBL: Bad NBLOCKS in QJACOBI!')
         ENDIF

         IF (IMAX.EQ.JMAX) THEN 
          IF (IPRINT.GE.7) THEN
            write(lupri,*) 'FINDSUMBL: Find IMAX=',IMAX,' JMAX=',JMAX
            write(lupri,*) 'FINDSUMBL: IBL(1:',NBL,';1:4) elements:'
            do i=1,NBL
             write(lupri,*)                                             &
       '>>> IBL(',i,';1)=',IBL(i,1),' IBL(',i,';2)=',IBL(i,2),          &
       ' IBL(',i,';3)=',IBL(i,3),' IBL(',i,';4)=',IBL(i,4)
            enddo 
           ENDIF
!           CALL QUIT('Error in FINDSUMBL: IMAX=JMAX ')
         ELSE IF(IMAX.GT.JMAX) THEN
            IMAX1 = IMAX
            IMAX = JMAX
            JMAX = IMAX1
         ENDIF

         IP = IMAX
         IQ = JMAX

         SUMSQ = EMAX

          IF (IPRINT.GE.10) THEN
            write(lupri,*)                                              &
     'FINDSUMBL: the whole investigated matrix'
            DO I=1,N
            DO J=1,I
            DO IZ=1,NZ
             IF (I.EQ.J) THEN
             write(lupri,*) '(',I,';',J,';',IZ,') :',A(I,J,IZ)
             ELSE
             write(lupri,*) '(',I,';',J,';',IZ,') :',A(I,J,IZ),A(J,I,IZ)
             ENDIF
            ENDDO
            ENDDO
            ENDDO
          ENDIF

         RETURN
         END

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&*
!  /* Deck testblocks  */
      SUBROUTINE TESTBLOCKS(N,NZ,NBL,IBL,IPRINT)
!****************************************************************************
!
!   Check specifix blocks of hermitian matrix to be diagonalized.
!
!   Called from QJACOBI
!
!   Written by Miro Ilias, Strasbourg, 2006
!
!****************************************************************************
#include"implicit.h"
#include"priunit.h"
      DIMENSION IBL(NBL,1:4)
      LOGICAL BAD
      DATA ZEROTHRS / 1.0D-11 /

      BAD = .FALSE.

!       ... print info on blocks first
      IF (IPRINT.GE.2) THEN
        write(lupri,*) 'QJACOBI(TESTBLOCKS): # of subblocks=',NBL
        IF (NBL.GT.0) THEN
         DO I=1, NBL
          write(lupri,*) I,'. block; left upper point',                 &
     IBL(I,1),IBL(I,2),' right bottom point:',IBL(I,3),IBL(I,4)
         ENDDO    
        ENDIF
        write(*,*)
      ENDIF
!     ... do own check
      DO I=1,NBL
       I1=IBL(I,1)
       J1=IBL(I,2)
       I2=IBL(I,3)
       J2=IBL(I,4)

       IF (I1.GT.I2.OR.J1.GT.J2) THEN
        write(lupri,*) 'I1,I2:',I1,I2,' J1,J2:',J1,J2,' * N=',N
        write(lupri,*) 'I1>I2 or J1>J2!'
        BAD=.TRUE.
       ENDIF

       IF (I1.LT.1.OR.I1.GT.N-1.OR.J1.LT.1.OR.J1.GT.N                   &
  .OR.I2.LT.1.OR.I2.GT.N-1.OR.J2.LT.1.OR.J2.GT.N) THEN
        write(lupri,*) 'I1.LT.1.OR.I1.GT.N-1.OR.J1.LT.1.OR.J1.GT.N!'
        write(lupri,*) 'OR.I2.LT.1.OR.I2.GT.N-1.OR.J2.LT.1.OR.J2.GT.N'
        BAD=.TRUE.
       ENDIF

       IF (.NOT.(I1.LT.J1.AND.I2.LT.J2)) THEN
        write(lupri,*) 'NOT(I1<I2 and I2<J2!'
        BAD=.TRUE.
       ENDIF

      ENDDO

      INUMEL = 0
      IF (BAD) THEN
       CALL QUIT('BAD BLOCKS !')
      ELSE
       IF (IPRINT.GE.2) THEN
        write(lupri,'(/2X,A)')                                          &
 'QJACOBI(TESTBLOCKS): The off-elements (indexes) to be'//              &
 ' diagonalized (zeroed) are OK:'
       ENDIF
        IF (IPRINT.GE.3) THEN
        DO I=1,NBL
        DO I1=IBL(I,1),IBL(I,3)
        DO J1=IBL(I,2),IBL(I,4)
          write(lupri,*) I1,J1,' <> ',J1,I1
          INUMEL = INUMEL + 1
        ENDDO
        ENDDO
        ENDDO
        ENDIF
      ENDIF

      IF (IPRINT.GE.3) THEN
        write(lupri,*)  & 
  'Number of elements (upper triangular) to be eliminated:',INUMEL    
        write(lupri,*)
      ENDIF

      RETURN
      END
!#endif /* MOD_UNRELEASED */

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&& 
module diagmod
!-------------------------------------------------------------------------------------------
!  Module for the standalone testing program of diagonalization routines,
! utils/diag.F
!
! Written by Miro Ilias, August 2014
!-------------------------------------------------------------------------------------------
  integer :: N=0,NZ=0, print_level
  integer :: lu=12
  real*8, allocatable :: Ar(:,:), ASr(:,:),Ai(:,:),ASi(:,:),  & 
                         Aj(:,:),ASj(:,:), Ak(:,:), ASk(:,:)
  character*21 :: inpfile="DIAGONALIZE_TESTS.INP"
  character*50 :: title_text
  character*50 :: matrix_file_name

  logical :: do_dsyevr = .false.
  logical :: do_dsyevr_ulf = .false.
  logical :: do_rs = .false.
  logical :: do_rsjaco = .false.
  logical :: do_qjaco = .false.
  logical :: do_householder_psb = .false.
  logical :: do_eigv_check = .false.
  
contains

  subroutine read_matrix
  ! allocates space and reads the matrix for the diagonalization testing from the file
  use include_dgroup_h, only: get_IPQTOQ
  implicit none
#include "priunit.h"
   integer :: i,j, ii,jj
   logical :: file_exists

   inquire(file = matrix_file_name, exist = file_exists )

   if (.not.file_exists) then
     write(LUPRI,"(2x,a,a)") "The matrix input file does not exist ! matrix_file_name=",trim(matrix_file_name)
     call quit("The matrix input file does not exist !")
   endif

   open (lu,file=trim(matrix_file_name), access="sequential", status="old")  
   read(lu,*) N,NZ
   write(LUPRI,"(a,i3,a,i1)") 'Matrix dimension, read N=',N,' and NZ=',NZ
  ! allocates space for matrices
   call allocate_space
  ! read matrix elements from the open file
   do i=1,N
   do j=1,N
     if (NZ==1) then
       read(lu,*) ii,jj,Ar(ii,jj)
       ASr(ii,jj)=Ar(ii,jj)
     else if (NZ==2) then
       read(lu,*) ii,jj,Ar(ii,jj),Ai(ii,jj)
       ASr(ii,jj)=Ar(ii,jj); ASi(ii,jj)=Ai(ii,jj)
     else if (NZ==4) then
       read(lu,*) ii,jj,Ar(ii,jj),Ai(ii,jj),Aj(ii,jj),Ak(ii,jj)
       ASr(ii,jj)=Ar(ii,jj); ASi(ii,jj)=Ai(ii,jj)
       ASj(ii,jj)=Aj(ii,jj); ASk(ii,jj)=Ak(ii,jj)
     endif
    enddo
    enddo
    close(lu,status="keep")
    write(LUPRI,"(2x,a)") "The matrix for diagonalization was read."
    if (print_level >=2) then
      if (NZ==1) then
      CALL PRQMAT(Ar,N,N,N,N,NZ,get_IPQTOQ(1,0),LUPRI) 
      endif
    endif
    ! print out read matrix
  end subroutine read_matrix

  subroutine allocate_space
    if (N<=1.or.NZ<=0) then 
       print *,'allocate_space: N=',N,' NZ=',NZ
       call quit('wrong allocation parameters N or NZ !')
    endif
    if (NZ==1) then
      allocate( Ar(N,N), stat=ierr )
      call alloc_stat_check(" Ar(N,N)",ierr)
      allocate( ASr(N,N), stat=ierr )
      call alloc_stat_check(" ASr(N,N)",ierr)
    else if (NZ==2) then
      allocate( Ar(N,N), Ai(N,N), stat=ierr )
      call alloc_stat_check(" Ar(N,N), Ai(N,N)",ierr)
      allocate( ASr(N,N), ASi(N,N), stat=ierr )
      call alloc_stat_check(" ASr(N,N), ASi(N,N)",ierr)
    else if (NZ==4) then
      allocate( Ar(N,N), Ai(N,N), Aj(N,N), Ak(N,N), stat=ierr )
      call alloc_stat_check(" Ar(N,N), Ai(N,N), Aj(N,N), Ak(N,N)",ierr)
      allocate( ASr(N,N), ASi(N,N), ASj(N,N), ASk(N,N), stat=ierr )
      call alloc_stat_check(" ASr(N,N), ASi(N,N), ASj(N,N), ASk(N,N)",ierr)
    else 
       print *,'wrong value of NZ=',NZ
       call quit('NZ must be either 1 or 2 or 4')
    endif
  end subroutine allocate_space

  subroutine generate_random_matrix
  ! allocates and generates random Hemitian matrix instead of reading it from a file
  end subroutine generate_random_matrix

  subroutine test_lapack_dsyevr_ulf
!-------------------------------------------------------------------------------------------
!    Test program for lapack symmetric eigensolver DSYEVR
!     By ulfek 2008, transferred by Miro Ilias, Aug 2014
!-------------------------------------------------------------------------------------------
   implicit none
#include "priunit.h"
   integer :: i, j, k, N2, N3, lwork, liwork, info, idum, M,  & 
           idum2,idum3,idum4
!     Store matrices in contigous memory blocks with (:,:) ?
   double precision, allocatable :: A(:,:),vec(:,:),eig(:)
   double precision, allocatable :: work(:)
   double precision ddum, workdum, maxerr, maxvecerr, pi, verr, ddum2
   integer, allocatable :: iwork(:), isuppz(:)
   character(len=32) :: arg
   logical :: issorted

      pi = 2*acos(0.0D0)
      write(LUPRI,"(/,2X,A,I2,A)") 'Testing DSYEVR diagonalization of',N,'**2 matrix (Ulf)'

#ifdef FORTRAN2003
      print *,'Using',sizeof(N)*8,' bit integers'
      print *,'Allocating matrix and vectors ',
     &     ((2*N**2 + 3*N)*sizeof(ddum))/2**10,' kB'
#endif
      call flush(LUPRI)
      allocate( A(N,N), vec(N,N), eig(N), isuppz(2*N) )
!     set A = 0.5 on the upper and lower second diagonals
      do i=1,N
         do j=1,N
            A(i,j) = 0.0d0
         enddo
      enddo
      do i=1,N-1
         A(i,i+1) = 0.5d0
         A(i+1,i) = 0.5d0
      enddo
 if (print_level >= 3) then
      print *,'Matrix to diagonalize:'
      do i=1,N
         write (LUPRI,'(100F4.1)') (A(i,j),j=1,N)
      enddo
 endif
      N2 = N
      N3 = N
      CALL DSYEVR('V','A','U',N,A,N2,DDUM,DDUM2,IDUM,IDUM2,0.0D0,  & 
           IDUM3,EIG,VEC,N3,IDUM4,WORKDUM,-1,liwork,-1,INFO)
      if (info.ne.0) then
         print *, 'DSYEVR memory query failed:',info
         call quit('DSYEVR memory query failed!')
      endif
      lwork = workdum
#ifdef FORTRAN2003
      print *,'Allocating work buffers',(lwork*sizeof(ddum)+ & 
            liwork*sizeof(idum))/2**10,' kB'
#endif
      call flush(LUPRI)
      allocate( work(lwork), iwork(liwork) )
      print *,'Diagonalizing..'
      call flush(LUPRI)
      CALL DSYEVR('V','A','U',N,A,N2,DDUM,DDUM2,IDUM,IDUM2,0.0D0,  & 
           M,EIG,VEC,N3,isuppz,work,lwork, iwork,liwork,info)
      if (info.ne.0.or.M.ne.N) then
         print *,'Diagonalization failed, info = ',info,' M =',M
         call quit('Diagonalization failed!')
      endif
      issorted = .true.
      do i=1,N-1
         if (eig(i).gt.eig(i+1)) then
            issorted = .false.
         endif
      enddo
      print *,'Eigenvalues sorted ascending:',issorted
      if (issorted) then
!     Check eigenvalues and vectors against analytical solution
!     Currently only implemented for ascending sorting
         maxerr = 0
         maxvecerr = 0
         do i=1,N
            j = (N-i+1)
            maxerr = max(maxerr,abs(eig(i) - cos((j*pi)/(N+1))))
            if (print_level >= 2) then
              print *,'vector ',i,' eigval:',eig(i)
              print *,'  <Calculated>             <Exact>'
            endif
            do k=1,N
!     Discard the phase of vectors with abs
               verr = abs(eig(i)*(abs(vec(k,i)) -  & 
                    abs(sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1)))))
               maxvecerr = max(maxvecerr, verr)
            if (print_level >= 3) then
               if (verr.gt.N*1.0D-14) then
                  print *,vec(k,i),sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1))   & 
                       ,'<<<<< WARNING'
               else
                  print *,vec(k,i),sqrt(2.0D0/(N+1))*sin((j*k*pi)/(N+1))
               endif
            endif
            enddo
        if (print_level >= 4) then
            print *,'maxvecerr so far:',maxvecerr
        endif
         enddo
         print *,'Maximum absolute eigenvalue error:',maxerr
         print *,'Maximum absolute vector element error:',maxvecerr
      else
         print *,'Eigenvalues not sorted so not testing them,'//  & 
              ' this can be ok.'
      endif
  end subroutine test_lapack_dsyevr_ulf

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
 subroutine test_lapack_dsyevr
 !
 ! tests lapack's routine DSYERV for real only symmetric matrices
 !
   use machine_parameters
   implicit none
#include "priunit.h"
   integer :: ierr, MATZ, LLWORK, KLWORK, LILWORK
   integer :: i
   integer :: IL, INFO, IU, LDA, LDZ, LIWORK, LWORK, M, LRA, LRV
   character*1 :: JOBZ, RANGE, UPLO
   real*8 :: ABSTOL, VL, VU
   integer, allocatable :: ISUPPZ(:), IWORK(:)
   real*8, allocatable ::  WORK(:), EIG(:),VEC(:,:)
   real*8 :: DDUM, WORKDUM
   integer :: IWORKDUM,IDUM
   real*8, parameter :: D0=0.0D0

   call get_machine_parameters
   ABSTOL=mp_sfmin; JOBZ='V'; RANGE = 'A'; UPLO = 'U'
   DDUM = D0; IDUM = 0; WORKDUM = D0; IWORKDUM = 0

   print *,'lapack_dsyerv ABSTOL=',ABSTOL

   MATZ = 1 
   IF(MATZ.EQ.1) THEN
        JOBZ='V'
   ELSEIF(MATZ.EQ.0) THEN
        JOBZ='N'
   ELSE
     WRITE(LUPRI,*) 'QDIAG: Illegal value of MATZ: ',MATZ
     CALL QUIT('QDIAG: Illegal value of MATZ !')
   ENDIF

   allocate(EIG(N),VEC(N,N), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv EIG(N),VEC(N,N)",ierr)
    
   ! 1.call to determine size 
   IERR =  0 
   ! FYI: LLWORK=-1; LIWORK=-1 on input means calculate temporary storage
   M = N ; LRA = N ;LRV = N
   CALL DSYEVR(JOBZ,'A','U',N,Ar,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  & 
               M,EIG,VEC,LRV,IDUM,EIG,-1,IWORKDUM,-1,IERR)

   LLWORK = NINT(EIG(1))
   LILWORK = IWORKDUM 

   allocate(WORK(LLWORK), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv WORK(KLWORK)",ierr)

   allocate(IWORK(LILWORK), stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv IWORK(LILWORK)",ierr)

   allocate(ISUPPZ(2*N),stat=ierr)
   call alloc_stat_check("test_lapack_dsyerv ISUPPZ",ierr)

   IERR = 0
   CALL DSYEVR(JOBZ,'A','U',N,Ar,LRA,DDUM,DDUM,IDUM,IDUM,ABSTOL,  &
               M,EIG,VEC,LRV,ISUPPZ,WORK,LLWORK,IWORK,LILWORK,IERR)     
   if (IERR /= 0) then
     print *,'The lapack_dsyerv routine ended with error ! ierr=',ierr
   endif

  ! ... print out eigenvalues
  if (print_level >= 1) then
     print *,"LAPACK DSYEVR eigenvalues:"
     do i=1,N
       print *,i,EIG(i)
     enddo
   endif

   if (do_eigv_check) call eigv_check(EIG,VEC,'**** LAPACK DSYEVR ****')
   deallocate(WORK,IWORK,EIG,VEC,ISUPPZ,stat=ierr)
  end subroutine test_lapack_dsyevr

  subroutine test_dirac_jacobi
  ! tests DIRAC's RSJACO for real only symmetric matrices
  end subroutine test_dirac_jacobi

  subroutine test_dirac_householder
  ! tests DIRAC's Houeholder for real, complex and quaternion Hermitian matrices
  end subroutine test_dirac_householder

  subroutine test_PSB_householder
  ! tests Paul Bagus own Houeholder routine for real symmetric matrices
  end subroutine test_PSB_householder


!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
  subroutine test_dirac_rs
! tests DIRAC's EISPACK routine RS real Hermitian matrix
#include "priunit.h"
    real*8, allocatable :: eigval(:),eigvect(:,:),temp1(:),temp2(:)
    integer :: ierr,i
! ... allocate arrays for the subroutine
    allocate(eigval(N),eigvect(N,N),temp1(N),temp2(N),stat=ierr)
    call alloc_stat_check("eigval(N),eigvect(N,N),temp1(N),temp2(N) ",ierr)
    print *,"arrays for RS subroutine allocated"
    matz=1 ! have both eigenvalues and eigenvectors
    call RS(N,N,Ar,eigval,matz,eigvect,temp1,temp2,ierr)
    if (ierr /= 0) then
      print *,'The RS routine ended with error !'
      !call quit('RS routine ended with error !')
    endif
! ... print out eigenvalues
    if (print_level >= 1) then
      print *,"Eispack's RS eigenvalues:"
      do i=1,N
        print *,i,eigval(i)
      enddo
    endif
    if (do_eigv_check) call eigv_check(eigval,eigvect,'**** EISPACK RS ****')
    deallocate(eigval,eigvect,temp1,temp2,stat=ierr)
  end subroutine test_dirac_rs

  subroutine test_dirac_ch
! tests DIRAC's EISPACK routine CH for complex Hermitian matrix
  end subroutine test_dirac_ch

  subroutine eigv_check(eigval,eigvect,routine_name)
  !-----------------------------------------------------------------------------------
  !
  ! Do the diagonalization checks using original Hermitian matrix, eigenvalues and
  ! eigenvectors
  ! 
  !           eigvect^+ * AS * eigvect - eigval  = 0
  !
  !           eigvect^+ * eigvect = I
  !
  !           eigvect * eigvect^+ = I
  !
  ! where ASr(:,:),ASi(:,:), ....
  ! 
  !-----------------------------------------------------------------------------------
   implicit none
#include "priunit.h"
#include "dgroup.h"
    real*8, intent(in) :: eigval(N),eigvect(N,N)
    character(*), intent(in) :: routine_name
    real*8, allocatable :: AXr(:,:),AYr(:,:),C(:,:)
    integer :: i,j, ierr, iprint
    real*8 :: C_norm_diag, C_norm_offdiag
    real*8 :: AXr_norm_diag, AXr_norm_offdiag
    real*8 :: AYr_norm_diag, AYr_norm_offdiag
    integer(kind=8) :: maxsize = -1
    
    allocate( AXr(N,N),AYr(N,N),C(N,N),stat=ierr )
    call alloc_stat_check("alloc AXr(N,N),AYr(N,N),C(N,N)",ierr)
 
   ! write(LUPRI,*) 'Eigenvectors matrix:'
   ! call prqmat(eigvect,N,N,N,N,1,IPQTOQ(1,0),LUPRI)

   ! Do  eigvect^{+)*ASr = AXr
    call dgemm('T','N',N,N,N,1.0d0,eigvect,N,ASr,N,0.0d0,AXr,N)

   ! Do  AXr * eigvect = C
    call dgemm('N','N',N,N,N,1.0d0,AXr,N,eigvect,N,0.0d0,C,N)

   ! Another check:  eigvect^{+)*eigvect = I = AXr
    call dgemm('T','N',N,N,N,1.0d0,eigvect,N,eigvect,N,0.0d0,AXr,N)
    if (print_level>=2) then
      write(LUPRI,*) "Printout of eigvect^{+)*eigvect =? I "
      call prqmat(AXr,N,N,N,N,1,IPQTOQ(1,0),LUPRI)
    endif

   ! Another check:  eigvect*eigvect^{+} = I = AYr
    call dgemm('N','T',N,N,N,1.0d0,eigvect,N,eigvect,N,0.0d0,AYr,N)
    if (print_level>=2) then
      write(LUPRI,*) "Printout of eigvect*eigvect^{+} =? I "
      call prqmat(AYr,N,N,N,N,1,IPQTOQ(1,0),LUPRI)
    endif

    ! substract eigenvalues ...
    C_norm_diag = 0.0d0; AXr_norm_diag = 0.0d0; AYr_norm_diag = 0.0d0
    do i=1,N
      C(i,i) = C(i,i) - eigval(i); 
      AXr(i,i)=AXr(i,i)-1.0d0; AYr(i,i)=AYr(i,i)-1.0d0
      C_norm_diag = C_norm_diag + DABS(C(i,i))
      AXr_norm_diag = AXr_norm_diag  + DABS(AXr(i,i))
      AYr_norm_diag = AYr_norm_diag  + DABS(AYr(i,i))
    enddo
    C_norm_diag   = C_norm_diag/DFLOAT(N)
    AXr_norm_diag = AXr_norm_diag/DFLOAT(N)
    AYr_norm_diag = AYr_norm_diag/DFLOAT(N)

    C_norm_offdiag = 0.0d0; AXr_norm_offdiag = 0.0d0;  AYr_norm_offdiag = 0.0d0
    do i = 1, N
    do j = 1, N
     if (i /= j) then
       C_norm_offdiag = C_norm_offdiag + DABS(C(i,j))
       AXr_norm_offdiag = AXr_norm_offdiag + DABS(AXr(i,j))
       AYr_norm_offdiag = AYr_norm_offdiag + DABS(AYr(i,j))
     endif
    enddo
    enddo
    C_norm_offdiag = C_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))
    AXr_norm_offdiag = AXr_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))
    AYr_norm_offdiag = AYr_norm_offdiag/((DFLOAT(N)*DFLOAT(N))-DFLOAT(N))

    ! ... Good printout
    write(LUPRI,"(2X,A)") routine_name
    write(LUPRI,"(1X,A,D10.4,1X,A,D10.4)") "U^{+}*A*U - eps ?= 0> norm/diag:",C_norm_diag,  "norm/offdiag:",C_norm_offdiag
    write(LUPRI,"(1X,A,D10.4,1X,A,D10.4)") "    U^{+}*U - I ?= 0> norm/diag:",AXr_norm_diag,"norm/offdiag:",AXr_norm_offdiag
    write(LUPRI,"(1X,A,D10.4,1X,A,D10.4)") "    U*U^{+} - I ?= 0> norm/diag:",AYr_norm_diag,"norm/offdiag:",AYr_norm_offdiag

    deallocate( AXr, AYr, C, stat=ierr )
    call alloc_stat_check("dealloc AXr(N,N), AYr(N,N) C(N,N)",ierr)

  end subroutine eigv_check

!> main testing routine of the QJACOBi quaternion diagonalization method
subroutine QJACOBITEST
use qjacobi_mod,only:A,EV,IPRINT,DOEVE,NBL,IBL
implicit none
#include"priunit.h"
integer :: istat,i,j
  write(LUPRI,"(/,2X,A)") "*** Quaternion Jacobi diagonalization test ***"
  ! ... get N,NZ from the file
  !CALL GET_N_NZ
  !> allocate matrix to be diagonalized, N,NZ,print_level from the "module diagmod"
  allocate(A(N,N,NZ),stat=istat)
  allocate(EV(N,N,NZ),stat=istat) !> space for eigenvectors (if desired)
  !> allocate space for blocks positions in the matrix to be diagonalized
  allocate(IBL(NBL,1:4),stat=istat)

  ! ... fill the hermitian matrix A (real, quaternion) matrix
  do i=1,N
  do j=1,N
   if      (NZ==1) then
     A(i,j,1) = ASr(i,j)
   else if (NZ==2) then
     A(i,j,1) = ASr(i,j)
     A(i,j,2) = ASi(i,j)
   else if (NZ==4) then
     A(i,j,1) = ASr(i,j)
     A(i,j,2) = ASi(i,j)
     A(i,j,3) = ASj(i,j)
     A(i,j,4) = ASk(i,j)
   else
    call quit('wrong NZ !')
   endif
  enddo
  enddo
  !CALL PREPMATRIX
  ! ... call of the external diagonalization routine ...
  IPRINT = print_level+1; NBL = 0 ; DOEVE = .true.
#ifdef MOD_UNRELEASED
  CALL QJACOBI(A,EV,N,NZ,NBL,IBL,DOEVE,IPRINT)
#else
  CALL QUIT('QJACOI is not in this DIRAC release !')
#endif


END subroutine

end module diagmod

#endif /* MOD_UNRELEASED */

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
Subroutine alloc_stat_check(alloc_label,istat)
!---------------------------------------------------------------------------
! very simple routine for checking status of dynamic allocation/deallocation 
! written by Miro Ilias, July 2014
!---------------------------------------------------------------------------
#include "priunit.h"
  character(*), intent(in) :: alloc_label
  integer, intent(in) :: istat
  if (istat /= 0) then
    write(LUPRI,*) 'Error in (de)allocation, variable :',alloc_label
    write(LUPRI,*) 'istat=',istat
    call quit('error in F90 dynamic allocation !')
  endif
End Subroutine alloc_stat_check
