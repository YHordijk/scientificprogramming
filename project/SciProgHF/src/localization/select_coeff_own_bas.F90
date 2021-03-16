module select_coeff_own_bas
implicit none
!================================================================================!
!!! This is not a true module. It depends on external subroutines and .h files !!!
!                                                                                !
!                This module is used in func_Pipek-Mezey module                  !
!          it will select MO coefficient for the projection analysis             !
!================================================================================!

public :: selown_drv

private :: SELOWN_mod, READT_MOD, READI_MOD, SELOWC_mod, SELOWE_mod, SELOWB_mod, SELOWI_mod

contains

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  subroutine selown_drv(iopt,cmo,eig,ibeig,inuc,nnuc,iunit,cmofil, &
                        csel,esel,ibsel,nvecs,ioff,kvec,           &
                        nstr)
  !------------------------------------------------------------------------------!
  implicit none

  ! external variables
  character*6, intent(in) :: cmofil
  integer, intent(in) :: iopt, inuc, nnuc, iunit

  integer, intent(in), dimension(  :  ) :: kvec
  integer, intent(in), dimension(2,0:2) :: nstr
  integer, intent(in), dimension(  2  ) :: nvecs, ioff
  integer, intent(in), dimension(  *  ) :: ibeig, ibsel
  real*8,  intent(in), dimension(  *  ) :: eig, cmo, csel, esel

  ! internal variables

  !------------------------------------------------------------------------------!
  call selown_mod(iopt,cmo,eig,ibeig,inuc,nnuc,iunit,cmofil, &
                  csel,esel,ibsel,nvecs,ioff,kvec,nstr)

  !------------------------------------------------------------------------------!

  end subroutine selown_drv

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE SELOWN_mod(IOPT,CMO,EIG,IBEIG,INUC,NNUC,IUNIT,CMOFIL, &
                    CSEL,ESEL,IBSEL,NVECS,IOFF,KVEC,NSTR)
!***********************************************************************
!
!     Select coeffiecients for symmetry independent nucleus INUC
!     calculated in its own basis
!     Read option is provided from bit-packed IOPT:
!       0001 - give restart info
!       0010 - read coefficients
!       0100 - read eigenvalues
!       1000 - read boson irrep identification
!
!     The array IDIM contains the following information:
!       IDIM(1,IFRP) = NPSH(IFRP)    : number of positronic solutions
!       IDIM(2,IFRP) = NESH(IFRP)    : number of electronic solutions
!       IDIM(3,IFRP) = NFBAS(IFRP,0) : number of AO-basis functions
!
!     modified for MCSCF coefficients, occupation numbers and mj-values 
!     by S. Knecht - April 2010
!
!     Written by T.Saue Sep 25 2000
!
!***********************************************************************
      USE MEMORY_ALLOCATOR
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "priunit.h"
#include "mxcent.h"
!
#include "dgroup.h"
#include "dcbbas.h"
#include "nuclei.h"
      DIMENSION EIG(*),CMO(*),IBEIG(*),NSTR(2,0:2),KVEC(*),NVECS(2), &
                CSEL(*),ESEL(*),IBSEL(*),IOFF(2)
!     Local variables
      LOGICAL LBIT,FNDLAB
      CHARACTER TEXT*74,CMOFIL*6,FMT*6,MXFORM*6
      DIMENSION IDIM(3,2),JDIM(2,0:3)
      REAL*8 CMCSCF_MAGIC
!
      INTEGER, ALLOCATABLE :: KBAS(:)
!
!     flag for MCSCF coefficients on file CMOFIL
      CMCSCF_MAGIC = 137.0D0
!
!     Check if the file is in formatted or unformatted form
!     =====================================================
!
!      KEYFORM = 2
!      OPEN(IUNIT,FILE=CMOFIL,STATUS='OLD',FORM='FORMATTED',ACCESS='SEQUENTIAL',IOSTAT=IER)
!      READ(IUNIT,'(A74)',IOSTAT=IER) TEXT
!      READ(IUNIT,*,IOSTAT=IER) IDUM

!      IF(IER /= 0)THEN
!        CLOSE(IUNIT,STATUS='KEEP')
        OPEN(IUNIT,FILE=CMOFIL,STATUS='OLD',FORM='UNFORMATTED',ACCESS='SEQUENTIAL',IOSTAT=IER)
        IF(IER /= 0) CALL QUIT('SELOWN_mod: ERROR while opening file')
        KEYFORM = 1
!      ENDIF
!
      REWIND IUNIT
!
!     Read title line
!     ===============
!
      IF(KEYFORM == 1)THEN
        IF(FNDLAB('INFO    ',IUNIT)) THEN      
          READ (IUNIT,END=10,ERR=20) TEXT,NSYM,NZBUF,((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
        ELSE
           CALL QUIT('SELOWN: Error reading info')
        ENDIF
!        READ(IUNIT,END=10,ERR=20) TEXT,NSYM,((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
      ELSEIF(KEYFORM == 2)THEN
        READ(IUNIT,'(A74)',IOSTAT=IER) TEXT
        IF(IER /= 0) CALL QUIT('SELOWN_mod: ERROR while reading TEXT')
        READ(IUNIT,*,IOSTAT=IER) NSYM,((IDIM(I,J),I = 1,3),J=1,NSYM)
        IF(IER /= 0) CALL QUIT('SELOWN_mod: ERROR while reading NSYM')
        READ(IUNIT,*,IOSTAT=IER) TOTERG
        IF(IER /= 0) CALL QUIT('SELOWN_mod: ERROR while reading TOTERG')
      ENDIF
!
!     set local IOPT_TMP
      IOPT_TMP = IOPT
!
!     found MCSCF coefficients on file CMOFIL?
      IF( TOTERG .eq. CMCSCF_MAGIC ) IOPT_TMP = 14
!
!     File info
!     =========
!
      IF(LBIT(IOPT_TMP,1)) THEN
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/A,A6,3X,A,'//FMT//')')           &
         'SELOWN: Coefficients read from file ',CMOFIL, &
         '- Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
        WRITE(LUPRI, '(2A, I8)') '- number of symmetry independent', &
         ' nuclei in this fragment : ', NNUC      
        WRITE(LUPRI,'(A,2I8)') &
         '- Positrons : ',(IDIM(1,I),I=1,NFSYM)
        WRITE(LUPRI,'(A,2I8)') &
         '- Electrons : ',(IDIM(2,I),I=1,NFSYM)
        WRITE(LUPRI,'(A,2I8)') &
         '- SO-basis  : ',(IDIM(3,I),I=1,NFSYM)
      ENDIF
!
!     Calculate dimensions
!     ====================
!
      NCDIM = 0
      NEDIM = 0
      DO IFRP = 1,NFSYM
        JDIM(IFRP,1) = IDIM(1,IFRP)
        JDIM(IFRP,2) = IDIM(2,IFRP)
        JDIM(IFRP,0) = IDIM(1,IFRP) + IDIM(2,IFRP)
        JDIM(IFRP,3) = IDIM(3,IFRP)                 
        NEDIM = NEDIM + (IDIM(1,IFRP)+IDIM(2,IFRP))
        NCDIM = NCDIM + (IDIM(1,IFRP)+IDIM(2,IFRP))*IDIM(3,IFRP)
      ENDDO
      NCDIM = NCDIM*NZ
!
!     Read coefficients
!     =================
!    
      IF(LBIT(IOPT_TMP,2)) THEN
        CALL REACMO_coefficients(IUNIT,CMO,NCDIM,NFSYM,NZ,    &
     &           JDIM(1,0),JDIM(1,3),JDIM(1,1),JDIM(1,2),IDIM)         
      ENDIF
!      IF(LBIT(IOPT_TMP,2)) THEN
!        CALL READT_MOD(IUNIT,NCDIM,CMO,KEYFORM)
!      ELSE
!        IF(KEYFORM == 1) READ(IUNIT)
!        IF(KEYFORM == 2) READ(IUNIT,*)
!      ENDIF
!
!     Read eigenvalues
!     ================
!
      IF(LBIT(IOPT_TMP,3)) THEN
        CALL REACMO_eigenvalues(IUNIT,EIG,NEDIM,NFSYM,NZ, &
     &           JDIM(1,0),JDIM(1,1),JDIM(1,2),IDIM)         
      ENDIF
!      IF(LBIT(IOPT_TMP,3)) THEN
!        CALL READT_MOD(IUNIT,NEDIM,EIG,KEYFORM)
!      ELSE
!        IF(KEYFORM == 1) READ(IUNIT)
!        IF(KEYFORM == 2) READ(IUNIT,*)
!      ENDIF
!
!     Read boson irrep info
!     =====================
!
      IF(LBIT(IOPT_TMP,4)) THEN
        CALL REACMO_supersymmetry(IUNIT,'SUPERSYM',IBEIG,NEDIM,NFSYM,NZ, &
     &           JDIM(1,0),JDIM(1,1),JDIM(1,2),IDIM)         
      ENDIF
!      IF(LBIT(IOPT_TMP,4)) THEN
!        CALL READI_MOD(IUNIT,NEDIM,IBEIG,KEYFORM)
!      ELSE
!        IF(KEYFORM == 1) READ(IUNIT)
!        IF(KEYFORM == 2) READ(IUNIT,*)
!      ENDIF
!
!     Select and adjust format to molecular basis
!     ===========================================
!
      ICMO  = 1
      IEIG  = 1
      ICOFF = 1
      IEOFF = 1
      DO IFRP = 1,NFSYM
      IF(NSTR(IFRP,0).GT.0) THEN
        NCMO = IDIM(1,IFRP)+IDIM(2,IFRP)
!       Generate pointer array for bases
        CALL ALLOC(KBAS,IDIM(3,IFRP),ID='KBAS in SELOWN_mod')
        CALL SELOWI_mod(INUC,NNUC,IFRP,KBAS,NSBAS)
        IF(NSBAS.NE.IDIM(3,IFRP)) THEN
          WRITE(LUPRI,'(A,I3/A,A4,A,I3)')                           &
            'SELOWN: Error in selection of coefficients in ircop ', &
            IFRP,                                                   &
            'for symmetry independent center ',NAMN(INUC),' no. ',INUC
          WRITE(LUPRI,'(A,I8)')                                     &
            'Number of basis functions for this center is:',NSBAS,  &
            'Number of basis functions in coefficient file:',       &
             IDIM(3,IFRP)
           CALL QUIT('SELOWN: Error in cf.selection !')
        ENDIF
!       Select coefficients
        IF(LBIT(IOPT_TMP,2)) THEN
          ICSEL = ICOFF + NFBAS(IFRP,0)*IOFF(IFRP)
          CALL SELOWC_mod(NZ,CMO(ICMO),IDIM(3,IFRP),NCMO,  &
                    IDIM(1,IFRP),IDIM(2,IFRP),             &
                    CSEL(ICSEL),NFBAS(IFRP,0),NVECS(IFRP), &
                    KVEC,KBAS,                             &
                    NSTR(IFRP,2),NSTR(IFRP,1))
          ICMO = ICMO + NCMO*IDIM(3,IFRP)*NZ
        ENDIF
        IESEL = IEOFF + IOFF(IFRP)
!       Select eigenvalues
        IF(LBIT(IOPT_TMP,3)) THEN
          CALL SELOWE_mod(EIG(IEIG),IDIM(1,IFRP),IDIM(2,IFRP), &
                      ESEL(IESEL),KVEC,                        &
                      NSTR(IFRP,2),NSTR(IFRP,1))
        ENDIF
!       Select boson irrep information
        IF(LBIT(IOPT_TMP,4)) THEN
          CALL SELOWB_mod(IBEIG(IEIG),IDIM(1,IFRP),IDIM(2,IFRP), &
                      IBSEL(IESEL),KVEC,                         &
                      NSTR(IFRP,2),NSTR(IFRP,1))
        ENDIF
        IEIG = IEIG + NCMO
        CALL DEALLOC(KBAS,ID='KBAS in SELOWN_mod')
        ICOFF = ICOFF + NFBAS(IFRP,0)*NVECS(IFRP)*NZ
        IEOFF = IEOFF + NVECS(IFRP)
      ENDIF
      ENDDO

!     Close file
      CLOSE(IUNIT,STATUS='KEEP')

      RETURN
 10   CONTINUE
      CALL QUIT('SELOWN: END reading TEXT')
 20   CONTINUE
      CALL QUIT('SELOWN: ERROR reading TEXT')

  END SUBROUTINE SELOWN_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE READT_MOD(IT,N,X,KEYFORM)
  IMPLICIT NONE
  INTEGER :: IER, IT, N, KEYFORM
  REAL*8  :: X(N)

  IF(KEYFORM == 1)THEN
    READ (IT,IOSTAT=IER) X
    IF(IER /= 0) CALL QUIT('READT_MOD: ERROR while reading real array in unformatted file')
  ELSEIF(KEYFORM == 2)THEN
    READ (IT,*,IOSTAT=IER) X
    IF(IER /= 0) CALL QUIT('READT_MOD: ERROR while reading real array in formatted file')
  ENDIF

  END SUBROUTINE READT_MOD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE READI_MOD(IT,N,INTX,KEYFORM)
  IMPLICIT NONE
  INTEGER :: IER, IT, N, KEYFORM
  INTEGER :: INTX(N)

  IF(KEYFORM == 1)THEN
    READ (IT,IOSTAT=IER) INTX
    IF(IER /= 0) CALL QUIT('READT_MOD: ERROR while reading integer array in unformatted file')
  ELSEIF(KEYFORM == 2)THEN
    READ (IT,*,IOSTAT=IER) INTX
    IF(IER /= 0) CALL QUIT('READT_MOD: ERROR while reading integer array in formatted file')
  ENDIF

  END SUBROUTINE READI_MOD

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE SELOWC_mod(NZ,CMO,KBAS,KORB,KPSH,KESH, &
                           CBUF,MBAS,MORB,      &
                           JVEC,JBAS,NPVEC,NEVEC)
!***********************************************************************
!     Pick out a set of vectors from CMO according to array JVEC
!
!     Written by T.Saue 1997
!     Last revision Jan 8 1997
!
!*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "priunit.h"
!
      DIMENSION CMO(KBAS,KORB,NZ),CBUF(MBAS,MORB,NZ),JVEC(*), &
                JBAS(*)
!
      NVEC = NPVEC + NEVEC
      NTOT = NVEC*MBAS
      DO IZ = 1,NZ
        CALL DZERO(CBUF(1,1,IZ),NTOT)
!
!       Positronic vectors
!
        DO I = 1,NPVEC
          II = KPSH+1+JVEC(I)
          DO J = 1,KBAS
            CBUF(JBAS(J),I,IZ) = CMO(J,II,IZ)
          ENDDO
        ENDDO
!
!       Electronic vectors
!
        DO I = NPVEC+1,NVEC
          II = KPSH+JVEC(I)
          DO J = 1,KBAS
            CBUF(JBAS(J),I,IZ) = CMO(J,II,IZ)
          ENDDO
        ENDDO
      ENDDO
!
      RETURN
  END SUBROUTINE SELOWC_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE SELOWE_mod(EIG,KPSH,KESH,EBUF,JVEC,NPVEC,NEVEC)
!***********************************************************************
!     Pick out a set of eigenvalues from EIG according to array JVEC
!
!     Written by L. Visscher 1997
!
!*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "priunit.h"
!
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION EIG(*),EBUF(NPVEC+NEVEC),JVEC(*)
!
      NVEC = NPVEC + NEVEC
!
!
!     Positronic eigenvalues
!
      DO I = 1,NPVEC
        II = KPSH+1+JVEC(I)
        EBUF(I) = EIG(II)
      ENDDO
!
!     Electronic eigenvalues
!
      DO I = NPVEC+1,NVEC
        II = KPSH+JVEC(I)
        EBUF(I) = EIG(II)
      ENDDO
!
      RETURN
  END SUBROUTINE SELOWE_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE SELOWB_mod(IBEIG,KPSH,KESH,IBBUF,JVEC,NPVEC,NEVEC)
!***********************************************************************
!     Pick out a set of from boson irrep info array IBEIG 
!     according to array JVEC
!
!     Written by L. Visscher 1997
!
!*****************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "priunit.h"
!
#include "dgroup.h"
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION IBEIG(*),IBBUF(NPVEC+NEVEC),JVEC(*)
!
      NVEC = NPVEC + NEVEC
!
!
!     Positronic eigenvalues
!
      DO I = 1,NPVEC
        II = KPSH+1+JVEC(I)
        IBBUF(I) = IBEIG(II)
      ENDDO
!
!     Electronic eigenvalues
!
      DO I = NPVEC+1,NVEC
        II = KPSH+JVEC(I)
        IBBUF(I) = IBEIG(II)
      ENDDO
!
      RETURN
  END SUBROUTINE SELOWB_mod

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&!

  SUBROUTINE SELOWI_mod(INUC,NNUC,IFRP,JBAS,NSBAS)
!***********************************************************************
!
!     Make index array pointing from basis functions of
!     symmetry independent centers INUC..INUC+NNUC-1 to the full
!     molecular basis
!
!     Written by T.Saue Sep 25 2000
!
!***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER (I-N)
#include "priunit.h"
#include "maxorb.h"
!
#include "dcbbas.h"
#include "dcblab.h"
      DIMENSION JBAS(*)
!
      NSBAS = 0
      DO I = 1,NFBAS(IFRP,0)
        II = IBAS(IFRP) + I
        ILAB = IPLAB(II,2)
        ICENT = JGET(IATTR(ILAB,2))
        IF ((ICENT.GE.INUC).AND.(ICENT.LT.(INUC+NNUC))) THEN
          NSBAS = NSBAS + 1
          JBAS(NSBAS) = I
        ENDIF
      ENDDO
!
      RETURN
  END SUBROUTINE SELOWI_mod

  FUNCTION JGET(I)
  INTEGER I, JGET
!
#if defined (SYS_CRAY) || defined (INT_STAR8)
      JGET = IAND(ISHFT(I,-32),65535)
#else
      JGET = IAND(ISHFT(I,-16),255)
#endif
!
  END FUNCTION JGET

!================================================================================!
!                             END OF THE MODULE                                  !
!================================================================================!
end module select_coeff_own_bas
