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

C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmoin */
      SUBROUTINE PCMOIN()
C*****************************************************************************
C
C     Read coefficients from formatted DFPCMO and
C     write to unformatted DFCOEF
C
C         IDIM(1,IFRP) - number of positronic solutions
C         IDIM(2,IFRP) - number of electronic solutions
C         IDIM(3,IFRP) - number of AO-basis functions
C
C     Written by T.Saue Feb 1997
C
C*****************************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbpsi.h"      
#include "dcbbas.h"
#include "dcborb.h"
#include "dgroup.h"
C
      real(8), allocatable :: CMO(:),EIG(:)
      integer, allocatable :: IBEIG(:,:)
      integer              :: IDIM(3,2)
      logical              :: TOBE,SSYM
      character            :: LABEL*8,TEXT*74

      CALL QENTER('PCMOIN')
C
C     Get coefficients from DFPCMO;
C     If file is not there, return quietly      
C
      INQUIRE(FILE='DFPCMO',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
        CALL QEXIT('PCMOIN')
        RETURN
      ENDIF
      OPEN(LUCOEF,FILE='DFPCMO',STATUS='OLD',FORM='FORMATTED',
     +     ACCESS='SEQUENTIAL')
C
      allocate(CMO(NCMOTQ*NKRBLK))
      NEIG=NORBT*NKRBLK
      allocate(EIG(NEIG))
      allocate(IBEIG(NEIG,2))
      SSYM=.FALSE.
      REWIND LUCOEF
 10   CONTINUE
      READ(LUCOEF,'(A8)',END=20) LABEL
      SELECT CASE (LABEL)
        CASE('INFO    ')
          READ(LUCOEF,'(A74)',iostat=ierr) TEXT
          IF (IERR /= 0) CALL QUIT('error reading DFPCMO: title line')
          READ(LUCOEF,*,iostat=ierr)
     &         NFSYM2,NDIM,((IDIM(I,J),I = 1,3),J=1,NFSYM2)
          IF (IERR /= 0) CALL 
     &         QUIT('error reading DFPCMO: symmetry information')
C         Check dimensions
          IF(NFSYM2.NE.NFSYM) THEN
             WRITE(LUPRI,'()')
     &         '* Expected number of fermion ircops: ',NFSYM,
     &         '   - from DFPCMO                   : ',NFSYM2
             CALL QUIT('PCMOIN: Error in symmetry information')
          ENDIF
          NCMOTQ2=0
          DO IFRP = 1,NFSYM
            NCMOTQ2 = NCMOTQ2 + (IDIM(1,IFRP)+IDIM(2,IFRP))*IDIM(3,IFRP)
          ENDDO
          NCMOTQ2 = NCMOTQ2*NDIM
          IF(NCMOTQ2.NE.NCMOTQ*NKRBLK) THEN
             WRITE(LUPRI,'(25X,A20)')       ' Expected /From file'       
             DO IFRP = 1,NFSYM
                WRITE(6,'(3X,A22,2I10)') 'Positronic solutions: ',
     &                     IDIM(1,IFRP),NPSH(IFRP)
                WRITE(6,'(3X,A22,2I10)') 'Electronic solutions: ',
     &                     IDIM(2,IFRP),NESH(IFRP)
                WRITE(6,'(3X,A22,2I10)') 'AO basis functions  : ',
     &                     IDIM(3,IFRP),NFBAS(IFRP,0)
             ENDDO   
             CALL QUIT('PCMOIN: Error in dimension')             
          ENDIF
          READ(LUCOEF,'(E24.16)',iostat=ierr) TOTERG
          IF (IERR /= 0) CALL QUIT('error reading DFPCMO: energy')
          GOTO 10
        CASE('COEFS   ')
          READ(LUCOEF,'(6F22.16)',iostat=ierr) CMO
          IF (IERR /= 0) CALL QUIT('error reading DFPCMO: coefficients')
          GOTO 10
        CASE('EVALS   ')
          READ(LUCOEF,'(6E22.12)',iostat=ierr) EIG
          IF (IERR /= 0) CALL QUIT('error reading DFPCMO: eigenvalues')
          GOTO 10
        CASE('SUPERSYM')
          READ(LUCOEF,*,iostat=ierr) (IBEIG(I,1),I=1,NEIG)
          IF (IERR /= 0) CALL QUIT('error reading DFPCMO: ibeig')
          GOTO 10
        CASE('KAPPA   ')
          SSYM=.TRUE.
          READ(LUCOEF,*,iostat=ierr) (IBEIG(I,2),I=1,NEIG)
          IF (IERR /= 0) CALL QUIT('error reading DFPCMO: ibeig2')
          GOTO 10
      END SELECT
 20   CONTINUE
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     Write coefficients to either KRMCOLD or DFCOEF
C
      IF(DOKRMC.or.DOKRCI.or.DOKRCC)THEN
        CALL OPNFIL(LUKRMC,'KRMCOLD','UNKNOWN','PCMOIN')
        CALL WRTKRMC(LUKRMC,'NEWORB  ',CMO,NCMOTQ*NKRBLK)
        CLOSE (LUKRMC,STATUS='KEEP')
        IWRTMO = 1
      ELSE
        CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','PCMOIN')        
        CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG,SSYM)
        CLOSE(LUCOEF,STATUS='KEEP')
        IWRTMO = 2
      END IF
      deallocate(CMO)
      deallocate(EIG)
      deallocate(IBEIG)
C
      CALL PRSYMB(LUPRI,'=',75,0)
      IF( IWRTMO .eq. 1 )THEN
        WRITE(LUPRI,'(A)')
     &      '* PCMOIN: Coefficients read from formatted DFPCMO ',
     &      '          and written to unformatted KRMCOLD'
      ELSE IF( IWRTMO .eq. 2 )THEN
        WRITE(LUPRI,'(A)')
     &      '* PCMOIN: Coefficients read from formatted DFPCMO ',
     &      '          and written to unformatted DFCOEF'
      END IF
      CALL PRSYMB(LUPRI,'=',75,0)

      CALL QEXIT('PCMOIN')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmout */
      SUBROUTINE PCMOUT()
C*****************************************************************************
C
C     Read coefficients from unformatted DFCOEF and
C     write to formatted DFPCMO
C
C     Written by T.Saue Feb 1997
C
C*****************************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbbas.h"
#include "dcborb.h"`
#include "dgroup.h"`      
C
      real(8), allocatable :: CMO(:),EIG(:)
      integer, allocatable :: IBEIG(:,:)
      logical              :: TOBE,TOBEK

      CALL QENTER('PCMOUT')

      NEIGB=NORBT*NKRBLK
      NCMOB=NCMOTQ*NKRBLK
      NDIMB=NZ*NKRBLK
      allocate(CMO(NCMOB))
      allocate(EIG(NEIGB))
      allocate(IBEIG(NEIGB,2))
C
C     Get coefficients from either KRMCSCF or DFCOEF 
C
      INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
      INQUIRE(FILE='DFCOEF' ,EXIST=TOBE)
      IF((.NOT.TOBE).and.(.NOT.TOBEK)) THEN
        WRITE(LUPRI,'(A)') 
     &     'PCMOUT: Unformatted file DFCOEF/KRMCSCF not found !'
        CALL QUIT('PCMOUT: Unformatted file DFCOEF/KRMCSCF not found !')
      ENDIF
      IF(TOBEK.and.(DOKRMC.or.DOKRCI.or.DOKRCC))THEN
        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','PCMOUT')
        JRDMO = -1
        CALL RREADMO(CMO,JRDMO,1,LUKRMC)
        CLOSE (LUKRMC,STATUS='KEEP')
C       ... set to zero
        TOTERG = 0.0D0
        CALL DZERO(EIG,NORBT)
        CALL IZERO(IBEIG,NORBT)
        IREADMO = 1
      ELSE IF(TOBEK.and.(.not.DOKRMC.or..not.DOKRCI.or..not.DOKRCC))THEN
C
C       read CI/MCSCF natural orbitals and occupancies
        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','ACMOU1')
        JRDMO = -1
        CALL RREADMO(CMO,JRDMO,4,LUKRMC)
        CALL REAKRMC(LUKRMC,'MCNATOCC',EIG,NORBT)
        CLOSE (LUKRMC,STATUS='KEEP')
C
C       temporary fix for use of occupancies in Mulliken population
C       analysis - set energy to 137.0D0 - SK Feb 2009
        TOTERG = 137.0D0
C       ... set to zero
        CALL IZERO(IBEIG,NORBT)
        IREADMO = 3
      ELSE
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PCMOUT')
        IOPT=14
        IF(ATOMIC) IOPT=IOPT+16
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)
        CLOSE(LUCOEF,STATUS='KEEP')
        IREADMO = 2
      END IF
C
C     Write coefficients to DFPCMO
C
      OPEN(LUCOEF,FILE='DFPCMO',STATUS='UNKNOWN',
     &            FORM='FORMATTED',ACCESS='SEQUENTIAL')
      CALL WRIPCMO(LUCOEF,CMO,EIG,IBEIG,TOTERG,ATOMIC)
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL PRSYMB(LUPRI,'=',75,0)
      IF( IREADMO .eq. 1 .or. IREADMO .eq. 3 )THEN
        WRITE(LUPRI,'(A)')
     &      '* PCMOUT: Coefficients read from unformatted KRMCSCF',
     &      '          and written to formatted DFPCMO'
      ELSE IF( IREADMO .eq. 2)THEN
        WRITE(LUPRI,'(A)')
     &      '* PCMOUT: Coefficients read from unformatted DFCOEF ',
     &      '          and written to formatted DFPCMO'
      END IF
      CALL PRSYMB(LUPRI,'=',75,0)

      deallocate(CMO)
      deallocate(EIG)
      deallocate(IBEIG)

      CALL QEXIT('PCMOUT')
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wricmo */
      SUBROUTINE WRICMO(IUNIT,CMO,EIG,IBEIG,TOTERG,SSYM)
C*****************************************************************************
C
C     Write to file:
C       CMO    - coefficients
C       EIG    - eigenvalues    
C       IBEIG  - supersymmetry information 
C       TOTERG - total energy
C
C     SSYM = T activates writing of second array of supersymmetry information
C
C     IDIM(1,J) - number of positronic solutions
C     IDIM(2,J) - number of electronic solutions
C     IDIM(3,J) - number of AO-basis functions
C
C     Written by T.Saue Sept 1 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ*NKRBLK),EIG(NORBT*NKRBLK),
     &          IBEIG(NORBT*NKRBLK,2)
      DIMENSION IDIM(3,2)
      LOGICAL   SSYM
C
      DO I = 1,NFSYM
        IDIM(1,I) = NPSH(I)
        IDIM(2,I) = NESH(I)
        IDIM(3,I) = NFBAS(I,0)
      ENDDO
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT     
      CALL NEWLAB('INFO    ',IUNIT,LUPRI)
      WRITE(IUNIT) TEXT,NFSYM,NZ*NKRBLK,((IDIM(I,J),I = 1,3),J=1,NFSYM),
     &     TOTERG
      CALL WRICMO_coefficients(IUNIT,CMO,NCMOTQ*NKRBLK)
      CALL WRICMO_eigenvalues(IUNIT,EIG,NORBT*NKRBLK)
      CALL WRICMO_supersymmetry(IUNIT,'SUPERSYM',IBEIG,NORBT*NKRBLK)
      IF(SSYM) CALL WRICMO_supersymmetry(IUNIT,'KAPPA   ',
     &  IBEIG(1,2),NORBT*NKRBLK)
!     Add additional information about basis set to create a simple interface file
!     for correlated calculations (like exacorr)
      CALL WRICMO_basisinfo(IUNIT)
      CALL NEWLAB('EOFLABEL',IUNIT,LUPRI)                              
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck WRICMO_coeffcients */
      SUBROUTINE WRICMO_coefficients(IUNIT,CMO,NCMO)
C***********************************************************************
C     Write SCF-coefficients to unformatted file
C     Written by Trond Saue Jan 2021      
C***********************************************************************
#include "priunit.h"
      INTEGER IUNIT,NCMO
      REAL*8 CMO(NCMO)
      CALL NEWLAB('COEFS   ',IUNIT,LUPRI)      
      WRITE(IUNIT) CMO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck WRICMO_eigenvalues */
      SUBROUTINE WRICMO_eigenvalues(IUNIT,EIG,NEIG)
C***********************************************************************
C     Write SCF-orbital energies to unformatted file
C     Written by Trond Saue Jan 2021      
C***********************************************************************
#include "priunit.h"
      INTEGER IUNIT,NEIG
      REAL*8 EIG(NEIG)
      CALL NEWLAB('EVALS   ',IUNIT,LUPRI)            
      WRITE(IUNIT) EIG
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck WRICMO_supersymmetry */
      SUBROUTINE WRICMO_supersymmetry(IUNIT,LABEL,IBEIG,NEIG)
C***********************************************************************
C     Write SCF-supersymmetry string to unformatted file
C     Written by Trond Saue Jan 2021      
C***********************************************************************
#include "priunit.h"
      INTEGER IUNIT,NEIG,IBEIG(NEIG)
      CHARACTER*8 LABEL
      CALL NEWLAB(LABEL,IUNIT,LUPRI)                  
      WRITE(IUNIT) IBEIG
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck WRICMO_basisinfo */
      SUBROUTINE WRICMO_basisinfo(IUNIT)
C***********************************************************************
C     Write SCF-coefficients to unformatted file
C     Written by Trond Saue Jan 2021      
C***********************************************************************
#include "priunit.h"
      INTEGER IUNIT
      CALL NEWLAB('BASINFO ',IUNIT,LUPRI)                        
      call append_basisinfo_to_dfcoef (IUNIT)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine append_basisinfo_to_dfcoef (IUNIT)

!      L. Visscher, October 2018
!      add basis set info to the end of dfcoef to not infere with normal usage
!      a better solution is to use a keyed-file system, but this is quite some
!      work to implement everywhere in DIRAC

!      IUNIT should be open and positioned at the right place (i.e. after IBEIG)

#include "implicit.h"
#include "maxorb.h"
#include "shells.h"
#include "aovec.h"
#include "primit.h"
#include "dcbbas.h"

       integer npriexp

!      First calculate and write dimensions so that we can allocate when reading this back
!      number of aos, number of shells (LC only)
       npriexp = jstrt(nlrgsh)+nuco(nlrgsh)
       write (IUNIT) ntbas(1),nlrgsh,npriexp,MXCONT

!      Now write arrays with basis set information
!      orb_momentum, atom_number, n_primives for all shells
       write (IUNIT) nhkt(1:nlrgsh),ncent(1:nlrgsh),nuco(1:nlrgsh)
!      start of exponent values, number of contracted for all shells
       write (IUNIT) jstrt(1:nlrgsh),numcf(1:nlrgsh)
!      coordinates (x, y, z) of center for all shells
       write (IUNIT) cent(1:nlrgsh,1:3,1)
!      values of the exponents for all shells (packed into a 1-d array)
       write (IUNIT) priexp(1:npriexp)
!      values of the contraction coefficients for all shells (packed into a 1-d array)
       write (IUNIT) priccf(1:npriexp,1:MXCONT)

       end subroutine append_basisinfo_to_dfcoef
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wripcmo */
      SUBROUTINE WRIPCMO(IUNIT,CMO,EIG,IBEIG,TOTERG,SSYM)
C*****************************************************************************
C
C     Write SCF-coefficients and eigenvalues to unformatted file
C
C     IDIM(1,J) - number of positronic solutions
C     IDIM(2,J) - number of electronic solutions
C     IDIM(3,J) - number of AO-basis functions

C     Written by BHP (adapted from TS WRICMO)
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "dgroup.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcbgen.h"
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ*NKRBLK),EIG(NORBT*NKRBLK),
     &          IBEIG(NORBT*NKRBLK,2)
      DIMENSION IDIM(3,2)
      LOGICAL   SSYM
C
      DO I = 1,NFSYM
        IDIM(1,I) = NPSH(I)
        IDIM(2,I) = NESH(I)
        IDIM(3,I) = NFBAS(I,0)
      ENDDO
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
C
      REWIND IUNIT
      WRITE(IUNIT,'(A8)') 'INFO    '
      WRITE(IUNIT,'(A74)') TEXT
      WRITE(IUNIT,'(8(X,I0))')
     &      NFSYM,NZ*NKRBLK,((IDIM(I,J),I = 1,3),J=1,NFSYM)
      WRITE(IUNIT,'(E24.16)') TOTERG
      WRITE(IUNIT,'(A8)') 'COEFS   '      
      WRITE(IUNIT,'(6F22.16)') CMO
      WRITE(IUNIT,'(A8)') 'EVALS   '            
      WRITE(IUNIT,'(6E22.12)') EIG
      WRITE(IUNIT,'(A8)') 'SUPERSYM'
      WRITE(IUNIT,'(66(X,I0))') (IBEIG(I,1),I=1,NORBT*NKRBLK)
      IF(SSYM) THEN
        WRITE(IUNIT,'(A8)') 'KAPPA   '
        WRITE(IUNIT,'(66(X,I0))') (IBEIG(I,2),I=1,NORBT*NKRBLK)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reacmo */
      SUBROUTINE REACMO(IUNIT,CMOFIL,CMO,EIG,IBEIG,TOTERG,IOPT)
C***********************************************************************
C
C     Read SCF-coefficients and eigenvalues from unformatted file
C     Read option is provided from bit-packed IOPT:
C      00001 - give restart info
C      00010 - read coefficients
C      00100 - read eigenvalues
C      01000 - read boson irrep identification
C      10000 - read second supersymmetry id (kappa atomic symmetry)
C     
C     IDIM(1,IFRP) - number of positronic solutions
C     IDIM(2,IFRP) - number of electronic solutions
C     IDIM(3,IFRP) - number of AO-basis functions
C
C     Written by T.Saue Sept 1 1995
C
C***********************************************************************
         use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL INFO,LBIT,TEST,FNDLAB
      CHARACTER TEXT*74,CMOFIL*6,FMT*6,MXFORM*6
      DIMENSION CMO(NCMOTQ*NKRBLK),EIG(NORBT*NKRBLK),
     &          IBEIG(NORBT*NKRBLK,2)
      DIMENSION IDIM(3,2)
C
      IF(IOPT.EQ.0) RETURN
C
C     Read title line - always done
C     =============================
C
      REWIND IUNIT
      CALL REACMO_info(IUNIT,TEXT,NSYM,NZBUF,IDIM,TOTERG)
      IF(IOPT.EQ.-1) RETURN
C
C     Compatibility check for NFSYM and NFBAS
C     =======================================
C
      IERR = 0
      IF(NSYM .NE. NFSYM) THEN
         IERR = IERR + 1
         INFO = .TRUE.
         WRITE(LUPRI,'(/A/A,2I5)')
     &      'REACMO FATAL ERROR: Incompatible number of fermion ircops',
     &      ' from file and in this run:',NSYM,NFSYM
      ELSEIF(NZBUF .NE. NZ*NKRBLK) THEN
         IERR = IERR + 1
         INFO = .TRUE.
         WRITE(LUPRI,'(/A/A,2I5)')
     &      'REACMO FATAL ERROR: Incompatible NZ parameter',
     &      ' from file and in this run:',NZBUF/NKRBLK,NZ
      ELSE
         DO IFRP = 1,NFSYM
           IF(IDIM(3,IFRP).NE.NFBAS(IFRP,0)) IERR = IERR + 1
         END DO
         IF (IERR.GT.0) THEN
           INFO = .TRUE.
           WRITE(LUPRI,'(/A/A,2(/A,I2,2I5))')
     &     'REACMO FATAL ERROR. Incompatible number of basis functions',
     &     ' from file and in this run:', (' Fermion ircop no.',IFRP,
     &     IDIM(3,IFRP),NFBAS(IFRP,0),IFRP=1,NFSYM)
         END IF
      END IF
C
C     File info
C     =========
C
      IF(LBIT(IOPT,1)) THEN
        FMT = MXFORM(TOTERG,20)
        WRITE(LUPRI,'(/3A,'//FMT//')')
     &   '* REACMO: Coefficients read from file ',CMOFIL,
     &   ' - Total energy: ',TOTERG
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
C
      IF (IERR.GT.0) CALL QUIT(
     &'REACMO: old DFCOEF file is not consistent with this calculation')
C
      IF(IOPT.LT.2) RETURN
C
C     Read coefficients
C     =================
C
      IF(LBIT(IOPT,2)) THEN
        NCMO=NCMOTQ*NKRBLK
        CALL REACMO_coefficients(IUNIT,CMO,NCMO,NFSYM,NZ,
     &           NORB,NFBAS(1,0),NPSH,NESH,IDIM)         
      ENDIF
C
C     Read eigenvalues
C     ================
C
      IF(LBIT(IOPT,3)) THEN
        NEIG=NORBT*NKRBLK
        CALL REACMO_eigenvalues(IUNIT,EIG,NEIG,NFSYM,NZ,
     &           NORB,NPSH,NESH,IDIM)         
      ENDIF
C
C     Read supersymmetry info
C     =======================
C
      IF(LBIT(IOPT,4)) THEN
        NEIG=NORBT*NKRBLK
        CALL REACMO_supersymmetry(IUNIT,'SUPERSYM',IBEIG,NEIG,NFSYM,NZ,
     &       NORB,NPSH,NESH,IDIM)
      ENDIF
C
C     Read second supersymmetry info (kappa atomic)
C     =============================================
C
      IF(LBIT(IOPT,5)) THEN
        NEIG=NORBT*NKRBLK
        CALL REACMO_supersymmetry(IUNIT,'KAPPA   ',IBEIG(1,2),
     &       NEIG,NFSYM,NZ,NORB,NPSH,NESH,IDIM)
      ENDIF
C
      RETURN
C
C     Error messages
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck reacmo_info */
      SUBROUTINE REACMO_info(IUNIT,TEXT,NFSYM,NZ,IDIM,TOTERG)
      implicit none
      integer, intent(in)    :: IUNIT
      integer, intent(out)   :: NFSYM,NZ,IDIM(3,2)
      character, intent(out) :: TEXT*74
      real*8, intent(out)    :: TOTERG
      integer                :: I,J
      LOGICAL FNDLAB
      IF(FNDLAB('INFO    ',IUNIT)) THEN
        READ (IUNIT,END=10,ERR=20) TEXT,NFSYM,NZ,
     &       ((IDIM(I,J),I = 1,3),J=1,NFSYM),TOTERG
      ELSE
        GOTO 20        
      ENDIF
      RETURN
 10   CONTINUE
      write(6,'(//A/A/A)')
     &   'REACMO: END OF FILE reading first record on DFCOEF.',
     &   '        Possible reason: DFCOEF was generated with'//
     &      ' 32-bit integer Dirac version.',
     &   '        This Dirac is compiled with 64 bit integers.'
      CALL QUIT('REACMO: END OF FILE reading first record')
 20   CONTINUE
      write(6,'(//A/A/A)')
     &   'REACMO: ERROR reading first record on DFCOEF.',
     &   '        Possible reason: DFCOEF was generated with'//
     &      ' 32-bit integer Dirac version.',
     &   '        This Dirac is compiled with 64 bit integers.'
      CALL QUIT('REACMO: ERROR reading first record')
      END
      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck reacmo_coefficients */
      SUBROUTINE REACMO_coefficients(IUNIT,CMO,NCMO,NFSYM,NZ,
     &           NORB,NFBAS,NPSH,NESH,IDIM)
C***********************************************************************
C
C     Read coefficients from unformatted file
C
C     IDIM(1,IFRP) - number of positronic solutions
C     IDIM(2,IFRP) - number of electronic solutions
C     IDIM(3,IFRP) - number of AO-basis functions
C
C     Written by T.Saue Jan 2021
C
C***********************************************************************
         use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION CMO(NCMO)
      DIMENSION NORB(NFSYM),NFBAS(NFSYM),NPSH(NFSYM),NESH(NFSYM)
      DIMENSION IDIM(3,NFSYM),NORB2(NFSYM)
      real(8), allocatable :: buf(:)
      LOGICAL TEST,FNDLAB
      IF(FNDLAB('COEFS   ',IUNIT)) THEN
C
C     Compatibility tests for number of orbitals
C
        TEST = .FALSE.
        NC = 0
        DO IFRP = 1,NFSYM
          NORB2(IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
          IF(NORB2(IFRP).NE.NORB(IFRP)) TEST=.TRUE.
          NC = NC + NORB2(IFRP)*NFBAS(IFRP)*NZ
       ENDDO
       IF(TEST) THEN
C
C         If dimensions do not fit, then the coefficients are read
C         into a buffer and adapted for the correct dimensions    
C
          call alloc(BUF,NC)
          CALL READT(IUNIT,NC,BUF)
          IC = 1
          JC = 1
          DO IFRP = 1,NFSYM
            CALL ADACMO(NFBAS(IFRP),NZ,
     &         CMO(IC),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &         BUF(JC),NORB2(IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
            IC = IC + NORB(IFRP) *NFBAS(IFRP)*NZ
            JC = JC + NORB2(IFRP)*NFBAS(IFRP)*NZ
          ENDDO
          call dealloc(buf)
        ELSE
C
C         Default case
C           
          READ (IUNIT,END=10,ERR=20) CMO
       ENDIF
      ELSE
        GOTO 20 
      ENDIF

      RETURN
 10   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading coefficients')
 20   CONTINUE
      CALL QUIT('REACMO: ERROR reading coefficients')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck reacmo_eigenvalues */
      SUBROUTINE REACMO_eigenvalues(IUNIT,EIG,NEIG,NFSYM,NZ,
     &           NORB,NPSH,NESH,IDIM)
C***********************************************************************
C
C     Read eigenvalues from unformatted file
C
C     IDIM(1,IFRP) - number of positronic solutions
C     IDIM(2,IFRP) - number of electronic solutions
C     IDIM(3,IFRP) - number of AO-basis functions
C
C     Written by T.Saue Jan 2021
C
C***********************************************************************
         use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION EIG(NEIG)
      DIMENSION NORB(NFSYM),NPSH(NFSYM),NESH(NFSYM)
      DIMENSION IDIM(3,NFSYM),NORB2(NFSYM)
      real(8), allocatable :: buf(:)
      LOGICAL TEST,FNDLAB
      IF(FNDLAB('EVALS   ',IUNIT)) THEN         
C
C     Compatibility tests for number of orbitals
C
        TEST = .FALSE.
        NO = 0
        DO IFRP = 1,NFSYM
          NORB2(IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
          IF(NORB2(IFRP).NE.NORB(IFRP)) TEST=.TRUE.
          NO = NO + NORB2(IFRP)
        ENDDO
        IF(TEST) THEN
C
C         If dimensions do not fit, then the eigenvalues are read
C         into a buffer and adapted for the correct dimensions    
C
            call alloc(BUF,NO)
            CALL READT(IUNIT,NO,BUF)
            IO = 1
            JO = 1
            DO IFRP = 1,NFSYM
              CALL ADAEIG(EIG(IO),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &              BUF(JO),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
              IO = IO + NORB(IFRP)
              JO = JO + NORB2(IFRP)
            ENDDO
            call dealloc(BUF)
        ELSE
C
C         Default case
C           
          READ (IUNIT,END=10,ERR=20) EIG
        ENDIF
      ELSE
        GOTO 20 
      ENDIF

      RETURN
 10   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading eigenvalues')
 20   CONTINUE
      CALL QUIT('REACMO: ERROR reading eigenvalues')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck reacmo_supersymmetry */
      SUBROUTINE REACMO_supersymmetry(IUNIT,LABEL,IBEIG,NEIG,NFSYM,NZ,
     &           NORB,NPSH,NESH,IDIM)
C***********************************************************************
C
C     Read supersymmetry information from unformatted file
C
C     IDIM(1,IFRP) - number of positronic solutions
C     IDIM(2,IFRP) - number of electronic solutions
C     IDIM(3,IFRP) - number of AO-basis functions
C
C     Written by T.Saue Jan 2021
C
C***********************************************************************
         use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
      DIMENSION IBEIG(NEIG)
      CHARACTER*8 LABEL   
      DIMENSION NORB(NFSYM),NPSH(NFSYM),NESH(NFSYM)
      DIMENSION IDIM(3,NFSYM),NORB2(NFSYM)
      integer, allocatable :: buf(:)
      LOGICAL TEST,FNDLAB
      IF(FNDLAB(LABEL,IUNIT)) THEN         
C
C     Compatibility tests for number of orbitals
C
        TEST = .FALSE.
        NO = 0
        DO IFRP = 1,NFSYM
          NORB2(IFRP) = IDIM(1,IFRP)+IDIM(2,IFRP)
          IF(NORB2(IFRP).NE.NORB(IFRP)) TEST=.TRUE.
          NO = NO + NORB2(IFRP)
        ENDDO
        IF(TEST) THEN
C
C         If dimensions do not fit, then the eigenvalues are read
C         into a buffer and adapted for the correct dimensions    
C
            call alloc(BUF,NO)
            CALL READI(IUNIT,NO,BUF)
            IO = 1
            JO = 1
            DO IFRP = 1,NFSYM
              CALL ADAIBO(IBEIG(IO),NORB(IFRP),NPSH(IFRP),NESH(IFRP),
     &              BUF(JO),IDIM(3,IFRP),IDIM(1,IFRP),IDIM(2,IFRP))
              IO = IO + NORB(IFRP)
              JO = JO + NORB2(IFRP)
            ENDDO
            call dealloc(BUF)
        ELSE
C
C         Default case
C           
          READ (IUNIT,END=10,ERR=20) IBEIG
        ENDIF
      ELSE
        GOTO 20 
      ENDIF

      RETURN
 10   CONTINUE
      CALL QUIT('REACMO: END OF FILE reading supersymmetry '
     &  //LABEL)
 20   CONTINUE
      CALL QUIT('REACMO: ERROR reading irrep supersymmetry '
     &  //LABEL)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck pcmout_1 */
      SUBROUTINE PCMOUT_1(CMO,EIG,IBEIG)
C*****************************************************************************
C
C     Write SCF/KRMCSCF-coefficients and eigenvalues to FORMATTED file
C
C     Written by T.Saue Feb 1997
C     extended for KRMC modules - S. Knecht Nov 2008
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL TOBE, TOBEK
      CHARACTER TEXT*74
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT)
      DIMENSION IDIM(3,2)
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wrifck */
      SUBROUTINE WRIFCK(IUNIT,FOCK,N)
C*****************************************************************************
C
C     Write N AO - Fock matrices to unformatted file
C
C     Written by T.Saue Sept 1 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbgen.h"
#include "dcbdhf.h"
      CHARACTER TEXT*74
      DIMENSION FOCK(N2BBASXQ*N)
C
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT
      WRITE(IUNIT) TEXT
      WRITE(IUNIT) N2BBASXQ,N
      WRITE(IUNIT) FOCK
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck reafck */
      SUBROUTINE REAFCK(IUNIT,FOCK,INFO,N)
C*****************************************************************************
C
C     Read AO Fock matrix from unformatted file
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbdhf.h"
      LOGICAL INFO, INFO_loc
      CHARACTER TEXT*74, FILNAM*60
      DIMENSION FOCK(N2BBASXQ*N)
C
      info_loc = info
!     infoc_loc = .true.
      REWIND IUNIT
      READ(IUNIT,END=10,ERR=20) TEXT
      READ(IUNIT) M2BBASXQ,M
      IF (N2BBASXQ .NE. M2BBASXQ .OR. N .GT. M) INFO_loc = .TRUE.
      IF(INFO_loc) THEN
        INQUIRE(UNIT=IUNIT, NAME=FILNAM)
        WRITE(LUPRI,'(/2A)')
     &   '* REAFCK: Fock matrix read from file ',FILNAM
        WRITE(LUPRI,'(2A)') '* Heading :',TEXT
      ENDIF
      IF (N2BBASXQ .NE. M2BBASXQ .OR. N .GT. M) THEN
         WRITE(LUPRI,'(/A,2(/A,2I10))') '*** FATAL ERROR ***',
     &      'Fock matrix size in this run and size on file: ',
     &      N2BBASXQ,M2BBASXQ,
     &      'Number of Fock matrices in this run and on file:',N,M
         CALL QUIT('REAFCK: inconsistent dimensions on file'//FILNAM)
      END IF
      READ (IUNIT,END=10,ERR=30) FOCK
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
C     End of file
      CALL QUIT('REAFCK: End of file')
 20   CONTINUE
C     Error reading file
      CALL QUIT('REAFCK: ERROR reading file - text part')
 30   CONTINUE
      CALL QUIT('REAFCK: ERROR reading file - fock part')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck wridns */
      SUBROUTINE WRIDNS(IUNIT,DMAT)
C*****************************************************************************
C
C     Write SO - density matrix to unformatted file
C
C     Written by T.Saue Oct 12 1995
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbgen.h"
#include "dcbdhf.h"
      CHARACTER TEXT*74
      DIMENSION DMAT(N2BBASXQ*NFMAT)
C
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      REWIND IUNIT
      WRITE(IUNIT) TEXT
      WRITE(IUNIT) DMAT
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readns */
      SUBROUTINE READNS(IUNIT,DMAT)
C*****************************************************************************
C
C     Read SO density matrix from unformatted file
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcbdhf.h"
      CHARACTER TEXT*74
      DIMENSION DMAT(N2BBASXQ*NFMAT)
C
      REWIND IUNIT
      READ(IUNIT,END=10,ERR=20) TEXT
      READ (IUNIT,END=30,ERR=40) DMAT
C
      RETURN
C
C     Error messages
C
 10   CONTINUE
      CALL QUIT('READNS: END reading TEXT')
 20   CONTINUE
      CALL QUIT('READNS: ERROR reading TEXT')
 30   CONTINUE
      CALL QUIT('READNS: END reading record 2')
 40   CONTINUE
      CALL QUIT('READNS: ERROR reading record 2')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numlst */
      SUBROUTINE NUMLST(A,ILIST,NLIST,IRNGE1,IRNGE2,IFRP,IP)
C*****************************************************************************
C
C     Original purpose was to read the numbers in the string A which has
C     the format : '3..6, 8, -27..2, 123, 14..15'.
C     
C     !! Number 0 (Zero) NOT counted or stored !!
C
C     Input: 
C           IF(IP.NE.0) : Give list of numbers in array
C                         ILIST(NLIST).
C           ELSE : Only IP as OUTPUT
C
C     Output:
C           IP : number of elements in the list A.
C                Duplicates only counted once.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     ITMPAR   : Temporary storage
C     If the string starts with 'all', all elements in the range enters.
C
C                           Jon K. Laerdahl 30.12.96
C     Minor polish T.Saue Jan 7 1997
C
C     Additional feature (LV, 15-5-2001) : Determine list based on energy or
C     other criteria. This is more user-friendly for standard cases.
C
C     This option is toggled by starting the string by energy, followed by the
C     lower treshold, upper treshold and minimal gap (to prevent cutting through
C     quasi-degenerate orbitals).
C
C*****************************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
#include "dcborb.h"
#include "dcbgen.h"
C
      CHARACTER A*(*)
      DIMENSION ILIST(NLIST)
      LOGICAL TOBE
  
      real(8), allocatable :: eig(:)
      integer, allocatable :: ibg(:)
      integer, allocatable :: itm(:)

      CALL QENTER('NUMLST')
C
C
C     Distinguish between the two branches
C
      NA = LEN(A)
!     IF (NA.GT.80) NA = 80
      IF (INDEX(A(1:NA),'energy') .EQ. 0) THEN
         allocate(itm(IRNGE2-IRNGE1+1))
         CALL NUMLS1(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,ITM)
         deallocate(itm)
      ELSE
         allocate(eig(norbt))
         allocate(ibg(norbt))
C
C        Get eigenvalues from DFCOEF
C
         INQUIRE(FILE='DFCOEF',EXIST=TOBE)
         IF(.NOT.TOBE.OR.I_DCBORB_SET.NE.1) THEN
C
C           DFCOEF does not (yet) exist, put IP = - 666, which makes 
C           NUMLS2 return with an empty list
C           Another showstopper occurs if we call this before the 
C           orbital common block is set (this info is currently needed
C           to correctly process the orbital energy list), this is
C           why we check for I_DCBORB_SET as well.
C
C           miro(14.3.'13): changed EIG in the NUMLS2 parameter list, which is of zero length, 
C           to DUMMY to satisfy run with runtime check flags
C
            IP = - 666
            CALL NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,DUMMY)
         ELSE
            CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PCMOUT')
            CALL REACMO(LUCOEF,'DFCOEF',DUM,EIG,IBG,DUM,4)
            CLOSE(LUCOEF,STATUS='KEEP')
            IOFF = IORB(IFRP)+NPSH(IFRP)
            CALL NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,EIG(1+IOFF))
         ENDIF
         deallocate(eig)
         deallocate(ibg)
      ENDIF

      CALL QEXIT('NUMLST')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls1 */
      SUBROUTINE NUMLS1(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,ITMPAR)
C*****************************************************************************
C
C     Reads the numbers in the string A which has
C     the format : '3..6, 8, -27..2, 123, 14..15'.
C     
C     !! Number 0 (Zero) NOT counted or stored !!
C
C     Input: 
C           IP > 0: Give list of numbers in array
C                         ILIST(NLIST).
C           IP = 0: Only IP as OUTPUT
C           IP < 0: Give list of numbers in array
C                         ILIST(NLIST) - but the order is maintained,
C                         ie, 1-8,10,9  .NE.  1-10
C
C     Output:
C           IP : number of elements in the list A.
C                Duplicates only counted once.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     ITMPAR   : Temporary storage
C     If the string starts with 'all', all elements in the range enters.w
C
C                           Jon K. Laerdahl 30.12.96
C     Minor polish T.Saue Jan 7 1997
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      PARAMETER (MAXLB = 400)
C
      CHARACTER A*(*),B*(MAXLB)
      CHARACTER*1 TEGN
      DIMENSION ITMPAR(IRNGE1:IRNGE2),ILIST(NLIST)

C
C     Initialize
C
      NUM1 = 0
      NUM2 = 0
      INO  = 0
      NA   = LEN(A)
!     IF (NA.GT.80) NA = 80
      IF(NA.GE.3) THEN
        IF( INDEX(A(1:NA),'all') .GT. 0) THEN
          IC = 0
          IF(IP.EQ.0) THEN
            DO 10 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 10
              IC = IC + 1
 10         CONTINUE
          ELSE
            DO 11 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 11
              IC = IC + 1
              IF (IC.GT.NLIST) GOTO 930
              ILIST(IC) = I
 11         CONTINUE
          ENDIF
          IP = IC
          RETURN
        ENDIF
      ENDIF
      IPOINT = ICHAR('.')
      ISPACE = ICHAR(' ')
      ISIGN  = ICHAR('-')
      ICOMMA = ICHAR(',')
      INULL  = ICHAR('0')
      ININE  = ICHAR('9')
      IOOH   = ICHAR('o')
      CALL IZERO(ITMPAR,IRNGE2-IRNGE1+1)
      B = ' '
C
C     Take out blanks
C
      NB = 0
      DO 20 I = 1,NA
         TEGN = A(I:I)
         IF (TEGN.EQ.' ') GOTO 20
         NB = NB + 1
         IF(NB.GT.MAXLB) GOTO 940
         B(NB:NB) = TEGN
 20   CONTINUE
      I = 1
 100  CONTINUE
      NEG1 = 1
      NEG2 = 1
      INT = ICHAR(B(I:I))
      IF (INT.EQ.ISIGN) THEN
         NEG1 = -1
         I = I+1
         INT = ICHAR(B(I:I))
      ENDIF
      IF ((INT.GE.INULL).AND.(INT.LE.ININE) .OR.
     &    (INT.EQ.IOOH) ) THEN
         IF (INT .EQ. IOOH) THEN
            I = I + 1
            TEGN = B(I:I)
            IF (TEGN .NE. 'o') GOTO 950
            IF (NEG1 .EQ. -1) THEN
               NUM1 = -IRNGE1
            ELSE
               NUM1 = IRNGE2
            END IF
            I = I + 1
         ELSE
C           ... found a number
            NUM1 = INT-INULL
            CALL READNM(I,NUM1,B)
            NUM1 = NUM1*NEG1
         END IF
         TEGN = B(I:I)
         IF (TEGN.EQ.'.') THEN
C        ... checking for range
            I = I+1
            TEGN = B(I:I)
            IF (TEGN.NE.'.') THEN
               GOTO 900
            ELSE
C           ... find end of range
               I = I+1
               INT = ICHAR(B(I:I))
               IF (INT.EQ.ISIGN) THEN
                  IF (NEG1.GT.0) GOTO 900
                  NEG2 = -1
                  I = I+1
                  INT = ICHAR(B(I:I))
               ENDIF
               IF ((INT.GE.INULL).AND.(INT.LE.ININE)) THEN
C              ... read end of range
                  NUM2 = INT-INULL
                  CALL READNM(I,NUM2,B)
                  NUM2 = NUM2*NEG2
               ELSE IF (INT .EQ. IOOH) THEN
                  I = I + 1
                  TEGN = B(I:I)
                  IF (TEGN .NE. 'o') GOTO 950
                  IF (NEG2 .EQ. -1) THEN
                     NUM2 = -IRNGE1
                  ELSE
                     NUM2 = IRNGE2
                  END IF
                  I = I + 1
               ELSE 
                  GOTO 900
               ENDIF
            ENDIF
            IF (NUM1.GT.NUM2) GOTO 900
            IF ((NUM1.GE.IRNGE1).AND.(NUM2.LE.IRNGE2)) THEN
               DO 200 J=NUM1,NUM2
                  INO = INO + 1
                  ITMPAR(J) = INO
 200           CONTINUE
            ELSE
               WRITE(LUPRI,'(A,I5,A,I5)') 'Range:', NUM1,'..',NUM2
               GOTO 920
            ENDIF
         ELSEIF ((TEGN.EQ.',').OR.(TEGN.EQ.' ')) THEN
C           ... single number
            IF ((NUM1.GE.IRNGE1).AND.(NUM1.LE.IRNGE2)) THEN
               INO = INO + 1
               ITMPAR(NUM1) = INO
            ELSE
               WRITE(LUPRI,'(A,I5)') 'Value:', NUM1
               GOTO 920
            ENDIF
         ELSE 
            GOTO 910
         ENDIF
      ELSEIF ((INT.NE.ISPACE).AND.(INT.NE.IPOINT)) THEN
         GOTO 910
      ENDIF
      I = I+1
      IF (I.LE.NB) GOTO 100
C     
C     Count number of strings
C
      ICNT = 0
      DO 300 I=IRNGE1,IRNGE2
C
C     Do NOT count ZERO...
C
         IF (I.EQ.0) GOTO 300
         IF (ITMPAR(I).NE.0) ICNT = ICNT+1
 300  CONTINUE
C
      IF(IP.NE.0) THEN
         CALL IZERO(ILIST,NLIST)
         IC = 0
         DO 400 I=IRNGE1,IRNGE2
C
C     Do NOT store ZERO...
C
            IF (I.EQ.0) GOTO 400
            IF (ITMPAR(I).NE.0) THEN
               IC = IC + 1
               IF (IC.GT.NLIST) GOTO 930
               IF (IP.GT.0) THEN
                  ILIST(IC) = I
               ELSE
                  ILIST(ITMPAR(I)) = I
               END IF
            ENDIF
 400     CONTINUE
      ENDIF
      IP = ICNT
C
      RETURN

 900  CONTINUE
      WRITE(LUPRI,'(//A,I3//3A/9X,200A1)')
     &     'NUMLS1: Wrong syntax in input at character ',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Wrong Syntax in Input')
 910  CONTINUE
      WRITE(LUPRI,'(/A,A1,A,I4,A,I3//3A/9X,200A1)')
     &   'NUMLS1: Wrong Character in Input: ',CHAR(INT),
     &   ' ( ASCII code',INT,') at character',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Wrong Character in Input')
 920  CONTINUE
      WRITE(LUPRI,'(/A,2I6//3A/9X,200A1)')
     &     'Out of Range in NUMLS1. Supported range:',IRNGE1,IRNGE2,
     &     'Entering string: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      WRITE(LUPRI,'(2X,A,/)')  
     & '...upper value of this range is defined by parameter'
      CALL FLSHFO(LUPRI) 
      CALL QUIT('NUMLS1: Out of Range!')
 930  CONTINUE
      WRITE(LUPRI,'(/A,I7//3A/9X,200A1)')
     &     'Out of Range in NUMLS1. NLIST:',NLIST,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Out of Range in ITMPAR') 
 940  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS1: Character string beyond MAXLB.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: Character string beyond MAXLB') 
 950  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS1: oo must be two consecutive characters.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS1: oo must be two consecutive characters.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls2 */
      SUBROUTINE NUMLS2(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,EIG)
C*****************************************************************************
C
C     Gives the list of orbitals that fulfill the energy criteria
C
C     Input: 
C           IP = 0: Only IP as OUTPUT
C           IP > 0: Give also list
C
C     Output:
C           IP : number of elements in the list A.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     EIG   : List of eigenvalues
C
C     L. Visscher, May 2001
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER A*(*)
      DIMENSION EIG(*)
      DIMENSION ILIST(NLIST)
      LOGICAL CHMIN,CHMAX
C     
C     Determine begin and end of the part of the string that contains
C     the energy values
C
      NA = LEN(A)
!     IF (NA.GT.80) NA = 80
      MA = INDEX(A(1:NA),'energy') + 6
C     
      READ (A(MA:NA),*,END=900,ERR=900) ELOW,EHIGH,EGAP
C
      IF (IP.EQ.-666) THEN
         IP = 0
         RETURN
      ENDIF
C
      WRITE (LUPRI,1000,IOSTAT=IOS) ELOW,EHIGH,EGAP
      IF (IOS.NE.0) THEN
        WRITE(LUPRI,*) "Error in formatting this WRITE statement !"
        WRITE(LUPRI,*) 
     & "Check numbers formats in order not to get *** in output !"
        WRITE(LUPRI,*) "number, ELOW=",ELOW
        WRITE(LUPRI,*) "number, EHIGH=",EHIGH
        WRITE(LUPRI,*) "number, EGAP=",EGAP
        CALL FLSHFO(LUPRI)
        !CALL QUIT('Error in formatting write statement !')
      ENDIF

C
C     Check that the specified ELOW and EHIGH are separated by a 
C     large enough gap from the next lower/higher energies
C
   1  EMIN = 1.E16
      EMAX = - 1.E16
C
C     Set start and end of search region
C
      IF (IRNGE1.LT.0) THEN
         ISTART = IRNGE1 + 1
      ELSEIF (IRNGE1.GT.0) THEN
         ISTART = IRNGE1
      ELSE
         WRITE(LUPRI,*) ' IRNGE1:',IRNGE1
         CALL QUIT ("Error in energy-based vector selection (NUMLST)")
      ENDIF
C
      IF (IRNGE2.LT.0) THEN
         IEND = IRNGE2 + 1
      ELSEIF (IRNGE2.GT.0) THEN
         IEND = IRNGE2
      ELSE
         WRITE(LUPRI,*) 'IRNGE2:',IRNGE2
         CALL QUIT ("Error in energy-based vector selection (NUMLST)")
      ENDIF
C
      DO I = ISTART, IEND
        EI = EIG(I)
        IF (EI.GE.ELOW.AND.EI.LE.EHIGH) THEN
           IF (EMIN.GT.EI) EMIN = EI
           IF (EMAX.LT.EI) EMAX = EI
        ENDIF
      ENDDO
C
      CHMIN = .FALSE.
      CHMAX = .FALSE.
C
      DO I = ISTART, IEND
        EI = EIG(I)
        IF (EI.LT.EMIN.AND.EI.GE.EMIN-EGAP) THEN
           IF (ELOW.GT.EI) ELOW = EI
           CHMIN = .TRUE.
        ENDIF
        IF (EI.GT.EMAX.AND.EI.LE.EMAX+EGAP) THEN
           IF (EHIGH.LT.EI) EHIGH = EI
           CHMAX = .TRUE.
        ENDIF
      ENDDO
C
      IF (CHMIN) WRITE (LUPRI,1001) 'minimum',ELOW
      IF (CHMAX) WRITE (LUPRI,1001) 'maximum',EHIGH
C
      IF (CHMIN.OR.CHMAX) GOTO 1
C
C     Determine the orbitals that fall within this energy range
C
      ICNT = 0
      DO I = ISTART, IEND
        EI = EIG(I)
        IF (EI.GE.ELOW.AND.EI.LE.EHIGH) THEN
           ICNT = ICNT + 1
           IF (IP.NE.0) ILIST(ICNT) = I
        ENDIF
      ENDDO
C
      IP = ICNT
      RETURN
C
 900  CONTINUE
      WRITE(LUPRI,'(//A)')
     & 'NUMLS2: Error reading energy criteria, the string is ',A(1:NA)
      CALL QUIT('NUMLS2: Error in Input')
 1000 FORMAT (5X,'Energy selection of active orbitals :',
     & F12.2,' < Eps. < ',F12.2,
     & ' with a mininum gap of ',F10.4,' au.')
 1001 FORMAT (/,' Changed ',A7,' energy for active orbital range to ',
     & F10.4,' because of near-degeneracies.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls4 */
      SUBROUTINE NUMLS4(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,EIG)
C*****************************************************************************
C
C     Gives the list of NOs that fulfill the NO occupancy selection 
C     criteria to use in subsequent correlation calculation
C
C     Input: 
C           IP = 0: Only IP as OUTPUT
C           IP > 0: Give also list
C           IP < 0: Exit after printing the criteria 
C
C     Output:
C           IP : number of elements in the list A.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C     EIG   : List of eigenvalues
C
C     S. Knecht, April 2008 
C
C     based on NUMLS2 by L. Visscher.
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
      CHARACTER A*(*)
      DIMENSION EIG(*)
      DIMENSION ILIST(NLIST)
      LOGICAL CHMIN,CHMAX
C     
C     Determine begin and end of the part of the string that contains
C     the energy values
C
      NA = LEN(A)
!     IF (NA.GT.80) NA = 80
      MA = INDEX(A(1:NA),'NO-occ') + 6
C     
      READ (A(MA:NA),*,END=900,ERR=900) OHIGH,OLOW,OGAP
C
      IF (IP.EQ.-666) THEN
         IP = 0
         RETURN
      ENDIF
C
      WRITE (LUPRI,1000) OHIGH,OLOW,OGAP*100
C
C     Check that the specified OLOW and OHIGH are separated by a 
C     large enough gap from the next lower/higher energies
C
   1  OMIN = -1.0D-01
      OMAX = 2.0D0
C
C     Set start and end of search region
C
      IF (IRNGE1.LT.0) THEN
         ISTART = IRNGE1 + 1
      ELSEIF (IRNGE1.GT.0) THEN
         ISTART = IRNGE1
      ELSE
         CALL QUIT ("Error in occupation-based vector 
     &               selection; IRNGE1 (NUMLST4)")
      ENDIF
C
      IF (IRNGE2.LT.0) THEN
         IEND = IRNGE2 + 1
      ELSEIF (IRNGE2.GT.0) THEN
         IEND = IRNGE2
      ELSE
         CALL QUIT ("Error in occupation-based vector 
     &               selection; IRNGE2 (NUMLST4)")
      ENDIF
C
      DO I = ISTART, IEND
        OI = EIG(I)
        IF (OI.GE.OLOW.AND.OI.LE.OHIGH) THEN
           IF (OMIN.GT.OI) OMIN = OI
           IF (OMAX.LT.OI) OMAX = OI
        ENDIF
      ENDDO
C
      CHMIN = .FALSE.
      CHMAX = .FALSE.
C
      DO I = ISTART, IEND
        OI = EIG(I)
        IF (OI.LT.OMIN.AND.OI.GE.OMIN-(OGAP*OLOW)) THEN
           IF (OLOW.GT.OI) OLOW = OI
           CHMIN = .TRUE.
        ENDIF
        IF (OI.GT.OMAX.AND.OI.LE.OMAX+(OGAP*OHIGH)) THEN
           IF (OHIGH.LT.OI) OHIGH = OI
           CHMAX = .TRUE.
        ENDIF
      ENDDO
C
      IF (CHMIN) WRITE (LUPRI,1001) 'minimum',OLOW
      IF (CHMAX) WRITE (LUPRI,1001) 'maximum',OHIGH
C
      IF (CHMIN.OR.CHMAX) GOTO 1
C
C     Determine the NOs that fall within this occupation range
C
      ICNT = 0
      DO I = ISTART, IEND
        OI = EIG(I)
        IF (OI.GE.OLOW.AND.OI.LE.OHIGH) THEN
           ICNT = ICNT + 1
           IF (IP.NE.0) ILIST(ICNT) = I
        ENDIF
      ENDDO
C
      IP = ICNT
      RETURN
C
 900  CONTINUE
      WRITE(LUPRI,'(//A)')
     &' NUMLS4: Error reading selection criteria, the string is ',
     &  A(1:NA)
      CALL Abend2('NUMLS4: Error in Input')
 1000 FORMAT (/3X,'Occ. sel. of active natural orbitals :',
     & F10.6,' < occ. < ',F10.6,' with a relative gap of ',F6.2,'%')
 1001 FORMAT (/' Changed ',A7,' occupation for active NO range to ',
     & F10.4,' because of near-degeneracies.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck NOSELPRI */
      SUBROUTINE NOSELPRI(NINDEX,I,NUMTOT)
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcborb.h"
C
      INTEGER NINDEX(*)
C
      WRITE(LUPRI,'(/3X,A,A3)')'* Fermion ircop ',FREP(I)
      WRITE(LUPRI,'(3X,A,I4/)')'* active natural spinors:',NUMTOT
      WRITE(LUPRI,'(3X,8I5)') (NINDEX(K),K=1,NUMTOT)
C
      END 
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck readnm */
      SUBROUTINE READNM(I,NUM,A)
      CHARACTER A*(*),TEGN*1
C
      INULL  = ICHAR('0')
      ININE  = ICHAR('9')
 10   I = I+1
      TEGN = A(I:I)
      INT = ICHAR(TEGN)
      IF ((INT.GE.INULL).AND.(INT.LE.ININE)) THEN
         NUM = NUM*10+(INT-INULL)
         GO TO 10
      ELSE
         RETURN
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C   /* Deck wrtiar */
      SUBROUTINE WRTIAR(ILIST,IP)
#include "implicit.h"
#include "priunit.h"
      DIMENSION ILIST(IP)
C
      WRITE(LUPRI,'(10I5)') (ILIST(I), I=1,IP)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck EIGEXT */
      SUBROUTINE EIGEXT(EIG,IORBCL,ITYP,WORK,KFREE,LFREE)
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
C
      DIMENSION EIG(*),IORBCL(4,NORBT),WORK(*)
      KFRSAV = KFREE
      CALL MEMGET('REAL',KCMO,NCMOTQ,WORK,KFREE,LFREE)
      CALL MEMGET('REAL',KEIG,NORBT ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIBEIG,NORBT ,WORK,KFREE,LFREE)
      LUBUF = 22
      CALL OPNFIL(LUBUF,'DFCOEF','OLD','EIGEXT')
      CALL REACMO(LUBUF,'DFCOEF',WORK(KCMO),WORK(KEIG),WORK(KIBEIG),
     &            TOTERG,14)
      CLOSE(LUBUF,STATUS='KEEP')
      DO I = 1,NORBT
      IF(IORBCL(1,I).EQ.ITYP) THEN
        IND      = IORBCL(2,I) 
        EIG(IND) = WORK(KEIG+I-1)
      ENDIF
      ENDDO
      CALL MEMREL('EIGEXT',WORK,KFRSAV,KFRSAV,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adacmo */
      SUBROUTINE ADACMO(NB,NZ,CF1,NO1,NP1,NE1,CF2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt coefficients CF2 to current dimensions CF1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CF1(NB,NO1,NZ),CF2(NB,NO2,NZ)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        NN = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)        
      ENDDO
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      NN = NB*NC
      DO IZ = 1,NZ
        CALL DCOPY(NN,CF2(1,IV2,IZ),1,CF1(1,IV1,IZ),1)        
      ENDDO
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        NN  = NB*NR
        DO IZ = 1,NZ
          CALL DZERO(CF1(1,IV1,IZ),NN)
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adaeig */
      SUBROUTINE ADAEIG(EIG1,NO1,NP1,NE1,EIG2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt eigenvalues EIG2 to current dimensions EIG1
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION EIG1(NO1),EIG2(NO2)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        CALL DZERO(EIG1(IV1),NR)
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)        
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL DCOPY(NC,EIG2(IV2),1,EIG1(IV1),1)        
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        CALL DZERO(EIG1(IV1),NR)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Adaibo */
      SUBROUTINE ADAIBO(IBO1,NO1,NP1,NE1,IBO2,NO2,NP2,NE2)
C***********************************************************************
C
C     Adapt boson irrep array IBO2 to current dimensions IBO1
C     If negative-energy solutions are missing,
C        then initial values are zeroed out.
C     If positive-energy solutions are missing,
C        then final values are zeroed out.      
C     Written by T.Saue
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION IBO1(NO1),IBO2(NO2)
C     Positronic part
      IV1 = 1
      IV2 = 1
      NC = MIN(NP1,NP2)
      NR = NP1-NP2
      IF(NR.GT.0) THEN
        CALL IZERO(IBO1(IV1),NR)
        IV1 = IV1 + NR
      ELSE
        IV2 = 1-NR
      ENDIF
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)        
C     Electronic part
      IV1 = NP1 + 1
      IV2 = NP2 + 1
      NC = MIN(NE1,NE2)
      NR = NE1-NE2
      CALL ICOPY(NC,IBO2(IV2),1,IBO1(IV1),1)        
      IF(NR.GT.0) THEN
        IV1 = IV1 + NC
        CALL IZERO(IBO1(IV1),NR)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Nucbas */
      SUBROUTINE NUCBAS(NUCORB,AA,IPRINT)
C***********************************************************************
C
C     Extract basis information for all centers
C
C     Written by T.Saue March 12 2001
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0)
C
#include "maxaqn.h"
#include "mxcent.h"
#include "maxorb.h"
#include "aovec.h"
#include "dummy.h"
C
#include "nuclei.h"
#include "ccom.h"
#include "shells.h"
#include "primit.h"
      CHARACTER SPDCAR*1
      DIMENSION NUCORB(NHTYP,2,NUCIND),AA(2,NHTYP,2,NUCIND)
C
C     Initialize
C
      NDIM = 2*NHTYP*NUCIND
      CALL IZERO(NUCORB,NDIM)
C
      JCENT = 0
      JPRIM = -1
      JC = -1
      DO ISHELL = 1,KMAX
        ICENT = NCENT(ISHELL)
        IF(ICENT.NE.JCENT) THEN
          JCENT = ICENT
          JLVAL = 0
        ENDIF
        IC = LCLASS(ISHELL)
        IF(IC.NE.JC) THEN
          JC    = IC
          JLVAL = 0
        ENDIF
        ILVAL = NHKT(ISHELL)
        IF(ILVAL.NE.JLVAL) THEN
          JLVAL = ILVAL
          NUCORB(ILVAL,IC,ICENT) = 0
          AA(1,ILVAL,IC,ICENT)=D0
          AA(2,ILVAL,IC,ICENT)=DUMMY
        ENDIF
        NUCORB(ILVAL,IC,ICENT)=NUCORB(ILVAL,IC,ICENT)+1
        IPRIM = JSTRT(ISHELL)
        IF(IPRIM.NE.JPRIM) THEN
          JPRIM = IPRIM
          NPRIM = NUCO(ISHELL)
          DO IEXP = 1,NPRIM
            A=PRIEXP(IPRIM+IEXP)
            AA(1,ILVAL,IC,ICENT)=MAX(AA(1,ILVAL,IC,ICENT),A)
            AA(2,ILVAL,IC,ICENT)=MIN(AA(2,ILVAL,IC,ICENT),A)
          ENDDO
        ENDIF
      ENDDO
C
      IF(IPRINT.GE.2) THEN
         CALL HEADER('NUCORB:Basis set information:',-1)
        DO I = 1,NUCIND
          WRITE(LUPRI,'(/A,A4,A/)') '*** Center: ',NAMN(I),' ***'
          WRITE(LUPRI,'(2X,A)') '* Large components:'
          IC = 1
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H :',AA(1,LL,IC,I),
     &      'Alpha_L :',AA(2,LL,IC,I)
          ENDIF
          ENDDO
          WRITE(LUPRI,'(2X,A)') '* Small components:'
          IC = 2
          DO LL = 1,NHTYP
          IF(NUCORB(LL,IC,I).GT.0) THEN
            L=LL-1
            WRITE(LUPRI,'(3X,A1,A,I6,2(3X,A,E12.5))')
     &      SPDCAR(L),'-orbitals:',NUCORB(LL,IC,I),
     &      'Alpha_H: ',AA(1,LL,IC,I),
     &      'Alpha_L: ',AA(2,LL,IC,I)
          ENDIF
          ENDDO
        ENDDO
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck acmout */
      SUBROUTINE ACMOUT()
C*****************************************************************************
C
C     Read coefficients from unformatted DFCOEF, 
C     take out symmetry and write to unformatted DFACMO 
C
C     There is a catch:
C       In HERMIT the basis set is generally symmetry-adapted. 
C       Routines like WTSOAO will transform from this SO-basis
C       to an AO-basis, where the symmetry adaption is taken out.
C       In this AO-basis looping is over
C        - atomic type
C         - symmetry independent center
C          - shell (exponents/contracted functions)
C           - degeneracy of center
C             - component
C       [NOTE that in DALTON the looping over degeneracy and components
C        is reversed !]
C
C       A consequence of this is that the ordering of the AO-basis may change
C       if the symmetry of the system is lowered, that is one does not use the
C       full symmetry.
C
C       In order to make a more universal AO-basis, ACMOUT will modify the
C       basis set ordering such that the looping goes as
C        - atomic type
C         - symmetry independent center
C           - degeneracy of center
C             - shell
C               - component
C
C       This is not completely foolproof, so the import of DFACMO to another
C       symmetry has to be handled with care !
C       
C     Written by T.Saue May 26 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "dcbpsi.h"
#include "dcbgen.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcblab.h"
#include "dgroup.h"
#include "nuclei.h"
      LOGICAL TOBE, TOBEK, FNDLAB, MCORB, MJAVA
      CHARACTER TEXT*74
      DIMENSION IDIM(3)
      real(8), allocatable :: CSO(:)
      real(8), allocatable :: EIG(:)
      integer, allocatable :: IBEIG(:)
      real(8), allocatable :: CAO(:)
C
C     First check if file KRMCSCF or DFCOEF is there
C
      MCORB = .FALSE.
      MJAVA = .FALSE.
      INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF((.NOT.TOBE).and.(.NOT.TOBEK)) THEN
        WRITE(LUPRI,'(A)') 
     &  'ACMOUT: Unformatted file DFCOEF/KRMCSCF not found !'
        CALL QUIT('ACMOUT: Unformatted file DFCOEF/KRMCSCF not found !')
      ENDIF
C     Allocate arrays
      allocate(CSO(NCMOTQ))
      allocate(EIG(NORBT))
      allocate(IBEIG(NORBT))
      allocate(CAO(N2BBASX*NZC1))
C
C     Get coefficients from KRMCSCF/DFCOEF --> stored in CSO
C     =======================================================
C
      IF(TOBEK)THEN
        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','ACMOU1')
        REWIND LUKRMC
C.......CI or MP2 natural orbitals
        MCORB = FNDLAB('MCCINATO',LUKRMC)
        REWIND LUKRMC
C.......boson/linear symmetry
        MJAVA = FNDLAB('MJVEC   ',LUKRMC)
        IF(MJAVA)THEN
          CALL IREAKRMC(LUKRMC,'MJVEC   ',IBEIG,NORBT)
        ELSE
          CALL IZERO(IBEIG,NORBT)
        END IF
        IF(.NOT.MCORB)THEN
          JRDMO = -1
          CALL RREADMO(CSO,JRDMO,1,LUKRMC)
!         ... set to zero
          TOTERG = 0.0D0
          CALL DZERO(EIG,NORBT)
          IREADMO = 1
        ELSE
!
!         read CI/MCSCF natural orbitals and occupancies
          JRDMO = -1
          CALL RREADMO(CSO,JRDMO,4,LUKRMC)
          CALL REAKRMC(LUKRMC,'MCNATOCC',EIG,NORBT)
!
!         temporary fix for use of occupancies in Mulliken population
!         analysis - set energy to 137.0D0 - SK Feb 2009
          TOTERG = 137.0D0
          IREADMO = 3
        END IF
        CLOSE (LUKRMC,STATUS='KEEP')
      ELSE
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','ACMOUT')
        CALL REACMO(LUCOEF,'DFCOEF',CSO,EIG,IBEIG,TOTERG,14)
        CLOSE(LUCOEF,STATUS='KEEP')
        IREADMO = 2
      END IF
C.....Transform SO-coefficients to C1
      CALL C1COEF(CSO,CAO,NORB,0,NORBT,1)
      deallocate(CSO)
C
C     Get correct column ordering
C     ===========================
C 
      NDIM = NTBAS(0)*NORBT*NZC1
C     ... only if coefficients are read from DFCOEF
      IF(IREADMO.EQ.2) THEN
C.......if more than one fermion ircop merge orbital classes
        CALL IZERO(IBEIG,NORBT)
        IF(NFSYM.EQ.2) CALL ACMOU3(EIG,CAO,NTBAS(0),NORBT)
      END IF
C
C     Write coefficients etc. to file
C     ===============================
C
      CALL OPNFIL(LUCOEF,'DFACMO','UNKNOWN','ACMOUT')
      DO I = 1,2
        IDIM(I) = 0
      ENDDO
      DO IFRP = 1,NFSYM
        IDIM(1) = IDIM(1) + NPSH(IFRP)  
        IDIM(2) = IDIM(2) + NESH(IFRP)  
      ENDDO
      IDIM(3) = NTBAS(0)
      TEXT(1:50) = TITLE
      CALL GTINFO(TEXT(51:74))
      CALL NEWLAB('INFO    ',LUCOEF,LUPRI)
      NSYM=1
      NZBUF=4
      WRITE(LUCOEF) TEXT,NSYM,NZBUF,(IDIM(I),I = 1,3),TOTERG
      CALL WRICMO_coefficients(LUCOEF,CAO,NTBAS(0)*NORBT*4)
      CALL WRICMO_eigenvalues(LUCOEF,EIG,NORBT)      
      CALL WRICMO_supersymmetry(LUCOEF,'SUPERSYM',IBEIG,NORBT)
!     Add additional information about basis set to create a simple interface file
!     for correlated calculations (like exacorr)
      CALL NEWLAB('BASINFO ',LUCOEF,LUPRI)                              
      call append_basisinfo_to_dfcoef (LUCOEF)
      CALL NEWLAB('EOFLABEL',LUCOEF,LUPRI)                                    
      CLOSE(LUCOEF,STATUS='KEEP')
C
      CALL PRSYMB(LUPRI,'=',75,0)
      IF( IREADMO .eq. 1 )THEN
        WRITE(LUPRI,'(A)')
     &      '* ACMOU1: MCSCF orb. coefficients read from KRMCSCF',
     &      '          and written to unformatted DFACMO (no symmetry)'
      ELSE IF( IREADMO .eq. 2)THEN
        WRITE(LUPRI,'(A)')
     &      '* ACMOU1: Coefficients read from unformatted DFCOEF ',
     &      '          and written to unformatted DFACMO (no symmetry)'
      ELSE IF( IREADMO .eq. 3)THEN
        WRITE(LUPRI,'(A)')
     &      '* ACMOU1: MC natural orb. coefficients read from KRMCSCF',
     &      '          and written to unformatted DFACMO (no symmetry)'
      END IF
      CALL PRSYMB(LUPRI,'=',75,0)
C
C     De-allocate arrays
      deallocate(EIG)
      deallocate(IBEIG)
      deallocate(CAO)
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck acmou2 */
      SUBROUTINE ACMOU2(CMO1,NMO1,CMO2,NMO2,NBAS,IBUF,NZ)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      DIMENSION CMO1(NBAS,NMO1,NZ),CMO2(NBAS,NMO2,NZ),IBUF(*)
C
      DO IZ = 1,NZ
        DO J = 1,NMO2
          DO I = 1,NBAS
            CMO1(I,J,IZ) = CMO2(IBUF(I),J,IZ)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck acmou3 */
      SUBROUTINE ACMOU3(EIG,CMO,NR,NC)
C***********************************************************************
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbdhf.h"
#include "dgroup.h"
      DIMENSION CMO1(NR,NC,NZC1),CMO2(NR,NC,NZC1),
     &          EIG1(NC),EIG2(NC),IVEC(2)     
      real(8), allocatable :: CBUF(:),EBUF(:)
      NDIM=NR*NC*NZC1
      allocate(CBUF(NDIM))
      allocate(EBUF(NC))
      CALL DCOPY(NDIM,CMO,1,CBUF,1)
      CALL DCOPY(NC  ,EIG,1,EBUF,1)
C     Initialize offsets; they will be updated inside ACMOU4
      IVEC(1) = 0
      IVEC(2) = NORB(1) 
      JOFF    = 0
C.....Negative-energy solutions
      CALL ACMOU4(CBUF,EBUF,CMO,EIG,NPSH,IVEC,
     &            JOFF,NFSYM,NR,NC,NZC1)
C.....Inactive orbitals      
      CALL ACMOU4(CBUF,EBUF,CMO,EIG,NISH,IVEC,JOFF,
     &            NFSYM,NR,NC,NZC1)
C.....Active orbitals
      DO IOPEN = 1,NOPEN
        CALL ACMOU4(CBUF,EBUF,CMO,EIG,NACSH(1,IOPEN),IVEC,JOFF,
     &              NFSYM,NR,NC,NZC1)
      ENDDO
C.....Virtual orbitals
      CALL ACMOU4(CBUF,EBUF,CMO,EIG,NSSH,IVEC,JOFF,
     &            NFSYM,NR,NC,NZC1)
C
      deallocate(CBUF)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck acmou4 */
      SUBROUTINE ACMOU4(CMO1,EIG1,CMO2,EIG2,NVEC,IVEC,JVEC,
     &                  NFSYM,NBAS,NORB,NZC1)
C***********************************************************************
C
C     Merge orbital coefficients/energies from different fermion ircops
C     of the same orbital class and order them by energy
C
C     Written by T. Saue in Nha Trang Aug 29 3014
C
C***********************************************************************
      implicit none
      integer, intent(in)   :: NFSYM,NBAS,NORB,NVEC(NFSYM),NZC1
      real(8), intent(in)   :: EIG1(NORB),CMO1(NBAS,NORB,NZC1)
      integer, intent(inout):: IVEC(NFSYM),JVEC
      real(8), intent(out)  :: EIG2(NORB),CMO2(NBAS,NORB,NZC1)
      integer, allocatable  :: IND(:),IBUF(:)
      real(8), allocatable  :: BUF(:)
      integer               :: NVECT,IFRP,I,II,JJ,IZ
C     Find total number of orbitals in this class
      NVECT = 0
      DO IFRP = 1,NFSYM
        NVECT = NVECT + NVEC(IFRP)
      ENDDO
      IF(NVECT.EQ.0) RETURN
C     Merge orbital eigenvalues into buffer
      allocate(IND(NVECT))
      allocate(IBUF(NVECT))
      allocate(BUF(NVECT))
      II = 0
      DO IFRP = 1,NFSYM
        DO I = 1,NVEC(IFRP)
          II = II + 1
          IND(II) = IVEC(IFRP)+I
          BUF(II) = EIG1(IVEC(IFRP)+I)
        ENDDO
      ENDDO
C     Merge/Sort on energy
      CALL INDEXX(NVECT,BUF,IBUF)
      DO I = 1, NVECT
        JJ = I + JVEC
        II = IND(IBUF(I))
        EIG2(JJ) = EIG1(II)
        DO IZ = 1,NZC1
          CALL DCOPY(NBAS,CMO1(1,II,IZ),1,CMO2(1,JJ,IZ),1)        
        ENDDO
      ENDDO
      DO IFRP = 1,NFSYM
        IVEC(IFRP) = IVEC(IFRP) + NVEC(IFRP)
      ENDDO       
      JVEC = JVEC + NVECT
      deallocate(IND)
      deallocate(IBUF)
      deallocate(BUF)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Reaprj */
      SUBROUTINE REAPRJ(LUPRJ,NP,NR,KRVEC,KBVEC,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Read DFPRJC, file from projection analysis
C     NP(IFRP,1) - number of fragment orbitals in fermion ircop IFRP
C     NP(IFRP,2) - number of MO (and polarization) orbitals 
C                  in fermion ircop IFRP
C     NP(IFRP,0) - nref + npol
C
C     Written by T. Saue, Oct 21 2003
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbprj.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION NR(NREFS,NFSYM),NP(2,0:2),WORK(*)
      NCDIM = 0
      NBDIM = 0
      READ(LUPRJ) NREFS
      DO IFRP = 1,NFSYM
C       Dimension....nref.......nmol......
        READ(LUPRJ) NP(IFRP,1),NP(IFRP,2),(NR(I,IFRP),I=1,NREFS)
        NP(IFRP,0) = NP(IFRP,1)+NP(IFRP,2)
        NCDIM = NCDIM + NFBAS(IFRP,0)*NP(IFRP,0)*NZ
        NBDIM = NBDIM + NP(IFRP,0)*NP(IFRP,2)*NZ
      ENDDO
      CALL MEMGET('REAL',KRVEC,NCDIM,WORK,KFREE,LFREE)       
      CALL MEMGET('REAL',KBVEC,NBDIM,WORK,KFREE,LFREE)       
      CALL READT(LUPRJ,NCDIM,WORK(KRVEC))
      CALL READT(LUPRJ,NBDIM,WORK(KBVEC))
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck numls3 */
      SUBROUTINE NUMLS3(A,ILIST,NLIST,IRNGE1,IRNGE2,IP,NZERO)
C*****************************************************************************
C
C     Reads the numbers in the string A which has
C     the format : '3..6, 8, -27..2, 123, 14..15'.
C     
C     !! Number 0 (Zero) IS counted AND stored, and returned as number
C
C     Input: 
C           IP.EQ.0: Only IP as OUTPUT
C           IP.NE.0: Give list of numbers in array ILIST(NLIST).
C
C     Output:
C           IP : number of elements in the list A.
C
C     IRNGE1..IRNGE2 is range of allowed values in A.
C
C     Based on the routine NUMLS1 written by Jon K. Laerdahl 30.12.96
C     Written by T.Saue 2003
C
C*****************************************************************************
#include "implicit.h"
#include "priunit.h"
C
C
      PARAMETER (MAXLB = 400)
C
      CHARACTER A*(*),B*(MAXLB)
      CHARACTER*1 TEGN
      DIMENSION ILIST(NLIST)
C
C     Initialize
C
      NZERO = 0
      NUM1  = 0
      NUM2  = 0
      INO   = 0
      NA = LEN(A)
      IF(NA.GE.3) THEN
        IF( INDEX(A(1:NA),'all') .GT. 0) THEN
          IC = 0
          IF(IP.EQ.0) THEN
            DO 10 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 10
              IC = IC + 1
 10         CONTINUE
          ELSE
            DO 11 I=IRNGE1,IRNGE2
              IF (I.EQ.0) GOTO 11
              IC = IC + 1
              IF (IC.GT.NLIST) GOTO 930
              ILIST(IC) = I
 11         CONTINUE
          ENDIF
          IP = IC
          RETURN
        ENDIF
      ENDIF
      IPOINT = ICHAR('.')
      ISPACE = ICHAR(' ')
      ISIGN  = ICHAR('-')
      ICOMMA = ICHAR(',')
      INULL  = ICHAR('0')
      ININE  = ICHAR('9')
      IOOH   = ICHAR('o')
      B = ' '
C
C     Take out blanks
C
      NB = 0
      DO 20 I = 1,NA
         TEGN = A(I:I)
         IF (TEGN.EQ.' ') GOTO 20
         NB = NB + 1
         IF(NB.GT.MAXLB) GOTO 940
         B(NB:NB) = TEGN
 20   CONTINUE
      I = 1
 100  CONTINUE
      NEG1 = 1
      NEG2 = 1
      INT = ICHAR(B(I:I))
      IF (INT.EQ.ISIGN) THEN
         NEG1 = -1
         I = I+1
         INT = ICHAR(B(I:I))
      ENDIF
      IF ((INT.GE.INULL).AND.(INT.LE.ININE) .OR.
     &    (INT.EQ.IOOH) ) THEN
         IF (INT .EQ. IOOH) THEN
C           ... "infinity"
            I = I + 1
            TEGN = B(I:I)
            IF (TEGN .NE. 'o') GOTO 950
            IF (NEG1 .EQ. -1) THEN
               NUM1 = -IRNGE1
            ELSE
               NUM1 = IRNGE2
            END IF
            I = I + 1
         ELSE
C           ... found a number
            NUM1 = INT-INULL
            CALL READNM(I,NUM1,B)
            NUM1 = NUM1*NEG1
         END IF
         TEGN = B(I:I)
         IF (TEGN.EQ.'.') THEN
C        ... checking for range
            I = I+1
            TEGN = B(I:I)
            IF (TEGN.NE.'.') THEN
               GOTO 900
            ELSE
C           ... find end of range
               I = I+1
               INT = ICHAR(B(I:I))
               IF (INT.EQ.ISIGN) THEN
                  IF (NEG1.GT.0) GOTO 900
                  NEG2 = -1
                  I = I+1
                  INT = ICHAR(B(I:I))
               ENDIF
               IF ((INT.GE.INULL).AND.(INT.LE.ININE)) THEN
C              ... read end of range
                  NUM2 = INT-INULL
                  CALL READNM(I,NUM2,B)
                  NUM2 = NUM2*NEG2
               ELSE IF (INT .EQ. IOOH) THEN
                  I = I + 1
                  TEGN = B(I:I)
                  IF (TEGN .NE. 'o') GOTO 950
                  IF (NEG2 .EQ. -1) THEN
                     NUM2 = -IRNGE1
                  ELSE
                     NUM2 = IRNGE2
                  END IF
                  I = I + 1
               ELSE 
                  GOTO 900
               ENDIF
            ENDIF
            IF (NUM1.GT.NUM2) GOTO 900
            IF ((NUM1.GE.IRNGE1).AND.(NUM2.LE.IRNGE2)) THEN
               DO 200 J=NUM1,NUM2
                  INO = INO + 1
                  IF(J.EQ.0) NZERO=NZERO+1
                  IF(IP.NE.0) ILIST(INO)=J
 200           CONTINUE
            ELSE
               WRITE(LUPRI,'(A,I5,A,I5)') 'Range:', NUM1,'..',NUM2
               GOTO 920
            ENDIF
         ELSEIF ((TEGN.EQ.',').OR.(TEGN.EQ.' ')) THEN
C           ... single number
            IF ((NUM1.GE.IRNGE1).AND.(NUM1.LE.IRNGE2)) THEN
               INO = INO + 1
               IF(NUM1.EQ.0) NZERO=NZERO+1
               IF(IP.NE.0) ILIST(INO)=NUM1
            ELSE
               WRITE(LUPRI,'(A,I5)') 'Value:', NUM1
               GOTO 920
            ENDIF
         ELSE 
            GOTO 910
         ENDIF
      ELSEIF ((INT.NE.ISPACE).AND.(INT.NE.IPOINT)) THEN
         GOTO 910
      ENDIF
      I = I+1
      IF (I.LE.NB) GOTO 100
      IP = INO
      RETURN
 900  CONTINUE
      WRITE(LUPRI,'(//A,I3//3A/9X,200A1)')
     &     'NUMLS3: Wrong syntax in input at character ',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Wrong Syntax in Input')
 910  CONTINUE
      WRITE(LUPRI,'(/A,A1,A,I4,A,I3//3A/9X,200A1)')
     &   'NUMLS3: Wrong Character in Input: ',CHAR(INT),
     &   ' ( ASCII code',INT,') at character',I,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Wrong Character in Input')
 920  CONTINUE
      WRITE(LUPRI,'(/A,2I6//3A/9X,200A1)')
     &     'Out of Range in NUMLS3. Range:',IRNGE1,IRNGE2,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Out of Range!')
 930  CONTINUE
      WRITE(LUPRI,'(/A,I7//3A/9X,200A1)')
     &     'Out of Range in NUMLS3. NLIST:',NLIST,
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Out of Range in ITMPAR') 
 940  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS3: Character string beyond MAXLB.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: Character string beyond MAXLB') 
 950  CONTINUE
      WRITE(LUPRI,'(/A//3A/9X,200A1)')
     &     'NUMLS3: oo must be two consecutive characters.',
     &     'String: "',B(1:NB),'"',(' ',J=1,I-1),'^'
      CALL QUIT('NUMLS3: oo must be two consecutive characters.')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck cmomi1 */
      SUBROUTINE CMOROT(CMO,EIG,IBEIG)
C***********************************************************************
C
C     Perform Jacobi rotations of user-specified angle on pairs
C     of orbitals
C
C     Written by T. Saue March 1 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D180=180.0D0)
#include "pi.h"
C
#include "dcbgen.h"
#include "dcbrot.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      DIMENSION CMO(NCMOTQ),EIG(NORBT),IBEIG(NORBT)

      integer, allocatable :: jvec(:)
      real(8), allocatable :: buf(:)
C
      X = ROTANG*PI/D180
      COSX = COS(X)
      SINX = SIN(X)      
      DO 10 IFRP = 1,NFSYM
        NVECS = 1
        allocate(jvec(norb(ifrp)))
        CALL NUMLST(VECROT(IFRP),JVEC,NORB(IFRP),
     &             -NPSH(IFRP),NESH(IFRP),IFRP,NVECS)
        CALL ORBCN2(JVEC,NVECS,IFRP,NPVECS,NEVECS)
        IF(MOD(NVECS,2).EQ.1) THEN
          CALL QUIT('CMOROT: Odd number of orbitals to rotate !')
        ENDIF
        IF(NVECS.EQ.0) GOTO 20
        allocate(buf(NFBAS(IFRP,0)))
        CALL JACORB(COSX,SINX,NVECS,JVEC,CMO(ICMOQ(IFRP)+1),
     &              NFBAS(IFRP,0),NORB(IFRP),BUF)
        deallocate(buf)
 20     CONTINUE
        deallocate(jvec)
 10   CONTINUE
C
      RETURN
      END      
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck jacorb */
      SUBROUTINE JACORB(COSX,SINX,NVECS,JVEC,CMO,NBAS,NORB,BUF)
C***********************************************************************
C
C     Perform Jacobi rotation of orbitals
C
C     Written by T. Saue March 1 2005
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
      DIMENSION CMO(NBAS,NORB,NZ),BUF(NBAS),JVEC(NORB)
      DO I=1,NVECS,2
        I1 = JVEC(I)
        I2 = JVEC(I+1)
        DO IZ = 1,NZ
          CALL DROT(NBAS,CMO(1,I1,IZ),1,CMO(1,I2,IZ),1,COSX,SINX)
        ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck preorb */
      SUBROUTINE PREORB()
C***********************************************************************
C
C     Manipulation of DFCOEF prior to wave function calculation
C
C***********************************************************************

      use memory_allocator

#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbreo.h"
#include "dgroup.h"
#include "dcbrot.h"
      LOGICAL TOBE
      real(8), allocatable :: cmo(:)
      real(8), allocatable :: eig(:)
      integer, allocatable :: ibeig(:)
C
C     Check if coefficients are on file
C
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(.NOT.TOBE) GOTO 1000

!     memory allocation
      allocate(cmo(ncmotq))
      allocate(eig(norbt))
      allocate(ibeig(norbt))
C
C     Read coefficients, eigenvalues and irrep identification 
C
      CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PREORB')
      IOPT=14
      IF(ATOMIC) THEN
        CALL ATOMIC_REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)           
      ELSE
        CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)
      ENDIF
C
C     Reorder orbitals
C
      IF (LMOORD) CALL REORD(CMO,EIG,IBEIG,
     &                       IMOORD,IREORD,MXREORD)
C
C     Rotate orbitals
C
      IF(ROTORB) CALL CMOROT(CMO,EIG,IBEIG)
C
C     Write new coefficients
C
      IF(ATOMIC) THEN
        CALL ATOMIC_WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
      ELSE
        CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG,.FALSE.)
      ENDIF   
      CLOSE(LUCOEF,STATUS='KEEP')

!     Memory deallocation
      deallocate(cmo)
      deallocate(eig)
      deallocate(ibeig)
      RETURN
C
 1000 CONTINUE
      WRITE(LUPRI,'(A)') 'PREORB: Coefficient file not found !'
      CALL QUIT('PREORB: Coefficients not found !')
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck pstorb */
      SUBROUTINE PSTORB()
C***********************************************************************
C
C     Manipulation of DFCOEF after a wave function calculation
C
C***********************************************************************
      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
      LOGICAL TOBE, TOBEK
      real(8), allocatable :: cmo(:)
      real(8), allocatable :: eig(:)
      integer, allocatable :: ibeig(:)
C

C
C     check if file KRMCSCF or DFCOEF is present
      INQUIRE(FILE='KRMCSCF',EXIST=TOBEK)
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF((.NOT.TOBE).and.(.NOT.TOBEK)) THEN
        WRITE(LUPRI,'(A)') 
     &  'PSTORB: Coefficient file DFCOEF/KRMCSCF not found !'
         CALL QUIT('PSTORB: Coefficient file not found !')
      ENDIF

!     assign mcscf file higher priority
      if(tobek .and. tobe) tobe = .false.
C
C     Memory allocation
      allocate(cmo(ncmotq))
      allocate(eig(norbt))
      allocate(ibeig(norbt))
C
      call pstor1(cmo,eig,ibeig,tobek)
C
C     Memory deallocation
      deallocate(cmo)
      deallocate(eig)
      deallocate(ibeig)
      RETURN
C
 1000 CONTINUE
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/*Deck Pstor1 */
      SUBROUTINE PSTOR1(CMO,EIG,IBEIG,read_mcscf)
C***********************************************************************
C
C     Post orbital modifications
C     Written by T. Saue Jan 27 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbpsi.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbreo.h"
      DIMENSION CMO(*),EIG(*),IBEIG(*)
      logical, intent(in) :: read_mcscf
      logical             :: mcorb, mcorb1, mcorb2, mjava
      logical             :: fndlab, mcci_natorb, boson_irrep_info

C     Get coefficients from KRMCSCF/DFCOEF
C     ====================================
      if(read_mcscf)then

        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','PSTOR1')
        REWIND LUKRMC
        MCORB1 = .not.FNDLAB('MCCINATO',LUKRMC)
        REWIND LUKRMC
        MCORB2 = .not.FNDLAB('NEWNATOB',LUKRMC)
        REWIND LUKRMC
        MJAVA = FNDLAB('MJVEC   ',LUKRMC)
        mcorb = .false.
        if(MCORB1 .and. MCORB2) mcorb = .true.

!       read mj-values / boson irrep info
        IF(MJAVA)THEN
          CALL IREAKRMC(LUKRMC,'MJVEC   ',IBEIG,NORBT)
        ELSE
          REWIND LUKRMC
          boson_irrep_info = FNDLAB('IBEIG   ',LUKRMC)
          if(boson_irrep_info)then
            CALL IREAKRMC(LUKRMC,'IBEIG   ',IBEIG,NORBT)
          else
            CALL IZERO(IBEIG,NORBT)
          end if
        END IF

!       read mcscf /natural orbitals
        IF(MCORB)THEN
          WRITE(LUPRI,'(/A)')
     &         ' (PSTOR1) Read orbitals from' //
     &                 ' label NEWORB on file KRMCSCF'
          JRDMO = -1
          CALL RREADMO(CMO,JRDMO,1,LUKRMC)
!         ... set to zero
          CALL DZERO(EIG,NORBT)
        ELSE
          REWIND LUKRMC
          mcci_natorb = .false.
          if(FNDLAB('MCCINATO',LUKRMC))then
            WRITE(LUPRI,'(/A)')
     &           ' (PSTOR1) Read orbitals from' //
     &                   ' label MCCINATO on file KRMCSCF'
!           read CI/MCSCF natural orbitals and occupancies
            REWIND LUKRMC
            JRDMO = -1
            CALL RREADMO(CMO,JRDMO,4,LUKRMC)
            CALL REAKRMC(LUKRMC,'MCNATOCC',EIG,NORBT)
            mcci_natorb = .true.
          else
            WRITE(LUPRI,'(/A)')
     &           ' (PSTOR1) Read orbitals from' //
     &                   ' label NEWNATOB on file KRMCSCF'
!           read MP2 natural orbitals and occupancies
            JRDMO = -1
            REWIND LUKRMC
            CALL RREADMO(CMO,JRDMO,3,LUKRMC)
            CALL REAKRMC(LUKRMC,'NEWNATOC',EIG,NORBT)
!           ... set to zero
!           CALL DZERO(EIG,NORBT)
          end if
        END IF

!       remove old MCSCF - we write a new one as this reordering step requires 
!       a recalculation of the reference state vector in any following
!       MCSCF/CI run.
        CLOSE (LUKRMC,STATUS='DELETE')
      ELSE
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','PSTOR1')
        IOPT=14
        IF(ATOMIC) THEN
          CALL ATOMIC_REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)           
        ELSE
          CALL REACMO(LUCOEF,'DFCOEF',CMO,EIG,IBEIG,TOTERG,IOPT)
        ENDIF
      END IF
C
C     Reorder orbitals
C
      IF (LMOORF) CALL REORD(CMO,EIG,IBEIG,IMOORF,IREORF,MXREORD)
C
C     Phase-adjustment of coefficients, if requested...
C
      IF(PHCOEF) THEN
        DO I = 1,NFSYM
          CALL PHATRA(CMO(ICMOQ(I)+1),NFBAS(I,0),NORB(I),NZ)
        ENDDO
      ENDIF
C
C     Write new coefficients
C
      if(read_mcscf)then
        CALL OPNFIL(LUKRMC,'KRMCSCF','UNKNOWN','PSTOR1')
        CALL NEWLAB('SODLABEL',LUKRMC,LUPRI)
        REWIND(LUKRMC)
!       mj-values / boson irrep info
        if(mjava)then
          CALL IWRTKRMC(LUKRMC,'MJVEC   ',IBEIG,NORBT)
        else
          if(boson_irrep_info)then
            CALL IWRTKRMC(LUKRMC,'IBEIG   ',IBEIG,NORBT)
          end if
        end if
!       reordered orbitals
        if(mcorb)then
          CALL WRTKRMC(LUKRMC,'NEWORB  ',CMO,NCMOTQ)
        else
          if(mcci_natorb)then
            CALL WRTKRMC(LUKRMC,'MCCINATO',CMO,NCMOTQ)
            CALL WRTKRMC(LUKRMC,'MCNATOCC',EIG,NORBT)
          else
            CALL WRTKRMC(LUKRMC,'NEWNATOB',CMO,NCMOTQ)
            CALL WRTKRMC(LUKRMC,'NEWNATOC',EIG,NORBT)
          end if
        end if
        close(LUKRMC,status='keep')
      else
        IF(ATOMIC) THEN
          CALL ATOMIC_WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG)
        ELSE
          CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,TOTERG,.FALSE.)
        ENDIF   
        CLOSE(LUCOEF,STATUS='KEEP')
      end if
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Orbnum */
      SUBROUTINE ORBNUM(VECREF,IFRP,KVEC,NSTR,WORK,KFREE,LFREE)
C***********************************************************************
C
C     From an orbital string VECREF for fermion ircop IFRP
C     find the total number of vectors NSTR(IFRP,0) and allocate
C     starting adress KVEC. 
C     
C     On input:
C     NSTR(IFRP,1) : max. number of vectors of positive energy
C     NSTR(IFRP,2) : max. number of vectors of negative energy
C
C
C     On output:
C     NSTR(IFRP,0) : total number of vectors
C     NSTR(IFRP,1) : vectors of positive energy
C     NSTR(IFRP,2) : vectors of negative energy
C     KVEC         : starting address for index array of orbital addresses
C     
C     Written by T. Saue Sep 1 2006
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcborb.h"
      CHARACTER VECREF*(*)
      DIMENSION NSTR(2,0:2),WORK(*)
C
      CALL MEMGET('INTE',KVEC,NORB(IFRP),WORK,KFREE,LFREE)
      NSTR(IFRP,0) = 1
      CALL NUMLST(VECREF,WORK(KVEC),NORB(IFRP),-NSTR(IFRP,2),
     &            NSTR(IFRP,1),IFRP,NSTR(IFRP,0))
      CALL ORBCNT(WORK(KVEC),NSTR(IFRP,0),NPSH(IFRP),NESH(IFRP),
     &            NSTR(IFRP,2),NSTR(IFRP,1))
      NSTR(IFRP,0) = NSTR(IFRP,1) + NSTR(IFRP,2)
      CALL MEMREL('ORBNUM',WORK,KVEC,KVEC,KFREE,LFREE)
      CALL MEMGET('INTE',KVEC,NSTR(IFRP,0),WORK,KFREE,LFREE)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck SelFrag */
      SUBROUTINE SELFRAG(CSEL,NVECS,ESEL,IBSEL,NFRAG,NSTR,FRAGFIL,
     &                   OWNBAS,NUCFRAG,IOPT,KVEC,CMO,EIG,IBEIG,
     &                   IOFF,KRMC_FLG,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Routine for extracting coefficients, eigenvalues and boson irrep
C     information for a set of fragments.
C
C     Returns:
C
C         CSEL - selected coefficient adapted to molecule basis dimensions
C         ESEL - corresponding orbital energies 
C         IBSEL - corresponding irrep information
C         NVECS(IFRP) - number of selected coefficients
C         NFRAG - number of fragments
C         NSTR(IFRP,I,NFRAG) - number of orbitals
C            I = 0: total
C            I = 1: positive-energy orbitals
C            I = 2: negative-energy orbitals
C         OWNBAS - logical variable; if true, then fragments are defined 
C                  in their own basis
C         NUCFRAG(NFRAG) - number of symmmetry-independent nuclei of this fragment
C         IOPT - flag for the coefficient reader CMO
C         KVEC - (IFRP,NFRAG) - start address for list of selected orbitals of fermion ircop
C                IFRP for each fragment
C         CMO, EIG,IBEIG - buffer space for all coefficients, energies, irrep info
C         IOFF - offset in list of orbitals; allows to fill CSEL in several calls to SELFRAG
C     For each fragments, this information is first read into:
C        CMO - coefficients
C        EIG - eigenvalues
C        IBEIG - boson irrep information
C     Then, using the list of orbitals with adresses given in KVEC,
C     the selected information is extracted into the corresponding
C     arrays CSEL,ESEL and IBSEL
C
C     modifications for analysis of MCSCF coefficients
C     by S. Knecht - April 2010
C
C     Written by Trond Saue Sep 6 2006
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
#include "dcbham.h"
      LOGICAL TOBE,OWNBAS,LBIT
      CHARACTER FRAGFIL(NFRAG)*6
      DIMENSION CSEL(*),ESEL(*),IBSEL(*),NVECS(*),CMO(*),KVEC(2,NFRAG),
     &          EIG(*),IBEIG(*),NSTR(2,0:2,NFRAG),NUCFRAG(NFRAG),WORK(*)
      DIMENSION IOFF(2)
      REAL*8 CMCSCF_MAGIC
!
!     flag for MCSCF coefficents present on fragment/molecular
!     coefficients file
      CMCSCF_MAGIC = 137.0D0
      KRMC_FLG = 0
C
C.....loop over fragments
      INUC  = 1
      DO 10 IFRAG = 1,NFRAG
        NVECT = 0
        DO IFRP = 1,NFSYM
          NVECT = NVECT + NSTR(IFRP,0,IFRAG)
        ENDDO
        IF(NVECT.EQ.0) THEN
          WRITE(LUPRI,'(A,A6)') 
     &       '* No coefficients selected for fragment ',FRAGFIL(IFRAG)
          GOTO 12
        ENDIF
        INQUIRE(FILE=FRAGFIL(IFRAG),EXIST=TOBE)
        IF (.NOT.TOBE) THEN
           WRITE(LUPRI,'(A,I2,A,I2,A,A6,A)') 
     &     'Fragment coefficient file no. ',
     &     IFRAG,'/',NFRAG,' named ',FRAGFIL(IFRAG),' not present'
           CALL QUIT('SELFRAG: No fragment coefficients !')
        ENDIF
C
        CALL OPNFIL(LUCOEF,FRAGFIL(IFRAG),'OLD','SELFRAG')
        IF(LBIT(IOPT,1)) THEN
          WRITE(LUPRI,'(A,A6)') 
     &       '* SELFRAG: Coefficients from file: ',FRAGFIL(IFRAG)
        ENDIF
C.......Fragment orbitals in own basis
        IF(OWNBAS) THEN
          CALL SELOWN(IOPT,CMO,EIG,IBEIG,
     &                INUC,NUCFRAG(IFRAG),LUCOEF,FRAGFIL(IFRAG),
     &                CSEL,ESEL,IBSEL,NVECS,IOFF,KVEC(1,IFRAG),
     &                NSTR(1,0,IFRAG),WORK,KFREE,LFREE)
C.......Fragment orbitals in molecular basis
        ELSE
!
!         check for MCSCF coefficients on fragment/molecule file
          CALL REACMO(LUCOEF,FRAGFIL(IFRAG),CMO,EIG,IBEIG,
     &                TOTERG,-1)
!
!         found MCSCF coefficients on file CMOFIL?
          IF( TOTERG .eq. CMCSCF_MAGIC ) KRMC_FLG = 1
          CALL REACMO(LUCOEF,FRAGFIL(IFRAG),CMO,EIG,IBEIG,
     &                TOTERG,IOPT)
!
          ICOFF = 1
          IEOFF = 1
          DO IFRP = 1, NFSYM
            IF(NSTR(IFRP,0,IFRAG).GT.0) THEN
C             Select coefficients
              IF(LBIT(IOPT,2)) THEN
                ICSEL = ICOFF + NFBAS(IFRP,0)*IOFF(IFRP)
                CALL SELCFS (CMO(1+ICMOQ(IFRP)),IFRP,CSEL(ICSEL),
     &                       NVECS(IFRP),WORK(KVEC(IFRP,IFRAG)),
     &                       NSTR(IFRP,2,IFRAG),NSTR(IFRP,1,IFRAG),
     &                       NFBAS(IFRP,0),NORB(IFRP))
              ENDIF
              IESEL = IEOFF + IOFF(IFRP)
C             Select eigenvalues
              IF(LBIT(IOPT,3)) THEN
                CALL SELEIG(EIG(1+IORB(IFRP)),IFRP,ESEL(IESEL),
     &               WORK(KVEC(IFRP,IFRAG)),
     &               NSTR(IFRP,2,IFRAG),NSTR(IFRP,1,IFRAG))
              ENDIF
C       Select boson irrep information
              IF(LBIT(IOPT,4)) THEN
                CALL SELIBEIG(IBEIG,IORB(IFRP),IFRP,
     &                 IBSEL(IESEL),WORK(KVEC(IFRP,IFRAG)),
     &                 NSTR(IFRP,2,IFRAG),NSTR(IFRP,1,IFRAG))
              ENDIF
              ICOFF = ICOFF + NFBAS(IFRP,0)*NVECS(IFRP)*NZ
              IEOFF = IEOFF + NVECS(IFRP)
            ENDIF
          ENDDO
        ENDIF
        DO IFRP = 1,NFSYM
          IOFF(IFRP) = IOFF(IFRP) + NSTR(IFRP,0,IFRAG)
        ENDDO
        CLOSE(LUCOEF,STATUS='KEEP')
 12     CONTINUE
        INUC = INUC + NUCFRAG(IFRAG)
 10   CONTINUE
      
      RETURN
C
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphcar */
      SUBROUTINE SPHCAR(SPHM,IREP,IC,NCAR,NSPH,
     &                  DONRM,IPRINT,WORK,LWORK)
C***********************************************************************
C
C***  Generate spherical transformation matrix for component IC
C***  and irrep IREP
C
C     Input:
C       IREP - irrep
C       IC   - L(1),S(2) component
C       DONRM - normalize column by column
C       
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
C
#include "ccom.h"
      LOGICAL DONRM
      DIMENSION SPHM(*),WORK(LWORK)
#include "memint.h"
C.....LMAX is L + 1 where L is top angular momentum quantum number
      LMAX  = NHTYP
C.....NFUN gives the number of Cartesian functions up to and including LMAX
      NFUN  = (LMAX*(LMAX+1)*(LMAX+2))/6
C.....NCRED gives the sum of squares of Cartesian functions up to and 
C     including LMAX
      NCRED = ((LMAX+2)*(LMAX+1)*(3*LMAX**2+6*LMAX+1)*LMAX)/60
      CALL MEMGET('REAL',KCRED,NCRED ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRED,3*LMAX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFUN,3*NFUN,WORK,KFREE,LFREE)
      CALL SPHCA1(SPHM,IREP,IC,NCAR,NSPH,LMAX,NFUN,
     &            DONRM,WORK(KCRED),WORK(KIRED),WORK(KIFUN),IPRINT)
      CALL MEMREL('SPHCAR',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck sphma1 */
      SUBROUTINE SPHCA1(SPHM,IREP,IC,NCAR,NSPH,LMAX,
     &                  NFUN,DONRM,CRED,IRED,IFUN,IPRINT)
C***********************************************************************
C
C***  Generate spherical transformation matrix for component IC
C***  and irrep IREP
C
C     Input:
C       IREP - irrep
C       IC   - L(1),S(2) component
C       
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
C
#include "ccom.h"
#include "sphtrm.h"
#include "symmet.h"
#include "pgroup.h"
#include "shells.h"
#include "nuclei.h"
#include "dgroup.h"
      LOGICAL DONRM
      CHARACTER COMP(2)*1
      DIMENSION IFUN(NFUN,3),
     &          SPHM(NCAR,NCAR),CRED(*),IRED(3,LMAX)
      DATA COMP/'L','S'/
#include "ibtfun.h"
C
C     Initialize
C     ==========
C
C.....Large components
      IF(IC.EQ.1) THEN
        JA1 = 1
        JA2 = NLRGSH
C.....Small components
      ELSEIF(IC.EQ.2) THEN
        JA1 = NLRGSH+1
        JA2 = KMAX
      ENDIF
C.....
      N2MAT = NCAR*NCAR
      CALL DZERO(SPHM,N2MAT)
C
C     Initial print section
C     =====================
C
C     Loop over shells
C     ================
C
      ISOFF = ICOS(IREP+1,IC)
      NAORB = NSTRT(JA1)
      INUC = 0
      IROW = 0
      ICOL = 0
      DO 10 JA = JA1,JA2
        JNUC = NCENT(JA)
C
C       New center; check out allowed symmetries
C       ----------------------------------------
        IF(JNUC.NE.INUC) THEN
          INUC = JNUC
          ILL  = 0
          IOFF = 0
        ENDIF
        LL  = NHKT(JA)
C
C       New L-value; extract relevant transformation matrix
C       ---------------------------------------------------
        IF(LL.NE.ILL) THEN
          ILL  = LL
          CALL EXTSPH(LL,IREP,JNUC,CRED,LMAX,NFUN,IFUN,IRED,
     &                DONRM,IOFF,IPRINT)
        ENDIF
C
C       Expand transformation matrix
C       ============================
C
        KLM  = IRED(2,LL)
        IF(KLM.EQ.0) GOTO 10
        KXYZ = IRED(1,LL)
        KOFF = IRED(3,LL)+1
        DO J = 1,KLM
          CALL DCOPY(KXYZ,CRED(KOFF),1,
     &                SPHM(IROW+1,ICOL+J),1)   
          KOFF = KOFF + KXYZ
        ENDDO
        IROW = IROW + KXYZ
        ICOL = ICOL + KLM
   10 CONTINUE
      NSPH = ICOL
C
C     Final print section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A,A3,2X,A1)') 
     &    'SPHCAR: Transformation matrix for ',REP(IREP),COMP(IC)
        CALL PRSYMB(LUPRI,'-',40,0)
        CALL OUTPUT(SPHM,1,NCAR,1,NSPH,NCAR,NPSH,-1,LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Extsph */
      SUBROUTINE EXTSPH(LL,IREP,JNUC,CRED,LMAX,NFUN,IFUN,IRED,
     &                  DONRM,IOFF,IPRINT)
C***********************************************************************
C
C     Extract spherical transformation matrix for symmetry adapted
C     basis functions of angular momentum L=LL-1 of irrep IREP
C     and centered on symmetry-independent center JNUC
C
C       IFUN(*,1) - points to Cartesian component in full list
C       IFUN(*,2) - points to spherical component in full list
C       IFUN(*,3) - points to Cartesian comp. in symmetry-reduced list
C     The symmetry-reduced spherical transformation matrix is
C     stored in CRED. For each L value
C       IRED(1,*) - gives number of symmetry-reduced Cartesians KXYZ
C       IRED(2,*) - gives number of symmetry-reduced sphericals KLM
C       IRED(3,*) - gives offsets to each L value
C     For each L value CRED has the dimension CRED(KXYZ,KLM)
C
C     A cut'n paste job. T. Saue Jun 23 2007
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0)
C
#include "mxcent.h"
#include "maxaqn.h"
#include "maxorb.h"
C
#include "ccom.h"
#include "symmet.h"
#include "sphtrm.h"
#include "nuclei.h"
      LOGICAl DONRM
      DIMENSION IFUN(NFUN,3),IRED(3,LMAX),CRED(*)
C
#include "ibtfun.h"
      L    = LL - 1
      NLM  = 2*L + 1
      NXYZ = (LL*(LL+1))/2
      LOFF = NHKOFF(LL)
      MULA = ISTBNU(JNUC)
C
C     Cartesian functions
C
      KXYZ = 0
      DO ICMP = 1,NXYZ
        IVARB = IBTXOR(IREP,ISYMAO(LL,ICMP))
        IF(IBTAND(MULA,IVARB).EQ.0) THEN
          KXYZ              = KXYZ + 1
          IFUN(LOFF+KXYZ,1) = ICMP
          IFUN(LOFF+ICMP,3) = KXYZ
        ENDIF
        IRED(1,LL) = KXYZ
      ENDDO
C
C     Solid harmonic functions
C
      KLM = 0
      DO ICMP = 1,NLM
        M = MDEF(L,ICMP)
        IVARB = IBTXOR(IREP,IREPLM(L,M))
        IF(IBTAND(MULA,IVARB).EQ.0) THEN
          KLM               = KLM + 1
          IFUN(LOFF+KLM,2)  = ICMP
        ENDIF
        IRED(2,LL) = KLM
      ENDDO
C
C     Extract spherical transformation matrix
C
      IRED(3,LL) = IOFF
      IF(KLM.GT.0) THEN
        CALL EXTELM('T',CSP(ISPADR(LL)),NLM,NXYZ,
     &              CRED(IOFF+1),KXYZ,KXYZ,KLM,
     &              IFUN(LOFF+1,1),IFUN(LOFF+1,2))
        IF(DONRM) THEN
          JOFF = IOFF
          DO J = 1,KLM
            FAC = D1/DNRM2(KXYZ,CRED(JOFF+1),1)
            CALL DSCAL(KXYZ,FAC,CRED(JOFF+1),1)
            JOFF = JOFF + KXYZ
          ENDDO
        ENDIF
        IOFF = IOFF + KXYZ*KLM
      ENDIF
      IF(IPRINT.GE.4) THEN
        KOFF = IRED(3,LL)
        WRITE(LUPRI,'(A)') '* EXTSPH:'
        WRITE(LUPRI,'(A,I3)')
     &       ' * Center : ',JNUC,
     &       ' * L-value: ',L
        IF(KLM.GT.0) THEN
          WRITE(lupri,*) 'KXYZ,KLM=',KXYZ,KLM
          CALL OUTPUT(CRED(KOFF+1),1,KXYZ,1,KLM,KXYZ,KLM,-1,LUPRI)
        ELSE
          WRITE(LUPRI,'(A)') ' - no functions - '
        ENDIF
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkbcar */
      SUBROUTINE RKBCAR(SPHM,IREP,IC,NCAR,NSPH,
     &                  DONRM,IRKB,IPRINT,WORK,LWORK)
C***********************************************************************
C
C     Generate spherical transformation matrix for component IC
C     and irrep IREP
C
C     For a given L value there is, in addition to the standard
C     solid harmonic Gaussians modified functions of the form
C
C     X(L-2,M) = N * sum(i+j+k=L-2) [c_{ijk}^{L-2,M}
C                *(G_{i+2,j,k} + G_{i,j+2,k} + G_{i,j,k+2})-2(2L-1)G_{ijk}]
C
C     where the normalization constant is
C
C     N = 1/{(2L+1)(2L-1)}
C     
C     Input:
C       IREP - irrep
C       IC   - L(1),S(2) component
C       DONRM - normalize column by column
C       NCAR - number of Cartesian Gaussians of this symmetry
C     Output:
C       SPHM - transformation matrix
C       
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "maxaqn.h"
C
#include "ccom.h"
      LOGICAL DONRM
      DIMENSION SPHM(*),IRKB(*),WORK(LWORK)
#include "memint.h"
      LMAX  = NHTYP
      NFUN  = (LMAX*(LMAX+1)*(LMAX+2))/6
      NCRED = ((LMAX+2)*(LMAX+1)*(3*LMAX**2+6*LMAX+1)*LMAX)/60
      IF(LMAX.GE.3) THEN
        L2    = LMAX - 2
        NCRED = NCRED + ((L2+2)*(L2+1)*(3*L2**2+6*L2+1)*L2)/60
      ENDIF
      CALL MEMGET('REAL',KCRED,NCRED ,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIRED,3*LMAX,WORK,KFREE,LFREE)
      CALL MEMGET('INTE',KIFUN,3*NFUN,WORK,KFREE,LFREE)
      CALL RKBCA1(SPHM,IREP,IC,NCAR,NSPH,LMAX,NFUN,DONRM,
     &            IRKB,WORK(KCRED),WORK(KIRED),WORK(KIFUN),IPRINT)
      CALL MEMREL('RKBCAR',WORK,KWORK,KWORK,KFREE,LFREE)
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck rkbca1 */
      SUBROUTINE RKBCA1(SPHM,IREP,IC,NCAR,NSPH,LMAX,
     &                  NFUN,DONRM,IRKB,CRED,IRED,IFUN,IPRINT)
C***********************************************************************
C
C***  Generate modified spherical transformation matrix for component IC
C***  and irrep IREP
C
C       IFUN(*,1) - points to Cartesian component in full list
C       IFUN(*,2) - points to spherical component in full list
C       IFUN(*,3) - points to Cartesian comp. in symmetry-reduced list
C     The symmetry-reduced spherical transformation matrix is
C     stored in CRED. For each L value
C       IRED(1,*) - gives number of symmetry-reduced Cartesians KXYZ
C       IRED(2,*) - gives number of symmetry-reduced sphericals KLM
C       IRED(3,*) - gives offsets to each L value
C     Input:
C       IREP - irrep
C       IC   - L(1),S(2) component
C       
C
CMI:  IPRINT is IPRHAM (from GMOTRA)
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0 = 0.0D0,D1 = 1.0D0,D2=2.0D0)
C
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "maxmom.h"
#include "aovec.h"
C
#include "ccom.h"
#include "xyzpow.h"
#include "sphtrm.h"
#include "primit.h"
#include "symmet.h"
#include "pgroup.h"
#include "shells.h"
#include "nuclei.h"
#include "dgroup.h"
      LOGICAL DONRM,LRKB
      CHARACTER COMP(2)*1
      DIMENSION IFUN(NFUN,3),SPHM(NCAR,NCAR),CRED(*),IRED(3,LMAX),
     &          IRKB(NPSHEL)
      DATA COMP/'L','S'/
#include "ibtfun.h"
C
C
C     Initialize
C     ==========
C
      IF(IC.EQ.1) THEN
        JA1 = 1
        JA2 = NLRGSH
      ELSEIF(IC.EQ.2) THEN
        JA1 = NLRGSH+1
        JA2 = KMAX
      ENDIF
      N2MAT = NCAR*NCAR
      CALL DZERO(SPHM,N2MAT)
C
C     Initial print section
C     =====================
C
      IF(IPRINT.GE.4) THEN
        WRITE(LUPRI,'(/,2X,A,A3,2X,A1)') 
     &  '*** RKBCA1: Transformation matrix for ',REP(IREP),COMP(IC)
        WRITE(LUPRI,*) 'JA1=',JA1,' JA2=',JA2
        CALL FLSHFO(LUPRI) 
      ENDIF
C
C     Loop over shells
C     ================
C     (Shell number will be different from primitive number
C      if large component is contracted, thus the double
C      loop structure)
C
C     Offset in SO-basis:
      ISOFF = ICOS(IREP+1,IC)
C     Offset in AO-basis:
      NAORB = NSTRT(JA1)
C
      INUC = 0
      ICOL = 0
      JP1_TEST = 0
C.....loop over shells
      DO 10 JA = JA1,JA2
        JP1  = JSTRT(JA) + 1
C
C       Skip primitives already processed in another shell.
C
        IF ( JP1_TEST .EQ. JP1 ) GO TO 10
        JP1_TEST = JP1
        JP2  = JSTRT(JA) + NUCO(JA)
C
        JNUC = NCENT(JA)
C
C       New center; check out allowed symmetries
C       ----------------------------------------
        IF(JNUC.NE.INUC) THEN
          INUC = JNUC
          ILL  = 0
          IOFF = 0
        ENDIF
        LL  = NHKT(JA)

        IF(IPRINT.GE.8) THEN
          write(LUPRI,*) 'JA, NHKT(JA)=',JA, NHKT(JA)
          write(LUPRI,*) 'LL=',LL
          CALL FLSHFO(LUPRI) 
        ENDIF

C   
C       New L-value; extract relevant transformation matrix
C       ---------------------------------------------------
        IF(LL.NE.ILL) THEN
          ILL  = LL
          JRKB = 0
          LOFF = NHKOFF(LL)
          CALL EXTSPH(LL,IREP,JNUC,CRED,LMAX,NFUN,IFUN,IRED,
     &                 DONRM,IOFF,IPRINT)
C.........these factors assure the correct normalized linear
C.........combination of L and (L-2) for the modified solid harmonic
          L = LL - 1
          IF(L.GT.0) THEN
            FACL=D1/SQRT((D2*dble(L)+D1)*(D2*dble(L)-D1))
          ELSE
            FACL   = D0
          ENDIF
          FACLM2 = -D2*(D2*dble(L)-D1)*FACL
        ENDIF
C
C       Expand transformation matrix
C       ============================
C
        KLM = IRED(2,LL)
        IF(KLM.EQ.0) GOTO 10
        KXYZ = IRED(1,LL)
C
C.......Loop over primitives
        DO 15 JP = JP1, JP2
C
C         A function can contribute to kinetic balance in two ways :
C         a. It balances a (L-1) large component function
C         b. It balances a (L+1) large component function
C
C         In case b it needs to be combined with a (L+2) function,
C         and we take this combination when we enter this (L+2) shell.
C         This means we can disgard the functions for the moment,
C         provided that case a. does not apply as well (as may happen
C         with family-type basis sets).
C
C         Look for (L-2) small function. If it does not exist, 
C         IRKB returns the pointer to (L-1) large component function
C         (case b). If this is zero it means that we have the
C         first function for case b and we may skip it.
C
          JB = IRKB(JP)
          IF (JB.EQ.0) GOTO 15
C
          IROW = NSTRT(JA + JP - JP1) + IFUN(LOFF+1,1)
C                       ^- to get correct shell number
C
C         JA is the pointer to the start of this series of uncontracted
C         shells, increase it by the distance from the start
C
          IROW = IPTSYM(IROW,IREP) - ISOFF
C
C         Fill in standard spherical transformation matrix
C
          KOFF = IRED(3,LL)+1
          DO J = 1,KLM
            CALL DCOPY(KXYZ,CRED(KOFF),1,SPHM(IROW,ICOL+J),1)   
            KOFF = KOFF + KXYZ
          ENDDO
          ICOL = ICOL + KLM
C
C         Possibly add RKB part
C
          LRKB = LL.GE.3
          IF(LRKB) THEN
            LLM2   = LL-2
            KLMM2  = IRED(2,LLM2)
            LRKB   = KLMM2.GT.0
          ENDIF
          IF(LRKB) THEN
C
C           Possibly generate RKB matrix
C
            IF(JRKB.EQ.0) THEN
              JRKB = 1
              NXYZM2 = (LLM2*(LLM2+1))/2
              NLMM2  = 2*LLM2-1
              KXYZM2 = IRED(1,LLM2)                            
              KRKB   = KXYZ + KXYZM2
              IM2    = IOFF + KXYZ
              LM2F   = NHKOFF(LLM2)
              CALL DZERO(CRED(IOFF+1),KRKB*KLMM2)
C
C             First extract spherical transformation matrix 
C             for L-2 partner
C
              CALL EXTELM('T',CSP(ISPADR(LLM2)),NLMM2,NXYZM2,
     &                   CRED(IM2+1),KRKB,KXYZM2,KLMM2,
     &                   IFUN(LM2F+1,1),IFUN(LM2F+1,2))
C
C              Then fill in corresponding number for r**2
C
C              A triplet (i,j,k) can be generated by looping over the
C              lower triangle of a square matrix of dimension L + 1
C              From an element (a,b) of the matrix the values of the triplet
C              are given by
C
C              i = L + 1 - a; j = a - b; k = b - 1
C
C              The inverse relations are:
C
C              b = k + 1; a = j + k + 1
C                 
C              Consider a triplet (i,j,k) corresponding to a matrix
C              element (a,b). The above relations show that
C                (i+2,j,k) --> (a,b)
C                (i,j+2,k) --> (a+2,b)
C                (i,j,k+2) --> (a+2,b+2)
C
              DO K = 1,KXYZM2
                IJ = IFUN(LM2F+K,1)
                IB = NVAL(IJ) + 1
                IA = MVAL(IJ) + IB
C**             x**2: (ia  ,ib  )
                IJX = IJ
                KOFF = IOFF + IFUN(LOFF+IJX,3)
                CALL DAXPY(KLMM2,FACL,CRED(IM2+K),KRKB,
     &                    CRED(KOFF),KRKB)
C**             y**2: (ia+2,ib  )
                IJY = IJ + IA + IA + 1
                KOFF = IOFF + IFUN(LOFF+IJY,3)
                CALL DAXPY(KLMM2,FACL,CRED(IM2+K),KRKB,
     &                     CRED(KOFF),KRKB)
C**             z**2: (ia+2,ib+2)
                IJZ = IJY + 2
                KOFF = IOFF + IFUN(LOFF+IJZ,3)
                CALL DAXPY(KLMM2,FACL,CRED(IM2+K),KRKB,
     &                    CRED(KOFF),KRKB)
              ENDDO
C
C             Scale matrix for L-2 partner
C
              KOFF = IM2+1
              DO K = 1,KLMM2
                CALL DSCAL(KXYZM2,FACLM2,CRED(KOFF),1)
                KOFF = KOFF + KRKB
              ENDDO
C
C             Normalize...
C
              IF(DONRM) THEN
                KOFF = IOFF+1
                DO J = 1,KLMM2
                  FAC = D1/DNRM2(KRKB,CRED(KOFF),1)
                  CALL DSCAL(KRKB,FAC,CRED(KOFF),1)
                  KOFF = KOFF + KRKB
                ENDDO
              ENDIF
              IF(IPRINT.GE.4) THEN
                WRITE(LUPRI,'(A)') '* RKB part:'
                CALL OUTPUT(CRED(IOFF+1),1,KRKB,1,KLMM2,
     &                      KRKB,KLMM2,-1,LUPRI)
              ENDIF
              IOFF = IOFF + KRKB*KLMM2
            ENDIF
C
C           Copy extra functions for L
C            
            IRM2 = NSTRT(JB) + IFUN(LM2F+1,1)
            IRM2 = IPTSYM(IRM2,IREP) - ISOFF
            KOFF = IRED(3,LL)+KXYZ*KLM+1
            DO J = 1,KLMM2
              CALL DCOPY(KXYZ,CRED(KOFF),1,
     &                   SPHM(IROW,ICOL+J),1)   
              KOFF = KOFF + KXYZ
              CALL DCOPY(KXYZM2,CRED(KOFF),1,
     &                   SPHM(IRM2,ICOL+J),1)   
              KOFF = KOFF + KXYZM2
            ENDDO
            IROW = IROW + KXYZ
            ICOL = ICOL + KLMM2
          ENDIF
 15     CONTINUE
 10   CONTINUE
      NSPH = ICOL
C
C     Final print section
C
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A,A3,2X,A1)') 
     &    'RKBCAR: Transformation matrix for ',REP(IREP),COMP(IC)
        WRITE(LUPRI,*) ' NCAR and NSPH:',NCAR,NSPH
        CALL PRSYMB(LUPRI,'-',40,0)
        CALL OUTPUT(SPHM,1,NCAR,1,NSPH,NCAR,NSPH,-1,LUPRI)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck acmoin */
      SUBROUTINE ACMOIN()
C*****************************************************************************
C
C     Read coefficients in C1 symmetry from unformatted DFACMO,
C     adapt to current symmetry and write to unformatted DFCOEF
C
C     Written by T. Saue Mar 3 2008
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
#include "mxcent.h"
#include "maxorb.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcbdhf.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dcblab.h"
#include "dgroup.h"
#include "nuclei.h"
      LOGICAL TOBE,FNDLAB
      CHARACTER TEXT*74
      DIMENSION IDIM(3),NBF(2)
      real(8), allocatable :: CMO1(:)
      real(8), allocatable :: CMO2(:)
      real(8), allocatable :: EIG1(:)
      real(8), allocatable :: EIG2(:)
      integer, allocatable :: IBEIG1(:)
      integer, allocatable :: IBEIG2(:)
      integer, allocatable :: IBUF(:)
C
      IPRINT = IPRGEN
C
C     First check if file DFACMO is there
C
      INQUIRE(FILE='DFACMO',EXIST=TOBE)
      IF(.NOT.TOBE) THEN
          WRITE(LUPRI,'(A)') 
     &      'ACMOIN: Unformatted file DFACMO not found !'
          CALL QUIT('ACMOIN: Unformatted file DFACMO not found !')
      ENDIF
C     Allocate arrays
      allocate(CMO1(N2BBASX*4))
      allocate(CMO2(N2BBASX*4))
      allocate(EIG1(NORBT))
      allocate(EIG2(NORBT))
      allocate(IBEIG1(NORBT))
      allocate(IBEIG2(NORBT))
      allocate(IBUF(NTBAS(0)*2))
C
C     Get coefficients from DFACMO --> stored in CMO2
C     ===============================================
C
      CALL OPNFIL(LUCOEF,'DFACMO','OLD','ACMOIN')
      IF(FNDLAB('INFO    ',LUCOEF)) THEN
         READ(LUCOEF) TEXT,NSYM,NZBUF,(IDIM(I),I = 1,3),TOTERG
      ELSE
         CALL QUIT('ACMOIN: Error reading info')
      ENDIF
      NDIM = NTBAS(0)*NORBT*4
      CALL REACMO_coefficients(LUCOEF,CMO2,NDIM,1,4,
     &           NORBT,NTBAS(0),NPSHT,NESHT,IDIM)         
      CALL REACMO_eigenvalues(LUCOEF,EIG2,NORBT,1,4,
     &           NORBT,NPSHT,NESHT,IDIM)         
C      CALL REACMO_supersymmetry(LUCOEF,'SUPERSYM',IBEIG2,NORBT,
C     &     1,4,NORBT,NPSHT,NESHT,IDIM)
      CALL IZERO(IBEIG2,NORBT)
      CLOSE(LUCOEF,STATUS='KEEP')

      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A)') 
     &    '* ACMOIN : Coefficients read from DFACMO'
          CALL PRQMAT(CMO2,NTBAS(0),NORBT,
     &                NTBAS(0),NORBT,4,IQDEF,LUPRI)
      ENDIF
C     Find correct ordering of nuclei
      CALL NUCORD(IBUF)
C
C
C     Get correct row ordering (NBAS) --> stored in CMO1     
C     ==================================================
C
      CALL ACMOU2(CMO1,NORBT,CMO2,NORBT,NTBAS(0),IBUF,4)
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A)') 
     &    '* ACMOIN : row-ordered coefficients'
          CALL PRQMAT(CMO1,NTBAS(0),NORBT,
     &                NTBAS(0),NORBT,4,IQDEF,LUPRI)
      ENDIF
C
C     Transform to symmetry-adapted (SO) basis (unsorted) -> CMO2
C
      NDIM = NORBT*4
      CALL VTAOSO(CMO1,CMO2,NTBAS(0),NDIM,0)
      IF(IPRINT.GE.3) THEN
        WRITE(LUPRI,'(A)') 
     &    '* ACMOIN : raw SO coefficients'
          CALL PRQMAT(CMO2,NTBAS(0),NORBT,
     &                NTBAS(0),NORBT,4,IQDEF,LUPRI)
      ENDIF
C
C     Analyze distribution of symmetry on quaternion units
C     and prepare for compression
C
        CALL QSYMANA(CMO2,NTBAS(0),NORBT,NBF,IBUF)
        DO IFRP = 1,NFSYM
          IF(NORB(IFRP).NE.NBF(IFRP)) THEN
            WRITE(LUPRI, *)
     &        '*** ERROR: ACMOIN Fermion ircop:',IFRP,
     &        '     Expected : ',NORB(IFRP),
     &        '     Found    : ',NBF(IFRP)
            CALL QUIT('ACMOIN: Wrong number of orbitals !')
          ENDIF
        ENDDO
C
C     Compact orbitals (and orbital energies) to sorted basis -> CMO1
C
      CALL MKQCMO(CMO2,EIG2,IBEIG2,CMO1,EIG1,IBEIG1,NBF,IBUF)
      IF(IPRINT.GE.3) THEN
        DO IFRP = 1,NFSYM
          WRITE(LUPRI,'(A,A,I2)') 
     &      '* ACMOIN : final SO coefficients.',
     &      ' Fermion ircop ',IFRP
          CALL PRQMAT(CMO1,NTBAS(0),NORBT,
     &                NTBAS(0),NORBT,NZ,IPQTOQ(1,0),LUPRI)
        ENDDO
      ENDIF
C
      CALL OPNFIL(LUCOEF,'DFCOEF','UNKNOWN','ACMOIN')
C     tsaue:
C     A possible source of error is that IBEIG1 may contain information
C     about linear or boson symmetry, and the reduction of the variational
C     space will not work in the default case. This needs some thinking; for 
C     now I simply set IBEIG1 to zero.
      CALL ICOPY (NORBT,0,0,IBEIG1,1) 
      CALL WRICMO(LUCOEF,CMO1,EIG1,IBEIG1,TOTERG,ATOMIC)
      CLOSE(LUCOEF,STATUS='KEEP')
C
C     De-allocate arrays
      deallocate(CMO1)
      deallocate(CMO2)
      deallocate(EIG1)
      deallocate(EIG2)
      deallocate(IBEIG1)
      deallocate(IBEIG2)
      deallocate(IBUF)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck qsymana */
      SUBROUTINE QSYMANA(CMO,NB,NO,NBF,IBUF)
C***********************************************************************
C
C     Input: Non-compressed quaternion MO-coefficients
C            in terms of symmetry-adapted basis functions
c            in unsorted basis (HERMIT order).
C     This routine will prepare for compression by checking if the 
C     coefficients have the right (quaternion) phase, that is that 
C     the symmetry-adapted functions are in the expected positions.
C
C     In the DIRAC the real and imaginary parts of a spinor have
C     fixed symmetries, that is
C
C     La  ( 0 , R_z)
C     Lb  (R_y, R_x)
C     Sa  (xyz, z  )
C     Sb  ( y , x  )
C
C     where R_q is rotation about axis q.
C
C     If the system has inversion symmetry such that there is a second
C     fermion ircop, the symmetry distribution in the second ircop is
C
C     La  (xyz, z  )
C     Lb  ( y , x  )
C     Sa  ( 0 , R_z)
C     Sb  (R_y, R_x)
C
C     This information is contained in the array JSPINR(IM,IC,IFRP)
C     where IM = 1,4 is the part, IC=1,2 refers to large or small component
C     and IFRP refers to fermion ircop.
C            
C     A quaternion number is given by 
C       q = a + bi + cj +dk
C     where (i,j,k) are the quaternion units. In general the quaternion
C     units i,j and k come in positions 2, 3, and 4, respectively.
C     There are two exceptions for complex group, having one totally symmetric
C     rotation (see routine DBLGRP):
C      - for a totally symmetric x-rotation the mapping is (x,y,z) --> (i,k,j)
C      - for a totally symmetric y-rotation the mapping is (x,y,z) --> (j,i,k)
C     contrary to the standard mapping
C      (x,y,z) --> (k,j,i)
C
C     Written by T. Saue Mar 4 2008
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "mxcent.h"
#include "maxaqn.h"
      PARAMETER(D0=0.0D0,DTOL=1.0D-1)
C
#include "dcbbas.h"
#include "dgroup.h"
#include "pgroup.h"
#include "symmet.h"
#include "dcbham.h"
      DIMENSION CMO(NB,NO,4),IBUF(NO,2)
      DIMENSION TMP(2),NBF(2),III(2),WB(4,8,2),IB(64,4)
#include "dcbibt.h"
C
C.....III is the totally symmetric irrep for the large components
      III(1) = 0
C.....III is the symmety of XYZ for the small components
      III(2) = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
C     Initialize
      II = 0
      DO IC = 1,MC
        DO IREP = 0,7
          DO IQ = 1,4
            II = II + 1
            IB(II,2) = IQ
            IB(II,3) = IREP
            IB(II,4) = IC
          ENDDO
        ENDDO
      ENDDO
      NBF(1) = 0
      NBF(2) = 0
      DO JORB = 1,NO
        CALL DZERO(WB,64)
C
C       Accumulate gerade and ungerade norms
C
        TMP(1)=D0
        TMP(2)=D0
        DO IC = 1,MC
          DO IREP = 0,MAXREP
            ISYM = IREP + 1
            IFRP = JBTOF(IREP,IC)
            IBB  = ICOS (ISYM,IC)
            NBB  = NBBAS(IREP,IC)
            IF(NBB.GT.0) THEN
C....Carefully note that the assignment of the quaternion units (i,j,k) with
C....respect to axes (x,y,z) may change; 
C....This is the case of complex groups with totally symmetric x- or y-rotation
C....and is the reason why the array IQMAP is used. It refers to the default
C....assigment used for the C1 coefficients.
              DO IQQ = 1,4
                IQ=IQMAP(IQQ)
                WB(IQQ,ISYM,IC) = DDOT(NBB,CMO(IBB+1,JORB,IQ),1,
     &                   CMO(IBB+1,JORB,IQ),1)
                TMP(IFRP) = TMP(IFRP) + WB(IQQ,ISYM,IC)
              ENDDO
C
C             Order contributions
C
              CALL INDEXX(4,WB(1,ISYM,IC),IB)
C       
C             Find the number of non-zero contributions;
C             it has to be less or equal to NZ,
C             otherwise the compression of coefficients will fail.
C
              MZ = 0
              DO IQ = 4,1,-1
                IF(WB(IB(IQ,1),ISYM,IC).GT.D0) MZ = MZ + 1
              ENDDO
              IF(MZ.GT.NZ) THEN
                WRITE(LUPRI,'(A,I5)')
     &          '*** ERROR (QSYMANA): MO-coefficient no. ',JORB
                IF(IC.EQ.1) THEN
                  WRITE(LUPRI,'(3X,A,A3,A,I3,A)')
     &          'Large component basis functions of symmetry ',
     &           REP(IREP),' has ',MZ,' non-zero components'
                ELSE 
                  WRITE(LUPRI,'(3X,A,A3,A,I3,A)')
     &          'Large component basis functions of symmetry ',
     &           REP(IREP),' has ',MZ,' non-zero components'
                ENDIF
                CALL QUIT('QSYMANA. Compression will fail !')
              ENDIF
            ENDIF
          ENDDO ! End loop over boson irreps
        ENDDO ! End loop over large and small components
C
C       Decide parity of MO
C
        IF(TMP(1).GT.TMP(2)) THEN
          IBUF(JORB,1) = 1
          NBF(1) = NBF(1) + 1
          IF(TMP(2).GT.DTOL) THEN
            WRITE(LUPRI,'(A,I5,3X,E7.2)')
     &        '*** WARNING from QSYMANA: Contamination of orbital:',
     &        JORB,TMP(2)
          ENDIF
        ELSE
          IBUF(JORB,1) = 2
          NBF(2) = NBF(2) + 1
          IF(TMP(1).GT.DTOL) THEN
            WRITE(LUPRI,'(A,I5,3X,E7.2)')
     &        '*** WARNING from QSYMANA: Contamination of orbital:',
     &        JORB,TMP(1)
          ENDIF
        ENDIF
C
C       Check if a quaternion phase transformation (IQ) is necessary
C
        CALL INDEXX(64,WB,IB)
C
C       For now only consider largest contribution;
C       a consistency check can be included later
C       by scanning through more contributions
C
        II   = IB(64,1)
        IREP = IB(II,3)
        IC   = IB(II,4)
C.......IQ is the position where the contribution from (IREP,IC) was found
        IQ   = IB(II,2)
C.......KQ is the position where this contribution would be expected
C.......(note that for complex/quaternion groups there are two/four possibilities, though)
        KQ   = JQBAS(IREP,IC)
C.......IBUF(JORB,2) indicates the quaternion phase transformation necessary to bring the
C.......contribution to the desired position
        IBUF(JORB,2) = IQMULT(IQ,KQ,1)
      ENDDO ! End loop over coefficients
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck mkqcmo */
      SUBROUTINE MKQCMO(CBU,EIGU,IBEIGU,CBS,EIGS,IBEIGS,NBF,IBUF)
C***********************************************************************
C
C     Starting from quaternion MO-coefficients in unsorted basis
C     make compressed coefficients in sorted basis
C
C     Written by T.Saue Mar 4 2008
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
C
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
C
#include "dcblab.h"
#include "dcbham.h"
#include "dgroup.h"
#include "symmet.h"
#include "dcbbas.h"
#include "dcborb.h"
C
      DIMENSION CBU(*),EIGU(*),IBEIGU(*),
     &          CBS(NCMOTQ),EIGS(*),IBEIGS(*),IBUF(NORBT,2),NBF(2)
C
      CALL DZERO(CBS,NCMOTQ)
      NBF(1) = 0
      NBF(2) = 0
      DO JU = 1,NORBT
        IFRP = IBUF(JU,1)
        JQ   = IBUF(JU,2)
        NBS = NFBAS(IFRP,0)
        NOS = NORB(IFRP)
        NBF(IFRP) = NBF(IFRP)+1
        JS = IORB(IFRP) + NBF(IFRP)
        EIGS(JS) = EIGU(JU)
        IBEIGS(JS) = IBEIGU(JU)
        IOFF = ICMOQ(IFRP)+1
        CALL MKQCMO1(CBU,NTBAS(0),NORBT,JU,IFRP,JQ,
     &               CBS(IOFF),NBS,NOS,NBF(IFRP))
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck mkqcmo1 */
      SUBROUTINE MKQCMO1(CBU,NBU,NOU,IVU,IFRP,JQ,CBS,NBS,NOS,IVS)
C***********************************************************************
C
C     Starting from quaternion MO-coefficient CBU(IVU) in unsorted basis
C     make compressed coefficient CBS(IVS) in sorted basis.
C
C     JQ is a quaternion unit used for phase transformation if needed.
C
C     Written by T.Saue Mar 4 2008
C
C***********************************************************************
      use quaternion_algebra
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "mxcent.h"
      PARAMETER(DM1=-1.0D0)
C
#include "symmet.h"
#include "dgroup.h"
#include "dcbham.h"
#include "dcbbas.h"
      DIMENSION CBU(NBU,NOU,4),CBS(NBS,NOS,NZ),III(2)
#include "dcbibt.h"
C
C.....III is the totally symmetric irrep for the large components
      III(1) = 0
C.....III is the symmety of XYZ for the small components
      III(2) = IBTXOR(ISYMAX(1,1),ISYMAX(1,2))
      NBRP = 4/NZ
C.....loop over components
      DO IC = 1,MC
C.......loop over symmetries associated with fermion ircop and component
        DO ISYM = 1,NBRP
          IP    = MOD(IFRP+IC,2) + 1
          IBSYM = JFSYM(ISYM,IP)
          IBRP  = IBSYM - 1
          IS    = IBBAS(IBRP,IC) - IBAS(IFRP) + 1
          IU    = ICOS(IBSYM,IC) + 1 
          IBRQ  = IBTXOR(IBRP,III(IC))
          NBB   = NBBAS(IBRP,IC)
C.........Note the use of IQMAP referring to the default asignment of
C.........quaternion units vectors to coordinate axes (x,y,z)
          DO IZ = 1,NZ
            IQ = IPQTOQ(IZ,IBRQ)
            KQ = IQMAP(IQMULT(IQ,JQ,1))
            CALL DCOPY(NBB,CBU(IU,IVU,KQ),1,CBS(IS,IVS,IZ),1)
            IFAC = IQPHASE(IQ,JQ,1)
            IF(IFAC.EQ.-1) CALL DSCAL(NBB,DM1,CBS(IS,IVS,IZ),1)
          ENDDO
        ENDDO
      ENDDO
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rkbimp */
      SUBROUTINE RKBIMP(IPRINT)
C***********************************************************************
C
C     Import RKB coefficients and add complementary UKB space
C
C     Written by T. Saue Feb 22 2010
C
C***********************************************************************

      use memory_allocator
#include "implicit.h"
#include "priunit.h"
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dcborb.h"
#include "dcbbas.h"
#include "dgroup.h"
      LOGICAL TOBE,FNDLAB
      CHARACTER TEXT*74
      DIMENSION IDIM(3,2)

      real(8), allocatable :: buf(:)
      real(8), allocatable :: cmo(:)
      real(8), allocatable :: eig(:)
      integer, allocatable :: ibeig(:)
      real(8), allocatable :: tmat(:)
      real(8), allocatable :: WORK(:)

      CALL QENTER('RKBIMP')      
      call legacy_lwork_get(LWORK)
#include "memint.h"
      call alloc(WORK,LWORK,id='WORK in RKBIMP')

C     Check for the presence of coefficient file
      INQUIRE(FILE='DFCOEF',EXIST=TOBE)
      IF(TOBE) THEN
        CALL OPNFIL(LUCOEF,'DFCOEF','OLD','RKBIMP')
        IF(FNDLAB('INFO    ',LUCOEF)) THEN         
C         IDIM(1,IFRP) - number of positronic solutions
C         IDIM(2,IFRP) - number of electronic solutions
C         IDIM(3,IFRP) - number of AO-basis functions
          READ (LUCOEF,END=10,ERR=20) TEXT,NSYM,NZBUF,
     &         ((IDIM(I,J),I = 1,3),J=1,NSYM),TOTERG
      ELSE
         CALL QUIT('RKBIMP: Error reading info')
      ENDIF
        
C.......Check that the number of negative energy solutions is smaller than UKB
C.......and calculate space needed to store UKB complement in orthonormal basis
        NBUF = 0
        DO IFRP = 1,NFSYM
          NUKBP = NPSH(IFRP)-IDIM(1,IFRP)
          IF(NUKBP.GT.0) THEN
            TOBE = .TRUE.
            NBUF = NBUF + NORB(IFRP)*NORB(IFRP)*NZ
          ENDIF
        ENDDO
        IF(TOBE) THEN
C.........Read coefficients, eigenvalues and boson irreps
          IOPT=14
          REWIND LUCOEF
          allocate(buf(nbuf))    ! buffer for storing UKB complement in orthonormal basis
          allocate(cmo(ncmotq))
          allocate(eig(norbt))
          allocate(ibeig(norbt))
          allocate(tmat(n2tmt))
C
          CALL REACMO(LUCOEF,'DFCOEF',cmo,eig,ibeig,TOTERG,IOPT)
C.........Zero IBEIG if no subblocks; this allows to treat linear systems
C         that are presently not handled in URKBAL
          IF(SUB_BL) THEN
            CALL QUIT('RKBIMP: SUB_BL must be programmed.')
          ELSE
            CALL ICOPY(NORBT,0,0,IBEIG,1)
          ENDIF
C.........At this point the first columns of the negative energy solutions
C         corresponding to the UKB complement has been zeroed out by REACMO
          CALL OPNFIL(LUTMAT,'AOMOMAT','OLD','RKBIMP')
C.........skip record with AO-to-MO transformation matrix
          READ(LUTMAT)
C.........Read MO-to-AO-transformation matrix
          CALL READT(LUTMAT,N2TMT,tmat)
          IOFF = 0
          DO IFRP = 1,NFSYM
            NUKBP=NPSH(IFRP)-IDIM(1,IFRP)
            IF(NUKBP.GT.0) THEN
              NRKB=IDIM(1,IFRP)+IDIM(2,IFRP)
              CALL RKBIMP_1(IFRP,NFBAS(IFRP,0),NORB(IFRP),NRKB,NUKBP,
     &             CMO(1+ICMOQ(IFRP)),BUF(1+IOFF),
     &             TMAT(1+I2TMT(IFRP)),NZT,WORK,KFREE,LFREE)
              IOFF = IOFF + NORB(IFRP)*NORB(IFRP)*NZ
            ENDIF
          ENDDO
C.........Backtransform UKB complement
          REWIND LUTMAT
          CALL READT(LUTMAT,N2TMT,TMAT)
          IOFF = 0
          DO IFRP = 1,NFSYM
            NUKBP=NPSH(IFRP)-IDIM(1,IFRP)
            IF(NUKBP.GT.0) THEN
              CALL RKBIMP_2(IFRP,NFBAS(IFRP,0),NORB(IFRP),NUKBP,
     &                      CMO(1+ICMOQ(IFRP)),EIG(1+IORB(IFRP)),
     &                      BUF(1+IOFF),TMAT(1+I2TMT(IFRP)),
     &                      NZT,IPRINT)
              IOFF = IOFF + NORB(IFRP)*NORB(IFRP)*NZ
            ENDIF
          ENDDO
          CALL WRICMO(LUCOEF,CMO,EIG,IBEIG,DHFERG,ATOMIC)

          deallocate(buf)
          deallocate(cmo)
          deallocate(eig)
          deallocate(ibeig)
          deallocate(tmat)

          CLOSE(LUTMAT,STATUS='KEEP')
          CLOSE(LUCOEF,STATUS='KEEP')
        ELSE
          CLOSE(LUCOEF,STATUS='KEEP')          
          WRITE(LUPRI,'(A)')
     &      '*** WARNING! RKBIMP: No need for UKB extension. Ignored.'
          RETURN
        ENDIF
      ELSE
        CALL QUIT('RKBIMP: Coefficient file not found !')
      ENDIF
C
      call dealloc(work)
      CALL QEXIT('RKBIMP')
      RETURN
 10   CONTINUE
      CALL QUIT('RKBIMP: END OF FILE reading TEXT')
 20   CONTINUE
      CALL QUIT('RKBIMP: ERROR reading TEXT')
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rkbimp_1 */
      SUBROUTINE RKBIMP_1(IFRP,NBAS,NORB,NRKB,NUKBP,CMO,CBUF,TINV,NZT,
     &                    WORK,KFREE,LFREE)
C***********************************************************************
C
C     Extend RKB coefficients of fermion ircop IFRP by UKB complement
C
C     Written by T. Saue Feb 23 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0)
C
#include "dcbham.h"
#include "dgroup.h"
      DIMENSION CMO(NBAS,NORB,NZ),TINV(NBAS,NORB,NZ),CBUF(NORB,NORB,NZ),
     &          WORK(*)
      real(8), allocatable :: cf(:)
C.....Transform RKB coefficients to orthonormal basis
      ICOL = 1 + NUKBP
      NDIM = NORB*NRKB*NZ
      allocate(cf(ndim))
      CALL QGEMM(NORB,NRKB,NBAS,D1,
     &            'H','N',IPQTOQ(1,0),TINV,NBAS,NORB,NZT,
     &            'N','N',IPQTOQ(1,0),CMO(1,ICOL,1),NBAS,NORB,NZ,
     &     D0,IPQTOQ(1,0),CF,NORB,NRKB,NZ)
C.....Generate complement (it fills the first NUKBP columns of CBUF
      CALL PRJORT(CBUF,CF,NORB,NRKB,NZ,WORK,KFREE,LFREE)
      deallocate (cf)
c     
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck rkbimp_2 */
      SUBROUTINE RKBIMP_2(IFRP,NBAS,NORB,NUKBP,CMO,EIG,
     &                    CBUF,TMAT,NZT,IPRINT)
C***********************************************************************
C
C     Extend RKB coefficients of fermion ircop IFRP by UKB complement
C
C     Written by T. Saue Feb 23 2010
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
      PARAMETER(D0=0.0D0,D1=1.0D0,D2=2.0D0)
C
#include "dcbgen.h"
#include "dcbham.h"
#include "dgroup.h"
      DIMENSION CMO(NBAS,NORB,NZ),EIG(NORB),
     &          TMAT(NBAS,NORB,NZ),CBUF(NORB,NORB,NZ)
      IF(IPRINT.GE.5) THEN
          WRITE(6,*) '* RKBIMP: CBUF on input',I
          CALL PRQMAT(CBUF(1,1,1),NORB,NORB,NORB,NORB,
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
C.....Approximate eigenvalue of UKB complement
      FAC=-D2*CVAL*CVAL
      CALL DCOPY(NUKBP,FAC,0,EIG,1)
C.....Backtransform UKB complement      
      CALL QGEMM(NBAS,NUKBP,NORB,D1,
     &          'N','N',IPQTOQ(1,0),TMAT,NBAS,NORB,NZT,
     &          'N','N',IPQTOQ(1,0),CBUF(1,1,1), NORB,NORB,NZ,
     &          D0,IPQTOQ(1,0),CMO(1,1,1),NBAS,NORB,NZ)
      IF(IPRINT.GE.5) THEN
        WRITE(6,*) '* RKBIMP: UKB complement..'
        CALL PRQMAT(CMO(1,1,1),NBAS,II,NBAS,NORB,
     &              NZ,IPQTOQ(1,0),LUPRI)
      ENDIF
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck getiopen */
      INTEGER FUNCTION GETIOPEN(IORBITAL,IFRP)
C***********************************************************************
C     
C     Which open shell does electron orbital IORBITAL belong to?
C
C     Written by J. Thyssen - Feb 11 1999
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
#include "dcbdhf.h"
#include "dcborb.h"
C
      IF (IORBITAL .LT. 1)
     &   CALL QUIT('illegal value in GETIOPEN (<1)')
      IF (IFRP .GT. NFSYM)
     &   CALL QUIT('illegal fermion symmetry in GETIOPEN')
C
      IF (IORBITAL .LE. NISH(IFRP)) THEN
         GETIOPEN = 0
         RETURN
      END IF
C
      ILOW = NISH(IFRP)+1
      DO IOPEN = 1,NOPEN
      IF(NACSH(IFRP,IOPEN).GT.0) THEN
         IHIGH = ILOW + NACSH(IFRP,IOPEN)
         IF ( (IORBITAL .GE. ILOW) .AND. (IORBITAL .LE. IHIGH) ) THEN
            GETIOPEN = IOPEN
            RETURN
         END IF
         ILOW = IHIGH+1
      ENDIF
      END DO
      CALL QUIT('illegal value in GETIOPEN (unoccupied orbital)')
      RETURN
      END

      SUBROUTINE Read_and_add_FCK(IUNIT,FOCK,INFO,N)
C*****************************************************************************
C
C     Add AO Fock matrix read from unformatted file to contents of Fock
C
C*****************************************************************************
C
#include "dcbbas.h"
      INTEGER IUNIT,N
      LOGICAL INFO
      Real*8 FOCK(N2BBASXQ*N)
      Real*8, Allocatable :: Buffer(:)

      Allocate (Buffer(N2BBASXQ*N))
      call REAFCK(IUNIT,Buffer,INFO,N)
      call DAXPY (N2BBASXQ*N,1.D0,Buffer,1,FOCK,1)
      DeAllocate (Buffer)

      End
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C  /* Deck getacord */
      SUBROUTINE GETACORD(COORA)
C*****************************************************************************
C
C    GETACORD : Make list atomic coordinates
C
C               Written oct.2001 by Jesper Kielberg Pedersen
C
C*****************************************************************************
#include "implicit.h"
#include "mxcent.h"
#include "maxorb.h"
#include "maxaqn.h"
#include "nuclei.h"
#include "symmet.h"
#include "pgroup.h"
      PARAMETER (D0 = 0.0D0 , DP5 = 0.5D0 , D2 = 2.0D0)
      DIMENSION COORA(3,*)
      CALL QENTER('GETACORD')
C
C     Make the full matrix of cartisian coordinates from CORD(NUCIND)
C
      JATOM = 0
      DO 100 ICENT = 1, NUCIND
         MULCNT = ISTBNU(ICENT)
         IF (MULT(MULCNT) .EQ. 1) THEN
            JATOM = JATOM + 1
            COORA(1,JATOM) = CORD(1,ICENT)
            COORA(2,JATOM) = CORD(2,ICENT)
            COORA(3,JATOM) = CORD(3,ICENT)
         ELSE
            DO 200 ISYMOP = 0, MAXOPR
               IF (IAND(ISYMOP,MULCNT) .EQ. 0) THEN
                  JATOM = JATOM + 1
                  COORA(1,JATOM) = 
     &                 PT(IAND(ISYMAX(1,1),ISYMOP))*CORD(1,ICENT)
                  COORA(2,JATOM) = 
     &                 PT(IAND(ISYMAX(2,1),ISYMOP))*CORD(2,ICENT)
                  COORA(3,JATOM) = 
     &                 PT(IAND(ISYMAX(3,1),ISYMOP))*CORD(3,ICENT)
               END IF
  200       CONTINUE
         END IF
  100 CONTINUE
      CALL QEXIT('GETACORD')
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C /* Deck Selcmo */
      SUBROUTINE SELCMO(ILABDF,IFRP,CMO,IVEC,NVEC,NPVEC,NEVEC,
     &                  CSEL,LRSEL,LCSEL,NCBAS,
     &                  IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     Select coefficients from CMO according to pointer index IVEC
C       ILABDF = 1: AO-basis
C       ILABDF = 2: SO-basis
C     The selected coefficients are stored in CSEL and returned
C     together with NCBAS which is its row-dimension.
C
C     Written by T. Saue May 30 2012
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
      DIMENSION CMO(*),CSEL(LRSEL,LCSEL,4),IVEC(NVEC),WORK(*)
      IF(ILABDF.EQ.1) THEN
C
C     Extract/expand coefficients in AO-basis
C     =======================================
C
        CALL SELCMO_AO(CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &                 CSEL,LRSEL,LCSEL,NTBAS(0),IFRP,IVEC,
     &                 NVEC,NPVEC,NEVEC,IPRINT,WORK,KFREE,LFREE)
      ELSE
C
C     Extract/expand coefficients in SO-basis
C     =======================================
C
        CALL SELCMO_SO(CMO(ICMOQ(IFRP)+1),NFBAS(IFRP,0),NORB(IFRP),
     &                 CSEL,LRSEL,LCSEL,IFRP,
     &                 IVEC,NVEC,NPVEC,NEVEC,IPRINT,WORK,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Selcmo_ao */
      SUBROUTINE SELCMO_AO(CMO,NBAS1,NORB1,
     &                     CSEL,LRSEL,LCSEL,NBAS2,IFRP,
     &                     IVEC,NVEC,NPVEC,NEVEC,
     &                     IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C     Select coefficients according to IVEC and
C     extract/expand in AO-basis
C
C     Written by T. Saue May 30 2012
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
      DIMENSION CMO(NBAS1,NORB1,NZ),CSEL(LRSEL,LCSEL,4),
     &          IVEC(NVEC),WORK(*)
      IF(NVEC.EQ.NORB1) THEN
C.....All coefficients selected
        IF(NBSYM.EQ.1) THEN
C.......C1 case: just copy
          CALL QEXPAND(CSEL,CMO,NVEC,IFRP,NBAS1)
        ELSE
C.......Higher symmetry: Expand to quaternion format and transform to AO
          NCBF = NBAS1*NVEC*4
          CALL MEMGET('REAL',KCBF,NCBF,WORK,KFREE,LFREE)
          CALL QEXPAND(WORK(KCBF),CMO,NVEC,IFRP,NBAS1)
          CALL WTSOAO(WORK(KCBF),CSEL,IFRP,NBAS1,NBAS2,NVEC*4,IPRINT)
          CALL MEMREL('SELCMO.AOa',WORK,1,KCBF,KFREE,LFREE)
        ENDIF
      ELSEIF(NVEC.LT.NORB1) THEN
C.....Not all coefficients selected
        NCBF = NBAS2*NVEC*4
        CALL MEMGET('REAL',KCBF,NCBF,WORK,KFREE,LFREE)
        LCB = NPVEC+NEVEC
        IF(NBSYM.EQ.1) THEN
C.........C1 case: select coefficients and copy
          CALL SELCFS(CMO,IFRP,WORK(KCBF),
     &                LCB,IVEC,NPVEC,NEVEC,
     &                NBAS1,NORB1)
          CALL QEXPAND(CSEL,WORK(KCBF),NVEC,IFRP,NBAS1)
        ELSE
C.........Higher symmetry: select, expand and transform
          CALL SELCFS(CMO,IFRP,CSEL,LCB,IVEC,NPVEC,NEVEC,
     &                NBAS1,NORB1)
          CALL QEXPAND(WORK(KCBF),CSEL,NVEC,IFRP,NBAS1)
          CALL WTSOAO(WORK(KCBF),CSEL,IFRP,
     &                NBAS1,NBAS2,NVEC*4,IPRINT)
        ENDIF
        CALL MEMREL('SELCMO.AOb',WORK,1,KCBF,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Selcmo_ao */
      SUBROUTINE SELCMO_SO(CMO,NBAS,NORB,CSEL,LRSEL,LCSEL,IFRP,
     &                     IVEC,NVEC,NPVEC,NEVEC,
     &                     IPRINT,WORK,KFREE,LFREE)
C***********************************************************************
C
C     Select coefficients according to IVEC and
C     extract/expand in SO-basis
C
C     Written by T. Saue May 30 2012
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"
      DIMENSION CMO(NBAS,NORB,NZ),CSEL(LRSEL,LCSEL,4),IVEC(NVEC),WORK(*)
      IF(NVEC.EQ.NORB) THEN
        CALL QEXPAND(CSEL,CMO,NVEC,IFRP,NBAS)
      ELSE
        NCBF  = NBAS*NVEC*NZ
        CALL MEMGET('REAL',KCBF,NCBF,WORK,KFREE,LFREE)
        LCB = NPVEC+NEVEC
        CALL SELCFS(CMO,IFRP,WORK(KCBF),LCB,IVEC,NPVEC,NEVEC,
     &              NBAS,NORB)
        CALL QEXPAND(CSEL,WORK(KCBF),NVEC,IFRP,NBAS)
        CALL MEMREL('SELCMO_SO',WORK,1,KCBF,KFREE,LFREE)
      ENDIF
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      subroutine resort_mdirac_vmat(tmat,vmat,nrow,ncol,ifsym,
     &                              nboson_orb,nz)
      implicit none
      real(8), intent(inout) :: tmat(nrow,ncol,nz)
      real(8), intent(inout) :: vmat(nrow,ncol,nz)
      integer, intent(in)    :: nrow
      integer, intent(in)    :: ncol
      integer, intent(in)    :: nz
      integer, intent(in)    :: nboson_orb(4,2,0:2)
      integer, intent(in)    :: ifsym
      integer                :: i, j, k, l, nv
      integer                :: nsym
      integer                :: nvec(2)
      integer                :: ioff(2)
      integer                :: isum

      call dzero(tmat,nrow*ncol*nz)
      call izero(nvec,2)
      call izero(ioff,2)

      nsym = 4 / nz ! # of boson irreps in gerade/ungerade

      ioff(1) = 0
      ioff(2) = isum(nsym,nboson_orb(1,ifsym,1),1)

      nv      = 0

      do j = 1, nsym
        do k = 2,1,-1
          do l = 1, nboson_orb(j,ifsym,k)
            nv = nv + 1
            do i = 1, nz
!                        copy from p1 p2 p3 p4 n1 n2 n3 n4 order (vmat)
!                        copy to   n1 p1 n2 p2 n3 p3 n4 p4 order (tmat)
!                        where n == negative energy shell
!                        where p == positive energy shell
!                        and   "1", "2", ... refer to boson irreps
              call dcopy(nrow,
     &                   vmat(1,ioff(k)+nvec(k)+l,i),1,
     &                   tmat(1,nv,i),1)
            end do
          end do
          nvec(k) = nvec(k) + nboson_orb(j,ifsym,k)
        end do
      end do
      end

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Setocc */
      SUBROUTINE SETOCC(NMOL,OCC,JMOL)
C***********************************************************************
C
C     
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbdhf.h"
#include "dcborb.h"
      DIMENSION OCC(NMOL),JMOL(3,NMOL)
C     Consistency check
      CALL DZERO(OCC,NMOL)
C     Negative-energy orbitals have zero occupation
C     Check positive-energy orbitals
      DO I = 1,NMOL
        IF(JMOL(2,I).EQ.1) THEN
          II    = JMOL(1,I)
          IFRP  = JMOL(3,I)
          NVORB = NISH(IFRP)
          IF(II.LE.NVORB) THEN
            OCC(I) = DF(0)
          ELSE
            DO IOPEN = 1,NOPEN
              NVORB =  NVORB + NACSH(IFRP,IOPEN)
              IF(II.LE.NVORB) THEN
                OCC(I) = DF(IOPEN)
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDDO
C
      RETURN
      END

C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck C1coef */
      SUBROUTINE C1COEF(CSO,CAO,NVEC,IOFF,NTOT,IOPT)
C***********************************************************************
C
C     Transform molecular coefficients CSO to CAO (C1 symmetry)
C     IOPT = 0 : transform to AO-basis
C     IOPT = 1 : transform to C1-basis
C
C     NOTE: 
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dcbbas.h"
#include "dcborb.h"
#include "dcbham.h"
#include "dgroup.h"
      DIMENSION CAO(*),CSO(*),NVEC(2)
      real(8), allocatable :: CBUF(:)
      integer, allocatable :: IBUF(:)
      allocate(IBUF(NTBAS(0)))
      IF(IOPT.EQ.1) THEN
        CALL NUCORD(IBUF)
      ELSE
        DO I = 1,NTBAS(0)
          IBUF(I) = I
        ENDDO
      ENDIF   
C
C     Transform to AO-basis --> stored in CAO
C     ========================================
C
C     NBSYM.EQ.1, implies C1 !!! Check how to handle that...
      ICAO  = 1+IOFF*NTBAS(0)
      ICSO  = 1
      IF(NZC1.EQ.1) THEN 
        DO IFRP = 1,NFSYM
        IF(NVEC(IFRP).GT.0) THEN
          IZOFF = ICAO
          CALL WTSOAO_mod(CSO(ICSO),CAO(IZOFF),IFRP,
     &                  NFBAS(IFRP,0),NTBAS(0),NVEC(IFRP),IBUF,0)
          ICAO  = ICAO  + NTBAS(0)*NVEC(IFRP)
          ICSO  = ICSO  + NFBAS(IFRP,0)*NVEC(IFRP)*NZ
        ENDIF
        ENDDO
      ELSE
        DO IFRP = 1,NFSYM
        IF(NVEC(IFRP).GT.0) THEN
          NCBF2 = NFBAS(IFRP,0)*NVEC(IFRP)*NZC1
          allocate(CBUF(NCBF2))
          CALL QEXPAND(CBUF,CSO(ICSO),NVEC(IFRP),IFRP,NFBAS(IFRP,0))
          IZOFF = ICAO
          ISOFF = 1
          DO IZ = 1,NZC1
            CALL WTSOAO_mod(CBUF(ISOFF),CAO(IZOFF),IFRP,
     &                  NFBAS(IFRP,0),NTBAS(0),NVEC(IFRP),IBUF,0)
            IZOFF = IZOFF + NTBAS(0)*NTOT
            ISOFF = ISOFF + NFBAS(IFRP,0)*NVEC(IFRP)
          ENDDO
          deallocate(CBUF)
          ICAO  = ICAO  + NTBAS(0)*NVEC(IFRP)
          ICSO  = ICSO  + NFBAS(IFRP,0)*NVEC(IFRP)*NZ
        ENDIF
        ENDDO
      ENDIF
      deallocate(IBUF)
C
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C/* Deck Nucord */
      SUBROUTINE NUCORD(IBUF)
C***********************************************************************
C
C     Find correct ordering of nuclei.
C     Written by Trond Saue
C
C***********************************************************************C
#include "implicit.h"
#include "mxcent.h"
      DIMENSION IBUF(*)
      integer, allocatable :: INUC(:)
#include "dcbham.h"
#include "dcbbas.h"
#include "nuclei.h"
C
      allocate(INUC(NTBAS(0)))
      CALL LSQMA1(INUC)
      I1 = 0
      I2 = 0
      DO IC = 1,MC
        DO ICENT = 1,NUCDEP
          DO IB = 1,NTBAS(IC)
          IF(INUC(I1+IB).EQ.ICENT) THEN
            I2 = I2 + 1
            IBUF(I1+IB) = I2
          ENDIF
          ENDDO
        ENDDO
        I1 = I1 + NTBAS(IC)
      ENDDO
      deallocate(INUC)
C
      RETURN
      END
      SUBROUTINE ZERO_MO(j,CMOBUF,NBAST,NMOT,NZ)
C
C     1-Oct-2015 Hans Joergen Aa. Jensen
C     Zero MO no. j in CMOBUF
C     (CMOBUF all of CMO, or a subset generated with SELCFS, or?)
C     In fact, this routine can zero column j in any quaternion matrix.
C
      integer j, nbast, nmot, nz
      real*8  cmobuf(nbast,nmot,nz)
      integer iz

      do iz = 1, nz
         CMOBUF(:,j,iz) = 0.0d0
      end do
      
      return
      end
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE FRAGINFO(INFOSEL,NVECT,KVEC,NFRAG,NSTR,
     &      WORK,KFREE,LFREE)
C***********************************************************************
C
C     Collect info about fragments in an integer array
C
C     Written by T. Saue Oct 5 2016
C
C***********************************************************************
#include "implicit.h"
#include "priunit.h"
C
#include "dgroup.h"      
      DIMENSION INFOSEL(2,NVECT),KVEC(2,NFRAG),NSTR(2,0:2,NFRAG),
     &     WORK(*)
      IOFF = 1
      DO IFRP = 1,NFSYM
         DO IFRAG = 1,NFRAG
            CALL IICOPY(NSTR(IFRP,0,IFRAG),WORK(KVEC(IFRP,IFRAG)),1,1,
     &           INFOSEL(1,IOFF),2,1)
            CALL ICOPY(NSTR(IFRP,0,IFRAG),IFRAG,0,
     &           INFOSEL(2,IOFF),2)
            IOFF = IOFF + NSTR(IFRP,0,IFRAG)
         ENDDO
      ENDDO
      RETURN
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck atomic_id */
      SUBROUTINE ATOMIC_ID(ID,KP,J,MJ,LL)
C***********************************************************************
C
C     Extract atomic quantum numbers from sub_block ID
C
C     m_j can be expressed as (2I-1)/2 * (-1)**(I-1)
C
C     The pair (I,J) is packed as a superindex IJ = J(J-1)/2 + I
C     where I = |kappa|
C     The final sub_block identification is then given as
C       ID = sgn(kappa)*IJ      
C
C     Written by Trond Saue Feb 25 2021
C
C***********************************************************************
      implicit none
      integer,intent(in)  :: ID
      integer,intent(out) :: KP,J,MJ,LL
      integer             :: INDI,INDJ
      INDJ = INT(SQRT(dble(2*ABS(ID)) + 0.25D0) + 0.4999D0)
      INDI = ABS(ID) - INDJ*(INDJ-1)/2
      KP   = INDJ * ID/ABS(ID)
      MJ   = (2*INDI - 1) * (-1)**(INDI+1)
      J    = 2*INDJ - 1
      IF(KP.GT.0)THEN
        LL = INDJ
      ELSE   
        LL = INDJ-1
      ENDIF
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck atomic_reacmo */
      SUBROUTINE ATOMIC_REACMO(IUNIT,FILNAM,CMO,EIG,IBEIG,TOTERG,IOPT)
C***********************************************************************
C
C     Wrapper routine to handle atomic supersymmetry; hope for better solution soon
C
C     Written by Trond Saue March 3 2021
C***********************************************************************
#include "implicit.h"
#include "dcborb.h"      
      DIMENSION CMO(*),EIG(*),IBEIG(NORBT)
      integer, allocatable :: IBEIG2(:,:)
      character            :: FILNAM*6
      allocate(IBEIG2(NORBT,2))
      IOPT2 = IOPT+16
      CALL REACMO(IUNIT,FILNAM,CMO,EIG,IBEIG2,TOTERG,IOPT2)
      CALL ATOMIC_ID_REDUCE(IBEIG2,IBEIG,NORBT)
      deallocate(IBEIG2)
      END     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck atomic_reacmo */
      SUBROUTINE ATOMIC_WRICMO(IUNIT,CMO,EIG,IBEIG,TOTERG)
C***********************************************************************
C
C     Wrapper routine to handle atomic supersymmetry; hope for better solution soon
C
C     Written by Trond Saue March 3 2021
C***********************************************************************
#include "implicit.h"
#include "dcborb.h"      
      DIMENSION CMO(*),EIG(*),IBEIG(NORBT)
      integer, allocatable :: IBEIG2(:,:)
      character            :: FILNAM*6
      allocate(IBEIG2(NORBT,2))
      DO I = 1,NORBT
         CALL ATOMIC_ID(IBEIG(I),IBEIG2(I,2),J,IBEIG2(I,1),LL)
      ENDDO
      CALL WRICMO(IUNIT,CMO,EIG,IBEIG2,TOTERG,.TRUE.)        
      deallocate(IBEIG2) 
      END     
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     /* Deck atomic_id_reduce */
      SUBROUTINE ATOMIC_ID_REDUCE(IBEIG2,IBEIG,N)
C***********************************************************************
C
C     Starting from a double atomic supersymmetry array, form
C     a single one
C
C     Written by Trond Saue Feb 25 2021      
C***********************************************************************
      implicit none
      integer,intent(in)  :: IBEIG2(N,2),N
      integer,intent(out) :: IBEIG(N)
      integer             :: I,KP,KPABS,KPOFF,MJABS,ISKP
      DO I = 1,N
        MJABS    = (ABS(IBEIG2(I,1))+1)/2
        KP       = IBEIG2(I,2)
        KPABS    = ABS(KP)
        ISKP     = KP/KPABS
        KPOFF    = KPABS*(KPABS-1)/2
        IBEIG(I) = ISKP*(KPOFF+MJABS)
      ENDDO
      END
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE CMOCLEAN(CMO,NBAS,NORB,NZ,FACTOR)
C***********************************************************************
C
C     Input is CMO-coefficients. 
C     The routine sets coefficients with norm below the 
C     threshold FACTOR to zero.
C
C     Written by Trond Saue Nov 13 2019
C
C***********************************************************************
      implicit none
      integer, intent(in)   :: NBAS,NORB,NZ
      real*8, intent(in)    :: FACTOR
      real*8, intent(inout) :: CMO(NBAS,NORB,NZ)
      real*8, allocatable   :: BUF(:)
      real*8                :: TMP
      integer               :: IORB,IBAS,IZ

      IF(NZ.EQ.1) THEN
        DO IORB = 1,NORB
          DO IBAS = 1,NBAS
          IF(ABS(CMO(IBAS,IORB,1)).LT.FACTOR) THEN
            CMO(IBAS,IORB,1) = 0.0D0
          ENDIF
          ENDDO
        ENDDO
      ELSE
        allocate(buf(nbas))
        DO IORB = 1,NORB
          BUF = 0.0D0
          DO IZ = 1,NZ
            DO IBAS = 1,NBAS
              BUF(IBAS) = BUF(IBAS) + CMO(IBAS,IORB,IZ)**2
            ENDDO
          ENDDO
          TMP = FACTOR*FACTOR
          DO IZ = 1,NZ
            DO IBAS = 1,NBAS
              IF(BUF(IBAS).LT.FACTOR) THEN
                CMO(IBAS,IORB,IZ)=0.0D0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
        deallocate(buf)
      ENDIF
      RETURN
      END
! -- end of dirgp.F --
